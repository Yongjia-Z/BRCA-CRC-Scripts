from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy import stats
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TCGA")
RESULTS = ROOT / "results"

MODULES = {
    "luminal": ["FOXA1", "GATA3", "ESR1", "SPDEF", "RARA", "XBP1", "TFAP2C"],
    "her2": ["TFAP2C", "TFAP2A", "EHF", "ELF3", "STAT1", "MYC", "ERBB2"],
    "basal": ["TP63", "SOX9", "FOXC1", "RUNX1", "RUNX2", "CEBPB", "RELA", "JUN", "FOS"],
}


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, pd.NA), axis=0).fillna(0)


def score_module(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in matrix.index]
    if not present:
        return pd.Series(0.0, index=matrix.columns)
    return zscore_rows(matrix.loc[present]).mean(axis=0)


def main() -> None:
    activity = pd.read_csv(RESULTS / "TCGA_BRCA_DoRothEA_like_TF_activity_matrix.tsv.gz", sep="\t", index_col=0)
    tf_expr = pd.read_csv(RESULTS / "TCGA_BRCA_TF_expression_matrix.tsv.gz", sep="\t", index_col=0)
    sample_ann = pd.read_csv(RESULTS / "TCGA_BRCA_expression_sample_annotation.tsv", sep="\t")
    high = pd.read_csv(RESULTS / "CRC_member_high_confidence_TFs.tsv", sep="\t")
    high_tfs = high.loc[high["CRC_evidence_score"] >= 7, "TF"].tolist()

    samples = sorted(set(activity.columns) & set(tf_expr.columns))
    activity = activity[samples]
    tf_expr = tf_expr[samples]

    scores = pd.DataFrame({"sample_id": samples})
    scores["pan_CRC_TF_activity_score"] = score_module(activity, high_tfs).values
    scores["pan_CRC_TF_expression_score"] = score_module(tf_expr, high_tfs).values
    scores["pan_CRC_score"] = scores[["pan_CRC_TF_activity_score", "pan_CRC_TF_expression_score"]].mean(axis=1)
    for module, genes in MODULES.items():
        scores[f"{module}_CRC_TF_activity_score"] = score_module(activity, genes).values
        scores[f"{module}_CRC_TF_expression_score"] = score_module(tf_expr, genes).values
        scores[f"{module}_CRC_score"] = scores[[f"{module}_CRC_TF_activity_score", f"{module}_CRC_TF_expression_score"]].mean(axis=1)

    scores = scores.merge(sample_ann, on="sample_id", how="left")
    scores.to_csv(RESULTS / "CRC_activity_scores_by_sample.tsv", sep="\t", index=False)

    score_cols = [c for c in scores.columns if c.endswith("_CRC_score")]
    rows = []
    for score_col in score_cols:
        if "SUBTYPE" not in scores.columns:
            continue
        groups = []
        for subtype, g in scores.groupby("SUBTYPE", dropna=True):
            vals = g[score_col].dropna()
            if len(vals) >= 5:
                groups.append(vals)
                rows.append({"score": score_col, "SUBTYPE": subtype, "mean_score": float(vals.mean()), "median_score": float(vals.median()), "n": int(len(vals))})
        if len(groups) >= 2:
            p = stats.kruskal(*groups).pvalue
            rows.append({"score": score_col, "SUBTYPE": "__Kruskal_p__", "mean_score": p, "median_score": p, "n": 0})
    subtype_summary = pd.DataFrame(rows)
    subtype_summary.to_csv(RESULTS / "CRC_activity_score_subtype_summary.tsv", sep="\t", index=False)

    for score_col in ["luminal_CRC_score", "her2_CRC_score", "basal_CRC_score", "pan_CRC_score"]:
        if score_col not in scores.columns or "SUBTYPE" not in scores.columns:
            continue
        order = scores.groupby("SUBTYPE")[score_col].median().sort_values(ascending=False).index.tolist()
        data = [scores.loc[scores["SUBTYPE"] == subtype, score_col].dropna() for subtype in order]

    log = {
        "n_samples_scored": int(scores.shape[0]),
        "n_pan_crc_tfs_score_ge7": int(len(high_tfs)),
        "score_columns": score_cols,
    }
    (RESULTS / "CRC_activity_score_summary.json").write_text(json.dumps(log, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(log, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
