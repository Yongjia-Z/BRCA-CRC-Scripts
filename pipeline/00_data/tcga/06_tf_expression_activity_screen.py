from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TCGA")
WORKSPACE = DATA_ROOT.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
EXPR = DATA / "expression"
RESULTS.mkdir(parents=True, exist_ok=True)

GENCODE_ANNO = DATA_ROOT / "GENCODE" / "results" / "gencode_v36_gene_annotation.tsv"
TF_LIST = DATA_ROOT / "TF_resources" / "results" / "Lambert_human_TF_list.tsv"
DOROTHEA = DATA_ROOT / "TF_resources" / "data" / "dorothea_hs_pancancer.tsv"


def sample_type_from_barcode(sample_id: str) -> str:
    code = sample_id[13:15] if len(sample_id) >= 15 else ""
    return {
        "01": "Primary Tumor",
        "02": "Recurrent Tumor",
        "06": "Metastatic",
        "11": "Solid Tissue Normal",
    }.get(code, f"Other_{code}")


def benjamini(pvalues: list[float]) -> np.ndarray:
    arr = np.array(pvalues, dtype=float)
    valid = np.isfinite(arr)
    out = np.full(arr.shape, np.nan)
    if valid.sum() > 0:
        out[valid] = multipletests(arr[valid], method="fdr_bh")[1]
    return out


def load_symbol_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm = pd.read_csv(EXPR / "TCGA-BRCA.star_tpm.tsv.gz", sep="\t")
    tpm["ensembl_gene_id"] = tpm["Ensembl_ID"].astype(str).str.replace(r"\..*$", "", regex=True)

    anno = pd.read_csv(GENCODE_ANNO, sep="\t")
    anno = anno[["ensembl_gene_id", "gene_name", "gene_type"]].drop_duplicates("ensembl_gene_id")

    merged = tpm.merge(anno, on="ensembl_gene_id", how="left")
    merged = merged.dropna(subset=["gene_name"])
    sample_cols = [c for c in tpm.columns if c not in ["Ensembl_ID", "ensembl_gene_id"]]

    expression = merged.groupby("gene_name", as_index=True)[sample_cols].mean()
    gene_anno = merged[["gene_name", "gene_type", "ensembl_gene_id"]].drop_duplicates("gene_name")
    return expression, gene_anno


def build_sample_annotation(samples: list[str]) -> pd.DataFrame:
    ann = pd.DataFrame(
        {
            "sample_id": samples,
            "patient_id": [s[:12] for s in samples],
            "sample_type": [sample_type_from_barcode(s) for s in samples],
        }
    )
    subtype_path = DATA / "cbioportal_brca_patient_clinical_wide.tsv"
    if subtype_path.exists():
        subtype = pd.read_csv(subtype_path, sep="\t")
        keep = ["patientId"] + [c for c in ["SUBTYPE", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS", "DSS_STATUS", "DSS_MONTHS", "PFS_STATUS", "PFS_MONTHS"] if c in subtype.columns]
        ann = ann.merge(subtype[keep], left_on="patient_id", right_on="patientId", how="left").drop(columns=["patientId"])
    return ann


def tf_expression_screen(expression: pd.DataFrame, sample_ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tf = pd.read_csv(TF_LIST, sep="\t")
    tf_symbols = sorted(set(tf["external_gene_name"].astype(str)) & set(expression.index.astype(str)))
    tf_expr = expression.loc[tf_symbols]

    primary_samples = sample_ann.loc[sample_ann["sample_type"] == "Primary Tumor", "sample_id"].tolist()
    normal_samples = sample_ann.loc[sample_ann["sample_type"] == "Solid Tissue Normal", "sample_id"].tolist()
    primary_samples = [s for s in primary_samples if s in tf_expr.columns]
    normal_samples = [s for s in normal_samples if s in tf_expr.columns]

    rows = []
    for gene, values in tf_expr.iterrows():
        primary = values[primary_samples].astype(float)
        normal = values[normal_samples].astype(float) if normal_samples else pd.Series(dtype=float)
        detected_fraction = float((primary > 1).mean()) if len(primary) else np.nan
        mean_primary = float(primary.mean()) if len(primary) else np.nan
        mean_normal = float(normal.mean()) if len(normal) else np.nan
        log2fc_tumor_normal = mean_primary - mean_normal if len(normal) else np.nan
        p_tn = np.nan
        if len(normal) >= 5 and len(primary) >= 5:
            p_tn = stats.mannwhitneyu(primary, normal, alternative="two-sided").pvalue

        subtype_p = np.nan
        subtype_means = {}
        if "SUBTYPE" in sample_ann.columns:
            groups = []
            for subtype, group_ann in sample_ann.loc[sample_ann["sample_type"] == "Primary Tumor"].groupby("SUBTYPE", dropna=True):
                cols = [s for s in group_ann["sample_id"].tolist() if s in tf_expr.columns]
                if len(cols) >= 5:
                    group_values = values[cols].astype(float)
                    groups.append(group_values)
                    subtype_means[f"mean_{subtype}"] = float(group_values.mean())
            if len(groups) >= 2:
                subtype_p = stats.kruskal(*groups).pvalue

        rows.append(
            {
                "TF": gene,
                "mean_primary_tpm_log": mean_primary,
                "mean_normal_tpm_log": mean_normal,
                "tumor_minus_normal": log2fc_tumor_normal,
                "primary_detected_fraction_expr_gt1": detected_fraction,
                "tumor_normal_p": p_tn,
                "subtype_kruskal_p": subtype_p,
                **subtype_means,
            }
        )

    summary = pd.DataFrame(rows)
    summary["tumor_normal_FDR"] = benjamini(summary["tumor_normal_p"].tolist())
    summary["subtype_kruskal_FDR"] = benjamini(summary["subtype_kruskal_p"].tolist())
    summary["detectable_in_BRCA"] = (summary["primary_detected_fraction_expr_gt1"] >= 0.2) & (summary["mean_primary_tpm_log"] > 1)
    summary = summary.sort_values(["detectable_in_BRCA", "subtype_kruskal_FDR", "mean_primary_tpm_log"], ascending=[False, True, False])
    return tf_expr, summary


def infer_tf_activity(expression: pd.DataFrame, sample_ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_samples = sample_ann.loc[sample_ann["sample_type"] == "Primary Tumor", "sample_id"].tolist()
    primary_samples = [s for s in primary_samples if s in expression.columns]
    expr = expression[primary_samples]

    dorothea = pd.read_csv(DOROTHEA, sep="\t")
    dorothea = dorothea[dorothea["target"].isin(expr.index)]
    dorothea = dorothea[dorothea["tf"].notna() & dorothea["target"].notna()]

    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0).fillna(0)

    activities = {}
    regulon_rows = []
    for tf, group in dorothea.groupby("tf"):
        group = group.drop_duplicates("target")
        targets = [g for g in group["target"].tolist() if g in z.index]
        if len(targets) < 10:
            continue
        mor = group.set_index("target").loc[targets, "mor"].astype(float)
        weighted = z.loc[targets].mul(mor, axis=0)
        activities[tf] = weighted.mean(axis=0)
        regulon_rows.append({"TF": tf, "n_targets_used": len(targets), "mean_mor": float(mor.mean())})

    activity = pd.DataFrame(activities).T
    regulon_summary = pd.DataFrame(regulon_rows).sort_values("n_targets_used", ascending=False)
    return activity, regulon_summary


def main() -> None:
    expression, gene_anno = load_symbol_expression()
    sample_ann = build_sample_annotation(expression.columns.tolist())
    sample_ann.to_csv(RESULTS / "TCGA_BRCA_expression_sample_annotation.tsv", sep="\t", index=False)
    gene_anno.to_csv(RESULTS / "TCGA_BRCA_expression_gene_annotation_matched.tsv", sep="\t", index=False)

    tf_expr, tf_summary = tf_expression_screen(expression, sample_ann)
    tf_expr.to_csv(RESULTS / "TCGA_BRCA_TF_expression_matrix.tsv.gz", sep="\t", compression="gzip")
    tf_summary.to_csv(RESULTS / "BRCA_detectable_TF_expression_summary.tsv", sep="\t", index=False)
    tf_summary.loc[tf_summary["detectable_in_BRCA"]].to_csv(RESULTS / "BRCA_detectable_TF_list.tsv", sep="\t", index=False)

    activity, regulon_summary = infer_tf_activity(expression, sample_ann)
    activity.to_csv(RESULTS / "TCGA_BRCA_DoRothEA_like_TF_activity_matrix.tsv.gz", sep="\t", compression="gzip")
    regulon_summary.to_csv(RESULTS / "TCGA_BRCA_DoRothEA_regulon_target_summary.tsv", sep="\t", index=False)

    activity_summary = []
    primary_ann = sample_ann[sample_ann["sample_id"].isin(activity.columns)].copy()
    if "SUBTYPE" in primary_ann.columns and not activity.empty:
        for tf, row in activity.iterrows():
            groups = []
            means = {}
            for subtype, g in primary_ann.groupby("SUBTYPE", dropna=True):
                cols = g["sample_id"].tolist()
                vals = row[cols].astype(float)
                if len(vals) >= 5:
                    groups.append(vals)
                    means[f"activity_mean_{subtype}"] = float(vals.mean())
            p = stats.kruskal(*groups).pvalue if len(groups) >= 2 else np.nan
            activity_summary.append({"TF": tf, "activity_subtype_kruskal_p": p, **means})
    activity_summary = pd.DataFrame(activity_summary)
    if not activity_summary.empty:
        activity_summary["activity_subtype_kruskal_FDR"] = benjamini(activity_summary["activity_subtype_kruskal_p"].tolist())
    activity_summary.to_csv(RESULTS / "TCGA_BRCA_TF_activity_subtype_summary.tsv", sep="\t", index=False)

    footprint_path = RESULTS / "ATAC_candidate_CRC_TF_footprint_expression_correlation.tsv"
    footprint = pd.read_csv(footprint_path, sep="\t") if footprint_path.exists() else pd.DataFrame()
    footprint_best = pd.DataFrame()
    if not footprint.empty:
        footprint_best = footprint.sort_values("flank_correlation", ascending=False).drop_duplicates("TF_Name")
        footprint_best = footprint_best.rename(columns={"TF_Name": "TF"})

    candidate = tf_summary.merge(activity_summary, on="TF", how="left") if not activity_summary.empty else tf_summary.copy()
    if not footprint_best.empty:
        candidate = candidate.merge(
            footprint_best[["TF", "TF_Motif_Name", "depth_correlation", "depth_FDR", "flank_correlation", "flank_FDR"]],
            on="TF",
            how="left",
        )
    candidate["expression_subtype_evidence"] = candidate["detectable_in_BRCA"] & (candidate["subtype_kruskal_FDR"] < 0.05)
    candidate["activity_subtype_evidence"] = candidate.get("activity_subtype_kruskal_FDR", pd.Series(np.nan, index=candidate.index)) < 0.05
    candidate["ATAC_footprint_expression_evidence"] = candidate.get("flank_FDR", pd.Series(np.nan, index=candidate.index)) < 0.05
    candidate["preliminary_CRC_evidence_count"] = (
        candidate["expression_subtype_evidence"].astype(int)
        + candidate["activity_subtype_evidence"].astype(int)
        + candidate["ATAC_footprint_expression_evidence"].astype(int)
    )
    candidate["candidate_CRC_entry"] = candidate["expression_subtype_evidence"] & (
        candidate["activity_subtype_evidence"] | candidate["ATAC_footprint_expression_evidence"]
    )
    candidate.to_csv(RESULTS / "BRCA_candidate_CRC_TF_initial_table.tsv", sep="\t", index=False)


    top_tfs = tf_summary.loc[tf_summary["detectable_in_BRCA"]].sort_values("subtype_kruskal_FDR").head(20)

    summary = {
        "n_expression_genes_mapped": int(expression.shape[0]),
        "n_expression_samples": int(expression.shape[1]),
        "n_lambert_tfs_in_expression": int(tf_expr.shape[0]),
        "n_detectable_tfs": int(tf_summary["detectable_in_BRCA"].sum()),
        "n_activity_tfs": int(activity.shape[0]),
        "n_candidate_crc_entry_tfs": int(candidate["candidate_CRC_entry"].sum()),
    }
    (RESULTS / "TF_expression_activity_screen_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "TF_expression_activity_screen_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
