from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "CistromeCancer")
WORKSPACE = DATA_ROOT.parent
DATA = ROOT / "data" / "CancerTarget"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

TCGA_RESULTS = DATA_ROOT / "TCGA" / "results"
TF_TABLE = TCGA_RESULTS / "BRCA_CRC_TF_multiomics_high_evidence_initial.tsv"


def read_target_file(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path)
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return df
    return df.rename(columns={df.columns[0]: "target_gene"})


def main() -> None:
    tf_table = pd.read_csv(TF_TABLE, sep="\t")
    high_tfs = tf_table["TF"].astype(str).head(80).tolist()
    tf_set = set(high_tfs)

    rows = []
    downloaded_tfs = sorted({p.name.split(".")[0] for p in DATA.glob("*.rp.csv")})
    for regulator in downloaded_tfs:
        rp = read_target_file(DATA / f"{regulator}.rp.csv")
        cor = read_target_file(DATA / f"{regulator}.cor.csv")
        if rp.empty:
            continue
        rp_cols = [c for c in ["BRCA_1", "BRCA_2"] if c in rp.columns]
        cor_cols = [c for c in ["BRCA_1", "BRCA_2"] if c in cor.columns]
        if not rp_cols:
            continue
        rp = rp[rp["target_gene"].astype(str).isin(tf_set)].copy()
        if rp.empty:
            continue
        cor_small = cor[["target_gene"] + cor_cols].copy() if (not cor.empty and cor_cols) else pd.DataFrame(columns=["target_gene"])
        merged = rp[["target_gene"] + rp_cols].merge(cor_small, on="target_gene", how="left", suffixes=("_rp", "_cor"))
        for _, row in merged.iterrows():
            target = str(row["target_gene"])
            if target == regulator:
                continue
            brca1_rp = float(row.get("BRCA_1_rp", row.get("BRCA_1", 0)) or 0)
            brca2_rp = float(row.get("BRCA_2_rp", row.get("BRCA_2", 0)) or 0)
            brca1_cor = float(row.get("BRCA_1_cor", 0) or 0)
            brca2_cor = float(row.get("BRCA_2_cor", 0) or 0)
            rows.append(
                {
                    "regulator_TF": regulator,
                    "target_TF": target,
                    "BRCA_1_regular_potential": brca1_rp,
                    "BRCA_2_regular_potential": brca2_rp,
                    "BRCA_1_expression_correlation": brca1_cor,
                    "BRCA_2_expression_correlation": brca2_cor,
                    "max_BRCA_regular_potential": max(brca1_rp, brca2_rp),
                    "max_BRCA_expression_correlation": max(brca1_cor, brca2_cor),
                }
            )

    edges = pd.DataFrame(rows)
    if edges.empty:
        edges.to_csv(RESULTS / "CancerTarget_candidate_TF_TF_edges.tsv", sep="\t", index=False)
        print(json.dumps({"n_edges": 0}, ensure_ascii=False))
        return

    edges["CancerTarget_TF_TF_edge_evidence"] = (edges["max_BRCA_regular_potential"] >= 0.5) & (edges["max_BRCA_expression_correlation"] >= 0.3)
    edges = edges.sort_values(["CancerTarget_TF_TF_edge_evidence", "max_BRCA_regular_potential", "max_BRCA_expression_correlation"], ascending=[False, False, False])
    edges.to_csv(RESULTS / "CancerTarget_candidate_TF_TF_edges.tsv", sep="\t", index=False)
    strong = edges[edges["CancerTarget_TF_TF_edge_evidence"]].copy()
    strong.to_csv(RESULTS / "CancerTarget_candidate_TF_TF_edges_strong.tsv", sep="\t", index=False)

    node_summary = pd.concat(
        [
            strong.groupby("regulator_TF").size().rename("out_degree"),
            strong.groupby("target_TF").size().rename("in_degree"),
        ],
        axis=1,
    ).fillna(0).reset_index().rename(columns={"index": "TF"})
    node_summary["total_degree"] = node_summary["out_degree"] + node_summary["in_degree"]
    node_summary = node_summary.sort_values("total_degree", ascending=False)
    node_summary.to_csv(RESULTS / "CancerTarget_TF_TF_edge_node_degree.tsv", sep="\t", index=False)

    plot_df = node_summary.head(25).iloc[::-1]

    summary = {
        "n_downloaded_regulator_tfs": len(downloaded_tfs),
        "n_candidate_tf_targets_considered": len(tf_set),
        "n_all_tf_tf_edges": int(edges.shape[0]),
        "n_strong_tf_tf_edges": int(strong.shape[0]),
        "n_nodes_in_strong_edges": int(node_summary.shape[0]),
    }
    (RESULTS / "CancerTarget_TF_TF_edge_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "CancerTarget_TF_TF_edge_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
