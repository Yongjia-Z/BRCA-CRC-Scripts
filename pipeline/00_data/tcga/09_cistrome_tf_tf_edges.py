from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TCGA")
WORKSPACE = DATA_ROOT.parent
RESULTS = ROOT / "results"
CISTROME_RESULTS = DATA_ROOT / "CistromeCancer" / "results"


def main() -> None:
    candidate = pd.read_csv(RESULTS / "BRCA_CRC_TF_multiomics_initial_evidence_table.tsv", sep="\t")
    edges = pd.read_csv(CISTROME_RESULTS / "CancerTarget_candidate_TF_TF_edges_strong.tsv", sep="\t")

    out_degree = edges.groupby("regulator_TF").size().rename("CancerTarget_out_degree").reset_index().rename(columns={"regulator_TF": "TF"})
    in_degree = edges.groupby("target_TF").size().rename("CancerTarget_in_degree").reset_index().rename(columns={"target_TF": "TF"})
    degree = out_degree.merge(in_degree, on="TF", how="outer").fillna(0)
    degree["CancerTarget_total_TF_TF_degree"] = degree["CancerTarget_out_degree"] + degree["CancerTarget_in_degree"]
    degree["CistromeCancer_TF_TF_edge_evidence"] = degree["CancerTarget_total_TF_TF_degree"] > 0

    integrated = candidate.merge(degree, on="TF", how="left")
    for col in ["CancerTarget_out_degree", "CancerTarget_in_degree", "CancerTarget_total_TF_TF_degree"]:
        integrated[col] = integrated[col].fillna(0).astype(int)
    integrated["CistromeCancer_TF_TF_edge_evidence"] = integrated["CistromeCancer_TF_TF_edge_evidence"].fillna(False)
    integrated["multiomics_CRC_evidence_count_with_TF_TF"] = integrated["multiomics_CRC_evidence_count"] + integrated["CistromeCancer_TF_TF_edge_evidence"].astype(int)
    integrated = integrated.sort_values(
        ["multiomics_CRC_evidence_count_with_TF_TF", "CancerTarget_total_TF_TF_degree", "mean_primary_tpm_log"],
        ascending=[False, False, False],
    )
    integrated.to_csv(RESULTS / "BRCA_CRC_TF_multiomics_with_CistromeCancer_TF_TF.tsv", sep="\t", index=False)

    high = integrated[integrated["multiomics_CRC_evidence_count_with_TF_TF"] >= 5].copy()
    high.to_csv(RESULTS / "BRCA_CRC_TF_high_evidence_with_TF_TF.tsv", sep="\t", index=False)

    plot_df = high.head(30).iloc[::-1]

    summary = {
        "n_candidate_tfs": int(integrated.shape[0]),
        "n_tfs_with_CistromeCancer_TF_TF_edge": int(integrated["CistromeCancer_TF_TF_edge_evidence"].sum()),
        "n_high_evidence_tfs_evidence_ge5": int(high.shape[0]),
        "n_strong_TF_TF_edges": int(edges.shape[0]),
    }
    (RESULTS / "CRC_TF_CistromeCancer_TF_TF_integration_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "CRC_TF_CistromeCancer_TF_TF_integration_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
