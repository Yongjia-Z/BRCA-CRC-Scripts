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
ENCODE_RESULTS = DATA_ROOT / "ENCODE" / "results"


def main() -> None:
    base = pd.read_csv(RESULTS / "BRCA_CRC_TF_multiomics_with_CistromeCancer_TF_TF.tsv", sep="\t")
    enhancer = pd.read_csv(CISTROME_RESULTS / "TCGA_enhancer_candidate_CRC_TF_gene_window_signal.tsv", sep="\t")
    h3k27ac = pd.read_csv(ENCODE_RESULTS / "ENCODE_MCF7_H3K27ac_candidate_TF_gene_window_summary.tsv", sep="\t")

    keep_enhancer = [
        "TF",
        "max_enhancer_signal_100kb",
        "sum_enhancer_signal_100kb",
        "CistromeCancer_enhancer_gene_window_evidence",
    ]
    keep_h3k27ac = [
        "TF",
        "ENCODE_MCF7_H3K27ac_files_with_overlap",
        "ENCODE_MCF7_H3K27ac_max_peaks_100kb",
        "ENCODE_MCF7_H3K27ac_max_signal",
        "ENCODE_MCF7_H3K27ac_gene_window_evidence",
    ]
    integrated = base.merge(enhancer[keep_enhancer], on="TF", how="left").merge(h3k27ac[keep_h3k27ac], on="TF", how="left")
    for col in ["CistromeCancer_enhancer_gene_window_evidence", "ENCODE_MCF7_H3K27ac_gene_window_evidence"]:
        integrated[col] = integrated[col].fillna(False).astype(bool)
    integrated["CRC_evidence_count_with_enhancer_H3K27ac"] = (
        integrated["multiomics_CRC_evidence_count_with_TF_TF"]
        + integrated["CistromeCancer_enhancer_gene_window_evidence"].astype(int)
        + integrated["ENCODE_MCF7_H3K27ac_gene_window_evidence"].astype(int)
    )
    integrated = integrated.sort_values(
        ["CRC_evidence_count_with_enhancer_H3K27ac", "CancerTarget_total_TF_TF_degree", "ENCODE_MCF7_H3K27ac_max_signal"],
        ascending=[False, False, False],
    )
    integrated.to_csv(RESULTS / "BRCA_CRC_TF_with_enhancer_H3K27ac_evidence.tsv", sep="\t", index=False)
    high = integrated[integrated["CRC_evidence_count_with_enhancer_H3K27ac"] >= 7].copy()
    high.to_csv(RESULTS / "BRCA_CRC_TF_high_evidence_with_enhancer_H3K27ac.tsv", sep="\t", index=False)

    plot_df = high.head(30).iloc[::-1]

    summary = {
        "n_candidate_tfs": int(integrated.shape[0]),
        "n_tfs_with_CistromeCancer_enhancer_signal": int(integrated["CistromeCancer_enhancer_gene_window_evidence"].sum()),
        "n_tfs_with_ENCODE_MCF7_H3K27ac_overlap": int(integrated["ENCODE_MCF7_H3K27ac_gene_window_evidence"].sum()),
        "n_high_evidence_tfs_evidence_ge7": int(high.shape[0]),
    }
    (RESULTS / "CRC_TF_enhancer_H3K27ac_integration_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "CRC_TF_enhancer_H3K27ac_integration_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
