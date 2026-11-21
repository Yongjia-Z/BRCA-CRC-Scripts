from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TCGA")
RESULTS = ROOT / "results"


EVIDENCE_COLS = [
    "expression_evidence",
    "TF_activity_evidence",
    "regulon_evidence",
    "ATAC_motif_evidence",
    "ATAC_footprint_evidence",
    "ATAC_accessible_gene_evidence",
    "H3K27ac_enhancer_evidence",
    "super_enhancer_proxy_evidence",
    "TF_TF_edge_evidence",
]


def main() -> None:
    df = pd.read_csv(RESULTS / "BRCA_CRC_TF_with_enhancer_H3K27ac_evidence.tsv", sep="\t")

    score = pd.DataFrame({"TF": df["TF"]})
    score["expression_evidence"] = df["expression_subtype_evidence"].fillna(False).astype(bool)
    score["TF_activity_evidence"] = df["activity_subtype_evidence"].fillna(False).astype(bool)
    score["regulon_evidence"] = df["TF"].isin(pd.read_csv(RESULTS / "TCGA_BRCA_DoRothEA_regulon_target_summary.tsv", sep="\t")["TF"])
    score["ATAC_motif_evidence"] = df["ATAC_motif_enrichment_evidence"].fillna(False).astype(bool)
    score["ATAC_footprint_evidence"] = df["ATAC_footprint_expression_evidence"].fillna(False).astype(bool)
    score["ATAC_accessible_gene_evidence"] = (
        df["ATAC_peak_to_TF_gene_evidence"].fillna(False).astype(bool)
        | df["ATAC_nearest_TF_peak_evidence"].fillna(False).astype(bool)
    )
    score["H3K27ac_enhancer_evidence"] = df["ENCODE_MCF7_H3K27ac_gene_window_evidence"].fillna(False).astype(bool)

    enhancer_signal = df["max_enhancer_signal_100kb"].fillna(0)
    score["super_enhancer_proxy_evidence"] = enhancer_signal >= enhancer_signal.quantile(0.75)
    score["TF_TF_edge_evidence"] = df["CistromeCancer_TF_TF_edge_evidence"].fillna(False).astype(bool)

    score["CRC_evidence_score"] = score[EVIDENCE_COLS].astype(int).sum(axis=1)
    score["CRC_confidence"] = pd.cut(
        score["CRC_evidence_score"],
        bins=[-1, 3, 5, 9],
        labels=["low", "medium", "high"],
    ).astype(str)
    score["high_confidence_CRC_TF"] = score["CRC_evidence_score"] >= 6

    keep_meta = [
        "TF",
        "mean_primary_tpm_log",
        "subtype_kruskal_FDR",
        "activity_subtype_kruskal_FDR",
        "CancerTarget_total_TF_TF_degree",
        "max_enhancer_signal_100kb",
        "ENCODE_MCF7_H3K27ac_max_signal",
    ]
    out = score.merge(df[[c for c in keep_meta if c in df.columns]], on="TF", how="left")
    out = out.sort_values(["CRC_evidence_score", "CancerTarget_total_TF_TF_degree", "mean_primary_tpm_log"], ascending=[False, False, False])
    out.to_csv(RESULTS / "CRC_member_evidence_score_table.tsv", sep="\t", index=False)
    out[out["high_confidence_CRC_TF"]].to_csv(RESULTS / "CRC_member_high_confidence_TFs.tsv", sep="\t", index=False)

    summary = out["CRC_confidence"].value_counts().rename_axis("CRC_confidence").reset_index(name="n_TFs")
    summary.to_csv(RESULTS / "CRC_member_evidence_score_summary.tsv", sep="\t", index=False)


    top = out.head(30).iloc[::-1]

    log = {
        "n_total_TFs_scored": int(out.shape[0]),
        "n_high_confidence_CRC_TFs_score_ge6": int(out["high_confidence_CRC_TF"].sum()),
        "n_medium_confidence": int((out["CRC_confidence"] == "medium").sum()),
        "n_low_confidence": int((out["CRC_confidence"] == "low").sum()),
        "note": "super_enhancer_proxy_evidence uses top quartile CistromeCancer enhancer signal; replace with ROSE/LILY/SEdb evidence when available.",
    }
    (RESULTS / "CRC_member_evidence_score_summary.json").write_text(json.dumps(log, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(log, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
