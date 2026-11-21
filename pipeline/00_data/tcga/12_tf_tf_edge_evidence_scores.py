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
    edges = pd.read_csv(CISTROME_RESULTS / "CancerTarget_candidate_TF_TF_edges.tsv", sep="\t")
    tf_score = pd.read_csv(RESULTS / "CRC_member_evidence_score_table.tsv", sep="\t")
    atac_peak_links = pd.read_csv(RESULTS / "ATAC_BRCA_peak_to_candidate_TF_gene_links.tsv", sep="\t")
    h3k27ac = pd.read_csv(DATA_ROOT / "ENCODE" / "results" / "ENCODE_MCF7_H3K27ac_candidate_TF_gene_window_summary.tsv", sep="\t")

    score_cols = ["TF", "ATAC_motif_evidence", "ATAC_footprint_evidence", "TF_activity_evidence", "H3K27ac_enhancer_evidence", "CRC_evidence_score"]
    reg_meta = tf_score[score_cols].add_prefix("regulator_").rename(columns={"regulator_TF": "regulator_TF"})
    tar_meta = tf_score[score_cols].add_prefix("target_").rename(columns={"target_TF": "target_TF"})
    scored = edges.merge(reg_meta, on="regulator_TF", how="left").merge(tar_meta, on="target_TF", how="left")

    target_peak_summary = (
        atac_peak_links.groupby("Linked_Gene")
        .agg(target_ATAC_peak_to_gene_links=("Peak_Name", "count"), target_ATAC_peak_min_FDR=("FDR_Post_CNA", "min"))
        .reset_index()
        .rename(columns={"Linked_Gene": "target_TF"})
    )
    scored = scored.merge(target_peak_summary, on="target_TF", how="left")
    scored["target_ATAC_peak_to_gene_links"] = scored["target_ATAC_peak_to_gene_links"].fillna(0).astype(int)

    scored["edge_CancerTarget_RP_evidence"] = scored["max_BRCA_regular_potential"] >= 0.5
    scored["edge_expression_correlation_evidence"] = scored["max_BRCA_expression_correlation"] >= 0.3
    scored["edge_target_ATAC_peak_to_gene_evidence"] = (scored["target_ATAC_peak_to_gene_links"] > 0) & (scored["target_ATAC_peak_min_FDR"].fillna(1) < 0.05)
    scored["edge_regulator_ATAC_motif_evidence"] = scored["regulator_ATAC_motif_evidence"].fillna(False).astype(bool)
    scored["edge_regulator_ATAC_footprint_evidence"] = scored["regulator_ATAC_footprint_evidence"].fillna(False).astype(bool)
    scored["edge_target_H3K27ac_evidence"] = scored["target_H3K27ac_enhancer_evidence"].fillna(False).astype(bool)
    scored["edge_activity_context_evidence"] = scored["regulator_TF_activity_evidence"].fillna(False).astype(bool) & scored["target_TF_activity_evidence"].fillna(False).astype(bool)

    edge_evidence_cols = [
        "edge_CancerTarget_RP_evidence",
        "edge_expression_correlation_evidence",
        "edge_target_ATAC_peak_to_gene_evidence",
        "edge_regulator_ATAC_motif_evidence",
        "edge_regulator_ATAC_footprint_evidence",
        "edge_target_H3K27ac_evidence",
        "edge_activity_context_evidence",
    ]
    scored["TF_TF_edge_evidence_score"] = scored[edge_evidence_cols].astype(int).sum(axis=1)
    scored["high_confidence_TF_TF_edge"] = scored["TF_TF_edge_evidence_score"] >= 5
    scored = scored.sort_values(["high_confidence_TF_TF_edge", "TF_TF_edge_evidence_score", "max_BRCA_regular_potential"], ascending=[False, False, False])
    scored.to_csv(RESULTS / "TF_TF_regulatory_edge_evidence_score_table.tsv", sep="\t", index=False)
    scored[scored["high_confidence_TF_TF_edge"]].to_csv(RESULTS / "TF_TF_regulatory_edge_high_confidence.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        {
            "edge_confidence": ["high_score_ge5", "medium_score_3_4", "low_score_lt3"],
            "n_edges": [
                int((scored["TF_TF_edge_evidence_score"] >= 5).sum()),
                int(((scored["TF_TF_edge_evidence_score"] >= 3) & (scored["TF_TF_edge_evidence_score"] < 5)).sum()),
                int((scored["TF_TF_edge_evidence_score"] < 3).sum()),
            ],
        }
    )
    summary.to_csv(RESULTS / "TF_TF_regulatory_edge_evidence_summary.tsv", sep="\t", index=False)


    log = {
        "n_edges_scored": int(scored.shape[0]),
        "n_high_confidence_edges_score_ge5": int((scored["TF_TF_edge_evidence_score"] >= 5).sum()),
    }
    (RESULTS / "TF_TF_regulatory_edge_evidence_summary.json").write_text(json.dumps(log, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(log, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
