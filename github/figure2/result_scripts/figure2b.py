from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import edge_score_column, load_crc_tfs, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure2b")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--source-table-label", default="figure1-4/figure2/figure2b.tsv")
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    tf_set = set(crc["TF"])
    edge_path = args.input_dir / "tf_tf_edge_evidence.csv"
    if not edge_path.exists():
        edge_path = args.input_dir / "figure2b.tsv"
    edges = require_table(edge_path, "tf_tf_edge_evidence")
    edges, score_col = edge_score_column(edges)
    evidence_cols = [c for c in edges.columns if c in ["TRRUST_curated_edge_evidence", "context_direct_binding", "edge_CancerTarget_RP_evidence", "edge_expression_correlation_evidence", "edge_target_ATAC_peak_to_gene_evidence", "edge_regulator_ATAC_motif_evidence", "edge_regulator_ATAC_footprint_evidence", "edge_target_H3K27ac_evidence", "edge_activity_context_evidence"]]
    edges = edges[edges["regulator_TF"].isin(tf_set) & edges["target_TF"].isin(tf_set)].copy()
    edges[score_col] = pd.to_numeric(edges[score_col], errors="coerce").fillna(0.0)
    panel_meta = [c for c in ["edge_context", "breast_context_peak_overlaps", "TNBC_peak_overlaps"] if c in edges.columns]
    panel_ready = {"edge_label", "edge_context", "breast_context_peak_overlaps", "TNBC_peak_overlaps"}.issubset(edges.columns)
    if panel_ready and edges.shape[0] <= args.top_n:
        out = edges.head(args.top_n).copy()
    else:
        out = edges.sort_values([score_col, "regulator_TF", "target_TF"], ascending=[False, True, True]).head(args.top_n).copy()
    out["edge_label"] = out["regulator_TF"].astype(str) + "->" + out["target_TF"].astype(str)
    out.insert(0, "source_table", out.get("source_table", args.source_table_label))
    columns = ["source_table", "regulator_TF", "target_TF", "edge_label", "edge_context", score_col, "breast_context_peak_overlaps", "TNBC_peak_overlaps"] + evidence_cols
    columns = list(dict.fromkeys([c for c in columns if c in out.columns]))
    write_table(out[columns], args.output_dir / "figure2b.csv")


if __name__ == "__main__":
    main()
