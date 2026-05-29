from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import edge_score_column, load_crc_tfs, read_manual_layout, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure2a")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--min-edge-score", type=float, default=6.0)
    return parser.parse_args()


def pagerank(nodes, edges, damping=0.85, n_iter=100):
    nodes = list(nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    if n == 0:
        return {}
    mat = np.zeros((n, n), dtype=float)
    for s, t, w in edges:
        if s in idx and t in idx and s != t:
            mat[idx[t], idx[s]] += float(w)
    colsum = mat.sum(axis=0)
    for j in range(n):
        mat[:, j] = 1.0 / n if colsum[j] == 0 else mat[:, j] / colsum[j]
    pr = np.ones(n) / n
    for _ in range(n_iter):
        pr = (1 - damping) / n + damping * (mat @ pr)
    return dict(zip(nodes, pr))


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    tf_set = set(crc["TF"])
    edges = require_table(args.input_dir / "tf_tf_edge_evidence.csv", "tf_tf_edge_evidence")
    edges, score_col = edge_score_column(edges)
    edges = edges[edges["regulator_TF"].isin(tf_set) & edges["target_TF"].isin(tf_set)].copy()
    edges[score_col] = pd.to_numeric(edges[score_col], errors="coerce").fillna(0.0)
    selected_edges = edges[edges[score_col].ge(args.min_edge_score)].copy()
    pr = pagerank(tf_set, selected_edges[["regulator_TF", "target_TF", score_col]].itertuples(index=False, name=None))
    nodes = crc.copy()
    nodes["PageRank"] = nodes["TF"].map(pr).fillna(0.0)
    layout = read_manual_layout(args.input_dir / "manual_layout.json")
    if not layout.empty:
        nodes = nodes.merge(layout[["TF", "x", "y"]], on="TF", how="left")
    out = pd.concat([nodes.assign(record_type="node"), selected_edges.assign(record_type="edge")], ignore_index=True, sort=False)
    write_table(out, args.output_dir / "figure2a.csv")
    (args.output_dir / "figure2a_run_audit.json").write_text(json.dumps({"analysis": "directed CRC TF topology with PageRank from context-supported TF-TF edges", "min_edge_score": args.min_edge_score, "n_edges": int(selected_edges.shape[0])}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
