from __future__ import annotations

import argparse
import sys
from pathlib import Path
from itertools import permutations

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import edge_score_column, load_crc_tfs, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure2c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--min-edge-score", type=float, default=6.0)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    tf_set = set(crc["TF"])
    
    # Load edges and determine score column
    edges = require_table(args.input_dir / "tf_tf_edge_evidence.csv", "tf_tf_edge_evidence")
    edges, score_col = edge_score_column(edges)
    
    # Filter edges between formal TFs with score >= min_edge_score
    edges = edges[edges["regulator_TF"].isin(tf_set) & edges["target_TF"].isin(tf_set)].copy()
    edges[score_col] = pd.to_numeric(edges[score_col], errors="coerce").fillna(0.0)
    edges = edges[edges[score_col].ge(args.min_edge_score)].copy()
    
    edge_map = {}
    for _, row in edges.sort_values(score_col, ascending=False).iterrows():
        edge_map.setdefault((row["regulator_TF"], row["target_TF"]), row)

    def edge_score(pair: tuple[str, str]) -> float:
        return float(edge_map[pair][score_col]) if pair in edge_map else 0.0

    member_cols = {
        "Luminal": "Luminal_member" if "Luminal_member" in crc.columns else "Luminal_state",
        "HER2": "HER2_member" if "HER2_member" in crc.columns else "HER2_state",
        "Basal/TNBC": "Basal/TNBC_member" if "Basal/TNBC_member" in crc.columns else "Basal/TNBC_state"
    }
    member = crc.set_index("TF")
    contexts = ["Luminal", "HER2", "Basal/TNBC"]

    def motif_context(nodes: list[str]) -> str:
        counts = {}
        for ctx, col in member_cols.items():
            if col in member.columns:
                # Support both 1/0 member flags and P/S/empty state strings
                vals = member.loc[nodes, col].fillna(0)
                if vals.dtype == object or str(vals.dtype).startswith("category"):
                    counts[ctx] = int(vals.isin(["P", "S", "1", "1.0", 1]).sum())
                else:
                    counts[ctx] = int(pd.to_numeric(vals, errors="coerce").fillna(0).sum())
            else:
                counts[ctx] = 0
        return max(contexts, key=lambda c: (counts[c], -contexts.index(c)))

    evidence_cols = [
        ("TRRUST", "TRRUST_curated_edge_evidence"),
        ("Direct", "context_direct_binding"),
        ("CancerTarget", "edge_CancerTarget_RP_evidence"),
        ("Expr corr.", "edge_expression_correlation_evidence"),
        ("Target ATAC", "edge_target_ATAC_peak_to_gene_evidence"),
        ("Reg. motif", "edge_regulator_ATAC_motif_evidence"),
        ("Reg. footprint", "edge_regulator_ATAC_footprint_evidence"),
        ("Target H3K27ac", "edge_target_H3K27ac_evidence"),
        ("Activity", "edge_activity_context_evidence"),
    ]

    def evidence_mean(pairs: list[tuple[str, str]], col: str) -> float:
        vals = []
        for pair in pairs:
            if pair not in edge_map or col not in edge_map[pair].index:
                continue
            vals.append(float(pd.to_numeric(pd.Series([edge_map[pair][col]]), errors="coerce").fillna(0).iloc[0]))
        return float(np.mean(vals)) if vals else 0.0

    rows = []
    
    # 1. Reciprocal pairs
    seen_rec = set()
    for a, b in sorted(edge_map):
        if (b, a) in edge_map:
            key = tuple(sorted([a, b]))
            if key in seen_rec:
                continue
            seen_rec.add(key)
            pairs = [(a, b), (b, a)]
            rec = {
                "type": "reciprocal_pair",
                "motif": f"{key[0]}<->{key[1]}",
                "context": motif_context(list(key)),
                "motif_score": sum(edge_score(p) for p in pairs)
            }
            for lab, col in evidence_cols:
                rec[lab] = evidence_mean(pairs, col)
            rec.update({"edge_1": f"{pairs[0][0]}->{pairs[0][1]}", "edge_2": f"{pairs[1][0]}->{pairs[1][1]}", "edge_3": ""})
            rows.append(rec)

    # 2. Feed-forward loops
    for a, b, c in permutations(sorted(tf_set), 3):
        pairs = [(a, b), (b, c), (a, c)]
        if all(p in edge_map for p in pairs):
            rec = {
                "type": "feed_forward_loop",
                "motif": f"{a}->{b}->{c}; {a}->{c}",
                "context": motif_context([a, b, c]),
                "motif_score": sum(edge_score(p) for p in pairs)
            }
            for lab, col in evidence_cols:
                rec[lab] = evidence_mean(pairs, col)
            rec.update({"edge_1": f"{a}->{b}", "edge_2": f"{b}->{c}", "edge_3": f"{a}->{c}"})
            rows.append(rec)

    # 3. Three-node cycles
    seen_cycle = set()
    for a, b, c in permutations(sorted(tf_set), 3):
        pairs = [(a, b), (b, c), (c, a)]
        if not all(p in edge_map for p in pairs):
            continue
        rotations = [(a, b, c), (b, c, a), (c, a, b)]
        key = min(rotations)
        if key in seen_cycle:
            continue
        seen_cycle.add(key)
        rec = {
            "type": "three_node_cycle",
            "motif": f"{key[0]}->{key[1]}->{key[2]}->{key[0]}",
            "context": motif_context(list(key)),
            "motif_score": sum(edge_score(p) for p in [(key[0], key[1]), (key[1], key[2]), (key[2], key[0])])
        }
        for lab, col in evidence_cols:
            rec[lab] = evidence_mean([(key[0], key[1]), (key[1], key[2]), (key[2], key[0])], col)
        rec.update({"edge_1": f"{key[0]}->{key[1]}", "edge_2": f"{key[1]}->{key[2]}", "edge_3": f"{key[2]}->{key[0]}"})
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values(["type", "motif_score"], ascending=[True, False]).reset_index(drop=True)
    write_table(out, args.output_dir / "figure2c.csv")


if __name__ == "__main__":
    main()
