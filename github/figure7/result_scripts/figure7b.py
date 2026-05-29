from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import read_gene_by_sample_matrix, require_table, standardize_tf_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure7b")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    candidates = standardize_tf_table(require_table(args.input_dir / "crc_and_supporting_tf_table.csv", "crc_and_supporting_tf_table"))
    expr = read_gene_by_sample_matrix(args.input_dir / "malignant_scrna_expression_matrix.csv.gz")
    meta = require_table(args.input_dir / "malignant_cell_metadata.csv", "malignant_cell_metadata")
    cell_col = "barcode" if "barcode" in meta.columns else meta.columns[0]
    cluster_col = next((c for c in ["cell_cluster", "leiden", "cluster"] if c in meta.columns), None)
    subtype_col = next((c for c in ["subtype_context", "subtype"] if c in meta.columns), None)
    meta = meta.rename(columns={cell_col: "barcode", cluster_col: "cell_cluster", subtype_col: "subtype_context"})
    common = expr.columns.intersection(meta["barcode"])
    expr = expr[common]
    meta = meta.set_index("barcode").reindex(common)
    rows = []
    for gene in candidates["TF"]:
        if gene not in expr.index:
            continue
        values = pd.to_numeric(expr.loc[gene], errors="coerce")
        for (subtype, cluster), cells in meta.groupby(["subtype_context", "cell_cluster"]).groups.items():
            vals = values.loc[list(cells)]
            rows.append({"subtype_context": subtype, "cell_cluster": cluster, "TF": gene, "mean_ln_cpm_plus_1": float(vals.mean()), "pct_expressing": float(vals.gt(0).mean()), "n_cells": int(vals.shape[0])})
    out = pd.DataFrame(rows).merge(candidates, on="TF", how="left")
    write_table(out, args.output_dir / "figure7b.csv")


if __name__ == "__main__":
    main()
