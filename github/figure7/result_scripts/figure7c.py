from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import module_score, read_gene_by_sample_matrix, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure7c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    expr = read_gene_by_sample_matrix(args.input_dir / "malignant_scrna_expression_matrix.csv.gz")
    meta = require_table(args.input_dir / "malignant_cell_metadata.csv", "malignant_cell_metadata")
    states = require_table(args.input_dir / "crc_state_signatures.csv", "crc_state_signatures")
    cell_col = "barcode" if "barcode" in meta.columns else meta.columns[0]
    meta = meta.rename(columns={cell_col: "barcode"}).set_index("barcode")
    common = expr.columns.intersection(meta.index)
    meta = meta.loc[common].copy()
    expr = expr[common]
    for state, group in states.groupby("state_id"):
        genes = group["gene"].astype(str).tolist()
        meta[state] = module_score(expr, genes, zscore_genes=False).reindex(common)
    score_cols = [c for c in states["state_id"].drop_duplicates() if c in meta.columns]
    meta["dominant_crc_state"] = meta[score_cols].idxmax(axis=1)
    meta["dominant_crc_state_score"] = meta[score_cols].max(axis=1)
    if "subtype" in meta.columns:
        cell_out = meta.reset_index().assign(record_type="cell_score")
        fraction = meta.groupby(["subtype", "dominant_crc_state"]).size().rename("n_cells").reset_index()
        fraction["fraction_within_subtype"] = fraction["n_cells"] / fraction.groupby("subtype")["n_cells"].transform("sum")
        out = pd.concat([cell_out, fraction.assign(record_type="subtype_fraction")], ignore_index=True, sort=False)
    else:
        out = meta.reset_index().assign(record_type="cell_score")
    write_table(out, args.output_dir / "figure7c.csv")


if __name__ == "__main__":
    main()
