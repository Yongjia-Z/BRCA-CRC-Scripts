from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import genes_for_context, load_crc_tfs, module_score, ordered_tfs, read_gene_by_sample_matrix, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure7a")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    expr = read_gene_by_sample_matrix(args.input_dir / "scrna_expression_matrix.csv.gz")
    meta = require_table(args.input_dir / "scrna_cell_metadata.csv", "scrna_cell_metadata")
    cell_col = "barcode" if "barcode" in meta.columns else meta.columns[0]
    meta = meta.rename(columns={cell_col: "barcode"})
    common = expr.columns.intersection(meta["barcode"])
    expr = expr[common]
    meta = meta.set_index("barcode").reindex(common).reset_index()
    modules = {"Pan": ordered_tfs(crc), "Luminal": genes_for_context(crc, "Luminal"), "HER2": genes_for_context(crc, "HER2"), "Basal/TNBC": genes_for_context(crc, "Basal/TNBC")}
    for name, genes in modules.items():
        meta[f"{name}_CRC_TF_score"] = module_score(expr, genes, zscore_genes=False).reindex(common).to_numpy()
    write_table(meta, args.output_dir / "figure7a.csv")


if __name__ == "__main__":
    main()
