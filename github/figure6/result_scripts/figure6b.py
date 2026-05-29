from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_or_compute_crc_tf_scores, module_score, read_gene_by_sample_matrix, require_table, spearman_rows, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure6b")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores = load_or_compute_crc_tf_scores(args.input_dir)
    expression = read_gene_by_sample_matrix(args.input_dir / "tcga_expression_matrix.csv.gz")
    signatures = require_table(args.input_dir / "tumor_state_signatures.csv", "tumor_state_signatures")
    rows = []
    for sig, group in signatures.groupby("signature"):
        genes = group["gene"].astype(str).tolist()
        rows.append(module_score(expression, genes, zscore_genes=True).rename(sig))
    sig_scores = pd.concat(rows, axis=1).T
    out = spearman_rows(scores.T, sig_scores, "CRC_score", "signature")
    write_table(out, args.output_dir / "figure6b.csv")


if __name__ == "__main__":
    main()
