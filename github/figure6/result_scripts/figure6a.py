from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import collapse_pam50_to_context, crc_tf_module_scores, load_crc_tfs, normalize_pam50, read_gene_by_sample_matrix, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure6a")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    expression = read_gene_by_sample_matrix(args.input_dir / "tcga_expression_matrix.csv.gz")
    activity = read_gene_by_sample_matrix(args.input_dir / "tcga_tf_activity_matrix.csv.gz")
    annotation = require_table(args.input_dir / "tcga_sample_annotation.csv", "tcga_sample_annotation")
    sample_col = "sample_id" if "sample_id" in annotation.columns else annotation.columns[0]
    subtype_col = next((c for c in ["SUBTYPE", "PAM50", "pam50", "subtype"] if c in annotation.columns), None)
    annotation = annotation.rename(columns={sample_col: "sample_id", subtype_col: "SUBTYPE"})
    annotation["PAM50"] = annotation["SUBTYPE"].map(normalize_pam50)
    annotation["SUBTYPE"] = annotation["PAM50"].map(collapse_pam50_to_context)
    common = expression.columns.intersection(activity.columns).intersection(annotation["sample_id"])
    scores = crc_tf_module_scores(crc, expression[common], activity[common])
    rows = [scores[col].rename("value").reset_index().rename(columns={"index": "sample_id"}).assign(score_label=col) for col in scores.columns]
    out = pd.concat(rows, ignore_index=True).merge(annotation, on="sample_id", how="left")
    write_table(scores.reset_index().rename(columns={"index": "sample_id"}), args.output_dir / "tcga_crc_tf_scores_by_sample.csv")
    write_table(out, args.output_dir / "figure6a.csv")


if __name__ == "__main__":
    main()
