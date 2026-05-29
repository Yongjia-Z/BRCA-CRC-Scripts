from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_crc_tfs, require_table, standardize_tf_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure3c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    matrix = require_table(args.input_dir / "motif_burden_matrix.csv", "motif_burden_matrix")
    matrix = standardize_tf_table(matrix, ["TF_motif", "motif", "TF", "gene"])
    out = matrix[matrix["TF"].isin(set(crc["TF"]))].copy()
    missing = sorted(set(crc["TF"]) - set(out["TF"]))
    if missing:
        out = pd.concat([out, pd.DataFrame({"TF": missing, "omission_reason": "motif not present in burden matrix"})], ignore_index=True, sort=False)
    write_table(out, args.output_dir / "figure3c.csv")


if __name__ == "__main__":
    main()
