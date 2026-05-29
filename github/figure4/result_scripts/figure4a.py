from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_crc_tfs, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure4a")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    scores = require_table(args.input_dir / "figure1a_all_tf_scores.csv", "figure1a_all_tf_scores")
    scores["H3K27ac_SE_MAX"] = scores[["Luminal_H3K27ac_SE", "HER2_H3K27ac_SE", "Basal_TNBC_H3K27ac_SE"]].apply(pd.to_numeric, errors="coerce").max(axis=1)
    scores["Direct_binding_MAX"] = scores[["Breast_context_direct_binding", "TNBC_direct_binding"]].apply(pd.to_numeric, errors="coerce").max(axis=1)
    cols = ["TF", "Luminal_H3K27ac_SE", "HER2_H3K27ac_SE", "Basal_TNBC_H3K27ac_SE", "H3K27ac_SE_MAX", "Breast_context_direct_binding", "TNBC_direct_binding", "Direct_binding_MAX"]
    out = crc[["TF"]].merge(scores[[c for c in cols if c in scores.columns]], on="TF", how="left")
    write_table(out, args.output_dir / "figure4a.csv")


if __name__ == "__main__":
    main()
