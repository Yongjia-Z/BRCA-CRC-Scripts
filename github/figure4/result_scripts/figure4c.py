from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_crc_tfs, ordered_tfs, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure4c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    order = ordered_tfs(crc)
    overlap = require_table(args.input_dir / "direct_binding_promoter_overlap_long.csv", "direct_binding_promoter_overlap_long")
    overlap = overlap[overlap["regulator_TF"].isin(order) & overlap["target_TF"].isin(order)].copy()
    overlap["weighted_self_overlap_ln"] = np.log1p(pd.to_numeric(overlap["n_peak_overlaps"], errors="coerce").fillna(0))
    mat = overlap.pivot_table(index="regulator_TF", columns="target_TF", values="weighted_self_overlap_ln", aggfunc="sum", fill_value=0).reindex(index=order, columns=order).reset_index()
    for tf in order:
        mat.loc[mat["regulator_TF"].eq(tf), tf] = np.nan
    write_table(mat, args.output_dir / "figure4c.csv")


if __name__ == "__main__":
    main()
