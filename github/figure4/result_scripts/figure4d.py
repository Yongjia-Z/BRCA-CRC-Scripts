from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_crc_tfs, log_scale01, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure4d")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    h3 = require_table(args.input_dir / "h3k27ac_se_proxy_summary.csv", "h3k27ac_se_proxy_summary")
    h3 = h3[h3["TF"].isin(set(crc["TF"]))].copy()
    for raw, scaled in [("n_H3K27ac_peaks_100kb", "MCF7_H3K27ac_scaled"), ("n_ROSE_like_SE_100kb", "ROSE_like_SE_scaled"), ("n_dbSUPER_breast_cell_lines", "dbSUPER_scaled")]:
        h3[scaled] = log_scale01(h3[raw]) if raw in h3.columns else 0.0
    h3["H3K27ac_dot_size"] = h3["MCF7_H3K27ac_scaled"]
    h3["ROSE_dot_size"] = h3["ROSE_like_SE_scaled"]
    h3["dbSUPER_dot_size"] = h3["dbSUPER_scaled"]
    out = crc[["TF"]].merge(h3, on="TF", how="left")
    write_table(out, args.output_dir / "figure4d.csv")


if __name__ == "__main__":
    main()
