from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_or_compute_crc_tf_scores, require_table, spearman_rows, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure6c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores = load_or_compute_crc_tf_scores(args.input_dir).T
    mcp = require_table(args.input_dir / "tcga_mcp_counter_scores.csv", "tcga_mcp_counter_scores").set_index("MCPcounter_population")
    out = spearman_rows(scores, mcp, "CRC_score", "MCPcounter_population")
    write_table(out, args.output_dir / "figure6c.csv")


if __name__ == "__main__":
    main()
