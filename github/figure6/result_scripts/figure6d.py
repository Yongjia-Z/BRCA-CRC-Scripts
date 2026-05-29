from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import collapse_pam50_to_context, load_or_compute_crc_tf_scores, normalize_pam50, require_table, residualize_by_category, spearman_rows, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure6d")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ann = require_table(args.input_dir / "tcga_sample_annotation.csv", "tcga_sample_annotation").set_index("sample_id")
    scores = load_or_compute_crc_tf_scores(args.input_dir)
    mcp = require_table(args.input_dir / "tcga_mcp_counter_scores.csv", "tcga_mcp_counter_scores").set_index("MCPcounter_population").T
    common = scores.index.intersection(mcp.index).intersection(ann.index)
    subtype = ann.loc[common, "SUBTYPE" if "SUBTYPE" in ann.columns else "PAM50"]
    context = subtype.map(normalize_pam50).map(collapse_pam50_to_context)
    score_res = scores.loc[common].apply(lambda s: residualize_by_category(s, context), axis=0).T
    mcp_res = mcp.loc[common].apply(lambda s: residualize_by_category(s, context), axis=0).T
    out = spearman_rows(score_res, mcp_res, "CRC_score", "MCPcounter_population")
    write_table(out, args.output_dir / "figure6d.csv")


if __name__ == "__main__":
    main()
