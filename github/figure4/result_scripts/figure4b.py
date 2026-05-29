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
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure4b")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    promoter = require_table(args.input_dir / "direct_binding_promoter_overlap_long.csv", "direct_binding_promoter_overlap_long")
    locus = require_table(args.input_dir / "direct_binding_locus_overlap_long.csv", "direct_binding_locus_overlap_long")
    prom_tf = promoter.groupby("regulator_TF").agg(Breast_ReMap_promoter_binding_TF_count=("target_TF", "nunique"), Breast_ReMap_peak_overlaps=("n_peak_overlaps", "sum")).reset_index().rename(columns={"regulator_TF": "TF"})
    if "context" in promoter.columns:
        tnbc_prom = promoter[promoter["context"].astype(str).str.contains("TNBC|Basal", case=False, na=False)].groupby("regulator_TF").agg(TNBC_ChIP_promoter_peak_count=("n_peak_overlaps", "sum")).reset_index().rename(columns={"regulator_TF": "TF"})
    else:
        tnbc_prom = pd.DataFrame({"TF": []})
    locus_tf = locus.groupby("regulator_TF").agg(TNBC_ChIP_locus_peak_count=("n_peak_overlaps", "sum")).reset_index().rename(columns={"regulator_TF": "TF"})
    out = crc[["TF"]].merge(prom_tf, on="TF", how="left").merge(tnbc_prom, on="TF", how="left").merge(locus_tf, on="TF", how="left").fillna(0)
    write_table(out, args.output_dir / "figure4b.csv")


if __name__ == "__main__":
    main()
