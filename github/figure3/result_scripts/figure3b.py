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
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure3b")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def normalize_peak_type(value):
    text = str(value).lower()
    if "promoter" in text or "tss" in text:
        return "Promoter"
    if "intron" in text:
        return "Intronic"
    if "exon" in text:
        return "Exonic"
    if "distal" in text or "intergenic" in text or "enhancer" in text:
        return "Distal"
    return "Other"


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    peaks = require_table(args.input_dir / "scatac_peak_annotation.csv", "scatac_peak_annotation")
    peaks = standardize_tf_table(peaks, ["TF", "nearestGene", "gene", "gene_name"])
    type_col = next((c for c in ["peakType", "annotation", "feature_type"] if c in peaks.columns), None)
    if type_col is None:
        raise ValueError("scatac_peak_annotation requires peakType/annotation column")
    peaks = peaks[peaks["TF"].isin(set(crc["TF"]))].copy()
    peaks["peakType"] = peaks[type_col].map(normalize_peak_type)
    long = peaks.groupby(["TF", "peakType"], as_index=False).size().rename(columns={"size": "n_peaks"})
    total = long.groupby("TF")["n_peaks"].sum().rename("total_n_peaks")
    long = long.merge(total, on="TF", how="left")
    long["percent"] = 100 * long["n_peaks"] / long["total_n_peaks"].replace(0, np.nan)
    wide = long.pivot_table(index="TF", columns="peakType", values="percent", fill_value=0).reset_index()
    out = long.merge(wide, on="TF", how="left")
    missing = sorted(set(crc["TF"]) - set(out["TF"].astype(str))) if "TF" in out.columns else crc["TF"].tolist()
    if missing:
        out = pd.concat([out, pd.DataFrame({"TF": missing, "omission_reason": "TF not present in scATAC peak annotation"})], ignore_index=True, sort=False)
    write_table(out, args.output_dir / "figure3b.csv")


if __name__ == "__main__":
    main()
