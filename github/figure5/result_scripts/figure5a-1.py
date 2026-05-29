from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import read_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure5a-1")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--panel", default="figure5a-1")
    parser.add_argument("--bin-size", type=int, default=200)
    return parser.parse_args()


def summarize_bigwig(path, chrom, start, end, bin_size):
    try:
        import pyBigWig
    except Exception:
        return {"positive_bins": np.nan, "scale_99_5pct": np.nan, "notes": "pyBigWig unavailable"}
    bw = pyBigWig.open(str(path))
    bins = list(range(int(start), int(end), int(bin_size)))
    vals = []
    for left in bins:
        right = min(left + bin_size, int(end))
        value = bw.stats(str(chrom), int(left), int(right), type="mean")[0]
        vals.append(0.0 if value is None else float(value))
    bw.close()
    arr = np.asarray(vals, dtype=float)
    return {"positive_bins": int(np.sum(arr > 0)), "scale_99_5pct": float(np.nanpercentile(arr, 99.5)) if arr.size else np.nan, "notes": "bigWig summarized by fixed bins"}


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = read_table(args.input_dir / "track_manifest.csv")
    if manifest.empty:
        raise FileNotFoundError("track_manifest.csv is required")
    if "panel" in manifest.columns:
        manifest = manifest[manifest["panel"].astype(str).eq(args.panel)].copy()
    rows = []
    for rec in manifest.to_dict("records"):
        file_name = str(rec.get("source_bigwig", rec.get("file", "")))
        path = args.input_dir / file_name
        chrom = rec.get("chrom", rec.get("chr"))
        start = rec.get("start")
        end = rec.get("end")
        summary = summarize_bigwig(path, chrom, start, end, args.bin_size) if path.exists() and pd.notna(chrom) else {"positive_bins": np.nan, "scale_99_5pct": np.nan, "notes": "bigWig missing"}
        row = dict(rec)
        row.update(summary)
        rows.append(row)
    write_table(pd.DataFrame(rows), args.output_dir / "figure5a-1.csv")


if __name__ == "__main__":
    main()
