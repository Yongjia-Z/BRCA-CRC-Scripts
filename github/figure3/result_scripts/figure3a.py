from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from concurrent.futures import ProcessPoolExecutor, as_completed
from crc_repro import count_region_overlaps, gene_windows, load_crc_tfs, read_bed, read_table, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure3a")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--flank-bp", type=int, default=100000)
    parser.add_argument("--max-workers", type=int, default=10)
    return parser.parse_args()


def worker(row, base, regions):
    peaks = read_bed(Path(base) / str(row["file"]))
    out = count_region_overlaps(peaks, regions)
    out["peakset"] = row.get("peakset", row.get("file"))
    out["peakset_label"] = row.get("peakset_label", row.get("file"))
    return out


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    genes = crc["TF"].tolist()
    summary = read_table(args.input_dir / "scatac_locus_peak_count_summary.csv")
    if not summary.empty:
        out = summary[summary["TF"].isin(genes)].copy()
    else:
        annotation = require_table(args.input_dir / "gene_annotation.csv", "gene_annotation")
        regions = gene_windows(annotation, genes, args.flank_bp)
        manifest = require_table(args.input_dir / "scatac_peakset_manifest.csv", "scatac_peakset_manifest")
        frames = []
        with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
            futures = [pool.submit(worker, rec, str(args.input_dir), regions) for rec in manifest.to_dict("records")]
            for fut in as_completed(futures):
                frames.append(fut.result())
        out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["TF", "peakset", "n_peak_overlaps"])
    if "n_peak_overlaps" in out.columns and "n_peaks" not in out.columns:
        out = out.rename(columns={"n_peak_overlaps": "n_peaks"})
    missing = sorted(set(genes) - set(out["TF"].astype(str))) if "TF" in out.columns else genes
    if missing:
        out = pd.concat([out, pd.DataFrame({"TF": missing, "n_peaks": np.nan, "omission_reason": "TF not present in scATAC peak-count input"})], ignore_index=True, sort=False)
    out["log1p_n_peaks"] = np.log1p(pd.to_numeric(out["n_peaks"], errors="coerce").fillna(0))
    write_table(out, args.output_dir / "figure3a.csv")


if __name__ == "__main__":
    main()
