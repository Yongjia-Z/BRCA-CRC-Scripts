from __future__ import annotations

import argparse
import gzip
import json
import re
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


SCRIPT = Path(__file__).resolve()
DATASET = (DATA_ROOT / "GSE85158_H3K27ac")
WORKSPACE = DATA_ROOT.parent
DATASETS = DATA_ROOT
DATA = DATASET / "data"
RESULTS = DATASET / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
PER_FILE_RESULTS = RESULTS / "GSE85158_file_level_TF_window_signal"
PER_FILE_RESULTS.mkdir(parents=True, exist_ok=True)

GTF_V19 = DATASETS / "GENCODE" / "data" / "gencode.v19.annotation.gtf.gz"
CHAIN_HG19_TO_HG38 = DATASETS / "UCSC_LiftOver" / "data" / "hg19ToHg38.over.chain.gz"
TF_SOURCE = DATASETS / "TCGA" / "results" / "BRCA_CRC_TF_with_enhancer_H3K27ac_evidence.tsv"

CELL_SUBTYPE = {
    "76NF2V": "normal_like",
    "MCF10A": "normal_like",
    "MCF7": "luminal",
    "ZR751": "luminal",
    "MB361": "HER2",
    "UACC812": "HER2",
    "SKBR3": "HER2",
    "AU565": "HER2",
    "HCC1954": "HER2",
    "MB231": "basal_TNBC",
    "MB436": "basal_TNBC",
    "MB468": "basal_TNBC",
    "HCC1937": "basal_TNBC",
}

WIG_RE = re.compile(r"GSM\d+_(?P<cell>.+)\.H3K27ac\.rep(?P<rep>\d+)\.wig\.gz$")
FIXED_STEP_RE = re.compile(r"fixedStep chrom=(\S+) start=(\d+) step=(\d+) span=(\d+)")


def read_tf_universe() -> list[str]:
    df = pd.read_csv(TF_SOURCE, sep="\t", usecols=["TF"])
    return sorted(df["TF"].astype(str).dropna().unique())


def parse_gtf_attributes(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in attr.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field or " " not in field:
            continue
        key, value = field.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def load_hg19_gene_windows(tf_list: list[str], flank: int = 100_000) -> pd.DataFrame:
    wanted = set(tf_list)
    rows = []
    with gzip.open(GTF_V19, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene = attrs.get("gene_name")
            if gene not in wanted:
                continue
            chrom = parts[0]
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            start0 = max(int(parts[3]) - 1 - flank, 0)
            end0 = int(parts[4]) + flank
            rows.append(
                {
                    "TF": gene,
                    "chrom_hg19": chrom,
                    "gene_start0_hg19": int(parts[3]) - 1,
                    "gene_end0_hg19": int(parts[4]),
                    "strand": parts[6],
                    "window_start0_hg19": start0,
                    "window_end0_hg19": end0,
                    "window_length_bp": end0 - start0,
                    "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
                }
            )
    windows = pd.DataFrame(rows).sort_values(["chrom_hg19", "window_start0_hg19", "TF"])
    windows = windows.drop_duplicates("TF", keep="first").reset_index(drop=True)
    missing = sorted(wanted.difference(windows["TF"].astype(str)))
    (RESULTS / "GSE85158_missing_TFs_in_GENCODE_v19.json").write_text(
        json.dumps({"n_missing": len(missing), "missing_TFs": missing}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return windows


def liftover_windows_hg19_to_hg38(windows: pd.DataFrame) -> pd.DataFrame:
    lo = LiftOver(str(CHAIN_HG19_TO_HG38))
    rows = []
    for rec in windows.itertuples(index=False):
        left = lo.convert_coordinate(rec.chrom_hg19, int(rec.window_start0_hg19))
        right = lo.convert_coordinate(rec.chrom_hg19, max(int(rec.window_end0_hg19) - 1, int(rec.window_start0_hg19)))
        if not left or not right:
            rows.append({**rec._asdict(), "chrom_hg38": "", "window_start0_hg38": np.nan, "window_end0_hg38": np.nan, "liftover_status": "unmapped"})
            continue
        l0, r0 = left[0], right[0]
        if l0[0] != r0[0]:
            rows.append({**rec._asdict(), "chrom_hg38": "", "window_start0_hg38": np.nan, "window_end0_hg38": np.nan, "liftover_status": "split_or_inconsistent"})
            continue
        start = min(int(l0[1]), int(r0[1]))
        end = max(int(l0[1]), int(r0[1])) + 1
        rows.append({**rec._asdict(), "chrom_hg38": l0[0], "window_start0_hg38": start, "window_end0_hg38": end, "liftover_status": "mapped"})
    return pd.DataFrame(rows)


def build_window_index(windows: pd.DataFrame) -> dict[str, dict[str, np.ndarray]]:
    index: dict[str, dict[str, np.ndarray]] = {}
    for chrom, sub in windows.groupby("chrom_hg19", sort=False):
        sub = sub.sort_values("window_start0_hg19").copy()
        index[chrom] = {
            "row_idx": sub.index.to_numpy(dtype=np.int64),
            "starts": sub["window_start0_hg19"].to_numpy(dtype=np.int64),
            "ends": sub["window_end0_hg19"].to_numpy(dtype=np.int64),
        }
    return index


def scan_one_wig(path: Path, windows: pd.DataFrame, index: dict[str, dict[str, np.ndarray]]) -> pd.DataFrame:
    match = WIG_RE.match(path.name)
    if not match:
        raise ValueError(f"Unexpected WIG filename: {path.name}")
    cell = match.group("cell")
    rep = int(match.group("rep"))
    subtype = CELL_SUBTYPE.get(cell, "unknown")

    n = windows.shape[0]
    signal_sum = np.zeros(n, dtype=np.float64)
    positive_bp = np.zeros(n, dtype=np.float64)
    max_signal = np.zeros(n, dtype=np.float64)
    n_signal_bins = np.zeros(n, dtype=np.int64)

    chrom = None
    pos1 = None
    step = None
    span = None
    chrom_index = None
    ptr = 0
    active: deque[int] = deque()
    n_lines = 0
    n_parsed_inside_windows = 0
    started = time.time()

    with gzip.open(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("track"):
                continue
            if line.startswith("fixedStep"):
                m = FIXED_STEP_RE.match(line)
                if not m:
                    raise ValueError(f"Malformed fixedStep line in {path.name}: {line}")
                chrom = m.group(1)
                pos1 = int(m.group(2))
                step = int(m.group(3))
                span = int(m.group(4))
                chrom_index = index.get(chrom)
                ptr = 0
                active = deque()
                continue

            if chrom_index is None or pos1 is None or step is None or span is None:
                continue
            start0 = pos1 - 1
            end0 = start0 + span
            pos1 += step
            n_lines += 1

            starts = chrom_index["starts"]
            ends = chrom_index["ends"]
            row_idx = chrom_index["row_idx"]
            while ptr < len(starts) and starts[ptr] < end0:
                if ends[ptr] > start0:
                    active.append(ptr)
                ptr += 1
            if active:
                while active and ends[active[0]] <= start0:
                    active.popleft()
            if not active:
                continue

            value = float(line)
            n_parsed_inside_windows += 1
            if value <= 0:
                continue
            for active_pos in active:
                global_row = int(row_idx[active_pos])
                overlap = min(end0, int(ends[active_pos])) - max(start0, int(starts[active_pos]))
                if overlap <= 0:
                    continue
                signal_sum[global_row] += value * overlap
                positive_bp[global_row] += overlap
                max_signal[global_row] = max(max_signal[global_row], value)
                n_signal_bins[global_row] += 1

    out = windows[["TF", "window_length_bp"]].copy()
    out["gsm_file"] = path.name
    out["cell_line"] = cell
    out["replicate"] = rep
    out["subtype_context"] = subtype
    out["signal_sum"] = signal_sum
    out["mean_signal_per_bp"] = signal_sum / out["window_length_bp"].replace(0, np.nan).astype(float)
    out["positive_bp"] = positive_bp
    out["positive_fraction"] = positive_bp / out["window_length_bp"].replace(0, np.nan).astype(float)
    out["max_signal"] = max_signal
    out["n_signal_bins"] = n_signal_bins
    out["elapsed_seconds"] = round(time.time() - started, 3)
    out["n_wig_bins_seen"] = n_lines
    out["n_wig_bins_inside_TF_windows"] = n_parsed_inside_windows
    return out


def quantify_or_load_one_file(args: tuple[str, pd.DataFrame, dict[str, dict[str, np.ndarray]], str]) -> dict[str, object]:
    path_str, windows, index, out_dir_str = args
    path = Path(path_str)
    out_dir = Path(out_dir_str)
    per_file_out = out_dir / f"{path.name.replace('.wig.gz', '')}.TF_window_signal.tsv"
    if per_file_out.exists():
        scanned = pd.read_csv(per_file_out, sep="\t")
        from_cache = True
    else:
        scanned = scan_one_wig(path, windows, index)
        scanned.to_csv(per_file_out, sep="\t", index=False)
        from_cache = False
    return {
        "file": path.name,
        "cell_line": str(scanned["cell_line"].iat[0]),
        "subtype_context": str(scanned["subtype_context"].iat[0]),
        "replicate": int(scanned["replicate"].iat[0]),
        "elapsed_seconds": float(scanned["elapsed_seconds"].iat[0]),
        "n_wig_bins_seen": int(scanned["n_wig_bins_seen"].iat[0]),
        "n_wig_bins_inside_TF_windows": int(scanned["n_wig_bins_inside_TF_windows"].iat[0]),
        "from_cache": from_cache,
        "output": str(per_file_out),
    }


def summarize(long_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_cell = (
        long_df.groupby(["TF", "cell_line", "subtype_context"], as_index=False)
        .agg(
            mean_signal_per_bp=("mean_signal_per_bp", "mean"),
            positive_fraction=("positive_fraction", "mean"),
            positive_bp=("positive_bp", "mean"),
            max_signal=("max_signal", "max"),
            n_replicates=("replicate", "nunique"),
        )
    )
    normal = (
        by_cell.loc[by_cell["subtype_context"].eq("normal_like")]
        .groupby("TF")["mean_signal_per_bp"]
        .max()
        .rename("normal_like_max_mean_signal")
    )
    by_subtype = (
        by_cell.groupby(["TF", "subtype_context"], as_index=False)
        .agg(
            subtype_max_mean_signal=("mean_signal_per_bp", "max"),
            subtype_max_positive_fraction=("positive_fraction", "max"),
            subtype_max_positive_bp=("positive_bp", "max"),
            subtype_max_signal=("max_signal", "max"),
            n_cell_lines=("cell_line", "nunique"),
        )
    )
    wide = by_subtype.pivot_table(index="TF", columns="subtype_context", values="subtype_max_mean_signal", fill_value=0)
    out = pd.DataFrame(index=wide.index)
    for context, label in [("luminal", "Luminal"), ("HER2", "HER2"), ("basal_TNBC", "Basal_TNBC"), ("normal_like", "Normal_like")]:
        out[f"GSE85158_{label}_H3K27ac_mean_signal"] = wide[context] if context in wide.columns else 0.0
    out = out.join(normal, how="left").fillna(0)
    for label in ["Luminal", "HER2", "Basal_TNBC"]:
        raw = out[f"GSE85158_{label}_H3K27ac_mean_signal"]
        out[f"GSE85158_{label}_minus_normal_mean_signal"] = np.maximum(raw - out["normal_like_max_mean_signal"], 0)
    return by_cell, out.reset_index()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit-files", type=int, default=0, help="Debug only: process the first N WIG files.")
    parser.add_argument("--max-workers", type=int, default=8, help="Number of WIG files to scan in parallel.")
    args = parser.parse_args()

    tf_list = read_tf_universe()
    windows = load_hg19_gene_windows(tf_list)
    windows_lifted = liftover_windows_hg19_to_hg38(windows)
    windows_lifted.to_csv(RESULTS / "GSE85158_TF_gene_windows_hg19_lifted_to_hg38.tsv", sep="\t", index=False)

    index = build_window_index(windows)
    wig_files = sorted(DATA.glob("GSM*_*.H3K27ac.rep*.wig.gz"))
    if args.limit_files and args.limit_files > 0:
        wig_files = wig_files[: args.limit_files]
    tasks = [(str(path), windows, index, str(PER_FILE_RESULTS)) for path in wig_files]
    manifest = []
    max_workers = max(1, int(args.max_workers))
    print(f"Processing {len(tasks)} WIG files with max_workers={max_workers}", flush=True)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(quantify_or_load_one_file, task): Path(task[0]).name for task in tasks}
        for done_i, future in enumerate(as_completed(futures), start=1):
            file_name = futures[future]
            result = future.result()
            manifest.append(result)
            cache_note = "cache" if result["from_cache"] else "new"
            print(f"[{done_i}/{len(tasks)}] finished {file_name} ({cache_note}, {result['elapsed_seconds']} sec)", flush=True)
    manifest = sorted(manifest, key=lambda x: x["file"])
    all_long = [pd.read_csv(row["output"], sep="\t") for row in manifest]
    long_df = pd.concat(all_long, ignore_index=True)
    long_df.to_csv(RESULTS / "GSE85158_TF_window_H3K27ac_signal_by_file.tsv", sep="\t", index=False)
    by_cell, subtype_summary = summarize(long_df)
    by_cell.to_csv(RESULTS / "GSE85158_TF_window_H3K27ac_signal_by_cell_line.tsv", sep="\t", index=False)
    subtype_summary.to_csv(RESULTS / "GSE85158_TF_window_H3K27ac_subtype_summary_hg19_to_hg38.tsv", sep="\t", index=False)
    pd.DataFrame(manifest).to_csv(RESULTS / "GSE85158_H3K27ac_quantification_manifest.tsv", sep="\t", index=False)
    (RESULTS / "GSE85158_H3K27ac_quantification_summary.json").write_text(
        json.dumps(
            {
                "n_tf_windows": int(windows.shape[0]),
                "n_wig_files_processed": len(wig_files),
                "coordinate_strategy": "GSE85158 WIG signal was quantified in native hg19 GENCODE v19 TF windows; TF windows were then lifted from hg19 to hg38 for audit and downstream integration.",
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(json.dumps({"n_files": len(wig_files), "n_tfs": int(windows.shape[0])}, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
