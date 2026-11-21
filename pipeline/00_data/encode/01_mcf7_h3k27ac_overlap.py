from __future__ import annotations

import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "ENCODE")
WORKSPACE = DATA_ROOT.parent
DATA = ROOT / "data" / "H3K27ac_MCF7"
RESULTS = ROOT / "results"
TCGA_RESULTS = DATA_ROOT / "TCGA" / "results"
CISTROME_RESULTS = DATA_ROOT / "CistromeCancer" / "results"
DATA.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)


def local_file(path: Path) -> dict:
    if not (path.exists() and path.stat().st_size > 0):
        raise FileNotFoundError(path)
    return {"status": "exists", "bytes": path.stat().st_size}


def parse_position(pos: str) -> tuple[str, int, int]:
    m = re.match(r"(chr[^:]+):(\d+)-(\d+)", str(pos))
    if not m:
        return "", -1, -1
    return m.group(1), int(m.group(2)), int(m.group(3))


def load_tf_windows() -> pd.DataFrame:
    pos = pd.read_csv(DATA_ROOT / "CistromeCancer" / "data" / "TCGA_enhancer" / "position.txt", sep="\t", header=None, names=["gene", "position"])
    candidate = pd.read_csv(TCGA_RESULTS / "BRCA_CRC_TF_high_evidence_with_TF_TF.tsv", sep="\t")
    pos = pos[pos["gene"].isin(candidate["TF"])].copy()
    parsed = pos["position"].apply(parse_position)
    pos["chrom"] = parsed.apply(lambda x: x[0])
    pos["start"] = parsed.apply(lambda x: x[1])
    pos["end"] = parsed.apply(lambda x: x[2])
    pos["window_start"] = (pos["start"] - 100000).clip(lower=0)
    pos["window_end"] = pos["end"] + 100000
    return pos


def count_overlaps(bed_gz: Path, windows: pd.DataFrame) -> pd.DataFrame:
    by_chrom = defaultdict(list)
    for row in windows.itertuples(index=False):
        by_chrom[row.chrom].append(row)
    stats = {row.gene: {"n_H3K27ac_peaks_100kb": 0, "max_peak_signal": 0.0, "sum_peak_signal": 0.0} for row in windows.itertuples(index=False)}
    with gzip.open(bed_gz, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            signal = float(parts[6]) if len(parts) > 6 else 0.0
            for row in by_chrom.get(chrom, []):
                if end < row.window_start or start > row.window_end:
                    continue
                rec = stats[row.gene]
                rec["n_H3K27ac_peaks_100kb"] += 1
                rec["max_peak_signal"] = max(rec["max_peak_signal"], signal)
                rec["sum_peak_signal"] += signal
    return pd.DataFrame([{"TF": gene, **rec} for gene, rec in stats.items()])


def main() -> None:
    files = pd.read_csv(REPO_ROOT / "data" / "dataset_lists" / "ENCODE_H3K27ac_breast_peak_files.tsv", sep="\t")
    chosen = files[
        (files["biosample_query"] == "MCF-7")
        & (files["assembly"] == "GRCh38")
        & (files["output_type"] == "replicated peaks")
    ].copy()
    chosen = chosen.sort_values(["preferred_default", "experiment_accession"], ascending=[False, True])
    rows = []
    for _, row in chosen.iterrows():
        filename = f"{row['file_accession']}.bed.gz"
        path = DATA / filename
        result = local_file(path)
        rows.append({**row.to_dict(), "local_path": str(path), **result})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(RESULTS / "ENCODE_MCF7_H3K27ac_peak_files_used.tsv", sep="\t", index=False)

    windows = load_tf_windows()
    all_stats = []
    for _, row in manifest.iterrows():
        stat = count_overlaps(Path(row["local_path"]), windows)
        stat["file_accession"] = row["file_accession"]
        stat["experiment_accession"] = row["experiment_accession"]
        all_stats.append(stat)
    combined = pd.concat(all_stats, ignore_index=True)
    combined.to_csv(RESULTS / "ENCODE_MCF7_H3K27ac_candidate_TF_gene_window_overlap_by_file.tsv", sep="\t", index=False)
    summary = (
        combined.groupby("TF")
        .agg(
            ENCODE_MCF7_H3K27ac_files_with_overlap=("n_H3K27ac_peaks_100kb", lambda x: int((x > 0).sum())),
            ENCODE_MCF7_H3K27ac_max_peaks_100kb=("n_H3K27ac_peaks_100kb", "max"),
            ENCODE_MCF7_H3K27ac_max_signal=("max_peak_signal", "max"),
            ENCODE_MCF7_H3K27ac_sum_signal=("sum_peak_signal", "sum"),
        )
        .reset_index()
    )
    summary["ENCODE_MCF7_H3K27ac_gene_window_evidence"] = summary["ENCODE_MCF7_H3K27ac_files_with_overlap"] > 0
    summary.to_csv(RESULTS / "ENCODE_MCF7_H3K27ac_candidate_TF_gene_window_summary.tsv", sep="\t", index=False)

    plot_df = summary.sort_values("ENCODE_MCF7_H3K27ac_max_signal", ascending=False).head(30).iloc[::-1]

    out = {
        "n_peak_files": int(manifest.shape[0]),
        "n_candidate_tfs": int(windows.shape[0]),
        "n_candidate_tfs_with_H3K27ac_overlap": int(summary["ENCODE_MCF7_H3K27ac_gene_window_evidence"].sum()),
    }
    (RESULTS / "ENCODE_MCF7_H3K27ac_integration_summary.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([out]).to_csv(RESULTS / "ENCODE_MCF7_H3K27ac_integration_summary.tsv", sep="\t", index=False)
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
