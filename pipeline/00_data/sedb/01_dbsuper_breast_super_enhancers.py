from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))

try:
    import pyranges as pr
except ImportError:
    pr = None


ROOT = (DATA_ROOT / "SEdb")
WORKSPACE = DATA_ROOT.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURE4 = WORKSPACE / "data" / "figure4"
for d in [DATA, RESULTS, FIGURE4]:
    d.mkdir(parents=True, exist_ok=True)

GENCODE = DATA_ROOT / "GENCODE" / "results" / "gencode_v36_gene_annotation.tsv"
CRC = DATA_ROOT / "TCGA" / "results" / "CRC_member_high_confidence_TFs.tsv"
CHAIN = DATA_ROOT / "UCSC_LiftOver" / "data" / "hg19ToHg38.over.chain.gz"

CELL_LINES = ["MCF-7", "HCC1954"]


def load_tf_windows() -> pd.DataFrame:
    genes = pd.read_csv(GENCODE, sep="\t")
    tfs = pd.read_csv(CRC, sep="\t")["TF"].dropna().unique()
    sub = genes[genes["gene_name"].isin(tfs)].copy()
    sub["window_start"] = (sub["start"] - 100000).clip(lower=0)
    sub["window_end"] = sub["end"] + 100000
    return sub[["gene_name", "chromosome", "start", "end", "window_start", "window_end"]]


def read_bed(path: Path) -> pd.DataFrame:
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")) or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                rows.append(parts[:3])
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["chromosome", "start", "end"])
    df = df.iloc[:, :3]
    df.columns = ["chromosome", "start", "end"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna()
    df[["start", "end"]] = df[["start", "end"]].astype(int)
    return df


def load_chain_blocks() -> dict[str, list[dict[str, object]]]:
    blocks: dict[str, list[dict[str, object]]] = {}
    with gzip.open(CHAIN, "rt") as fh:
        current = None
        for raw in fh:
            line = raw.strip()
            if not line:
                current = None
                continue
            parts = line.split()
            if parts[0] == "chain":
                current = {
                    "tName": parts[2],
                    "tStart": int(parts[5]),
                    "qName": parts[7],
                    "qSize": int(parts[8]),
                    "qStrand": parts[9],
                    "qStart": int(parts[10]),
                    "tPos": int(parts[5]),
                    "qPos": int(parts[10]),
                }
                blocks.setdefault(current["tName"], [])
                continue
            if current is None:
                continue
            size = int(parts[0])
            t0 = current["tPos"]
            q0 = current["qPos"]
            blocks[current["tName"]].append(
                {
                    "t_start": t0,
                    "t_end": t0 + size,
                    "q_chrom": current["qName"],
                    "q_start": q0,
                    "q_size": current["qSize"],
                    "q_strand": current["qStrand"],
                }
            )
            if len(parts) >= 3:
                current["tPos"] += size + int(parts[1])
                current["qPos"] += size + int(parts[2])
    return blocks


def lift_bed_midpoint(bed: pd.DataFrame, blocks: dict[str, list[dict[str, object]]]) -> pd.DataFrame:
    rows = []
    for row in bed.itertuples(index=False):
        chrom = row.chromosome
        midpoint = int((int(row.start) + int(row.end)) / 2)
        half_width = max(1, int((int(row.end) - int(row.start)) / 2))
        for b in blocks.get(chrom, []):
            if int(b["t_start"]) <= midpoint < int(b["t_end"]):
                offset = midpoint - int(b["t_start"])
                if b["q_strand"] == "+":
                    mid = int(b["q_start"]) + offset
                else:
                    mid = int(b["q_size"]) - (int(b["q_start"]) + offset)
                rows.append({"chromosome": b["q_chrom"], "start": max(0, mid - half_width), "end": mid + half_width, "liftover_method": "midpoint"})
                break
    return pd.DataFrame(rows)


def overlap(se: pd.DataFrame, windows: pd.DataFrame, cell: str) -> pd.DataFrame:
    if se.empty or windows.empty:
        return pd.DataFrame()

    if pr is not None:
        se_pr = se.reset_index().rename(
            columns={"index": "se_index", "chromosome": "Chromosome", "start": "Start", "end": "End"}
        )[["Chromosome", "Start", "End", "se_index"]]
        win_pr = windows.rename(
            columns={"chromosome": "Chromosome", "window_start": "Start", "window_end": "End"}
        )[["Chromosome", "Start", "End", "gene_name", "start", "end"]].rename(
            columns={"start": "gene_start", "end": "gene_end"}
        )
        joined = pr.PyRanges(win_pr).join(pr.PyRanges(se_pr)).df
        if joined.empty:
            return pd.DataFrame()
        grouped = (
            joined.groupby(["gene_name", "Chromosome", "gene_start", "gene_end"])
            .agg(n_dbSUPER_SE_100kb=("se_index", "nunique"), nearest_overlap_start=("Start_b", "min"), nearest_overlap_end=("End_b", "min"))
            .reset_index()
        )
        grouped.insert(0, "cell_line", cell)
        return grouped.rename(columns={"gene_name": "TF", "Chromosome": "chromosome"})

    rows = []
    for chrom, win in windows.groupby("chromosome", sort=False):
        bed = se[se["chromosome"] == chrom].sort_values("start")
        if bed.empty:
            continue
        starts = bed["start"].to_numpy(dtype=np.int64)
        ends = bed["end"].to_numpy(dtype=np.int64)
        for w in win.itertuples(index=False):
            w_start = int(w.window_start)
            w_end = int(w.window_end)
            limit = np.searchsorted(starts, w_end, side="left")
            if limit == 0:
                continue
            hit_mask = ends[:limit] > w_start
            if not hit_mask.any():
                continue
            hit_idx = np.nonzero(hit_mask)[0]
            rows.append(
                {
                    "cell_line": cell,
                    "TF": w.gene_name,
                    "chromosome": chrom,
                    "gene_start": int(w.start),
                    "gene_end": int(w.end),
                    "n_dbSUPER_SE_100kb": int(hit_idx.size),
                    "nearest_overlap_start": int(starts[hit_idx][0]),
                    "nearest_overlap_end": int(ends[hit_idx][0]),
                    "interval_overlap_engine": "chromosome-grouped numpy.searchsorted fallback",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    windows = load_tf_windows()
    blocks = load_chain_blocks()
    all_rows = []
    summaries = []
    for cell in CELL_LINES:
        out = DATA / f"dbSUPER_hg19_{cell}.bed"
        ok = out.exists() and out.stat().st_size > 0
        if not ok:
            summaries.append({"cell_line": cell, "downloaded": False})
            continue
        bed_hg19 = read_bed(out)
        bed = lift_bed_midpoint(bed_hg19, blocks)
        bed.to_csv(RESULTS / f"dbSUPER_hg38_lifted_{cell}.bed.tsv", sep="\t", index=False)
        ov = overlap(bed, windows, cell)
        ov.to_csv(RESULTS / f"dbSUPER_{cell}_CRC_TF_100kb_overlap.tsv", sep="\t", index=False)
        if not ov.empty:
            all_rows.append(ov)
        summaries.append({
            "cell_line": cell,
            "downloaded": True,
            "n_SE_hg19": int(bed_hg19.shape[0]),
            "n_SE_lifted_hg38": int(bed.shape[0]),
            "n_TF_with_SE_100kb": int(ov["TF"].nunique()) if not ov.empty else 0,
            "interval_overlap_engine": "PyRanges" if pr is not None else "chromosome-grouped numpy.searchsorted fallback",
        })
    if all_rows:
        all_ov = pd.concat(all_rows, ignore_index=True)
        all_ov.to_csv(RESULTS / "dbSUPER_breast_CRC_TF_100kb_overlap_all.tsv", sep="\t", index=False)
        counts = all_ov.groupby("TF")["cell_line"].nunique().sort_values(ascending=False).head(30)
    (RESULTS / "dbSUPER_breast_SE_summary.json").write_text(json.dumps(summaries, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(summaries).to_csv(RESULTS / "dbSUPER_breast_SE_summary.tsv", sep="\t", index=False)
    print(json.dumps(summaries, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
