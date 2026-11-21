from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = DATA_ROOT.parent
ENCODE_DIR = DATA_ROOT / "ENCODE"
DATA_DIR = ENCODE_DIR / "data" / "H3K27ac_MCF7"
RESULTS_DIR = ENCODE_DIR / "results"
GENCODE_PATH = DATA_ROOT / "GENCODE" / "results" / "gencode_v36_gene_annotation.tsv"
CRC_TF_PATH = DATA_ROOT / "TCGA" / "results" / "CRC_member_evidence_score_table.tsv"
TCGA_RESULTS_DIR = DATA_ROOT / "TCGA" / "results"

STITCH_DISTANCE = 12_500
PROMOTER_EXCLUSION = 2_500
GENE_WINDOW = 100_000


def read_narrowpeak(path: Path) -> pd.DataFrame:
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak",
    ]
    with gzip.open(path, "rt") as handle:
        df = pd.read_csv(handle, sep="\t", header=None, names=cols)
    df["source_file"] = path.name
    return df


def load_all_peaks() -> pd.DataFrame:
    files = sorted(DATA_DIR.glob("*.bed.gz"))
    if not files:
        raise FileNotFoundError(f"No H3K27ac narrowPeak files found in {DATA_DIR}")
    peaks = pd.concat([read_narrowpeak(path) for path in files], ignore_index=True)
    peaks = peaks[peaks["chrom"].str.match(r"^chr([0-9]+|X|Y)$", na=False)].copy()
    peaks["start"] = peaks["start"].astype(int)
    peaks["end"] = peaks["end"].astype(int)
    peaks["signalValue"] = pd.to_numeric(peaks["signalValue"], errors="coerce").fillna(0.0)
    return peaks


def promoter_intervals() -> pd.DataFrame:
    genes = pd.read_csv(GENCODE_PATH, sep="\t")
    genes = genes[genes["gene_type"].eq("protein_coding")].copy()
    genes["tss"] = genes.apply(lambda r: r["start"] if r["strand"] == "+" else r["end"], axis=1)
    genes["promoter_start"] = (genes["tss"] - PROMOTER_EXCLUSION).clip(lower=0).astype(int)
    genes["promoter_end"] = (genes["tss"] + PROMOTER_EXCLUSION).astype(int)
    return genes[["chromosome", "promoter_start", "promoter_end"]].rename(
        columns={"chromosome": "chrom"}
    )


def overlaps_any(region: pd.Series, intervals: pd.DataFrame) -> bool:
    sub = intervals[intervals["chrom"].eq(region["chrom"])]
    if sub.empty:
        return False
    return bool(((sub["promoter_start"] <= region["end"]) & (sub["promoter_end"] >= region["start"])).any())


def stitch_peaks(peaks: pd.DataFrame, promoters: pd.DataFrame) -> pd.DataFrame:
    stitched_rows = []
    cluster_id = 0
    for chrom, grp in peaks.sort_values(["chrom", "start", "end"]).groupby("chrom", sort=False):
        current = None
        for row in grp.itertuples(index=False):
            if current is None:
                cluster_id += 1
                current = {
                    "enhancer_id": f"MCF7_H3K27ac_stitched_{cluster_id}",
                    "chrom": chrom,
                    "start": int(row.start),
                    "end": int(row.end),
                    "n_peaks": 1,
                    "total_signal": float(row.signalValue),
                    "max_signal": float(row.signalValue),
                    "source_files": {row.source_file},
                }
                continue
            if int(row.start) - current["end"] <= STITCH_DISTANCE:
                current["end"] = max(current["end"], int(row.end))
                current["n_peaks"] += 1
                current["total_signal"] += float(row.signalValue)
                current["max_signal"] = max(current["max_signal"], float(row.signalValue))
                current["source_files"].add(row.source_file)
            else:
                stitched_rows.append(current)
                cluster_id += 1
                current = {
                    "enhancer_id": f"MCF7_H3K27ac_stitched_{cluster_id}",
                    "chrom": chrom,
                    "start": int(row.start),
                    "end": int(row.end),
                    "n_peaks": 1,
                    "total_signal": float(row.signalValue),
                    "max_signal": float(row.signalValue),
                    "source_files": {row.source_file},
                }
        if current is not None:
            stitched_rows.append(current)

    stitched = pd.DataFrame(stitched_rows)
    stitched["source_files"] = stitched["source_files"].apply(lambda x: ",".join(sorted(x)))
    stitched["length"] = stitched["end"] - stitched["start"]
    stitched["promoter_overlap"] = stitched.apply(lambda r: overlaps_any(r, promoters), axis=1)
    stitched = stitched[~stitched["promoter_overlap"]].copy()
    stitched["rank_signal_desc"] = stitched["total_signal"].rank(method="first", ascending=False).astype(int)
    return stitched.sort_values("rank_signal_desc")


def call_super_enhancers(stitched: pd.DataFrame) -> tuple[pd.DataFrame, float, int]:
    ranked = stitched.sort_values("total_signal", ascending=True).reset_index(drop=True).copy()
    if ranked.shape[0] < 3:
        ranked["is_ROSE_like_super_enhancer"] = False
        return ranked, float("inf"), 0
    x = pd.Series(range(ranked.shape[0]), dtype=float)
    x = x / x.max()
    y = ranked["total_signal"].astype(float)
    y = (y - y.min()) / (y.max() - y.min())
    idx = int((y - x).idxmax())
    cutoff = float(ranked.loc[idx, "total_signal"])
    ranked["ROSE_like_rank_low_to_high"] = ranked.index + 1
    ranked["is_ROSE_like_super_enhancer"] = ranked["total_signal"] >= cutoff
    if ranked["is_ROSE_like_super_enhancer"].mean() > 0.20:
        cutoff = float(ranked["total_signal"].quantile(0.95))
        idx = int((ranked["total_signal"] >= cutoff).idxmax())
        ranked["is_ROSE_like_super_enhancer"] = ranked["total_signal"] >= cutoff
    return ranked, cutoff, idx + 1


def link_se_to_candidate_tfs(super_enhancers: pd.DataFrame) -> pd.DataFrame:
    crc = pd.read_csv(CRC_TF_PATH, sep="\t")
    genes = pd.read_csv(GENCODE_PATH, sep="\t")
    tf_genes = genes[genes["gene_name"].isin(crc["TF"])].copy()
    tf_genes = tf_genes.sort_values(["gene_name", "gene_type"]).drop_duplicates("gene_name")
    rows = []
    for gene in tf_genes.itertuples(index=False):
        window_start = max(0, int(gene.start) - GENE_WINDOW)
        window_end = int(gene.end) + GENE_WINDOW
        sub = super_enhancers[
            super_enhancers["chrom"].eq(gene.chromosome)
            & (super_enhancers["start"] <= window_end)
            & (super_enhancers["end"] >= window_start)
        ]
        if sub.empty:
            rows.append(
                {
                    "TF": gene.gene_name,
                    "chromosome": gene.chromosome,
                    "gene_start": int(gene.start),
                    "gene_end": int(gene.end),
                    "n_ROSE_like_SE_100kb": 0,
                    "max_ROSE_like_SE_signal": 0.0,
                    "ROSE_like_SE_ids": "",
                    "MCF7_ROSE_like_SE_evidence": False,
                }
            )
        else:
            rows.append(
                {
                    "TF": gene.gene_name,
                    "chromosome": gene.chromosome,
                    "gene_start": int(gene.start),
                    "gene_end": int(gene.end),
                    "n_ROSE_like_SE_100kb": int(sub.shape[0]),
                    "max_ROSE_like_SE_signal": float(sub["total_signal"].max()),
                    "ROSE_like_SE_ids": ";".join(sub["enhancer_id"].astype(str).head(20)),
                    "MCF7_ROSE_like_SE_evidence": True,
                }
            )
    out = pd.DataFrame(rows)
    out = out.merge(crc[["TF", "CRC_evidence_score", "CRC_confidence"]], on="TF", how="left")
    return out.sort_values(["MCF7_ROSE_like_SE_evidence", "max_ROSE_like_SE_signal"], ascending=[False, False])


def main() -> None:
    peaks = load_all_peaks()
    promoters = promoter_intervals()
    stitched = stitch_peaks(peaks, promoters)
    ranked, cutoff, cutoff_rank = call_super_enhancers(stitched)
    super_enhancers = ranked[ranked["is_ROSE_like_super_enhancer"]].copy()
    tf_overlap = link_se_to_candidate_tfs(super_enhancers)
    crc = pd.read_csv(CRC_TF_PATH, sep="\t")
    integrated = crc.merge(
        tf_overlap[
            [
                "TF",
                "n_ROSE_like_SE_100kb",
                "max_ROSE_like_SE_signal",
                "ROSE_like_SE_ids",
                "MCF7_ROSE_like_SE_evidence",
            ]
        ],
        on="TF",
        how="left",
    )
    integrated["n_ROSE_like_SE_100kb"] = integrated["n_ROSE_like_SE_100kb"].fillna(0).astype(int)
    integrated["max_ROSE_like_SE_signal"] = integrated["max_ROSE_like_SE_signal"].fillna(0.0)
    integrated["ROSE_like_SE_ids"] = integrated["ROSE_like_SE_ids"].fillna("")
    integrated["MCF7_ROSE_like_SE_evidence"] = integrated["MCF7_ROSE_like_SE_evidence"].fillna(False)
    integrated["CRC_evidence_score_with_MCF7_ROSE_like_SE_proxy"] = (
        integrated["CRC_evidence_score"].astype(int)
        + integrated["MCF7_ROSE_like_SE_evidence"].astype(int)
    )

    stitched.to_csv(RESULTS_DIR / "ENCODE_MCF7_ROSE_like_stitched_enhancers.tsv", sep="\t", index=False)
    super_enhancers.to_csv(RESULTS_DIR / "ENCODE_MCF7_ROSE_like_super_enhancers.tsv", sep="\t", index=False)
    tf_overlap.to_csv(RESULTS_DIR / "ENCODE_MCF7_ROSE_like_SE_candidate_TF_overlap.tsv", sep="\t", index=False)
    integrated.to_csv(
        TCGA_RESULTS_DIR / "CRC_member_evidence_score_table_with_MCF7_ROSE_like_SE_proxy.tsv",
        sep="\t",
        index=False,
    )

    summary = {
        "n_input_H3K27ac_peaks": int(peaks.shape[0]),
        "n_stitched_non_promoter_enhancers": int(stitched.shape[0]),
        "stitch_distance_bp": STITCH_DISTANCE,
        "promoter_exclusion_bp": PROMOTER_EXCLUSION,
        "ROSE_like_signal_cutoff": cutoff,
        "ROSE_like_cutoff_rank_low_to_high": cutoff_rank,
        "n_ROSE_like_super_enhancers": int(super_enhancers.shape[0]),
        "n_candidate_TFs_with_ROSE_like_SE_100kb": int(tf_overlap["MCF7_ROSE_like_SE_evidence"].sum()),
        "note": "ROSE-like stitching/ranking from ENCODE MCF-7 H3K27ac narrowPeak files; if tangent cutoff is too permissive, top 5% signal is used as a conservative fallback. Not an official ROSE/LILY run.",
    }
    pd.DataFrame([summary]).to_csv(
        RESULTS_DIR / "ENCODE_MCF7_ROSE_like_SE_summary.tsv", sep="\t", index=False
    )
    (RESULTS_DIR / "ENCODE_MCF7_ROSE_like_SE_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
