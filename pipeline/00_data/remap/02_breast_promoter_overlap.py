from __future__ import annotations

import math
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


ROOT = (CORE_ROOT / "remap_reanalysis")
CORE_DATA = ROOT.parent
MANUSCRIPT = CORE_DATA.parent
WORKSPACE = MANUSCRIPT.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
PEAK_DIR = DATA / "breast_tier1_dataset_peaks"

GENCODE = DATA_ROOT / "GENCODE" / "results" / "gencode_v36_gene_annotation.tsv"
MANIFEST = REPO_ROOT / "data" / "dataset_lists" / "remap_breast_core16_tier1_download_manifest.tsv"

CORE16 = [
    "ESR1",
    "FOXA1",
    "SPDEF",
    "GATA3",
    "PBX1",
    "AR",
    "CREB3L4",
    "TFAP2C",
    "EHF",
    "MYBL2",
    "NR2F6",
    "STAT1",
    "ELF3",
    "RUNX1",
    "TP63",
    "FOXC1",
]

MODULE = {
    "ESR1": "luminal",
    "FOXA1": "luminal",
    "SPDEF": "luminal",
    "GATA3": "luminal",
    "PBX1": "luminal",
    "AR": "luminal",
    "CREB3L4": "luminal",
    "TFAP2C": "her2",
    "EHF": "her2",
    "MYBL2": "her2",
    "NR2F6": "her2",
    "STAT1": "her2",
    "ELF3": "her2",
    "RUNX1": "basal",
    "TP63": "basal",
    "FOXC1": "basal",
}

PALETTE = {
    "red": "#c35655",
    "blue": "#5e90b8",
    "light_green": "#b1c5ba",
    "light_orange": "#f5c19b",
    "dark_green": "#6d817d",
    "pink": "#eea5d1",
    "dark_gray": "#3e3c3d",
    "light_gray": "#e8e8e8",
    "white": "#ffffff",
}
MODULE_COLORS = {"luminal": PALETTE["blue"], "her2": PALETTE["pink"], "basal": PALETTE["red"]}


def bool_num(x) -> float:
    if pd.isna(x):
        return 0.0
    if isinstance(x, str):
        return 1.0 if x.strip().lower() in {"true", "1", "yes", "y"} else 0.0
    return 1.0 if bool(x) else 0.0


def load_promoters() -> pd.DataFrame:
    genes = pd.read_csv(GENCODE, sep="\t")
    genes = genes[genes["gene_name"].isin(CORE16)].copy()
    genes["tss"] = np.where(genes["strand"].eq("+"), genes["start"], genes["end"]).astype(int)
    genes["promoter_start"] = (genes["tss"] - 2000).clip(lower=0).astype(int)
    genes["promoter_end"] = (genes["tss"] + 2000).astype(int)
    genes["module"] = genes["gene_name"].map(MODULE)
    return genes[["gene_name", "module", "chromosome", "strand", "start", "end", "tss", "promoter_start", "promoter_end"]]


def read_bed(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2],
        names=["chromosome", "start", "end"],
        dtype={"chromosome": str, "start": np.int64, "end": np.int64},
    )


def overlap_summary(peaks: pd.DataFrame, promoters: pd.DataFrame, regulator: str) -> pd.DataFrame:
    if peaks.empty:
        return pd.DataFrame()
    peaks = peaks.dropna(subset=["chromosome", "start", "end"]).copy()
    peaks = peaks[peaks["end"].gt(peaks["start"])]
    peaks["tf"] = regulator

    if pr is not None:
        peak_pr = pr.PyRanges(
            peaks.rename(columns={"chromosome": "Chromosome", "start": "Start", "end": "End"})[
                ["Chromosome", "Start", "End", "tf"]
            ]
        )
        prom_pr = pr.PyRanges(
            promoters.rename(
                columns={"chromosome": "Chromosome", "promoter_start": "Start", "promoter_end": "End"}
            )[["Chromosome", "Start", "End", "gene_name", "module", "tss"]]
        )
        joined = peak_pr.join(prom_pr).df
        if joined.empty:
            return pd.DataFrame()
        joined["overlap_bp"] = np.minimum(joined["End"], joined["End_b"]) - np.maximum(joined["Start"], joined["Start_b"])
        return (
            joined.groupby(["tf", "gene_name", "module"], as_index=False)
            .agg(n_peak_overlaps=("overlap_bp", "size"), max_overlap_bp=("overlap_bp", "max"), median_overlap_bp=("overlap_bp", "median"))
        )

    rows: list[dict[str, object]] = []
    for chrom, peak_chr in peaks.groupby("chromosome", sort=False):
        prom_chr = promoters[promoters["chromosome"].eq(chrom)]
        if prom_chr.empty:
            continue
        peak_chr = peak_chr.sort_values("start").reset_index(drop=True)
        starts = peak_chr["start"].to_numpy(np.int64)
        ends = peak_chr["end"].to_numpy(np.int64)
        for prom in prom_chr.itertuples(index=False):
            p_start = int(prom.promoter_start)
            p_end = int(prom.promoter_end)
            limit = np.searchsorted(starts, p_end, side="left")
            if limit == 0:
                continue
            idx = np.nonzero(ends[:limit] > p_start)[0]
            if idx.size == 0:
                continue
            ov = np.minimum(ends[idx], p_end) - np.maximum(starts[idx], p_start)
            rows.append(
                {
                    "tf": regulator,
                    "gene_name": prom.gene_name,
                    "module": prom.module,
                    "n_peak_overlaps": int(idx.size),
                    "max_overlap_bp": int(ov.max()),
                    "median_overlap_bp": float(np.median(ov)),
                }
            )
    return pd.DataFrame(rows)


def compute_breast_remap_overlaps() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    promoters = load_promoters()
    manifest = pd.read_csv(MANIFEST, sep="\t")
    manifest = manifest.drop_duplicates("dataset_file").copy()
    manifest["local_peak_file"] = manifest["dataset_file"].map(lambda x: str(PEAK_DIR / str(x)))
    manifest["local_peak_file_present"] = manifest["dataset_file"].map(lambda x: (PEAK_DIR / str(x)).exists())
    manifest.to_csv(RESULTS / "remap_breast_tier1_download_status.tsv", sep="\t", index=False)
    missing = manifest[~manifest["local_peak_file_present"]].copy()
    missing.to_csv(RESULTS / "remap_breast_tier1_missing_downloads.tsv", sep="\t", index=False)

    all_rows: list[pd.DataFrame] = []
    errors: list[dict[str, str]] = []
    for row in manifest[manifest["local_peak_file_present"]].itertuples(index=False):
        path = PEAK_DIR / str(row.dataset_file)
        try:
            peaks = read_bed(path)
            one = overlap_summary(peaks, promoters, str(row.target_name))
        except Exception as exc:
            errors.append({"dataset_file": str(row.dataset_file), "target_name": str(row.target_name), "error": str(exc)})
            continue
        if one.empty:
            continue
        one["dataset_file"] = str(row.dataset_file)
        one["biotype"] = str(row.biotype)
        one["breast_context_group"] = str(row.breast_context_group)
        one["experiment"] = str(row.experiment)
        one["source"] = str(row.source)
        one["download_url"] = str(row.download_url)
        one["n_peaks_reported"] = row.n_peaks_reported
        all_rows.append(one)

    if errors:
        pd.DataFrame(errors).to_csv(RESULTS / "remap_breast_tier1_overlap_read_errors.tsv", sep="\t", index=False)

    if all_rows:
        dataset_overlap = pd.concat(all_rows, ignore_index=True)
    else:
        dataset_overlap = pd.DataFrame(
            columns=[
                "tf",
                "gene_name",
                "module",
                "n_peak_overlaps",
                "max_overlap_bp",
                "median_overlap_bp",
                "dataset_file",
                "biotype",
                "breast_context_group",
                "experiment",
                "source",
                "download_url",
                "n_peaks_reported",
            ]
        )
    dataset_overlap = dataset_overlap.rename(columns={"tf": "regulator_TF", "gene_name": "target_TF", "module": "target_module"})
    dataset_overlap.to_csv(RESULTS / "remap_breast_tier1_promoter_overlap_by_dataset.tsv", sep="\t", index=False)

    edge = (
        dataset_overlap.groupby(["regulator_TF", "target_TF", "target_module"], as_index=False)
        .agg(
            n_breast_remap_dataset_overlaps=("dataset_file", "nunique"),
            n_breast_remap_peak_overlaps=("n_peak_overlaps", "sum"),
            max_breast_remap_peak_overlaps_per_dataset=("n_peak_overlaps", "max"),
            max_breast_remap_overlap_bp=("max_overlap_bp", "max"),
            n_breast_biotypes=("biotype", "nunique"),
            breast_biotypes=("biotype", lambda x: ";".join(sorted(set(map(str, x))))),
            breast_experiments=("experiment", lambda x: ";".join(sorted(set(map(str, x))))),
            breast_dataset_files=("dataset_file", lambda x: ";".join(sorted(set(map(str, x))))),
        )
    )
    edge["breast_remap_promoter_binding_evidence"] = edge["n_breast_remap_dataset_overlaps"].gt(0)
    edge["breast_remap_support_level"] = np.select(
        [edge["n_breast_biotypes"].ge(2), edge["n_breast_remap_dataset_overlaps"].ge(2)],
        ["strong_multi_biotype", "moderate_multi_dataset"],
        default="single_dataset",
    )
    edge.to_csv(RESULTS / "remap_breast_tier1_promoter_overlap_edges.tsv", sep="\t", index=False)

    tf_support = (
        edge.groupby("regulator_TF", as_index=False)
        .agg(
            n_target_promoters_with_breast_remap=("target_TF", "nunique"),
            n_breast_remap_dataset_overlaps=("n_breast_remap_dataset_overlaps", "sum"),
            n_breast_remap_peak_overlaps=("n_breast_remap_peak_overlaps", "sum"),
            n_breast_biotypes=("n_breast_biotypes", "max"),
        )
    )
    tf_support["breast_remap_TF_edge_evidence"] = tf_support["n_target_promoters_with_breast_remap"].gt(0)
    tf_support.to_csv(RESULTS / "remap_breast_tier1_TF_support_summary.tsv", sep="\t", index=False)
    return dataset_overlap, edge, tf_support


def main() -> None:
    dataset_overlap, edge, tf_support = compute_breast_remap_overlaps()
    summary = {
        "n_dataset_promoter_overlap_rows": int(dataset_overlap.shape[0]),
        "n_breast_remap_supported_edges": int(edge.shape[0]),
        "n_regulator_TFs_with_breast_remap_support": int(tf_support.shape[0]),
        "interval_engine": "PyRanges" if pr is not None else "chromosome-grouped numpy.searchsorted fallback",
    }
    pd.DataFrame([summary]).to_csv(RESULTS / "remap_breast_redraw_summary.tsv", sep="\t", index=False)
    print(summary)


if __name__ == "__main__":
    main()
