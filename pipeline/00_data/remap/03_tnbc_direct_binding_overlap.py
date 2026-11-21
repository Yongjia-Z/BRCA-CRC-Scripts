from __future__ import annotations

import gzip
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
WORKSPACE = DATA_ROOT.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
TNBC_RESULTS = RESULTS / "TNBC_direct_binding"
PEAK_DIR = DATA / "TNBC_direct_binding_peaks"


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


def find_gencode() -> Path:
    hits = list(WORKSPACE.glob("*/GENCODE/results/gencode_v36_gene_annotation.tsv"))
    if not hits:
        raise FileNotFoundError("Could not locate gencode_v36_gene_annotation.tsv under the workspace.")
    return hits[0]


def load_gene_windows(flank: int) -> pd.DataFrame:
    genes = pd.read_csv(find_gencode(), sep="\t")
    genes = genes[genes["gene_name"].isin(CORE16)].copy()
    genes["tss"] = np.where(genes["strand"].eq("+"), genes["start"], genes["end"]).astype(np.int64)
    if flank == 2000:
        genes["window_start"] = (genes["tss"] - flank).clip(lower=0).astype(np.int64)
        genes["window_end"] = (genes["tss"] + flank).astype(np.int64)
        window_type = "TSS_pm_2kb_promoter"
    else:
        genes["window_start"] = (genes["start"] - flank).clip(lower=0).astype(np.int64)
        genes["window_end"] = (genes["end"] + flank).astype(np.int64)
        window_type = f"gene_body_pm_{flank // 1000}kb_locus"
    genes["module"] = genes["gene_name"].map(MODULE)
    genes["window_type"] = window_type
    return genes[
        [
            "gene_name",
            "module",
            "chromosome",
            "strand",
            "start",
            "end",
            "tss",
            "window_start",
            "window_end",
            "window_type",
        ]
    ]


def build_manifest() -> pd.DataFrame:
    rows = [
        {
            "dataset": "GSE213840",
            "regulator_TF": "FOXC1",
            "cell_line": "BT549",
            "subtype_context": "TNBC/basal-like",
            "assay": "ChIP-seq",
            "source": "GEO-derived processed narrowPeak",
            "assembly_status": "used_as_hg38_in_current_project",
            "usable_for_hg38_overlap": True,
            "file_name": "GSE213840_FOXC1_BT549_hg38.narrowPeak.gz",
        },
        {
            "dataset": "GSE213840",
            "regulator_TF": "FOXC1",
            "cell_line": "Hs578T",
            "subtype_context": "TNBC/basal-like",
            "assay": "ChIP-seq",
            "source": "GEO-derived processed narrowPeak",
            "assembly_status": "used_as_hg38_in_current_project",
            "usable_for_hg38_overlap": True,
            "file_name": "GSE213840_FOXC1_Hs578T_hg38.narrowPeak.gz",
        },
        {
            "dataset": "GSE213840",
            "regulator_TF": "FOXC1",
            "cell_line": "MDA-MB-231",
            "subtype_context": "TNBC/basal-like",
            "assay": "ChIP-seq",
            "source": "GEO-derived processed narrowPeak",
            "assembly_status": "used_as_hg38_in_current_project",
            "usable_for_hg38_overlap": True,
            "file_name": "GSE213840_FOXC1_MDA-MB-231_hg38.narrowPeak.gz",
        },
        {
            "dataset": "GSE213840",
            "regulator_TF": "FOXC1",
            "cell_line": "MDA-MB-468",
            "subtype_context": "TNBC/basal-like",
            "assay": "ChIP-seq",
            "source": "GEO-derived processed narrowPeak",
            "assembly_status": "used_as_hg38_in_current_project",
            "usable_for_hg38_overlap": True,
            "file_name": "GSE213840_FOXC1_MDA-MB-468_hg38.narrowPeak.gz",
        },
        {
            "dataset": "GSE72009",
            "regulator_TF": "TP63",
            "cell_line": "MDA-MB-231",
            "subtype_context": "TNBC/basal-like",
            "assay": "ChIP-seq",
            "source": "GEO original supplementary narrowPeak",
            "assembly_status": "original_GEO_coordinates_kept_for_traceability_not_used_in_hg38_overlap",
            "usable_for_hg38_overlap": False,
            "file_name": "GSE72009_TP63_MDA-MB-231_original_narrowPeak.gz",
        },
        {
            "dataset": "GSE72009",
            "regulator_TF": "TP63",
            "cell_line": "MDA-MB-231",
            "subtype_context": "TNBC/basal-like",
            "assay": "ChIP-seq",
            "source": "GEO GSE72009 reprocessed by ReMap2022 dataset-level hg38 peak",
            "assembly_status": "hg38_ReMap2022_reprocessed_dataset_level",
            "usable_for_hg38_overlap": True,
            "file_name": "GSE72009_TP63_MDA-MB-231_ReMap2022_hg38.bed.gz",
        },
    ]
    manifest = pd.DataFrame(rows)
    manifest["local_peak_file"] = manifest["file_name"].map(lambda x: str(PEAK_DIR / x))
    manifest["file_exists"] = manifest["file_name"].map(lambda x: (PEAK_DIR / x).exists())
    peak_counts = []
    for file_name in manifest["file_name"]:
        path = PEAK_DIR / file_name
        if not path.exists():
            peak_counts.append(np.nan)
            continue
        with gzip.open(path, "rt") as handle:
            peak_counts.append(sum(1 for _ in handle if _.strip() and not _.startswith("#")))
    manifest["n_peaks_in_file"] = peak_counts
    return manifest


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


def interval_overlap(peaks: pd.DataFrame, windows: pd.DataFrame, regulator_tf: str) -> pd.DataFrame:
    peaks = peaks.dropna(subset=["chromosome", "start", "end"]).copy()
    peaks = peaks[peaks["end"].gt(peaks["start"])]
    if peaks.empty:
        return pd.DataFrame()

    if pr is not None:
        peak_pr = pr.PyRanges(
            peaks.rename(columns={"chromosome": "Chromosome", "start": "Start", "end": "End"})[
                ["Chromosome", "Start", "End"]
            ]
        )
        window_pr = pr.PyRanges(
            windows.rename(columns={"chromosome": "Chromosome", "window_start": "Start", "window_end": "End"})[
                ["Chromosome", "Start", "End", "gene_name", "module", "window_type"]
            ]
        )
        joined = peak_pr.join(window_pr).df
        if joined.empty:
            return pd.DataFrame()
        joined["overlap_bp"] = np.minimum(joined["End"], joined["End_b"]) - np.maximum(joined["Start"], joined["Start_b"])
        out = (
            joined.groupby(["gene_name", "module", "window_type"], as_index=False)
            .agg(
                n_peak_overlaps=("overlap_bp", "size"),
                max_overlap_bp=("overlap_bp", "max"),
                median_overlap_bp=("overlap_bp", "median"),
            )
        )
        out["regulator_TF"] = regulator_tf
        return out

    rows: list[dict[str, object]] = []
    for chrom, peak_chr in peaks.groupby("chromosome", sort=False):
        window_chr = windows[windows["chromosome"].eq(chrom)]
        if window_chr.empty:
            continue
        peak_chr = peak_chr.sort_values("start").reset_index(drop=True)
        starts = peak_chr["start"].to_numpy(np.int64)
        ends = peak_chr["end"].to_numpy(np.int64)
        for window in window_chr.itertuples(index=False):
            start = int(window.window_start)
            end = int(window.window_end)
            limit = np.searchsorted(starts, end, side="left")
            if limit == 0:
                continue
            idx = np.nonzero(ends[:limit] > start)[0]
            if idx.size == 0:
                continue
            overlap_bp = np.minimum(ends[idx], end) - np.maximum(starts[idx], start)
            rows.append(
                {
                    "gene_name": window.gene_name,
                    "module": window.module,
                    "window_type": window.window_type,
                    "n_peak_overlaps": int(idx.size),
                    "max_overlap_bp": int(overlap_bp.max()),
                    "median_overlap_bp": float(np.median(overlap_bp)),
                    "regulator_TF": regulator_tf,
                }
            )
    return pd.DataFrame(rows)


def compute_tnbc_overlaps(manifest: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    promoter_windows = load_gene_windows(2000)
    locus_windows = load_gene_windows(100000)
    windows = pd.concat([promoter_windows, locus_windows], ignore_index=True)

    all_rows: list[pd.DataFrame] = []
    usable = manifest[manifest["file_exists"] & manifest["usable_for_hg38_overlap"]].copy()
    for rec in usable.itertuples(index=False):
        peaks = read_bed(Path(rec.local_peak_file))
        one = interval_overlap(peaks, windows, rec.regulator_TF)
        if one.empty:
            continue
        one["dataset"] = rec.dataset
        one["cell_line"] = rec.cell_line
        one["subtype_context"] = rec.subtype_context
        one["assay"] = rec.assay
        one["source"] = rec.source
        one["assembly_status"] = rec.assembly_status
        one["file_name"] = rec.file_name
        one["n_peaks_in_file"] = rec.n_peaks_in_file
        all_rows.append(one)

    if all_rows:
        by_dataset = pd.concat(all_rows, ignore_index=True)
    else:
        by_dataset = pd.DataFrame(
            columns=[
                "gene_name",
                "module",
                "window_type",
                "n_peak_overlaps",
                "max_overlap_bp",
                "median_overlap_bp",
                "regulator_TF",
                "dataset",
                "cell_line",
                "subtype_context",
                "assay",
                "source",
                "assembly_status",
                "file_name",
                "n_peaks_in_file",
            ]
        )

    by_dataset = by_dataset.rename(columns={"gene_name": "target_TF", "module": "target_module"})
    edge = (
        by_dataset.groupby(["regulator_TF", "target_TF", "target_module", "window_type"], as_index=False)
        .agg(
            n_tnbc_datasets=("file_name", "nunique"),
            n_tnbc_cell_lines=("cell_line", "nunique"),
            n_tnbc_peak_overlaps=("n_peak_overlaps", "sum"),
            max_tnbc_peak_overlaps_per_dataset=("n_peak_overlaps", "max"),
            max_tnbc_overlap_bp=("max_overlap_bp", "max"),
            tnbc_cell_lines=("cell_line", lambda x: ";".join(sorted(set(map(str, x))))),
            tnbc_datasets=("dataset", lambda x: ";".join(sorted(set(map(str, x))))),
            tnbc_files=("file_name", lambda x: ";".join(sorted(set(map(str, x))))),
        )
    )
    edge["TNBC_direct_binding_evidence"] = edge["n_tnbc_datasets"].gt(0)
    return by_dataset, edge


def wide_edge(edge: pd.DataFrame, window_type: str, value_col: str = "n_tnbc_peak_overlaps") -> pd.DataFrame:
    subset = edge[edge["window_type"].eq(window_type)].copy()
    if subset.empty:
        return pd.DataFrame(0, index=[], columns=CORE16)
    mat = subset.pivot_table(index="regulator_TF", columns="target_TF", values=value_col, aggfunc="sum", fill_value=0)
    mat = mat.reindex(index=["FOXC1", "TP63"], columns=CORE16, fill_value=0)
    return mat


def main() -> None:
    manifest = build_manifest()
    manifest.to_csv(TNBC_RESULTS / "tnbc_direct_binding_dataset_manifest.tsv", sep="\t", index=False)

    by_dataset, edge = compute_tnbc_overlaps(manifest)
    by_dataset.to_csv(TNBC_RESULTS / "tnbc_direct_binding_overlap_by_dataset.tsv", sep="\t", index=False)
    edge.to_csv(TNBC_RESULTS / "tnbc_direct_binding_overlap_edges.tsv", sep="\t", index=False)

    promoter_mat = wide_edge(edge, "TSS_pm_2kb_promoter")
    locus_mat = wide_edge(edge, "gene_body_pm_100kb_locus")
    promoter_mat.to_csv(TNBC_RESULTS / "figure4c_tnbc_direct_promoter_overlap_matrix.tsv", sep="\t")
    locus_mat.to_csv(TNBC_RESULTS / "figure4c_tnbc_direct_locus_100kb_overlap_matrix.tsv", sep="\t")


    summary = pd.DataFrame(
        [
            {
                "n_manifest_rows": len(manifest),
                "n_hg38_overlap_usable_rows": int(manifest["usable_for_hg38_overlap"].sum()),
                "n_dataset_overlap_rows": len(by_dataset),
                "n_edge_window_rows": len(edge),
                "n_promoter_edges": int(edge["window_type"].eq("TSS_pm_2kb_promoter").sum()),
                "n_locus_edges": int(edge["window_type"].str.contains("gene_body_pm_100kb", regex=False).sum()),
                "interval_engine": "pyranges" if pr is not None else "chromosome-grouped numpy.searchsorted fallback",
            }
        ]
    )
    summary.to_csv(TNBC_RESULTS / "tnbc_direct_binding_overlap_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
