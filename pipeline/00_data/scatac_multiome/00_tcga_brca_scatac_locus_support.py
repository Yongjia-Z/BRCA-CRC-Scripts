from __future__ import annotations

import json
import shutil
import tarfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))

ROOT = DATA_ROOT.parent
BASE = DATA_ROOT / "scATAC_multiome"
DATA = BASE / "data"
RESULTS = BASE / "results"
FIG8 = ROOT / "data" / "figure8"

RESULTS.mkdir(parents=True, exist_ok=True)
FIG8.mkdir(parents=True, exist_ok=True)


MODULES = {
    "luminal_CRC": ["FOXA1", "GATA3", "ESR1", "SPDEF", "PBX1", "AR", "CREB3L4"],
    "basal_CRC": ["RUNX1", "TP63", "FOXC1"],
    "her2_CRC": ["TFAP2C", "MYBL2", "EHF", "STAT1", "ELF3", "NR2F6"],
}
CORE_TFS = sorted({g for genes in MODULES.values() for g in genes})
MOTIF_ALIASES = {
    "FOXA1": ["FOX/"],
    "FOXC1": ["FOX/"],
    "GATA3": ["GATA"],
    "ESR1": ["NR/"],
    "AR": ["NR/"],
    "NR2F6": ["NR/"],
    "SPDEF": ["SPDEF/"],
    "PBX1": ["HD/"],
    "CREB3L4": ["CREB3/XBP1", "CREB/ATF/"],
    "RUNX1": ["RUNX/"],
    "TP63": ["P53-like/"],
    "TFAP2C": ["TFAP2/"],
    "MYBL2": ["MYB/"],
    "EHF": ["ETS/"],
    "ELF3": ["ETS/"],
    "STAT1": ["STAT/"],
}


def copy_fig(path: Path) -> None:
    shutil.copy2(path, FIG8 / path.name)


def read_tar_csv(tar_path: Path, member_suffix: str, **kwargs) -> pd.DataFrame:
    with tarfile.open(tar_path, "r:gz") as tar:
        member = next(m for m in tar.getmembers() if m.isfile() and m.name.endswith(member_suffix))
        with tar.extractfile(member) as fh:
            return pd.read_csv(fh, **kwargs)


def analyze_peaksets() -> tuple[pd.DataFrame, pd.DataFrame]:
    tar_path = DATA / "BRCA_peaksets.tar.gz"
    peak_members = {
        "BRCA_cancer_control": "BRCA_cancer_control_peakset.csv",
        "Breast_BASAL_healthy": "Breast_BASAL_healthy.csv",
        "Breast_LuminalHR_healthy": "Breast_LuminalHR_healthy.csv",
        "Breast_LuminalSEC_healthy": "Breast_LuminalSEC_healthy.csv",
    }
    summaries = []
    tf_rows = []
    usecols = ["seqnames", "start", "end", "score", "replicateScoreQuantile", "groupScoreQuantile", "Reproducibility", "distToGeneStart", "nearestGene", "peakType", "distToTSS", "GC"]
    for label, suffix in peak_members.items():
        df = read_tar_csv(tar_path, suffix, usecols=lambda c: c in usecols)
        df["peakset"] = label
        df["nearestGene"] = df["nearestGene"].astype(str)
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
        df["groupScoreQuantile"] = pd.to_numeric(df["groupScoreQuantile"], errors="coerce")
        tf = df[df["nearestGene"].isin(CORE_TFS)].copy()
        tf_rows.append(tf)
        by_gene = (
            tf.groupby(["peakset", "nearestGene", "peakType"], dropna=False)
            .agg(
                n_peaks=("nearestGene", "count"),
                mean_score=("score", "mean"),
                median_score=("score", "median"),
                mean_group_score_quantile=("groupScoreQuantile", "mean"),
                min_abs_dist_to_tss=("distToTSS", lambda x: float(np.nanmin(np.abs(pd.to_numeric(x, errors="coerce")))) if len(x) else np.nan),
            )
            .reset_index()
        )
        summaries.append(by_gene)
        del df
    summary = pd.concat(summaries, ignore_index=True)
    tf_peaks = pd.concat(tf_rows, ignore_index=True)
    summary.to_csv(RESULTS / "TCGA_BRCA_scATAC_CRC_TF_locus_peak_summary.tsv", sep="\t", index=False)
    tf_peaks.to_csv(RESULTS / "TCGA_BRCA_scATAC_CRC_TF_locus_peaks.tsv.gz", sep="\t", index=False, compression="gzip")
    return summary, tf_peaks


def analyze_motifs() -> tuple[pd.DataFrame, pd.DataFrame]:
    tar_path = DATA / "Breast_cleaned_Cancer_peakset.tar.gz"
    rows = []
    total_rows = []
    with tarfile.open(tar_path, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(".csv")]
        for member in members:
            sample = Path(member.name).name.replace("_cleaned_motifs.csv", "")
            total = 0
            counter = Counter()
            with tar.extractfile(member) as fh:
                for chunk in pd.read_csv(fh, usecols=["group_name"], chunksize=1_000_000):
                    names = chunk["group_name"].astype(str)
                    total += len(names)
                    counter.update(names)
            total_rows.append({"sample": sample, "total_motif_instances": total})
            for module, genes in MODULES.items():
                for gene in genes:
                    aliases = MOTIF_ALIASES.get(gene, [gene])
                    matched = [name for name in counter if any(name == alias or name.startswith(alias) for alias in aliases)]
                    rows.append(
                        {
                            "sample": sample,
                            "module": module,
                            "TF_motif": gene,
                            "matched_motif_groups": ";".join(sorted(matched)),
                            "motif_instances": int(sum(counter[name] for name in matched)),
                        }
                    )
    counts = pd.DataFrame(rows).merge(pd.DataFrame(total_rows), on="sample", how="left")
    counts["motif_instances_per_million"] = counts["motif_instances"] / counts["total_motif_instances"] * 1e6
    counts.to_csv(RESULTS / "TCGA_BRCA_scATAC_CRC_TF_motif_counts_by_sample.tsv", sep="\t", index=False)
    module = (
        counts.groupby(["sample", "module"])
        .agg(
            total_crc_motif_instances=("motif_instances", "sum"),
            mean_motif_instances_per_million=("motif_instances_per_million", "mean"),
            total_motif_instances=("total_motif_instances", "first"),
        )
        .reset_index()
    )
    module.to_csv(RESULTS / "TCGA_BRCA_scATAC_CRC_module_motif_counts_by_sample.tsv", sep="\t", index=False)
    return counts, module


def stats_tables(peak_summary: pd.DataFrame, motif_counts: pd.DataFrame, module_counts: pd.DataFrame) -> dict:
    motif_stats = (
        motif_counts.groupby(["module", "TF_motif"])
        .agg(
            n_samples=("sample", "nunique"),
            total_motif_instances=("motif_instances", "sum"),
            median_motif_instances_per_million=("motif_instances_per_million", "median"),
            mean_motif_instances_per_million=("motif_instances_per_million", "mean"),
        )
        .reset_index()
        .sort_values("total_motif_instances", ascending=False)
    )
    motif_stats.to_csv(RESULTS / "TCGA_BRCA_scATAC_CRC_TF_motif_summary.tsv", sep="\t", index=False)
    module_stats = (
        module_counts.groupby("module")
        .agg(
            n_samples=("sample", "nunique"),
            median_module_motif_per_million=("mean_motif_instances_per_million", "median"),
            mean_module_motif_per_million=("mean_motif_instances_per_million", "mean"),
        )
        .reset_index()
    )
    module_stats.to_csv(RESULTS / "TCGA_BRCA_scATAC_CRC_module_motif_summary.tsv", sep="\t", index=False)
    locus = peak_summary[peak_summary["peakset"].eq("BRCA_cancer_control")].groupby("nearestGene")["n_peaks"].sum().sort_values(ascending=False)
    return {
        "crc_tfs_with_cancer_locus_peaks": int((locus > 0).sum()),
        "top_locus_peak_genes": locus.head(8).to_dict(),
        "motif_samples": int(motif_counts["sample"].nunique()),
        "top_motif_tfs": motif_stats.head(8).to_dict(orient="records"),
    }


def main() -> None:
    peak_summary, tf_peaks = analyze_peaksets()
    motif_counts, module_counts = analyze_motifs()
    summary = stats_tables(peak_summary, motif_counts, module_counts)
    summary.update({"tf_locus_peak_rows": int(tf_peaks.shape[0]), "motif_count_rows": int(motif_counts.shape[0])})
    with open(RESULTS / "TCGA_BRCA_scATAC_CRC_open_chromatin_validation_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
