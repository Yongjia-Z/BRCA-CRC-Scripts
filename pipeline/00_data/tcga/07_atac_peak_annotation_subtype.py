from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TCGA")
DATA = ROOT / "data"
ATAC = DATA / "ATAC"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def read_brca_peak_annotation() -> pd.DataFrame:
    s2 = ATAC / "TCGA_ATAC_Data_S2_peak_calls.xlsx"
    peak = pd.read_excel(s2, sheet_name="BRCA-Only_PeakCalls", header=21)
    peak = peak.dropna(subset=["Peak_Name"])
    peak.to_csv(RESULTS / "ATAC_BRCA_peak_annotation.tsv", sep="\t", index=False)
    return peak


def main() -> None:
    peak_anno = read_brca_peak_annotation()
    matrix = pd.read_csv(ATAC / "brca_peak_Log2Counts_dedup.tsv", sep="\t").rename(columns={"sample": "Peak_Name"})
    sample_ann = pd.read_csv(RESULTS / "ATAC_sample_accessibility_summary.tsv", sep="\t")
    matrix = matrix[matrix["Peak_Name"].isin(set(peak_anno["Peak_Name"]))]

    merged_peak = peak_anno.merge(pd.read_csv(RESULTS / "ATAC_peak_accessibility_summary.tsv", sep="\t"), left_on="Peak_Name", right_on="peak_id", how="left")
    merged_peak.to_csv(RESULTS / "ATAC_BRCA_peak_annotation_with_accessibility.tsv", sep="\t", index=False)

    annotation_summary = peak_anno["Annotation"].fillna("missing").value_counts().rename_axis("Annotation").reset_index(name="n_peaks")
    annotation_summary.to_csv(RESULTS / "ATAC_BRCA_peak_genomic_annotation_summary.tsv", sep="\t", index=False)


    if "Distance To Nearest TSS" in peak_anno.columns:
        dist = pd.to_numeric(peak_anno["Distance To Nearest TSS"], errors="coerce").dropna()

    subtype_groups = {}
    for subtype in ["BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal"]:
        cols = sample_ann.loc[sample_ann.get("SUBTYPE") == subtype, "sample_id"].tolist()
        cols = [c for c in cols if c in matrix.columns]
        if len(cols) >= 5:
            subtype_groups[subtype] = cols

    subtype_result = matrix[["Peak_Name"]].copy()
    arrays = []
    for subtype, cols in subtype_groups.items():
        vals = matrix[cols].astype(float)
        subtype_result[f"mean_{subtype}"] = vals.mean(axis=1)
        arrays.append(vals.to_numpy())

    if len(arrays) >= 2:
        f_stat, pvals = stats.f_oneway(*arrays, axis=1)
        subtype_result["subtype_ANOVA_F"] = f_stat
        subtype_result["subtype_ANOVA_p"] = pvals
        subtype_result["subtype_ANOVA_FDR"] = multipletests(np.nan_to_num(pvals, nan=1.0), method="fdr_bh")[1]
        mean_cols = [c for c in subtype_result.columns if c.startswith("mean_BRCA_")]
        subtype_result["best_subtype"] = subtype_result[mean_cols].idxmax(axis=1).str.replace("mean_", "", regex=False)
        subtype_result["max_minus_min_subtype_mean"] = subtype_result[mean_cols].max(axis=1) - subtype_result[mean_cols].min(axis=1)
    subtype_result = subtype_result.merge(peak_anno, on="Peak_Name", how="left")
    subtype_result.to_csv(RESULTS / "ATAC_BRCA_subtype_specific_peak_analysis.tsv", sep="\t", index=False)

    sig = subtype_result[(subtype_result["subtype_ANOVA_FDR"] < 0.05) & (subtype_result["max_minus_min_subtype_mean"] >= 0.5)].copy()
    sig.to_csv(RESULTS / "ATAC_BRCA_subtype_specific_peaks_FDR005_delta05.tsv", sep="\t", index=False)
    subtype_counts = sig["best_subtype"].value_counts().rename_axis("best_subtype").reset_index(name="n_subtype_specific_peaks")
    subtype_counts.to_csv(RESULTS / "ATAC_BRCA_subtype_specific_peak_counts.tsv", sep="\t", index=False)


    summary = {
        "n_brca_peaks_with_annotation": int(peak_anno.shape[0]),
        "n_peaks_in_matrix_after_annotation_match": int(matrix.shape[0]),
        "n_subtype_groups_used": len(subtype_groups),
        "n_subtype_specific_peaks_FDR005_delta05": int(sig.shape[0]),
    }
    (RESULTS / "ATAC_peak_annotation_subtype_analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "ATAC_peak_annotation_subtype_analysis_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
