from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TCGA")
DATA = ROOT / "data"
ATAC = DATA / "ATAC"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def main() -> None:
    matrix_path = ATAC / "brca_peak_Log2Counts_dedup.tsv"
    subtype_path = DATA / "cbioportal_brca_patient_clinical_wide.tsv"

    matrix = pd.read_csv(matrix_path, sep="\t")
    matrix = matrix.rename(columns={"sample": "peak_id"})
    sample_cols = [c for c in matrix.columns if c != "peak_id"]

    numeric = matrix[sample_cols]
    peak_summary = pd.DataFrame(
        {
            "peak_id": matrix["peak_id"],
            "peak_mean_accessibility": numeric.mean(axis=1),
            "peak_sd_accessibility": numeric.std(axis=1),
            "peak_min_accessibility": numeric.min(axis=1),
            "peak_max_accessibility": numeric.max(axis=1),
        }
    )
    peak_summary.to_csv(RESULTS / "ATAC_peak_accessibility_summary.tsv", sep="\t", index=False)

    sample_summary = pd.DataFrame(
        {
            "sample_id": sample_cols,
            "patient_id": [sample[:12] for sample in sample_cols],
            "sample_mean_accessibility": numeric.mean(axis=0).values,
            "sample_total_accessibility": numeric.sum(axis=0).values,
            "sample_sd_accessibility": numeric.std(axis=0).values,
            "n_peaks_positive": (numeric > 0).sum(axis=0).values,
        }
    )

    if subtype_path.exists():
        subtype = pd.read_csv(subtype_path, sep="\t")
        if "SUBTYPE" in subtype.columns:
            sample_summary = sample_summary.merge(
                subtype[["patientId", "SUBTYPE"]],
                left_on="patient_id",
                right_on="patientId",
                how="left",
            ).drop(columns=["patientId"])

    sample_summary.to_csv(RESULTS / "ATAC_sample_accessibility_summary.tsv", sep="\t", index=False)

    top_variable = peak_summary.sort_values("peak_sd_accessibility", ascending=False).head(1000)
    top_variable.to_csv(RESULTS / "ATAC_top1000_variable_peaks.tsv", sep="\t", index=False)

    matrix_summary = {
        "n_peaks": int(matrix.shape[0]),
        "n_atac_samples": int(len(sample_cols)),
        "n_samples_with_subtype": int(sample_summary.get("SUBTYPE", pd.Series(dtype=object)).notna().sum())
        if "SUBTYPE" in sample_summary.columns
        else 0,
        "mean_peak_accessibility": float(peak_summary["peak_mean_accessibility"].mean()),
        "mean_sample_total_accessibility": float(sample_summary["sample_total_accessibility"].mean()),
    }
    (RESULTS / "ATAC_peak_matrix_summary.json").write_text(
        json.dumps(matrix_summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    pd.DataFrame([matrix_summary]).to_csv(RESULTS / "ATAC_peak_matrix_summary.tsv", sep="\t", index=False)


    plot_df = sample_summary.sort_values("sample_total_accessibility", ascending=False)
    colors = "#F28E2B"

    if "SUBTYPE" in sample_summary.columns and sample_summary["SUBTYPE"].notna().any():
        subtype_order = sample_summary.groupby("SUBTYPE")["sample_total_accessibility"].median().sort_values(ascending=False).index
        grouped = [sample_summary.loc[sample_summary["SUBTYPE"] == subtype, "sample_total_accessibility"] for subtype in subtype_order]

    print(json.dumps(matrix_summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
