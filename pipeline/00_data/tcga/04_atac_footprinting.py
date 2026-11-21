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

BREAST_CRC_CANDIDATE_TFS = [
    "ESR1",
    "FOXA1",
    "GATA3",
    "XBP1",
    "TFAP2A",
    "TFAP2C",
    "TP63",
    "FOXC1",
    "SOX9",
    "RUNX1",
    "RUNX2",
    "JUN",
    "FOS",
    "NFKB1",
    "RELA",
    "STAT3",
    "TEAD1",
    "MYC",
]


def read_excel_table(path: Path, sheet_name: str, header_keyword: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=sheet_name, header=None)
    header_row = None
    for idx, row in raw.iterrows():
        values = [str(value) for value in row.dropna().tolist()]
        if header_keyword in values:
            header_row = idx
            break
    if header_row is None:
        raise ValueError(f"Cannot find header keyword {header_keyword} in {path.name}:{sheet_name}")
    table = pd.read_excel(path, sheet_name=sheet_name, header=header_row)
    table = table.dropna(how="all")
    return table


def main() -> None:
    matrix = pd.read_csv(ATAC / "brca_peak_Log2Counts_dedup.tsv", sep="\t", nrows=1)
    brca_patient_ids = sorted({col[:12] for col in matrix.columns if col != "sample"})

    s4_path = ATAC / "TCGA_ATAC_Data_S4_motif_enrichment.xlsx"
    s6_path = ATAC / "TCGA_ATAC_Data_S6_footprinting.xlsx"

    sheet_rows = []
    for workbook in [s4_path, s6_path]:
        xls = pd.ExcelFile(workbook)
        for sheet in xls.sheet_names:
            sheet_rows.append({"workbook": workbook.name, "sheet": sheet})
    pd.DataFrame(sheet_rows).to_csv(RESULTS / "ATAC_footprinting_workbook_sheets.tsv", sep="\t", index=False)

    motif = read_excel_table(s4_path, "Motif_Enrichment", "TF_Motif_Name")
    cluster_cols = [c for c in motif.columns if str(c).startswith("Cluster_")]
    for col in cluster_cols:
        motif[col] = pd.to_numeric(motif[col], errors="coerce")
    motif["max_cluster_motif_enrichment"] = motif[cluster_cols].max(axis=1)
    motif["best_cluster"] = motif[cluster_cols].idxmax(axis=1)
    motif_top = motif.sort_values("max_cluster_motif_enrichment", ascending=False).head(50)
    motif_top.to_csv(RESULTS / "ATAC_motif_enrichment_top50.tsv", sep="\t", index=False)

    footprint_depth = read_excel_table(s6_path, "Footprint_Depth", "TF_Motif_Name")
    sample_cols = [c for c in footprint_depth.columns if str(c).startswith("TCGA-")]
    brca_cols = [c for c in sample_cols if c[:12] in brca_patient_ids]
    for col in brca_cols:
        footprint_depth[col] = pd.to_numeric(footprint_depth[col], errors="coerce")

    footprint_summary = footprint_depth[["TF_Motif_Name", "CIS-BP_ID", "TF_Name", "Family_Name"]].copy()
    footprint_summary["n_brca_atac_samples"] = len(brca_cols)
    footprint_summary["brca_mean_footprint_depth"] = footprint_depth[brca_cols].mean(axis=1)
    footprint_summary["brca_median_footprint_depth"] = footprint_depth[brca_cols].median(axis=1)
    footprint_summary["brca_sd_footprint_depth"] = footprint_depth[brca_cols].std(axis=1)
    footprint_summary.to_csv(RESULTS / "ATAC_BRCA_footprint_depth_summary.tsv", sep="\t", index=False)

    deepest = footprint_summary.sort_values("brca_mean_footprint_depth", ascending=True).head(50)
    deepest.to_csv(RESULTS / "ATAC_BRCA_deepest_footprint_top50.tsv", sep="\t", index=False)

    correlation = read_excel_table(s6_path, "Correlation Depth Flank to RNA", "TF_Motif_Name")
    for col in ["depth_correlation", "depth_FDR", "flank_correlation", "flank_FDR"]:
        if col in correlation.columns:
            correlation[col] = pd.to_numeric(correlation[col], errors="coerce")
    correlation.to_csv(RESULTS / "ATAC_footprint_expression_correlation.tsv", sep="\t", index=False)

    candidate_footprint = footprint_summary[
        footprint_summary["TF_Name"].astype(str).str.upper().isin(BREAST_CRC_CANDIDATE_TFS)
    ].copy()
    candidate_correlation = correlation[
        correlation["TF_Name"].astype(str).str.upper().isin(BREAST_CRC_CANDIDATE_TFS)
    ].copy()
    candidate_footprint.to_csv(RESULTS / "ATAC_BRCA_candidate_CRC_TF_footprint.tsv", sep="\t", index=False)
    candidate_correlation.to_csv(RESULTS / "ATAC_candidate_CRC_TF_footprint_expression_correlation.tsv", sep="\t", index=False)

    plot_df = deepest.head(20).iloc[::-1]

    candidate_plot = candidate_correlation.sort_values("flank_correlation", ascending=False).head(20)

    summary = {
        "n_brca_atac_samples_in_matrix": len(brca_patient_ids),
        "n_brca_columns_in_footprint_depth": len(brca_cols),
        "n_motifs_in_motif_enrichment": int(motif.shape[0]),
        "n_motifs_in_footprint_depth": int(footprint_summary.shape[0]),
        "n_candidate_crc_tf_motifs_with_footprint": int(candidate_footprint.shape[0]),
        "n_candidate_crc_tf_motifs_with_correlation": int(candidate_correlation.shape[0]),
    }
    (RESULTS / "ATAC_footprinting_analysis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    pd.DataFrame([summary]).to_csv(RESULTS / "ATAC_footprinting_analysis_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
