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


def read_excel_table(path: Path, sheet_name: str, header_keyword: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=sheet_name, header=None)
    header_row = None
    for idx, row in raw.iterrows():
        if header_keyword in [str(x) for x in row.dropna().tolist()]:
            header_row = idx
            break
    if header_row is None:
        raise ValueError(f"Cannot find {header_keyword} in {path.name}:{sheet_name}")
    return pd.read_excel(path, sheet_name=sheet_name, header=header_row).dropna(how="all")


def main() -> None:
    candidate = pd.read_csv(RESULTS / "BRCA_candidate_CRC_TF_initial_table.tsv", sep="\t")
    peak_anno = pd.read_csv(RESULTS / "ATAC_BRCA_peak_annotation.tsv", sep="\t")
    subtype_peaks = pd.read_csv(RESULTS / "ATAC_BRCA_subtype_specific_peak_analysis.tsv", sep="\t")

    motif = read_excel_table(ATAC / "TCGA_ATAC_Data_S4_motif_enrichment.xlsx", "Motif_Enrichment", "TF_Motif_Name")
    cluster_cols = [c for c in motif.columns if str(c).startswith("Cluster_")]
    for col in cluster_cols:
        motif[col] = pd.to_numeric(motif[col], errors="coerce")
    motif["max_cluster_motif_enrichment"] = motif[cluster_cols].max(axis=1)
    motif["best_motif_cluster"] = motif[cluster_cols].idxmax(axis=1)
    motif_best = motif.sort_values("max_cluster_motif_enrichment", ascending=False).drop_duplicates("TF_Name")
    motif_best = motif_best.rename(columns={"TF_Name": "TF"})

    links = read_excel_table(ATAC / "TCGA_ATAC_Data_S7_peak_to_gene_links.xlsx", "BRCA_Links", "Peak_Name")
    links.to_csv(RESULTS / "ATAC_BRCA_peak_to_gene_links.tsv", sep="\t", index=False)

    linked_tf = links[links["Linked_Gene"].isin(candidate["TF"])].copy()
    linked_tf.to_csv(RESULTS / "ATAC_BRCA_peak_to_candidate_TF_gene_links.tsv", sep="\t", index=False)
    linked_summary = (
        linked_tf.groupby("Linked_Gene")
        .agg(
            n_peak_to_TF_gene_links=("Peak_Name", "count"),
            min_peak_to_TF_gene_FDR=("FDR_Post_CNA", "min"),
            max_peak_to_TF_gene_correlation=("Correlation_Post_CNA", "max"),
        )
        .reset_index()
        .rename(columns={"Linked_Gene": "TF"})
    )

    nearest_tf = peak_anno[peak_anno["Gene Symbol"].isin(candidate["TF"])].copy() if "Gene Symbol" in peak_anno.columns else pd.DataFrame()
    if not nearest_tf.empty:
        nearest_summary = (
            nearest_tf.groupby("Gene Symbol")
            .agg(n_nearest_TF_peaks=("Peak_Name", "count"))
            .reset_index()
            .rename(columns={"Gene Symbol": "TF"})
        )
    else:
        nearest_summary = pd.DataFrame(columns=["TF", "n_nearest_TF_peaks"])

    subtype_sig = subtype_peaks[(subtype_peaks["subtype_ANOVA_FDR"] < 0.05) & (subtype_peaks["max_minus_min_subtype_mean"] >= 0.5)]
    subtype_tf_links = linked_tf.merge(subtype_sig[["Peak_Name", "best_subtype", "subtype_ANOVA_FDR", "max_minus_min_subtype_mean"]], on="Peak_Name", how="inner")
    subtype_tf_links.to_csv(RESULTS / "ATAC_subtype_specific_peak_to_candidate_TF_links.tsv", sep="\t", index=False)
    subtype_tf_summary = (
        subtype_tf_links.groupby(["Linked_Gene", "best_subtype"])
        .size()
        .reset_index(name="n_subtype_specific_peak_to_TF_links")
        .rename(columns={"Linked_Gene": "TF"})
    )
    subtype_tf_pivot = subtype_tf_summary.pivot_table(
        index="TF", columns="best_subtype", values="n_subtype_specific_peak_to_TF_links", fill_value=0
    ).reset_index()
    subtype_tf_pivot.columns = ["TF"] + [f"n_peak_links_{c}" for c in subtype_tf_pivot.columns[1:]]

    integrated = candidate.merge(
        motif_best[["TF", "TF_Motif_Name", "Family_Name", "max_cluster_motif_enrichment", "best_motif_cluster"]],
        on="TF",
        how="left",
        suffixes=("", "_motif_enrichment"),
    )
    integrated = integrated.merge(linked_summary, on="TF", how="left")
    integrated = integrated.merge(nearest_summary, on="TF", how="left")
    integrated = integrated.merge(subtype_tf_pivot, on="TF", how="left")

    for col in ["n_peak_to_TF_gene_links", "n_nearest_TF_peaks"]:
        integrated[col] = integrated[col].fillna(0).astype(int)
    peak_link_cols = [c for c in integrated.columns if c.startswith("n_peak_links_BRCA_")]
    for col in peak_link_cols:
        integrated[col] = integrated[col].fillna(0).astype(int)

    integrated["ATAC_motif_enrichment_evidence"] = integrated["max_cluster_motif_enrichment"].fillna(0) >= 10
    integrated["ATAC_peak_to_TF_gene_evidence"] = (integrated["n_peak_to_TF_gene_links"] > 0) & (integrated["min_peak_to_TF_gene_FDR"].fillna(1) < 0.05)
    integrated["ATAC_nearest_TF_peak_evidence"] = integrated["n_nearest_TF_peaks"] > 0
    evidence_cols = [
        "expression_subtype_evidence",
        "activity_subtype_evidence",
        "ATAC_footprint_expression_evidence",
        "ATAC_motif_enrichment_evidence",
        "ATAC_peak_to_TF_gene_evidence",
        "ATAC_nearest_TF_peak_evidence",
    ]
    for col in evidence_cols:
        if col not in integrated.columns:
            integrated[col] = False
    integrated["multiomics_CRC_evidence_count"] = integrated[evidence_cols].astype(bool).sum(axis=1)
    integrated = integrated.sort_values(["multiomics_CRC_evidence_count", "mean_primary_tpm_log"], ascending=[False, False])
    integrated.to_csv(RESULTS / "BRCA_CRC_TF_multiomics_initial_evidence_table.tsv", sep="\t", index=False)

    high = integrated[integrated["multiomics_CRC_evidence_count"] >= 4].copy()
    high.to_csv(RESULTS / "BRCA_CRC_TF_multiomics_high_evidence_initial.tsv", sep="\t", index=False)

    plot_df = integrated.head(30).iloc[::-1]

    summary = {
        "n_candidate_tfs": int(integrated.shape[0]),
        "n_peak_to_gene_links": int(links.shape[0]),
        "n_peak_to_candidate_TF_gene_links": int(linked_tf.shape[0]),
        "n_candidate_tfs_with_peak_to_gene_links": int((integrated["n_peak_to_TF_gene_links"] > 0).sum()),
        "n_high_evidence_initial_tfs_evidence_ge4": int(high.shape[0]),
    }
    (RESULTS / "CRC_TF_ATAC_integration_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "CRC_TF_ATAC_integration_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
