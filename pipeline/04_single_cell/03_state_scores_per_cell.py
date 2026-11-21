from __future__ import annotations

import json
from pathlib import Path


import numpy as np
import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


COLORS = {
    "red": "#c35655",
    "blue": "#5e90b8",
    "light_green": "#b1c5ba",
    "light_orange": "#f5c19b",
    "dark_green": "#6d817d",
    "pink": "#eea5d1",
    "dark_gray": "#3e3c3d",
    "luminal": "#6192b9",
    "her2": "#ee77d2",
    "basal": "#c45958",
}

CELLTYPE_COLORS = {
    "Cancer Epithelial": COLORS["red"],
    "T cells": COLORS["blue"],
    "B cells": COLORS["light_green"],
    "Myeloid": COLORS["dark_green"],
    "Fibroblasts": COLORS["light_orange"],
    "Endothelial": COLORS["pink"],
    "Plasmablasts": "#9e9e9e",
    "Normal Epithelial": "#7f7f7f",
}

SUBTYPE_COLORS = {"ER+": COLORS["luminal"], "HER2+": COLORS["her2"], "TNBC": COLORS["basal"]}

STATE_COLORS = {
    "Luminal_S1": "#6192b9",
    "Luminal_S3": "#b1dfdd",
    "Luminal_S4": "#6d817d",
    "HER2_S3": "#ee77d2",
    "Basal_TNBC_S2": "#e2baad",
    "Basal_TNBC_S4": "#c45958",
    "Low/none": "#d0d0d0",
}

SUBTYPE_TO_CONTEXT = {"ER+": "Luminal", "HER2+": "HER2", "TNBC": "Basal/TNBC"}


ROOT = DATA_ROOT.parent
CORE = CORE_ROOT
OUT = next(CORE.glob("Remap*/figure5-8/figure8"))


def read_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    atlas = pd.read_csv(OUT / "figure8a_gse176078_crc13_scanpy_umap_scores.tsv.gz", sep="\t")
    cell = pd.read_csv(OUT / "figure8i_scrna_guided_crc_cell_expression_by_subtype.tsv.gz", sep="\t")
    sig = pd.read_csv(OUT / "figure8i_scrna_guided_crc_state_signature_table.tsv", sep="\t")
    sig = sig[sig["recommended_for_bulk_projection"].astype(bool)].copy()
    sig["signature_gene_list"] = sig["signature_genes"].fillna("").apply(lambda x: [g for g in str(x).split(";") if g])
    return atlas, cell, sig


def add_state_scores(malignant: pd.DataFrame, cell: pd.DataFrame, sig: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    state_rows = []
    cell = cell.copy()
    malignant = malignant.merge(cell.drop(columns=["orig.ident"], errors="ignore"), on="barcode", how="left", suffixes=("", "_state"))
    for _, row in sig.iterrows():
        genes = [gene for gene in row["signature_gene_list"] if gene in malignant.columns]
        state_col = row["state_id"]
        context = row["subtype_context"]
        values = pd.Series(np.nan, index=malignant.index, dtype=float)
        mask = malignant["subtype"].map(SUBTYPE_TO_CONTEXT).eq(context)
        if genes:
            values.loc[mask] = malignant.loc[mask, genes].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        malignant[state_col] = values
        state_rows.append(
            {
                "state_id": state_col,
                "state_label": row["state_label"],
                "subtype_context": context,
                "n_genes_defined": len(row["signature_gene_list"]),
                "n_genes_used_for_cell_score": len(genes),
                "genes_used_for_cell_score": ";".join(genes),
            }
        )
    state_cols = sig["state_id"].tolist()
    malignant["dominant_crc_state"] = "Low/none"
    malignant["dominant_crc_state_score"] = np.nan
    for idx, row in malignant.iterrows():
        context = SUBTYPE_TO_CONTEXT.get(row["subtype"], None)
        if context is None:
            continue
        context_states = sig.loc[sig["subtype_context"].eq(context), "state_id"].tolist()
        vals = pd.to_numeric(row[context_states], errors="coerce")
        if vals.notna().any():
            best = vals.idxmax()
            best_score = float(vals.max())
            malignant.at[idx, "dominant_crc_state"] = best if best_score > 0 else "Low/none"
            malignant.at[idx, "dominant_crc_state_score"] = best_score
    return malignant, pd.DataFrame(state_rows)


def main() -> None:
    atlas, cell, sig = read_data()
    malignant = atlas[atlas["celltype_major"].astype(str).eq("Cancer Epithelial")].copy()
    malignant, state_coverage = add_state_scores(malignant, cell, sig)
    malignant.to_csv(OUT / "figure8j_malignant_crc_state_umap_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    state_coverage.to_csv(OUT / "figure8j_crc_state_cell_score_gene_coverage.tsv", sep="\t", index=False)
    summary = (
        malignant.groupby(["subtype", "dominant_crc_state"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    summary["fraction_within_subtype"] = summary["n_cells"] / summary.groupby("subtype")["n_cells"].transform("sum")
    summary.to_csv(OUT / "figure8j_dominant_crc_state_cell_counts.tsv", sep="\t", index=False)
    payload = {
        "purpose": "Candidate Figure 7a: retain all-cell atlas context and then filter to malignant epithelial cells to display sample subtype, dominant scRNA-guided CRC state and dominant state score.",
        "important_note": "Malignant-only panels use the existing all-cell UMAP coordinates with non-malignant cells removed; they are malignant-filtered global UMAP views, not a freshly recomputed malignant-only embedding.",
        "outputs": [
            "figure8j_allcell_and_malignant_crc_state_umap.png/pdf/svg/tiff",
            "figure8j_malignant_crc_state_umap_scores.tsv.gz",
            "figure8j_crc_state_cell_score_gene_coverage.tsv",
            "figure8j_dominant_crc_state_cell_counts.tsv",
        ],
    }
    (OUT / "figure8j_allcell_malignant_crc_state_umap_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
