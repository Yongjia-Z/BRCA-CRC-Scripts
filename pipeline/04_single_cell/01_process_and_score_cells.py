from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))

os.environ.setdefault("OMP_NUM_THREADS", "10")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "10")
os.environ.setdefault("MKL_NUM_THREADS", "10")
os.environ.setdefault("NUMBA_NUM_THREADS", "10")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


CRC13 = ["ESR1", "GATA3", "PGR", "PBX1", "FOXA1", "SPDEF", "ELF3", "EHF", "AR", "TFAP2C", "RUNX1", "TP63", "FOXC1"]
MODULES = {
    "Pan": CRC13,
    "Luminal": ["ESR1", "GATA3", "PGR", "PBX1", "FOXA1", "SPDEF", "ELF3", "EHF", "AR", "TFAP2C", "RUNX1"],
    "HER2": ["FOXA1", "SPDEF", "ELF3", "EHF", "AR", "TFAP2C", "RUNX1"],
    "Basal/TNBC": ["RUNX1", "TP63", "FOXC1"],
}
COLORS = {
    "red": "#c35655",
    "blue": "#5e90b8",
    "light_orange": "#f5c19b",
    "dark_green": "#6d817d",
    "pink": "#eea5d1",
    "dark_gray": "#3e3c3d",
    "luminal": "#6192b9",
    "her2": "#ee77d2",
    "basal": "#c45958",
}


ROOT = DATA_ROOT.parent
DATA = DATA_ROOT / "GSE176078_scRNA" / "data" / "Wu_etal_2021_BRCA_scRNASeq"
H5AD = DATA / "GSE176078_Wu_etal_2021_BRCA_scRNASeq.h5ad"
META = DATA / "metadata.csv"
OUT = CORE_ROOT / "remap_reanalysis" / "figure5-8" / "figure8"
OUT.mkdir(parents=True, exist_ok=True)


if hasattr(sc.settings, "n_jobs"):
    sc.settings.n_jobs = 10
sc.settings.verbosity = 2


def attach_metadata(adata: ad.AnnData) -> ad.AnnData:
    meta = pd.read_csv(META).rename(columns={"Unnamed: 0": "barcode"})
    for col in ["barcode", "orig.ident", "subtype", "celltype_major", "celltype_minor"]:
        if col not in meta.columns:
            meta[col] = "Unknown"
    meta = meta[["barcode", "orig.ident", "subtype", "celltype_major", "celltype_minor"]].copy()
    if set(meta["barcode"]).intersection(set(adata.obs_names)):
        meta = meta.set_index("barcode")
        overlap = [col for col in ["orig.ident", "subtype", "celltype_major", "celltype_minor"] if col in adata.obs.columns]
        if overlap:
            adata.obs = adata.obs.drop(columns=overlap)
        adata.obs = adata.obs.join(meta, how="left")
        adata.obs["barcode"] = adata.obs_names
    else:
        meta = meta.iloc[: adata.n_obs]
        for col in meta.columns:
            adata.obs[col] = meta[col].to_numpy()
    for col in ["orig.ident", "subtype", "celltype_major", "celltype_minor"]:
        adata.obs[col] = adata.obs[col].fillna("Unknown").astype(str)
    return adata


def row_mean(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.mean(axis=1)).ravel()
    return np.asarray(x).mean(axis=1)


def add_crc_scores(adata: ad.AnnData) -> pd.DataFrame:
    rows = [{"TF": gene, "present_in_h5ad": gene in adata.var_names} for gene in CRC13]
    presence = pd.DataFrame(rows)
    presence.to_csv(OUT / "figure8a_gse176078_scanpy_umap_crc13_gene_presence.tsv", sep="\t", index=False)
    module_rows = []
    for module, genes in MODULES.items():
        available = [gene for gene in genes if gene in adata.var_names]
        adata.obs[f"{module}_CRC_TF_score"] = row_mean(adata[:, available].X) if available else np.nan
        module_rows.append({"module": module, "available_genes": ";".join(available), "n_available": len(available)})
    return pd.DataFrame(module_rows)


def compute_umap(adata: ad.AnnData) -> str:
    if "X_umap" in adata.obsm and adata.obsm["X_umap"].shape[1] >= 2:
        return "existing X_umap"
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="cell_ranger", subset=False)
    sc.tl.pca(adata, n_comps=30, svd_solver="randomized", use_highly_variable=True, random_state=20260515)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=20260515)
    sc.tl.umap(adata, min_dist=0.35, spread=1.0, random_state=20260515)
    return "scanpy log-normalized HVG/PCA/neighbors/UMAP"


def main() -> None:
    adata = sc.read_h5ad(H5AD)
    if sparse.issparse(adata.X) and not sparse.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
    adata = attach_metadata(adata)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    module_df = add_crc_scores(adata)
    backend = compute_umap(adata)
    obs = adata.obs[["barcode", "orig.ident", "subtype", "celltype_major", "celltype_minor"]].copy()
    obs["UMAP1"] = adata.obsm["X_umap"][:, 0]
    obs["UMAP2"] = adata.obsm["X_umap"][:, 1]
    for module in MODULES:
        obs[f"{module}_CRC_TF_score"] = adata.obs[f"{module}_CRC_TF_score"].to_numpy()
    obs.to_csv(OUT / "figure8a_gse176078_crc13_scanpy_umap_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    module_df.to_csv(OUT / "figure8a_gse176078_crc13_scanpy_umap_modules.tsv", sep="\t", index=False)
    summary = {
        "dataset": "GSE176078",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "embedding_backend": backend,
        "max_plot_cells": 60000,
        "scatter_rasterized_in_pdf_svg": True,
        "crc_tfs": CRC13,
    }
    (OUT / "figure8a_gse176078_crc13_scanpy_umap_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
