from __future__ import annotations

import gc
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))

N_JOBS = os.cpu_count() or 8
os.environ["OMP_NUM_THREADS"] = str(N_JOBS)
os.environ["OPENBLAS_NUM_THREADS"] = str(N_JOBS)
os.environ["MKL_NUM_THREADS"] = str(N_JOBS)
os.environ["NUMEXPR_NUM_THREADS"] = str(N_JOBS)

import anndata as ad
import pandas as pd
from scipy.io import mmread

ROOT = DATA_ROOT.parent
DATA = DATA_ROOT / "GSE176078_scRNA" / "data" / "Wu_etal_2021_BRCA_scRNASeq"
H5AD = DATA / "GSE176078_Wu_etal_2021_BRCA_scRNASeq.h5ad"


def main() -> None:
    if H5AD.exists() and H5AD.stat().st_size > 0:
        print(f"exists\t{H5AD}\t{H5AD.stat().st_size}")
        return
    genes = pd.read_csv(DATA / "count_matrix_genes.tsv", header=None)[0].astype(str).tolist()
    meta = pd.read_csv(DATA / "metadata.csv").rename(columns={"Unnamed: 0": "barcode"})
    with open(DATA / "count_matrix_sparse.mtx", "rb") as fh:
        x = mmread(fh).tocsr()
    x = x.T.tocsr()
    adata = ad.AnnData(X=x)
    adata.var_names = pd.Index(genes)
    adata.obs_names = pd.Index(meta["barcode"].astype(str))
    for col in meta.columns:
        if col != "barcode":
            adata.obs[col] = meta[col].values
    adata.write_h5ad(H5AD, compression="gzip")
    print(f"created\t{H5AD}\t{adata.n_obs}\t{adata.n_vars}\t{H5AD.stat().st_size}")
    del adata, x, genes, meta
    gc.collect()


if __name__ == "__main__":
    main()
