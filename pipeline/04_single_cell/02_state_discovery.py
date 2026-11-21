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
from scipy import sparse, stats
from scipy.cluster import hierarchy
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


SUBTYPE_CONFIGS = {
    "Luminal": {
        "obs_subtype": "ER+",
        "formal": ["GATA3", "FOXA1", "ESR1", "RUNX1", "ELF3", "SPDEF", "PBX1", "PGR", "EHF", "AR", "TFAP2C"],
        "supporting": ["SREBF1", "EGR3", "RARA", "ARNT2", "CREB3L4", "XBP1", "HSF1"],
        "color": "#6192b9",
    },
    "HER2": {
        "obs_subtype": "HER2+",
        "formal": ["FOXA1", "TFAP2C", "ELF3", "SPDEF", "EHF", "AR", "RUNX1"],
        "supporting": ["SREBF1", "RARA", "HSF1", "HOXB2", "ESRRA", "FOXP4", "SMAD3", "ETV6", "HMGA1", "STAT3"],
        "color": "#ee77d2",
    },
    "Basal/TNBC": {
        "obs_subtype": "TNBC",
        "formal": ["TP63", "FOXC1", "RUNX1"],
        "supporting": ["HSF1", "FOXP4", "SMAD3", "CEBPB", "NFIX", "ETV6", "HMGA1", "NFE2L3", "EN1", "ELF5", "SOX9", "E2F3", "FOSL1", "STAT3"],
        "extra_context": ["KLF5", "NFIB", "TFCP2L1", "E2F2", "SOX8", "RELA", "TFCP2"],
        "color": "#c45958",
    },
}

COLORS = {
    "red": "#c35655",
    "blue": "#5e90b8",
    "light_green": "#b1c5ba",
    "light_orange": "#f5c19b",
    "dark_green": "#6d817d",
    "pink": "#eea5d1",
    "dark_gray": "#3e3c3d",
    "gray": "#d9d9d9",
}


ROOT = DATA_ROOT.parent
CORE = CORE_ROOT
OUT = next(CORE.glob("Remap*/figure5-8/figure8"))
H5AD = next(ROOT.rglob("GSE176078_Wu_etal_2021_BRCA_scRNASeq.h5ad"))


def ordered_unique(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def candidate_source_table() -> pd.DataFrame:
    rows = []
    for subtype, cfg in SUBTYPE_CONFIGS.items():
        for gene in cfg["formal"]:
            rows.append({"subtype_context": subtype, "TF": gene, "candidate_class": "formal_crc_tf"})
        for gene in cfg["supporting"]:
            rows.append({"subtype_context": subtype, "TF": gene, "candidate_class": "borderline_supporting"})
        for gene in cfg.get("extra_context", []):
            rows.append({"subtype_context": subtype, "TF": gene, "candidate_class": "prior_scrna_context"})
    return pd.DataFrame(rows).drop_duplicates(["subtype_context", "TF", "candidate_class"])


def extract_logcpm(adata: ad.AnnData, genes: list[str]) -> tuple[pd.DataFrame, list[str], list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    missing = [gene for gene in genes if gene not in adata.var_names]
    if not present:
        return pd.DataFrame(index=adata.obs_names), present, missing
    x = adata[:, present].X
    if sparse.issparse(x):
        x = x.toarray()
    x = np.asarray(x, dtype=np.float32)
    total = pd.to_numeric(adata.obs["nCount_RNA"], errors="coerce").fillna(0).to_numpy(dtype=np.float32)
    total = np.maximum(total, 1.0)
    x = np.log1p(x / total[:, None] * 1e4)
    return pd.DataFrame(x, index=adata.obs_names, columns=present), present, missing


def choose_cell_k(x: np.ndarray, n_cells: int) -> tuple[int, pd.DataFrame]:
    if n_cells < 500:
        k_values = [2, 3]
        min_size = 20
    elif n_cells < 2500:
        k_values = [3, 4, 5]
        min_size = 50
    else:
        k_values = [3, 4, 5, 6]
        min_size = max(80, int(n_cells * 0.035))
    rng = np.random.default_rng(20260515)
    rows = []
    best = None
    if n_cells > 5000:
        sample_idx = rng.choice(n_cells, size=5000, replace=False)
    else:
        sample_idx = np.arange(n_cells)
    for k in k_values:
        labels = MiniBatchKMeans(n_clusters=k, random_state=20260515, batch_size=2048, n_init=30).fit_predict(x)
        counts = pd.Series(labels).value_counts()
        try:
            sil = float(silhouette_score(x[sample_idx], labels[sample_idx], metric="euclidean"))
        except Exception:
            sil = np.nan
        row = {"k": k, "silhouette": sil, "min_cluster_size": int(counts.min()), "max_cluster_size": int(counts.max())}
        rows.append(row)
        if counts.min() >= min_size and np.isfinite(sil):
            if best is None or sil > best["silhouette"]:
                best = row
    table = pd.DataFrame(rows)
    if best is None:
        best = table.sort_values(["silhouette", "min_cluster_size"], ascending=False).iloc[0].to_dict()
    return int(best["k"]), table


def cluster_cells(expr: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    genes = expr.columns[(expr > 0).mean(axis=0) >= 0.01].tolist()
    if len(genes) < 2:
        genes = expr.columns.tolist()
    x = expr[genes].to_numpy(dtype=np.float32)
    x = StandardScaler(with_mean=True, with_std=True).fit_transform(x)
    n_components = max(2, min(8, x.shape[1], x.shape[0] - 1))
    x_pca = PCA(n_components=n_components, random_state=20260515).fit_transform(x)
    k, k_table = choose_cell_k(x_pca, x_pca.shape[0])
    labels = MiniBatchKMeans(n_clusters=k, random_state=20260515, batch_size=2048, n_init=60).fit_predict(x_pca)
    return pd.Series(labels.astype(str), index=expr.index), k_table


def cluster_genes(expr: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame, list[str]]:
    genes = expr.columns[(expr > 0).mean(axis=0) >= 0.01].tolist()
    if len(genes) < 2:
        return pd.Series({gene: 1 for gene in expr.columns}), pd.DataFrame(), expr.columns.tolist()
    corr = expr[genes].corr(method="spearman").fillna(0.0)
    dist = 1 - corr
    condensed = dist.where(np.triu(np.ones(dist.shape), k=1).astype(bool)).stack().to_numpy()
    z = hierarchy.linkage(condensed, method="average")
    maxclust = min(5, max(2, len(genes) // 4))
    labels = hierarchy.fcluster(z, t=maxclust, criterion="maxclust")
    gene_cluster = pd.Series(labels, index=genes)
    leaves = hierarchy.leaves_list(z)
    ordered = [genes[i] for i in leaves]
    corr_long = corr.reset_index(names="TF1").melt(id_vars="TF1", var_name="TF2", value_name="spearman_r")
    return gene_cluster, corr_long, ordered


def summarize_gene_metrics(expr: pd.DataFrame, cell_clusters: pd.Series, source: pd.DataFrame, subtype: str, gene_cluster: pd.Series) -> pd.DataFrame:
    rows = []
    for gene in expr.columns:
        by_cluster = []
        for cluster, idx in cell_clusters.groupby(cell_clusters).groups.items():
            vals = expr.loc[idx, gene]
            by_cluster.append(
                {
                    "cell_cluster": str(cluster),
                    "mean_ln_cpm_plus_1": float(vals.mean()),
                    "pct_expressing": float((vals > 0).mean() * 100),
                }
            )
        cluster_df = pd.DataFrame(by_cluster)
        peak = cluster_df.sort_values(["mean_ln_cpm_plus_1", "pct_expressing"], ascending=False).iloc[0]
        class_values = source.loc[(source["subtype_context"].eq(subtype)) & (source["TF"].eq(gene)), "candidate_class"].tolist()
        rows.append(
            {
                "subtype_context": subtype,
                "TF": gene,
                "candidate_classes": ";".join(class_values) if class_values else "not_configured",
                "formal_crc_tf": "formal_crc_tf" in class_values,
                "borderline_supporting": "borderline_supporting" in class_values,
                "prior_scrna_context": "prior_scrna_context" in class_values,
                "gene_correlation_cluster": int(gene_cluster.get(gene, -1)),
                "peak_cell_cluster": str(peak["cell_cluster"]),
                "peak_mean_ln_cpm_plus_1": float(peak["mean_ln_cpm_plus_1"]),
                "peak_pct_expressing": float(peak["pct_expressing"]),
                "overall_mean_ln_cpm_plus_1": float(expr[gene].mean()),
                "overall_pct_expressing": float((expr[gene] > 0).mean() * 100),
            }
        )
    return pd.DataFrame(rows)


def build_signature_table(gene_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (subtype, cluster), part in gene_table.groupby(["subtype_context", "gene_correlation_cluster"], observed=True):
        if cluster < 0:
            continue
        formal = part.loc[part["formal_crc_tf"], "TF"].tolist()
        supporting = part.loc[part["borderline_supporting"], "TF"].tolist()
        extra = part.loc[part["prior_scrna_context"], "TF"].tolist()
        if not formal:
            continue
        all_genes = part["TF"].tolist()
        dominant_cluster = part["peak_cell_cluster"].mode().iloc[0] if not part["peak_cell_cluster"].mode().empty else "NA"
        label_genes = formal[:3] + [g for g in supporting[:2] if g not in formal]
        state_id = f"{subtype.replace('/','_')}_S{int(cluster)}"
        rows.append(
            {
                "state_id": state_id,
                "subtype_context": subtype,
                "gene_correlation_cluster": int(cluster),
                "state_label": f"{subtype} " + "/".join(label_genes) + " axis",
                "dominant_cell_cluster": dominant_cluster,
                "n_genes": len(all_genes),
                "n_formal_crc_tfs": len(formal),
                "formal_crc_tfs": ";".join(formal),
                "borderline_supporting_tfs": ";".join(supporting),
                "prior_scrna_context_tfs": ";".join(extra),
                "signature_genes": ";".join(all_genes),
                "recommended_for_bulk_projection": len(formal) >= 1 and len(all_genes) >= 2,
            }
        )
    return pd.DataFrame(rows).sort_values(["subtype_context", "gene_correlation_cluster"])


def dotplot_table(expr: pd.DataFrame, cell_clusters: pd.Series, subtype: str, gene_order: list[str], gene_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in gene_order:
        if gene not in expr.columns:
            continue
        gt = gene_table[(gene_table["subtype_context"].eq(subtype)) & (gene_table["TF"].eq(gene))].iloc[0]
        for cluster, idx in cell_clusters.groupby(cell_clusters).groups.items():
            vals = expr.loc[idx, gene]
            rows.append(
                {
                    "subtype_context": subtype,
                    "cell_cluster": str(cluster),
                    "TF": gene,
                    "gene_correlation_cluster": int(gt["gene_correlation_cluster"]),
                    "formal_crc_tf": bool(gt["formal_crc_tf"]),
                    "candidate_classes": gt["candidate_classes"],
                    "mean_ln_cpm_plus_1": float(vals.mean()),
                    "pct_expressing": float((vals > 0).mean() * 100),
                    "n_cells": int(len(vals)),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    source = candidate_source_table()
    all_candidates = ordered_unique(source["TF"].tolist())
    adata = ad.read_h5ad(H5AD, backed="r")
    expr_all, present, missing = extract_logcpm(adata, all_candidates)
    obs = adata.obs[["orig.ident", "subtype", "celltype_major", "celltype_minor", "nCount_RNA", "nFeature_RNA"]].copy()
    adata.file.close()
    source["present_in_gse176078_h5ad"] = source["TF"].isin(present)
    source["missing_from_gse176078_h5ad"] = source["TF"].isin(missing)
    source.to_csv(OUT / "figure8i_scrna_guided_crc_candidate_source.tsv", sep="\t", index=False)

    gene_tables = []
    dot_tables = []
    k_tables = []
    correlation_tables = []
    cell_score_tables = []
    for subtype, cfg in SUBTYPE_CONFIGS.items():
        genes = ordered_unique(cfg["formal"] + cfg["supporting"] + cfg.get("extra_context", []))
        genes = [gene for gene in genes if gene in expr_all.columns]
        mask = obs["celltype_major"].astype(str).eq("Cancer Epithelial") & obs["subtype"].astype(str).eq(cfg["obs_subtype"])
        subtype_expr = expr_all.loc[mask, genes].copy()
        subtype_obs = obs.loc[mask].copy()
        if subtype_expr.shape[0] < 100 or subtype_expr.shape[1] < 2:
            continue
        cell_clusters, k_table = cluster_cells(subtype_expr)
        k_table.insert(0, "subtype_context", subtype)
        k_tables.append(k_table)
        gene_cluster, corr_long, gene_order = cluster_genes(subtype_expr)
        if not corr_long.empty:
            corr_long.insert(0, "subtype_context", subtype)
            correlation_tables.append(corr_long)
        gene_table = summarize_gene_metrics(subtype_expr, cell_clusters, source, subtype, gene_cluster)
        gene_tables.append(gene_table)
        dot_tables.append(dotplot_table(subtype_expr, cell_clusters, subtype, gene_order, gene_table))
        cell_scores = pd.DataFrame(
            {
                "barcode": subtype_expr.index,
                "subtype_context": subtype,
                "orig.ident": subtype_obs["orig.ident"].to_numpy(),
                "cell_cluster": cell_clusters.reindex(subtype_expr.index).to_numpy(),
            }
        )
        for gene in genes:
            cell_scores[gene] = subtype_expr[gene].to_numpy(dtype=np.float32)
        cell_score_tables.append(cell_scores)

    gene_table = pd.concat(gene_tables, ignore_index=True)
    signature = build_signature_table(gene_table)
    dot = pd.concat(dot_tables, ignore_index=True)
    gene_table.to_csv(OUT / "figure8i_scrna_guided_crc_state_gene_table.tsv", sep="\t", index=False)
    signature.to_csv(OUT / "figure8i_scrna_guided_crc_state_signature_table.tsv", sep="\t", index=False)
    dot.to_csv(OUT / "figure8i_scrna_guided_crc_state_dotplot.tsv", sep="\t", index=False)
    pd.concat(k_tables, ignore_index=True).to_csv(OUT / "figure8i_scrna_guided_crc_cell_cluster_k_selection.tsv", sep="\t", index=False)
    pd.concat(cell_score_tables, ignore_index=True).to_csv(OUT / "figure8i_scrna_guided_crc_cell_expression_by_subtype.tsv.gz", sep="\t", index=False, compression="gzip")
    if correlation_tables:
        pd.concat(correlation_tables, ignore_index=True).to_csv(OUT / "figure8i_scrna_guided_crc_tf_spearman_correlations.tsv", sep="\t", index=False)
    summary = {
        "dataset": "GSE176078 Wu et al. BRCA scRNA-seq",
        "analysis": "malignant epithelial cells split by clinical subtype; candidate TF expression normalized as ln(CPM+1); cell clusters from candidate TF PCA/KMeans; gene modules from Spearman correlation hierarchical clustering",
        "subtype_configs": SUBTYPE_CONFIGS,
        "n_signatures": int(signature.shape[0]),
        "outputs": [
            "figure8i_scrna_guided_crc_tf_state_dotplot.png/pdf/svg/tiff",
            "figure8i_scrna_guided_crc_state_signature_table.tsv",
            "figure8i_scrna_guided_crc_state_gene_table.tsv",
        ],
    }
    (OUT / "figure8i_scrna_guided_crc_state_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
