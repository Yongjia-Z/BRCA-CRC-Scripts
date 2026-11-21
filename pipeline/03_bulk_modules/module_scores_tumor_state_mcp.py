from __future__ import annotations

import json
from pathlib import Path


import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


CRC_ORDER = ["ESR1", "GATA3", "PGR", "PBX1", "FOXA1", "SPDEF", "ELF3", "EHF", "AR", "TFAP2C", "RUNX1", "TP63", "FOXC1"]
SCORE_ORDER = ["pan_CRC13_score", "luminal_CRC13_score", "her2_CRC13_score", "basal_TNBC_CRC13_score"]
PAM50_ORDER = ["BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal", "BRCA_Normal"]
PAM50_LABELS = {"BRCA_LumA": "LumA", "BRCA_LumB": "LumB", "BRCA_Her2": "HER2", "BRCA_Basal": "Basal", "BRCA_Normal": "Normal"}
PAM50_COLORS = {
    "BRCA_LumA": "#2166AC",
    "BRCA_LumB": "#67A9CF",
    "BRCA_Her2": "#EE77D2",
    "BRCA_Basal": "#C45958",
    "BRCA_Normal": "#6D817D",
}
MODULE_COLORS = {
    "pan_CRC13_score": "#3e3c3d",
    "luminal_CRC13_score": "#6192b9",
    "her2_CRC13_score": "#ee77d2",
    "basal_TNBC_CRC13_score": "#c45958",
}

SIGNATURES = {
    "ER signaling": ["ESR1", "PGR", "GREB1", "TFF1", "XBP1", "FOXA1", "GATA3"],
    "HER2 amplicon": ["ERBB2", "GRB7", "PGAP3", "MIEN1", "STARD3", "TCAP"],
    "EMT": ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "CDH2"],
    "Cell cycle": ["MKI67", "TOP2A", "CCNB1", "CCNE1", "CDK1", "AURKA", "BIRC5"],
    "MYC targets": ["MYC", "NPM1", "ODC1", "NCL", "LDHA", "MCM2", "MCM5"],
    "DNA repair": ["BRCA1", "BRCA2", "RAD51", "MRE11", "CHEK1", "FANCD2"],
    "Hypoxia": ["HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "BNIP3"],
    "Inflammatory response": ["IL6", "CXCL8", "TNF", "NFKB1", "RELA", "STAT1", "CXCL10"],
    "T cell signature": ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "CD8B", "GZMB", "PRF1"],
    "B cell signature": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Macrophage signature": ["CD68", "CD163", "CSF1R", "MRC1", "MSR1"],
    "Stromal signature": ["COL1A1", "COL1A2", "ACTA2", "FAP", "DCN", "LUM"],
    "Immune checkpoint": ["PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2"],
}


def find_one(root: Path, name: str, must: bool = True) -> Path | None:
    matches = list(root.rglob(name))
    if not matches:
        if must:
            raise FileNotFoundError(name)
        return None
    return matches[0]


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    values = df.astype(float)
    mean = values.mean(axis=1)
    std = values.std(axis=1).replace(0, np.nan)
    return values.sub(mean, axis=0).div(std, axis=0).fillna(0)


def zscore_columns(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std()
    if pd.isna(sd) or sd == 0:
        return x * 0
    return (x - x.mean()) / sd


def load_crc_membership(fig1a_path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    table = pd.read_csv(fig1a_path, sep="\t")
    table = table[table["TF"].isin(CRC_ORDER)].copy()
    table["TF"] = pd.Categorical(table["TF"], categories=CRC_ORDER, ordered=True)
    table = table.sort_values("TF")
    modules = {
        "luminal": table.loc[table["Luminal_member"].astype(int) == 1, "TF"].tolist(),
        "her2": table.loc[table["HER2_member"].astype(int) == 1, "TF"].tolist(),
        "basal_TNBC": table.loc[table["Basal/TNBC_member"].astype(int) == 1, "TF"].tolist(),
        "pan": CRC_ORDER,
    }
    return table, modules


def score_module(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in matrix.index]
    if not present:
        return pd.Series(np.nan, index=matrix.columns)
    return matrix.loc[present].mean(axis=0)


def build_crc13_scores(expr_path: Path, activity_path: Path, annotation_path: Path, modules: dict[str, list[str]]) -> pd.DataFrame:
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    activity = pd.read_csv(activity_path, sep="\t", index_col=0)
    expr = expr.loc[[g for g in CRC_ORDER if g in expr.index]]
    activity = activity.loc[[g for g in CRC_ORDER if g in activity.index]]
    common_samples = [s for s in expr.columns if s in activity.columns]
    expr_z = zscore_rows(expr[common_samples])
    activity_z = zscore_rows(activity[common_samples])
    scores = pd.DataFrame({"sample_id": common_samples})
    score_map = {
        "pan": "pan_CRC13",
        "luminal": "luminal_CRC13",
        "her2": "her2_CRC13",
        "basal_TNBC": "basal_TNBC_CRC13",
    }
    for module, genes in modules.items():
        prefix = score_map[module]
        expr_score = score_module(expr_z, genes)
        activity_score = score_module(activity_z, genes)
        scores[f"{prefix}_TF_expression_score"] = scores["sample_id"].map(expr_score.to_dict())
        scores[f"{prefix}_TF_activity_score"] = scores["sample_id"].map(activity_score.to_dict())
        scores[f"{prefix}_score"] = scores[[f"{prefix}_TF_expression_score", f"{prefix}_TF_activity_score"]].mean(axis=1)
    annotation = pd.read_csv(annotation_path, sep="\t")
    scores = scores.merge(annotation, on="sample_id", how="left")
    return scores


def load_signature_expression(root: Path, samples: list[str], signatures: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm_path = find_one(root, "TCGA-BRCA.star_tpm.tsv.gz")
    gene_anno_path = find_one(root, "TCGA_BRCA_expression_gene_annotation_matched.tsv")
    gene_anno = pd.read_csv(gene_anno_path, sep="\t")
    wanted_genes = sorted({g for genes in signatures.values() for g in genes})
    wanted_anno = gene_anno[gene_anno["gene_name"].isin(wanted_genes)].copy()
    wanted_ids = set(wanted_anno["ensembl_gene_id"].astype(str))
    id_to_symbol = dict(zip(wanted_anno["ensembl_gene_id"].astype(str), wanted_anno["gene_name"].astype(str)))
    chunks = []
    usecols = ["Ensembl_ID"] + samples
    for chunk in pd.read_csv(tpm_path, sep="\t", usecols=lambda c: c == "Ensembl_ID" or c in samples, chunksize=5000):
        chunk["ensembl_gene_id"] = chunk["Ensembl_ID"].astype(str).str.replace(r"\..*$", "", regex=True)
        chunk = chunk[chunk["ensembl_gene_id"].isin(wanted_ids)]
        if not chunk.empty:
            chunk["gene_name"] = chunk["ensembl_gene_id"].map(id_to_symbol)
            chunks.append(chunk[["gene_name"] + samples])
    if not chunks:
        raise RuntimeError("No signature genes were found in TCGA TPM matrix.")
    expr = pd.concat(chunks, axis=0).groupby("gene_name")[samples].mean()
    coverage_rows = []
    for name, genes in signatures.items():
        present = [g for g in genes if g in expr.index]
        coverage_rows.append({"signature": name, "n_genes_defined": len(genes), "n_genes_present": len(present), "genes_present": ";".join(present)})
    return expr, pd.DataFrame(coverage_rows)


def signature_scores(expr: pd.DataFrame, signatures: dict[str, list[str]]) -> pd.DataFrame:
    z = zscore_rows(expr)
    out = {}
    for name, genes in signatures.items():
        present = [g for g in genes if g in z.index]
        if present:
            out[name] = z.loc[present].mean(axis=0)
    return pd.DataFrame(out).reset_index().rename(columns={"index": "sample_id"})


def spearman_table(df: pd.DataFrame, x_cols: list[str], y_cols: list[str], x_label: str, y_label: str) -> pd.DataFrame:
    rows = []
    for x_col in x_cols:
        for y_col in y_cols:
            sub = df[[x_col, y_col]].dropna()
            if len(sub) < 20 or sub[x_col].nunique() < 3 or sub[y_col].nunique() < 3:
                rows.append({x_label: x_col, y_label: y_col, "spearman_rho": np.nan, "pvalue": np.nan, "n": len(sub)})
                continue
            rho, p = stats.spearmanr(sub[x_col], sub[y_col])
            rows.append({x_label: x_col, y_label: y_col, "spearman_rho": float(rho), "pvalue": float(p), "n": int(len(sub))})
    table = pd.DataFrame(rows)
    ok = table["pvalue"].notna()
    table["fdr"] = np.nan
    if ok.any():
        table.loc[ok, "fdr"] = multipletests(table.loc[ok, "pvalue"], method="fdr_bh")[1]
    return table


def residualize_against_subtype(df: pd.DataFrame, cols: list[str], subtype_col: str = "SUBTYPE") -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    dummies = pd.get_dummies(df[subtype_col].fillna("missing"), drop_first=False, dtype=float)
    x = np.column_stack([np.ones(len(dummies)), dummies.to_numpy()])
    for col in cols:
        y = pd.to_numeric(df[col], errors="coerce")
        ok = y.notna()
        res = pd.Series(np.nan, index=df.index, dtype=float)
        if ok.sum() >= 20:
            beta = np.linalg.lstsq(x[ok.to_numpy()], y[ok].to_numpy(), rcond=None)[0]
            res.loc[ok] = y.loc[ok].to_numpy() - x[ok.to_numpy()].dot(beta)
        out[col] = res
    return out


def clean_score_label(label: str) -> str:
    return {
        "pan_CRC13_score": "Pan CRC TFs",
        "luminal_CRC13_score": "Luminal CRC TFs",
        "her2_CRC13_score": "HER2 CRC TFs",
        "basal_TNBC_CRC13_score": "Basal/TNBC CRC TFs",
    }.get(label, label)


def plot_score_by_subtype(scores: pd.DataFrame, out_dir: Path) -> None:
    long = scores.melt(id_vars=["sample_id", "SUBTYPE"], value_vars=SCORE_ORDER, var_name="score", value_name="value").dropna()
    long = long[long["SUBTYPE"].isin(PAM50_ORDER)].copy()
    long["score_label"] = long["score"].map(clean_score_label)
    long.to_csv(out_dir / "figure6a_crc13_module_scores_by_pam50.tsv", sep="\t", index=False)


def main() -> None:
    root = Path.cwd()
    out_dir = CORE_ROOT / "remap_reanalysis" / "figure5-8" / "figure6"
    raw_dir = out_dir / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)


    fig1a_path = REPO_ROOT / "data" / "crc_tf_table_modules.tsv"
    membership, modules = load_crc_membership(fig1a_path)
    membership.to_csv(out_dir / "figure6_crc13_module_membership.tsv", sep="\t", index=False)

    expr_path = find_one(root, "TCGA_BRCA_TF_expression_matrix.tsv.gz")
    activity_path = find_one(root, "TCGA_BRCA_DoRothEA_like_TF_activity_matrix.tsv.gz")
    annotation_path = find_one(root, "TCGA_BRCA_expression_sample_annotation.tsv")
    scores = build_crc13_scores(expr_path, activity_path, annotation_path, modules)

    signature_expr, coverage = load_signature_expression(root, scores["sample_id"].tolist(), SIGNATURES)
    coverage.to_csv(out_dir / "figure6_signature_gene_coverage.tsv", sep="\t", index=False)
    sig_scores = signature_scores(signature_expr, SIGNATURES)
    combined = scores.merge(sig_scores, on="sample_id", how="left")

    mcp_path = find_one(root, "TCGA_BRCA_MCPcounter_scores.tsv")
    mcp_raw = pd.read_csv(mcp_path, sep="\t")
    mcp = mcp_raw.set_index("cell_population").T.reset_index().rename(columns={"index": "sample_id"})
    combined = combined.merge(mcp, on="sample_id", how="left")
    combined.to_csv(out_dir / "figure6_crc13_scores_with_tumor_state_mcp.tsv", sep="\t", index=False)

    plot_score_by_subtype(scores, out_dir)

    tumor_corr = spearman_table(combined, SCORE_ORDER, list(SIGNATURES.keys()), "CRC_score", "signature")
    tumor_corr.to_csv(out_dir / "figure6b_crc13_tumor_state_signature_correlation.tsv", sep="\t", index=False)

    mcp_cols = [c for c in mcp.columns if c != "sample_id"]
    mcp_corr = spearman_table(combined, SCORE_ORDER, mcp_cols, "CRC_score", "MCPcounter_population")
    mcp_corr.to_csv(out_dir / "figure6c_crc13_mcp_counter_correlation.tsv", sep="\t", index=False)

    residual_cols = SCORE_ORDER + mcp_cols
    residuals = residualize_against_subtype(combined, residual_cols)
    residuals["sample_id"] = combined["sample_id"]
    residuals["SUBTYPE"] = combined["SUBTYPE"]
    adj = residuals[["sample_id", "SUBTYPE"] + residual_cols].copy()
    adj_corr = spearman_table(adj, SCORE_ORDER, mcp_cols, "CRC_score", "MCPcounter_population")
    adj_corr.to_csv(out_dir / "figure6d_crc13_mcp_counter_pam50_adjusted_correlation.tsv", sep="\t", index=False)

    summary = {
        "crc_tf_order": CRC_ORDER,
        "module_membership": modules,
        "n_samples_scored": int(scores.shape[0]),
        "n_samples_with_pam50": int(scores["SUBTYPE"].notna().sum()),
        "expression_matrix": str(expr_path),
        "activity_matrix": str(activity_path),
        "sample_annotation": str(annotation_path),
        "mcp_counter_scores": str(mcp_path),
        "tumor_state_signature_method": "Mean row z-score over compact curated marker genes from TCGA-BRCA TPM; used for interpretive tumor-state context, not formal Hallmark ssGSEA.",
        "mcp_counter_method": "Formal MCP-counter score matrix previously generated by MCPcounter.estimate; correlations recomputed against CRC13 scores.",
        "pam50_adjustment_method": "Linear residualization of CRC13 scores and MCP-counter populations against PAM50 subtype dummies, followed by Spearman correlation of residuals.",
    }
    (out_dir / "figure6_crc13_recompute_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
