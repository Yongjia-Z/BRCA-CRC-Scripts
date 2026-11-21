from __future__ import annotations

import json
import re
from pathlib import Path


import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = DATA_ROOT.parent
CORE = CORE_ROOT
REM = next(CORE.glob("Remap*/figure5-8"))
FIG7 = REM / "figure7"
FIG8 = REM / "figure8"
DEPMAP = next(ROOT.glob("*/DepMap"))
GDSC_CTRP = next(ROOT.glob("*/GDSC_CTRP"))

COLORS = {
    "red": "#c35655",
    "blue": "#5e90b8",
    "dark_gray": "#3e3c3d",
    "grid": "#d9d9d9",
}


def gene_symbol(col: str) -> str:
    m = re.match(r"^(.+?) \(\d+\)$", col)
    return m.group(1) if m else col


def zscore_columns(df: pd.DataFrame) -> pd.DataFrame:
    x = df.astype(float)
    return x.sub(x.mean(axis=0), axis=1).div(x.std(axis=0, ddof=0).replace(0, np.nan), axis=1).fillna(0)


def load_signatures() -> pd.DataFrame:
    sig = pd.read_csv(FIG8 / "figure8i_scrna_guided_crc_state_signature_table.tsv", sep="\t")
    sig = sig[sig["recommended_for_bulk_projection"].astype(bool)].copy()
    sig["signature_gene_list"] = sig["signature_genes"].fillna("").apply(lambda x: [g for g in str(x).split(";") if g])
    return sig


def depmap_gene_columns(path: Path, genes: list[str]) -> dict[str, str]:
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    mapping = {}
    for col in cols:
        sym = gene_symbol(col)
        if sym in genes and sym not in mapping:
            mapping[sym] = col
    return mapping


def load_breast_models() -> pd.DataFrame:
    model = pd.read_csv(DEPMAP / "data" / "Model.csv")
    breast = model[model["OncotreeLineage"].astype(str).eq("Breast")].copy()
    keep = ["ModelID", "CellLineName", "OncotreeLineage", "OncotreeSubtype", "OncotreeCode", "ModelSubtypeFeatures", "PrimaryOrMetastasis"]
    return breast[keep].drop_duplicates("ModelID")


def load_depmap_expression(signatures: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = sorted({g for gene_list in signatures["signature_gene_list"] for g in gene_list})
    path = DEPMAP / "data" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    mapping = depmap_gene_columns(path, genes)
    usecols = ["ModelID", "IsDefaultEntryForModel", *mapping.values()]
    expr = pd.read_csv(path, usecols=lambda c: c in usecols)
    if "IsDefaultEntryForModel" in expr.columns:
        expr = expr[expr["IsDefaultEntryForModel"].astype(str).str.lower().isin(["true", "yes", "1"])].copy()
    expr = expr.drop(columns=[c for c in ["IsDefaultEntryForModel"] if c in expr.columns])
    expr = expr.rename(columns={v: k for k, v in mapping.items()})
    expr = expr.groupby("ModelID", as_index=False).mean(numeric_only=True)
    breast = load_breast_models()
    expr = breast.merge(expr, on="ModelID", how="inner")
    coverage = []
    score_df = expr[["ModelID", "CellLineName", "OncotreeSubtype", "OncotreeCode", "ModelSubtypeFeatures", "PrimaryOrMetastasis"]].copy()
    gene_expr = expr.set_index("ModelID")[[g for g in genes if g in expr.columns]]
    gene_z = zscore_columns(gene_expr)
    for _, row in signatures.iterrows():
        present = [g for g in row["signature_gene_list"] if g in gene_z.columns]
        missing = [g for g in row["signature_gene_list"] if g not in gene_z.columns]
        score_df[row["state_id"]] = gene_z[present].mean(axis=1).reindex(score_df["ModelID"]).to_numpy() if present else np.nan
        coverage.append({"state_id": row["state_id"], "n_defined": len(row["signature_gene_list"]), "n_present_depmap_expression": len(present), "present_genes": ";".join(present), "missing_genes": ";".join(missing)})
    return score_df, pd.DataFrame(coverage)


def load_depmap_crispr(signatures: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = sorted({g for gene_list in signatures["signature_gene_list"] for g in gene_list})
    path = DEPMAP / "data" / "CRISPRGeneEffect.csv"
    mapping = depmap_gene_columns(path, genes)
    usecols = ["Unnamed: 0", *mapping.values()]
    crispr = pd.read_csv(path, usecols=lambda c: c in usecols)
    crispr = crispr.rename(columns={"Unnamed: 0": "ModelID", **{v: k for k, v in mapping.items()}})
    breast = load_breast_models()
    crispr = breast.merge(crispr, on="ModelID", how="inner")
    coverage = []
    for _, row in signatures.iterrows():
        present = [g for g in row["signature_gene_list"] if g in crispr.columns]
        missing = [g for g in row["signature_gene_list"] if g not in crispr.columns]
        coverage.append({"state_id": row["state_id"], "n_defined": len(row["signature_gene_list"]), "n_present_crispr": len(present), "present_crispr_genes": ";".join(present), "missing_crispr_genes": ";".join(missing)})
        crispr[f"{row['state_id']}_mean_gene_effect"] = crispr[present].mean(axis=1) if present else np.nan
    return crispr, pd.DataFrame(coverage)


def spearman_pair(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    sub = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if sub.shape[0] < 8 or sub["x"].nunique() < 4 or sub["y"].nunique() < 4:
        return int(sub.shape[0]), np.nan, np.nan
    rho, p = stats.spearmanr(sub["x"], sub["y"])
    return int(sub.shape[0]), float(rho), float(p)


def crispr_correlations(scores: pd.DataFrame, crispr: pd.DataFrame, signatures: pd.DataFrame) -> pd.DataFrame:
    merged = scores.merge(crispr[["ModelID", *[f"{sid}_mean_gene_effect" for sid in signatures["state_id"]]]], on="ModelID", how="inner")
    rows = []
    for _, state in signatures.iterrows():
        sid = state["state_id"]
        n, rho, p = spearman_pair(merged[sid], merged[f"{sid}_mean_gene_effect"])
        rows.append(
            {
                "state_id": sid,
                "state_label": state["state_label"],
                "subtype_context": state["subtype_context"],
                "n_breast_cell_lines": n,
                "spearman_r_state_expr_vs_mean_gene_effect": rho,
                "pvalue": p,
                "interpretation": "negative rho means higher state expression associates with stronger CRISPR dependency because lower GeneEffect is more essential",
            }
        )
    out = pd.DataFrame(rows)
    mask = out["pvalue"].notna()
    out.loc[mask, "fdr"] = multipletests(out.loc[mask, "pvalue"], method="fdr_bh")[1]
    return out


def main() -> None:
    signatures = load_signatures()
    scores, expr_cov = load_depmap_expression(signatures)
    crispr, crispr_cov = load_depmap_crispr(signatures)
    scores.to_csv(FIG7 / "figure7f_scrna_state_depmap_breast_expression_scores.tsv", sep="\t", index=False)
    pd.merge(expr_cov, crispr_cov, on="state_id", how="outer").to_csv(FIG7 / "figure7f_scrna_state_depmap_gene_coverage.tsv", sep="\t", index=False)
    crispr_corr = crispr_correlations(scores, crispr, signatures)
    crispr_corr.to_csv(FIG7 / "figure7f_scrna_state_depmap_crispr_correlation.tsv", sep="\t", index=False)
    summary = {
        "method": "DepMap breast cell lines were scored for scRNA-guided CRC state expression using mean z-scored logTPM. CRISPR uses mean DepMap GeneEffect for genes in each state; lower GeneEffect means stronger dependency. Drug response uses GDSC LN_IC50/AUC and CTRP area_under_curve/log2_ec50.",
        "n_breast_models_with_expression": int(scores.shape[0]),
        "outputs": [
            "figure7f_scrna_state_depmap_crispr_correlation.png/pdf/svg/tiff",
            "figure7f_scrna_state_depmap_crispr_correlation.tsv",
            "figure7g_scrna_state_drug_response_correlation_heatmap.png/pdf/svg/tiff",
            "figure7g_scrna_state_drug_response_correlations.tsv",
        ],
    }
    (FIG7 / "figure7f_g_scrna_state_functional_projection_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
