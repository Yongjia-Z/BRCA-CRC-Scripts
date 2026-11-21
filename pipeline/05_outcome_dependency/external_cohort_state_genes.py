from __future__ import annotations

import json
import warnings
from pathlib import Path


import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


COLORS = {
    "Luminal": "#6192b9",
    "HER2": "#ee77d2",
    "Basal/TNBC": "#c45958",
    "axis": "#3e3c3d",
    "grid": "#d9d9d9",
}

ROOT = DATA_ROOT.parent
CORE = CORE_ROOT
REM = next(CORE.glob("Remap*/figure5-8"))
FIG7 = REM / "figure7"
FIG8 = REM / "figure8"
RAW = FIG7 / "raw"
RAW.mkdir(parents=True, exist_ok=True)

TCGA_ENDPOINTS = [
    ("OS", "OS_MONTHS", "OS_STATUS"),
    ("DFS", "DFS_MONTHS", "DFS_STATUS"),
    ("DSS", "DSS_MONTHS", "DSS_STATUS"),
    ("PFS", "PFS_MONTHS", "PFS_STATUS"),
]


def clean_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def parse_event(series: pd.Series) -> pd.Series:
    return series.astype(str).str.startswith("1:").astype(float)


def row_zscore(expr: pd.DataFrame) -> pd.DataFrame:
    x = expr.astype(float)
    return x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)


def load_signatures() -> pd.DataFrame:
    sig = pd.read_csv(FIG8 / "figure8i_scrna_guided_crc_state_signature_table.tsv", sep="\t")
    sig = sig[sig["recommended_for_bulk_projection"].astype(bool)].copy()
    sig["signature_gene_list"] = sig["signature_genes"].fillna("").apply(lambda x: [g for g in str(x).split(";") if g])
    return sig


def score_signatures(expr: pd.DataFrame, signatures: pd.DataFrame, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    z = row_zscore(expr)
    scores = pd.DataFrame({"sample_id": z.columns})
    coverage_rows = []
    for _, row in signatures.iterrows():
        genes = row["signature_gene_list"]
        present = [g for g in genes if g in z.index]
        missing = [g for g in genes if g not in z.index]
        score_col = row["state_id"]
        scores[score_col] = z.loc[present].mean(axis=0).to_numpy() if present else np.nan
        coverage_rows.append(
            {
                "cohort": cohort,
                "state_id": row["state_id"],
                "subtype_context": row["subtype_context"],
                "n_defined": len(genes),
                "n_present": len(present),
                "present_genes": ";".join(present),
                "missing_genes": ";".join(missing),
            }
        )
    return scores, pd.DataFrame(coverage_rows)


def load_tcga_expr(genes: list[str]) -> pd.DataFrame:
    path = next(ROOT.rglob("TCGA_BRCA_TF_expression_matrix.tsv.gz"))
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=1000):
        sub = chunk[chunk["gene_name"].isin(genes)].copy()
        if not sub.empty:
            chunks.append(sub)
    expr = pd.concat(chunks, ignore_index=True).drop_duplicates("gene_name").set_index("gene_name")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr


def load_tcga_clinical() -> pd.DataFrame:
    path = REM / "figure6" / "figure6_crc13_scores_with_tumor_state_mcp.tsv"
    clin = pd.read_csv(path, sep="\t")
    clin = clin[clin["sample_type"].astype(str).eq("Primary Tumor")].copy()
    for endpoint, time_col, status_col in TCGA_ENDPOINTS:
        clin[f"{endpoint}_event"] = parse_event(clin[status_col])
        clin[time_col] = clean_numeric(clin[time_col])
    return clin


def load_metabric_expr(genes: list[str]) -> tuple[pd.DataFrame, str]:
    return pd.read_csv(RAW / "METABRIC_scrna_state_signature_genes_zscores_wide.tsv", sep="\t", index_col=0), "cached_cBioPortal"


def load_metabric_clinical() -> pd.DataFrame:
    clin = pd.read_csv(FIG7 / "figure7a_metabric_crc13_expression_scores.tsv", sep="\t")
    clin["pam50_like"] = clin["CLAUDIN_SUBTYPE"].replace({"Her2": "HER2"})
    clin["OS_event"] = parse_event(clin["OS_STATUS"])
    clin["RFS_event"] = parse_event(clin["RFS_STATUS"])
    clin["OS_MONTHS"] = clean_numeric(clin["OS_MONTHS"])
    clin["RFS_MONTHS"] = clean_numeric(clin["RFS_MONTHS"])
    return clin


def load_gse_expr(genes: list[str]) -> pd.DataFrame:
    out = RAW / "GSE96058_scrna_state_signature_genes_expression.tsv"
    if out.exists() and out.stat().st_size > 0:
        return pd.read_csv(out, sep="\t", index_col=0)
    raw = next(ROOT.rglob("GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"))
    rows = []
    for chunk in pd.read_csv(raw, chunksize=2000):
        gene_col = chunk.columns[0]
        sub = chunk[chunk[gene_col].isin(genes)].copy()
        if not sub.empty:
            rows.append(sub)
    expr = pd.concat(rows, ignore_index=True).rename(columns={pd.concat(rows, ignore_index=True).columns[0]: "gene"}).drop_duplicates("gene").set_index("gene")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.to_csv(out, sep="\t")
    return expr


def load_gse_clinical() -> pd.DataFrame:
    clin = pd.read_csv(FIG7 / "figure7b_gse96058_crc13_expression_scores.tsv", sep="\t")
    if "isRepl" in clin.columns:
        clin = clin[clin["isRepl"].astype(str).str.lower().isin(["false", "0", "nan"])].copy()
    clin["pam50_like"] = clin["pam50_subtype"].replace({"Her2": "HER2"})
    clin["OS_MONTHS"] = clean_numeric(clin["ovrallSurvDays"]) / 30.4375
    clin["OS_event"] = clean_numeric(clin["ovrallSurvEvent"])
    return clin


def subset_for_state(data: pd.DataFrame, state: pd.Series, cohort: str) -> pd.DataFrame:
    subtype = state["subtype_context"]
    if cohort == "TCGA-BRCA":
        if subtype == "Luminal":
            return data[data["SUBTYPE"].isin(["BRCA_LumA", "BRCA_LumB"])].copy()
        if subtype == "HER2":
            return data[data["SUBTYPE"].eq("BRCA_Her2")].copy()
        return data[data["SUBTYPE"].eq("BRCA_Basal")].copy()
    if subtype == "Luminal":
        return data[data["pam50_like"].isin(["LumA", "LumB"])].copy()
    if subtype == "HER2":
        return data[data["pam50_like"].isin(["HER2"])].copy()
    return data[data["pam50_like"].eq("Basal")].copy()


def run_cox(data: pd.DataFrame, time_col: str, event_col: str, score_col: str) -> dict[str, float | str]:
    sub = data[[time_col, event_col, score_col]].copy()
    sub[time_col] = clean_numeric(sub[time_col])
    sub[event_col] = clean_numeric(sub[event_col])
    sub[score_col] = clean_numeric(sub[score_col])
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
    sub = sub[sub[time_col] > 0].copy()
    n = int(sub.shape[0])
    events = int(sub[event_col].sum()) if n else 0
    if n < 30 or events < 5 or sub[score_col].std(ddof=0) <= 0:
        return {"n": n, "events": events, "hr_per_1sd": np.nan, "ci_low": np.nan, "ci_high": np.nan, "pvalue": np.nan, "status": "insufficient"}
    score = ((sub[score_col] - sub[score_col].mean()) / sub[score_col].std(ddof=0)).to_numpy(float)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = PHReg(sub[time_col].to_numpy(float), score[:, None], status=sub[event_col].to_numpy(float), ties="breslow").fit(disp=False)
        beta = float(fit.params[0])
        ci = fit.conf_int()
        return {
            "n": n,
            "events": events,
            "hr_per_1sd": float(np.exp(beta)),
            "ci_low": float(np.exp(float(ci[0, 0]))),
            "ci_high": float(np.exp(float(ci[0, 1]))),
            "pvalue": float(fit.pvalues[0]),
            "status": "ok",
        }
    except Exception as exc:
        return {"n": n, "events": events, "hr_per_1sd": np.nan, "ci_low": np.nan, "ci_high": np.nan, "pvalue": np.nan, "status": f"cox_failed:{type(exc).__name__}"}


def p_text(p: float) -> str:
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return f"{p:.1e}"
    return f"{p:.3f}"


def main() -> None:
    signatures = load_signatures()
    all_genes = sorted({gene for genes in signatures["signature_gene_list"] for gene in genes})

    coverage_all = []
    results = []

    tcga_expr = load_tcga_expr(all_genes)
    tcga_scores, cov = score_signatures(tcga_expr, signatures, "TCGA-BRCA")
    coverage_all.append(cov)
    tcga = load_tcga_clinical().merge(tcga_scores, left_on="sample_id", right_on="sample_id", how="left")
    for _, state in signatures.iterrows():
        sub = subset_for_state(tcga, state, "TCGA-BRCA")
        for endpoint, time_col, status_col in TCGA_ENDPOINTS:
            row = run_cox(sub, time_col, f"{endpoint}_event", state["state_id"])
            results.append({"cohort": "TCGA-BRCA", "endpoint": endpoint, "state_id": state["state_id"], "state_label": state["state_label"], "subtype_context": state["subtype_context"], **row})

    met_expr, met_source = load_metabric_expr(all_genes)
    met_scores, cov = score_signatures(met_expr, signatures, "METABRIC")
    cov["metabric_expression_source"] = met_source
    coverage_all.append(cov)
    met = load_metabric_clinical().merge(met_scores, left_on="sample_id", right_on="sample_id", how="left")
    for _, state in signatures.iterrows():
        sub = subset_for_state(met, state, "METABRIC")
        for endpoint, time_col, event_col in [("OS", "OS_MONTHS", "OS_event"), ("RFS", "RFS_MONTHS", "RFS_event")]:
            row = run_cox(sub, time_col, event_col, state["state_id"])
            results.append({"cohort": "METABRIC", "endpoint": endpoint, "state_id": state["state_id"], "state_label": state["state_label"], "subtype_context": state["subtype_context"], **row})

    gse_expr = load_gse_expr(all_genes)
    gse_scores, cov = score_signatures(gse_expr, signatures, "GSE96058")
    coverage_all.append(cov)
    gse = load_gse_clinical().merge(gse_scores, left_on="sample_id", right_on="sample_id", how="left")
    for _, state in signatures.iterrows():
        sub = subset_for_state(gse, state, "GSE96058")
        row = run_cox(sub, "OS_MONTHS", "OS_event", state["state_id"])
        results.append({"cohort": "GSE96058", "endpoint": "OS", "state_id": state["state_id"], "state_label": state["state_label"], "subtype_context": state["subtype_context"], **row})

    results_df = pd.DataFrame(results)
    coverage_df = pd.concat(coverage_all, ignore_index=True)
    results_df.to_csv(FIG7 / "figure7e_scrna_state_continuous_cox_results.tsv", sep="\t", index=False)
    coverage_df.to_csv(FIG7 / "figure7e_scrna_state_bulk_gene_coverage.tsv", sep="\t", index=False)
    payload = {
        "method": "Project scRNA-guided CRC state signatures into bulk cohorts using mean row-z-scored expression; Cox HR is per 1 SD higher state score.",
        "metabric_expression_source": met_source,
        "n_states": int(signatures.shape[0]),
        "outputs": [
            "figure7e_scrna_state_continuous_cox_forest.png/pdf/svg/tiff",
            "figure7e_scrna_state_continuous_cox_results.tsv",
            "figure7e_scrna_state_bulk_gene_coverage.tsv",
        ],
    }
    (FIG7 / "figure7e_scrna_state_survival_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
