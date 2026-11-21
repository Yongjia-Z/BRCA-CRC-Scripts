from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))
REM58 = CORE_ROOT / "remap_reanalysis" / "figure5-8"
EXPR = DATA_ROOT / "TCGA" / "results" / "TCGA_BRCA_TF_expression_matrix.tsv.gz"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / "source_tables" / "figure8b"


def collapse_tcga_subtype_to_context(value: object) -> str:
    text = str(value).strip().lower().replace("_", "").replace("-", "").replace(" ", "")
    if text in {"brcaluma", "brcalumb", "luma", "lumb", "luminal"}:
        return "Luminal"
    if "her2" in text:
        return "HER2"
    if "basal" in text or "tnbc" in text:
        return "Basal/TNBC"
    return str(value).strip()


def row_zscore(expr: pd.DataFrame) -> pd.DataFrame:
    x = expr.astype(float)
    return x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0.0)


def load_signatures() -> pd.DataFrame:
    sig = pd.read_csv(REM58 / "figure8" / "figure8i_scrna_guided_crc_state_signature_table.tsv", sep="\t")
    sig = sig[sig["recommended_for_bulk_projection"].astype(bool)].copy()
    sig["signature_gene_list"] = sig["signature_genes"].fillna("").map(lambda x: [g for g in str(x).split(";") if g])
    order = {"Luminal": 0, "HER2": 1, "Basal/TNBC": 2}
    sig["context_order"] = sig["subtype_context"].map(order).fillna(9)
    return sig.sort_values(["context_order", "state_id"]).reset_index(drop=True)


def load_tcga_state_model_input(truncation_months: float = 120.0) -> pd.DataFrame:
    sig = load_signatures()
    all_genes = sorted({g for genes in sig["signature_gene_list"] for g in genes})
    chunks = []
    for chunk in pd.read_csv(EXPR, sep="\t", chunksize=1000):
        sub = chunk[chunk["gene_name"].isin(all_genes)].copy()
        if not sub.empty:
            chunks.append(sub)
    expr = pd.concat(chunks, ignore_index=True).drop_duplicates("gene_name").set_index("gene_name")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    z = row_zscore(expr)
    scores = pd.DataFrame({"sample_id": z.columns})
    for _, row in sig.iterrows():
        genes = [g for g in row["signature_gene_list"] if g in z.index]
        scores[row["state_id"]] = z.loc[genes].mean(axis=0).to_numpy() if genes else np.nan
    clin = pd.read_csv(REM58 / "figure6" / "figure6_crc13_scores_with_tumor_state_mcp.tsv", sep="\t")
    clin = clin[clin["sample_type"].astype(str).eq("Primary Tumor")].copy()
    data = clin.merge(scores, on="sample_id", how="left")
    data["SUBTYPE_CONTEXT"] = data["SUBTYPE"].map(collapse_tcga_subtype_to_context)

    def event_from_status(s: pd.Series) -> pd.Series:
        return s.astype(str).str.startswith("1:").astype(int)

    rows = []
    endpoints = [("OS", "OS_MONTHS", "OS_STATUS"), ("PFS", "PFS_MONTHS", "PFS_STATUS")]
    for _, state in sig.iterrows():
        subset = data[data["SUBTYPE_CONTEXT"].eq(state["subtype_context"])].copy()
        for endpoint, time_col, status_col in endpoints:
            temp = subset[["sample_id", state["state_id"], time_col, status_col]].copy()
            temp = temp.rename(columns={state["state_id"]: "state_score", time_col: "time_months", status_col: "status_raw"})
            temp["event"] = event_from_status(temp["status_raw"])
            temp["state_score"] = pd.to_numeric(temp["state_score"], errors="coerce")
            temp["time_months"] = pd.to_numeric(temp["time_months"], errors="coerce")
            temp = temp.dropna(subset=["state_score", "time_months", "event"])
            temp = temp[temp["time_months"].gt(0)].copy()
            temp["score_z"] = (temp["state_score"] - temp["state_score"].mean()) / temp["state_score"].std(ddof=0)
            temp["time_truncated"] = temp["time_months"].clip(upper=truncation_months)
            temp["event_truncated"] = np.where(temp["time_months"].gt(truncation_months), 0, temp["event"])
            temp["cohort"] = "TCGA-BRCA"
            temp["endpoint"] = endpoint
            temp["state_id"] = state["state_id"]
            temp["state_label"] = state["state_label"]
            temp["subtype_context"] = state["subtype_context"]
            rows.append(temp)
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    model_input = load_tcga_state_model_input(120.0)
    inp = OUT / "figure8b_tcga_cox_model_input.tsv"
    res = OUT / "figure8b_tcga_os_pfs_cox_ph_truncated_results.tsv"
    model_input.to_csv(inp, sep="\t", index=False)
    subprocess.run(["Rscript", str(Path(__file__).with_name("tcga_state_cox.R")), str(inp), str(res)], check=True)
    print("rows", pd.read_csv(res, sep="\t").shape[0], "->", res)


if __name__ == "__main__":
    main()
