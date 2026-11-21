from __future__ import annotations

from pathlib import Path


import numpy as np
import pandas as pd
from scipy.stats import kruskal
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = DATA_ROOT.parent
REM = CORE_ROOT / "remap_reanalysis" / "figure5-8"
FIG6 = REM / "figure6"
FIG7 = REM / "figure7"
FIG8 = REM / "figure8"
SIG = FIG8 / "figure8i_scrna_guided_crc_state_signature_table.tsv"
OUT = FIG7 / "figure8a_tcga_state_score_by_pam50"

PAM50_ORDER = ["LumA", "LumB", "HER2", "Basal", "Normal"]
PAM50_COLORS = {
    "LumA": "#2166AC",
    "LumB": "#67A9CF",
    "HER2": "#ee77d2",
    "Basal": "#c45958",
    "Normal": "#6d817d",
}
STATE_ORDER = ["Luminal_S1", "Luminal_S3", "Luminal_S4", "HER2_S3", "Basal_TNBC_S2", "Basal_TNBC_S4"]
STATE_LABELS = {
    "Luminal_S1": "Luminal_S1\nGATA3/FOXA1/ESR1",
    "Luminal_S3": "Luminal_S3\nSPDEF/EHF/AR",
    "Luminal_S4": "Luminal_S4\nRUNX1/ELF3/TFAP2C",
    "HER2_S3": "HER2_S3\nFOXA1/TFAP2C/ELF3",
    "Basal_TNBC_S2": "Basal_S2\nTP63 axis",
    "Basal_TNBC_S4": "Basal_S4\nFOXC1/RUNX1 axis",
}


def p_label(p: float) -> str:
    if not np.isfinite(p):
        return "Kruskal p=NA"
    if p < 0.001:
        return "Kruskal p<0.001"
    return f"Kruskal p={p:.3f}"


def row_zscore(expr: pd.DataFrame) -> pd.DataFrame:
    x = expr.astype(float)
    return x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)


def load_signatures() -> pd.DataFrame:
    sig = pd.read_csv(SIG, sep="\t")
    sig = sig[sig["recommended_for_bulk_projection"].astype(bool)].copy()
    sig = sig[sig["state_id"].isin(STATE_ORDER)].copy()
    sig["state_id"] = pd.Categorical(sig["state_id"], STATE_ORDER, ordered=True)
    sig = sig.sort_values("state_id")
    sig["signature_gene_list"] = sig["signature_genes"].fillna("").apply(lambda x: [g for g in str(x).split(";") if g])
    return sig


def load_tcga_expr(genes: list[str]) -> pd.DataFrame:
    path = next(ROOT.rglob("TCGA_BRCA_TF_expression_matrix.tsv.gz"))
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=1000):
        sub = chunk[chunk["gene_name"].isin(genes)].copy()
        if not sub.empty:
            chunks.append(sub)
    expr = pd.concat(chunks, ignore_index=True).drop_duplicates("gene_name").set_index("gene_name")
    return expr.apply(pd.to_numeric, errors="coerce")


def main() -> None:
    sig = load_signatures()
    genes = sorted({g for xs in sig["signature_gene_list"] for g in xs})
    expr = load_tcga_expr(genes)
    z = row_zscore(expr)
    scores = pd.DataFrame({"sample_id": z.columns})
    coverage = []
    for _, row in sig.iterrows():
        present = [g for g in row["signature_gene_list"] if g in z.index]
        missing = [g for g in row["signature_gene_list"] if g not in z.index]
        scores[row["state_id"]] = z.loc[present].mean(axis=0).to_numpy() if present else np.nan
        coverage.append(
            {
                "state_id": row["state_id"],
                "state_label": row["state_label"],
                "n_defined": len(row["signature_gene_list"]),
                "n_present": len(present),
                "present_genes": ";".join(present),
                "missing_genes": ";".join(missing),
            }
        )
    clin = pd.read_csv(FIG6 / "figure6_crc13_scores_with_tumor_state_mcp.tsv", sep="\t")
    clin = clin[clin["sample_type"].astype(str).eq("Primary Tumor")].copy()
    clin["pam50"] = clin["SUBTYPE"].astype(str).str.replace("BRCA_", "", regex=False).replace({"Her2": "HER2"})
    data = clin[["sample_id", "pam50"]].merge(scores, on="sample_id", how="inner")
    data = data[data["pam50"].isin(PAM50_ORDER)].copy()
    long = data.melt(id_vars=["sample_id", "pam50"], value_vars=STATE_ORDER, var_name="state_id", value_name="state_score")
    long["state_label"] = long["state_id"].map(STATE_LABELS)
    long.to_csv(OUT.with_name(OUT.name + "_long.tsv"), sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(OUT.with_name(OUT.name + "_gene_coverage.tsv"), sep="\t", index=False)

    p_rows = []
    for state in STATE_ORDER:
        sub = long[long["state_id"].eq(state)]
        groups = [sub.loc[sub["pam50"].eq(p), "state_score"].dropna().to_numpy() for p in PAM50_ORDER]
        stat, p = kruskal(*groups)
        p_rows.append({"state_id": state, "kruskal_statistic": stat, "pvalue": p})
    ptab = pd.DataFrame(p_rows)
    ptab.to_csv(OUT.with_name(OUT.name + "_kruskal_tests.tsv"), sep="\t", index=False)

    rng = np.random.default_rng(20260515)


if __name__ == "__main__":
    main()
