from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import os

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


FIG_DIR = CORE_ROOT / "figure3"
ROOT = DATA_ROOT.parent
DATA_FIG_DIR = ROOT / "data" / "figure3"
RAW = FIG_DIR / "raw"
TCGA_DIR = DATA_ROOT / "TCGA"
TCGA_RESULTS = TCGA_DIR / "results"
TCGA_ATAC = TCGA_DIR / "data" / "ATAC"

PALETTE = {
    "red": "#c35655",
    "blue": "#5e90b8",
    "light_green": "#b1c5ba",
    "light_orange": "#f5c19b",
    "dark_green": "#6d817d",
    "pink": "#eea5d1",
    "dark_gray": "#3e3c3d",
    "mid_gray": "#bdbdbd",
    "light_gray": "#f2f2f2",
    "white": "#ffffff",
}

MODULE_COLORS = {
    "luminal": PALETTE["blue"],
    "her2": PALETTE["pink"],
    "basal": PALETTE["red"],
}

CORE_ROWS = [
    ("ESR1", "luminal"),
    ("FOXA1", "luminal"),
    ("SPDEF", "luminal"),
    ("GATA3", "luminal"),
    ("PBX1", "luminal"),
    ("AR", "luminal"),
    ("CREB3L4", "luminal"),
    ("TFAP2C", "her2"),
    ("EHF", "her2"),
    ("MYBL2", "her2"),
    ("NR2F6", "her2"),
    ("STAT1", "her2"),
    ("ELF3", "her2"),
    ("RUNX1", "basal"),
    ("TP63", "basal"),
    ("FOXC1", "basal"),
]
CORE_TFS = [x[0] for x in CORE_ROWS]
MODULE_BY_TF = dict(CORE_ROWS)

PEAKSET_ORDER = [
    "BRCA_cancer_control",
    "Breast_BASAL_healthy",
    "Breast_LuminalHR_healthy",
    "Breast_LuminalSEC_healthy",
]
PEAKSET_LABELS = {
    "BRCA_cancer_control": "BRCA cancer/control",
    "Breast_BASAL_healthy": "Healthy basal",
    "Breast_LuminalHR_healthy": "Healthy LumHR",
    "Breast_LuminalSEC_healthy": "Healthy LumSEC",
}
PEAKTYPE_ORDER = ["Promoter", "Intronic", "Distal", "Exonic", "Other"]
PEAKTYPE_COLORS = {
    "Promoter": PALETTE["dark_gray"],
    "Intronic": PALETTE["light_green"],
    "Distal": PALETTE["light_orange"],
    "Exonic": "#dedede",
    "Other": PALETTE["mid_gray"],
}


def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[valid] = restored
    return out


def spearman_pair(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 8 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
        return np.nan, np.nan
    r, p = spearmanr(x[ok], y[ok])
    return float(r), float(p)


def build_figure3d():
    window_bp = 50000
    peak_anno = pd.read_csv(
        TCGA_RESULTS / "ATAC_BRCA_peak_annotation_with_accessibility.tsv",
        sep="\t",
        usecols=["Peak_Name", "Gene Symbol", "Annotation", "Distance To Nearest TSS"],
    )
    peak_anno = peak_anno[
        peak_anno["Gene Symbol"].isin(CORE_TFS)
        & (pd.to_numeric(peak_anno["Distance To Nearest TSS"], errors="coerce").abs() <= window_bp)
    ].copy()
    peak_anno["Distance To Nearest TSS"] = pd.to_numeric(peak_anno["Distance To Nearest TSS"], errors="coerce")
    needed_peaks = set(peak_anno["Peak_Name"].astype(str))

    atac_header = pd.read_csv(TCGA_ATAC / "brca_peak_Log2Counts_dedup.tsv", sep="\t", nrows=0).columns.tolist()
    peak_col = atac_header[0]
    atac_samples = atac_header[1:]
    expr_header = pd.read_csv(TCGA_RESULTS / "TCGA_BRCA_TF_expression_matrix.tsv.gz", sep="\t", nrows=0).columns.tolist()
    expr_samples = expr_header[1:]
    expr_by_patient = {}
    for col in expr_samples:
        if col[13:15] == "01":
            expr_by_patient.setdefault(col[:12], col)
    sample_pairs = []
    for atac_sample in atac_samples:
        expr_sample = atac_sample if atac_sample in expr_samples else expr_by_patient.get(atac_sample[:12])
        if expr_sample is not None:
            sample_pairs.append((atac_sample, expr_sample, atac_sample[:12]))
    atac_cols = [x[0] for x in sample_pairs]
    expr_cols = [x[1] for x in sample_pairs]

    chunks = []
    for chunk in pd.read_csv(TCGA_ATAC / "brca_peak_Log2Counts_dedup.tsv", sep="\t", chunksize=50000):
        keep = chunk[peak_col].astype(str).isin(needed_peaks)
        if keep.any():
            chunks.append(chunk.loc[keep, [peak_col] + atac_cols].copy())
    if not chunks:
        raise ValueError("No local ATAC peaks were found for the 16 core CRC TFs.")
    atac = pd.concat(chunks, ignore_index=True).rename(columns={peak_col: "Peak_Name"})
    atac["Peak_Name"] = atac["Peak_Name"].astype(str)
    atac = atac.set_index("Peak_Name")
    atac[atac_cols] = atac[atac_cols].apply(pd.to_numeric, errors="coerce")

    expr = pd.read_csv(TCGA_RESULTS / "TCGA_BRCA_TF_expression_matrix.tsv.gz", sep="\t")
    expr = expr[expr["gene_name"].isin(CORE_TFS)].set_index("gene_name")
    expr[expr_cols] = expr[expr_cols].apply(pd.to_numeric, errors="coerce")

    peak_tests = []
    rows = []
    local_by_tf = {tf: peak_anno[peak_anno["Gene Symbol"] == tf].copy() for tf in CORE_TFS}
    for tf in CORE_TFS:
        local = local_by_tf[tf]
        local_peak_ids = [p for p in local["Peak_Name"].astype(str).tolist() if p in atac.index]
        expr_vec = expr.loc[tf, expr_cols].to_numpy(dtype=float)
        local_matrix = atac.loc[local_peak_ids, atac_cols] if local_peak_ids else pd.DataFrame(columns=atac_cols)
        mean_vec = local_matrix.mean(axis=0).to_numpy(dtype=float) if local_peak_ids else np.full(len(atac_cols), np.nan)
        mean_r, mean_p = spearman_pair(mean_vec, expr_vec)
        row = {
            "TF": tf,
            "module": MODULE_BY_TF[tf],
            "window_bp": window_bp,
            "n_local_peaks": len(local_peak_ids),
            "n_matched_TCGA_ATAC_RNA_samples": len(sample_pairs),
            "mean_local_ATAC_RNA_spearman": mean_r,
            "mean_local_ATAC_RNA_p": mean_p,
            "best_local_peak": np.nan,
            "best_local_peak_type": np.nan,
            "best_local_peak_distance_to_TSS": np.nan,
            "best_local_peak_RNA_spearman": np.nan,
            "best_local_peak_RNA_p": np.nan,
        }
        per_peak = []
        for peak_id in local_peak_ids:
            r, p = spearman_pair(atac.loc[peak_id, atac_cols].to_numpy(dtype=float), expr_vec)
            annot = local[local["Peak_Name"].astype(str) == peak_id].iloc[0]
            test = {
                "TF": tf,
                "Peak_Name": peak_id,
                "Peak_Type": annot["Annotation"],
                "Distance_To_Nearest_TSS": annot["Distance To Nearest TSS"],
                "peak_RNA_spearman": r,
                "peak_RNA_p": p,
            }
            per_peak.append(test)
            peak_tests.append(test)
        valid_peak = [x for x in per_peak if np.isfinite(x["peak_RNA_spearman"])]
        if valid_peak:
            best = max(valid_peak, key=lambda x: abs(x["peak_RNA_spearman"]))
            row.update(
                {
                    "best_local_peak": best["Peak_Name"],
                    "best_local_peak_type": best["Peak_Type"],
                    "best_local_peak_distance_to_TSS": best["Distance_To_Nearest_TSS"],
                    "best_local_peak_RNA_spearman": best["peak_RNA_spearman"],
                    "best_local_peak_RNA_p": best["peak_RNA_p"],
                }
            )
        rows.append(row)

    out = pd.DataFrame(rows)
    out["mean_local_ATAC_RNA_FDR"] = bh_fdr(out["mean_local_ATAC_RNA_p"].to_numpy())
    peak_test_df = pd.DataFrame(peak_tests)
    if not peak_test_df.empty:
        peak_test_df["peak_RNA_FDR_global"] = bh_fdr(peak_test_df["peak_RNA_p"].to_numpy())
        peak_test_df.to_csv(RAW / "figure3d__all_local_peak_RNA_correlations.tsv", sep="\t", index=False)
        fdr_map = peak_test_df.set_index(["TF", "Peak_Name"])["peak_RNA_FDR_global"].to_dict()
        out["best_local_peak_RNA_FDR_global"] = [
            fdr_map.get((row.TF, row.best_local_peak), np.nan) for row in out.itertuples(index=False)
        ]
    else:
        out["best_local_peak_RNA_FDR_global"] = np.nan
    out = out[
        [
            "TF",
            "module",
            "window_bp",
            "n_local_peaks",
            "n_matched_TCGA_ATAC_RNA_samples",
            "mean_local_ATAC_RNA_spearman",
            "mean_local_ATAC_RNA_p",
            "mean_local_ATAC_RNA_FDR",
            "best_local_peak",
            "best_local_peak_type",
            "best_local_peak_distance_to_TSS",
            "best_local_peak_RNA_spearman",
            "best_local_peak_RNA_p",
            "best_local_peak_RNA_FDR_global",
        ]
    ]
    out.to_csv(FIG_DIR / "figure3d.tsv", sep="\t", index=False)
    out.to_csv(RAW / "figure3d__core_TF_local_ATAC_RNA_coupling_50kb.tsv", sep="\t", index=False)
    shutil.copy2(FIG_DIR / "figure3d.tsv", DATA_FIG_DIR / "figure3d.tsv")

    columns = ["mean_local_ATAC_RNA_spearman", "best_local_peak_RNA_spearman"]
    labels = ["Mean local ATAC\nvs RNA", "Strongest local peak\nvs RNA"]
    fdr_cols = ["mean_local_ATAC_RNA_FDR", "best_local_peak_RNA_FDR_global"]
    matrix = out.set_index("TF")[columns].astype(float).reindex(CORE_TFS)
    fdr = out.set_index("TF")[fdr_cols].astype(float).reindex(CORE_TFS)
    n_rows, n_cols = matrix.shape

    for i, gene in enumerate(CORE_TFS):
        y = n_rows - i - 1
        for j, col in enumerate(columns):
            val = matrix.iloc[i, j]
            if pd.isna(val):
                color = PALETTE["light_gray"]
                label = "NA"
                txt_color = PALETTE["dark_gray"]
            else:
                label = f"{val:.2f}"
                if fdr.iloc[i, j] < 0.05:
                    label += "*"
                txt_color = PALETTE["white"] if abs(val) > 0.42 else PALETTE["dark_gray"]


def main():
    RAW.mkdir(exist_ok=True)
    DATA_FIG_DIR.mkdir(parents=True, exist_ok=True)
    build_figure3d()
    print(f"Rebuilt Figure 3A-D in {FIG_DIR}")


if __name__ == "__main__":
    main()
