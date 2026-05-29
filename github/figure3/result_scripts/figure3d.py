from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import bh_fdr, gene_windows, load_crc_tfs, read_gene_by_sample_matrix, require_table, write_table


CONTEXTS = ["Luminal", "HER2", "Basal/TNBC"]
OUTPUT_COLUMNS = [
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
    "mean_significant",
    "best_significant",
    "mean_p_is_descriptive_only",
    "robust_peak_significant_FDR05",
    "figure3d_significance_rule",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure3d")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--flank-bp", type=int, default=50000)
    return parser.parse_args()


def first_existing(columns: pd.Index, names: list[str]) -> str | None:
    lookup = {str(c).lower(): c for c in columns}
    for name in names:
        if name.lower() in lookup:
            return lookup[name.lower()]
    return None


def tf_module(row: pd.Series) -> str:
    for context in CONTEXTS:
        for suffix in ["member", "state"]:
            col = f"{context}_{suffix}"
            if col not in row.index:
                continue
            value = row[col]
            if pd.isna(value):
                continue
            if str(value).strip() in {"1", "1.0", "P", "S", "True", "true"}:
                return context.replace("/TNBC", "").lower()
    return "pan"


def peak_metadata(peaks: pd.DataFrame, peak_id: str, peak: str, rec: dict) -> tuple[object, object]:
    hit = peaks.loc[peaks["peak_id"].astype(str).eq(str(peak))]
    if hit.empty:
        return np.nan, np.nan
    row = hit.iloc[0]
    type_col = first_existing(peaks.columns, ["peak_type", "feature_type", "genomic_context", "annotation", "peak_annotation"])
    dist_col = first_existing(peaks.columns, ["distance_to_TSS", "distance_to_tss", "distance", "nearest_tss_distance"])
    peak_type = row[type_col] if type_col else np.nan
    if dist_col:
        distance = pd.to_numeric(pd.Series([row[dist_col]]), errors="coerce").iloc[0]
    else:
        start_col = first_existing(peaks.columns, ["start", "peak_start"])
        end_col = first_existing(peaks.columns, ["end", "peak_end"])
        if start_col and end_col and "tss" in rec:
            midpoint = (float(row[start_col]) + float(row[end_col])) / 2.0
            distance = midpoint - float(rec["tss"])
        else:
            distance = np.nan
    return peak_type, distance


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    crc = load_crc_tfs(args.input_dir)
    crc = crc.copy()
    crc["module"] = crc.apply(tf_module, axis=1)
    module_map = crc.set_index("TF")["module"].to_dict()
    genes = crc["TF"].tolist()

    rna = read_gene_by_sample_matrix(args.input_dir / "tcga_rna_expression_matrix.csv.gz")
    atac = read_gene_by_sample_matrix(args.input_dir / "tcga_atac_peak_matrix.csv.gz")
    peaks = require_table(args.input_dir / "atac_peak_annotation.csv", "atac_peak_annotation")
    annotation = require_table(args.input_dir / "gene_annotation.csv", "gene_annotation")

    peak_id = first_existing(peaks.columns, ["peak_id", "peak", "id"]) or peaks.columns[0]
    peaks = peaks.rename(columns={peak_id: "peak_id"})
    chrom_col = first_existing(peaks.columns, ["chrom", "chromosome", "chr"])
    start_col = first_existing(peaks.columns, ["start", "peak_start"])
    end_col = first_existing(peaks.columns, ["end", "peak_end"])
    if not chrom_col or not start_col or not end_col:
        raise ValueError("atac_peak_annotation requires chromosome, start, and end columns")

    windows = gene_windows(annotation, genes, args.flank_bp)
    common = rna.columns.intersection(atac.columns)
    rows = []

    for rec in windows.to_dict("records"):
        tf = rec["TF"]
        if tf not in rna.index:
            rows.append({"TF": tf, "module": module_map.get(tf, "pan"), "omission_reason": "TF RNA not present"})
            continue

        sub = peaks[
            (peaks[chrom_col].astype(str).eq(str(rec["chrom"])))
            & (pd.to_numeric(peaks[end_col], errors="coerce") > rec["region_start"])
            & (pd.to_numeric(peaks[start_col], errors="coerce") < rec["region_end"])
        ]
        peak_ids = [p for p in sub["peak_id"].astype(str) if p in atac.index]
        corrs = []

        for peak in peak_ids:
            x = pd.to_numeric(atac.loc[peak, common], errors="coerce")
            y = pd.to_numeric(rna.loc[tf, common], errors="coerce")
            mask = x.notna() & y.notna()
            if mask.sum() >= 4:
                rho, pvalue = stats.spearmanr(x[mask], y[mask])
                corrs.append((peak, float(rho), float(pvalue)))

        if not corrs:
            rows.append(
                {
                    "TF": tf,
                    "module": module_map.get(tf, "pan"),
                    "window_bp": args.flank_bp,
                    "n_local_peaks": len(peak_ids),
                    "omission_reason": "no matched local ATAC-RNA peak",
                }
            )
            continue

        best = max(corrs, key=lambda item: abs(item[1]))
        best_type, best_distance = peak_metadata(sub, "peak_id", best[0], rec)
        rows.append(
            {
                "TF": tf,
                "module": module_map.get(tf, "pan"),
                "window_bp": args.flank_bp,
                "n_local_peaks": len(peak_ids),
                "n_matched_TCGA_ATAC_RNA_samples": int(len(common)),
                "mean_local_ATAC_RNA_spearman": float(np.mean([item[1] for item in corrs])),
                "mean_local_ATAC_RNA_p": float(np.mean([item[2] for item in corrs])),
                "best_local_peak": best[0],
                "best_local_peak_type": best_type,
                "best_local_peak_distance_to_TSS": best_distance,
                "best_local_peak_RNA_spearman": best[1],
                "best_local_peak_RNA_p": best[2],
            }
        )

    out = pd.DataFrame(rows)
    out["mean_local_ATAC_RNA_FDR"] = bh_fdr(out["mean_local_ATAC_RNA_p"]) if "mean_local_ATAC_RNA_p" in out else np.nan
    out["best_local_peak_RNA_FDR_global"] = (
        bh_fdr(out["best_local_peak_RNA_p"]) if "best_local_peak_RNA_p" in out else np.nan
    )
    out["mean_significant"] = pd.to_numeric(out.get("mean_local_ATAC_RNA_p"), errors="coerce").lt(0.05)
    out["best_significant"] = pd.to_numeric(out.get("best_local_peak_RNA_p"), errors="coerce").lt(0.05)
    out["mean_p_is_descriptive_only"] = True
    out["robust_peak_significant_FDR05"] = pd.to_numeric(out["best_local_peak_RNA_FDR_global"], errors="coerce").lt(0.05)
    out["figure3d_significance_rule"] = (
        "strongest local peak RNA correlation tested by global BH-FDR; mean p-value is not used as a formal test"
    )
    out = out[[c for c in OUTPUT_COLUMNS if c in out.columns]]
    write_table(out, args.output_dir / "figure3d.csv")


if __name__ == "__main__":
    main()
