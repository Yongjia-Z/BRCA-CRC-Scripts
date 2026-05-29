from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import bh_fdr, load_crc_tfs, require_table, row_zscore, write_table


OUTPUT_COLUMNS = [
    "source_table",
    "state_id",
    "state_label",
    "subtype_context",
    "n_breast_cell_lines",
    "spearman_r_state_expr_vs_dependency_strength",
    "spearman_r_state_expr_vs_mean_gene_effect",
    "pvalue",
    "interpretation",
    "fdr",
    "gene",
    "used_in_crc_tf_heatmap",
    "present_in_depmap_crispr",
    "depmap_column",
    "row_label",
    "Luminal_state",
    "HER2_state",
    "Basal/TNBC_state",
    "Luminal_member",
    "HER2_member",
    "Basal/TNBC_member",
    "ER+/luminal",
    "HER2+",
    "TNBC/basal",
    "Other/unknown",
]
CORRELATION_SOURCE = "figure8c_state_expression_vs_crispr_dependency.csv"
COVERAGE_SOURCE = "figure8c_depmap_crispr_gene_coverage.csv"
DEPENDENCY_SOURCE = "figure8c_depmap_crispr_gene_dependency_strength.csv"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure8c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def gene_symbol(col: str) -> str:
    match = re.match(r"^(.+?) \(\d+\)$", str(col))
    return match.group(1) if match else str(col)


def depmap_gene_columns(cols: list[str], genes: list[str]) -> dict[str, str]:
    mapping = {}
    for col in cols:
        symbol = gene_symbol(col)
        if symbol in genes and symbol not in mapping:
            mapping[symbol] = col
    return mapping


def load_signatures(input_dir: Path) -> pd.DataFrame:
    sig = require_table(input_dir / "crc_state_signatures.csv", "crc_state_signatures")
    rows = []
    for state_id, group in sig.groupby("state_id", sort=False):
        rows.append(
            {
                "state_id": state_id,
                "subtype_context": group["subtype_context"].dropna().iloc[0]
                if "subtype_context" in group.columns and group["subtype_context"].notna().any()
                else "Luminal",
                "state_label": group["state_label"].dropna().iloc[0]
                if "state_label" in group.columns and group["state_label"].notna().any()
                else state_id,
                "signature_gene_list": group["gene"].dropna().astype(str).tolist(),
            }
        )
    return pd.DataFrame(rows)


def load_breast_models(meta: pd.DataFrame) -> pd.DataFrame:
    lineage_col = "OncotreeLineage" if "OncotreeLineage" in meta.columns else None
    breast = meta[meta[lineage_col].astype(str).eq("Breast")].copy() if lineage_col else meta.copy()
    keep = [
        c
        for c in ["ModelID", "CellLineName", "OncotreeSubtype", "OncotreeCode", "ModelSubtypeFeatures", "lineage_subtype"]
        if c in breast.columns
    ]
    breast = breast[keep].drop_duplicates("ModelID").copy()
    subtype_text = breast.drop(columns=["ModelID"]).fillna("").astype(str).agg(" ".join, axis=1).str.upper()
    group = np.full(breast.shape[0], "Other/unknown", dtype=object)
    group[subtype_text.str.contains("TNBC|TRIPLE|BASAL", regex=True)] = "TNBC/basal"
    group[subtype_text.str.contains("ER\\+|LUMINAL", regex=True)] = "ER+/luminal"
    group[subtype_text.str.contains("HER2", regex=True)] = "HER2+"
    breast["cell_line_group"] = group
    return breast[["ModelID", "cell_line_group"]].drop_duplicates("ModelID")


def add_crc_membership(df: pd.DataFrame, crc: pd.DataFrame) -> pd.DataFrame:
    cols = ["TF", "Luminal_state", "HER2_state", "Basal/TNBC_state", "Luminal_member", "HER2_member", "Basal/TNBC_member"]
    keep = [c for c in cols if c in crc.columns]
    annot = crc[keep].rename(columns={"TF": "gene"}).drop_duplicates("gene") if "TF" in keep else pd.DataFrame()
    return df.merge(annot, on="gene", how="left") if not annot.empty else df


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    crc = load_crc_tfs(args.input_dir)
    crc_tfs = crc["TF"].dropna().astype(str).tolist()
    sig = load_signatures(args.input_dir)
    sig_genes = sorted({gene for genes in sig["signature_gene_list"] for gene in genes})
    all_genes = sorted(set(crc_tfs) | set(sig_genes))

    dep_path = args.input_dir / "depmap_crispr_gene_effect_matrix.csv.gz"
    expr_path = args.input_dir / "depmap_expression_matrix.csv.gz"
    meta = require_table(args.input_dir / "depmap_cell_line_metadata.csv", "depmap_cell_line_metadata")

    dep_cols = pd.read_csv(dep_path, nrows=0).columns.tolist()
    expr_cols = pd.read_csv(expr_path, nrows=0).columns.tolist()
    dep_map = depmap_gene_columns(dep_cols, all_genes)
    expr_map = depmap_gene_columns(expr_cols, sig_genes)

    dep = pd.read_csv(dep_path, usecols=lambda c: c in set(["Unnamed: 0", "ModelID"] + list(dep_map.values())))
    dep_id = "ModelID" if "ModelID" in dep.columns else "Unnamed: 0"
    dep = dep.rename(columns={dep_id: "ModelID", **{v: k for k, v in dep_map.items()}})

    expr = pd.read_csv(expr_path, usecols=lambda c: c in set(["Unnamed: 0", "ModelID"] + list(expr_map.values())))
    expr_id = "ModelID" if "ModelID" in expr.columns else "Unnamed: 0"
    expr = expr.rename(columns={expr_id: "ModelID", **{v: k for k, v in expr_map.items()}})

    breast = load_breast_models(meta)
    dep = breast.merge(dep, on="ModelID", how="inner").drop_duplicates("ModelID")
    expr = breast[["ModelID", "cell_line_group"]].merge(expr, on="ModelID", how="inner").drop_duplicates("ModelID")
    common = sorted(set(dep["ModelID"]).intersection(expr["ModelID"]))
    dep = dep[dep["ModelID"].isin(common)].copy()
    expr = expr[expr["ModelID"].isin(common)].copy()

    dep_genes = [gene for gene in all_genes if gene in dep.columns]
    expr_genes = [gene for gene in sig_genes if gene in expr.columns]
    dep_strength = dep.set_index("ModelID")[dep_genes].apply(pd.to_numeric, errors="coerce") * -1
    expr_matrix = expr.set_index("ModelID")[expr_genes].apply(pd.to_numeric, errors="coerce").T

    z = row_zscore(expr_matrix)
    state_scores = pd.DataFrame(index=z.columns)
    state_meta = sig.set_index("state_id")[["state_label", "subtype_context"]].to_dict("index")
    for _, row in sig.iterrows():
        present = [gene for gene in row["signature_gene_list"] if gene in z.index]
        state_scores[row["state_id"]] = z.loc[present].mean(axis=0) if present else np.nan
    state_scores = state_scores.loc[common]

    state_corr_rows = []
    state_mean_cols = []
    for _, row in sig.iterrows():
        sid = row["state_id"]
        present = [gene for gene in row["signature_gene_list"] if gene in dep.columns]
        mean_col = f"{sid}_mean_gene_effect"
        state_mean_cols.append(mean_col)
        dep[mean_col] = dep[present].mean(axis=1) if present else np.nan
    state_merged = state_scores.reset_index().merge(dep[["ModelID", *state_mean_cols]], on="ModelID", how="inner")
    for _, row in sig.iterrows():
        sid = row["state_id"]
        mean_col = f"{sid}_mean_gene_effect"
        sub = pd.DataFrame(
            {
                "x": pd.to_numeric(state_merged[sid], errors="coerce"),
                "y": pd.to_numeric(state_merged[mean_col], errors="coerce"),
            }
        ).dropna()
        if sub.shape[0] < 8 or sub["x"].nunique() < 4 or sub["y"].nunique() < 4:
            rho, pvalue = np.nan, np.nan
        else:
            rho, pvalue = stats.spearmanr(sub["x"], sub["y"])
        state_corr_rows.append(
            {
                "source_table": "figure8c_state_expression_vs_mean_gene_effect.csv",
                "state_id": sid,
                "state_label": row["state_label"],
                "subtype_context": row["subtype_context"],
                "n_breast_cell_lines": int(sub.shape[0]),
                "spearman_r_state_expr_vs_mean_gene_effect": float(rho) if pd.notna(rho) else np.nan,
                "pvalue": float(pvalue) if pd.notna(pvalue) else np.nan,
                "interpretation": "negative rho means higher state expression associates with stronger CRISPR dependency because lower GeneEffect is more essential",
            }
        )
    state_corr = pd.DataFrame(state_corr_rows)
    state_corr["fdr"] = bh_fdr(state_corr["pvalue"]) if not state_corr.empty else []
    state_corr.to_csv(args.output_dir / "figure8c_state_expression_vs_mean_gene_effect.csv", index=False)

    corr_rows = []
    for state in state_scores.columns:
        for gene in [gene for gene in crc_tfs if gene in dep_strength.columns]:
            tmp = pd.concat([state_scores[state], dep_strength[gene]], axis=1).dropna()
            if tmp.shape[0] < 5:
                continue
            rho, pvalue = stats.spearmanr(tmp.iloc[:, 0], tmp.iloc[:, 1])
            meta_row = state_meta.get(state, {})
            corr_rows.append(
                {
                    "source_table": CORRELATION_SOURCE,
                    "state_id": state,
                    "state_label": meta_row.get("state_label", state),
                    "subtype_context": meta_row.get("subtype_context", np.nan),
                    "n_breast_cell_lines": int(tmp.shape[0]),
                    "spearman_r_state_expr_vs_dependency_strength": float(rho),
                    "pvalue": float(pvalue),
                    "interpretation": "positive rho means higher state expression associates with stronger CRISPR dependency because dependency_strength is -GeneEffect",
                    "gene": gene,
                }
            )
    corr = pd.DataFrame(corr_rows)
    corr["fdr"] = bh_fdr(corr["pvalue"]) if not corr.empty else []

    coverage = pd.DataFrame(
        [
            {
                "source_table": COVERAGE_SOURCE,
                "gene": gene,
                "used_in_crc_tf_heatmap": gene in set(crc_tfs),
                "present_in_depmap_crispr": gene in dep_map,
                "depmap_column": dep_map.get(gene, np.nan),
            }
            for gene in all_genes
        ]
    )

    group_summary = (
        dep[["ModelID", "cell_line_group"]]
        .join(dep_strength[[gene for gene in crc_tfs if gene in dep_strength.columns]], on="ModelID")
        .groupby("cell_line_group", dropna=False)
        .median(numeric_only=True)
        .T
        .reset_index()
        .rename(columns={"index": "gene"})
    )
    for group in ["ER+/luminal", "HER2+", "TNBC/basal", "Other/unknown"]:
        if group not in group_summary.columns:
            group_summary[group] = np.nan
    group_summary.insert(0, "source_table", DEPENDENCY_SOURCE)
    group_summary["row_label"] = group_summary["gene"]
    group_summary = add_crc_membership(group_summary, crc)

    out = pd.concat([corr, state_corr, coverage, group_summary], ignore_index=True, sort=False)
    for col in OUTPUT_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    out = out[OUTPUT_COLUMNS]
    write_table(out, args.output_dir / "figure8c.csv")


if __name__ == "__main__":
    main()
