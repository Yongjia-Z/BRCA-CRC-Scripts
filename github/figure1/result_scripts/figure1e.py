from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import PAM50_CONTEXT_ORDER, bh_fdr, collapse_pam50_to_context, load_crc_tfs, normalize_pam50, read_gene_by_sample_matrix, require_table, write_table
from scipy import stats


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure1e")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--top-n", type=int, default=20)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    expression = read_gene_by_sample_matrix(args.input_dir / "tcga_tf_expression_matrix.csv.gz")
    annotation = require_table(args.input_dir / "tcga_sample_annotation.csv", "tcga_sample_annotation")
    sample_col = "sample_id" if "sample_id" in annotation.columns else annotation.columns[0]
    subtype_col = next((c for c in ["SUBTYPE", "PAM50", "pam50", "subtype"] if c in annotation.columns), None)
    if subtype_col is None:
        raise ValueError("sample annotation requires subtype")
    annotation = annotation.rename(columns={sample_col: "sample_id", subtype_col: "SUBTYPE"})
    annotation["PAM50"] = annotation["SUBTYPE"].map(normalize_pam50)
    annotation["SUBTYPE"] = annotation["PAM50"].map(collapse_pam50_to_context)
    common = [s for s in expression.columns if s in set(annotation["sample_id"])]
    expression = expression[common]
    ann = annotation.set_index("sample_id").reindex(common)
    rows = []
    for gene, values in expression.iterrows():
        present_subtypes = [st for st in PAM50_CONTEXT_ORDER if ann["SUBTYPE"].eq(st).any()]
        arrays = [pd.to_numeric(values[ann["SUBTYPE"].eq(st).to_numpy()], errors="coerce").dropna().to_numpy(float) for st in present_subtypes]
        arrays = [a for a in arrays if len(a) > 1]
        if len(arrays) < 2:
            stat = 0.0
            p = 1.0
        else:
            res = stats.kruskal(*arrays)
            stat = float(res.statistic)
            p = float(res.pvalue)
        rec = {"TF": gene, "kruskal_statistic": stat, "pvalue": p, "mean_primary_tpm_log": float(pd.to_numeric(values, errors="coerce").mean()), "primary_detected_fraction_expr_gt1": float(pd.to_numeric(values, errors="coerce").gt(1).mean())}
        for st in PAM50_CONTEXT_ORDER:
            mask = ann["SUBTYPE"].eq(st).to_numpy()
            col_st = st.replace("/", "_").replace(" ", "_")
            rec[f"mean_BRCA_{col_st}"] = float(pd.to_numeric(values[mask], errors="coerce").mean()) if mask.sum() else np.nan
        rows.append(rec)
    out = pd.DataFrame(rows)
    out["subtype_kruskal_FDR"] = bh_fdr(out["pvalue"])
    out["minus_log10_subtype_FDR"] = -np.log10(out["subtype_kruskal_FDR"].clip(lower=1e-300))
    out["is_formal_CRC_TF"] = out["TF"].isin(set(crc["TF"]))
    selected = out.sort_values("subtype_kruskal_FDR").head(args.top_n).sort_values("minus_log10_subtype_FDR", ascending=False).copy()
    write_table(selected, args.output_dir / "figure1e.csv")
    write_table(out, args.output_dir / "figure1e_all_tf_subtype_association.csv")


if __name__ == "__main__":
    main()
