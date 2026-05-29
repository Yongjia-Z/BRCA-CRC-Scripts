from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import collapse_pam50_to_context, load_crc_tfs, module_score, normalize_pam50, ordered_tfs, read_gene_by_sample_matrix, require_table, row_zscore, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure1b")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--clip-z", type=float, default=2.5)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    genes = ordered_tfs(crc)
    expression = read_gene_by_sample_matrix(args.input_dir / "tcga_expression_matrix.csv.gz")
    annotation = require_table(args.input_dir / "tcga_sample_annotation.csv", "tcga_sample_annotation")
    sample_col = "sample_id" if "sample_id" in annotation.columns else annotation.columns[0]
    subtype_col = next((c for c in ["SUBTYPE", "PAM50", "pam50", "subtype"] if c in annotation.columns), None)
    if subtype_col is None:
        raise ValueError("sample annotation requires PAM50/SUBTYPE column")
    annotation = annotation.rename(columns={sample_col: "sample_id", subtype_col: "SUBTYPE"})
    annotation["PAM50_original"] = annotation["SUBTYPE"].map(normalize_pam50)
    annotation["SUBTYPE_CONTEXT"] = annotation["PAM50_original"].map(collapse_pam50_to_context)
    common_samples = [s for s in expression.columns if s in set(annotation["sample_id"])]
    expression = expression[common_samples]
    present = [g for g in genes if g in expression.index]
    missing = [g for g in genes if g not in expression.index]
    z = row_zscore(expression.loc[present], clip=args.clip_z)
    long = z.reset_index().rename(columns={"gene": "TF"}).melt(id_vars="TF", var_name="sample_id", value_name="z_score")
    ann = annotation.set_index("sample_id").reindex(common_samples).reset_index()
    score = module_score(expression, present, zscore_genes=True).rename("Formal CRC expr").reset_index().rename(columns={"index": "sample_id"})
    out = long.merge(ann, on="sample_id", how="left").merge(score, on="sample_id", how="left")
    if missing:
        out = pd.concat([out, pd.DataFrame({"TF": missing, "sample_id": pd.NA, "z_score": pd.NA, "omission_reason": "not present in expression matrix"})], ignore_index=True, sort=False)
    write_table(out, args.output_dir / "figure1b.csv")
    (args.output_dir / "figure1b_run_audit.json").write_text(json.dumps({"n_samples": len(common_samples), "n_crc_tfs_present": len(present), "missing_tfs": missing, "subtype_grouping": "raw subtype kept in SUBTYPE; collapsed Luminal context stored in SUBTYPE_CONTEXT", "analysis": "row-wise z-score expression heatmap table and CRC TF expression module score"}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
