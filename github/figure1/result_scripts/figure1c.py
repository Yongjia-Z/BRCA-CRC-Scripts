from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import collapse_pam50_to_context, load_crc_tfs, normalize_pam50, ordered_tfs, read_gene_by_sample_matrix, require_table, row_zscore, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure1c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--clip-z", type=float, default=2.5)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    genes = ordered_tfs(crc)
    expression = read_gene_by_sample_matrix(args.input_dir / "metabric_expression_matrix.csv.gz")
    annotation = require_table(args.input_dir / "metabric_sample_annotation.csv", "metabric_sample_annotation")
    sample_col = "sample_id" if "sample_id" in annotation.columns else annotation.columns[0]
    subtype_col = next((c for c in ["SUBTYPE", "PAM50", "pam50", "claudin_subtype", "subtype"] if c in annotation.columns), None)
    annotation = annotation.rename(columns={sample_col: "sample_id"})
    if subtype_col is not None:
        annotation = annotation.rename(columns={subtype_col: "SUBTYPE"})
        annotation["PAM50_original"] = annotation["SUBTYPE"].map(normalize_pam50)
        annotation["SUBTYPE_CONTEXT"] = annotation["SUBTYPE"].map(collapse_pam50_to_context)
    common_samples = [s for s in expression.columns if s in set(annotation["sample_id"])]
    expression = expression[common_samples]
    present = [g for g in genes if g in expression.index]
    missing = [g for g in genes if g not in expression.index]
    z = row_zscore(expression.loc[present], clip=args.clip_z)
    out = z.reset_index().rename(columns={"gene": "TF"}).melt(id_vars="TF", var_name="sample_id", value_name="z_score")
    out = out.merge(annotation, on="sample_id", how="left")
    if missing:
        out = pd.concat([out, pd.DataFrame({"TF": missing, "record_type": "unavailable", "omission_reason": "not present in external expression matrix"})], ignore_index=True, sort=False)
    write_table(out, args.output_dir / "figure1c.csv")
    (args.output_dir / "figure1c_run_audit.json").write_text(json.dumps({"n_samples": len(common_samples), "n_crc_tfs_present": len(present), "missing_tfs": missing, "subtype_grouping": "raw subtype kept in SUBTYPE; LumA/LumB collapsed to Luminal only in SUBTYPE_CONTEXT", "analysis": "external cohort CRC TF row-wise z-score expression heatmap table"}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
