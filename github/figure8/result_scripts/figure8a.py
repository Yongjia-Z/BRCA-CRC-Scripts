from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import PAM50_CONTEXT_ORDER, collapse_pam50_to_context, module_score, normalize_pam50, read_gene_by_sample_matrix, require_table, write_table
from scipy import stats


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure8a")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    expr = read_gene_by_sample_matrix(args.input_dir / "bulk_expression_matrix.csv.gz")
    ann = require_table(args.input_dir / "bulk_sample_annotation.csv", "bulk_sample_annotation")
    states = require_table(args.input_dir / "crc_state_signatures.csv", "crc_state_signatures")
    sample_col = "sample_id" if "sample_id" in ann.columns else ann.columns[0]
    subtype_col = next((c for c in ["pam50", "PAM50", "SUBTYPE", "subtype"] if c in ann.columns), None)
    ann = ann.rename(columns={sample_col: "sample_id", subtype_col: "pam50"})
    ann["PAM50"] = ann["pam50"].map(normalize_pam50)
    ann["pam50"] = ann["PAM50"].map(collapse_pam50_to_context)
    common = expr.columns.intersection(ann["sample_id"])
    rows = []
    for state, group in states.groupby("state_id"):
        score = module_score(expr[common], group["gene"].astype(str).tolist(), zscore_genes=True)
        temp = pd.DataFrame({"sample_id": common, "state_id": state, "state_score": score.reindex(common).to_numpy()}).merge(ann, on="sample_id", how="left")
        arrays = [temp.loc[temp["pam50"].eq(st), "state_score"].dropna().to_numpy(float) for st in PAM50_CONTEXT_ORDER if temp["pam50"].eq(st).sum() > 1]
        stat, p = stats.kruskal(*arrays) if len(arrays) >= 2 else (np.nan, np.nan)
        temp["kruskal_statistic"] = stat
        temp["pvalue"] = p
        rows.append(temp)
    write_table(pd.concat(rows, ignore_index=True), args.output_dir / "figure8a.csv")


if __name__ == "__main__":
    main()
