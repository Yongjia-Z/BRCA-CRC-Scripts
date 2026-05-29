from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import PAM50_CONTEXT_ORDER, collapse_pam50_to_context, load_crc_tfs, normalize_pam50, read_gene_by_sample_matrix, require_table, row_zscore, write_table

CRC_ORDER = ["ESR1", "GATA3", "PGR", "PBX1", "FOXA1", "SPDEF", "ELF3", "EHF", "AR", "TFAP2C", "RUNX1", "TP63", "FOXC1"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure1f")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    genes = [gene for gene in CRC_ORDER if gene in set(crc["TF"].astype(str))]
    activity = read_gene_by_sample_matrix(args.input_dir / "tcga_tf_activity_matrix.csv.gz")
    annotation = require_table(args.input_dir / "tcga_sample_annotation.csv", "tcga_sample_annotation")
    sample_col = "sample_id" if "sample_id" in annotation.columns else annotation.columns[0]
    subtype_col = next((c for c in ["SUBTYPE", "PAM50", "pam50", "subtype"] if c in annotation.columns), None)
    annotation = annotation.rename(columns={sample_col: "sample_id", subtype_col: "SUBTYPE"})
    annotation["PAM50"] = annotation["SUBTYPE"].map(normalize_pam50)
    annotation["SUBTYPE"] = annotation["PAM50"].map(collapse_pam50_to_context)
    common = [s for s in activity.columns if s in set(annotation["sample_id"])]
    activity = row_zscore(activity[common])
    ann = annotation.set_index("sample_id").reindex(common)
    rows = []
    for gene in genes:
        if gene not in activity.index:
            rows.append({"TF": gene, "SUBTYPE": "Unavailable", "mean_TF_activity_zscore": np.nan, "omission_reason": "not present in activity matrix"})
            continue
        for subtype in [st for st in PAM50_CONTEXT_ORDER if ann["SUBTYPE"].eq(st).any()]:
            idx = ann.index[ann["SUBTYPE"].eq(subtype)]
            rows.append({"TF": gene, "SUBTYPE": subtype, "mean_TF_activity_zscore": float(activity.loc[gene, list(idx)].mean())})
    write_table(pd.DataFrame(rows), args.output_dir / "figure1f.csv")


if __name__ == "__main__":
    main()
