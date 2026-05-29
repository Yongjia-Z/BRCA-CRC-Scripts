from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_crc_tfs, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure1d")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--min-layer-score", type=float, default=0.0)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores = require_table(args.input_dir / "figure1a_all_tf_scores.csv", "figure1a_all_tf_scores")
    crc = load_crc_tfs(args.input_dir)
    layers = [c for c in ["RNA", "TF_activity_regulon", "ATAC_motif", "ATAC_footprint", "ATAC_peak_to_gene_accessibility", "context_matched_H3K27ac_SE_MAX", "context_matched_direct_binding_MAX", "scRNA_malignant", "scATAC_cis_accessibility", "DepMap_breast_dependency"] if c in scores.columns]
    data = scores[["TF"] + layers].copy()
    for col in layers:
        data[col] = pd.to_numeric(data[col], errors="coerce").fillna(0).gt(args.min_layer_score)
    groups = data.groupby(layers, dropna=False)["TF"].agg(lambda x: ";".join(sorted(map(str, x)))).reset_index()
    groups["n_TF"] = groups["TF"].str.split(";").map(len)
    crc_set = set(crc["TF"])
    groups["formal_CRC_TFs"] = groups["TF"].map(lambda text: ";".join([g for g in text.split(";") if g in crc_set]))
    groups["n_formal_CRC_TF"] = groups["formal_CRC_TFs"].map(lambda x: 0 if x == "" else len(x.split(";")))
    groups["n_active_layers"] = groups[layers].sum(axis=1)
    groups = groups.sort_values(["n_active_layers", "n_TF"], ascending=[False, False])
    layer_sizes = pd.DataFrame({"evidence_layer": layers, "n_supported_TF_in_universe": [int(data[c].sum()) for c in layers], "n_total_TF_universe": int(data.shape[0])})
    formal_long = data[data["TF"].isin(crc_set)].melt(id_vars="TF", value_vars=layers, var_name="evidence_layer", value_name="is_supported")
    formal_long["is_formal_CRC_TF"] = True
    out = pd.concat([groups, layer_sizes, formal_long], ignore_index=True, sort=False)
    write_table(out, args.output_dir / "figure1d.csv")
    (args.output_dir / "figure1d_run_audit.json").write_text(json.dumps({"analysis": "UpSet-style evidence-layer intersection from full TF universe", "layers": layers, "n_tf_universe": int(data.shape[0])}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
