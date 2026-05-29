from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import load_crc_tfs, read_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure5c")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--distal-exclusion-bp", type=int, default=20000)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    crc = load_crc_tfs(args.input_dir)
    contacts = read_table(args.input_dir / "hic_promoter_contacts_long.csv")
    if contacts.empty:
        raise FileNotFoundError("hic_promoter_contacts_long.csv is required")
    contacts = contacts[contacts["TF"].isin(set(crc["TF"]))].copy()
    contacts["distance_to_tss"] = pd.to_numeric(contacts["distance_to_tss"], errors="coerce").abs()
    contacts["contact"] = pd.to_numeric(contacts["contact"], errors="coerce").fillna(0.0)
    contacts["is_distal"] = contacts["distance_to_tss"].ge(args.distal_exclusion_bp)
    group_cols = [c for c in ["TF", "context", "cell_line", "dataset", "matrix_type"] if c in contacts.columns]
    rows = []
    for key, group in contacts.groupby(group_cols, dropna=False):
        key_values = key if isinstance(key, tuple) else (key,)
        rec = dict(zip(group_cols, key_values))
        distal = group[group["is_distal"]]
        rec.update({"n_window_bins": int(group.shape[0]), "n_nonzero_distal_bins": int(distal["contact"].gt(0).sum()), "distal_contact_sum": float(distal["contact"].sum()), "max_distal_contact": float(distal["contact"].max()) if not distal.empty else 0.0})
        rows.append(rec)
    out = pd.DataFrame(rows)
    out["ln_distal_contact_sum"] = np.log1p(out["distal_contact_sum"])
    out["status"] = np.where(out["n_nonzero_distal_bins"].gt(0), "available", "no_distal_contact")
    write_table(out, args.output_dir / "figure5c.csv")


if __name__ == "__main__":
    main()
