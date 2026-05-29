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
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure5b")
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
    distal = contacts[contacts["distance_to_tss"].ge(args.distal_exclusion_bp)].copy()
    group_cols = [c for c in ["TF", "context", "cell_line"] if c in distal.columns]
    out = distal.groupby(group_cols, as_index=False).agg(distal_contact_sum=("contact", "sum"), max_distal_contact=("contact", "max"), n_distal_bins=("contact", "size"))
    out["ln_distal_contact_sum"] = np.log1p(out["distal_contact_sum"])
    out["ln_max_distal_contact"] = np.log1p(out["max_distal_contact"])
    write_table(out, args.output_dir / "figure5b.csv")


if __name__ == "__main__":
    main()
