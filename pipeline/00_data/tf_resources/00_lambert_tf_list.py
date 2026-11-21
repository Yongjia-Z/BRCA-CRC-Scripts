from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TF_resources")
DATA = ROOT / "data"
RESULTS = ROOT / "results"
DATA.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

def find_column(columns: list[str], candidates: list[str]) -> str | None:
    lookup = {col.lower().strip(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    for col in columns:
        lowered = col.lower().strip()
        if any(candidate.lower() in lowered for candidate in candidates):
            return col
    return None


def process_lambert(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    ensembl_col = find_column(raw.columns.tolist(), ["Ensembl ID", "ensembl_gene_id", "ensembl"])
    symbol_col = find_column(raw.columns.tolist(), ["HGNC symbol", "external_gene_name", "Symbol"])
    is_tf_col = find_column(raw.columns.tolist(), ["Is TF?", "is_tf"])
    family_col = find_column(raw.columns.tolist(), ["DBD", "Family", "TF family"])

    if symbol_col is None or is_tf_col is None:
        raise ValueError(f"Cannot identify symbol/is_tf columns in Lambert table: {raw.columns.tolist()[:10]}")

    tf = raw.loc[raw[is_tf_col].astype(str).str.lower().eq("yes")].copy()
    out = pd.DataFrame(
        {
            "ensembl_gene_id": tf[ensembl_col] if ensembl_col else pd.NA,
            "external_gene_name": tf[symbol_col].astype(str),
            "is_tf": tf[is_tf_col].astype(str),
            "tf_category": tf[family_col].astype(str) if family_col else pd.NA,
        }
    )
    out = out.drop_duplicates("external_gene_name").sort_values("external_gene_name")
    return out


def main() -> None:
    tf_list = process_lambert(DATA / "Lambert_human_TF_database.csv")
    tf_list.to_csv(RESULTS / "Lambert_human_TF_list.tsv", sep="\t", index=False)
    summary = pd.DataFrame(
        [
            {
                "resource": "Lambert_human_TF_census",
                "n_TF": int(tf_list.shape[0]),
                "n_unique_symbols": int(tf_list["external_gene_name"].nunique()),
            }
        ]
    )
    summary.to_csv(RESULTS / "TF_resource_summary.tsv", sep="\t", index=False)
    print(json.dumps({"n_lambert_tfs": int(tf_list.shape[0])}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
