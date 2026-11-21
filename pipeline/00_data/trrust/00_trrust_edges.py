from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = DATA_ROOT.parent
DATASET_DIR = DATA_ROOT / "TRRUST"
DATA_DIR = DATASET_DIR / "data"
RESULTS_DIR = DATASET_DIR / "results"
TCGA_RESULTS_DIR = DATA_ROOT / "TCGA" / "results"

DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def load_trrust(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["regulator_TF", "target_gene", "mode_of_regulation", "pmid"],
        dtype=str,
    )
    df["regulator_TF"] = df["regulator_TF"].str.upper()
    df["target_gene"] = df["target_gene"].str.upper()
    df = df.dropna(subset=["regulator_TF", "target_gene"]).drop_duplicates()
    return df


def integrate_with_crc_edges(trrust: pd.DataFrame) -> dict:
    edge_path = TCGA_RESULTS_DIR / "TF_TF_regulatory_edge_evidence_score_table.tsv"
    if not edge_path.exists():
        return {"has_tcga_edge_table": False}

    edges = pd.read_csv(edge_path, sep="\t")
    edges["regulator_key"] = edges["regulator_TF"].astype(str).str.upper()
    edges["target_key"] = edges["target_TF"].astype(str).str.upper()

    trrust_pairs = trrust[["regulator_TF", "target_gene"]].drop_duplicates().rename(
        columns={"regulator_TF": "regulator_key", "target_gene": "target_key"}
    )
    trrust_pairs["TRRUST_curated_edge_evidence"] = True

    merged = edges.merge(trrust_pairs, on=["regulator_key", "target_key"], how="left")
    merged["TRRUST_curated_edge_evidence"] = merged["TRRUST_curated_edge_evidence"].fillna(False)
    merged = merged.drop(columns=["regulator_key", "target_key"])
    merged["TF_TF_edge_evidence_score_with_TRRUST"] = (
        merged["TF_TF_edge_evidence_score"].astype(int)
        + merged["TRRUST_curated_edge_evidence"].astype(int)
    )
    merged["high_confidence_TF_TF_edge_with_TRRUST"] = (
        merged["TF_TF_edge_evidence_score_with_TRRUST"] >= 5
    )

    overlap = merged[merged["TRRUST_curated_edge_evidence"]].copy()
    merged.to_csv(
        TCGA_RESULTS_DIR / "TF_TF_regulatory_edge_evidence_score_table_with_TRRUST.tsv",
        sep="\t",
        index=False,
    )
    overlap.to_csv(RESULTS_DIR / "TRRUST_overlap_with_CRC_TF_TF_edges.tsv", sep="\t", index=False)

    top = (
        overlap.groupby("regulator_TF")
        .size()
        .reset_index(name="n_TRRUST_supported_CRC_edges")
        .sort_values("n_TRRUST_supported_CRC_edges", ascending=False)
    )
    top.to_csv(RESULTS_DIR / "TRRUST_supported_CRC_edge_regulator_counts.tsv", sep="\t", index=False)

    if not top.empty:
        plot_df = top.head(25).iloc[::-1]

    summary = {
        "has_tcga_edge_table": True,
        "n_crc_edges": int(edges.shape[0]),
        "n_trrust_edges": int(trrust.shape[0]),
        "n_trrust_unique_pairs": int(trrust_pairs.shape[0]),
        "n_crc_edges_supported_by_TRRUST": int(overlap.shape[0]),
        "n_regulators_with_TRRUST_support": int(overlap["regulator_TF"].nunique()) if not overlap.empty else 0,
    }
    return summary


def main() -> None:
    data_path = DATA_DIR / "trrust_rawdata.human.tsv"
    trrust = load_trrust(data_path)

    trrust.to_csv(RESULTS_DIR / "TRRUST_human_TF_target_edges.tsv", sep="\t", index=False)

    summary = {
        "input_file": str(data_path),
        "n_rows": int(trrust.shape[0]),
        "n_unique_regulator_TFs": int(trrust["regulator_TF"].nunique()),
        "n_unique_targets": int(trrust["target_gene"].nunique()),
    }
    integration = integrate_with_crc_edges(trrust)
    summary.update(integration)

    pd.DataFrame([summary]).to_csv(RESULTS_DIR / "TRRUST_human_summary.tsv", sep="\t", index=False)
    (RESULTS_DIR / "TRRUST_human_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
