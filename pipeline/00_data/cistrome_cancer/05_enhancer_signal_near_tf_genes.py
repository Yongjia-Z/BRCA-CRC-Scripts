from __future__ import annotations

import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "CistromeCancer")
WORKSPACE = DATA_ROOT.parent
DATA = ROOT / "data" / "TCGA_enhancer"
RESULTS = ROOT / "results"
TCGA_RESULTS = DATA_ROOT / "TCGA" / "results"


def parse_position(pos: str) -> tuple[str, int, int]:
    match = re.match(r"(chr[^:]+):(\d+)-(\d+)", str(pos))
    if not match:
        return "", -1, -1
    return match.group(1), int(match.group(2)), int(match.group(3))


def load_positions() -> pd.DataFrame:
    pos = pd.read_csv(DATA / "position.txt", sep="\t", header=None, names=["gene", "position"])
    parsed = pos["position"].apply(parse_position)
    pos["chrom"] = parsed.apply(lambda x: x[0])
    pos["start"] = parsed.apply(lambda x: x[1])
    pos["end"] = parsed.apply(lambda x: x[2])
    pos["tss_window_start"] = (pos["start"] - 100000).clip(lower=0)
    pos["tss_window_end"] = pos["end"] + 100000
    return pos


def summarize_bedgraph(path: Path, genes: pd.DataFrame) -> pd.DataFrame:
    by_chrom = defaultdict(list)
    for row in genes.itertuples(index=False):
        by_chrom[row.chrom].append(row)

    stats = {
        row.gene: {
            "n_segments": 0,
            "sum_signal": 0.0,
            "max_signal": 0.0,
            "mean_signal": 0.0,
            "n_strong_segments_signal_ge10": 0,
        }
        for row in genes.itertuples(index=False)
    }

    with gzip.open(path, "rt") as handle:
        for line in handle:
            chrom, start, end, score = line.rstrip("\n").split("\t")[:4]
            start_i = int(start)
            end_i = int(end)
            score_f = float(score)
            for gene_row in by_chrom.get(chrom, []):
                if end_i < gene_row.tss_window_start or start_i > gene_row.tss_window_end:
                    continue
                rec = stats[gene_row.gene]
                rec["n_segments"] += 1
                rec["sum_signal"] += score_f
                rec["max_signal"] = max(rec["max_signal"], score_f)
                if score_f >= 10:
                    rec["n_strong_segments_signal_ge10"] += 1

    rows = []
    for gene, rec in stats.items():
        rec["mean_signal"] = rec["sum_signal"] / rec["n_segments"] if rec["n_segments"] else 0.0
        rows.append({"TF": gene, **rec})
    return pd.DataFrame(rows)


def main() -> None:
    candidate = pd.read_csv(TCGA_RESULTS / "BRCA_CRC_TF_high_evidence_with_TF_TF.tsv", sep="\t")
    genes = load_positions()
    tf_genes = genes[genes["gene"].isin(candidate["TF"])].copy()
    tf_genes.to_csv(RESULTS / "TCGA_enhancer_candidate_TF_gene_positions.tsv", sep="\t", index=False)

    b1 = summarize_bedgraph(DATA / "BRCA_1_enhancer.bedGraph.gz", tf_genes)
    b2 = summarize_bedgraph(DATA / "BRCA_2_enhancer.bedGraph.gz", tf_genes)
    b1 = b1.add_prefix("BRCA_1_").rename(columns={"BRCA_1_TF": "TF"})
    b2 = b2.add_prefix("BRCA_2_").rename(columns={"BRCA_2_TF": "TF"})
    merged = b1.merge(b2, on="TF", how="outer")
    merged["max_enhancer_signal_100kb"] = merged[["BRCA_1_max_signal", "BRCA_2_max_signal"]].max(axis=1)
    merged["sum_enhancer_signal_100kb"] = merged[["BRCA_1_sum_signal", "BRCA_2_sum_signal"]].sum(axis=1)
    merged["CistromeCancer_enhancer_gene_window_evidence"] = merged["max_enhancer_signal_100kb"] >= 10
    merged = candidate.merge(merged, on="TF", how="left")
    merged.to_csv(RESULTS / "TCGA_enhancer_candidate_CRC_TF_gene_window_signal.tsv", sep="\t", index=False)

    plot_df = merged.sort_values("max_enhancer_signal_100kb", ascending=False).head(30).iloc[::-1]

    summary = {
        "n_candidate_tfs": int(candidate.shape[0]),
        "n_candidate_tfs_with_gene_position": int(tf_genes.shape[0]),
        "n_candidate_tfs_with_enhancer_signal_ge10": int(merged["CistromeCancer_enhancer_gene_window_evidence"].fillna(False).sum()),
    }
    (RESULTS / "TCGA_enhancer_candidate_TF_signal_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "TCGA_enhancer_candidate_TF_signal_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
