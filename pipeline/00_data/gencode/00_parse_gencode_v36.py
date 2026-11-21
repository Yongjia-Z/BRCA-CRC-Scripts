from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "GENCODE")
DATA = ROOT / "data"
RESULTS = ROOT / "results"
DATA.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

def parse_attr(attr: str, key: str) -> str | None:
    match = re.search(fr'{key} "([^"]+)"', attr)
    return match.group(1) if match else None


def parse_gene_annotation(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attr = parts[8]
            gene_id = parse_attr(attr, "gene_id")
            gene_name = parse_attr(attr, "gene_name")
            gene_type = parse_attr(attr, "gene_type")
            rows.append(
                {
                    "gene_id": gene_id,
                    "ensembl_gene_id": gene_id.split(".")[0] if gene_id else None,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                    "chromosome": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    gtf_path = DATA / "gencode.v36.annotation.gtf.gz"
    anno = parse_gene_annotation(gtf_path)
    anno.to_csv(RESULTS / "gencode_v36_gene_annotation.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "gtf_bytes": gtf_path.stat().st_size,
                "n_genes": int(anno.shape[0]),
                "n_gene_symbols": int(anno["gene_name"].nunique()),
                "n_protein_coding": int((anno["gene_type"] == "protein_coding").sum()),
            }
        ]
    )
    summary.to_csv(RESULTS / "gencode_v36_gene_annotation_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary.iloc[0].to_dict(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
