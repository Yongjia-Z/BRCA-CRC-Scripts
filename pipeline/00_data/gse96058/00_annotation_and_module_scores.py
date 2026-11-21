from __future__ import annotations

import csv
import gzip
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "GSE96058")
WORKSPACE = DATA_ROOT.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
for d in [DATA, RESULTS]:
    d.mkdir(parents=True, exist_ok=True)

MODULES = {
    "luminal_CRC": ["FOXA1", "GATA3", "ESR1", "TFAP2C", "SPDEF", "RARA", "CREB3L4", "PBX1", "AR"],
    "basal_CRC": ["TP63", "SOX9", "RUNX1", "RUNX2", "JUN", "FOS", "CEBPB", "RELA", "FOXC1"],
    "her2_CRC": ["EHF", "EGR3", "STAT1", "TFAP2B", "ELF3", "MYBL2", "NR2F6", "TFAP2C"],
}


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            key = parts[0].replace("!", "")
            values = [v.strip().strip('"') for v in parts[1:]]
            fields.setdefault(key, values)
    if not fields:
        return pd.DataFrame()
    n = max(len(v) for v in fields.values())
    rows = []
    for i in range(n):
        row = {}
        for key, values in fields.items():
            if i < len(values):
                row[key] = values[i]
        rows.append(row)
    ann = pd.DataFrame(rows)
    if "Sample_geo_accession" in ann.columns:
        ann = ann.rename(columns={"Sample_geo_accession": "sample_id"})

    char_cols = [c for c in ann.columns if c == "Sample_characteristics_ch1"]
    if "Sample_characteristics_ch1" in ann.columns:
        pass
    ann.to_csv(RESULTS / "GSE96058_series_matrix_sample_annotation_raw.tsv", sep="\t", index=False)
    return ann


def parse_characteristics(path: Path) -> pd.DataFrame:
    sample_ids = None
    char_rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                sample_ids = [v.strip().strip('"') for v in next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                char_rows.append([v.strip().strip('"') for v in next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]])
    if sample_ids is None:
        return pd.DataFrame()
    rows = [{"sample_id": sid} for sid in sample_ids]
    for values in char_rows:
        for i, value in enumerate(values):
            if i >= len(rows) or ": " not in value:
                continue
            key, val = value.split(": ", 1)
            rows[i][key.strip().replace(" ", "_").replace("-", "_")] = val.strip()
    ann = pd.DataFrame(rows)
    ann.to_csv(RESULTS / "GSE96058_sample_annotation.tsv", sep="\t", index=False)
    return ann


def load_or_convert_sampdesc() -> pd.DataFrame:
    rda = DATA / "GSE96058_sampDesc.rda"
    tsv = RESULTS / "GSE96058_sampDesc.tsv"
    if rda.exists() and not tsv.exists():
        r_code = (
            f"load({json.dumps(str(rda.as_posix()))}); "
            f"x <- GSE96058_sampDesc; "
            f"x$sample_id <- as.character(x$title); "
            f"write.table(x, file={json.dumps(str(tsv.as_posix()))}, sep='\\t', quote=FALSE, row.names=FALSE)"
        )
        try:
            subprocess.run(["Rscript", "-e", r_code], check=True, timeout=120)
        except Exception as exc:
            (RESULTS / "GSE96058_sampDesc_convert_error.txt").write_text(str(exc), encoding="utf-8")
    if tsv.exists():
        ann = pd.read_csv(tsv, sep="\t")
        ann.to_csv(RESULTS / "GSE96058_sample_annotation.tsv", sep="\t", index=False)
        return ann
    return pd.DataFrame()


def extract_crc_expression(path: Path) -> pd.DataFrame:
    genes = sorted({g for module in MODULES.values() for g in module})
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        gene_col = 0
        for row in reader:
            if not row:
                continue
            symbol = row[gene_col].strip().strip('"')
            if symbol in genes:
                vals = [pd.to_numeric(x, errors="coerce") for x in row[1:]]
                rows.append([symbol] + vals)
    if not rows:
        return pd.DataFrame()
    expr = pd.DataFrame(rows, columns=["gene"] + header[1:])
    expr.to_csv(RESULTS / "GSE96058_core_CRC_TF_expression_matrix.tsv", sep="\t", index=False)
    return expr


def score_modules(expr: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    if expr.empty:
        return pd.DataFrame()
    mat = expr.set_index("gene")
    mat = mat.apply(pd.to_numeric, errors="coerce")
    mat = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0)
    scores = pd.DataFrame({"sample_id": mat.columns})
    for module, genes in MODULES.items():
        present = [g for g in genes if g in mat.index]
        scores[module] = mat.loc[present].mean(axis=0).values if present else np.nan
    scores["pan_CRC"] = scores[list(MODULES)].mean(axis=1)
    merged = scores.merge(ann, on="sample_id", how="left")
    merged.to_csv(RESULTS / "GSE96058_CRC_module_scores.tsv", sep="\t", index=False)
    return merged


def find_group_columns(df: pd.DataFrame) -> list[str]:
    candidates = []
    tokens = ["subtype", "pam50", "er_status", "her2_status", "tnbc", "grade", "size", "node", "survival", "rfs"]
    for col in df.columns:
        low = col.lower()
        if any(t in low for t in tokens) and df[col].nunique(dropna=True) >= 2:
            candidates.append(col)
    return candidates[:20]


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in find_group_columns(scores):
        for module in list(MODULES) + ["pan_CRC"]:
            sub = scores[[module, col]].dropna()
            if sub.empty:
                continue
            groups = [g[module].dropna().values for _, g in sub.groupby(col) if len(g[module].dropna()) >= 5]
            if len(groups) >= 2:
                stat, p = stats.kruskal(*groups)
                rows.append({"comparison": col, "module": module, "test": "Kruskal", "statistic": float(stat), "pvalue": float(p), "n": int(sub.shape[0])})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = stats.false_discovery_control(out["pvalue"].fillna(1.0), method="bh")
    out.to_csv(RESULTS / "GSE96058_CRC_module_clinical_association_summary.tsv", sep="\t", index=False)
    return out


def main() -> None:
    series_path = DATA / "GSE96058_series_matrix.txt.gz"
    expr_path = DATA / "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
    series_ok = series_path.exists()
    sampdesc_ok = (DATA / "GSE96058_sampDesc.rda").exists()
    expr_ok = expr_path.exists()

    ann = load_or_convert_sampdesc()
    if ann.empty and series_ok:
        ann = parse_characteristics(series_path)
    expr = extract_crc_expression(expr_path) if expr_ok else pd.DataFrame()
    scores = score_modules(expr, ann) if not expr.empty else pd.DataFrame()
    assoc = summarize(scores) if not scores.empty else pd.DataFrame()

    summary = {
        "series_matrix_found": series_ok,
        "sampdesc_found": sampdesc_ok,
        "expression_found": expr_ok,
        "n_annotation_samples": int(ann.shape[0]) if not ann.empty else 0,
        "n_crc_genes_found": int(expr.shape[0]) if not expr.empty else 0,
        "n_scored_samples": int(scores.shape[0]) if not scores.empty else 0,
        "n_associations": int(assoc.shape[0]) if not assoc.empty else 0,
    }
    (RESULTS / "GSE96058_CRC_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([summary]).to_csv(RESULTS / "GSE96058_CRC_validation_summary.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
