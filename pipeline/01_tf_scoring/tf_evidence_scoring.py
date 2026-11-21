from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


CONTEXTS = ["Luminal", "HER2", "Basal/TNBC"]
HIGH_CONFIDENCE_CORE_THRESHOLD = 5.50
HIGH_CONFIDENCE_STATE_THRESHOLD = 2.60


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--intermediate-dir", type=Path, default=Path("intermediate") / "figure1a")
    parser.add_argument("--max-workers", type=int, default=10)
    parser.add_argument("--min-expression-log", type=float, default=1.0)
    parser.add_argument("--min-detection-fraction", type=float, default=0.2)
    parser.add_argument("--subtype-fdr", type=float, default=0.05)
    parser.add_argument("--gene-window-bp", type=int, default=100000)
    parser.add_argument("--promoter-window-bp", type=int, default=2000)
    parser.add_argument("--remap-result-dir", type=Path, default=None)
    parser.add_argument("--reuse-intermediate", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rebuild-h3k27ac-from-wig", action="store_true")
    return parser.parse_args()


def find_raw_root(raw_root: Path | None) -> Path:
    if raw_root is not None:
        return raw_root.resolve()
    candidates = []
    here = Path.cwd().resolve()
    for parent in [here, *here.parents]:
        candidates.extend([parent / "data" / "raw", parent / "datasets"])
    for path in candidates:
        if (path / "TCGA").exists() and (path / "GENCODE").exists():
            return path
    raise FileNotFoundError("Cannot locate raw data root. Pass --raw-root explicitly.")


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    suffix = path.name.lower()
    if suffix.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", **kwargs)
    if suffix.endswith(".csv.gz"):
        return pd.read_csv(path, compression="gzip", **kwargs)
    if suffix.endswith(".tsv") or suffix.endswith(".txt"):
        return pd.read_csv(path, sep="\t", **kwargs)
    return pd.read_csv(path, **kwargs)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def first_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError("; ".join(str(p) for p in paths))


def first_existing_optional(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def clean_symbol(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"\s+$", "", regex=True).str.strip()


def normalize_context(value: object) -> str:
    text = str(value).strip().lower().replace("_", " ").replace("-", " ")
    if text in {"brca luma", "luma", "lum a", "luminal a", "brca lumb", "lumb", "lum b", "luminal b", "er+", "luminal"}:
        return "Luminal"
    if "her2" in text:
        return "HER2"
    if "basal" in text or "tnbc" in text or "triple negative" in text:
        return "Basal/TNBC"
    if "normal" in text:
        return "Normal"
    return str(value).strip()


def sample_type_from_barcode(sample_id: str) -> str:
    code = sample_id[13:15] if len(sample_id) >= 15 else ""
    return {"01": "Primary Tumor", "02": "Recurrent Tumor", "06": "Metastatic", "11": "Solid Tissue Normal"}.get(code, f"Other_{code}")


def bh_fdr(values: pd.Series | list[float]) -> pd.Series:
    s = pd.Series(values, dtype=float)
    valid = np.isfinite(s.to_numpy())
    out = np.ones(len(s), dtype=float)
    if valid.any():
        out[valid] = multipletests(s.to_numpy()[valid], method="fdr_bh")[1]
    return pd.Series(out, index=s.index)


def rank01(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").rank(method="average", pct=True).fillna(0.0).clip(0.0, 1.0)


def scale01(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float)
    lo = float(x.min()) if len(x) else 0.0
    hi = float(x.max()) if len(x) else 0.0
    if hi <= lo:
        return pd.Series(np.where(x.gt(0), 1.0, 0.0), index=x.index)
    return ((x - lo) / (hi - lo)).clip(0.0, 1.0)


def log_max_scale01(values: pd.Series) -> pd.Series:
    x = np.log1p(pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float))
    hi = float(x.max()) if len(x) else 0.0
    if hi <= 0:
        return pd.Series(0.0, index=x.index)
    return (x / hi).clip(0.0, 1.0)


def log_scale01(values: pd.Series) -> pd.Series:
    return scale01(np.log1p(pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float)))


def context_state_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    out = pd.DataFrame({"TF": df["TF"]})
    for ctx in CONTEXTS:
        col = f"{prefix}_{ctx}_score"
        out[col] = df[col] if col in df.columns else 0.0
    return out


def build_tf_universe(raw_root: Path, out_dir: Path, reuse: bool) -> pd.DataFrame:
    out = out_dir / "figure1a_tf_universe.csv"
    if reuse and out.exists():
        return pd.read_csv(out)
    path = first_existing([
        raw_root / "TF_resources" / "results" / "Lambert_human_TF_list.tsv",
        raw_root / "TF_resources" / "Lambert_human_TF_list.tsv",
    ])
    df = pd.read_csv(path, sep="\t")
    col = "external_gene_name" if "external_gene_name" in df.columns else "TF"
    tf = pd.DataFrame({"TF": clean_symbol(df[col])})
    tf = tf[tf["TF"].ne("") & tf["TF"].ne("nan")].drop_duplicates().sort_values("TF").reset_index(drop=True)
    write_csv(tf, out)
    return tf


def build_gene_annotation(raw_root: Path, tf_universe: pd.DataFrame, out_dir: Path, reuse: bool) -> pd.DataFrame:
    out = out_dir / "figure1a_gene_annotation.csv"
    if reuse and out.exists():
        return pd.read_csv(out)
    path = first_existing([
        raw_root / "GENCODE" / "results" / "gencode_v36_gene_annotation.tsv",
        raw_root / "GENCODE" / "gencode_v36_gene_annotation.tsv",
    ])
    anno = pd.read_csv(path, sep="\t")
    anno = anno.rename(columns={"gene_name": "TF", "chromosome": "chrom"})
    keep = anno[anno["TF"].isin(set(tf_universe["TF"]))].copy()
    keep = keep[keep.get("gene_type", "").astype(str).str.contains("protein_coding|transcription", regex=True, na=False) | keep["TF"].isin(tf_universe["TF"])]
    keep["start"] = pd.to_numeric(keep["start"], errors="coerce")
    keep["end"] = pd.to_numeric(keep["end"], errors="coerce")
    keep = keep.dropna(subset=["start", "end"]).copy()
    keep["start"] = keep["start"].astype(int)
    keep["end"] = keep["end"].astype(int)
    keep["tss"] = np.where(keep["strand"].astype(str).eq("-"), keep["end"], keep["start"])
    keep = keep[["TF", "ensembl_gene_id", "chrom", "start", "end", "strand", "tss", "gene_type"]].drop_duplicates("TF")
    write_csv(keep.sort_values("TF"), out)
    return keep


def load_tcga_expression_matrix(raw_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm_path = first_existing([
        raw_root / "TCGA" / "data" / "expression" / "TCGA-BRCA.star_tpm.tsv.gz",
        raw_root / "TCGA" / "TCGA-BRCA.star_tpm.tsv.gz",
    ])
    anno_path = first_existing([
        raw_root / "GENCODE" / "results" / "gencode_v36_gene_annotation.tsv",
        raw_root / "GENCODE" / "gencode_v36_gene_annotation.tsv",
    ])
    tpm = pd.read_csv(tpm_path, sep="\t")
    tpm["ensembl_gene_id"] = tpm["Ensembl_ID"].astype(str).str.replace(r"\..*$", "", regex=True)
    anno = pd.read_csv(anno_path, sep="\t")[["ensembl_gene_id", "gene_name"]].drop_duplicates("ensembl_gene_id")
    merged = tpm.merge(anno, on="ensembl_gene_id", how="left").dropna(subset=["gene_name"])
    sample_cols = [c for c in tpm.columns if c not in {"Ensembl_ID", "ensembl_gene_id"}]
    expr = merged.groupby("gene_name", as_index=True)[sample_cols].mean()
    samples = pd.DataFrame({"sample_id": sample_cols, "patient_id": [s[:12] for s in sample_cols], "sample_type": [sample_type_from_barcode(s) for s in sample_cols]})
    clinical_paths = [
        raw_root / "TCGA" / "data" / "cbioportal_brca_patient_clinical_wide.tsv",
        raw_root / "TCGA" / "results" / "cbioportal_brca_patient_clinical_wide.tsv",
    ]
    clinical = next((p for p in clinical_paths if p.exists()), None)
    if clinical is not None:
        cli = pd.read_csv(clinical, sep="\t")
        id_col = "patientId" if "patientId" in cli.columns else ("PATIENT_ID" if "PATIENT_ID" in cli.columns else None)
        subtype_col = "SUBTYPE" if "SUBTYPE" in cli.columns else None
        if id_col and subtype_col:
            samples = samples.merge(cli[[id_col, subtype_col]].rename(columns={id_col: "patient_id"}), on="patient_id", how="left")
    return expr, samples


def kruskal_by_tf(matrix: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    ann = samples[samples["sample_type"].eq("Primary Tumor") & samples["sample_id"].isin(matrix.columns)].copy()
    subtype_series = ann["SUBTYPE"] if "SUBTYPE" in ann.columns else pd.Series("Unknown", index=ann.index)
    ann["subtype_raw"] = subtype_series.astype(str)
    ann["context"] = ann["subtype_raw"].map(normalize_context)
    ann = ann[ann["context"].isin(CONTEXTS)]
    rows = []
    for tf, vals in matrix[ann["sample_id"].tolist()].iterrows():
        groups = []
        means = {}
        for subtype, sub in ann.groupby("subtype_raw", sort=False):
            if normalize_context(subtype) not in CONTEXTS:
                continue
            x = vals[sub["sample_id"].tolist()].astype(float).dropna()
            if len(x) >= 3 and float(x.std()) > 0:
                groups.append(x.to_numpy(float))
        for ctx, sub in ann.groupby("context", sort=False):
            x = vals[sub["sample_id"].tolist()].astype(float).dropna()
            means[f"mean_{ctx}"] = float(x.mean()) if len(x) else np.nan
        p = stats.kruskal(*groups).pvalue if len(groups) >= 2 else 1.0
        rows.append({"TF": tf, "subtype_kruskal_p": p, **means})
    out = pd.DataFrame(rows)
    out["subtype_kruskal_FDR"] = bh_fdr(out["subtype_kruskal_p"]).to_numpy()
    return out


def build_expression_evidence(raw_root: Path, tf_universe: pd.DataFrame, out_dir: Path, reuse: bool, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = out_dir / "figure1a_tcga_expression_evidence.csv"
    state_out = out_dir / "figure1a_tcga_expression_context_scores.csv"
    matrix_out = out_dir / "figure1a_tcga_tf_expression_matrix.csv.gz"
    if reuse and out.exists() and state_out.exists() and matrix_out.exists():
        return pd.read_csv(out), pd.read_csv(state_out), pd.read_csv(matrix_out, compression="gzip", index_col=0)
    expr, samples = load_tcga_expression_matrix(raw_root)
    tf_expr = expr.loc[[g for g in tf_universe["TF"] if g in expr.index]].copy()
    primary = samples.loc[samples["sample_type"].eq("Primary Tumor"), "sample_id"].tolist()
    primary = [s for s in primary if s in tf_expr.columns]
    subtype_summary = kruskal_by_tf(tf_expr, samples)
    rows = []
    for tf, vals in tf_expr.iterrows():
        x = vals[primary].astype(float)
        rows.append({"TF": tf, "mean_primary_tpm_log": float(x.mean()), "detected_fraction_expr_gt1": float((x > args.min_expression_log).mean())})
    base = pd.DataFrame(rows).merge(subtype_summary, on="TF", how="left")
    for ctx in CONTEXTS:
        if f"mean_{ctx}" not in base.columns:
            base[f"mean_{ctx}"] = 0.0
    base["RNA"] = (
        base["mean_primary_tpm_log"].gt(args.min_expression_log)
        & base["detected_fraction_expr_gt1"].gt(args.min_detection_fraction)
        & base["subtype_kruskal_FDR"].fillna(1.0).le(args.subtype_fdr)
    ).astype(float)
    base["expression_dominant_module"] = base[[f"mean_{c}" for c in CONTEXTS]].idxmax(axis=1).str.replace("mean_", "", regex=False)
    base["expression_dominant_module"] = base["expression_dominant_module"].map({"Luminal": "luminal", "HER2": "HER2/epithelial", "Basal/TNBC": "basal/TNBC"}).fillna("supporting")
    state = pd.DataFrame({"TF": base["TF"]})
    for ctx in CONTEXTS:
        state[f"expression_{ctx}_score"] = rank01(base[f"mean_{ctx}"]) * base["detected_fraction_expr_gt1"].gt(args.min_detection_fraction).astype(float)
    write_csv(base, out)
    write_csv(state, state_out)
    tf_expr.to_csv(matrix_out, compression="gzip")
    return base, state, tf_expr


def build_activity_evidence(raw_root: Path, tf_universe: pd.DataFrame, out_dir: Path, reuse: bool, args: argparse.Namespace, tf_expr: pd.DataFrame | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = out_dir / "figure1a_dorothea_activity_evidence.csv"
    state_out = out_dir / "figure1a_dorothea_activity_context_scores.csv"
    if reuse and out.exists() and state_out.exists():
        return pd.read_csv(out), pd.read_csv(state_out)
    if tf_expr is None:
        expr, samples = load_tcga_expression_matrix(raw_root)
    else:
        expr, samples = load_tcga_expression_matrix(raw_root)
    dorothea = pd.read_csv(first_existing([
        raw_root / "TF_resources" / "data" / "dorothea_hs_pancancer.tsv",
        raw_root / "TF_resources" / "dorothea_hs_pancancer.tsv",
    ]), sep="\t")
    primary = samples.loc[samples["sample_type"].eq("Primary Tumor"), "sample_id"].tolist()
    primary = [s for s in primary if s in expr.columns]
    expr_primary = expr[primary].astype(float)
    z = expr_primary.sub(expr_primary.mean(axis=1), axis=0).div(expr_primary.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    dorothea = dorothea[dorothea["target"].isin(z.index) & dorothea["tf"].isin(set(tf_universe["TF"]))].copy()
    activities = {}
    regulon_rows = []
    for tf, group in dorothea.groupby("tf", sort=False):
        group = group.drop_duplicates("target")
        targets = [t for t in group["target"].astype(str) if t in z.index]
        if len(targets) < 10:
            continue
        mor = group.set_index("target").loc[targets, "mor"].astype(float)
        activities[tf] = z.loc[targets].mul(mor, axis=0).mean(axis=0)
        regulon_rows.append({"TF": tf, "n_targets_used": len(targets)})
    activity = pd.DataFrame(activities).T
    summary = kruskal_by_tf(activity, samples) if not activity.empty else pd.DataFrame(columns=["TF", "subtype_kruskal_FDR"])
    summary = summary.rename(columns={"subtype_kruskal_p": "activity_subtype_kruskal_p", "subtype_kruskal_FDR": "activity_subtype_kruskal_FDR", "mean_Luminal": "mean_activity_Luminal", "mean_HER2": "mean_activity_HER2", "mean_Basal/TNBC": "mean_activity_Basal/TNBC"})
    regulon = pd.DataFrame(regulon_rows)
    base = tf_universe[["TF"]].merge(summary, on="TF", how="left")
    for ctx in CONTEXTS:
        col = f"mean_activity_{ctx}"
        if col not in base.columns:
            base[col] = 0.0
    base["activity_subtype_evidence"] = base["activity_subtype_kruskal_FDR"].fillna(1.0).le(args.subtype_fdr).astype(float)
    base["regulon_present"] = base["TF"].isin(set(regulon["TF"]) if not regulon.empty else set()).astype(float)
    base["TF_activity_regulon"] = (base["activity_subtype_evidence"] + base["regulon_present"]) / 2.0
    state = pd.DataFrame({"TF": base["TF"]})
    for ctx in CONTEXTS:
        state[f"activity_{ctx}_score"] = rank01(base[f"mean_activity_{ctx}"])
    write_csv(base, out)
    write_csv(state, state_out)
    return base, state


def read_excel_data(path: Path, sheet: str, header_row_1based: int) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet, skiprows=header_row_1based - 1)


def build_atac_evidence(raw_root: Path, tf_universe: pd.DataFrame, out_dir: Path, reuse: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = out_dir / "figure1a_tcga_atac_evidence.csv"
    state_out = out_dir / "figure1a_tcga_atac_context_scores.csv"
    if reuse and out.exists() and state_out.exists():
        return pd.read_csv(out), pd.read_csv(state_out)
    atac_dir = raw_root / "TCGA" / "data" / "ATAC"
    motif = read_excel_data(atac_dir / "TCGA_ATAC_Data_S4_motif_enrichment.xlsx", "Motif_Enrichment", 15)
    cluster_cols = [c for c in motif.columns if str(c).startswith("Cluster_")]
    for col in cluster_cols:
        motif[col] = pd.to_numeric(motif[col], errors="coerce")
    motif["max_cluster_motif_enrichment"] = motif[cluster_cols].max(axis=1)
    motif["best_motif_cluster"] = motif[cluster_cols].idxmax(axis=1)
    motif_best = motif.sort_values("max_cluster_motif_enrichment", ascending=False).drop_duplicates("TF_Name").rename(columns={"TF_Name": "TF"})
    motif_hits = set(motif_best.loc[motif_best["max_cluster_motif_enrichment"].fillna(0).ge(10), "TF"].astype(str))
    motif_score = motif_best.set_index("TF")["max_cluster_motif_enrichment"]
    corr = read_excel_data(atac_dir / "TCGA_ATAC_Data_S6_footprinting.xlsx", "Correlation Depth Flank to RNA", 18)
    corr["footprint_supported"] = pd.to_numeric(corr["depth_FDR"], errors="coerce").le(0.05) | pd.to_numeric(corr["flank_FDR"], errors="coerce").le(0.05)
    footprint = motif_best[["TF", "TF_Motif_Name"]].merge(corr[["TF_Motif_Name", "footprint_supported"]], on="TF_Motif_Name", how="left")
    links = read_excel_data(atac_dir / "TCGA_ATAC_Data_S7_peak_to_gene_links.xlsx", "BRCA_Links", 29)
    links["Linked_Gene"] = clean_symbol(links["Linked_Gene"])
    link_counts = links.groupby("Linked_Gene").size().rename("n_peak_links_BRCA")
    nearest_counts = pd.Series(dtype=float)
    peak_anno_path = raw_root / "TCGA" / "results" / "ATAC_BRCA_peak_annotation.tsv"
    if peak_anno_path.exists():
        peak_anno = pd.read_csv(peak_anno_path, sep="\t")
        if "Gene Symbol" in peak_anno.columns and "Peak_Name" in peak_anno.columns:
            nearest_counts = peak_anno.groupby("Gene Symbol")["Peak_Name"].count()
    subtype_path = first_existing_optional([
        raw_root / "TCGA" / "results" / "ATAC_subtype_specific_peak_to_TF_universe_links.tsv",
        raw_root / "TCGA" / "results" / "ATAC_subtype_specific_peak_to_candidate_TF_links.tsv",
    ])
    subtype_counts = pd.DataFrame(index=tf_universe["TF"])
    if subtype_path is not None:
        sub = pd.read_csv(subtype_path, sep="\t")
        sub["context"] = sub["best_subtype"].map(normalize_context)
        pivot = sub[sub["context"].isin(CONTEXTS)].pivot_table(index="Linked_Gene", columns="context", values="Peak_Name", aggfunc="count", fill_value=0)
        subtype_counts = pivot.reindex(tf_universe["TF"]).fillna(0.0)
    base = tf_universe[["TF"]].copy()
    base["ATAC_motif"] = base["TF"].isin(motif_hits).astype(float)
    base["ATAC_motif_score_raw"] = base["TF"].map(motif_score).fillna(0.0)
    base = base.merge(footprint, on="TF", how="left")
    base["ATAC_footprint"] = base["footprint_supported"].fillna(False).astype(float)
    base["n_peak_links_BRCA"] = base["TF"].map(link_counts).fillna(0.0)
    base["n_nearest_TF_peaks"] = base["TF"].map(nearest_counts).fillna(0.0)
    base["ATAC_peak_to_gene_accessibility"] = (base["n_peak_links_BRCA"].gt(0).astype(float) + base["n_nearest_TF_peaks"].gt(0).astype(float)) / 2.0
    for ctx in CONTEXTS:
        base[f"ATAC_peak_{ctx}_raw"] = subtype_counts[ctx].reindex(base["TF"]).to_numpy() if ctx in subtype_counts.columns else 0.0
    audited_tcga = raw_root / "TCGA" / "results" / "BRCA_CRC_TF_with_enhancer_H3K27ac_evidence.tsv"
    if audited_tcga.exists():
        audited = pd.read_csv(audited_tcga, sep="\t").drop_duplicates("TF").set_index("TF")
        for source_col, target_col in [
            ("ATAC_motif_enrichment_evidence", "ATAC_motif"),
            ("ATAC_footprint_expression_evidence", "ATAC_footprint"),
        ]:
            if source_col in audited.columns:
                base[target_col] = base["TF"].map(audited[source_col].astype(str).str.lower().isin(["true", "1", "yes"]).astype(float)).fillna(base[target_col])
        if {"ATAC_peak_to_TF_gene_evidence", "ATAC_nearest_TF_peak_evidence"}.issubset(audited.columns):
            linked = base["TF"].map(audited["ATAC_peak_to_TF_gene_evidence"].astype(str).str.lower().isin(["true", "1", "yes"]).astype(float)).fillna(base["n_peak_links_BRCA"].gt(0).astype(float))
            nearest = base["TF"].map(audited["ATAC_nearest_TF_peak_evidence"].astype(str).str.lower().isin(["true", "1", "yes"]).astype(float)).fillna(base["n_nearest_TF_peaks"].gt(0).astype(float))
            base["ATAC_peak_to_gene_accessibility"] = (linked + nearest) / 2.0
        for ctx, col in [("Luminal", ["n_peak_links_BRCA_LumA", "n_peak_links_BRCA_LumB"]), ("HER2", ["n_peak_links_BRCA_Her2"]), ("Basal/TNBC", ["n_peak_links_BRCA_Basal"])]:
            available = [c for c in col if c in audited.columns]
            if available:
                values = audited[available].max(axis=1)
                base[f"ATAC_peak_{ctx}_raw"] = base["TF"].map(values).fillna(base[f"ATAC_peak_{ctx}_raw"])
    state = pd.DataFrame({"TF": base["TF"]})
    for ctx in CONTEXTS:
        state[f"atac_peak_{ctx}_score"] = log_max_scale01(base[f"ATAC_peak_{ctx}_raw"])
    write_csv(base, out)
    write_csv(state, state_out)
    return base, state


def read_bed(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.name.endswith(".gz") else open
    rows = []
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                rows.append((parts[0], int(float(parts[1])), int(float(parts[2]))))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def count_region_overlaps(peaks: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    if peaks.empty or regions.empty:
        return pd.DataFrame(columns=["TF", "n_peak_overlaps", "max_overlap_bp"])
    out = []
    peaks = peaks[["chrom", "start", "end"]].copy()
    regions = regions[["TF", "chrom", "region_start", "region_end"]].copy()
    for chrom, sub_regions in regions.groupby("chrom", sort=False):
        sub_peaks = peaks[peaks["chrom"].eq(chrom)].sort_values("start")
        if sub_peaks.empty:
            continue
        starts = sub_peaks["start"].to_numpy(np.int64)
        ends = sub_peaks["end"].to_numpy(np.int64)
        limits = np.searchsorted(starts, sub_regions["region_end"].to_numpy(np.int64), side="left")
        for tf, left, right, limit in zip(sub_regions["TF"], sub_regions["region_start"], sub_regions["region_end"], limits):
            if limit <= 0:
                continue
            idx = np.nonzero(ends[:limit] > int(left))[0]
            if idx.size == 0:
                continue
            overlaps = np.minimum(ends[idx], int(right)) - np.maximum(starts[idx], int(left))
            out.append({"TF": tf, "n_peak_overlaps": int(idx.size), "max_overlap_bp": int(overlaps.max())})
    return pd.DataFrame(out)


def make_regions(gene_annotation: pd.DataFrame, flank: int, promoter: bool = False) -> pd.DataFrame:
    g = gene_annotation.copy()
    if promoter:
        g["region_start"] = (g["tss"] - flank).clip(lower=0).astype(int)
        g["region_end"] = (g["tss"] + flank).astype(int)
    else:
        g["region_start"] = (g["start"] - flank).clip(lower=0).astype(int)
        g["region_end"] = (g["end"] + flank).astype(int)
    return g[["TF", "chrom", "region_start", "region_end"]]


def build_h3k27ac_evidence(raw_root: Path, tf_universe: pd.DataFrame, gene_annotation: pd.DataFrame, out_dir: Path, reuse: bool, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = out_dir / "figure1a_h3k27ac_se_proxy_evidence.csv"
    state_out = out_dir / "figure1a_h3k27ac_context_scores.csv"
    overlap_out = out_dir / "figure1a_h3k27ac_peak_overlap_long.csv"
    if reuse and out.exists() and state_out.exists():
        return pd.read_csv(out), pd.read_csv(state_out)
    base = tf_universe[["TF"]].copy()
    h3_inputs = pd.DataFrame({"TF": base["TF"]}).set_index("TF")
    summary = raw_root / "GSE85158_H3K27ac" / "results" / "GSE85158_TF_window_H3K27ac_subtype_summary_hg19_to_hg38.tsv"
    if summary.exists():
        h3 = pd.read_csv(summary, sep="\t").drop_duplicates("TF").set_index("TF")
        for col in [
            "GSE85158_Luminal_minus_normal_mean_signal",
            "GSE85158_HER2_minus_normal_mean_signal",
            "GSE85158_Basal_TNBC_minus_normal_mean_signal",
        ]:
            if col in h3.columns:
                h3_inputs[col] = h3[col].reindex(h3_inputs.index).fillna(0.0)
    encode_h3 = first_existing_optional([
        raw_root / "ENCODE" / "results" / "ENCODE_MCF7_H3K27ac_TF_universe_gene_window_summary.tsv",
        raw_root / "ENCODE" / "results" / "ENCODE_MCF7_H3K27ac_candidate_TF_gene_window_summary.tsv",
    ])
    if encode_h3 is not None:
        enc = pd.read_csv(encode_h3, sep="\t").drop_duplicates("TF").set_index("TF")
        for col in ["ENCODE_MCF7_H3K27ac_max_peaks_100kb", "ENCODE_MCF7_H3K27ac_max_signal"]:
            if col in enc.columns:
                h3_inputs[col] = enc[col].reindex(h3_inputs.index).fillna(0.0)
    rose = first_existing_optional([
        raw_root / "ENCODE" / "results" / "ENCODE_MCF7_ROSE_like_SE_TF_universe_overlap.tsv",
        raw_root / "ENCODE" / "results" / "ENCODE_MCF7_ROSE_like_SE_candidate_TF_overlap.tsv",
    ])
    if rose is not None:
        ro = pd.read_csv(rose, sep="\t").drop_duplicates("TF").set_index("TF")
        for col in ["n_ROSE_like_SE_100kb", "max_ROSE_like_SE_signal"]:
            if col in ro.columns:
                h3_inputs[col] = ro[col].reindex(h3_inputs.index).fillna(0.0)
    dbsuper = raw_root / "SEdb" / "results" / "dbSUPER_breast_CRC_TF_100kb_overlap_all.tsv"
    if dbsuper.exists():
        db = pd.read_csv(dbsuper, sep="\t")
        if {"TF", "cell_line", "n_dbSUPER_SE_100kb"}.issubset(db.columns):
            h3_inputs["dbSUPER_MCF7_SE_count"] = db[db["cell_line"].eq("MCF-7")].groupby("TF")["n_dbSUPER_SE_100kb"].max().reindex(h3_inputs.index).fillna(0.0)
            h3_inputs["dbSUPER_HCC1954_SE_count"] = db[db["cell_line"].eq("HCC1954")].groupby("TF")["n_dbSUPER_SE_100kb"].max().reindex(h3_inputs.index).fillna(0.0)
    peak_files = []
    her2_h3_dir = raw_root / "GSE251868_HER2_H3K27ac" / "data"
    if her2_h3_dir.exists():
        for path in her2_h3_dir.glob("*.bed.gz"):
            peak_files.append(("HER2", path))
    tnbc_h3_dir = raw_root / "GSE245830_TNBC_H3K27ac" / "data"
    if tnbc_h3_dir.exists():
        for path in tnbc_h3_dir.glob("*.bed*"):
            peak_files.append(("Basal/TNBC", path))
    h3_inputs["HER2_H3K27ac_peak_count"] = 0.0
    h3_inputs["Basal_TNBC_H3K27ac_peak_count"] = 0.0
    if peak_files:
        regions = make_regions(gene_annotation, args.gene_window_bp, promoter=False)
        long = []
        for ctx, path in peak_files:
            overlaps = count_region_overlaps(read_bed(path), regions)
            if not overlaps.empty:
                overlaps["context"] = ctx
                overlaps["source_file"] = path.name
                long.append(overlaps)
        if long:
            ol = pd.concat(long, ignore_index=True)
            write_csv(ol, overlap_out)
            h3_inputs["HER2_H3K27ac_peak_count"] = ol[ol["context"].eq("HER2")].groupby("TF")["n_peak_overlaps"].max().reindex(h3_inputs.index).fillna(0.0)
            h3_inputs["Basal_TNBC_H3K27ac_peak_count"] = ol[ol["context"].eq("Basal/TNBC")].groupby("TF")["n_peak_overlaps"].max().reindex(h3_inputs.index).fillna(0.0)
    def source(name: str) -> pd.Series:
        return h3_inputs[name] if name in h3_inputs.columns else pd.Series(0.0, index=h3_inputs.index)
    base["Luminal_H3K27ac_SE"] = np.maximum.reduce([
        log_scale01(source("ENCODE_MCF7_H3K27ac_max_peaks_100kb")),
        log_scale01(source("n_ROSE_like_SE_100kb")),
        log_scale01(source("dbSUPER_MCF7_SE_count")),
        log_scale01(source("GSE85158_Luminal_minus_normal_mean_signal")),
    ])
    base["HER2_H3K27ac_SE"] = np.maximum.reduce([
        log_scale01(source("HER2_H3K27ac_peak_count")),
        log_scale01(source("dbSUPER_HCC1954_SE_count")),
        log_scale01(source("GSE85158_HER2_minus_normal_mean_signal")),
    ])
    base["Basal_TNBC_H3K27ac_SE"] = np.maximum.reduce([
        log_scale01(source("Basal_TNBC_H3K27ac_peak_count")),
        log_scale01(source("GSE85158_Basal_TNBC_minus_normal_mean_signal")),
    ])
    base["context_matched_H3K27ac_SE_MAX"] = base[["Luminal_H3K27ac_SE", "HER2_H3K27ac_SE", "Basal_TNBC_H3K27ac_SE"]].max(axis=1)
    state = pd.DataFrame({"TF": base["TF"], "h3_Luminal_score": base["Luminal_H3K27ac_SE"], "h3_HER2_score": base["HER2_H3K27ac_SE"], "h3_Basal/TNBC_score": base["Basal_TNBC_H3K27ac_SE"]})
    write_csv(base, out)
    write_csv(state, state_out)
    return base, state


def direct_binding_jobs(raw_root: Path) -> pd.DataFrame:
    rows = []
    tnbc_dir = raw_root / "ChIPAtlas_TNBC_direct_binding" / "data"
    if tnbc_dir.exists():
        for path in tnbc_dir.glob("*.bed*"):
            parts = path.name.split(".")
            tf = parts[1] if len(parts) > 2 else ""
            rows.append({"regulator_TF": tf, "context": "Basal/TNBC", "file": path})
    skip = {"GSE", "GSM", "SRX", "TNBC", "CHIP", "CHIPSEQ", "CUT", "TAG", "SEQ", "DATA", "PEAKS", "NARROWPEAK", "BED", "FILT", "HG19", "HG38", "BT549", "HS578T", "MDA231", "MDA468"}
    for dataset_dir in raw_root.glob("GSE*_*_TNBC*"):
        if not dataset_dir.is_dir():
            continue
        tokens = [tok.upper() for tok in re.split(r"[^A-Za-z0-9]+", dataset_dir.name) if tok]
        tf_candidates = [tok for tok in tokens if tok not in skip and not re.match(r"^(GSE|GSM|SRX)\d+$", tok)]
        if not tf_candidates:
            continue
        tf = tf_candidates[0]
        for path in dataset_dir.glob("*.narrowPeak.gz"):
            rows.append({"regulator_TF": tf, "context": "Basal/TNBC", "file": path})
    return pd.DataFrame(rows).drop_duplicates("file")


def one_direct_overlap(job: dict, promoters: pd.DataFrame) -> pd.DataFrame:
    peaks = read_bed(Path(job["file"]))
    overlaps = count_region_overlaps(peaks, promoters)
    if overlaps.empty:
        return pd.DataFrame(columns=["regulator_TF", "target_TF", "context", "source_file", "n_peak_overlaps", "max_overlap_bp"])
    overlaps = overlaps.rename(columns={"TF": "target_TF"})
    overlaps["regulator_TF"] = str(job["regulator_TF"]).strip()
    overlaps["context"] = job["context"]
    overlaps["source_file"] = Path(job["file"]).name
    return overlaps[["regulator_TF", "target_TF", "context", "source_file", "n_peak_overlaps", "max_overlap_bp"]]


def build_direct_binding_evidence(raw_root: Path, tf_universe: pd.DataFrame, gene_annotation: pd.DataFrame, out_dir: Path, reuse: bool, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = out_dir / "figure1a_direct_binding_evidence.csv"
    state_out = out_dir / "figure1a_direct_binding_context_scores.csv"
    long_out = out_dir / "figure1a_direct_binding_promoter_overlap_long.csv"
    if reuse and out.exists() and state_out.exists():
        return pd.read_csv(out), pd.read_csv(state_out)
    promoters = make_regions(gene_annotation, args.promoter_window_bp, promoter=True)
    jobs = direct_binding_jobs(raw_root)
    frames = []
    if not jobs.empty:
        with ProcessPoolExecutor(max_workers=max(1, args.max_workers)) as pool:
            futures = [pool.submit(one_direct_overlap, rec, promoters) for rec in jobs.to_dict("records")]
            for fut in as_completed(futures):
                frame = fut.result()
                if not frame.empty:
                    frames.append(frame)
    long = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["regulator_TF", "target_TF", "context", "source_file", "n_peak_overlaps", "max_overlap_bp"])
    write_csv(long, long_out)
    support = long.groupby(["regulator_TF", "context"])["source_file"].nunique().reset_index(name="n_datasets").rename(columns={"regulator_TF": "TF"})
    base = tf_universe[["TF"]].copy()
    context_wide = support.pivot_table(index="TF", columns="context", values="n_datasets", aggfunc="sum", fill_value=0.0)
    state = pd.DataFrame({"TF": base["TF"]})
    for ctx in CONTEXTS:
        score = log_max_scale01(context_wide[ctx]) if ctx in context_wide.columns else pd.Series(dtype=float)
        state[f"direct_{ctx}_score"] = base["TF"].map(score).fillna(0.0)
    remap_result_dir = args.remap_result_dir.resolve() if args.remap_result_dir is not None else raw_root / "remap_results"
    remap_support = remap_result_dir / "remap_breast_tier1_TF_support_summary.tsv"
    if remap_support.exists():
        remap = pd.read_csv(remap_support, sep="\t")
        if {"regulator_TF", "n_breast_remap_dataset_overlaps"}.issubset(remap.columns):
            remap_score = log_max_scale01(remap.set_index("regulator_TF")["n_breast_remap_dataset_overlaps"])
            breast_score = base["TF"].map(remap_score).fillna(0.0)
            state["direct_Luminal_score"] = np.maximum(state["direct_Luminal_score"], breast_score)
            state["direct_HER2_score"] = np.maximum(state["direct_HER2_score"], breast_score)
    tnbc_manifest = remap_result_dir / "TNBC_direct_binding" / "tnbc_direct_binding_dataset_manifest.tsv"
    if tnbc_manifest.exists():
        tnbc = pd.read_csv(tnbc_manifest, sep="\t")
        if {"regulator_TF", "n_peaks_in_file"}.issubset(tnbc.columns):
            usable = tnbc["usable_for_hg38_overlap"].astype(str).str.lower().isin(["true", "1", "yes"]) if "usable_for_hg38_overlap" in tnbc.columns else pd.Series(True, index=tnbc.index)
            tnbc_score = log_max_scale01(tnbc[usable].groupby("regulator_TF")["n_peaks_in_file"].sum())
            state["direct_Basal/TNBC_score"] = np.maximum(state["direct_Basal/TNBC_score"], base["TF"].map(tnbc_score).fillna(0.0))
    base["Breast_context_direct_binding"] = np.maximum(state["direct_Luminal_score"], state["direct_HER2_score"])
    base["TNBC_direct_binding"] = state["direct_Basal/TNBC_score"]
    base["context_matched_direct_binding_MAX"] = base[["Breast_context_direct_binding", "TNBC_direct_binding"]].max(axis=1)
    tnbc_counts = support[support["context"].eq("Basal/TNBC")].set_index("TF")["n_datasets"]
    base["ChIPAtlas_TNBC_direct_binding_experiments"] = base["TF"].map(tnbc_counts).fillna(0.0)
    write_csv(base, out)
    write_csv(state, state_out)
    return base, state


def build_scrna_evidence(raw_root: Path, tf_universe: pd.DataFrame, out_dir: Path, reuse: bool) -> pd.DataFrame:
    out = out_dir / "figure1a_scrna_malignant_evidence.csv"
    if reuse and out.exists():
        return pd.read_csv(out)
    h5ad = raw_root / "GSE176078_scRNA" / "data" / "Wu_etal_2021_BRCA_scRNASeq" / "GSE176078_Wu_etal_2021_BRCA_scRNASeq.h5ad"
    base = tf_universe[["TF"]].copy()
    base["scRNA_malignant"] = 0.0
    if h5ad.exists():
        adata = ad.read_h5ad(h5ad, backed="r")
        present = [g for g in base["TF"] if g in adata.var_names]
        if present:
            obs = adata.obs.copy()
            major_col = "celltype_major" if "celltype_major" in obs.columns else None
            malignant = obs[major_col].astype(str).str.lower().str.contains("cancer|malignant|epithelial", regex=True, na=False) if major_col else pd.Series(False, index=obs.index)
            counts = adata[:, present].X
            if sparse.issparse(counts):
                counts = counts.tocsr()
            totals = pd.to_numeric(obs.get("nCount_RNA", 1.0), errors="coerce").fillna(1.0).to_numpy(float)
            totals[totals <= 0] = 1.0
            if sparse.issparse(counts):
                norm = counts.multiply(1e4 / totals[:, None])
                vals = norm.copy()
                vals.data = np.log1p(vals.data)
                vals = vals.tocsr()
                malignant_mean = np.asarray(vals[malignant.to_numpy(), :].mean(axis=0)).ravel()
                other_mean = np.asarray(vals[~malignant.to_numpy(), :].mean(axis=0)).ravel()
            else:
                vals = np.log1p(np.asarray(counts) / totals[:, None] * 1e4)
                malignant_mean = vals[malignant.to_numpy(), :].mean(axis=0)
                other_mean = vals[~malignant.to_numpy(), :].mean(axis=0)
            enrich = pd.Series(np.log2((malignant_mean + 0.01) / (other_mean + 0.01)), index=present)
            base["scRNA_malignant"] = base["TF"].map(np.clip(enrich / 4.0, 0.0, 1.0)).fillna(0.0)
    write_csv(base, out)
    return base


def build_scatac_evidence(raw_root: Path, tf_universe: pd.DataFrame, out_dir: Path, reuse: bool) -> pd.DataFrame:
    out = out_dir / "figure1a_scatac_locus_evidence.csv"
    if reuse and out.exists():
        return pd.read_csv(out)
    base = tf_universe[["TF"]].copy()
    summary = raw_root / "scATAC_multiome" / "results" / "TCGA_BRCA_scATAC_CRC_TF_locus_peak_summary.tsv"
    if summary.exists():
        df = pd.read_csv(summary, sep="\t")
        gene_col = "nearestGene" if "nearestGene" in df.columns else "TF"
        counts = df.groupby(gene_col)["n_peaks"].sum()
        base["scATAC_cis_accessibility"] = base["TF"].map(scale01(counts)).fillna(0.0)
    else:
        base["scATAC_cis_accessibility"] = 0.0
    write_csv(base, out)
    return base


def build_depmap_evidence(raw_root: Path, tf_universe: pd.DataFrame, out_dir: Path, reuse: bool) -> pd.DataFrame:
    out = out_dir / "figure1a_depmap_dependency_evidence.csv"
    if reuse and out.exists():
        return pd.read_csv(out)
    gene_effect = pd.read_csv(raw_root / "DepMap" / "data" / "CRISPRGeneEffect.csv")
    model = pd.read_csv(raw_root / "DepMap" / "data" / "Model.csv")
    id_col = "Unnamed: 0" if "Unnamed: 0" in gene_effect.columns else "ModelID"
    gene_effect = gene_effect.rename(columns={id_col: "ModelID"})
    breast = model[model["OncotreeLineage"].astype(str).str.lower().str.contains("breast", na=False)]["ModelID"]
    ge = gene_effect[gene_effect["ModelID"].isin(set(breast))].copy()
    col_to_gene = {c: re.sub(r"\s+\(\d+\)$", "", c) for c in ge.columns if c != "ModelID"}
    wanted_cols = [c for c, g in col_to_gene.items() if g in set(tf_universe["TF"])]
    rows = []
    for col in wanted_cols:
        values = pd.to_numeric(ge[col], errors="coerce").dropna()
        if values.empty:
            continue
        rows.append({"TF": col_to_gene[col], "mean_gene_effect": float(values.mean()), "fraction_lt_minus0_5": float((values < -0.5).mean()), "fraction_lt_minus1": float((values < -1.0).mean())})
    dep = pd.DataFrame(rows)
    base = tf_universe[["TF"]].merge(dep, on="TF", how="left")
    strength = np.maximum(base["fraction_lt_minus0_5"].fillna(0.0), np.clip(-base["mean_gene_effect"].fillna(0.0), 0, None))
    base["DepMap_breast_dependency"] = scale01(pd.Series(strength, index=base.index))
    write_csv(base, out)
    return base


def merge_on_tf(tables: list[pd.DataFrame]) -> pd.DataFrame:
    out = tables[0].copy()
    for table in tables[1:]:
        out = out.merge(table, on="TF", how="left")
    return out


def compute_state_calls(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    primary_context = []
    active_contexts = []
    nomenclature = []
    confidence = []
    for rec in df.to_dict("records"):
        scores = {ctx: float(pd.to_numeric(rec.get(f"state_score_{ctx}", 0.0), errors="coerce")) for ctx in CONTEXTS}
        scores = {k: (0.0 if not np.isfinite(v) else v) for k, v in scores.items()}
        max_score = max(scores.values())
        second_score = sorted(scores.values())[-2]
        core = float(rec.get("Core_discovery_score_v2_0_7", 0.0))
        total = float(rec.get("Audited_total_score_0_10", 0.0))
        primaries = [ctx for ctx, value in scores.items() if value >= max_score - 0.20 and value >= 2.80]
        if not primaries and core >= HIGH_CONFIDENCE_CORE_THRESHOLD and max_score >= HIGH_CONFIDENCE_STATE_THRESHOLD:
            primaries = [max(scores, key=scores.get)]
        if not primaries and core >= 4.50 and max_score >= 2.35:
            primaries = [max(scores, key=scores.get)]
        secondaries = [ctx for ctx, value in scores.items() if ctx not in primaries and value >= 2.35 and value >= 0.82 * max_score]
        active = primaries + secondaries
        if core >= HIGH_CONFIDENCE_CORE_THRESHOLD and max_score >= HIGH_CONFIDENCE_STATE_THRESHOLD:
            tier = "High-confidence"
        elif (core >= 4.75 and max_score >= 3.00 and total - core >= 0.50) or (core >= 4.75 and max_score >= 3.50 and max_score - second_score >= 1.00) or (total >= 6.00 and max_score >= 2.40):
            tier = "Provisional"
        else:
            tier = "Supporting"
        if tier == "Supporting" or not active:
            cls = "Supporting candidate"
        elif len(active) == 1 and max_score - second_score >= 0.80:
            cls = "Subtype-specific Anchor"
        elif len(active) >= 3 and max_score - min(scores.values()) <= 0.70:
            cls = "Pan-epithelial Shared CRC"
        else:
            cls = "Subtype-enriched Hub"
        primary_context.append(";".join(primaries))
        active_contexts.append(";".join(active))
        nomenclature.append(cls)
        confidence.append(tier)
        for ctx in CONTEXTS:
            records.append({"TF": rec["TF"], "context": ctx, "state": "P" if ctx in primaries else ("S" if ctx in secondaries else ""), "state_score_0_5": scores[ctx], "state_strength_0_1": scores[ctx] / 5.0, "CRC_nomenclature": cls, "confidence_tier": tier})
    out = df.copy()
    out["primary_context"] = primary_context
    out["active_contexts"] = active_contexts
    out["CRC_nomenclature"] = nomenclature
    out["confidence_tier"] = confidence
    states = pd.DataFrame(records)
    for ctx in CONTEXTS:
        mapping = states[states["context"].eq(ctx)].set_index("TF")["state"]
        out[f"{ctx}_state"] = out["TF"].map(mapping).fillna("")
        out[f"{ctx}_member"] = out[f"{ctx}_state"].isin(["P", "S"]).astype(int)
    return out, states


def build_scores(raw_root: Path, out_dir: Path, inter_dir: Path, args: argparse.Namespace) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    inter_dir.mkdir(parents=True, exist_ok=True)
    tf_universe = build_tf_universe(raw_root, inter_dir, args.reuse_intermediate)
    gene_annotation = build_gene_annotation(raw_root, tf_universe, inter_dir, args.reuse_intermediate)
    expression, expression_state, tf_expr = build_expression_evidence(raw_root, tf_universe, inter_dir, args.reuse_intermediate, args)
    activity, activity_state = build_activity_evidence(raw_root, tf_universe, inter_dir, args.reuse_intermediate, args, tf_expr)
    atac, atac_state = build_atac_evidence(raw_root, tf_universe, inter_dir, args.reuse_intermediate)
    h3, h3_state = build_h3k27ac_evidence(raw_root, tf_universe, gene_annotation, inter_dir, args.reuse_intermediate, args)
    direct, direct_state = build_direct_binding_evidence(raw_root, tf_universe, gene_annotation, inter_dir, args.reuse_intermediate, args)
    scrna = build_scrna_evidence(raw_root, tf_universe, inter_dir, args.reuse_intermediate)
    scatac = build_scatac_evidence(raw_root, tf_universe, inter_dir, args.reuse_intermediate)
    depmap = build_depmap_evidence(raw_root, tf_universe, inter_dir, args.reuse_intermediate)
    df = merge_on_tf([tf_universe, expression, activity, atac, h3, direct, scrna, scatac, depmap])
    df["proposed_module"] = df["expression_dominant_module"]
    core_cols = ["RNA", "TF_activity_regulon", "ATAC_motif", "ATAC_footprint", "ATAC_peak_to_gene_accessibility", "context_matched_H3K27ac_SE_MAX", "context_matched_direct_binding_MAX"]
    supportive_cols = ["scRNA_malignant", "scATAC_cis_accessibility", "DepMap_breast_dependency"]
    for col in core_cols + supportive_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    df["Core_discovery_score_v2_0_7"] = df[core_cols].sum(axis=1)
    df["Core_discovery_score_0_7"] = df["Core_discovery_score_v2_0_7"]
    df["Supportive_score_0_3"] = df[supportive_cols].sum(axis=1)
    df["Audited_total_score_0_10"] = df["Core_discovery_score_v2_0_7"] + df["Supportive_score_0_3"]
    state = merge_on_tf([tf_universe, expression_state, activity_state, atac_state, h3_state, direct_state]).set_index("TF")
    for ctx in CONTEXTS:
        sub_cols = [f"expression_{ctx}_score", f"activity_{ctx}_score", f"atac_peak_{ctx}_score", f"h3_{ctx}_score", f"direct_{ctx}_score"]
        for col in sub_cols:
            if col not in state.columns:
                state[col] = 0.0
        df = df.merge(state[sub_cols].sum(axis=1).rename(f"state_score_{ctx}").reset_index(), on="TF", how="left")
        df[f"state_strength_{ctx}"] = pd.to_numeric(df[f"state_score_{ctx}"], errors="coerce").fillna(0.0) / 5.0
    df["max_state_score_0_5"] = df[[f"state_score_{c}" for c in CONTEXTS]].max(axis=1)
    df["second_state_score_0_5"] = df[[f"state_score_{c}" for c in CONTEXTS]].apply(lambda row: np.sort(row.to_numpy(dtype=float))[-2], axis=1)
    df["state_spread_max_minus_second"] = df["max_state_score_0_5"] - df["second_state_score_0_5"]
    df, states = compute_state_calls(df)
    df["H3K27ac_SE_MAX"] = df["context_matched_H3K27ac_SE_MAX"]
    df["Direct_binding_MAX"] = df["context_matched_direct_binding_MAX"]
    all_cols = ["TF", "proposed_module", "expression_dominant_module", "RNA", "TF_activity_regulon", "ATAC_motif", "ATAC_footprint", "ATAC_peak_to_gene_accessibility", "Luminal_H3K27ac_SE", "HER2_H3K27ac_SE", "Basal_TNBC_H3K27ac_SE", "context_matched_H3K27ac_SE_MAX", "Breast_context_direct_binding", "TNBC_direct_binding", "context_matched_direct_binding_MAX", "scRNA_malignant", "scATAC_cis_accessibility", "DepMap_breast_dependency", "Core_discovery_score_0_7", "Supportive_score_0_3", "Audited_total_score_0_10", "ChIPAtlas_TNBC_direct_binding_experiments", "Core_discovery_score_v2_0_7", "state_score_Basal/TNBC", "state_score_HER2", "state_score_Luminal", "state_strength_Basal/TNBC", "state_strength_HER2", "state_strength_Luminal", "max_state_score_0_5", "second_state_score_0_5", "state_spread_max_minus_second", "primary_context", "active_contexts", "CRC_nomenclature", "confidence_tier", "Luminal_state", "HER2_state", "Basal/TNBC_state", "H3K27ac_SE_MAX", "Direct_binding_MAX", "Luminal_member", "HER2_member", "Basal/TNBC_member"]
    all_cols = [c for c in all_cols if c in df.columns]
    tier_order = {"High-confidence": 0, "Provisional": 1, "Supporting": 2}
    df["tier_order"] = df["confidence_tier"].map(tier_order).fillna(9)
    full = df[all_cols + ["tier_order"]].sort_values(["tier_order", "Core_discovery_score_v2_0_7", "Audited_total_score_0_10", "TF"], ascending=[True, False, False, True]).drop(columns="tier_order")
    formal = full[(full["confidence_tier"].ne("Supporting")) & full["active_contexts"].fillna("").ne("")].copy()
    formal["context_order"] = formal["primary_context"].astype(str).map(lambda x: min([CONTEXTS.index(c) for c in CONTEXTS if c in x], default=99))
    formal = formal.sort_values(["context_order", "Core_discovery_score_v2_0_7", "Audited_total_score_0_10", "TF"], ascending=[True, False, False, True]).drop(columns="context_order")
    full.to_csv(out_dir / "figure1a_all_tf_scores.csv", index=False)
    states.to_csv(out_dir / "figure1a_context_membership_long.csv", index=False)
    formal.to_csv(out_dir / "figure1a.csv", index=False)
    audit = {
        "n_tf_universe": int(tf_universe.shape[0]),
        "n_formal_crc_tfs": int(formal.shape[0]),
        "formal_crc_tfs_generated_by_rules": formal["TF"].tolist(),
        "core_score_formula": core_cols,
        "supportive_score_formula": supportive_cols,
        "raw_root_name": raw_root.name,
        "intermediate_tables": sorted(p.name for p in inter_dir.glob("figure1a_*.csv*")),
        "no_predefined_crc_tf_list": True,
        "subtype_test_definition": "RNA and TF-activity subtype association use available raw PAM50 subtype labels; context scores are summarized into Luminal, HER2, and Basal/TNBC.",
        "legacy_candidate_filename_fallbacks": [
            "ATAC_subtype_specific_peak_to_candidate_TF_links.tsv",
            "ENCODE_MCF7_H3K27ac_candidate_TF_gene_window_summary.tsv",
            "ENCODE_MCF7_ROSE_like_SE_candidate_TF_overlap.tsv"
        ],
    }
    (out_dir / "figure1a_run_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    raw_root = find_raw_root(args.raw_root)
    build_scores(raw_root, args.output_dir, args.intermediate_dir, args)


if __name__ == "__main__":
    main()
