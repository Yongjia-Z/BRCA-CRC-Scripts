from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CONTEXTS = ["Luminal", "HER2", "Basal/TNBC"]
PAM50_ORDER = ["LumA", "LumB", "HER2", "Basal", "Normal"]
PAM50_CONTEXT_ORDER = ["Luminal", "HER2", "Basal/TNBC", "Normal"]

def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    if path.name.endswith(".csv.gz"):
        return pd.read_csv(path, compression="gzip")
    if path.name.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", compression="gzip")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() == ".json":
        return pd.json_normalize(json.loads(path.read_text(encoding="utf-8")))
    return pd.read_csv(path, sep=None, engine="python")

def require_table(path: Path, label: str) -> pd.DataFrame:
    df = read_table(path)
    if df.empty:
        raise FileNotFoundError(f"required input missing or empty: {label} ({path})")
    return df

def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".tsv":
        df.to_csv(path, sep="\t", index=False)
    else:
        df.to_csv(path, index=False)

def first_col(df: pd.DataFrame, names: list[str]) -> str | None:
    lookup = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lookup:
            return lookup[name.lower()]
    return None

def numeric(df: pd.DataFrame, names: list[str], default: float = 0.0) -> pd.Series:
    col = first_col(df, names)
    if col is None:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default).astype(float)

def as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().isin(["true", "1", "yes", "y"]).astype(float)

def normalize_context(value: object) -> str:
    text = str(value).strip().lower().replace("_", " ").replace("-", " ")
    if text in {"luma", "lum a", "luminal a", "lumb", "lum b", "luminal b", "er+", "er positive", "luminal"}:
        return "Luminal"
    if "her2" in text:
        return "HER2"
    if "basal" in text or "tnbc" in text or "triple negative" in text:
        return "Basal/TNBC"
    if "normal" in text:
        return "Normal"
    return str(value).strip()

def normalize_pam50(value: object) -> str:
    text = str(value).strip().lower().replace("_", "").replace("-", "")
    if text in {"luma", "luminala"}:
        return "LumA"
    if text in {"lumb", "luminalb"}:
        return "LumB"
    if "her2" in text:
        return "HER2"
    if "basal" in text or "tnbc" in text:
        return "Basal"
    if "normal" in text:
        return "Normal"
    return str(value).strip()

def collapse_pam50_to_context(value: object) -> str:
    text = str(value).strip().lower().replace("_", "").replace("-", "").replace(" ", "")
    if text in {
        "luma", "brcaluma", "luminala", "brcaluminala",
        "lumb", "brcalumb", "luminalb", "brcaluminalb",
        "luminal", "brcaluminal", "er+", "erpositive",
    }:
        return "Luminal"
    if "her2" in text:
        return "HER2"
    if "basal" in text or "tnbc" in text or "triplenegative" in text:
        return "Basal/TNBC"
    if "normal" in text:
        return "Normal"
    return str(value).strip()

def standardize_tf_table(df: pd.DataFrame, names: list[str] = ["TF", "gene", "gene_name", "symbol"]) -> pd.DataFrame:
    tf_col = first_col(df, names)
    if tf_col is None:
        raise ValueError("TF column not found")
    out = df.rename(columns={tf_col: "TF"}).copy()
    out["TF"] = out["TF"].astype(str).str.strip()
    return out

def load_crc_tfs(input_dir: Path, filename: str = "crc_tf_table.csv") -> pd.DataFrame:
    candidates = [
        input_dir / filename,
        input_dir / "figure1a.csv",
        Path(filename),
        Path("figure1a.csv"),
        Path("figure1") / "Result & Script" / "figure1a.csv",
        Path("figure1") / "Result & Script" / filename,
    ]
    path = next((p for p in candidates if p.exists()), candidates[0])
    df = require_table(path, path.name)
    df = standardize_tf_table(df)
    df["_input_order"] = np.arange(df.shape[0])
    order_col = first_col(df, ["display_order", "plot_order", "order", "TF_order"])
    if order_col is not None:
        df["_sort_order"] = pd.to_numeric(df[order_col], errors="coerce").fillna(df["_input_order"])
    else:
        df["_sort_order"] = df["_input_order"]
    out = df.sort_values(["_sort_order", "_input_order"]).drop_duplicates("TF")
    return out.drop(columns=["_input_order", "_sort_order"])

def ordered_tfs(crc: pd.DataFrame) -> list[str]:
    return standardize_tf_table(crc)["TF"].dropna().astype(str).drop_duplicates().tolist()

def genes_for_context(crc: pd.DataFrame, context: str) -> list[str]:
    crc = standardize_tf_table(crc)
    col = f"{context}_member"
    if col not in crc.columns:
        return ordered_tfs(crc)
    keep = pd.to_numeric(crc[col], errors="coerce").fillna(0).astype(int).eq(1)
    return crc.loc[keep, "TF"].dropna().astype(str).drop_duplicates().tolist()

def crc_tf_module_scores(crc: pd.DataFrame, expression: pd.DataFrame, activity: pd.DataFrame | None = None) -> pd.DataFrame:
    common = expression.columns if activity is None else expression.columns.intersection(activity.columns)
    modules = {
        "Pan CRC TFs": ordered_tfs(crc),
        "Luminal CRC TFs": genes_for_context(crc, "Luminal"),
        "HER2 CRC TFs": genes_for_context(crc, "HER2"),
        "Basal/TNBC CRC TFs": genes_for_context(crc, "Basal/TNBC"),
    }
    rows = []
    for label, genes in modules.items():
        expr_score = module_score(expression[common], genes, zscore_genes=True)
        if activity is not None:
            act_score = module_score(activity[common], genes, zscore_genes=True)
            score = pd.concat([expr_score.rename("expression"), act_score.rename("activity")], axis=1).mean(axis=1)
        else:
            score = expr_score
        rows.append(score.rename(label))
    return pd.concat(rows, axis=1)

def load_or_compute_crc_tf_scores(input_dir: Path) -> pd.DataFrame:
    score_path = input_dir / "tcga_crc_tf_scores_by_sample.csv"
    if score_path.exists():
        scores = require_table(score_path, score_path.name)
        sample_col = first_col(scores, ["sample_id", "Sample", "sample"])
        if sample_col is None:
            sample_col = scores.columns[0]
        return scores.rename(columns={sample_col: "sample_id"}).set_index("sample_id")
    crc = load_crc_tfs(input_dir)
    expression = read_gene_by_sample_matrix(input_dir / "tcga_expression_matrix.csv.gz")
    activity_path = input_dir / "tcga_tf_activity_matrix.csv.gz"
    activity = read_gene_by_sample_matrix(activity_path) if activity_path.exists() else None
    return crc_tf_module_scores(crc, expression, activity)

def rank01(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    return x.rank(method="average", pct=True).fillna(0.0).clip(0.0, 1.0)

def scale01(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float)
    if len(x) == 0:
        return x
    lo = float(x.min())
    hi = float(x.max())
    if hi <= lo:
        return pd.Series(np.where(x.gt(0), 1.0, 0.0), index=values.index)
    return ((x - lo) / (hi - lo)).clip(0.0, 1.0)

def log_scale01(values: pd.Series) -> pd.Series:
    return scale01(np.log1p(pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float)))

def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").fillna(1.0).clip(0.0, 1.0)
    n = len(p)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(out, index=p.index)

def row_zscore(matrix: pd.DataFrame, clip: float | None = None) -> pd.DataFrame:
    x = matrix.apply(pd.to_numeric, errors="coerce")
    mean = x.mean(axis=1)
    sd = x.std(axis=1).replace(0, np.nan)
    z = x.sub(mean, axis=0).div(sd, axis=0).fillna(0.0)
    if clip is not None:
        z = z.clip(-clip, clip)
    return z

def col_zscore(matrix: pd.DataFrame) -> pd.DataFrame:
    x = matrix.apply(pd.to_numeric, errors="coerce")
    mean = x.mean(axis=0)
    sd = x.std(axis=0).replace(0, np.nan)
    return x.sub(mean, axis=1).div(sd, axis=1).fillna(0.0)

def read_gene_by_sample_matrix(path: Path) -> pd.DataFrame:
    df = require_table(path, path.name)
    gene_col = first_col(df, ["TF", "gene", "gene_name", "symbol"])
    if gene_col is None:
        gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene"}).drop_duplicates("gene")
    return df.set_index("gene")

def module_score(matrix: pd.DataFrame, genes: list[str], zscore_genes: bool = True) -> pd.Series:
    present = [g for g in genes if g in matrix.index]
    if not present:
        return pd.Series(np.nan, index=matrix.columns)
    x = matrix.loc[present]
    if zscore_genes:
        x = row_zscore(x)
    return x.mean(axis=0)

def kruskal_by_gene(expr_long: pd.DataFrame, gene_col: str, value_col: str, group_col: str) -> pd.DataFrame:
    rows = []
    for gene, group in expr_long.groupby(gene_col, sort=False):
        arrays = [g[value_col].dropna().to_numpy(dtype=float) for _, g in group.groupby(group_col) if g[value_col].dropna().shape[0] > 1]
        if len(arrays) < 2 or sum(np.nanstd(a) > 0 for a in arrays) == 0:
            p = 1.0
            stat = 0.0
        else:
            try:
                res = stats.kruskal(*arrays)
                p = float(res.pvalue)
                stat = float(res.statistic)
            except ValueError:
                p = 1.0
                stat = 0.0
        rows.append({"TF": gene, "kruskal_statistic": stat, "pvalue": p})
    out = pd.DataFrame(rows)
    out["fdr"] = bh_fdr(out["pvalue"])
    return out

def spearman_rows(left: pd.DataFrame, right: pd.DataFrame, left_name: str, right_name: str) -> pd.DataFrame:
    rows = []
    common = left.columns.intersection(right.columns)
    for a in left.index:
        for b in right.index:
            x = pd.to_numeric(left.loc[a, common], errors="coerce")
            y = pd.to_numeric(right.loc[b, common], errors="coerce")
            mask = x.notna() & y.notna()
            if mask.sum() < 4:
                rho = np.nan
                p = np.nan
            else:
                rho, p = stats.spearmanr(x[mask], y[mask])
            rows.append({left_name: a, right_name: b, "spearman_rho": rho, "pvalue": p, "n": int(mask.sum())})
    out = pd.DataFrame(rows)
    out["fdr"] = bh_fdr(out["pvalue"].fillna(1.0))
    return out

def residualize_by_category(values: pd.Series, category: pd.Series) -> pd.Series:
    y = pd.to_numeric(values, errors="coerce")
    design = pd.get_dummies(category.fillna("Unknown").astype(str), drop_first=False, dtype=float)
    design.insert(0, "intercept", 1.0)
    mask = y.notna()
    if mask.sum() <= design.shape[1]:
        return y * np.nan
    beta = np.linalg.lstsq(design.loc[mask].to_numpy(float), y.loc[mask].to_numpy(float), rcond=None)[0]
    fitted = pd.Series(design.to_numpy(float) @ beta, index=y.index)
    return y - fitted

def read_bed(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.name.endswith(".gz") else open
    chroms = []
    starts = []
    ends = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                chroms.append(parts[0])
                starts.append(int(float(parts[1])))
                ends.append(int(float(parts[2])))
            except ValueError:
                continue
    return pd.DataFrame({"chrom": chroms, "start": starts, "end": ends})

def edge_score_column(edges: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    if "context_edge_score" in edges.columns:
        return edges, "context_edge_score"
    if "TF_TF_edge_evidence_score" in edges.columns:
        return edges, "TF_TF_edge_evidence_score"
    evidence_cols = [c for c in edges.columns if str(c).endswith("evidence") or str(c).endswith("_evidence")]
    if not evidence_cols:
        raise ValueError("edge evidence table requires context_edge_score, TF_TF_edge_evidence_score, or evidence-layer columns")
    out = edges.copy()
    out["context_edge_score"] = out[evidence_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum(axis=1)
    return out, "context_edge_score"

def read_manual_layout(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["TF", "x", "y"])
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        payload = None
    if isinstance(payload, dict) and isinstance(payload.get("nodes"), dict):
        rows = []
        for tf, rec in payload["nodes"].items():
            if isinstance(rec, dict):
                rows.append({"TF": tf, "x": rec.get("x"), "y": rec.get("y")})
        return pd.DataFrame(rows)
    table = read_table(path)
    if {"TF", "x", "y"}.issubset(table.columns):
        return table[["TF", "x", "y"]]
    return pd.DataFrame(columns=["TF", "x", "y"])

def count_region_overlaps(peaks: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if peaks.empty or regions.empty:
        return pd.DataFrame(columns=["TF", "n_peak_overlaps", "max_overlap_bp"])
    for chrom, sub_regions in regions.groupby("chrom", sort=False):
        sub_peaks = peaks[peaks["chrom"].eq(chrom)]
        if sub_peaks.empty:
            continue
        sub_peaks = sub_peaks.sort_values("start")
        starts = sub_peaks["start"].to_numpy(np.int64)
        ends = sub_peaks["end"].to_numpy(np.int64)
        r_start = sub_regions["region_start"].to_numpy(np.int64)
        r_end = sub_regions["region_end"].to_numpy(np.int64)
        limits = np.searchsorted(starts, r_end, side="left")
        for tf, left, right, limit in zip(sub_regions["TF"].to_numpy(), r_start, r_end, limits):
            if limit == 0:
                continue
            idx = np.nonzero(ends[:limit] > left)[0]
            if idx.size == 0:
                continue
            ov = np.minimum(ends[idx], right) - np.maximum(starts[idx], left)
            rows.append({"TF": tf, "n_peak_overlaps": int(idx.size), "max_overlap_bp": int(ov.max())})
    return pd.DataFrame(rows)

def gene_windows(gene_annotation: pd.DataFrame, genes: list[str], flank: int) -> pd.DataFrame:
    anno = standardize_tf_table(gene_annotation)
    chrom_col = first_col(anno, ["chrom", "chromosome", "chr"])
    start_col = first_col(anno, ["start", "gene_start", "tx_start"])
    end_col = first_col(anno, ["end", "gene_end", "tx_end"])
    strand_col = first_col(anno, ["strand"])
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError("gene annotation needs chromosome, start, and end")
    anno = anno[anno["TF"].isin(genes)].drop_duplicates("TF").copy()
    anno["chrom"] = anno[chrom_col].astype(str)
    anno["start"] = pd.to_numeric(anno[start_col], errors="coerce").astype(int)
    anno["end"] = pd.to_numeric(anno[end_col], errors="coerce").astype(int)
    anno["strand"] = anno[strand_col].astype(str) if strand_col else "+"
    anno["tss"] = np.where(anno["strand"].eq("-"), anno["end"], anno["start"])
    out = anno[["TF", "chrom", "start", "end", "strand", "tss"]].copy()
    out["region_start"] = (out["start"] - flank).clip(lower=0).astype(int)
    out["region_end"] = (out["end"] + flank).astype(int)
    return out

def promoter_windows(gene_annotation: pd.DataFrame, genes: list[str], flank: int) -> pd.DataFrame:
    windows = gene_windows(gene_annotation, genes, 0)
    windows["region_start"] = (windows["tss"] - flank).clip(lower=0).astype(int)
    windows["region_end"] = (windows["tss"] + flank).astype(int)
    return windows

def continuous_cox(df: pd.DataFrame, time_col: str, event_col: str, score_col: str, covariates: list[str] | None = None) -> dict:
    try:
        from lifelines import CoxPHFitter
    except Exception as exc:
        raise RuntimeError("lifelines is required for Cox analysis") from exc
    cols = [time_col, event_col, score_col] + list(covariates or [])
    data = df[cols].copy()
    data[score_col] = pd.to_numeric(data[score_col], errors="coerce")
    sd = data[score_col].std()
    if not np.isfinite(sd) or sd == 0:
        return {"n": int(data.dropna().shape[0]), "events": int(pd.to_numeric(data[event_col], errors="coerce").fillna(0).sum()), "hr_per_1sd": np.nan, "ci_low": np.nan, "ci_high": np.nan, "pvalue": np.nan, "status": "constant_score"}
    data[score_col] = (data[score_col] - data[score_col].mean()) / sd
    for col in covariates or []:
        if data[col].dtype == object or str(data[col].dtype).startswith("category"):
            dummies = pd.get_dummies(data[col].astype(str), prefix=col, drop_first=True, dtype=float)
            data = pd.concat([data.drop(columns=[col]), dummies], axis=1)
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < 20 or pd.to_numeric(data[event_col], errors="coerce").sum() < 5:
        return {"n": int(data.shape[0]), "events": int(pd.to_numeric(data[event_col], errors="coerce").fillna(0).sum()), "hr_per_1sd": np.nan, "ci_low": np.nan, "ci_high": np.nan, "pvalue": np.nan, "status": "insufficient_events"}
    model = CoxPHFitter()
    model.fit(data, duration_col=time_col, event_col=event_col)
    row = model.summary.loc[score_col]
    return {"n": int(data.shape[0]), "events": int(data[event_col].sum()), "hr_per_1sd": float(row["exp(coef)"]), "ci_low": float(row["exp(coef) lower 95%"]), "ci_high": float(row["exp(coef) upper 95%"]), "pvalue": float(row["p"]), "status": "ok"}

def km_groups(scores: pd.Series, method: str = "quartile") -> pd.Series:
    x = pd.to_numeric(scores, errors="coerce")
    out = pd.Series("Middle", index=x.index, dtype=object)
    if method == "median":
        med = x.median()
        out[x <= med] = "Low"
        out[x > med] = "High"
    else:
        q1 = x.quantile(0.25)
        q3 = x.quantile(0.75)
        out[x <= q1] = "Low"
        out[x >= q3] = "High"
    out[x.isna()] = np.nan
    return out
