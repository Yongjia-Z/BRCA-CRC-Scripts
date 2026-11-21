from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg
import os

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(os.environ.get("CRC_DATA_ROOT", REPO_ROOT / "datasets"))
CORE_ROOT = Path(os.environ.get("CRC_CORE_ROOT", REPO_ROOT / "core_data"))


ROOT = (DATA_ROOT / "TCGA")
WORKSPACE = DATA_ROOT.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SCORE_COLS = [
    "pan_CRC_score",
    "luminal_CRC_score",
    "her2_CRC_score",
    "basal_CRC_score",
]

ENDPOINTS = {
    "OS": ("OS_MONTHS", "OS_STATUS"),
    "DFS": ("DFS_MONTHS", "DFS_STATUS"),
    "DSS": ("DSS_MONTHS", "DSS_STATUS"),
    "PFS": ("PFS_MONTHS", "PFS_STATUS"),
}


def parse_event(value: object) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value)
    if text.startswith("1:"):
        return 1.0
    if text.startswith("0:"):
        return 0.0
    return np.nan


def safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def clean_stage(value: object) -> str | np.nan:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().upper().replace("STAGE ", "").replace("STAGE", "")
    text = text.replace("X", "X")
    return text


def stage_to_order(stage: str) -> int:
    mapping = {
        "0": 0,
        "I": 1,
        "IA": 1,
        "IB": 1,
        "II": 2,
        "IIA": 2,
        "IIB": 2,
        "III": 3,
        "IIIA": 3,
        "IIIB": 3,
        "IIIC": 3,
        "IV": 4,
        "IVA": 4,
        "IVB": 4,
        "X": np.nan,
    }
    if pd.isna(stage):
        return np.nan
    return mapping.get(str(stage), np.nan)


def logrank_test(time: pd.Series, event: pd.Series, group: pd.Series) -> float:
    df = pd.DataFrame({"time": time, "event": event, "group": group}).dropna()
    if df["group"].nunique() < 2:
        return np.nan
    groups = sorted(df["group"].unique())
    if len(groups) != 2:
        return np.nan
    g0, g1 = groups
    times = np.sort(df.loc[df["event"] == 1, "time"].unique())
    o0 = e0 = v0 = 0.0
    for t in times:
        risk = df["time"] >= t
        events = (df["time"] == t) & (df["event"] == 1)
        n = int(risk.sum())
        d = int(events.sum())
        if n <= 1 or d == 0:
            continue
        n0 = int((risk & (df["group"] == g0)).sum())
        d0 = int((events & (df["group"] == g0)).sum())
        exp0 = d * n0 / n
        var0 = (n0 * (n - n0) * d * (n - d)) / (n**2 * (n - 1))
        o0 += d0
        e0 += exp0
        v0 += var0
    if v0 <= 0:
        return np.nan
    chi2 = (o0 - e0) ** 2 / v0
    return float(stats.chi2.sf(chi2, 1))


def zscore_series(s: pd.Series) -> pd.Series:
    s = safe_numeric(s)
    std = s.std(ddof=0)
    if pd.isna(std) or std == 0:
        return s * 0
    return (s - s.mean()) / std


def prepare_data() -> pd.DataFrame:
    score = pd.read_csv(RESULTS / "CRC_activity_scores_by_sample.tsv", sep="\t")
    patient = pd.read_csv(ROOT / "data" / "cbioportal_brca_patient_clinical_wide.tsv", sep="\t")
    patient = patient[
        [
            "patientId",
            "AGE",
            "AJCC_PATHOLOGIC_TUMOR_STAGE",
            "PATH_N_STAGE",
            "PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT",
        ]
    ].copy()
    patient = patient.rename(columns={"patientId": "patient_id"})

    df = score.merge(patient, left_on="patient_id", right_on="patient_id", how="left")
    df["OS_event"] = df["OS_STATUS"].map(parse_event)
    df["DFS_event"] = df["DFS_STATUS"].map(parse_event)
    df["DSS_event"] = df["DSS_STATUS"].map(parse_event)
    df["PFS_event"] = df["PFS_STATUS"].map(parse_event)
    df["AJCC_PATHOLOGIC_TUMOR_STAGE_clean"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(clean_stage)
    df["AJCC_PATHOLOGIC_TUMOR_STAGE_order"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE_clean"].map(stage_to_order)
    df["PATH_N_STAGE_clean"] = df["PATH_N_STAGE"].astype(str).str.upper().replace({"NAN": np.nan})
    df["LYMPH_NODE_PRESENTATION"] = df["PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT"].astype(str)
    df["LYMPH_NODE_POSITIVE"] = np.where(
        df["PATH_N_STAGE_clean"].isin(["N1", "N1A", "N1B", "N1C", "N2", "N2A", "N2B", "N3", "N3A", "N3B", "N3C"]),
        1,
        np.where(df["PATH_N_STAGE_clean"].eq("N0"), 0, np.nan),
    )
    df["AGE"] = safe_numeric(df["AGE"])
    return df


def association_table(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    clinical_vars = {
        "AJCC_PATHOLOGIC_TUMOR_STAGE_order": "ordered_stage",
        "PATH_N_STAGE_clean": "nodal_stage",
        "LYMPH_NODE_POSITIVE": "lymph_node_positive",
        "AGE": "continuous",
    }
    for score in SCORE_COLS:
        for col, typ in clinical_vars.items():
            sub = df[[score, col]].dropna()
            if sub.empty:
                continue
            if typ == "continuous":
                rho, p = stats.spearmanr(sub[score], sub[col])
                rows.append(
                    {
                        "score": score,
                        "clinical_var": col,
                        "test": "Spearman",
                        "statistic": float(rho),
                        "pvalue": float(p),
                        "n": int(sub.shape[0]),
                    }
                )
            elif typ == "ordered_stage":
                rho, p = stats.spearmanr(sub[score], sub[col])
                rows.append(
                    {
                        "score": score,
                        "clinical_var": col,
                        "test": "Spearman",
                        "statistic": float(rho),
                        "pvalue": float(p),
                        "n": int(sub.shape[0]),
                    }
                )
            else:
                groups = [g[score].dropna().values for _, g in sub.groupby(col) if len(g[score].dropna()) >= 5]
                if len(groups) >= 2:
                    stat, p = stats.kruskal(*groups)
                    rows.append(
                        {
                            "score": score,
                            "clinical_var": col,
                            "test": "Kruskal",
                            "statistic": float(stat),
                            "pvalue": float(p),
                            "n": int(sub.shape[0]),
                        }
                    )
    out = pd.DataFrame(rows)
    out["fdr"] = stats.false_discovery_control(out["pvalue"].fillna(1.0), method="bh") if not out.empty else []
    out.to_csv(RESULTS / "TCGA_CRC_clinical_association_summary.tsv", sep="\t", index=False)
    return out


def cox_summary(df: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    time_col, status_col = ENDPOINTS[endpoint]
    base = df[["sample_id"] + SCORE_COLS + [time_col, status_col, "AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE_order", "SUBTYPE", "PATH_N_STAGE_clean"]].copy()
    base[status_col] = base[status_col].map(parse_event)
    base[time_col] = safe_numeric(base[time_col])
    base = base.dropna(subset=[time_col, status_col])
    rows = []
    for score in SCORE_COLS:
        sub = base[[time_col, status_col, score]].dropna().copy()
        if sub[score].nunique() < 2 or sub.shape[0] < 20:
            continue
        sub[score] = zscore_series(sub[score])
        try:
            model = PHReg(sub[time_col], sub[[score]], status=sub[status_col])
            fit = model.fit(disp=False)
            coef = float(fit.params[0])
            se = float(fit.bse[0])
            hr = float(np.exp(coef))
            lcl = float(np.exp(coef - 1.96 * se))
            ucl = float(np.exp(coef + 1.96 * se))
            p = float(fit.pvalues[0])
        except Exception:
            coef = se = hr = lcl = ucl = p = np.nan
        rows.append(
            {
                "endpoint": endpoint,
                "model": "univariate",
                "covariate": score,
                "n": int(sub.shape[0]),
                "hazard_ratio": hr,
                "ci_lower": lcl,
                "ci_upper": ucl,
                "pvalue": p,
            }
        )

    covars = [c for c in ["pan_CRC_score", "AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE_order", "SUBTYPE", "PATH_N_STAGE_clean"] if c in base.columns]
    mv = base[[time_col, status_col, "pan_CRC_score", "AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE_order", "SUBTYPE"]].copy().dropna()
    if not mv.empty and mv["pan_CRC_score"].nunique() > 1:
        mv["pan_CRC_score"] = zscore_series(mv["pan_CRC_score"])
        mv["AGE"] = zscore_series(mv["AGE"])
        mv = pd.get_dummies(mv, columns=["SUBTYPE"], drop_first=True)
        X = mv.drop(columns=[time_col, status_col]).apply(pd.to_numeric, errors="coerce").astype(float)
        X_cols = list(X.columns)
        try:
            fit = PHReg(mv[time_col].astype(float).values, X.values, status=mv[status_col].astype(float).values).fit(disp=False)
            idx = X_cols.index("pan_CRC_score")
            coef = float(fit.params[idx])
            se = float(fit.bse[idx])
            rows.append(
                {
                    "endpoint": endpoint,
                    "model": "multivariable",
                    "covariate": "pan_CRC_score",
                    "n": int(mv.shape[0]),
                    "hazard_ratio": float(np.exp(coef)),
                    "ci_lower": float(np.exp(coef - 1.96 * se)),
                    "ci_upper": float(np.exp(coef + 1.96 * se)),
                    "pvalue": float(fit.pvalues[idx]),
                }
            )
        except Exception:
            pass
    return pd.DataFrame(rows)


def main() -> None:
    df = prepare_data()
    df.to_csv(RESULTS / "TCGA_CRC_clinical_analysis_dataset.tsv", sep="\t", index=False)

    assoc = association_table(df)

    cox_rows = []
    for endpoint in ENDPOINTS:
        cox_rows.append(cox_summary(df, endpoint))
    cox_df = pd.concat(cox_rows, ignore_index=True) if cox_rows else pd.DataFrame()
    cox_df.to_csv(RESULTS / "TCGA_CRC_cox_summary.tsv", sep="\t", index=False)

    if cox_df.empty:
        logrank = pd.DataFrame()
    else:
        logrank = cox_df[cox_df["model"] == "univariate"][["endpoint", "covariate", "hazard_ratio", "pvalue", "n"]]
    logrank.to_csv(RESULTS / "TCGA_CRC_survival_logrank_summary.tsv", sep="\t", index=False)

    summary = {
        "n_samples": int(df.shape[0]),
        "n_assoc_rows": int(assoc.shape[0]) if not assoc.empty else 0,
        "n_cox_rows": int(cox_df.shape[0]) if not cox_df.empty else 0,
    }
    (RESULTS / "TCGA_CRC_clinical_analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
