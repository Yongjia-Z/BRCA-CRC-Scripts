from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "_repro_common"))
from crc_repro import read_gene_by_sample_matrix, require_table, write_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("data") / "figure8b")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--truncation-months", type=float, default=120.0)
    return parser.parse_args()


def row_zscore(expr: pd.DataFrame) -> pd.DataFrame:
    x = expr.apply(pd.to_numeric, errors="coerce")
    return x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0.0)


def load_signatures(input_dir: Path) -> pd.DataFrame:
    sig = require_table(input_dir / "crc_state_signatures.csv", "crc_state_signatures")
    rows = []
    for state_id, group in sig.groupby("state_id", sort=False):
        context = group["subtype_context"].dropna().iloc[0] if "subtype_context" in group.columns else "Luminal"
        label = group["state_label"].dropna().iloc[0] if "state_label" in group.columns else state_id
        rows.append(
            {
                "state_id": state_id,
                "subtype_context": context,
                "state_label": label,
                "signature_gene_list": group["gene"].dropna().astype(str).tolist(),
            }
        )
    return pd.DataFrame(rows)


def event_from_status(values: pd.Series) -> pd.Series:
    text = values.astype(str)
    if text.str.contains(":").any():
        return text.str.startswith("1:").astype(int)
    return pd.to_numeric(values, errors="coerce").fillna(0).astype(int)


def first_present(columns: pd.Index, candidates: list[str]) -> str | None:
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def build_model_input(input_dir: Path, truncation_months: float) -> pd.DataFrame:
    sig = load_signatures(input_dir)
    all_genes = sorted({gene for genes in sig["signature_gene_list"] for gene in genes})

    expr = read_gene_by_sample_matrix(input_dir / "bulk_expression_matrix.csv.gz")
    z = row_zscore(expr[expr.index.isin(all_genes)])

    scores = pd.DataFrame({"sample_id": z.columns})
    for _, row in sig.iterrows():
        genes = [gene for gene in row["signature_gene_list"] if gene in z.index]
        scores[row["state_id"]] = z.loc[genes].mean(axis=0).to_numpy() if genes else np.nan

    clinical = require_table(input_dir / "bulk_clinical_survival.csv", "bulk_clinical_survival")
    sample_col = first_present(clinical.columns, ["sample_id", "sample", "Sample"]) or clinical.columns[0]
    clinical = clinical.rename(columns={sample_col: "sample_id"})
    if "sample_type" in clinical.columns:
        clinical = clinical[clinical["sample_type"].astype(str).eq("Primary Tumor")].copy()

    data = clinical.merge(scores, on="sample_id", how="inner")
    endpoints = [
        ("OS", ["OS_MONTHS", "OS_time", "time_months"], ["OS_STATUS", "OS_event", "status_raw"]),
        ("PFS", ["PFS_MONTHS", "PFS_time"], ["PFS_STATUS", "PFS_event"]),
    ]

    rows = []
    for _, state in sig.iterrows():
        subset = data.copy()
        if "SUBTYPE" in subset.columns:
            if state["subtype_context"] == "Luminal":
                subset = subset[subset["SUBTYPE"].isin(["BRCA_LumA", "BRCA_LumB"])].copy()
            elif state["subtype_context"] == "HER2":
                subset = subset[subset["SUBTYPE"].eq("BRCA_Her2")].copy()
            else:
                subset = subset[subset["SUBTYPE"].eq("BRCA_Basal")].copy()

        for endpoint, time_candidates, status_candidates in endpoints:
            time_col = first_present(subset.columns, time_candidates)
            status_col = first_present(subset.columns, status_candidates)
            if not time_col or not status_col:
                continue

            temp = subset[["sample_id", state["state_id"], time_col, status_col]].copy()
            temp = temp.rename(columns={state["state_id"]: "state_score", time_col: "time_months", status_col: "status_raw"})
            temp["event"] = event_from_status(temp["status_raw"])
            temp["state_score"] = pd.to_numeric(temp["state_score"], errors="coerce")
            temp["time_months"] = pd.to_numeric(temp["time_months"], errors="coerce")
            temp = temp.dropna(subset=["state_score", "time_months", "event"])
            temp = temp[temp["time_months"].gt(0)].copy()
            if temp.shape[0] < 5:
                continue

            sd = temp["state_score"].std(ddof=0)
            temp["score_z"] = (temp["state_score"] - temp["state_score"].mean()) / sd if sd > 0 else np.nan
            temp["time_truncated"] = temp["time_months"].clip(upper=truncation_months)
            temp["event_truncated"] = np.where(temp["time_months"].gt(truncation_months), 0, temp["event"])
            temp["cohort"] = "TCGA-BRCA"
            temp["endpoint"] = endpoint
            temp["state_id"] = state["state_id"]
            temp["state_label"] = state["state_label"]
            temp["subtype_context"] = state["subtype_context"]
            rows.append(temp)

    if not rows:
        raise ValueError("No survival rows generated; check clinical column names and signature overlap")
    return pd.concat(rows, ignore_index=True)


def write_cox_script(path: Path) -> None:
    path.write_text(
        r'''
suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly=TRUE)
inp <- args[1]
out <- args[2]
trunc_months <- as.numeric(args[3])
d <- read.delim(inp, check.names=FALSE)
keys <- unique(d[, c("cohort", "endpoint", "state_id", "state_label", "subtype_context")])
rows <- list()
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  sub <- d[d$cohort == k$cohort & d$endpoint == k$endpoint & d$state_id == k$state_id, ]
  sub <- sub[is.finite(sub$score_z) & is.finite(sub$time_truncated) & is.finite(sub$event_truncated), ]
  n <- nrow(sub)
  events <- sum(sub$event_truncated, na.rm=TRUE)
  rec <- data.frame(cohort=k$cohort, endpoint=k$endpoint, state_id=k$state_id, state_label=k$state_label,
                    subtype_context=k$subtype_context, n=n, events=events, hr_per_1sd=NA_real_,
                    ci_low=NA_real_, ci_high=NA_real_, pvalue=NA_real_, ph_test_p=NA_real_,
                    followup_truncation_months=trunc_months, status="insufficient_events")
  if (n >= 30 && events >= 5 && sd(sub$score_z, na.rm=TRUE) > 0) {
    fit <- try(coxph(Surv(time_truncated, event_truncated) ~ score_z, data=sub, ties="breslow"), silent=TRUE)
    if (!inherits(fit, "try-error")) {
      sm <- summary(fit)
      zph <- try(cox.zph(fit), silent=TRUE)
      rec$hr_per_1sd <- sm$coefficients["score_z", "exp(coef)"]
      rec$ci_low <- sm$conf.int["score_z", "lower .95"]
      rec$ci_high <- sm$conf.int["score_z", "upper .95"]
      rec$pvalue <- sm$coefficients["score_z", "Pr(>|z|)"]
      if (!inherits(zph, "try-error")) {
        rec$ph_test_p <- zph$table["score_z", "p"]
      }
      rec$status <- "ok"
    } else {
      rec$status <- "cox_failed"
    }
  }
  rows[[length(rows) + 1]] <- rec
}
res <- do.call(rbind, rows)
write.table(res, out, sep="\t", quote=FALSE, row.names=FALSE)
''',
        encoding="utf-8",
    )


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    model_input = build_model_input(args.input_dir, args.truncation_months)
    input_path = args.output_dir / "figure8b_tcga_cox_model_input.tsv"
    result_path = args.output_dir / "figure8b_tcga_os_pfs_cox_ph_truncated_results.tsv"
    r_script = args.output_dir / "figure8b_cox_ph_truncated.R"

    model_input.to_csv(input_path, sep="\t", index=False)
    write_cox_script(r_script)
    subprocess.run(["Rscript", str(r_script), str(input_path), str(result_path), str(args.truncation_months)], check=True)

    res_df = pd.read_csv(result_path, sep="\t")
    write_table(res_df, args.output_dir / "figure8b.csv")


if __name__ == "__main__":
    main()
