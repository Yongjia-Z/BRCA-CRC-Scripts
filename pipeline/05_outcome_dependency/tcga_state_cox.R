suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly=TRUE)
inp <- args[1]
out <- args[2]
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
                    followup_truncation_months=120, status="insufficient_events")
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
