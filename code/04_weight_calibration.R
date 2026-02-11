# 04_weight_calibration.R — robust w-hat from simulation percentiles (fixed)
# Inputs:  out/sim_recovery_percentiles.csv  (from 01_fake_data_recovery.R)
# Outputs: out/weight_calibration.txt, out/weight_calibration.png

suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
dir.create("out", showWarnings=FALSE)

# ---- 1) Load percentiles ----------------------------------------------------
perc <- read.csv("out/sim_recovery_percentiles.csv", stringsAsFactors=FALSE)
stopifnot(all(c("w","metric","empirical","percentile") %in% names(perc)))
targets <- c("P_each","P_onean","Delta_mean")
perc <- perc |>
  filter(metric %in% targets) |>
  mutate(w = as.numeric(w),
         percentile = as.numeric(percentile)) |>
  arrange(metric, w)

# sanity echo (optional)
write.csv(perc, "out/weight_calibration_input.csv", row.names=FALSE)

# ---- 2) Invert percentile -> w at 0.5 (with explicit sorting by percentile) --
w_at_0p5_extrap <- function(w, p){
  ok <- is.finite(w) & is.finite(p)
  w <- w[ok]; p <- p[ok]
  if(length(w) < 2) return(NA_real_)
  # sort by percentile (x) so approx() has increasing x
  o <- order(p)           # <---- critical fix: sort by p, not by w
  p <- p[o]; w <- w[o]
  # handle duplicate (p,w) pairs
  dup <- !duplicated(cbind(p,w))
  p <- p[dup]; w <- w[dup]
  if(length(p) < 2) return(NA_real_)
  # interpolate if 0.5 inside; otherwise extrapolate from nearest edge
  if(min(p) <= 0.5 && max(p) >= 0.5){
    return(approx(x = p, y = w, xout = 0.5, ties = "ordered")$y)
  }
  # choose two closest p's to 0.5 for a stable extrapolation
  idx <- order(abs(p - 0.5))[1:2]
  p1 <- p[idx[1]]; p2 <- p[idx[2]]; w1 <- w[idx[1]]; w2 <- w[idx[2]]
  w1 + (0.5 - p1) * (w2 - w1) / (p2 - p1)
}

w_table <- perc |>
  group_by(metric) |>
  summarise(
    w_hat = w_at_0p5_extrap(w, percentile),
    method = case_when(
      min(percentile, na.rm=TRUE) <= 0.5 & max(percentile, na.rm=TRUE) >= 0.5 ~ "interpolation",
      min(percentile, na.rm=TRUE) >  0.5                                      ~ "extrapolation (low)",
      max(percentile, na.rm=TRUE) <  0.5                                      ~ "extrapolation (high)",
      TRUE                                                                    ~ "unknown"
    ),
    .groups = "drop"
  )

# ---- 3) Two consensuses -----------------------------------------------------
consensus_median <- median(w_table$w_hat, na.rm=TRUE)

# SSE-minimising consensus: piecewise-linear per metric, optimise w in [0,1]
interp_fun <- function(rows){
  function(wq){
    wr <- rows$w; pr <- rows$percentile
    ok <- is.finite(wr) & is.finite(pr); wr <- wr[ok]; pr <- pr[ok]
    if(length(wr) < 2) return(NA_real_)
    o <- order(wr); wr <- wr[o]; pr <- pr[o]
    if(min(wr) <= wq && wq <= max(wr)){
      approx(x=wr, y=pr, xout=wq, ties="ordered")$y
    } else if(wq < min(wr)){
      i <- order(wr)[1:2]; p1 <- pr[i[1]]; p2 <- pr[i[2]]; w1 <- wr[i[1]]; w2 <- wr[i[2]]
      p1 + (wq - w1) * (p2 - p1) / (w2 - w1)
    } else {
      i <- order(wr, decreasing=TRUE)[1:2]; p1 <- pr[i[1]]; p2 <- pr[i[2]]; w1 <- wr[i[1]]; w2 <- wr[i[2]]
      p1 + (wq - w1) * (p2 - p1) / (w2 - w1)
    }
  }
}
per_metric_fun <- purrr::map(split(perc, perc$metric), interp_fun)
obj <- function(wq){
  vals <- vapply(per_metric_fun, function(f) f(wq), numeric(1))
  mean((vals - 0.5)^2)
}
consensus_sse <- optimize(obj, interval=c(0,1))$minimum

# ---- 4) Write results -------------------------------------------------------
sink("out/weight_calibration.txt")
cat("w-hat by metric (target percentile = 0.50)\n")
for(i in seq_len(nrow(w_table))){
  cat(sprintf("  %-12s : %.3f   [%s]\n",
              w_table$metric[i], w_table$w_hat[i], w_table$method[i]))
}
cat(sprintf("\nConsensus w-hat (median of per-metric hats): %.3f\n", consensus_median))
cat(sprintf("Consensus w-hat (min SSE across metrics):   %.3f\n", consensus_sse))
sink()

# ---- 5) Plot with per-metric hats + both consensuses ------------------------
plt <- perc |>
  ggplot(aes(w, percentile, color=metric)) +
  geom_hline(yintercept=0.5, linetype=2) +
  geom_point(size=3) +
  geom_line() +
  # per-metric hats
  geom_vline(data=w_table, aes(xintercept=w_hat, color=metric), linetype=3, show.legend=FALSE) +
  # consensuses
  geom_vline(xintercept=consensus_median, linetype=3) +
  geom_vline(xintercept=consensus_sse,  linetype=3, alpha=.6) +
  scale_x_continuous(breaks=sort(unique(perc$w))) +
  labs(title="Calibration to percentiles across mixture weights",
       x="Mixture weight w (pronoun share)", y="Empirical percentile\nwithin simulated distribution") +
  theme_minimal()
ggsave("out/weight_calibration.png", plt, width=7, height=4.8, dpi=300, bg="white")

# Console echo
print(w_table)
cat(sprintf("Consensus (median): %.3f\n", consensus_median))
cat(sprintf("Consensus (min SSE): %.3f\n", consensus_sse))
