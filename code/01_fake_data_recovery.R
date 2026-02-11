# 01_fake_data_recovery.R — simulation-recovery with hard-coded mixture sweep
# Mixture weights (pronoun share): 0.35, 0.50, 0.65
# Inputs:  matrix_clean.csv  [required]
#          feature_blocks_by_order.csv  [optional; not required for this script]
# Outputs: out/sim_recovery_metrics_w035.csv, ... w050.csv, ... w065.csv
#          out/sim_recovery_metrics_all.csv
#          out/sim_recovery_summary_w035.csv, ... w050.csv, ... w065.csv
#          out/sim_recovery_summary_all.csv
#          out/sim_recovery_percentiles.csv
#          out/sim_recovery_plots_w035.png, ... w050.png, ... w065.png

suppressPackageStartupMessages({
  library(tidyverse); library(Matrix); library(proxy); library(vegan); library(glmnet); library(ggplot2)
})

dir.create("out", showWarnings=FALSE)

# ---------- Load observed matrix ----------
dat <- read.csv("matrix_clean.csv", stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("lemma","class") %in% names(dat)))
items <- dat$lemma
y_raw <- factor(dat$class, levels=c("determinative","pronoun"))
X_obs <- as.matrix(dat[, setdiff(names(dat), c("lemma","class"))])
storage.mode(X_obs) <- "numeric"; X_obs[is.na(X_obs)] <- 0
rownames(X_obs) <- items

recips <- intersect(c("each_other","one_another"), items)
if(length(recips)!=2) stop("Could not find both reciprocals in lemma column.")

# Sanity check: reciprocals' rows are not identical; report how many cells differ
row_diff_count <- sum(X_obs["each_other", ] != X_obs["one_another", ])
identical_rows <- row_diff_count == 0
writeLines(sprintf("Reciprocal rows identical: %s; differing feature count: %d", identical_rows, row_diff_count))
if(identical_rows) warning("each_other and one_another are identical in the matrix — unexpected.")

# Label sets
pron_set <- items[y_raw=="pronoun"]
det_set  <- items[y_raw=="determinative"]
fused_col <- names(dat)[grepl("fused", names(dat), ignore.case=TRUE)][1]
is_fused <- if(!is.na(fused_col)) dat[[fused_col]]==1 else rep(FALSE, nrow(dat))
fused_set <- items[y_raw=="determinative" & is_fused]
if(length(fused_set)<2) stop("Fused-determinative set not detected (need at least 2).")

# Helper functions
jacc <- function(M) proxy::dist(M, method="Jaccard")
as_mat <- function(D){ M <- as.matrix(D); rownames(M) <- colnames(M) <- attr(D,"Labels"); M }
mean_to <- function(D, src, tgt){ M <- as_mat(D); mean(M[src, tgt, drop=FALSE]) }
delta_to_fused <- function(D, recip, pron_targets, fused_targets){ mean_to(D, recip, pron_targets) - mean_to(D, recip, fused_targets) }
nearest_of <- function(D, item, k=5){ M <- as_mat(D); d <- M[item, ]; setdiff(names(sort(d)), item)[1:k] }

# Empirical reference quantities for overlays
D_emp <- jacc(X_obs)
emp_delta_each  <- delta_to_fused(D_emp, "each_other", pron_set, fused_set)
emp_delta_onean <- delta_to_fused(D_emp, "one_another", pron_set, fused_set)
emp_delta_mean  <- mean(c(emp_delta_each, emp_delta_onean))

# Empirical ridge probabilities for reciprocals (read if available; otherwise compute)
emp_ridge_path <- "recip_pred_probs.csv"
if(file.exists(emp_ridge_path)){
  emp <- read.csv(emp_ridge_path, stringsAsFactors=FALSE)
  # ridge alpha name may be "ridge" or alpha==0; handle both
  if("alpha" %in% names(emp)){
    emp_ridge <- emp |> dplyr::filter(alpha == 0 | grepl("ridge", as.character(alpha))) |> dplyr::select(lemma, P_pron)
  } else {
    emp_ridge <- emp |> dplyr::select(lemma, P_pron)
  }
  emp_e <- emp_ridge$P_pron[match("each_other", emp_ridge$lemma)]
  emp_o <- emp_ridge$P_pron[match("one_another", emp_ridge$lemma)]
} else {
  nonrec <- !(items %in% recips)
  fit_emp <- glmnet::cv.glmnet(x = X_obs[nonrec,], y = y_raw[nonrec], family="binomial", alpha=0, nfolds=10, type.measure="deviance")
  pr_emp  <- as.numeric(predict(fit_emp, newx=X_obs[items %in% recips,], s="lambda.min", type="response"))
  names(pr_emp) <- recips
  emp_e <- pr_emp["each_other"]; emp_o <- pr_emp["one_another"]
}
writeLines(sprintf("Empirical ridge P(pronoun): each_other=%.3f ; one_another=%.3f", emp_e, emp_o))

# ---------- Empirical per-feature probabilities by class (for simulation) ----------
p_pron <- colMeans(X_obs[y_raw=="pronoun", , drop=FALSE])
p_det  <- colMeans(X_obs[y_raw=="determinative", , drop=FALSE])
# shrink away from 0/1 with a simple Laplace smoothing
J <- ncol(X_obs)
shrink <- function(p, n=length(p)) (p*n + 1) / (n + 2)
p_pron <- shrink(p_pron, n=sum(y_raw=="pronoun"))
p_det  <- shrink(p_det,  n=sum(y_raw=="determinative"))

# ---------- Simulation design ----------
set.seed(2025)
R <- 500            # replicates per weight
k <- 5              # k-NN size
weights <- c(0.35, 0.50, 0.65)

simulate_matrix <- function(w_boundary){
  I <- nrow(X_obs); J <- ncol(X_obs)
  P <- matrix(0.0, I, J); rownames(P) <- items; colnames(P) <- colnames(X_obs)
  P[y_raw=="pronoun", ]       <- matrix(rep(p_pron, each=sum(y_raw=="pronoun")),       ncol=J)
  P[y_raw=="determinative", ] <- matrix(rep(p_det,  each=sum(y_raw=="determinative")), ncol=J)
  # mixture for reciprocals: w * pron + (1-w) * det
  P[recips, ] <- w_boundary*matrix(rep(p_pron, each=length(recips)), ncol=J) +
                 (1-w_boundary)*matrix(rep(p_det,  each=length(recips)), ncol=J)
  Xsim <- matrix(rbinom(I*J, size=1, prob=as.vector(P)), nrow=I, ncol=J)
  rownames(Xsim) <- items; colnames(Xsim) <- colnames(X_obs)
  Xsim
}

one_run <- function(Xs){
  D  <- jacc(Xs)
  d_each  <- delta_to_fused(D, "each_other", pron_set, fused_set)
  d_onean <- delta_to_fused(D, "one_another", pron_set, fused_set)
  knn_each  <- nearest_of(D, "each_other", k=k)
  knn_onean <- nearest_of(D, "one_another", k=k)
  share_pron_each  <- mean(knn_each  %in% pron_set)
  share_pron_onean <- mean(knn_onean %in% pron_set)
  # supervised P(pronoun) for reciprocals (CV on non-reciprocals)
  nonrec <- !(items %in% recips)
  fit <- glmnet::cv.glmnet(x = Xs[nonrec,], y = y_raw[nonrec], family="binomial", alpha=0, nfolds=10, type.measure="deviance")
  p_rec <- as.numeric(predict(fit, newx=Xs[items %in% recips,], s="lambda.min", type="response"))
  tibble(delta_each=d_each, delta_onean=d_onean,
         share_pron_each=share_pron_each, share_pron_onean=share_pron_onean,
         P_each=p_rec[match("each_other", recips)], P_onean=p_rec[match("one_another", recips)])
}

summarise_long <- function(df){
  df |>
    summarise(across(
      everything(),
      list(mean = mean, q05 = ~quantile(.x, .05), q50 = median, q95 = ~quantile(.x, .95))
    )) |>
    pivot_longer(cols = everything(),
                 names_to = c("metric", ".value"),
                 names_pattern = "^(.*)_(mean|q05|q50|q95)$")
}

percentile_of <- function(x, v) mean(v <= x, na.rm=TRUE)

# ---------- Run the sweep ----------
all_metrics <- list()
all_summaries <- list()
percentile_rows <- list()

for(w in weights){
  message(sprintf("Running w=%.2f ...", w))
  res <- purrr::map_dfr(1:R, ~one_run(simulate_matrix(w)), .id="replicate")
  res$w <- w
  # save per-w metrics
  tag <- sprintf("w%03d", as.integer(round(w*100)))
  write.csv(res, sprintf("out/sim_recovery_metrics_%s.csv", tag), row.names=FALSE)

  # per-w summary
  summ <- res |> select(-replicate) |> group_by(w) |> summarise_long()
  write.csv(summ, sprintf("out/sim_recovery_summary_%s.csv", tag), row.names=FALSE)

  all_metrics[[tag]]  <- res
  all_summaries[[tag]] <- summ

  # percentiles for overlays
  pct <- tibble(
    w = w,
    metric = c("P_each","P_onean","Delta_mean"),
    empirical = c(emp_e, emp_o, emp_delta_mean),
    percentile = c(percentile_of(emp_e, res$P_each),
                   percentile_of(emp_o, res$P_onean),
                   percentile_of(emp_delta_mean, (res$delta_each + res$delta_onean)/2))
  )
  percentile_rows[[tag]] <- pct

  # plots with empirical overlays
  p1 <- ggplot(res, aes(delta_each))  + geom_histogram(bins=40) + ggtitle("Δdistance: each_other")
  p2 <- ggplot(res, aes(delta_onean)) + geom_histogram(bins=40) + ggtitle("Δdistance: one_another")
  p3 <- ggplot(res, aes((delta_each + delta_onean)/2)) +
        geom_histogram(bins=40) + geom_vline(xintercept = emp_delta_mean, linetype = 2) +
        ggtitle("Δdistance: mean of reciprocals (dashed = empirical)")
  p4 <- ggplot(res, aes(share_pron_each))  + geom_histogram(bins=20) + ggtitle("kNN pronoun share: each_other")
  p5 <- ggplot(res, aes(share_pron_onean)) + geom_histogram(bins=20) + ggtitle("kNN pronoun share: one_another")
  p6 <- ggplot(res, aes(P_each))  + geom_histogram(bins=20) + geom_vline(xintercept = emp_e, linetype = 2) +
        ggtitle("P(pronoun): each_other (dashed = empirical ridge)")
  p7 <- ggplot(res, aes(P_onean)) + geom_histogram(bins=20) + geom_vline(xintercept = emp_o, linetype = 2) +
        ggtitle("P(pronoun): one_another (dashed = empirical ridge)")
  png(sprintf("out/sim_recovery_plots_%s.png", tag), width=1600, height=1400, res=170)
  print(p1); print(p2); print(p3); print(p4); print(p5); print(p6); print(p7)
  dev.off()
}

# ---------- Combine and save ----------
metrics_all  <- bind_rows(all_metrics)
summary_all  <- bind_rows(all_summaries)
percentiles  <- bind_rows(percentile_rows)

write.csv(metrics_all,  "out/sim_recovery_metrics_all.csv",  row.names=FALSE)
write.csv(summary_all,  "out/sim_recovery_summary_all.csv",  row.names=FALSE)
write.csv(percentiles,  "out/sim_recovery_percentiles.csv",  row.names=FALSE)

# Console recap
print(percentiles)
