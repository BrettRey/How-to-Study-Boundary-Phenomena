# 02_ppc_cmdstan.R — posterior predictive checks targeted to reported summaries
# PPCs: (a) Procrustes correlation of PCoA geometries; (b) k-NN overlap for anchors, for reciprocals,
# and their difference. Saves CSV summaries and three PNGs.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(vegan)
  library(cmdstanr)
  library(posterior)
})

set.seed(2025)
dir.create("out", showWarnings = FALSE)

# ---------- CmdStan availability ----------
cmdstan_ok <- tryCatch({
  v <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  if (inherits(v, "try-error") || is.null(v)) {
    p <- try(cmdstanr::cmdstan_path(), silent = TRUE)
    !(inherits(p, "try-error") || is.null(p)) && dir.exists(p)
  } else TRUE
}, error = function(e) FALSE)
if (!cmdstan_ok) stop("CmdStan not available. Run: cmdstanr::install_cmdstan()")

# ---------- Load ----------
dat <- read.csv("matrix_clean.csv", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("lemma","class") %in% names(dat)))
items <- dat$lemma
y_lab <- factor(dat$class, levels = c("determinative","pronoun"))

X <- as.matrix(dat[, setdiff(names(dat), c("lemma","class"))])
storage.mode(X) <- "numeric"; X[is.na(X)] <- 0
rownames(X) <- items
I <- nrow(X); J <- ncol(X)

# blocks (optional)
blocks <- rep(1L, J)
if (file.exists("feature_blocks_by_order.csv")) {
  fb <- read.csv("feature_blocks_by_order.csv", stringsAsFactors = FALSE)
  names(fb) <- tolower(names(fb))
  bl <- factor(fb$block[match(colnames(X), fb$feature)], levels = c("morph","synt","sem","phon"))
  blocks <- as.integer(bl)
  if (anyNA(blocks)) {
    mx <- max(blocks, na.rm = TRUE)
    blocks[is.na(blocks)] <- mx + 1L
  }
}
B <- as.integer(max(blocks))

# ---------- Long format for Stan ----------
grid <- expand.grid(i = seq_len(I), j = seq_len(J))
y <- as.integer(as.vector(X))
item <- as.integer(grid$i)
feat <- as.integer(grid$j)
block <- as.integer(blocks[feat])
class_i <- as.integer(y_lab == "pronoun")

# ---------- Stan model ----------
stan_code <- "
data{
  int<lower=1> N; int<lower=1> I; int<lower=1> J; int<lower=1> B;
  array[N] int<lower=1, upper=I> item;
  array[N] int<lower=1, upper=J> feat;
  array[N] int<lower=1, upper=B> block;
  array[N] int<lower=0, upper=1> y;
  array[I] int<lower=0, upper=1> class_i; // 1=pronoun
}
parameters{
  real alpha; real beta_class;
  vector[I] u_item; vector[J] v_feat; vector[B] w_block;
  real<lower=0> sigma_item; real<lower=0> sigma_feat; real<lower=0> sigma_block;
}
model{
  alpha ~ normal(0, 1.5);
  beta_class ~ normal(0, 1);
  u_item ~ normal(0, sigma_item);
  v_feat ~ normal(0, sigma_feat);
  w_block ~ normal(0, sigma_block);
  sigma_item ~ exponential(1);
  sigma_feat ~ exponential(1);
  sigma_block ~ exponential(1);
  for(n in 1:N){
    real eta = alpha + beta_class * class_i[item[n]] + u_item[item[n]] + v_feat[feat[n]] + w_block[block[n]];
    y[n] ~ bernoulli_logit(eta);
  }
}
generated quantities{
  array[N] int y_rep;        // one posterior predictive replicate
  for(n in 1:N){
    real eta = alpha + beta_class * class_i[item[n]] + u_item[item[n]] + v_feat[feat[n]] + w_block[block[n]];
    y_rep[n] = bernoulli_logit_rng(eta);
  }
}
"
tmp <- write_stan_file(stan_code)

data_list <- list(
  N = length(y), I = I, J = J, B = B,
  item = item, feat = feat, block = block, y = y,
  class_i = as.array(class_i)
)

mod <- cmdstan_model(tmp)
fit <- mod$sample(data = data_list, seed = 2025,
                  chains = 4, parallel_chains = 4,
                  iter_warmup = 500, iter_sampling = 500, refresh = 0)
sink("out/stan_fit_summary.txt")
print(fit$summary(c("alpha","beta_class","sigma_item","sigma_feat","sigma_block")))
sink()

# ---------- Helpers ----------
# Safe Jaccard for binary matrices; d(empty,empty)=0
jacc_safedist <- function(M){
  M <- (M > 0L) * 1L
  n <- nrow(M); D <- matrix(0, n, n)
  for(i in seq_len(n-1)){
    xi <- M[i, ]
    for(j in (i+1):n){
      xj <- M[j, ]
      a <- sum(xi & xj); b <- sum(xi & (!xj)); c <- sum((!xi) & xj)
      denom <- a + b + c
      dij <- if (denom == 0) 0 else 1 - a/denom
      D[i, j] <- D[j, i] <- dij
    }
  }
  rownames(D) <- rownames(M); colnames(D) <- rownames(M)
  stats::as.dist(D)
}
as_mat <- function(D){ M <- as.matrix(D); lab <- attr(D, "Labels"); if (!is.null(lab)) rownames(M) <- colnames(M) <- lab; M }
nearest_of <- function(D, it, k = 5){
  M <- as_mat(D); d <- M[it, ]; d[it] <- Inf
  res <- names(sort(d, decreasing = FALSE)); res <- res[res != it]
  res[seq_len(min(k, length(res)))]
}
pcoa2 <- function(D, k = 2){
  fit <- cmdscale(D, k = k, add = TRUE, eig = FALSE)
  if (is.list(fit)) fit$points else fit
}

# Robust Procrustes correlation:
# - use symmetric scaling;
# - clamp 1 - ss into [0,1] to avoid tiny negative roundoff;
# - on any error or degeneracy, fall back to protest()'s t0; otherwise NA.
procrust_R <- function(D_obs, D_rep, k = 2){
  Mo <- as_mat(D_obs); Mr <- as_mat(D_rep)
  keep <- rowSums(!is.finite(Mo)) == 0 & rowSums(!is.finite(Mr)) == 0
  if (sum(keep) < (k + 1)) return(NA_real_)
  Xo <- pcoa2(stats::as.dist(Mo[keep, keep]), k = k)
  Xr <- pcoa2(stats::as.dist(Mr[keep, keep]), k = k)
  if (!is.matrix(Xo) || !is.matrix(Xr) || anyNA(Xo) || anyNA(Xr)) return(NA_real_)
  out <- tryCatch({
    pr <- vegan::procrustes(Xo, Xr, scale = TRUE, symmetric = TRUE)
    R <- sqrt(pmin(pmax(1 - pr$ss, 0), 1))
    as.numeric(R)
  }, error = function(e) {
    t0 <- try(vegan::protest(Xo, Xr, permutations = 0)$t0, silent = TRUE)
    if (inherits(t0, "try-error") || !is.finite(t0)) NA_real_ else as.numeric(t0)
  })
  out
}

# ---------- Anchor sets ----------
recips <- intersect(c("each_other","one_another"), rownames(X))
pron_set <- rownames(X)[y_lab=="pronoun"]
det_set  <- rownames(X)[y_lab=="determinative"]
anchor_pron <- intersect(c("he","she","it_plain","I","me","you_pron_sing","you_pron_plur","we_pron","they_plur","them_plur","him","her_acc","herself","himself","themselves"),
                         pron_set)
if (length(anchor_pron) < 8) anchor_pron <- head(setdiff(pron_set, recips), 10)
anchor_det  <- intersect(c("the","a","this","that","these","those","some","no","all","more","most","many","few","less","least"),
                         det_set)
if (length(anchor_det) < 8) anchor_det <- head(det_set, 10)

# ---------- Observed dissimilarities ----------
D_obs <- jacc_safedist(X)

# ---------- PPC extraction ----------
draws <- fit$draws("y_rep", format = "matrix")
S <- min(200, nrow(draws))  # number of PPC replicates to score

ppc_procr        <- numeric(S)
ppc_knn_anchor   <- numeric(S)
ppc_knn_recip    <- numeric(S)

for (s in seq_len(S)){
  yrep <- matrix(draws[s, ], nrow = I, ncol = J, byrow = FALSE,
                 dimnames = list(rownames(X), colnames(X)))
  D_rep <- jacc_safedist(yrep)

  # (a) Procrustes correlation (robust)
  ppc_procr[s] <- procrust_R(D_obs, D_rep, k = 2)

  # (b) k-NN overlap (Jaccard over neighbour sets)
  k <- 5
  targets <- c(anchor_pron, anchor_det, recips)
  knn_obs <- setNames(lapply(targets, function(it) nearest_of(D_obs, it, k)), targets)
  knn_rep <- setNames(lapply(targets, function(it) nearest_of(D_rep, it, k)), targets)
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

  anchor_ids <- c(anchor_pron, anchor_det)
  ppc_knn_anchor[s] <- mean(mapply(jaccard, knn_obs[anchor_ids], knn_rep[anchor_ids]), na.rm = TRUE)
  ppc_knn_recip[s]  <- if (length(recips)) mean(mapply(jaccard, knn_obs[recips],     knn_rep[recips]),     na.rm = TRUE) else NA_real_
}

# ---------- Save CSVs ----------
write.csv(tibble(draw = 1:S, procrustes_R = ppc_procr),
          "out/ppc_procrustes.csv", row.names = FALSE)
write.csv(tibble(draw = 1:S,
                 knn_overlap_anchors = ppc_knn_anchor,
                 knn_overlap_reciprocals = ppc_knn_recip),
          "out/ppc_knn_overlap.csv", row.names = FALSE)

# ---------- PPC summaries (for appendix) ----------
ppc_sum <- tibble(
  stat = c("procrustes_R", "knn_overlap_anchors", "knn_overlap_reciprocals", "knn_overlap_diff"),
  q05  = c(quantile(ppc_procr, .05, na.rm=TRUE),
           quantile(ppc_knn_anchor, .05, na.rm=TRUE),
           quantile(ppc_knn_recip,  .05, na.rm=TRUE),
           quantile(ppc_knn_anchor - ppc_knn_recip, .05, na.rm=TRUE)),
  q50  = c(median(ppc_procr, na.rm=TRUE),
           median(ppc_knn_anchor, na.rm=TRUE),
           median(ppc_knn_recip,  na.rm=TRUE),
           median(ppc_knn_anchor - ppc_knn_recip, na.rm=TRUE)),
  q95  = c(quantile(ppc_procr, .95, na.rm=TRUE),
           quantile(ppc_knn_anchor, .95, na.rm=TRUE),
           quantile(ppc_knn_recip,  .95, na.rm=TRUE),
           quantile(ppc_knn_anchor - ppc_knn_recip, .95, na.rm=TRUE))
)
write.csv(ppc_sum, "out/ppc_summary.csv", row.names = FALSE)

pr_anchor_gt_recip <- mean(ppc_knn_anchor > ppc_knn_recip, na.rm = TRUE)
writeLines(sprintf("PPC: Pr(anchors overlap > reciprocals) = %.3f", pr_anchor_gt_recip),
           "out/ppc_summary.txt")

# ---------- Plots ----------
png("out/ppc_procrustes_hist.png", width = 900, height = 600, res = 140)
pp <- ppc_procr[is.finite(ppc_procr)]
if (length(pp) >= 5) {
  hist(pp, breaks = 30, main = "PPC: Procrustes correlation (PCoA Jaccard)", xlab = "R")
  abline(v = median(pp), lty = 2)
} else {
  plot.new(); title(main = "PPC: Procrustes correlation (PCoA Jaccard)")
  mtext("No finite values produced; check inputs.", side = 3, line = 0.5)
}
dev.off()

png("out/ppc_knn_overlap_hist.png", width = 900, height = 600, res = 140)
aa <- ppc_knn_anchor[is.finite(ppc_knn_anchor)]
hist(aa, breaks = 30, main = "PPC: k-NN overlap (anchors)", xlab = "Jaccard")
abline(v = median(aa), lty = 2)
dev.off()

png("out/ppc_knn_diff_hist.png", width = 900, height = 600, res = 140)
dd <- (ppc_knn_anchor - ppc_knn_recip)[is.finite(ppc_knn_anchor - ppc_knn_recip)]
hist(dd, breaks = 30, main = "PPC: anchor minus reciprocal k-NN overlap", xlab = "Jaccard difference")
abline(v = median(dd), lty = 2)
dev.off()
