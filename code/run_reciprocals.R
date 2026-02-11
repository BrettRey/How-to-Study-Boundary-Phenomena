# run_reciprocals.R — one-shot runner (R + cmdstanr)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(readr)
  library(stringr)
  library(posterior)
})

# ---- write Stan model (new array syntax) ----
stan_txt <- '
data {
  int<lower=1> N_items;
  int<lower=1> N_features;
  int<lower=1> N_blocks;
  array[N_features] int<lower=1, upper=N_blocks> block_id;
  array[N_items, N_features] int<lower=0, upper=1> y_obs;
  array[N_items, N_features] int<lower=0, upper=1> is_obs;
  array[N_items, N_features] int<lower=0, upper=1> y_holdout;
  int<lower=0> K_pron_anchors;
  array[K_pron_anchors] int<lower=1, upper=N_items> pron_anchor_idx;
  int<lower=0> K_det_anchors;
  array[K_det_anchors] int<lower=1, upper=N_items> det_anchor_idx;
  real<lower=0> eta_anchor_sd;
}
parameters {
  vector<lower=0, upper=1>[N_items] eta;
  vector[N_blocks] alpha;
  vector[N_blocks] beta;
  vector[N_features] u;
  vector[N_items]    v_raw;
  real<lower=0> tau_f;
  real<lower=0> tau_i;
  vector<lower=0, upper=1>[N_blocks] gamma_pos;
  vector<lower=0, upper=1>[N_blocks] gamma_neg;
}
transformed parameters {
  vector[N_items] v = tau_i * v_raw;
}
model {
  alpha ~ normal(0, 2);
  beta  ~ normal(0, 2);
  tau_f ~ normal(0, 1);
  tau_i ~ normal(0, 1);
  u     ~ normal(0, tau_f);
  v_raw ~ normal(0, 1);
  eta   ~ beta(1, 1);

  for (k in 1:K_pron_anchors)
    target += normal_lpdf(eta[pron_anchor_idx[k]] | 1, eta_anchor_sd);
  for (k in 1:K_det_anchors)
    target += normal_lpdf(eta[det_anchor_idx[k]] | 0, eta_anchor_sd);

  // slightly tighter measurement-error priors
  gamma_pos ~ beta(2, 50);  // mean ~ 0.038
  gamma_neg ~ beta(2, 50);

  for (i in 1:N_items)
    for (j in 1:N_features)
      if (is_obs[i, j] == 1) {
        int b = block_id[j];
        real lin = alpha[b] + beta[b] * eta[i] + u[j] + v[i];
        real p   = inv_logit(lin);                          // P(x=1)
        real q   = (1 - gamma_neg[b]) * p + gamma_pos[b] * (1 - p); // P(y=1)
        y_obs[i, j] ~ bernoulli(q);
      }
}
generated quantities {
  array[N_items, N_features] int y_rep;
  real lpd_mask = 0;
  for (i in 1:N_items)
    for (j in 1:N_features) {
      int b = block_id[j];
      real lin = alpha[b] + beta[b] * eta[i] + u[j] + v[i];
      real p   = inv_logit(lin);
      real q   = (1 - gamma_neg[b]) * p + gamma_pos[b] * (1 - p);
      y_rep[i, j] = bernoulli_rng(q);
      if (is_obs[i, j] == 0)
        lpd_mask += bernoulli_lpmf(y_holdout[i, j] | q);
    }
}
'
writeLines(stan_txt, "reciprocals_hpc.stan")

# ---- load matrix + blocks ----
X_df <- read_csv("matrix_clean.csv", show_col_types = FALSE)
if (!"lemma" %in% names(X_df)) stop("matrix_clean.csv must have a 'lemma' column.")
# drop a gold label if present
if ("class" %in% names(X_df)) X_df <- X_df[, setdiff(names(X_df), "class")]

item_names    <- X_df$lemma
feature_names <- setdiff(names(X_df), "lemma")
X_mat <- as.matrix(X_df[, feature_names, drop = FALSE])
storage.mode(X_mat) <- "integer"
if (!all(X_mat %in% c(0L,1L))) stop("All feature cells must be 0/1.")

fb <- read_csv("feature_blocks_by_order.csv", show_col_types = FALSE)
if (!all(c("feature","block") %in% names(fb))) stop("feature_blocks_by_order.csv needs columns: feature, block")
fb <- fb[ fb$block %in% c("morph","phon","sem","synt"), ]
fb <- fb[ match(feature_names, fb$feature), ]
if (any(is.na(fb$feature))) stop("Block mapping missing for some feature names.")
block_levels <- c("morph","phon","sem","synt")
block_id_vec <- match(fb$block, block_levels)

# ---- mask some cells for predictive scoring ----
set.seed(20250820)
mask   <- matrix(rbinom(nrow(X_mat) * ncol(X_mat), 1, 0.15), nrow(X_mat), ncol(X_mat))
is_obs <- 1 - mask

# ---- anchors (only those that actually exist in your lemma list) ----
canon <- function(x) tolower(trimws(x))
iname <- canon(item_names)

want_pron <- canon(c("i","you_pron_sing","it","he","we","they_plur","she","me","him","us_pron","them_plur"))  
pron_anchor_idx <- match(want_pron, iname); pron_anchor_idx <- pron_anchor_idx[!is.na(pron_anchor_idx)]

want_det <- canon(c("the","a","this","that","an","all","some","no","these","more"))
det_anchor_idx <- match(want_det, iname); det_anchor_idx <- det_anchor_idx[!is.na(det_anchor_idx)]

message("Pronoun anchors: ", if (length(pron_anchor_idx)) paste(item_names[pron_anchor_idx], collapse=", ") else "(none)")
message("Determiner anchors: ", if (length(det_anchor_idx)) paste(item_names[det_anchor_idx], collapse=", ") else "(none)")

# ---- Stan data ----
stan_data <- list(
  N_items         = nrow(X_mat),
  N_features      = ncol(X_mat),
  N_blocks        = length(block_levels),
  block_id        = as.integer(block_id_vec),
  y_obs           = X_mat,
  is_obs          = is_obs,
  y_holdout       = X_mat,
  K_pron_anchors  = length(pron_anchor_idx),
  pron_anchor_idx = as.integer(pron_anchor_idx),
  K_det_anchors   = length(det_anchor_idx),
  det_anchor_idx  = as.integer(det_anchor_idx),
  eta_anchor_sd   = 0.03   # slightly tighter anchors
)

# ---- compile + sample ----
mod <- cmdstan_model("reciprocals_hpc.stan")
fit <- mod$sample(
  data = stan_data,
  seed = 2025,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.9, max_treedepth = 12
)

# ---- results ----
eta_draws <- fit$draws("eta")
lpd_mask  <- fit$draws("lpd_mask")

qfun <- function(x) { qs <- as.numeric(stats::quantile(x, c(0.05,0.50,0.95))); names(qs) <- c("q5","q50","q95"); qs }

recip_names <- c("each_other","one_another")
recip_idx   <- match(recip_names, iname)

cat("\nReciprocals (η on pronoun↔determinative continuum):\n")
if (all(!is.na(recip_idx))) {
  out_rec <- posterior::summarise_draws(eta_draws[, , recip_idx, drop = FALSE], "mean", qfun)
  out_rec <- cbind(item = recip_names, as.data.frame(out_rec))
  print(out_rec, row.names = FALSE)
} else {
  cat("Could not find 'each_other' and/or 'one_another' in lemma.\n")
}

# anchors (quick sanity)
anchor_idx <- c(pron_anchor_idx, det_anchor_idx)
if (length(anchor_idx)) {
  cat("\nAnchors (η summaries):\n")
  out_anc <- posterior::summarise_draws(eta_draws[, , anchor_idx, drop = FALSE], "mean", qfun)
  out_anc <- cbind(item = item_names[anchor_idx], as.data.frame(out_anc))
  print(out_anc, row.names = FALSE)
}

# predictive score
n_masked <- sum(mask)
mean_lpd_per_cell <- mean(as.array(lpd_mask)) / n_masked
cat(sprintf("\nMean masked log predictive density per cell: %.3f  (masked cells = %d)\n",
            mean_lpd_per_cell, n_masked))
