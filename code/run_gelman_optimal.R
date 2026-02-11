suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(loo)
  library(vegan)
})

# ---- paths ----
matrix_path <- "matrix_clean.csv"
blocks_path <- "feature_blocks_by_order.csv"
stan_file   <- "reciprocals_hpc_feature_misclass.stan"

set.seed(2025)

# ---- load matrix robustly: first column is names; blanks -> 0; TRUE/FALSE -> 1/0 ----
Xraw <- read_csv(matrix_path, na = c("", "NA", "NaN", "."), show_col_types = FALSE)

# identify name column
name_col <- if ("lemma" %in% names(Xraw)) "lemma" else if ("item" %in% names(Xraw)) "item" else names(Xraw)[1]
item_names <- as.character(Xraw[[name_col]])
X <- Xraw %>% select(-all_of(name_col))

# normalize to 0/1 integers
to01 <- function(x){
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x))
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("TRUE","T","true","t","Yes","yes","Y","y")] <- "1"
  x[x %in% c("FALSE","F","false","f","No","no","N","n","")] <- "0"
  as.integer(x)
}
X <- X %>% mutate(across(everything(), to01))
X[is.na(X)] <- 0L

# sanity: binary check
bad_cols <- names(X)[!vapply(X, function(col) all(col %in% c(0L,1L)), logical(1))]
if (length(bad_cols)) stop(sprintf("Non-binary values in columns: %s", paste(bad_cols, collapse=", ")))

X <- as.matrix(X); storage.mode(X) <- "integer"
N_items <- nrow(X); N_features <- ncol(X)

# ---- blocks ----
blk <- read_csv(blocks_path, show_col_types = FALSE, col_names = "block")
if (nrow(blk) != N_features) stop(sprintf("Block file rows (%d) != #features (%d)", nrow(blk), N_features))
block_labels <- tolower(blk$block)
ulab <- unique(block_labels); block_map <- setNames(seq_along(ulab), ulab)
block_id <- as.integer(block_map[block_labels]); N_blocks <- length(ulab)
message(sprintf("Blocks: %s", paste(sprintf("%s=%d", names(block_map), block_map), collapse=", ")))

# ---- anchors ----
pron_anchor_candidates <- c("I","you_pron_sing","it_plain","he","we","they_plur","she","me","him","us_pron","them_plur")
det_anchor_candidates  <- c("the","a","an","this","that","these","those","some","no","all","more")

norm <- function(s) tolower(gsub("[^a-z_]", "", s))
iname_norm <- norm(item_names)
match_any <- function(cands, pool) { as.integer(unique(stats::na.omit(match(cands, pool)))) }

pron_anchor_idx <- match_any(norm(pron_anchor_candidates), iname_norm)
det_anchor_idx  <- match_any(norm(det_anchor_candidates),  iname_norm)

message(sprintf("Pronoun anchors matched: %s", paste(item_names[pron_anchor_idx], collapse=", ")))
message(sprintf("Determiner anchors matched: %s", paste(item_names[det_anchor_idx], collapse=", ")))
if (length(pron_anchor_idx) < 4 || length(det_anchor_idx) < 4) {
  stop("Anchor match too thin. Tell me the exact lemmas in your CSV for pronoun and determiner anchors.")
}

# ---- anchor priors (logit scale on eta_raw) ----
mu_eta_raw <- rep(0, N_items); sd_eta_raw <- rep(5, N_items)
mu_eta_raw[pron_anchor_idx] <-  2.5; sd_eta_raw[pron_anchor_idx] <- 0.5
mu_eta_raw[det_anchor_idx]  <- -2.5; sd_eta_raw[det_anchor_idx]  <- 0.5

# ---- reciprocals ----
recip_names <- c("each_other","one_another")
recip_idx   <- match(norm(recip_names), iname_norm)

# ---- 15% mask and indices (pass to Stan) ----
n_mask <- round(0.15 * N_items * N_features)
lin_idx <- sample.int(N_items * N_features, n_mask, replace = FALSE)
idx_i <- as.integer(((lin_idx - 1) %% N_items) + 1)
idx_j <- as.integer(((lin_idx - 1) %/% N_items) + 1)
mask <- matrix(0L, N_items, N_features)
mask[cbind(idx_i, idx_j)] <- 1L

stan_data <- list(
  N_items    = N_items,
  N_features = N_features,
  N_blocks   = N_blocks,
  block_id   = as.array(block_id),
  y          = X,
  mask       = mask,
  mu_eta_raw = as.vector(mu_eta_raw),
  sd_eta_raw = as.vector(sd_eta_raw),
  N_mask     = length(idx_i),
  idx_i      = idx_i,
  idx_j      = idx_j
)

# ---- compile & sample ----
mod <- cmdstan_model(stan_file, stanc_options = list("O1"))
fit <- mod$sample(
  data = stan_data,
  seed = 2025,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.9, max_treedepth = 12,
  refresh = 200
)

# ---- eta summaries ----
eta_draws <- fit$draws("eta")
qfun <- function(x) { qs <- as.numeric(stats::quantile(x, c(0.05,0.50,0.95))); names(qs) <- c("q5","q50","q95"); qs }
summ_eta <- posterior::summarise_draws(eta_draws, "mean", qfun) %>%
  as_tibble() %>% mutate(item = item_names) %>% select(item, mean, q5, q50, q95)

cat("\nReciprocals (eta):\n")
if (all(!is.na(recip_idx))) print(summ_eta[recip_idx, ], row.names = FALSE) else cat("No match for reciprocals.\n")

cat("\nAnchors (quick sanity):\n")
print(summ_eta[c(pron_anchor_idx, det_anchor_idx), ] %>% arrange(desc(mean)), row.names = FALSE)

# ---- Δelpd vs simple baselines ----
ll_mask <- fit$draws("log_lik_mask", format = "draws_matrix")  # draws x N_mask
log_mean_exp <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
lppd_main <- apply(ll_mask, 2, log_mean_exp)

eps <- 1e-6
p_j <- pmin(pmax(colMeans(X), eps), 1 - eps)
p_i <- pmin(pmax(rowMeans(X), eps), 1 - eps)
y_mask_vec <- X[cbind(idx_i, idx_j)]
lpd_baseA <- ifelse(y_mask_vec == 1L, log(p_j[idx_j]), log(1 - p_j[idx_j]))
lpd_baseB <- ifelse(y_mask_vec == 1L, log(p_i[idx_i]), log(1 - p_i[idx_i]))

deltaA <- lppd_main - lpd_baseA
deltaB <- lppd_main - lpd_baseB
cat(sprintf("\nΔelpd per masked cell (Main − Baseline A[feature]): mean=%.3f, SE=%.3f\n",
            mean(deltaA), sd(deltaA)/sqrt(length(deltaA))))
cat(sprintf("Δelpd per masked cell (Main − Baseline B[item]):   mean=%.3f, SE=%.3f\n",
            mean(deltaB), sd(deltaB)/sqrt(length(deltaB))))

# ---- PPC summaries (row/col/block) ----
row_exp <- fit$draws("row_sum_exp", format = "draws_matrix")
col_exp <- fit$draws("col_sum_exp", format = "draws_matrix")
row_obs <- rowSums(X); col_obs <- colSums(X)

row_ppc <- tibble(
  item = item_names,
  obs  = row_obs,
  mean = colMeans(row_exp),
  q5   = apply(row_exp, 2, function(v) quantile(v, 0.05)),
  q95  = apply(row_exp, 2, function(v) quantile(v, 0.95))
)
col_ppc <- tibble(
  feature = colnames(X),
  obs  = col_obs,
  mean = colMeans(col_exp),
  q5   = apply(col_exp, 2, function(v) quantile(v, 0.05)),
  q95  = apply(col_exp, 2, function(v) quantile(v, 0.95))
)
cat("\nRow-sum PPC (first 10 rows):\n"); print(head(row_ppc, 10), row.names = FALSE)
cat("\nCol-sum PPC (first 10 cols):\n"); print(head(col_ppc, 10), row.names = FALSE)

# ---- NN overlap & ordination geometry ----
j_obs <- vegdist(X, method = "jaccard", binary = TRUE)
pcoa_obs <- cmdscale(j_obs, k = 2, eig = FALSE)
nn_obs_k <- 5
nn_obs <- apply(as.matrix(j_obs), 1, function(row) order(row, decreasing = FALSE)[2:(nn_obs_k+1)])

pars <- c("alpha_b","beta_b","u_feat","v_item","eta","gplus","gminus")
drw  <- fit$draws(pars, format = "draws_df") %>% dplyr::slice_sample(n = 200)

make_q <- function(dd) {
  getv <- function(prefix, n) as.numeric(dd[, stringr::str_detect(names(dd), paste0("^", prefix, "\\["))][1:n])
  alpha <- getv("alpha_b", N_blocks); beta <- getv("beta_b", N_blocks)
  u     <- getv("u_feat",  N_features); v  <- getv("v_item",  N_items)
  eta   <- getv("eta",     N_items)
  gp    <- getv("gplus",   N_features); gm <- getv("gminus", N_features)
  Q <- matrix(NA_real_, N_items, N_features)
  for (j in 1:N_features) {
    b <- block_id[j]
    for (i in 1:N_items) {
      lp <- alpha[b] + beta[b] * eta[i] + u[j] + v[i]
      p  <- 1/(1+exp(-lp))
      Q[i, j] <- (1 - gm[j]) * p + gp[j] * (1 - p)
    }
  }
  Q
}

nn_overlap <- c(); proc_cor <- c()
for (r in 1:nrow(drw)) {
  Q <- make_q(drw[r, , drop = FALSE])
  Xrep <- (Q > 0.5) * 1L
  j_rep <- vegdist(Xrep, method = "jaccard", binary = TRUE)
  nn_rep <- apply(as.matrix(j_rep), 1, function(row) order(row, decreasing = FALSE)[2:(nn_obs_k+1)])
  overlap <- mean(sapply(1:N_items, function(i) length(intersect(nn_obs[, i], nn_rep[, i]))/nn_obs_k))
  nn_overlap <- c(nn_overlap, overlap)

  pcoa_rep <- cmdscale(j_rep, k = 2, eig = FALSE)
  pr <- vegan::procrustes(pcoa_obs, pcoa_rep, scale = TRUE)
  proc_cor <- c(proc_cor, stats::cor(c(pr$Yrot), c(pr$X)))
}
cat(sprintf("\nNN@5 overlap (y_rep vs y): mean=%.3f, q5=%.3f, q95=%.3f\n",
            mean(nn_overlap), quantile(nn_overlap, 0.05), quantile(nn_overlap, 0.95)))
cat(sprintf("PCoA Procrustes corr:       mean=%.3f, q5=%.3f, q95=%.3f\n",
            mean(proc_cor), quantile(proc_cor, 0.05), quantile(proc_cor, 0.95)))

cat("\nDone.\n")
