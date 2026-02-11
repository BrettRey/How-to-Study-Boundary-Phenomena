# run_reciprocals_multiverse.R  — FULL REPLACEMENT
# Fixes: adds y_obs and sd_anchor to stan_data (required by your Stan models).
# Also keeps block labels tolerant + template writer for permanent fix.

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(stringr)
  library(vegan)
})

# ---------------------- paths ----------------------
path_matrix <- "matrix_clean.csv"
path_blocks <- "feature_blocks_by_order.csv"

# Stan model files (adjust paths if yours differ)
stan_m0 <- "m0_nomisclass.stan"    # no misclass; per-block alpha,beta; uses anchors
stan_m1 <- "m1_common_slope.stan"  # variant; still declares sd_anchor in data section
stan_m2 <- "m2_full.stan"          # full model: misclass + anchors

set.seed(2025)

# ---------------------- load matrix ----------------------
raw <- suppressMessages(readr::read_csv(path_matrix, show_col_types = FALSE))

to01 <- function(z) {
  if (is.logical(z)) return(as.integer(z))
  if (is.numeric(z)) return(as.integer(z %in% c(1, 0)))
  if (is.character(z)) return(as.integer(tolower(z) %in% c("1","true","t","yes","y")))
  as.integer(rep(0L, length(z)))
}

is_binary_like <- function(x) {
  all(na.omit(as.character(x)) %in%
        c("0","1","TRUE","FALSE","T","F","true","false","Yes","No","yes","no"))
}

first_col <- raw[[1]]
if (!is_binary_like(first_col)) {
  item_names <- as.character(first_col)
  X <- raw[,-1, drop = FALSE]
} else {
  item_names <- paste0("item_", seq_len(nrow(raw)))
  X <- raw
}

X <- X %>% mutate(across(everything(), to01))
X[is.na(X)] <- 0L
X <- as.matrix(X) %>% unname()
storage.mode(X) <- "integer"

N_items    <- nrow(X)
N_features <- ncol(X)

feat_names <- colnames(raw)
if (!is_binary_like(first_col)) feat_names <- feat_names[-1]
if (is.null(feat_names)) feat_names <- paste0("f", seq_len(N_features))

cat(sprintf("Loaded matrix: %d items × %d features\n", N_items, N_features))

# ---------------------- blocks ----------------------
read_blocks <- function(p) {
  df <- tryCatch(read.csv(p, stringsAsFactors = FALSE, header = TRUE),
                 error = function(e) NULL)
  if (is.null(df)) return(NULL)
  if (ncol(df) == 1) return(tibble(block = tolower(trimws(df[[1]]))))
  pick <- intersect(names(df), c("block","label","group","Block","Label","Group",
                                 "BLOCK","LABEL","GROUP"))
  if (length(pick) == 0) tibble(block = tolower(trimws(df[[1]])))
  else tibble(block = tolower(trimws(df[[pick[1]]])))
}

blk <- read_blocks(path_blocks)

if (is.null(blk)) {
  warning("Could not read feature_blocks_by_order.csv. Defaulting ALL features to 'morph'.")
  blk <- tibble(block = rep("morph", N_features))
} else if (nrow(blk) != N_features) {
  tmpl <- tibble(index = seq_len(N_features), feature = feat_names,
                 block  = if (nrow(blk) == N_features) blk$block else NA_character_)
  out_tmpl <- "feature_blocks_TEMPLATE.csv"
  suppressWarnings(readr::write_csv(tmpl, out_tmpl))
  cat("\n[NOTE] Wrote correction template:", out_tmpl,
      "\n      Fill 'block' with morph/phon/sem/synt; save as",
      basename(path_blocks), "to make permanent.\n\n")
  if (nrow(blk) < N_features) {
    pad_lab <- if (nrow(blk) > 0) blk$block[nrow(blk)] else "morph"
    blk <- tibble(block = c(blk$block, rep(pad_lab, N_features - nrow(blk))))
    cat(sprintf("[AUTO-FIX] Padded %d label(s) with '%s'.\n",
                N_features - nrow(blk), pad_lab))
  } else if (nrow(blk) > N_features) {
    blk <- tibble(block = blk$block[seq_len(N_features)])
    cat(sprintf("[AUTO-FIX] Trimmed labels from %d to %d.\n",
                nrow(blk), N_features))
  }
}

map_label <- function(s) {
  if (str_starts(s, "morph")) return("morph")
  if (str_starts(s, "phon"))  return("phon")
  if (str_starts(s, "sem"))   return("sem")
  if (str_starts(s, "synt"))  return("synt")
  return(NA_character_)
}
lab_std <- vapply(blk$block, map_label, character(1))
if (anyNA(lab_std)) {
  bad <- unique(blk$block[is.na(lab_std)])
  stop("Unknown block labels (expect morph/phon/sem/synt prefixes): ", paste(bad, collapse=", "))
}

block_map <- c(morph = 1L, phon = 2L, sem = 3L, synt = 4L)
block_id  <- unname(as.integer(block_map[lab_std]))
N_blocks  <- 4L

cat("Block label summary:\n"); print(table(lab_std))
cat("Block ID range/check: min=", min(block_id), " max=", max(block_id), "\n\n", sep="")

# ---------------------- anchors ----------------------
pron_anchors <- c("I","you_pron_sing","he","she","me","him","us_pron",
                  "they_plur","them_plur","it_plain")
det_anchors  <- c("the","a","this","that","all","some","no","these","those","more")

iname <- item_names
pron_anchor_idx <- match(pron_anchors, iname) %>% discard(is.na)
det_anchor_idx  <- match(det_anchors,  iname) %>% discard(is.na)

cat("Anchors matched:\n")
cat("  Pronoun:      ", paste(iname[pron_anchor_idx], collapse=", "), "\n", sep="")
cat("  Determinative:", paste(iname[det_anchor_idx],  collapse=", "), "\n\n", sep="")

# ---------------------- mask ----------------------
mask_frac <- 0.15
mask <- matrix(0L, N_items, N_features)
N_mask <- as.integer(round(mask_frac * N_items * N_features))
choices <- sample.int(N_items * N_features, N_mask, replace = FALSE)
idx_i <- ((choices - 1L) %% N_items) + 1L
idx_j <- ((choices - 1L) %/% N_items) + 1L
mask[cbind(idx_i, idx_j)] <- 1L

# ---------------------- stan data ----------------------
# IMPORTANT: your Stan models expect y_obs and sd_anchor in the data.
anchor_sd <- 1.0  # adjust if you want tighter/wider anchor priors

stan_data <- list(
  N_items    = N_items,
  N_features = N_features,
  N_blocks   = N_blocks,
  block_id   = block_id,
  # both names provided for safety (some models use y, some y_obs)
  y          = X,
  y_obs      = X,
  mask       = mask,
  N_mask     = N_mask,
  idx_i      = as.array(idx_i),
  idx_j      = as.array(idx_j),
  N_pron_anchor   = length(pron_anchor_idx),
  pron_anchor_idx = if (length(pron_anchor_idx)) as.array(pron_anchor_idx) else array(0L, 0),
  N_det_anchor    = length(det_anchor_idx),
  det_anchor_idx  = if (length(det_anchor_idx))  as.array(det_anchor_idx)  else array(0L, 0),
  sd_anchor       = anchor_sd
)

# ---------------------- compile & sample ----------------------
m0 <- cmdstan_model(stan_m0)
m1 <- cmdstan_model(stan_m1)
m2 <- cmdstan_model(stan_m2)

sample_args <- list(
  seed = 2025,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.9, max_treedepth = 12,
  refresh = 200
)
sample_with <- function(model, data, args) do.call(model$sample, c(list(data = data), args))

cat("Fitting m0 (no misclassification)…\n")
fit0 <- sample_with(m0, stan_data, sample_args)

cat("Fitting m1 (no anchors / common slope)…\n")
fit1 <- sample_with(m1, stan_data, sample_args)

cat("Fitting m2 (full model: misclassification + anchors)…\n")
fit2 <- sample_with(m2, stan_data, sample_args)

# ---------------------- results (m2) ----------------------
recips <- c("each_other","one_another")
recip_idx <- match(recips, iname)

eta2 <- fit2$draws("eta")
qfun <- function(x) { qs <- as.numeric(quantile(x, c(0.05,0.50,0.95))); names(qs) <- c("q5","q50","q95"); qs }

cat("\nReciprocals (η, m2):\n")
if (all(!is.na(recip_idx))) {
  out <- posterior::summarise_draws(eta2[ , , recip_idx, drop = FALSE], "mean", qfun) |>
    as_tibble() |> mutate(item = recips) |> select(item, mean, q5, q50, q95)
  print(out, n = nrow(out), row.names = FALSE)
} else {
  cat("Could not find 'each_other' and/or 'one_another' in item names.\n")
}

anc_idx <- c(pron_anchor_idx, det_anchor_idx)
if (length(anc_idx)) {
  cat("\nAnchors (η, m2) quick check:\n")
  outa <- posterior::summarise_draws(eta2[ , , anc_idx, drop = FALSE], "mean", qfun) |>
    as_tibble() |> mutate(item = iname[anc_idx]) |> select(item, mean, q5, q50, q95) |>
    arrange(desc(mean))
  print(outa, n = nrow(outa), row.names = FALSE)
}

# ---------------------- Δelpd vs baselines (m2) ----------------------
ll2 <- fit2$draws("log_lik_mask", format = "draws_matrix")
log_mean_exp <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
lppd2 <- apply(ll2, 2, log_mean_exp)

eps <- 1e-6
p_feat <- pmin(pmax(colMeans(X), eps), 1 - eps)
p_item <- pmin(pmax(rowMeans(X), eps), 1 - eps)
y_mask <- X[cbind(idx_i, idx_j)]

lpd_feat <- ifelse(y_mask == 1L, log(p_feat[idx_j]), log(1 - p_feat[idx_j]))
lpd_item <- ifelse(y_mask == 1L, log(p_item[idx_i]), log(1 - p_item[idx_i]))

deltaA <- lppd2 - lpd_feat
deltaB <- lppd2 - lpd_item
cat(sprintf("\nΔelpd per masked cell (m2 − per-feature): mean=%.3f, SE=%.3f\n",
            mean(deltaA), sd(deltaA)/sqrt(length(deltaA))))
cat(sprintf("Δelpd per masked cell (m2 − per-item):   mean=%.3f, SE=%.3f\n",
            mean(deltaB), sd(deltaB)/sqrt(length(deltaB))))

# ---------------------- ordination geometry (m2) ----------------------
pars <- c("alpha_b","beta_b","u_feat","v_item","eta","gplus","gminus")
if (all(pars %in% fit2$metadata()$model_params)) {
  drw <- fit2$draws(pars, format = "draws_df") %>% dplyr::slice_sample(n = 100)
  make_q <- function(dd) {
    getv <- function(prefix, n) as.numeric(dd[, stringr::str_detect(names(dd), paste0("^", prefix, "\\["))][1:n])
    alpha <- getv("alpha_b", 4); beta <- getv("beta_b", 4)
    u     <- getv("u_feat",  N_features); v  <- getv("v_item", N_items)
    eta   <- getv("eta",     N_items)
    gp    <- if ("gplus[1]" %in% names(dd))  getv("gplus",  N_features) else rep(0, N_features)
    gm    <- if ("gminus[1]" %in% names(dd)) getv("gminus", N_features) else rep(0, N_features)
    Q <- matrix(NA_real_, N_items, N_features)
    for (j in 1:N_features) {
      b <- block_id[j]
      for (i in 1:N_items) {
        lp <- alpha[b] + beta[b]*eta[i] + u[j] + v[i]
        p  <- 1/(1+exp(-lp))
        Q[i,j] <- (1 - gm[j]) * p + gp[j] * (1 - p)
      }
    }
    Q
  }
  j_obs <- vegdist(X, method = "jaccard", binary = TRUE)
  pcoa_obs <- cmdscale(j_obs, k = 2, eig = FALSE)
  nn_k <- 5
  nn_obs <- apply(as.matrix(j_obs), 1, function(row) order(row)[2:(nn_k+1)])

  nn_overlap <- c(); proc_cor <- c()
  for (r in 1:nrow(drw)) {
    Q <- make_q(drw[r, , drop = FALSE])
    Xrep <- (Q > 0.5) * 1L
    j_rep <- vegdist(Xrep, method = "jaccard", binary = TRUE)
    nn_rep <- apply(as.matrix(j_rep), 1, function(row) order(row)[2:(nn_k+1)])
    ov <- mean(sapply(1:N_items, function(i) length(intersect(nn_obs[, i], nn_rep[, i]))/nn_k))
    nn_overlap <- c(nn_overlap, ov)
    pcoa_rep <- cmdscale(j_rep, k = 2, eig = FALSE)
    pr <- vegan::procrustes(pcoa_obs, pcoa_rep, scale = TRUE)
    proc_cor <- c(proc_cor, stats::cor(c(pr$Yrot), c(pr$X)))
  }
  cat(sprintf("\nNN@5 overlap (rep vs obs): mean=%.3f, q5=%.3f, q95=%.3f\n",
              mean(nn_overlap), quantile(nn_overlap, 0.05), quantile(nn_overlap, 0.95)))
  cat(sprintf("PCoA Procrustes corr:     mean=%.3f, q5=%.3f, q95=%.3f\n",
              mean(proc_cor), quantile(proc_cor, 0.05), quantile(proc_cor, 0.95)))
} else {
  cat("\n[Note] Skipping geometry check: not all parameters found in m2 output.\n")
}

cat("\nDone.\n")
