# reciprocals_jaccard_permutation.R
# Run: Rscript reciprocals_jaccard_permutation.R
# Purpose: Primary Jaccard-based permutation test that reciprocals are closer to fused determinatives than to pronouns.
# Notes:
#   - Uses robust column detection: looks for 'lemma' and 'class' (falls back to 'category').
#   - Binarizes all feature columns deterministically; tokens like "true"/"yes"/"1" -> 1; blanks/NA -> 0.
#   - Holds reciprocals fixed; shuffles labels among non-reciprocals while preserving fused/pronoun group sizes.

# -------------------- Configuration --------------------
csv_path      <- "matrix_clean.csv"
set.seed(20250819)
n_perm        <- 10000        # set 1000 for a quick check
report_euclid <- FALSE        # TRUE to add Euclidean robustness check

# Target lexeme sets (match the original Python analysis)
fused_names       <- c("someone","anyone","anything","everything","somebody","anybody")
pronoun_names     <- c("he","him","himself","she","her","herself","they","them","themselves")
reciprocal_names  <- c("each_other","one_another")

# -------------------- Helpers --------------------
key_us <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

normalize_binary <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x > 0))
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    y <- trimws(tolower(x))
    y[y %in% c("", "na", "nan")] <- NA
    out <- rep(NA_integer_, length(y))
    out[y %in% c("1","true","t","yes","y")] <- 1L
    out[y %in% c("0","false","f","no","n")] <- 0L
    suppressWarnings(num <- as.numeric(y))
    out[is.na(out) & !is.na(num)] <- as.integer(num > 0)
    out[is.na(out)] <- 0L
    return(out)
  }
  suppressWarnings(num <- as.numeric(x))
  out <- as.integer(num > 0)
  out[is.na(out)] <- 0L
  out
}

jaccard_dist_vec <- function(x, y) {
  # x, y binary integer vectors
  inter <- sum((x == 1L) & (y == 1L))
  uni   <- sum((x == 1L) | (y == 1L))
  if (uni == 0L) return(0)
  1 - inter / uni
}

jaccard_to_all <- function(v, Mbin) {
  # distances from vector v (binary) to every row of Mbin (binary)
  apply(Mbin, 1, function(r) jaccard_dist_vec(v, r))
}

mean_group_dist <- function(dists, idx) {
  if (length(idx) == 0) return(NA_real_)
  mean(dists[idx])
}

perm_pvalue <- function(stat_obs, stat_null, direction=c("greater","less")) {
  direction <- match.arg(direction)
  m <- length(stat_null)
  b <- if (direction == "greater") sum(stat_null >= stat_obs) else sum(stat_null <= stat_obs)
  p <- (b + 1) / (m + 1)
  se <- sqrt(p * (1 - p) / m)
  ci <- p + c(-1,1) * 1.96 * se
  list(p=p, se=se, ci=ci)
}

# -------------------- Load & prepare data --------------------
DF <- read.csv(csv_path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

# Detect core columns
cols_lower <- tolower(names(DF))
lem_col <- if ("lemma" %in% cols_lower) names(DF)[match("lemma", cols_lower)] else
           if ("lexeme" %in% cols_lower) names(DF)[match("lexeme", cols_lower)] else
           if ("word" %in% cols_lower) names(DF)[match("word", cols_lower)] else
           stop("No lemma/lexeme/word column found. Columns: ", paste(names(DF), collapse=", "))

cat_col <- if ("class" %in% cols_lower) names(DF)[match("class", cols_lower)] else
           if ("category" %in% cols_lower) names(DF)[match("category", cols_lower)] else
           stop("No class/category column found. Columns: ", paste(names(DF), collapse=", "))

# Canonicalise lemma keys
DF$lemma_key <- key_us(DF[[lem_col]])
cat_key      <- tolower(trimws(DF[[cat_col]]))

# Build feature matrix by excluding lemma + class/category
feature_cols <- setdiff(names(DF), c(lem_col, cat_col, "lemma_key"))
feat <- DF[, feature_cols, drop=FALSE]
feat[] <- lapply(feat, normalize_binary)
Xbin <- as.matrix(feat)
storage.mode(Xbin) <- "integer"

# Optional: numeric version for Euclidean robustness
Xnum <- NULL
if (isTRUE(report_euclid)) {
  Xnum <- as.matrix(DF[, feature_cols, drop=FALSE])
  storage.mode(Xnum) <- "double"
  Xnum[is.na(Xnum)] <- 0
}

# -------------------- Index sets --------------------
fused_keys       <- key_us(fused_names)
pronoun_keys     <- key_us(pronoun_names)
reciprocal_keys  <- key_us(reciprocal_names)

idx_fused <- which(DF$lemma_key %in% fused_keys)
idx_pron  <- which(DF$lemma_key %in% pronoun_keys & grepl("^pron", cat_key))
idx_recip <- which(DF$lemma_key %in% reciprocal_keys)

msg <- function(v, lab) paste0(lab, ": ", length(v), " [", paste(DF[[lem_col]][v], collapse=", "), "]")
cat("Found\n",
    msg(idx_fused, "  fused determinatives"), "\n",
    msg(idx_pron,  "  pronouns"), "\n",
    msg(idx_recip, "  reciprocals (held-out test items)"), "\n\n", sep="")

if (length(idx_recip) != 2)
  stop("Expected 2 reciprocals; found ", length(idx_recip), ". Check lemma spellings (", lem_col, ").")
if (length(idx_fused) == 0 || length(idx_pron) == 0)
  stop("Empty fused or pronoun group after matching lemmas. Category column used: ", cat_col, ".")

# Extract matrices
XF <- Xbin[idx_fused, , drop=FALSE]
XP <- Xbin[idx_pron,  , drop=FALSE]
XR <- Xbin[idx_recip, , drop=FALSE]

# -------------------- Precompute reciprocal-to-all Jaccard distances --------------------
recip_to_all <- lapply(1:nrow(XR), function(i) jaccard_to_all(XR[i,], Xbin))

# -------------------- Observed statistics (Jaccard) --------------------
delta_r <- vapply(1:nrow(XR), function(i) {
  dists <- recip_to_all[[i]]
  dP <- mean_group_dist(dists, idx_pron)
  dF <- mean_group_dist(dists, idx_fused)
  dP - dF
}, numeric(1))

delta_obs <- mean(delta_r)          # positive => closer to fused than to pronouns
count_obs <- sum(delta_r > 0)

cat(sprintf("Observed Jaccard results:\n  mean Δ (d_to_pron - d_to_fused): %.4f\n  count closer to fused (out of %d): %d\n\n",
            delta_obs, length(delta_r), count_obs))

# -------------------- Permutation test (shuffle labels among non-reciprocals) --------------------
pool_idx <- setdiff(seq_len(nrow(DF)), idx_recip)
# Preserve group sizes as in the observed data
nF <- length(idx_fused)
nP <- length(idx_pron)

null_delta  <- numeric(n_perm)
null_counts <- integer(n_perm)

for (b in 1:n_perm) {
  perm <- sample(pool_idx, length(pool_idx), replace=FALSE)
  fused_b <- perm[1:nF]
  pron_b  <- perm[(nF+1):(nF+nP)]

  d_r_b <- vapply(1:length(recip_to_all), function(i) {
    dists <- recip_to_all[[i]]
    dP <- mean_group_dist(dists, pron_b)
    dF <- mean_group_dist(dists, fused_b)
    dP - dF
  }, numeric(1))

  null_delta[b]  <- mean(d_r_b)
  null_counts[b] <- sum(d_r_b > 0)
}

p_delta <- perm_pvalue(delta_obs,  null_delta,  "greater")
p_count <- perm_pvalue(count_obs,  null_counts, "greater")

cat(sprintf("Permutation p (Jaccard, mean Δ): %.4f  [SE %.4f, 95%% CI %.4f..%.4f]\n",
            p_delta$p, p_delta$se, p_delta$ci[1], p_delta$ci[2]))
cat(sprintf("Permutation p (Jaccard, count):  %.4f  [SE %.4f, 95%% CI %.4f..%.4f]\n\n",
            p_count$p, p_count$se, p_count$ci[1], p_count$ci[2]))

# -------------------- Optional robustness: Euclidean-to-centroid --------------------
if (isTRUE(report_euclid)) {
  euclid_centroid <- function(M) colMeans(M)
  EF <- euclid_centroid(Xnum[idx_fused, , drop=FALSE])
  EP <- euclid_centroid(Xnum[idx_pron,  , drop=FALSE])
  dist_e <- function(v, c) sqrt(sum((v - c)^2))

  delta_r_e <- vapply(idx_recip, function(i) {
    dP <- dist_e(Xnum[i,], EP)
    dF <- dist_e(Xnum[i,], EF)
    dP - dF
  }, numeric(1))
  delta_obs_e <- mean(delta_r_e)

  null_delta_e <- numeric(n_perm)
  for (b in 1:n_perm) {
    perm <- sample(pool_idx, length(pool_idx), replace=FALSE)
    fused_b <- perm[1:nF]
    pron_b  <- perm[(nF+1):(nF+nP)]
    EF_b <- colMeans(Xnum[fused_b, , drop=FALSE])
    EP_b <- colMeans(Xnum[pron_b,  , drop=FALSE])
    d_r_b <- vapply(idx_recip, function(i) {
      dP <- sqrt(sum((Xnum[i,]-EP_b)^2))
      dF <- sqrt(sum((Xnum[i,]-EF_b)^2))
      dP - dF
    }, numeric(1))
    null_delta_e[b] <- mean(d_r_b)
  }
  p_delta_e <- perm_pvalue(delta_obs_e, null_delta_e, "greater")
  cat(sprintf("Permutation p (Euclidean centroid, mean Δ): %.4f  [SE %.4f, 95%% CI %.4f..%.4f]\n",
              p_delta_e$p, p_delta_e$se, p_delta_e$ci[1], p_delta_e$ci[2]))
}

# -------------------- Summary --------------------
cat("\n=== SUMMARY (Jaccard primary) ===\n")
cat(sprintf("Observed mean Δ: %.4f (positive ⇒ reciprocals closer to fused determinatives)\n", delta_obs))
cat(sprintf("Observed count closer to fused: %d of %d\n", count_obs, length(delta_r)))
cat(sprintf("p (mean Δ): %.4f; SE %.4f; 95%% CI [%.4f, %.4f]\n", p_delta$p, p_delta$se, p_delta$ci[1], p_delta$ci[2]))
cat(sprintf("p (count):  %.4f; SE %.4f; 95%% CI [%.4f, %.4f]\n", p_count$p, p_count$se, p_count$ci[1], p_count$ci[2]))
