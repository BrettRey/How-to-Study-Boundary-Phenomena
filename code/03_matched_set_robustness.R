# 03_matched_set_robustness.R — rotate matched subsets with shared permutations
# Canonical subset + K rotations; reuse the same permuted matrices for every subset;
# summarize the observed Delta distribution across rotations and save a one-line appendix summary.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(Matrix)
  library(proxy)
  library(vegan)
})
set.seed(2025)

root_dir <- if (file.exists("data/matrix_clean.csv")) "." else if (file.exists("../data/matrix_clean.csv")) ".." else stop("Run from repo root or code/.")
data_path <- function(...) file.path(root_dir, "data", ...)
plot_path <- function(...) file.path(root_dir, "plots", ...)
root_path <- function(...) file.path(root_dir, ...)

dir.create(data_path(), showWarnings = FALSE, recursive = TRUE)
dir.create(plot_path(), showWarnings = FALSE, recursive = TRUE)

house_blue <- "#6AADE4"
house_red <- "#E85D4C"

# ---------- Load ----------
dat <- read.csv(data_path("matrix_clean.csv"), stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("lemma","class") %in% names(dat)))
items <- dat$lemma
y_raw <- factor(dat$class, levels=c("determinative","pronoun"))
X <- as.matrix(dat[, setdiff(names(dat), c("lemma","class"))]); storage.mode(X) <- "numeric"; X[is.na(X)] <- 0
rownames(X) <- items

pron_all <- items[y_raw=="pronoun"]
det_all  <- items[y_raw=="determinative"]
recips <- intersect(c("each_other","one_another"), items)
stopifnot(length(recips)==2)

# fused-det set
fused_col <- names(dat)[grepl("fused", names(dat), ignore.case=TRUE)][1]
is_fused <- if(!is.na(fused_col)) dat[[fused_col]]==1 else rep(FALSE, nrow(dat))
fused_all <- items[y_raw=="determinative" & is_fused]
stopifnot(length(fused_all) >= 6)

# ---------- Helpers ----------
as_mat <- function(D){ M <- as.matrix(D); rn <- attr(D, "Labels"); rownames(M) <- colnames(M) <- rn; M }
mean_delta_for_sets <- function(Mdist, recip_ids, pron_ids, fused_ids){
  delta_item <- function(item_id){
    m_pron <- mean(Mdist[item_id, pron_ids], na.rm = TRUE)
    m_fuse <- mean(Mdist[item_id, fused_ids], na.rm = TRUE)
    m_pron - m_fuse
  }
  mean(vapply(recip_ids, delta_item, numeric(1)))
}
p_two_sided <- function(T_obs, T_null){
  T0 <- mean(T_null)
  p  <- (sum(abs(T_null - T0) >= abs(T_obs - T0)) + 1) / (length(T_null) + 1)
  list(T0 = T0, p = p)
}
upper_tail_summary <- function(T_obs, T_null){
  (sum(T_null >= T_obs) + 1) / (length(T_null) + 1)
}
mcse_binomial <- function(p, n){
  sqrt(p * (1 - p) / n)
}

# ---------- Subset design ----------
n_per_class <- min(6, length(setdiff(fused_all, recips)), length(setdiff(pron_all, recips)))
canon_fused <- head(sort(setdiff(fused_all, recips)), n_per_class)
canon_pron  <- head(sort(setdiff(pron_all,  recips)), n_per_class)

K <- 100
rot_fused <- replicate(K, sample(setdiff(fused_all, recips), n_per_class), simplify = FALSE)
rot_pron  <- replicate(K, sample(setdiff(pron_all,  recips), n_per_class), simplify = FALSE)

# ---------- Restrict to used rows and draw permutations once ----------
all_items <- rownames(X)
idx <- function(v) match(v, all_items)
recip_ids <- idx(recips)
canon_fused_ids <- idx(canon_fused); canon_pron_ids <- idx(canon_pron)
rot_fused_ids  <- lapply(rot_fused, idx); rot_pron_ids <- lapply(rot_pron, idx)

rows_used <- unique(c(recip_ids, canon_fused_ids, canon_pron_ids, unlist(rot_fused_ids), unlist(rot_pron_ids)))
X_used <- X[rows_used, , drop = FALSE]; lab_used <- rownames(X_used)

B <- 5000; burn <- 1500   # adjust up for the paper
message(sprintf("Permuting a %d x %d presence/absence matrix once for B=%d permutations...",
                nrow(X_used), ncol(X_used), B))
perm <- vegan::permatswap(X_used, method = "quasiswap", mtype = "prab", fixedmar = "both", times = B, burnin = burn)

# index maps within X_used
map_id <- function(orig_ids){ match(all_items[orig_ids], lab_used) }
recip_u <- map_id(recip_ids)
canon_fused_u <- map_id(canon_fused_ids); canon_pron_u <- map_id(canon_pron_ids)
rot_fused_u <- lapply(rot_fused_ids, map_id); rot_pron_u <- lapply(rot_pron_ids, map_id)

# ---------- Observed stats on X_used ----------
D_obs <- proxy::dist(X_used, method = "Jaccard"); M_obs <- as_mat(D_obs)
T_obs_canon <- mean_delta_for_sets(M_obs, recip_u, canon_pron_u, canon_fused_u)
T_obs_rot   <- vapply(seq_len(K), function(k) mean_delta_for_sets(M_obs, recip_u, rot_pron_u[[k]], rot_fused_u[[k]]), numeric(1))

# ---------- Null stats for all subsets (share the permutations) ----------
T_null_canon <- numeric(B)
T_null_rot   <- matrix(NA_real_, nrow = B, ncol = K)

pb <- txtProgressBar(min = 0, max = B, style = 3)
for (b in seq_len(B)) {
  Zu <- perm$perm[[b]]
  D_b <- proxy::dist(Zu, method = "Jaccard"); Mb <- as_mat(D_b)
  T_null_canon[b] <- mean_delta_for_sets(Mb, recip_u, canon_pron_u, canon_fused_u)
  for (k in seq_len(K)) {
    T_null_rot[b, k] <- mean_delta_for_sets(Mb, recip_u, rot_pron_u[[k]], rot_fused_u[[k]])
  }
  setTxtProgressBar(pb, b)
}
close(pb)

canon_p <- p_two_sided(T_obs_canon, T_null_canon)
rot_p   <- lapply(seq_len(K), function(k) p_two_sided(T_obs_rot[k], T_null_rot[, k]))
canon_tail <- upper_tail_summary(T_obs_canon, T_null_canon)
canon_tail_mcse <- mcse_binomial(canon_tail, B + 1)

# ---------- Assemble tidy results ----------
res_canon <- tibble(
  fused = paste(canon_fused, collapse = ";"),
  pron  = paste(canon_pron,  collapse = ";"),
  T_obs = T_obs_canon,
  T0 = canon_p$T0,
  p_two_sided = canon_p$p,
  tail_area_upper = canon_tail,
  tail_area_mcse = canon_tail_mcse
)
res_rot <- tibble(
  fused = vapply(rot_fused, paste, character(1), collapse = ";"),
  pron  = vapply(rot_pron,  paste, character(1), collapse = ";"),
  T_obs = T_obs_rot,
  T0    = vapply(rot_p, `[[`, numeric(1), "T0"),
  p_two_sided = vapply(rot_p, `[[`, numeric(1), "p")
)
all_res <- bind_rows(mutate(res_canon, subset="canonical"),
                     mutate(res_rot,   subset="rotation"))
write.csv(all_res, data_path("matched_subset_robustness.csv"), row.names = FALSE)
write.csv(
  tibble(draw = seq_len(B), delta = T_null_canon),
  data_path("quasiswap_reference_delta.csv"),
  row.names = FALSE
)

# ---------- Manifest for appendix ----------
sink(root_path("matched_subset_manifest.txt"))
cat("Canonical fused-determinatives:\n"); cat(paste0("  - ", canon_fused, collapse = "\n")); cat("\n\n")
cat("Canonical pronouns:\n");           cat(paste0("  - ", canon_pron,  collapse = "\n"));  cat("\n")
sink()

# ---------- Plots ----------
rot_effects <- filter(all_res, subset == "rotation")$T_obs
effect_lines <- tibble(
  x = c(quantile(rot_effects, 0.25), median(rot_effects), quantile(rot_effects, 0.75)),
  kind = c("iqr", "median", "iqr")
)

p_hist <- ggplot(tibble(delta = rot_effects), aes(delta)) +
  geom_histogram(bins = 30, fill = house_blue, colour = "white", linewidth = 0.2, alpha = 0.8) +
  geom_vline(data = effect_lines,
             aes(xintercept = x, linetype = kind),
             linewidth = 0.7, colour = "grey20", show.legend = FALSE) +
  geom_vline(xintercept = res_canon$T_obs, linetype = 2, linewidth = 0.9, colour = house_red) +
  scale_linetype_manual(values = c("median" = 1, "iqr" = 3)) +
  labs(
    x = expression(Delta~"(pronoun" - "compound determinative)"),
    y = "count",
    title = "Distribution across rotations; dashed = prespecified set"
  ) +
  theme_minimal(base_size = 12)

ggsave(plot_path("matched_subset_robustness.png"), p_hist, width = 7, height = 4.5, dpi = 300, bg = "white")

# ---------- One-line appendix summary ----------
rot <- filter(all_res, subset=="rotation")$T_obs

# Canonical effect as a percentile of rotation effects
canon_percentile <- ecdf(rot)(res_canon$T_obs)
extra <- sprintf("Canonical Δ percentile among rotations: %.1f%%", 100*canon_percentile)

txt <- sprintf(
  paste(
    "Canonical Δ = %.3f",
    "Canonical upper-tail summary = %.3f (MC SE ≈ %.3f)",
    "Reference model: quasiswap with preserved row/column totals; B = %d; burn-in = %d; seed = 2025",
    "Rotations: median Δ = %.3f; q25 = %.3f; q75 = %.3f; Pr(Δ>0) = %.3f",
    "%s",
    sep = "\n"
  ),
  res_canon$T_obs,
  canon_tail,
  canon_tail_mcse,
  B,
  burn,
  median(rot), quantile(rot,.25), quantile(rot,.75),
  mean(rot > 0),
  extra
)

writeLines(txt)
writeLines(txt, root_path("matched_subset_robustness.txt"))

# --- Canonical reference distribution plot (Δ under quasiswap) ---
null_df <- tibble(delta = T_null_canon)
p <- ggplot(null_df, aes(delta)) +
  geom_histogram(bins = 40, fill = house_blue, colour = "white", linewidth = 0.2, alpha = 0.8) +
  geom_vline(xintercept = mean(T_null_canon), linetype = 3, colour = "grey50", linewidth = 0.5) +   # permutation mean E0[Δ]
  geom_vline(xintercept = T_obs_canon,  linetype = 2, colour = house_red, linewidth = 0.8) +         # observed Δ
  labs(title = "Quasiswap reference distribution for Δ (B = 5,000)",
       x = "Δ under row-/column-preserving reference model",
       y = "count") +
  theme_minimal(base_size = 12)

ggsave(plot_path("permutation_null.png"), p, width = 7, height = 4.5, dpi = 300, bg = "white")
