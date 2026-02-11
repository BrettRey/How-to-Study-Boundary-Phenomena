# 03_matched_set_robustness.R — rotate matched subsets with shared permutations
# Canonical subset + K rotations; reuse the same permuted matrices for every subset;
# report per-subset two-sided p around the perm-mean; save a one-line appendix summary.

suppressPackageStartupMessages({
  library(tidyverse); library(Matrix); library(proxy); library(vegan)
})
set.seed(2025)
dir.create("out", showWarnings=FALSE)

# ---------- Load ----------
dat <- read.csv("matrix_clean.csv", stringsAsFactors=FALSE, check.names=FALSE)
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
  if (b %% 200 == 0) saveRDS(list(T_null_canon = T_null_canon[1:b],
                                  T_null_rot = T_null_rot[1:b, , drop = FALSE]),
                             file = sprintf("out/pp_null_partial_b%05d.rds", b))
  setTxtProgressBar(pb, b)
}
close(pb)

canon_p <- p_two_sided(T_obs_canon, T_null_canon)
rot_p   <- lapply(seq_len(K), function(k) p_two_sided(T_obs_rot[k], T_null_rot[, k]))

# ---------- Assemble tidy results ----------
res_canon <- tibble(
  fused = paste(canon_fused, collapse = ";"),
  pron  = paste(canon_pron,  collapse = ";"),
  T_obs = T_obs_canon, T0 = canon_p$T0, p_two_sided = canon_p$p
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
write.csv(all_res, "out/matched_subset_robustness.csv", row.names = FALSE)

# ---------- Manifest for appendix ----------
sink("out/matched_subset_manifest.txt")
cat("Canonical fused-determinatives:\n"); cat(paste0("  - ", canon_fused, collapse = "\n")); cat("\n\n")
cat("Canonical pronouns:\n");           cat(paste0("  - ", canon_pron,  collapse = "\n"));  cat("\n")
sink()

# ---------- Plots ----------
p_scatter <- ggplot(all_res, aes(T_obs, p_two_sided, color = subset)) +
  geom_point(alpha = .7) + geom_hline(yintercept = .05, linetype = 2) +
  labs(x = "Observed Δ (pron − fused)", y = "Two-sided p (around null mean)", title = "Matched-subset robustness")
p_hist <- ggplot(filter(all_res, subset == "rotation"), aes(p_two_sided)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = res_canon$p_two_sided, linetype = 2) +
  labs(x = "p(two-sided)", title = "Distribution across rotations; dashed = canonical")

png("out/matched_subset_robustness.png", width = 1400, height = 600, res = 150)
print(p_scatter); print(p_hist); dev.off()

# ---------- One-line appendix summary ----------
rot <- filter(all_res, subset=="rotation")$p_two_sided

# Canonical p as a percentile of rotation p's
canon_percentile <- ecdf(rot)(res_canon$p_two_sided)
extra <- sprintf("Canonical p percentile among rotations: %.1f%%", 100*canon_percentile)

txt <- sprintf(
  "Canonical two-sided p = %.3f\nRotations: median p = %.3f; q05 = %.3f; q95 = %.3f; Pr(p<.05) = %.3f; Pr(p<.10) = %.3f\n%s",
  res_canon$p_two_sided,
  median(rot), quantile(rot,.05), quantile(rot,.95),
  mean(rot < .05), mean(rot < .10),
  extra
)

writeLines(txt)
writeLines(txt, "out/matched_subset_robustness.txt")

# --- Canonical permutation null plot (Δ under quasiswap) ---
suppressPackageStartupMessages(library(tidyverse))

null_df <- tibble(delta = T_null_canon)
p <- ggplot(null_df, aes(delta)) +
  geom_histogram(bins = 40, linewidth = 0.2) +
  geom_vline(xintercept = mean(T_null_canon), linetype = 3) +   # permutation mean E0[Δ]
  geom_vline(xintercept = T_obs_canon,  linetype = 2) +         # observed Δ
  labs(title = "Permutation null for Δ (quasiswap, B = 5,000)",
       x = expression(Delta~"(under null)"),
       y = "count") +
  theme_minimal(base_size = 12)

ggsave("out/permutation_null.png", p, width = 7, height = 4.5, dpi = 300, bg = "white")
