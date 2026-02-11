# reciprocals_multiverse.R — end-to-end analysis for reciprocals
# Requires: tidyverse, Matrix, proxy, vegan, glmnet, ggplot2, scales, loo
# Optional: brms, bayesplot  (for PPCs)
# Input files: matrix_clean.csv, feature_blocks_by_order.csv (optional)
# Outputs to: out/

# --------------------------- 0) Packages ---------------------------
req <- c("tidyverse","Matrix","proxy","vegan","glmnet","ggplot2","scales","loo")
opt <- c("brms","bayesplot")
for(p in req) if(!requireNamespace(p, quietly=TRUE)) stop(sprintf("Package %s not installed.", p))
invisible(lapply(c(req, intersect(opt, installed.packages()[,1])), require,
                 character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))
dir.create("out", showWarnings=FALSE)

# --------------------------- 1) Load and prepare ---------------------------
dat <- read.csv("matrix_clean.csv", stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("lemma","class") %in% names(dat)))

items <- dat$lemma
y_raw <- factor(dat$class, levels=c("determinative","pronoun"))
X <- as.matrix(dat[, setdiff(names(dat), c("lemma","class"))])
storage.mode(X) <- "numeric"; X[is.na(X)] <- 0
rownames(X) <- items

# Target & comparison sets
recips <- intersect(c("each_other","one_another"), items)
if(length(recips) != 2) stop("Could not find both reciprocals in lemma column.")

# fused determinatives (try to detect; fall back to empty if not present)
fused_col <- names(dat)[grepl("fused", names(dat), ignore.case=TRUE)][1]
is_fused <- if(!is.na(fused_col)) dat[[fused_col]]==1 else rep(FALSE, nrow(dat))
fused_set <- items[y_raw=="determinative" & is_fused]
pron_set  <- items[y_raw=="pronoun"]

# Anchors (conservative defaults)
anchor_pron <- intersect(c("he","she","it_plain","I","me","you_pron_sing","you_pron_plur","we_pron","they_plur","them_plur","him","her_acc","herself","himself","themselves"),
                         pron_set)
if(length(anchor_pron) < 8) anchor_pron <- head(setdiff(pron_set, recips), 10)
anchor_det  <- intersect(c("the","a","this","that","these","those","some","no","all","more","most","many","few","less","least"),
                         items[y_raw=="determinative"])
if(length(anchor_det) < 8) {
  central_det <- setdiff(items[y_raw=="determinative" & !is_fused], fused_set)
  anchor_det <- head(central_det, 10)
}

# Optional feature blocks (by file name you mentioned)
blocks <- NULL
fb_paths <- c("feature_blocks_by_order.csv","feature_blocks.csv")
fb_path <- fb_paths[file.exists(fb_paths)][1]
if(!is.na(fb_path)) {
  fb <- read.csv(fb_path, stringsAsFactors=FALSE)
  # Expect columns: feature, block  (block in {morph, synt, sem, phon})
  names(fb) <- tolower(names(fb))
  stopifnot("feature" %in% names(fb))
  # try to find the block column robustly
  block_col <- if("block" %in% names(fb)) "block" else {
    cand <- names(fb)[sapply(fb, function(v) all(tolower(unique(v)) %in% c("morph","synt","sem","phon")))]
    if(length(cand)==0) stop("Couldn't locate a block column in the feature-blocks file.")
    cand[1]
  }
  blocks <- fb[[block_col]][match(colnames(X), fb$feature)]
  if(any(is.na(blocks))) warning("Some features missing block labels; ablations will ignore those.")
} else {
  message("feature_blocks_by_order.csv not found; block ablations will be skipped.")
}

# --------------------------- 2) Distances & ordinations ---------------------------
# Weighted Jaccard for binary with IDF-like weights
wjaccard_dist <- function(X, w) {
  n <- nrow(X); labs <- rownames(X)
  D <- matrix(0, n, n)
  for(i in 1:n){
    xi <- X[i,]*w
    for(j in i:n){
      xj <- X[j,]*w
      num <- sum(pmin(xi, xj)); den <- sum(pmax(xi, xj))
      d <- if(den>0) 1 - num/den else 0
      D[i,j] <- D[j,i] <- d
    }
  }
  out <- as.dist(D)
  attr(out, "Labels") <- labs
  out
}
idf <- log(nrow(X) / pmax(1, colSums(X)))

# Hamming via Manhattan on binary matrix (scaled to [0,1])
hamming_dist <- function(X) {
  out <- stats::dist(X, method="manhattan")/ncol(X)
  attr(out, "Labels") <- rownames(X)
  out
}

D_list <- list(
  jaccard = proxy::dist(X, method="Jaccard"),
  dice    = proxy::dist(X, method="Dice"),
  hamming = hamming_dist(X),
  wj_idf  = wjaccard_dist(X, idf)
)

# Ordination plot helper (now extracts labels from cmdscale result)
ord_plot <- function(D, title="PCoA"){
  pts <- cmdscale(D, k=2, eig=TRUE)
  lab <- rownames(pts$points)
  coords <- tibble(lemma=lab, Dim1=pts$points[,1], Dim2=pts$points[,2],
                   Class=as.character(y_raw[match(lab, items)]),
                   isRecip=lab %in% recips)
  ggplot(coords, aes(Dim1, Dim2, shape=isRecip, color=Class)) +
    geom_point(size=3, alpha=.85) +
    scale_shape_manual(values=c(`FALSE`=16,`TRUE`=17)) +
    labs(title=title, x="PCo 1", y="PCo 2") +
    theme_minimal()
}
g1 <- ord_plot(D_list$jaccard, "PCoA (Jaccard)")
ggsave("out/pcoa_jaccard.png", g1, width=6.5, height=5.0, dpi=300, bg="white")

# --------------------------- 3) k-NN & Δdistance diagnostics ---------------------------
as_mat <- function(D) { M <- as.matrix(D); rownames(M) <- colnames(M) <- attr(D,"Labels"); M }
mean_dist_to <- function(D, src, tgt) {
  M <- as_mat(D)
  mean(M[src, tgt, drop=FALSE])
}
delta_to_fused <- function(D, recip, pron_targets=pron_set, fused_targets=fused_set) {
  m_pron <- mean_dist_to(D, recip, pron_targets)
  m_fuse <- mean_dist_to(D, recip, fused_targets)
  m_pron - m_fuse  # >0 => closer to fused determinatives
}
nearest_of <- function(D, item, k=5) {
  M <- as_mat(D); d <- M[item, ]
  setdiff(names(sort(d)), item)[1:k]
}
knn_each  <- nearest_of(D_list$jaccard, "each_other", k=5)
knn_onean <- nearest_of(D_list$jaccard, "one_another", k=5)
writeLines(c("kNN(each_other):", paste("  ", knn_each),
             "kNN(one_another):", paste("  ", knn_onean)))

delta_tbl <- tibble(metric=names(D_list),
                    each_other = sapply(D_list, delta_to_fused, recip="each_other"),
                    one_another = sapply(D_list, delta_to_fused, recip="one_another"))
write.csv(delta_tbl, "out/delta_distance_by_metric.csv", row.names=FALSE)

# --------------------------- 4) Supervised classification + calibration ---------------------------
nonrec_idx <- which(!(items %in% recips))
recip_idx  <- which(items %in% recips)
X_nr <- X[nonrec_idx, , drop=FALSE]; y_nr <- y_raw[nonrec_idx]

K <- 10; set.seed(2025)
fold_id <- sample(rep(1:K, length.out=length(nonrec_idx)))
alphas <- c(ridge=0, elastic05=0.5)
calib_out <- list(); recip_pred <- list(); ll_pairs <- list()

for(a_nm in names(alphas)){
  a <- alphas[[a_nm]]
  # Out-of-fold probs for calibration
  oof <- rep(NA_real_, length(nonrec_idx))
  for(k in 1:K){
    tr <- which(fold_id!=k); te <- which(fold_id==k)
    fit <- glmnet::cv.glmnet(x = X_nr[tr,], y = y_nr[tr], family="binomial", alpha=a, nfolds=5, type.measure="deviance")
    p  <- as.numeric(predict(fit, newx=X_nr[te,], s="lambda.min", type="response"))
    oof[te] <- p
  }
  calib_df <- tibble(p=oof, y=as.integer(y_nr=="pronoun")) |>
    mutate(bin=cut(p, breaks=seq(0,1,by=.1), include.lowest=TRUE)) |>
    group_by(bin) |> summarise(pred=mean(p), obs=mean(y), n=dplyr::n(), .groups="drop")
  calib_out[[a_nm]] <- calib_df

  # Fit on all non-reciprocals; predict reciprocals
  fit_all <- glmnet::cv.glmnet(x = X_nr, y = y_nr, family="binomial", alpha=a, nfolds=10, type.measure="deviance")
  p_rec  <- as.numeric(predict(fit_all, newx=X[recip_idx,], s="lambda.min", type="response"))
  recip_pred[[a_nm]] <- tibble(lemma=items[recip_idx], P_pron=p_rec, alpha=a)

  # Bayes factor from recip log-likelihoods under alternative labels
  ll_pron <- sum(log(pmax(1e-12, p_rec)))   # labelling reciprocals as pronouns
  ll_det  <- sum(log(pmax(1e-12, 1-p_rec))) # labelling reciprocals as determinatives
  BF_det_pron <- exp(ll_det - ll_pron)
  ll_pairs[[a_nm]] <- c(ll_pron=ll_pron, ll_det=ll_det, BF_det_pron=BF_det_pron)
}

# Out-of-fold reliability with equal-frequency bins and Brier decomposition
calib_df_full <- tibble(p=oof, y=as.integer(y_nr=="pronoun"))
nbin <- 10
qs <- quantile(calib_df_full$p, probs=seq(0,1,length.out=nbin+1), na.rm=TRUE)
qs[1] <- 0; qs[length(qs)] <- 1
calib_bins <- calib_df_full |>
  mutate(bin = cut(p, breaks=unique(qs), include.lowest=TRUE, dig.lab=3)) |>
  group_by(bin) |>
  summarise(pred=mean(p), obs=mean(y), n=dplyr::n(), .groups="drop")

brier <- mean((calib_df_full$p - calib_df_full$y)^2)
refinement <- var(calib_df_full$p)
calibration_term <- brier - (mean(calib_df_full$y)*(1-mean(calib_df_full$y)) - refinement)
writeLines(sprintf("Brier=%.3f; Calibration=%.3f; Refinement=%.3f", brier, calibration_term, refinement))

calib_plot <- ggplot(calib_bins, aes(pred, obs, size=n)) +
  geom_abline(linetype=2) + geom_point() +
  scale_size_continuous(range=c(2.5,6)) +
  labs(x="Predicted P(pronoun)", y="Observed frequency",
       title="Reliability diagram (equal-frequency bins)") +
  theme_minimal()
ggsave("out/calibration.png", calib_plot, width=7, height=4.5, dpi=300, bg="white")


recip_preds <- bind_rows(recip_pred)
write.csv(recip_preds, "out/recip_pred_probs.csv", row.names=FALSE)
write.csv(as.data.frame(do.call(rbind, ll_pairs)), "out/recip_ll_bayesfactor.csv")

# Posterior from prior via Bayes factor (Gelman-style odds update)
odds_update <- function(prior_p=.85, BF_det_pron){
  prior_odds <- (1-prior_p)/prior_p   # odds for determinative vs pronoun
  post_odds  <- prior_odds * BF_det_pron
  post_p_pron <- 1/(1+post_odds)
  c(posterior_pronoun=post_p_pron)
}
BF_used <- as.data.frame(do.call(rbind, ll_pairs))$BF_det_pron[1]
post85  <- odds_update(.85, BF_used)
writeLines(sprintf("Bayes factor (det:pron) = %.3f; posterior P(pronoun | prior .85) = %.3f",
                   BF_used, post85))

# --------------------------- 5) Specification curve ---------------------------
drop_block <- function(X, block_name){
  keep <- ifelse(is.na(blocks), TRUE, blocks!=block_name)
  X[, keep, drop=FALSE]
}
grid <- expand.grid(dist_metric=c("jaccard","dice","hamming","wj_idf"),
                    alpha_name=names(alphas),
                    block_drop=c("none", if(is.null(blocks)) character(0) else c("morph","synt","sem","phon")),
                    stringsAsFactors=FALSE)

metric_get <- function(metric, X){
  if(metric=="wj_idf") return(wjaccard_dist(X, log(nrow(X) / pmax(1, colSums(X)))))
  if(metric=="jaccard") return(proxy::dist(X, method="Jaccard"))
  if(metric=="dice")    return(proxy::dist(X, method="Dice"))
  if(metric=="hamming") return(hamming_dist(X))
  stop("Unknown metric")
}

spec_row <- function(row){
  Xm <- if(is.null(blocks) || row$block_drop=="none") X else drop_block(X, row$block_drop)
  Dm <- metric_get(row$dist_metric, Xm)
  d_each  <- delta_to_fused(Dm, "each_other")
  d_onean <- delta_to_fused(Dm, "one_another")
  # classifier on reduced feature set (non-reciprocal items only)
  X_nr_m <- Xm[nonrec_idx, , drop=FALSE]; y_nr_m <- y_raw[nonrec_idx]
  fit <- glmnet::cv.glmnet(X_nr_m, y_nr_m, family="binomial", alpha=alphas[[row$alpha_name]], nfolds=10, type.measure="deviance")
  p_rec <- as.numeric(predict(fit, newx=Xm[recip_idx,], s="lambda.min", type="response"))
  ll_pron <- sum(log(pmax(1e-12, p_rec))); ll_det <- sum(log(pmax(1e-12, 1-p_rec)))
  BF <- exp(ll_det - ll_pron)
  tibble(dist=row$dist_metric, alpha=row$alpha_name, block=row$block_drop,
         delta_each=d_each, delta_onean=d_onean,
         P_each=p_rec[match("each_other", items[recip_idx])],
         P_onean=p_rec[match("one_another", items[recip_idx])],
         BF_det_pron=BF)
}
spec_res <- purrr::map_dfr(seq_len(nrow(grid)), ~spec_row(grid[.x, , drop=FALSE]))
write.csv(spec_res, "out/spec_curve.csv", row.names=FALSE)
spec_plot <- spec_res |>
  mutate(idx=row_number(), mean_delta=(delta_each+delta_onean)/2) |>
  ggplot(aes(idx, mean_delta, shape=alpha, color=dist)) + geom_hline(yintercept=0, linetype=2) +
  geom_point() + facet_wrap(~block, scales="free_x") +
  labs(y="Δdistance (pron − fused) ; >0 => closer to fused", x="specification index", title="Specification curve") +
  theme_minimal()
ggsave("out/spec_curve.png", spec_plot, width=8, height=4.8, dpi=300, bg="white")

# --------------------------- 6) Margin-preserving permutation (swap) ---------------------------
set.seed(2025)

# choose a balanced matched subset (fix the seed and the selection for reproducibility)
n_per_class <- min(6, length(setdiff(fused_set, recips)), length(setdiff(pron_set, recips)))
if(n_per_class < 3) warning("Small matched subset; permutation will be very low power.")
fused_n <- head(setdiff(fused_set, recips), n_per_class)
pron_n  <- head(setdiff(pron_set,  recips), n_per_class)
subset_items <- c(recips, fused_n, pron_n)
Xs <- X[subset_items, , drop=FALSE]

stat_delta_subset <- function(M, recips, pron_targets, fused_targets) {
  D <- proxy::dist(M, method="Jaccard")
  Mmat <- as.matrix(D); rownames(Mmat) <- colnames(Mmat) <- attr(D,"Labels")
  delta_item <- function(item) {
    m_pron <- mean(Mmat[item, pron_targets, drop=FALSE])
    m_fuse <- mean(Mmat[item, fused_targets, drop=FALSE])
    m_pron - m_fuse   # >0 => closer to fused determinatives
  }
  mean(sapply(recips, delta_item))
}

# observed statistic
T_obs <- stat_delta_subset(Xs, recips, pron_n, fused_n)

# permutation null via quasiswap (preserve both margins)
B <- 50000  # bump for stability in the paper
perm <- vegan::permatswap(Xs, method="quasiswap", mtype="prab", fixedmar="both",
                          times=B, burnin=2000)
T_null <- vapply(perm$perm, function(M) stat_delta_subset(M, recips, pron_n, fused_n), numeric(1))
T0 <- mean(T_null)  # reference point

# two-sided p-value around T0
p_two <- (sum(abs(T_null - T0) >= abs(T_obs - T0)) + 1) / (B + 1)
perm_out <- tibble(observed=T_obs, center=T0, p_two_sided=p_two, B=B)
write.csv(perm_out, "out/permutation_swap.csv", row.names=FALSE)
saveRDS(list(null=T_null, obs=T_obs, center=T0, subset=list(recips=recips, pron=pron_n, fused=fused_n)),
        "out/permutation_null.rds")

# --------------------------- 7) Design analysis (Type S/M with coherent decision rule) ---------------------------
design_analysis <- function(T_obs, T_null, T0, alpha=.05, R=5000){
  # repeat-experiment estimates under an effect equal to T_obs
  noise <- sample(T_null, R, replace=TRUE) - T0
  That <- T_obs + noise
  # same decision rule as the permutation test
  thr <- quantile(abs(T_null - T0), 1 - alpha)
  sig <- abs(That - T0) >= thr
  typeS <- if(any(sig)) mean(sign(That[sig]) != sign(T_obs)) else NA_real_
  typeM <- if(any(sig)) mean(abs(That[sig]) / abs(T_obs)) else NA_real_
  list(significant = (abs(T_obs - T0) >= thr),
       typeS = typeS, typeM = typeM, thr = thr, alpha = alpha)
}
da <- design_analysis(T_obs, T_null, T0, alpha=.05, R=5000)
writeLines(sprintf("Permutation test p(two-sided)=%.3f; significant=%s; Type-S≈%.3f; Type-M≈%s",
                   p_two, da$significant, da$typeS,
                   ifelse(is.na(da$typeM),"NA", sprintf("%.2f×", da$typeM))))


# --------------------------- 8) Optional: Generative PPCs ---------------------------
if(requireNamespace("brms", quietly=TRUE) && !is.null(blocks)){
  long <- as_tibble(X, rownames="lemma") |>
    pivot_longer(-lemma, names_to="feature", values_to="x") |>
    mutate(class = y_raw[match(lemma, items)],
           block = blocks[match(feature, colnames(X))])
  f1 <- brms::bf(x ~ 1 + class + (1|lemma) + (1|feature) + (1|block), family=bernoulli())
  fit1 <- brms::brm(f1, data=long, iter=2000, chains=4, cores=4, seed=2025, refresh=0)
  f0 <- brms::bf(x ~ 1 + (1|lemma) + (1|feature) + (1|block), family=bernoulli())
  fit0 <- brms::brm(f0, data=long, iter=2000, chains=4, cores=4, seed=2025, refresh=0)
  comp <- loo::loo_compare(loo::loo(fit1), loo::loo(fit0))
  saveRDS(list(fit1=fit1, fit0=fit0, loo=comp), file="out/brms_models.rds")
  # quick PPC: mask 15% cells
  set.seed(2025)
  mask_idx <- sample(nrow(long), size=round(0.15*nrow(long)))
  long_mask <- long; long_mask$x[mask_idx] <- NA
  fit_mask <- brms::brm(x ~ 1 + class + (1|lemma) + (1|feature) + (1|block), family=bernoulli(),
                        data=long_mask, iter=2000, chains=4, cores=4, seed=2025, refresh=0)
  yrep <- posterior_predict(fit_mask, newdata=long_mask, draws=200)
  obs_mask <- long$x[mask_idx]
  p_hat <- colMeans(yrep[, mask_idx, drop=FALSE])
  elpd <- mean(obs_mask*log(pmax(1e-6,p_hat)) + (1-obs_mask)*log(pmax(1e-6,1-p_hat)))
  writeLines(sprintf("Masked-cell mean log predictive density (approx): %.3f", elpd))
}

# --------------------------- 9) Session info snapshot ---------------------------
sink("out/sessionInfo.txt"); print(sessionInfo()); sink()
