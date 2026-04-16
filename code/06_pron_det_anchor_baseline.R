#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(glmnet))

root <- getwd()
if (!file.exists(file.path(root, "data", "matrix_clean.csv")) &&
    basename(root) == "code") {
  root <- dirname(root)
}
if (!file.exists(file.path(root, "data", "matrix_clean.csv"))) {
  stop("Run from the repo root or the code/ directory.")
}

source(file.path(root, "code", "helpers.R"))

matrix_path <- file.path(root, "data", "matrix_clean.csv")
inventory_path <- file.path(root, "data", "inventory_annotations.csv")
anchor_sets_path <- file.path(root, "data", "anchor_sets.csv")
pron_path <- file.path(root, "data", "clear_pronoun_anchors.txt")
det_path <- file.path(root, "data", "clear_determinative_anchors.txt")
review_path <- file.path(root, "data", "anchor_review_items.txt")

score_path <- file.path(root, "data", "pron_det_anchor_scores.csv")
summary_path <- file.path(root, "pron_det_anchor_summary.txt")
plot_path <- file.path(root, "plots", "pron_det_anchor_mds.png")
heatmap_path <- file.path(root, "plots", "pron_det_anchor_heatmap.png")
clustermap_path <- file.path(root, "plots", "pron_det_anchor_clustermap.png")
focus_clustermap_path <- file.path(root, "plots", "pron_det_review_focus_clustermap.png")

read_manifest <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

nearest_anchor <- function(i, idx, items, D) {
  candidates <- setdiff(idx, i)
  if (!length(candidates)) {
    return(list(lemma = NA_character_, distance = NA_real_))
  }
  d <- D[i, candidates]
  j <- candidates[which.min(d)]
  list(lemma = items$lemma[j], distance = D[i, j])
}

mean_anchor_distance <- function(i, idx, D) {
  candidates <- setdiff(idx, i)
  if (!length(candidates)) {
    return(NA_real_)
  }
  mean(D[i, candidates])
}

matrix_df <- read.csv(
  matrix_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

inventory <- read.csv(
  inventory_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

anchor_sets <- read.csv(
  anchor_sets_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pron_anchors <- read_manifest(pron_path)
det_anchors <- read_manifest(det_path)
review_items <- read_manifest(review_path)

feature_cols <- setdiff(names(matrix_df), c("lemma", "class"))
feature_df <- matrix_df[, feature_cols, drop = FALSE]
feature_df[] <- lapply(feature_df, normalize_binary)
Xbin <- as.matrix(feature_df)
storage.mode(Xbin) <- "integer"

items <- data.frame(
  lemma = matrix_df$lemma,
  lemma_key = key_us(matrix_df$lemma),
  seed_class = matrix_df$class,
  stringsAsFactors = FALSE
)

inventory <- inventory[
  match(items$lemma, inventory$lemma),
  c(
    "lemma",
    "revised_class",
    "category_source",
    "personhood_profile",
    "anchor_recommendation",
    "notes"
  )
]

anchor_sets <- anchor_sets[
  match(items$lemma, anchor_sets$lemma),
  c("lemma", "anchor_bucket")
]

items$revised_class <- inventory$revised_class
items$category_source <- inventory$category_source
items$personhood_profile <- inventory$personhood_profile
items$anchor_recommendation <- inventory$anchor_recommendation
items$notes <- inventory$notes
items$anchor_bucket <- anchor_sets$anchor_bucket

if (!setequal(pron_anchors, items$lemma[items$anchor_bucket == "clear_pronoun_anchor"])) {
  stop("Pronoun anchor manifest and anchor_sets.csv disagree.")
}
if (!setequal(det_anchors, items$lemma[items$anchor_bucket == "clear_determinative_anchor"])) {
  stop("Determinative anchor manifest and anchor_sets.csv disagree.")
}
if (!setequal(review_items, items$lemma[items$anchor_bucket == "review_item"])) {
  stop("Review manifest and anchor_sets.csv disagree.")
}

pron_idx <- which(items$anchor_bucket == "clear_pronoun_anchor")
det_idx <- which(items$anchor_bucket == "clear_determinative_anchor")
train_idx <- c(pron_idx, det_idx)
y_train <- as.integer(train_idx %in% pron_idx)

nfolds <- max(3L, min(10L, min(sum(y_train == 0L), sum(y_train == 1L))))

set.seed(20260413)
cv_model <- cv.glmnet(
  Xbin[train_idx, , drop = FALSE],
  y_train,
  family = "binomial",
  alpha = 0,
  nfolds = nfolds,
  keep = TRUE,
  type.measure = "class"
)

lambda_idx <- which.min(abs(cv_model$lambda - cv_model$lambda.min))
train_cv_prob <- as.numeric(cv_model$fit.preval[, lambda_idx])
train_cv_pred <- as.integer(train_cv_prob >= 0.5)
cv_accuracy <- mean(train_cv_pred == y_train)
cv_pron_recall <- mean(train_cv_pred[y_train == 1L] == 1L)
cv_det_recall <- mean(train_cv_pred[y_train == 0L] == 0L)

prob_pronoun <- as.numeric(
  predict(cv_model, newx = Xbin, s = "lambda.min", type = "response")
)

D <- calculate_jaccard_matrix(Xbin)
mds <- cmdscale(as.dist(D), k = 2)

mean_pron_dist <- numeric(nrow(Xbin))
mean_det_dist <- numeric(nrow(Xbin))
nearest_pron_lemma <- character(nrow(Xbin))
nearest_det_lemma <- character(nrow(Xbin))
nearest_pron_dist <- numeric(nrow(Xbin))
nearest_det_dist <- numeric(nrow(Xbin))

for (i in seq_len(nrow(Xbin))) {
  pron_nn <- nearest_anchor(i, pron_idx, items, D)
  det_nn <- nearest_anchor(i, det_idx, items, D)
  mean_pron_dist[i] <- mean_anchor_distance(i, pron_idx, D)
  mean_det_dist[i] <- mean_anchor_distance(i, det_idx, D)
  nearest_pron_lemma[i] <- pron_nn$lemma
  nearest_det_lemma[i] <- det_nn$lemma
  nearest_pron_dist[i] <- pron_nn$distance
  nearest_det_dist[i] <- det_nn$distance
}

scores <- items
scores$prob_pronoun <- prob_pronoun
scores$prob_determinative <- 1 - prob_pronoun
scores$predicted_class <- ifelse(prob_pronoun >= 0.5, "pronoun", "determinative")
scores$boundary_uncertainty <- 1 - 2 * abs(prob_pronoun - 0.5)
scores$mean_jaccard_to_pronoun_anchors <- mean_pron_dist
scores$mean_jaccard_to_determinative_anchors <- mean_det_dist
scores$distance_advantage_pronoun <- mean_det_dist - mean_pron_dist
scores$nearest_pronoun_anchor <- nearest_pron_lemma
scores$nearest_pronoun_distance <- nearest_pron_dist
scores$nearest_determinative_anchor <- nearest_det_lemma
scores$nearest_determinative_distance <- nearest_det_dist
scores$mds1 <- mds[, 1]
scores$mds2 <- mds[, 2]

write.csv(scores, score_path, row.names = FALSE)

included_idx <- which(scores$anchor_bucket != "excluded_item")
heatmap_order <- included_idx[order(scores$prob_pronoun[included_idx], decreasing = TRUE)]
D_heatmap <- D[heatmap_order, heatmap_order, drop = FALSE]
rownames(D_heatmap) <- scores$lemma[heatmap_order]
colnames(D_heatmap) <- scores$lemma[heatmap_order]

heatmap_ann <- data.frame(
  Bucket = factor(
    scores$anchor_bucket[heatmap_order],
    levels = c("clear_pronoun_anchor", "review_item", "clear_determinative_anchor")
  ),
  RevisedClass = factor(
    scores$revised_class[heatmap_order],
    levels = c("pronoun", "determinative")
  ),
  row.names = scores$lemma[heatmap_order]
)

bucket_colors <- c(
  clear_pronoun_anchor = "#2F66B0",
  review_item = "#1A1A1A",
  clear_determinative_anchor = "#B84435"
)

heat_cols <- colorRampPalette(c("#0B3C5D", "#F7F7F7", "#8C1D18"))(100)
breaks <- seq(0, 1, length.out = length(heat_cols) + 1)
n_heat <- nrow(D_heatmap)
z_heat <- t(D_heatmap[n_heat:1, , drop = FALSE])

png(heatmap_path, width = 2600, height = 2600, res = 220)
par(mar = c(3, 1, 3, 10), xpd = NA)
image(
  x = seq_len(n_heat),
  y = seq_len(n_heat),
  z = z_heat,
  col = heat_cols,
  breaks = breaks,
  axes = FALSE,
  xlab = "",
  ylab = ""
)
box()
title(main = "Pronoun vs Determinative Anchor Distances")
mtext("Ordered by model P(pronoun); excluded items omitted", side = 1, line = 1.8)

legend_x <- n_heat + max(4, ceiling(0.08 * n_heat))
legend_y <- seq(1, n_heat, length.out = length(heat_cols))
image(
  x = legend_x,
  y = legend_y,
  z = matrix(seq_along(heat_cols), nrow = 1),
  col = heat_cols,
  breaks = seq(0.5, length(heat_cols) + 0.5, by = 1),
  add = TRUE
)
rect(
  xleft = legend_x - 0.4,
  ybottom = min(legend_y),
  xright = legend_x + 0.4,
  ytop = max(legend_y),
  border = "black"
)
text(legend_x + 1.2, min(legend_y), labels = "0", adj = c(0, 0.5), cex = 0.85)
text(legend_x + 1.2, n_heat / 2, labels = "0.5", adj = c(0, 0.5), cex = 0.85)
text(legend_x + 1.2, max(legend_y), labels = "1", adj = c(0, 0.5), cex = 0.85)
text(
  legend_x - 0.3,
  max(legend_y) + max(4, ceiling(0.04 * n_heat)),
  labels = "Jaccard\ndistance",
  adj = c(0, 0),
  cex = 0.85
)
legend(
  x = legend_x - 1.0,
  y = max(legend_y) * 0.82,
  legend = c("Pronoun anchor zone", "Review items mixed in order", "Determinative anchor zone"),
  fill = bucket_colors[c("clear_pronoun_anchor", "review_item", "clear_determinative_anchor")],
  border = NA,
  bty = "n",
  cex = 0.9
)
dev.off()

bucket_side_colors <- bucket_colors[as.character(heatmap_ann$Bucket)]
hc_heatmap <- hclust(as.dist(D_heatmap), method = "average")

png(clustermap_path, width = 2600, height = 2600, res = 220)
par(xpd = NA)
heatmap(
  D_heatmap,
  Rowv = as.dendrogram(hc_heatmap),
  Colv = as.dendrogram(hc_heatmap),
  scale = "none",
  col = heat_cols,
  breaks = breaks,
  symm = TRUE,
  labRow = rep("", n_heat),
  labCol = rep("", n_heat),
  margins = c(6, 6),
  RowSideColors = bucket_side_colors,
  ColSideColors = bucket_side_colors,
  main = "Pronoun vs Determinative Anchor Cluster Map"
)
legend(
  "topright",
  inset = c(-0.22, 0),
  legend = c("Clear pronoun anchor", "Review item", "Clear determinative anchor"),
  fill = bucket_colors[c("clear_pronoun_anchor", "review_item", "clear_determinative_anchor")],
  border = NA,
  bty = "n",
  cex = 0.95
)
dev.off()

focus_review_scores <- scores[scores$anchor_bucket == "review_item", ]
focus_lemmas <- unique(c(
  focus_review_scores$lemma,
  focus_review_scores$nearest_pronoun_anchor,
  focus_review_scores$nearest_determinative_anchor
))
focus_idx <- match(focus_lemmas, scores$lemma)
focus_idx <- focus_idx[!is.na(focus_idx)]
focus_scores <- scores[focus_idx, ]
D_focus <- D[focus_idx, focus_idx, drop = FALSE]
rownames(D_focus) <- focus_scores$lemma
colnames(D_focus) <- focus_scores$lemma
hc_focus <- hclust(as.dist(D_focus), method = "average")
focus_labels <- ifelse(
  focus_scores$anchor_bucket == "review_item",
  paste0("*", focus_scores$lemma),
  focus_scores$lemma
)
focus_side_colors <- c(
  clear_pronoun_anchor = "#2F66B0",
  clear_determinative_anchor = "#B84435",
  review_item = "#1A1A1A"
)[focus_scores$anchor_bucket]

png(focus_clustermap_path, width = 2600, height = 2600, res = 220)
par(xpd = NA)
heatmap(
  D_focus,
  Rowv = as.dendrogram(hc_focus),
  Colv = as.dendrogram(hc_focus),
  scale = "none",
  col = heat_cols,
  breaks = breaks,
  symm = TRUE,
  labRow = focus_labels,
  labCol = focus_labels,
  margins = c(12, 12),
  cexRow = 0.9,
  cexCol = 0.9,
  RowSideColors = focus_side_colors,
  ColSideColors = focus_side_colors,
  main = "Review-Focused Pronoun vs Determinative Cluster Map"
)
legend(
  "topright",
  inset = c(-0.24, 0),
  legend = c(
    "Clear pronoun anchor",
    "Clear determinative anchor",
    "Review item"
  ),
  fill = c("#2F66B0", "#B84435", "#1A1A1A"),
  border = NA,
  bty = "n",
  cex = 0.95
)
text(
  par("usr")[2] * 1.18,
  par("usr")[4] * 0.9,
  labels = "* = review item",
  adj = c(0, 0.5),
  cex = 0.95
)
dev.off()

plot_colors <- c(
  clear_pronoun_anchor = "#2F66B0",
  clear_determinative_anchor = "#B84435",
  review_item = "#1A1A1A",
  excluded_item = "#BDBDBD"
)
plot_pch <- c(
  clear_pronoun_anchor = 19,
  clear_determinative_anchor = 17,
  review_item = 15,
  excluded_item = 1
)
plot_cex <- c(
  clear_pronoun_anchor = 0.9,
  clear_determinative_anchor = 0.9,
  review_item = 1.0,
  excluded_item = 0.8
)

png(plot_path, width = 1800, height = 1400, res = 180)
par(mar = c(5, 5, 3, 8), xpd = NA)
plot(
  scores$mds1,
  scores$mds2,
  col = plot_colors[scores$anchor_bucket],
  pch = plot_pch[scores$anchor_bucket],
  cex = plot_cex[scores$anchor_bucket],
  xlab = "MDS1 (Jaccard distance)",
  ylab = "MDS2 (Jaccard distance)",
  main = "Pronoun vs Determinative Anchor Baseline"
)
review_idx <- which(scores$anchor_bucket == "review_item")
text(
  scores$mds1[review_idx],
  scores$mds2[review_idx],
  labels = scores$lemma[review_idx],
  pos = 3,
  cex = 0.65,
  xpd = TRUE
)
legend(
  x = par("usr")[2] + 0.06 * diff(par("usr")[1:2]),
  y = par("usr")[4],
  legend = c(
    "Clear pronoun anchor",
    "Clear determinative anchor",
    "Review item",
    "Excluded item"
  ),
  col = plot_colors[c(
    "clear_pronoun_anchor",
    "clear_determinative_anchor",
    "review_item",
    "excluded_item"
  )],
  pch = plot_pch[c(
    "clear_pronoun_anchor",
    "clear_determinative_anchor",
    "review_item",
    "excluded_item"
  )],
  pt.cex = plot_cex[c(
    "clear_pronoun_anchor",
    "clear_determinative_anchor",
    "review_item",
    "excluded_item"
  )],
  bty = "o",
  bg = "white",
  xjust = 0,
  yjust = 1
)
dev.off()

review_scores <- scores[scores$anchor_bucket == "review_item", ]
review_scores <- review_scores[order(review_scores$prob_pronoun, decreasing = TRUE), ]
ambiguous_review <- review_scores[
  order(abs(review_scores$prob_pronoun - 0.5)),
  c("lemma", "prob_pronoun", "predicted_class", "boundary_uncertainty")
]

summary_lines <- c(
  "Pronoun vs Determinative Anchor Baseline",
  sprintf("Date: %s", Sys.Date()),
  "",
  sprintf("Clear pronoun anchors: %d", length(pron_idx)),
  sprintf("Clear determinative anchors: %d", length(det_idx)),
  sprintf("Review items scored out of sample: %d", nrow(review_scores)),
  sprintf("Excluded items carried through only for reference: %d", sum(scores$anchor_bucket == "excluded_item")),
  "",
  sprintf("Cross-validated anchor accuracy: %.3f", cv_accuracy),
  sprintf("Cross-validated pronoun recall: %.3f", cv_pron_recall),
  sprintf("Cross-validated determinative recall: %.3f", cv_det_recall),
  sprintf("Ridge lambda.min: %.6f", cv_model$lambda.min),
  "",
  "Review items by P(pronoun):"
)

for (i in seq_len(nrow(review_scores))) {
  row <- review_scores[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "  %-18s p=%.3f  pred=%-13s  near_pron=%s(%.3f)  near_det=%s(%.3f)",
      row$lemma,
      row$prob_pronoun,
      row$predicted_class,
      row$nearest_pronoun_anchor,
      row$nearest_pronoun_distance,
      row$nearest_determinative_anchor,
      row$nearest_determinative_distance
    )
  )
}

summary_lines <- c(summary_lines, "", "Most ambiguous review items:")
top_n <- min(5L, nrow(ambiguous_review))
for (i in seq_len(top_n)) {
  row <- ambiguous_review[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "  %-18s p=%.3f  pred=%-13s  uncertainty=%.3f",
      row$lemma,
      row$prob_pronoun,
      row$predicted_class,
      row$boundary_uncertainty
    )
  )
}

writeLines(summary_lines, summary_path)

cat(sprintf("Wrote %s\n", score_path))
cat(sprintf("Wrote %s\n", summary_path))
cat(sprintf("Wrote %s\n", plot_path))
cat(sprintf("Wrote %s\n", heatmap_path))
cat(sprintf("Wrote %s\n", clustermap_path))
cat(sprintf("Wrote %s\n", focus_clustermap_path))
