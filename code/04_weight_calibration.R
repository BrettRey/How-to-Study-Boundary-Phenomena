# 04_weight_calibration.R — predictive-fit curves across mixture weights
# Inputs:  data/matrix_clean.csv, optionally clear-anchor manifests;
#          falls back to data/matched_subset_robustness.csv if needed
# Outputs: plots/weight_calibration.png, data/weight_calibration_curve.csv,
#          weight_calibration.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
})

root_dir <- if (file.exists("data/matrix_clean.csv")) "." else if (file.exists("../data/matrix_clean.csv")) ".." else stop("Run from repo root or code/.")
data_path <- function(...) file.path(root_dir, "data", ...)
plot_path <- function(...) file.path(root_dir, "plots", ...)
root_path <- function(...) file.path(root_dir, ...)

dir.create(data_path(), showWarnings = FALSE, recursive = TRUE)
dir.create(plot_path(), showWarnings = FALSE, recursive = TRUE)

house_blue <- "#6AADE4"
house_red <- "#E85D4C"

dat <- read.csv(data_path("matrix_clean.csv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("lemma", "class") %in% names(dat)))

recips <- c("each_other", "one_another")

items <- dat$lemma
X <- as.matrix(dat[, setdiff(names(dat), c("lemma", "class"))])
storage.mode(X) <- "numeric"
X[is.na(X)] <- 0
rownames(X) <- items

stopifnot(all(recips %in% rownames(X)))

read_anchor_file <- function(path) {
  if (!file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  trimws(lines[nzchar(trimws(lines))])
}

pron_anchors <- read_anchor_file(data_path("clear_pronoun_anchors.txt"))
det_anchors  <- read_anchor_file(data_path("clear_determinative_anchors.txt"))

if (length(pron_anchors) == 0 || length(det_anchors) == 0) {
  matched <- read.csv(data_path("matched_subset_robustness.csv"), stringsAsFactors = FALSE)
  canon <- matched |> filter(subset == "canonical") |> slice(1)
  if (nrow(canon) != 1) stop("Could not locate the canonical matched subset.")
  split_items <- function(x) strsplit(x, ";", fixed = TRUE)[[1]]
  pron_anchors <- split_items(canon$pron)
  det_anchors  <- split_items(canon$fused)
}

stopifnot(all(pron_anchors %in% rownames(X)))
stopifnot(all(det_anchors %in% rownames(X)))

# Jeffreys smoothing avoids zero/one feature probabilities.
smoothed_rate <- function(M) {
  (colSums(M) + 0.5) / (nrow(M) + 1)
}

p_pron <- smoothed_rate(X[pron_anchors, , drop = FALSE])
p_det  <- smoothed_rate(X[det_anchors,  , drop = FALSE])

log_loss <- function(x, p) {
  p <- pmin(pmax(p, 1e-9), 1 - 1e-9)
  -mean(x * log(p) + (1 - x) * log1p(-p))
}

w_grid <- seq(0, 1, by = 0.001)

curve_df <- map_dfr(recips, function(recip) {
  x <- X[recip, ]
  tibble(
    lemma = recip,
    w = w_grid,
    mean_log_loss = vapply(w_grid, function(w) log_loss(x, w * p_pron + (1 - w) * p_det), numeric(1))
  )
})

best_df <- curve_df |>
  group_by(lemma) |>
  slice_min(order_by = mean_log_loss, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(label = c("each other", "one another")[match(lemma, recips)])

write.csv(curve_df, data_path("weight_calibration_curve.csv"), row.names = FALSE)

summary_lines <- c(
  "Predictive blend minima (mean log loss; lower is better):",
  sprintf("  each_other   : w = %.3f, mean log loss = %.4f",
          best_df$w[best_df$lemma == "each_other"],
          best_df$mean_log_loss[best_df$lemma == "each_other"]),
  sprintf("  one_another  : w = %.3f, mean log loss = %.4f",
          best_df$w[best_df$lemma == "one_another"],
          best_df$mean_log_loss[best_df$lemma == "one_another"])
)
writeLines(summary_lines, root_path("weight_calibration.txt"))

lemma_labels <- c("each_other" = "each other", "one_another" = "one another")
curve_plot <- curve_df |>
  mutate(label = recode(lemma, !!!lemma_labels)) |>
  ggplot(aes(w, mean_log_loss, colour = label)) +
  geom_line(linewidth = 1) +
  geom_vline(data = best_df, aes(xintercept = w, colour = label), linetype = 2, linewidth = 0.8, show.legend = FALSE) +
  scale_colour_manual(values = c("each other" = house_blue, "one another" = house_red)) +
  labs(
    title = "Predictive fit across mixture weights",
    x = "Mixture weight w (pronoun share)",
    y = "Mean log loss (lower is better)",
    colour = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(plot_path("weight_calibration.png"), curve_plot, width = 7, height = 4.8, dpi = 300, bg = "white")

print(best_df)
