# --- eta_strip.R : make figs/eta_strip.pdf (and png) --------------------------
# Inputs expected in the environment:
#   - fit2            : cmdstanr fit with parameter 'eta'
#   - item_names OR iname : character vector of item names (length N_items)
#   - pron_anchor_idx : integer indices of pronoun anchors (optional but recommended)
#   - det_anchor_idx  : integer indices of determinative anchors (optional but recommended)
#
# Output:
#   figs/eta_strip.pdf  (anchor-calibrated η, tall enough for labels)
#   figs/eta_strip.png

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# 0) Names ---------------------------------------------------------------
item_names <- if (exists("item_names")) item_names else if (exists("iname")) iname else NULL
stopifnot(!is.null(item_names))
N_items <- length(item_names)

# 1) Pull η draws as a matrix: rows=draws, cols=eta[1]..eta[N_items] ----
eta_dm <- fit2$draws("eta", format = "draws_matrix")
# Keep only eta[...] columns in correct order
eta_cols <- grep("^eta\\[\\d+\\]$", colnames(eta_dm), value = TRUE)
eta_dm   <- eta_dm[, eta_cols, drop = FALSE]
stopifnot(ncol(eta_dm) == N_items)

# 2) Anchor-calibrate η to [0,1] using anchor means per draw --------------
#    (0 = determinative anchors' mean, 1 = pronoun anchors' mean)
if (!exists("pron_anchor_idx")) pron_anchor_idx <- integer(0)
if (!exists("det_anchor_idx"))  det_anchor_idx  <- integer(0)

anchor_ok <- length(pron_anchor_idx) > 0 && length(det_anchor_idx) > 0

if (anchor_ok) {
  # Per-draw anchor means
  mu_pron <- rowMeans(eta_dm[, pron_anchor_idx, drop = FALSE])
  mu_det  <- rowMeans(eta_dm[, det_anchor_idx,  drop = FALSE])
  denom   <- (mu_pron - mu_det)

  # Guard against degenerate denominator
  # If a draw has near-zero separation, skip calibration for that draw (leave original)
  near_zero <- abs(denom) < 1e-6
  denom[near_zero] <- 1  # avoid division by ~0; those rows won't change much anyway

  # Linear transform per draw: (η - mean_det) / (mean_pron - mean_det), then clamp to [0,1]
  # Do this efficiently in matrix form
  eta_centered <- sweep(eta_dm, 1, mu_det, FUN = "-")        # subtract det mean (per row)
  eta_scaled   <- sweep(eta_centered, 1, denom, FUN = "/")   # divide by separation (per row)
  eta_scaled[eta_scaled < 0] <- 0
  eta_scaled[eta_scaled > 1] <- 1

  eta_for_plot <- eta_scaled
  scale_note   <- " (anchor-calibrated)"
} else {
  # No anchors available: plot raw η on the model’s native scale
  eta_for_plot <- eta_dm
  scale_note   <- " (model scale)"
  warning("No anchor indices found; plotting raw η without anchor calibration.")
}

# 3) Summaries for η[i] --------------------------------------------------
qfun <- function(x) {
  qs <- as.numeric(stats::quantile(x, c(0.05, 0.50, 0.95)))
  names(qs) <- c("q5","q50","q95")
  qs
}
eta_summ <- apply(eta_for_plot, 2, qfun) |>
  t() |>
  as.data.frame() |>
  tibble::as_tibble(.name_repair = "minimal") |>
  mutate(idx = as.integer(str_extract(colnames(eta_dm), "(?<=\\[)\\d+(?=\\])"))) |>
  arrange(idx)

stopifnot(nrow(eta_summ) == N_items)

# 4) Label groups: anchors, reciprocals, other ---------------------------
recips <- c("each_other","one_another")
recip_idx <- match(recips, item_names)

grp <- rep("other", N_items)
if (length(pron_anchor_idx)) grp[pron_anchor_idx] <- "pron_anchor"
if (length(det_anchor_idx))  grp[det_anchor_idx]  <- "det_anchor"
grp[recip_idx[!is.na(recip_idx)]] <- "reciprocal"

plot_df <- eta_summ |>
  transmute(
    item = item_names[idx],
    q5  = q5,
    q50 = q50,
    q95 = q95,
    group = factor(grp, levels = c("det_anchor","reciprocal","other","pron_anchor"))
  ) |>
  arrange(q50) |>
  mutate(item = factor(item, levels = item))

# 5) Plot (dot + 90% interval), highlight reciprocals/anchors -------------
# Unicode arrow via expression to avoid mbcsToSbcs warning; Cairo fixes glyphs in PDF.
ylab_expr <- expression(paste("Posterior ", eta, .(scale_note), " (pronoun \u2194 determinative)"))

p <- ggplot(plot_df, aes(x = item, y = q50, ymin = q5, ymax = q95, color = group, shape = group)) +
  geom_pointrange(size = 0.4, linewidth = 0.4, alpha = 0.95) +
  coord_flip() +
  labs(x = NULL, y = ylab_expr) +
  scale_shape_manual(
    values = c(det_anchor = 16, reciprocal = 17, other = 16, pron_anchor = 16),
    drop = FALSE
  ) +
  scale_color_manual(
    values = c(det_anchor = "#D62728", reciprocal = "#9467BD", other = "gray35", pron_anchor = "#1F77B4"),
    drop = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_size = 10) +
  theme(
    legend.position   = "top",
    legend.title      = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text.y       = element_text(hjust = 1),
    plot.margin       = margin(10, 15, 10, 15)
  )

# 6) Save to figs/ with height scaled to #items ---------------------------
dir.create("figs", showWarnings = FALSE)

# Height rule of thumb: ~0.10 in per item + 2 in top/bottom padding, cap at 18 in.
height_in <- min(18, 0.10 * N_items + 2.0)

use_cairo <- isTRUE(capabilities("cairo"))
pdf_device <- if (use_cairo) grDevices::cairo_pdf else grDevices::pdf

ggsave(filename = "figs/eta_strip.pdf", plot = p,
       width = 7.5, height = height_in, device = pdf_device)

# Optional raster for slides
ggsave(filename = "figs/eta_strip.png", plot = p,
       width = 7.5, height = height_in, dpi = 300)
# ------------------------------------------------------------------------------
