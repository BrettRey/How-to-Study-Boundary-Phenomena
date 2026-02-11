# ==============================================================================
# File: eta_strip.R
# Goal: Make a clean, publication-ready strip figure of posterior eta by item
# Output: figs/eta_strip.pdf  (vector, good for LaTeX) and figs/eta_strip.png
#
# Requirements in workspace:
#   - fit2        : a cmdstanr fit containing parameter "eta"
#   - item_names  : character vector of item names (or 'iname')
# Optional:
#   - pron_anchor_idx, det_anchor_idx : integer indices of anchors
#
# This version:
#   - clearer color/shape/size hierarchy (Okabe–Ito palette)
#   - subtle intervals for “other”, strong emphasis for anchors + reciprocals
#   - dashed reference at 0.5 (continuum center)
#   - wrapped, human-readable labels (“each_other” -> “each other”)
#   - Cairo PDF for the ↔ label if available
# ==============================================================================

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# ------------------------------ inputs ----------------------------------------
item_names <- if (exists("item_names")) item_names else if (exists("iname")) iname else NULL
stopifnot(!is.null(item_names))
N_items <- length(item_names)

# Optional indices default to empty
if (!exists("pron_anchor_idx")) pron_anchor_idx <- integer(0)
if (!exists("det_anchor_idx"))  det_anchor_idx  <- integer(0)

# Names for reciprocals (as in your matrix/fit parameter order)
recips <- c("each_other","one_another")
recip_idx <- match(recips, item_names)

# -------------------------- summarise eta -------------------------------------
qfun <- function(x) {
  qs <- as.numeric(stats::quantile(x, c(0.05, 0.50, 0.95)))
  names(qs) <- c("q5","q50","q95"); qs
}

eta_raw <- fit2$draws("eta")              # draws array
eta_summ <- posterior::summarise_draws(eta_raw, "mean", qfun) |>
  as_tibble() |>
  mutate(idx = as.integer(str_extract(variable, "(?<=\\[)\\d+(?=\\])"))) |>
  arrange(idx)

stopifnot(nrow(eta_summ) == N_items)

# ------------------------ grouping + display labels ---------------------------
grp <- rep("other", N_items)
grp[pron_anchor_idx] <- "pron_anchor"
grp[det_anchor_idx]  <- "det_anchor"
grp[recip_idx[!is.na(recip_idx)]] <- "reciprocal"

# Human-readable labels for plotting only
label_pretty <- item_names |>
  str_replace_all("_", " ") |>
  stringr::str_trim()

# If labels are long, wrap for readability (affects only axis text)
wrap_lab <- function(x, width = 22L) str_replace_all(stringr::str_wrap(x, width), "\n", "\n")
label_wrapped <- wrap_lab(label_pretty, width = 26L)

plot_df <- eta_summ |>
  transmute(
    item      = factor(label_wrapped, levels = label_wrapped),  # keep original order
    q5, q50, q95,
    group     = factor(grp, levels = c("det_anchor","reciprocal","other","pron_anchor")),
    is_other  = (grp == "other")
  ) |>
  arrange(q50) |>
  mutate(item = factor(item, levels = item))  # reorder by median

# --------------------------- aesthetics ---------------------------------------
# Okabe–Ito colorblind-safe palette
cols <- c(
  det_anchor  = "#D55E00",   # vermilion
  reciprocal  = "#CC79A7",   # reddish purple
  other       = "grey70",
  pron_anchor = "#0072B2"    # blue
)

shapes <- c(det_anchor = 16, reciprocal = 17, other = 16, pron_anchor = 16)

# Split data: faint “other” layer vs highlighted anchors/recips layer
df_other <- filter(plot_df, group == "other")
df_hi    <- filter(plot_df, group != "other")

# ------------------------------ plot ------------------------------------------
p <- ggplot() +
  # reference: center of the continuum
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4, color = "grey55") +

  # faint layer: all "other" items
  geom_linerange(
    data = df_other,
    aes(x = item, ymin = q5, ymax = q95),
    linewidth = 0.35, alpha = 0.35, color = cols["other"]
  ) +
  geom_point(
    data = df_other,
    aes(x = item, y = q50),
    size = 0.9, alpha = 0.6, color = cols["other"], shape = shapes["other"]
  ) +

  # highlight layer: anchors + reciprocals
  geom_linerange(
    data = df_hi,
    aes(x = item, ymin = q5, ymax = q95, color = group),
    linewidth = 0.7, alpha = 0.9
  ) +
  geom_point(
    data = df_hi,
    aes(x = item, y = q50, color = group, shape = group),
    size = 2.2, alpha = 0.95, stroke = 0.2
  ) +

  coord_flip() +
  scale_color_manual(values = cols, drop = FALSE) +
  scale_shape_manual(values = shapes, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    x = NULL,
    y = expression(paste("Posterior ", eta, " (pronoun \u2194 determinative)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    plot.margin = margin(8, 10, 8, 10)
  )

# ------------------------------ save ------------------------------------------
dir.create("figs", showWarnings = FALSE)

use_cairo <- isTRUE(capabilities("cairo"))
pdf_dev   <- if (use_cairo) grDevices::cairo_pdf else grDevices::pdf

ggsave("figs/eta_strip.pdf", plot = p, width = 7.5, height = 5, device = pdf_dev)
ggsave("figs/eta_strip.png", plot = p, width = 7.5, height = 5, dpi = 300)

message("Saved: figs/eta_strip.pdf and figs/eta_strip.png")
# ==============================================================================
