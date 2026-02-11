# restyle_plots.R — Regenerate presentation plots with house-style typography & colours
# Run from the code/ directory (or adjust paths)
# Produces: plots/mca_reciprocals.png, plots/spec_curve.png, plots/permutation_null.png
# NOTE: eta_strip.pdf requires Stan posteriors and must be rerun separately.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(showtext)
  library(sysfonts)
})

# ──────────────────── House palette & theme ────────────────────
# From .house-style/plot_style.py
primary   <- "#2E5077"
secondary <- "#E85D4C"
tertiary  <- "#4DA375"
quaternary <- "#9B6B9E"
quinary   <- "#D4A03E"
light     <- "#E8E8E8"
dark      <- "#2D2D2D"
accent    <- "#6AADE4"

# Register EB Garamond (local install)
font_add("EB Garamond",
         regular = "/Users/brettreynolds/Library/Fonts/EBGaramond12-Regular.otf",
         italic  = "/Users/brettreynolds/Library/Fonts/EBGaramond12-Italic.otf")
showtext_auto()

house_theme <- function(base_size = 13) {
  theme_bw(base_size = base_size, base_family = "EB Garamond") %+replace%
    theme(
      # No grid
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Remove top/right spines
      axis.line.x = element_line(colour = dark, linewidth = 0.4),
      axis.line.y = element_line(colour = dark, linewidth = 0.4),
      panel.border = element_blank(),
      # Frameless legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = rel(0.9)),
      # Axis text
      axis.text = element_text(colour = dark),
      axis.title = element_text(colour = dark),
      # Title
      plot.title = element_text(colour = dark, face = "plain", size = rel(1.1)),
      plot.subtitle = element_text(colour = "grey40", size = rel(0.9)),
      # Margins
      plot.margin = margin(8, 12, 8, 8)
    )
}

# ──────────────────── 1) MCA plot ────────────────────
message(">>> MCA plot")
library(FactoMineR)
library(ggrepel)
source("helpers.R")

DF <- load_linguistic_matrix("../data/matrix_clean.csv")
Xbin <- DF$features_binary
labels <- DF$labels

fused_names <- c("someone","anyone","anything","everything","somebody","anybody")
pronoun_names <- c("he","him","himself","she","her","herself","they","them","themselves")
reciprocal_names <- c("each_other","one_another")

mca_data <- as.data.frame(Xbin)
mca_data[] <- lapply(mca_data, as.factor)
mca_res <- MCA(mca_data, graph = FALSE)

mca_plot_data <- data.frame(
  Dim1 = mca_res$ind$coord[,1],
  Dim2 = mca_res$ind$coord[,2],
  Category = labels$category,
  IsReciprocal = labels$lemma_key %in% key_us(reciprocal_names),
  IsFused = labels$lemma_key %in% key_us(fused_names),
  Lemma = labels$lemma
)

# Map categories to house colours
cat_colours <- c(
  "determinative" = secondary,   # red
  "pronoun"       = accent       # blue
)

p_mca <- ggplot(mca_plot_data, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = Category, shape = IsReciprocal, size = IsReciprocal),
             alpha = 0.8) +
  scale_colour_manual(values = cat_colours) +
  scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3.5), guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "Other", "TRUE" = "Reciprocal")) +
  geom_text_repel(
    data = filter(mca_plot_data, IsReciprocal),
    aes(label = gsub("_", " ", Lemma)),
    family = "EB Garamond", size = 6, fontface = "italic",
    nudge_y = 0.15, segment.color = "grey60"
  ) +
  stat_ellipse(aes(colour = Category), level = 0.68, linewidth = 0.5,
               linetype = "dashed") +
  house_theme(base_size = 16) +
  labs(
    x = sprintf("Dimension 1 (%.1f%%)", mca_res$eig[1, "percentage of variance"]),
    y = sprintf("Dimension 2 (%.1f%%)", mca_res$eig[2, "percentage of variance"]),
    colour = NULL, shape = NULL
  ) +
  theme(legend.position = "top")

ggsave("../plots/mca_reciprocals.png", p_mca,
       width = 7, height = 5.5, dpi = 300, bg = "white")
message("  Saved plots/mca_reciprocals.png")


# ──────────────────── 2) Spec curve plot ────────────────────
message(">>> Spec curve plot")

spec_res <- read.csv("../data/spec_curve.csv", stringsAsFactors = FALSE)

# Block labels for facets
block_labels <- c(
  "none" = "All properties",
  "morph" = "Without morphology",
  "synt"  = "Without syntax",
  "sem"   = "Without semantics",
  "phon"  = "Without phonology"
)

dist_colours <- c(
  "jaccard" = primary,
  "dice"    = secondary,
  "hamming" = tertiary,
  "wj_idf"  = quinary
)

dist_labels <- c(
  "jaccard" = "Jaccard",
  "dice"    = "Dice",
  "hamming" = "Hamming",
  "wj_idf"  = "IDF-weighted Jaccard"
)

p_spec <- spec_res |>
  mutate(
    idx = row_number(),
    mean_delta = (delta_each + delta_onean) / 2,
    block = factor(block, levels = names(block_labels), labels = block_labels)
  ) |>
  ggplot(aes(idx, mean_delta, shape = alpha, colour = dist)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~block, scales = "free_x") +
  scale_colour_manual(values = dist_colours, labels = dist_labels, name = "Metric:") +
  scale_shape_manual(values = c("ridge" = 16, "elastic05" = 17),
                     labels = c("ridge" = "Ridge", "elastic05" = "Elastic net"),
                     name = "Regularization:") +
  house_theme(base_size = 16) +
  labs(
    y = expression(paste(Delta, "distance (pronoun \u2212 determinative)")),
    x = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(family = "EB Garamond", size = 14),
    legend.position = "none"
  )

ggsave("../plots/spec_curve.png", p_spec,
       width = 8, height = 4.8, dpi = 300, bg = "white")
message("  Saved plots/spec_curve.png")


# ──────────────────── 3) Permutation null plot ────────────────────
message(">>> Permutation null plot (running quasiswap, B = 5000...)")
library(vegan)
library(proxy)

dat <- read.csv("../data/matrix_clean.csv", stringsAsFactors = FALSE, check.names = FALSE)
items <- dat$lemma
y_raw <- factor(dat$class, levels = c("determinative", "pronoun"))
X <- as.matrix(dat[, setdiff(names(dat), c("lemma", "class"))])
storage.mode(X) <- "numeric"; X[is.na(X)] <- 0
rownames(X) <- items

recips <- c("each_other", "one_another")
pron_set <- items[y_raw == "pronoun"]
fused_col <- names(dat)[grepl("fused", names(dat), ignore.case = TRUE)][1]
is_fused <- if (!is.na(fused_col)) dat[[fused_col]] == 1 else rep(FALSE, nrow(dat))
fused_set <- items[y_raw == "determinative" & is_fused]

n_per_class <- min(6, length(setdiff(fused_set, recips)), length(setdiff(pron_set, recips)))
fused_n <- head(setdiff(fused_set, recips), n_per_class)
pron_n  <- head(setdiff(pron_set, recips), n_per_class)
subset_items <- c(recips, fused_n, pron_n)
Xs <- X[subset_items, , drop = FALSE]

as_mat <- function(D) {
  M <- as.matrix(D)
  rn <- attr(D, "Labels")
  rownames(M) <- colnames(M) <- rn
  M
}

stat_delta_subset <- function(M, recips, pron_targets, fused_targets) {
  D <- proxy::dist(M, method = "Jaccard")
  Mmat <- as.matrix(D)
  rownames(Mmat) <- colnames(Mmat) <- attr(D, "Labels")
  delta_item <- function(item) {
    m_pron <- mean(Mmat[item, pron_targets, drop = FALSE])
    m_fuse <- mean(Mmat[item, fused_targets, drop = FALSE])
    m_pron - m_fuse
  }
  mean(sapply(recips, delta_item))
}

set.seed(2025)
T_obs <- stat_delta_subset(Xs, recips, pron_n, fused_n)

B <- 5000
perm <- vegan::permatswap(Xs, method = "quasiswap", mtype = "prab",
                          fixedmar = "both", times = B, burnin = 2000)
T_null <- vapply(perm$perm,
                 function(M) stat_delta_subset(M, recips, pron_n, fused_n),
                 numeric(1))
T0 <- mean(T_null)

null_df <- data.frame(delta = T_null)

p_perm <- ggplot(null_df, aes(delta)) +
  geom_histogram(bins = 40, fill = accent, colour = "white", linewidth = 0.2,
                 alpha = 0.8) +
  geom_vline(xintercept = T0, linetype = 3, colour = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = T_obs, linetype = 2, colour = secondary, linewidth = 0.8) +
  house_theme(base_size = 16) +
  labs(
    x = expression(paste(Delta, " (under null)")),
    y = "Count"
  )

ggsave("../plots/permutation_null.png", p_perm,
       width = 7, height = 4.5, dpi = 300, bg = "white")
message("  Saved plots/permutation_null.png")

message("\n>>> Done. eta_strip.pdf requires Stan posteriors -- rerun from 02_ppc_cmdstan.R + eta_strip.R")
