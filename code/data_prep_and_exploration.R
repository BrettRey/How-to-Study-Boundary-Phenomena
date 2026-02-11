# 01_explore_reciprocals.R
# Exploratory visualization and distance analysis

library(tidyverse)
library(FactoMineR)  # for MCA
library(factoextra)  # for visualization
library(ape)         # for PCoA
library(pheatmap)    # for heatmaps

# -------------------- Configuration --------------------
csv_path <- "matrix_clean.csv"
set.seed(20250819)

# Target lexeme sets
fused_names <- c("someone","anyone","anything","everything","somebody","anybody")
pronoun_names <- c("he","him","himself","she","her","herself","they","them","themselves")
reciprocal_names <- c("each_other","one_another")

# -------------------- Load and prepare data --------------------
source("helpers.R")  # Contains the helper functions from your script

DF <- load_linguistic_matrix(csv_path)
Xbin <- DF$features_binary
labels <- DF$labels

# -------------------- MCA --------------------
# Prepare data for MCA (needs factors)
mca_data <- as.data.frame(Xbin)
mca_data[] <- lapply(mca_data, as.factor)

# Run MCA
mca_res <- MCA(mca_data, graph = FALSE)

# Create custom visualization
mca_plot_data <- data.frame(
  Dim1 = mca_res$ind$coord[,1],
  Dim2 = mca_res$ind$coord[,2],
  Category = labels$category,
  IsReciprocal = labels$lemma_key %in% key_us(reciprocal_names),
  IsFused = labels$lemma_key %in% key_us(fused_names),
  Lemma = labels$lemma
)

p_mca <- ggplot(mca_plot_data, aes(x = Dim1, y = Dim2)) +
  # Map color, shape, AND size to your data
  geom_point(aes(color = Category, shape = IsReciprocal, size = IsReciprocal)) +
  
  # Manually set the size for dots (FALSE) and triangles (TRUE)
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3)) +
  
  geom_text(data = filter(mca_plot_data, IsReciprocal), 
            aes(label = Lemma), vjust = -1) +
  stat_ellipse(aes(color = Category), level = 0.68) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "MCA: Reciprocals in Feature Space",
       subtitle = sprintf("Dimensions 1-2 explain %.1f%% of variance", 
                          sum(mca_res$eig[1:2, "percentage of variance"])))

ggsave("mca_reciprocals.png", p_mca, width = 10, height = 8, bg = "white")

# -------------------- PCoA with Jaccard --------------------
# Calculate Jaccard distance matrix
D_jaccard <- calculate_jaccard_matrix(Xbin)

# Run PCoA
pcoa_res <- pcoa(as.dist(D_jaccard))

# Define items to label
items_to_label <- c("each_other", "one_another",  # reciprocals
                    "someone", "anyone", "no_one",   # compound pronouns
                    "somebody", "anybody", "nobody",
                    "something", "anything", "nothing", 
                    "everything", "everybody",
                    "once", "twice", "thrice")  # multiplicatives

# Extract coordinates and create plot
pcoa_plot_data <- data.frame(
  PCo1 = pcoa_res$vectors[,1],
  PCo2 = pcoa_res$vectors[,2],
  Category = labels$category,
  IsReciprocal = labels$lemma_key %in% key_us(reciprocal_names),
  IsFused = labels$lemma_key %in% key_us(fused_names),
  ShouldLabel = labels$lemma_key %in% key_us(items_to_label),
  Lemma = labels$lemma
)

# Load the ggrepel library for better labels (optional but recommended)
# install.packages("ggrepel")
library(ggrepel)

p_pcoa <- ggplot(pcoa_plot_data, aes(x = PCo1, y = PCo2)) +
  # This part is correct
  geom_point(aes(color = Category, shape = IsReciprocal, size = IsReciprocal)) +
  scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3)) + # Adjusted dot size slightly
  
  # --- FIX IS HERE ---
  # Use geom_text_repel instead of geom_text and REMOVE the size = 3
  geom_text_repel(data = filter(pcoa_plot_data, ShouldLabel), 
                  aes(label = Lemma),
                  max.overlaps = 15) + # ggrepel is better at avoiding overlap
  
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.key = element_blank()) + # Hide size guide
  labs(title = "PCoA with Jaccard Distances",
       subtitle = sprintf("Axes 1-2 explain %.1f%% of variance", 
                          100 * sum(pcoa_res$values$Relative_eig[1:2]))) +
  guides(size = "none") # This hides the unnecessary size legend

ggsave("pcoa_reciprocals.png", p_pcoa, width = 10, height = 8, bg = "white")

# -------------------- k-Nearest Neighbors --------------------
knn_analysis <- find_nearest_neighbors(D_jaccard, labels, k = 5)
print(knn_analysis$reciprocal_neighbors)

# -------------------- Distance Heatmap (FIXED) --------------------
# Subset to key comparisons
key_items <- c(reciprocal_names, fused_names, pronoun_names)
key_idx <- which(labels$lemma %in% key_items)
D_subset <- D_jaccard[key_idx, key_idx]

# Need to ensure rownames and colnames are set
rownames(D_subset) <- labels$lemma[key_idx]
colnames(D_subset) <- labels$lemma[key_idx]

# Create annotation dataframe
annotation_row <- data.frame(
  Type = factor(case_when(
    labels$lemma[key_idx] %in% reciprocal_names ~ "Reciprocal",
    labels$lemma[key_idx] %in% fused_names ~ "Fused Det",
    labels$lemma[key_idx] %in% pronoun_names ~ "Pronoun",
    TRUE ~ "Other"
  ), levels = c("Reciprocal", "Fused Det", "Pronoun")),
  row.names = labels$lemma[key_idx]
)

# Define colors for annotation
ann_colors <- list(
  Type = c("Reciprocal" = "orange", "Fused Det" = "blue", "Pronoun" = "green")
)

pheatmap(D_subset, 
         annotation_row = annotation_row,
         annotation_col = annotation_row,
         annotation_colors = ann_colors,
         main = "Jaccard Distances: Key Comparisons",
         filename = "distance_heatmap.png",
         color = colorRampPalette(c("darkblue", "white", "darkred"))(50),
         breaks = seq(0, 1, length.out = 51))

# -------------------- Summary Statistics --------------------
# Mean distances from reciprocals to each category
recip_idx <- which(labels$lemma_key %in% key_us(reciprocal_names))
fused_idx <- which(labels$lemma_key %in% key_us(fused_names))
pron_idx <- which(labels$lemma_key %in% key_us(pronoun_names) & 
                   grepl("pronoun", labels$category, ignore.case = TRUE))

for (i in recip_idx) {
  d_to_fused <- mean(D_jaccard[i, fused_idx])
  d_to_pron <- mean(D_jaccard[i, pron_idx])
  cat(sprintf("%s: d_to_fused = %.3f, d_to_pron = %.3f, difference = %.3f\n",
              labels$lemma[i], d_to_fused, d_to_pron, d_to_pron - d_to_fused))
}