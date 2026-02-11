# 05_core_core_permutation.R
# Permutation test for Δ_core using COCA-based cores (hard-coded).
# Δ_core = mean_r [ d(r, det-core medoid) - d(r, pron-core medoid) ]

suppressPackageStartupMessages({
  library(tidyverse); library(Matrix); library(proxy); library(vegan); library(ggplot2)
})
set.seed(2025); dir.create("out", showWarnings = FALSE)

# --- load matrix --------------------------------------------------------------
dat <- read.csv("matrix_clean.csv", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("lemma","class") %in% names(dat)))
items <- dat$lemma
y_raw <- factor(dat$class, levels = c("determinative","pronoun"))
X <- as.matrix(dat[, setdiff(names(dat), c("lemma","class"))]); storage.mode(X) <- "numeric"; X[is.na(X)] <- 0
rownames(X) <- items

# --- hard-coded COCA cores (present in your matrix) ---------------------------
PRON_CORE <- c("I","you_pron_sing","it_plain","he","we_pron","they_plur","she","me","who_int","him")
DET_CORE  <- c("the","a","this","that","what_det","all","some","these","those","no")

stopifnot(all(PRON_CORE %in% items), all(DET_CORE %in% items))

recips <- intersect(c("each_other","one_another"), items); stopifnot(length(recips)==2)

# --- helpers -----------------------------------------------------------------
as_mat <- function(D){ M <- as.matrix(D); lab <- attr(D,"Labels"); if (!is.null(lab)) rownames(M) <- colnames(M) <- lab; M }
medoid_of <- function(D, ids){ M <- as_mat(D)[ids, ids, drop=FALSE]; ids[ which.min(rowSums(M)) ] }
delta_core <- function(D, recips, pron_core, det_core){
  M <- as_mat(D)
  m_pron <- medoid_of(D, pron_core)
  m_det  <- medoid_of(D, det_core)
  mean(M[recips, m_det] - M[recips, m_pron])
}

# restrict to rows actually used (speeds up)
rows_used <- unique(c(recips, PRON_CORE, DET_CORE))
Xu <- X[rows_used, , drop = FALSE]

# --- observed ----------------------------------------------------------------
D_obs <- proxy::dist(Xu, method = "Jaccard")
T_obs <- delta_core(D_obs, recips, PRON_CORE, DET_CORE)

# --- quasiswap null (recompute medoids each draw) ----------------------------
B <- 5000; burn <- 1500
perm <- vegan::permatswap(Xu, method="quasiswap", mtype="prab",
                          fixedmar="both", times=B, burnin=burn)

T_null <- numeric(B)
for (b in seq_len(B)) {
  Db <- proxy::dist(perm$perm[[b]], method="Jaccard")
  T_null[b] <- delta_core(Db, recips, PRON_CORE, DET_CORE)
}

# two-sided p around permutation mean
T0 <- mean(T_null)
p_two_sided <- (sum(abs(T_null - T0) >= abs(T_obs - T0)) + 1) / (B + 1)

# --- outputs -----------------------------------------------------------------
write.csv(data.frame(delta_core_null = T_null), "out/core_core_perm_null.csv", row.names = FALSE)
writeLines(sprintf("Observed Δ_core = %.4f; E0[Δ_core] = %.4f; two-sided p = %.3f",
                   T_obs, T0, p_two_sided),
           "out/core_core_perm_summary.txt")

p <- ggplot(data.frame(delta = T_null), aes(delta)) +
  geom_histogram(bins = 40, linewidth = 0.2) +
  geom_vline(xintercept = T0,  linetype = 3) +
  geom_vline(xintercept = T_obs, linetype = 2) +
  labs(title = "Permutation null for Δ_core (quasiswap; cores only)",
       x = expression(Delta[core]~"= d(recips, det-core medoid) - d(recips, pron-core medoid)"),
       y = "count") +
  theme_minimal(base_size = 12)
ggsave("out/core_core_perm_null.png", p, width = 7, height = 4.5, dpi = 300, bg = "white")
