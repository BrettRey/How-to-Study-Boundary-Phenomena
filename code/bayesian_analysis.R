# 04_bayesian_analysis.R
# Bayesian synthesis using actual classifier predictions

library(tidyverse)
library(glmnet)
library(ggplot2)

source("00_helpers.R")

# -------------------- Get classifier predictions --------------------
DF <- load_linguistic_matrix("matrix_clean.csv")
Xbin <- DF$features_binary
labels <- DF$labels

# Identify groups
recip_idx <- which(labels$lemma_key %in% key_us(c("each_other", "one_another")))
all_pron_idx <- which(grepl("pronoun", labels$category, ignore.case = TRUE))
all_det_idx <- which(grepl("determinative", labels$category, ignore.case = TRUE))

# Train classifier excluding reciprocals
train_idx <- setdiff(1:nrow(Xbin), recip_idx)
X_train <- Xbin[train_idx, ]
y_train <- ifelse(train_idx %in% all_pron_idx, "pronoun", "determinative")

# Get predictions for reciprocals
cv_model <- cv.glmnet(X_train, y_train, family = "binomial", 
                       alpha = 0, nfolds = 10)
X_test <- Xbin[recip_idx, ]
pred_probs <- predict(cv_model, newx = X_test, 
                       s = "lambda.min", type = "response")[,1]

cat("Classifier predictions:\n")
for (i in 1:length(recip_idx)) {
  cat(sprintf("  %s: P(pronoun) = %.3f\n", 
              labels$lemma[recip_idx[i]], pred_probs[i]))
}

# -------------------- Bayesian update using binary classification --------------------
cat("\n=== Bayesian Analysis ===\n")

# Prior: 0.85 probability that reciprocals are pronouns
prior_prob <- 0.85
prior_alpha <- prior_prob * 100  # effective sample size of 100
prior_beta <- (1 - prior_prob) * 100

# Count how many reciprocals are classified as pronouns (P > 0.5)
observed_pronoun <- sum(pred_probs > 0.5)
observed_total <- length(pred_probs)

# Calculate posterior using conjugate beta-binomial
posterior_alpha <- prior_alpha + observed_pronoun
posterior_beta <- prior_beta + (observed_total - observed_pronoun)

prior_mean <- prior_alpha / (prior_alpha + prior_beta)
posterior_mean <- posterior_alpha / (posterior_alpha + posterior_beta)

# Calculate Bayes Factor using odds ratios (CORRECTED)
prior_odds_pron <- prior_alpha / prior_beta  # 85/15 = 5.667
posterior_odds_pron <- posterior_alpha / posterior_beta  

# BF for pronoun vs determinative
bf_pron_det <- posterior_odds_pron / prior_odds_pron
bf_det_pron <- 1 / bf_pron_det

cat(sprintf("Prior P(pronoun): %.3f\n", prior_mean))
cat(sprintf("Observed: %d/%d classified as pronouns\n", observed_pronoun, observed_total))
cat(sprintf("Posterior P(pronoun): %.3f\n", posterior_mean))
cat(sprintf("\nOdds ratio calculation:\n"))
cat(sprintf("  Prior odds (pronoun:det): %.3f\n", prior_odds_pron))
cat(sprintf("  Posterior odds (pronoun:det): %.3f\n", posterior_odds_pron))
cat(sprintf("  BF (pronoun:det): %.3f\n", bf_pron_det))
cat(sprintf("  BF (det:pronoun): %.3f\n", bf_det_pron))

# Interpret using Kass-Raftery scale
if (bf_det_pron < 1) {
  bf_to_interpret <- bf_pron_det
  direction <- "pronoun"
} else {
  bf_to_interpret <- bf_det_pron
  direction <- "determinative"
}

if (bf_to_interpret < 1.0) {
  # This shouldn't happen given our setup
  bf_interpretation <- "evidence against the hypothesis"
} else if (bf_to_interpret <= 3.2) {
  bf_interpretation <- "not worth more than a bare mention"
} else if (bf_to_interpret <= 10) {
  bf_interpretation <- "substantial"
} else if (bf_to_interpret <= 20) {
  bf_interpretation <- "strong"
} else if (bf_to_interpret <= 150) {
  bf_interpretation <- "very strong"
} else {
  bf_interpretation <- "decisive"
}

cat(sprintf("\nKass-Raftery interpretation: Evidence for %s is %s (BF = %.2f)\n", 
            direction, bf_interpretation, bf_to_interpret))

# -------------------- Sensitivity analysis --------------------
cat("\n=== Sensitivity Analysis ===\n")

prior_probs <- seq(0.25, 0.95, by = 0.05)
posteriors <- numeric(length(prior_probs))
bayes_factors <- numeric(length(prior_probs))

for (i in 1:length(prior_probs)) {
  p <- prior_probs[i]
  a <- p * 100
  b <- (1 - p) * 100
  
  # Update with observed data
  post_a <- a + observed_pronoun
  post_b <- b + (observed_total - observed_pronoun)
  posteriors[i] <- post_a / (post_a + post_b)
  
  # Calculate BF
  prior_odds <- a / b
  post_odds <- post_a / post_b
  bayes_factors[i] <- (1 / (post_odds / prior_odds))  # BF for determinative
}

# Create sensitivity plot
sensitivity_data <- data.frame(
  Prior = prior_probs,
  Posterior = posteriors,
  BF_det = bayes_factors
)

sensitivity_plot <- ggplot(sensitivity_data, aes(x = Prior, y = Posterior)) +
  geom_line(linewidth = 1, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0.85, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dotted") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Bayesian Sensitivity Analysis",
       subtitle = sprintf("Observed: %d/%d classified as pronouns", 
                         observed_pronoun, observed_total),
       x = "Prior P(reciprocals are pronouns)",
       y = "Posterior P(reciprocals are pronouns)")

ggsave("bayesian_sensitivity.png", sensitivity_plot, 
       width = 8, height = 6, 
       bg = "white", 
       dpi = 300)

# Plot Bayes factors across priors
bf_plot <- ggplot(sensitivity_data, aes(x = Prior, y = BF_det)) +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 3.2, linetype = "dotted", color = "orange", alpha = 0.7) +
  geom_hline(yintercept = 20, linetype = "dotted", color = "red", alpha = 0.7) +
  geom_vline(xintercept = 0.85, color = "red", linetype = "dashed") +
  scale_y_log10() +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Bayes Factor Sensitivity",
       subtitle = "BF(det:pron) across different priors",
       x = "Prior P(reciprocals are pronouns)",
       y = "Bayes Factor (det:pron) [log scale]") +
  annotate("text", x = 0.3, y = 3.2, label = "Substantial (3.2)", 
           hjust = 0, vjust = -0.5, size = 3, color = "orange") +
  annotate("text", x = 0.3, y = 20, label = "Strong (20)", 
           hjust = 0, vjust = -0.5, size = 3, color = "red")

ggsave("bayes_factor_sensitivity.png", bf_plot, 
       width = 8, height = 6, 
       bg = "white", 
       dpi = 300)

# -------------------- Summary --------------------
cat("\n=== SUMMARY ===\n")
cat("Evidence from classifier:\n")
for (i in 1:length(pred_probs)) {
  classification <- ifelse(pred_probs[i] > 0.5, "pronoun", "determinative")
  cat(sprintf("  %s: P(pronoun)=%.3f → classified as %s\n", 
              labels$lemma[recip_idx[i]], pred_probs[i], classification))
}

cat(sprintf("\nWith prior P(pronoun) = %.2f:\n", prior_prob))
cat(sprintf("  Posterior P(pronoun) = %.3f\n", posterior_mean))
cat(sprintf("  BF(det:pron) = %.2f\n", bf_det_pron))
cat(sprintf("  Interpretation: %s\n", bf_interpretation))

if (posterior_mean < prior_prob) {
  cat("→ Evidence shifts belief toward determinative (but weakly)\n")
} else if (posterior_mean > prior_prob) {
  cat("→ Evidence shifts belief toward pronoun (but weakly)\n")
} else {
  cat("→ Evidence does not change prior belief\n")
}

cat("\nKass-Raftery scale reference:\n")
cat("  BF 1-3.2:    Not worth more than a bare mention\n")
cat("  BF 3.2-10:   Substantial\n") 
cat("  BF 10-20:    Strong\n")
cat("  BF 20-150:   Very strong\n")
cat("  BF >150:     Decisive\n")