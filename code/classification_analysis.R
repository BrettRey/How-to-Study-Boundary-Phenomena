# 02_classification_analysis.R (REVISED)
# Classification using the FULL dataset, not just matched subset

library(glmnet)
library(caret)
library(tidyverse)

source("helpers.R")

# -------------------- Load data --------------------
DF <- load_linguistic_matrix("matrix_clean.csv")
Xbin <- DF$features_binary
labels <- DF$labels

# Identify reciprocals (to be held out)
recip_idx <- which(labels$lemma_key %in% key_us(c("each_other", "one_another")))

# Use ALL pronouns and determinatives, not just the matched subset
all_pron_idx <- which(grepl("pronoun", labels$category, ignore.case = TRUE))
all_det_idx <- which(grepl("determinative", labels$category, ignore.case = TRUE))

cat(sprintf("Using full dataset:\n"))
cat(sprintf("  Pronouns: %d\n", length(all_pron_idx)))
cat(sprintf("  Determinatives: %d\n", length(all_det_idx)))
cat(sprintf("  Reciprocals (held out): %d\n", length(recip_idx)))

# -------------------- Train classifier excluding reciprocals --------------------
train_idx <- setdiff(1:nrow(Xbin), recip_idx)
X_train <- Xbin[train_idx, ]
y_train <- ifelse(train_idx %in% all_pron_idx, "pronoun", "determinative")

# Check class balance
cat(sprintf("\nTraining set composition:\n"))
cat(sprintf("  Pronouns: %d\n", sum(y_train == "pronoun")))
cat(sprintf("  Determinatives: %d\n", sum(y_train == "determinative")))

# Use cross-validation to select regularization
cv_model <- cv.glmnet(X_train, y_train, family = "binomial", 
                       alpha = 0, nfolds = 10)  # More folds with larger dataset

# Get predictions for reciprocals
X_test <- Xbin[recip_idx, ]
pred_probs <- predict(cv_model, newx = X_test, 
                       s = "lambda.min", type = "response")[,1]

cat("\nClassifier predictions for reciprocals:\n")
for (i in 1:length(recip_idx)) {
  cat(sprintf("  %s: P(pronoun) = %.3f\n", 
              labels$lemma[recip_idx[i]], pred_probs[i]))
}

# -------------------- Model comparison --------------------
# Model 1: Reciprocals as pronouns
y_all_pron <- ifelse(1:nrow(Xbin) %in% c(all_pron_idx, recip_idx), 
                      "pronoun", "determinative")
cv_model_pron <- cv.glmnet(Xbin, y_all_pron, family = "binomial", 
                            alpha = 0, nfolds = 10)
ll_pron <- calculate_log_likelihood(cv_model_pron, Xbin, y_all_pron)

# Model 2: Reciprocals as determinatives  
y_all_det <- ifelse(1:nrow(Xbin) %in% all_pron_idx, "pronoun", "determinative")
cv_model_det <- cv.glmnet(Xbin, y_all_det, family = "binomial", 
                           alpha = 0, nfolds = 10)
ll_det <- calculate_log_likelihood(cv_model_det, Xbin, y_all_det)

cat(sprintf("\nModel comparison:\n"))
cat(sprintf("  Log-likelihood (reciprocals as pronouns): %.3f\n", ll_pron))
cat(sprintf("  Log-likelihood (reciprocals as determinatives): %.3f\n", ll_det))
cat(sprintf("  Likelihood ratio: %.3f\n", exp(ll_pron - ll_det)))

if (ll_pron > ll_det) {
  cat("  → Better fit when treating reciprocals as pronouns\n")
} else {
  cat("  → Better fit when treating reciprocals as determinatives\n")
}

# -------------------- Calibration analysis --------------------
calibration_data <- assess_calibration(cv_model, X_train, y_train)
cat(sprintf("\nCalibration error: %.3f\n", calibration_data$error))