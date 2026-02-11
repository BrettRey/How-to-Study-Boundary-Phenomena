# helpers.R
# Shared helper functions for reciprocal analysis

# -------------------- Data cleaning helpers --------------------
key_us <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

normalize_binary <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x > 0))
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    y <- trimws(tolower(x))
    y[y %in% c("", "na", "nan")] <- NA
    out <- rep(NA_integer_, length(y))
    out[y %in% c("1","true","t","yes","y")] <- 1L
    out[y %in% c("0","false","f","no","n")] <- 0L
    suppressWarnings(num <- as.numeric(y))
    out[is.na(out) & !is.na(num)] <- as.integer(num > 0)
    out[is.na(out)] <- 0L
    return(out)
  }
  suppressWarnings(num <- as.numeric(x))
  out <- as.integer(num > 0)
  out[is.na(out)] <- 0L
  out
}

# -------------------- Data loading --------------------
load_linguistic_matrix <- function(csv_path) {
  DF <- read.csv(csv_path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  
  # Detect core columns
  cols_lower <- tolower(names(DF))
  lem_col <- if ("lemma" %in% cols_lower) names(DF)[match("lemma", cols_lower)] else
             if ("lexeme" %in% cols_lower) names(DF)[match("lexeme", cols_lower)] else
             if ("word" %in% cols_lower) names(DF)[match("word", cols_lower)] else
             stop("No lemma/lexeme/word column found")
  
  cat_col <- if ("class" %in% cols_lower) names(DF)[match("class", cols_lower)] else
             if ("category" %in% cols_lower) names(DF)[match("category", cols_lower)] else
             stop("No class/category column found")
  
  # Prepare labels
  labels <- data.frame(
    lemma = DF[[lem_col]],
    lemma_key = key_us(DF[[lem_col]]),
    category = tolower(trimws(DF[[cat_col]])),
    stringsAsFactors = FALSE
  )
  
  # Build binary feature matrix
  feature_cols <- setdiff(names(DF), c(lem_col, cat_col))
  feat <- DF[, feature_cols, drop=FALSE]
  feat[] <- lapply(feat, normalize_binary)
  Xbin <- as.matrix(feat)
  storage.mode(Xbin) <- "integer"
  
  list(
    features_binary = Xbin,
    labels = labels,
    raw_data = DF
  )
}

# -------------------- Distance calculations --------------------
jaccard_dist_vec <- function(x, y) {
  inter <- sum((x == 1L) & (y == 1L))
  uni <- sum((x == 1L) | (y == 1L))
  if (uni == 0L) return(0)
  1 - inter / uni
}

calculate_jaccard_matrix <- function(Xbin) {
  n <- nrow(Xbin)
  D <- matrix(0, n, n)
  rownames(D) <- rownames(Xbin)
  colnames(D) <- rownames(Xbin)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- jaccard_dist_vec(Xbin[i,], Xbin[j,])
      D[i,j] <- D[j,i] <- d
    }
  }
  D
}

# -------------------- Nearest neighbor analysis --------------------
find_nearest_neighbors <- function(D, labels, k = 5) {
  recip_idx <- which(labels$lemma_key %in% key_us(c("each_other", "one_another")))
  
  results <- list()
  for (i in recip_idx) {
    distances <- D[i, -i]  # exclude self
    sorted_idx <- order(distances)[1:k]
    
    # Map back to original indices
    original_idx <- setdiff(1:nrow(D), i)[sorted_idx]
    
    neighbor_info <- data.frame(
      lemma = labels$lemma[original_idx],
      category = labels$category[original_idx],
      distance = distances[sorted_idx]
    )
    
    results[[labels$lemma[i]]] <- neighbor_info
  }
  
  list(reciprocal_neighbors = results)
}

# -------------------- Classification helpers --------------------
calculate_log_likelihood <- function(cv_model, X, y_true) {
  # Get predicted probabilities
  pred_probs <- predict(cv_model, newx = X, s = "lambda.min", type = "response")[,1]
  
  # Convert y_true to binary (1 for pronoun, 0 for determinative)
  y_binary <- as.integer(y_true == "pronoun")
  
  # Calculate log-likelihood
  ll <- sum(y_binary * log(pred_probs + 1e-10) + 
            (1 - y_binary) * log(1 - pred_probs + 1e-10))
  ll
}

assess_calibration <- function(model, X, y_true, n_bins = 10) {
  pred_probs <- predict(model, newx = X, s = "lambda.min", type = "response")[,1]
  y_binary <- as.integer(y_true == "pronoun")
  
  # Bin predictions
  bins <- cut(pred_probs, breaks = seq(0, 1, length.out = n_bins + 1), 
              include.lowest = TRUE)
  
  calibration <- data.frame(
    bin = levels(bins),
    predicted = tapply(pred_probs, bins, mean),
    observed = tapply(y_binary, bins, mean),
    count = tapply(y_binary, bins, length)
  )
  
  # Remove NA bins
  calibration <- na.omit(calibration)
  
  # Calculate calibration error (mean absolute difference)
  error <- mean(abs(calibration$predicted - calibration$observed))
  
  list(
    calibration_table = calibration,
    error = error
  )
}

# -------------------- Bayesian helpers --------------------
calculate_bayes_factor <- function(prior_alpha, prior_beta, observed_success, observed_total) {
  # Bayes factor for pronoun hypothesis vs determinative hypothesis
  # Using the Savage-Dickey density ratio approximation
  
  # Prior probability at 0.5 (null hypothesis of no difference)
  prior_at_half <- dbeta(0.5, prior_alpha, prior_beta)
  
  # Posterior parameters
  post_alpha <- prior_alpha + observed_success
  post_beta <- prior_beta + (observed_total - observed_success)
  
  # Posterior probability at 0.5
  posterior_at_half <- dbeta(0.5, post_alpha, post_beta)
  
  # Bayes factor
  bf <- posterior_at_half / prior_at_half
  
  # Return BF in favor of pronoun hypothesis (>0.5) vs determinative (<0.5)
  1 / bf  # Invert because we want evidence for pronoun (away from 0.5)
}

# -------------------- Feature block analysis --------------------
identify_feature_blocks <- function(feature_names) {
  # Categorize features into morphological, syntactic, semantic, phonological
  # This is a simplified version - you'd need to customize based on your actual features
  
  blocks <- list(
    morphological = grep("morph|inflect|deriv|compound|affix", feature_names, ignore.case = TRUE),
    syntactic = grep("synt|position|order|complement|modifier|function", feature_names, ignore.case = TRUE),
    semantic = grep("sem|meaning|animate|definite|refer|plural|singular", feature_names, ignore.case = TRUE),
    phonological = grep("phon|sound|initial|final|stress", feature_names, ignore.case = TRUE)
  )
  
  # Assign remaining features
  all_assigned <- unlist(blocks)
  remaining <- setdiff(1:length(feature_names), all_assigned)
  blocks$other <- remaining
  
  blocks
}

# -------------------- Reporting helpers --------------------
format_p_value <- function(p) {
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

report_summary <- function(results_list) {
  cat("\n=== ANALYSIS SUMMARY ===\n")
  cat(sprintf("Date: %s\n", Sys.Date()))
  cat(sprintf("R version: %s\n", R.version.string))
  cat("\n")
  
  # Add specific results reporting here
  invisible(results_list)
}