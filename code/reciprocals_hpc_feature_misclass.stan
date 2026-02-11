data {
  int<lower=1> N_items;
  int<lower=1> N_features;
  int<lower=1> N_blocks;
  array[N_features] int<lower=1, upper=N_blocks> block_id;

  array[N_items, N_features] int<lower=0, upper=1> y;     // observed 0/1
  array[N_items, N_features] int<lower=0, upper=1> mask;  // 1 = held-out

  // Anchor priors for eta_raw (logit scale)
  vector[N_items] mu_eta_raw;
  vector<lower=0>[N_items] sd_eta_raw;

  // Masked indices (sizes must come from data)
  int<lower=0> N_mask;
  array[N_mask] int<lower=1, upper=N_items> idx_i;
  array[N_mask] int<lower=1, upper=N_features> idx_j;
}

parameters {
  // Latent item coordinates (logit scale); eta = inv_logit(eta_raw)
  vector[N_items] eta_raw;

  // Block-level intercepts/slopes
  vector[N_blocks] alpha_b;
  vector[N_blocks] beta_b;

  // Random intercepts
  vector[N_features] u_feat;
  vector[N_items]   v_item;
  real<lower=0> tau_feat;
  real<lower=0> tau_item;

  // Feature-level misclassification with partial pooling (by block)
  vector[N_blocks] mu_gplus_b;
  vector<lower=0>[N_blocks] sd_gplus_b;
  vector[N_blocks] mu_gminus_b;
  vector<lower=0>[N_blocks] sd_gminus_b;

  vector[N_features] z_gplus;
  vector[N_features] z_gminus;
}

transformed parameters {
  vector[N_items] eta = inv_logit(eta_raw);
  vector<lower=0, upper=1>[N_features] gplus;
  vector<lower=0, upper=1>[N_features] gminus;

  for (j in 1:N_features) {
    int b = block_id[j];
    gplus[j]  = inv_logit(mu_gplus_b[b]  + sd_gplus_b[b]  * z_gplus[j]);
    gminus[j] = inv_logit(mu_gminus_b[b] + sd_gminus_b[b] * z_gminus[j]);
  }
}

model {
  // Priors
  eta_raw ~ normal(mu_eta_raw, sd_eta_raw);

  alpha_b ~ normal(0, 1.5);
  beta_b  ~ normal(0, 1.5);

  u_feat  ~ normal(0, tau_feat);
  v_item  ~ normal(0, tau_item);
  tau_feat ~ normal(0, 1);
  tau_item ~ normal(0, 1);

  // Misclassification hyperpriors: mild skepticism (~2%)
  mu_gplus_b  ~ normal(logit(0.02), 1.0);
  mu_gminus_b ~ normal(logit(0.02), 1.0);
  sd_gplus_b  ~ normal(0, 0.5);
  sd_gminus_b ~ normal(0, 0.5);

  z_gplus  ~ normal(0, 1);
  z_gminus ~ normal(0, 1);

  // Likelihood (marginal over true x_ij)
  for (i in 1:N_items) {
    for (j in 1:N_features) {
      if (mask[i, j] == 0) {
        int b = block_id[j];
        real lp   = alpha_b[b] + beta_b[b] * eta[i] + u_feat[j] + v_item[i];
        real p_ij = inv_logit(lp);                                 // Pr(true=1)
        real q_ij = (1 - gminus[j]) * p_ij + gplus[j] * (1 - p_ij); // Pr(y=1)
        y[i, j] ~ bernoulli(q_ij);
      }
    }
  }
}

generated quantities {
  // Masked-cell log-lik for Δelpd
  vector[N_mask] log_lik_mask;

  // Posterior predictive expectations for row/col sums
  vector[N_items] row_sum_exp;
  vector[N_features] col_sum_exp;

  // Expectations
  for (i in 1:N_items) row_sum_exp[i] = 0;
  for (j in 1:N_features) col_sum_exp[j] = 0;

  for (i in 1:N_items) {
    for (j in 1:N_features) {
      int b = block_id[j];
      real lp   = alpha_b[b] + beta_b[b] * eta[i] + u_feat[j] + v_item[i];
      real p_ij = inv_logit(lp);
      real q_ij = (1 - gminus[j]) * p_ij + gplus[j] * (1 - p_ij);
      row_sum_exp[i] += q_ij;
      col_sum_exp[j] += q_ij;
    }
  }

  // Masked log-lik
  for (k in 1:N_mask) {
    int i = idx_i[k];
    int j = idx_j[k];
    int b = block_id[j];
    real lp   = alpha_b[b] + beta_b[b] * eta[i] + u_feat[j] + v_item[i];
    real p_ij = inv_logit(lp);
    real q_ij = (1 - gminus[j]) * p_ij + gplus[j] * (1 - p_ij);
    log_lik_mask[k] = bernoulli_lpmf(y[i, j] | q_ij);
  }
}
