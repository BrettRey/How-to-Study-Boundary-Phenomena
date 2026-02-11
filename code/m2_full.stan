data {
  int<lower=1> N_items;
  int<lower=1> N_features;
  int<lower=1> N_blocks;
  array[N_features] int<lower=1, upper=N_blocks> block_id;
  array[N_items, N_features] int<lower=0, upper=1> y_obs;

  int<lower=0> N_mask;
  array[N_items, N_features] int<lower=0, upper=1> mask;
  array[N_mask] int<lower=1, upper=N_items> idx_i;
  array[N_mask] int<lower=1, upper=N_features> idx_j;

  int<lower=0> N_pron_anchor;
  int<lower=0> N_det_anchor;
  array[N_pron_anchor] int<lower=1, upper=N_items> pron_anchor_idx;
  array[N_det_anchor] int<lower=1, upper=N_items> det_anchor_idx;
  real<lower=0> sd_anchor;
}
parameters {
  vector[N_blocks] alpha_b;
  vector[N_blocks] beta_b;

  vector[N_features] u_feat;
  real<lower=0> sigma_f;
  vector[N_items] v_item;
  real<lower=0> sigma_i;

  vector[N_blocks] mu_gplus_b;
  vector<lower=0>[N_blocks] sd_gplus_b;
  vector[N_blocks] mu_gminus_b;
  vector<lower=0>[N_blocks] sd_gminus_b;
  vector[N_features] gplus_raw;
  vector[N_features] gminus_raw;

  vector[N_items] z;
}
transformed parameters {
  vector[N_features] gplus;
  vector[N_features] gminus;
  for (j in 1:N_features) {
    int b = block_id[j];
    gplus[j]  = inv_logit(mu_gplus_b[b]  + sd_gplus_b[b]  * gplus_raw[j]);
    gminus[j] = inv_logit(mu_gminus_b[b] + sd_gminus_b[b] * gminus_raw[j]);
  }
}
model {
  alpha_b ~ normal(0, 1.5);
  beta_b  ~ normal(0, 1.5);
  u_feat  ~ normal(0, sigma_f);
  sigma_f ~ normal(0, 1) T[0,];
  v_item  ~ normal(0, sigma_i);
  sigma_i ~ normal(0, 1) T[0,];

  mu_gplus_b  ~ normal(logit(0.02), 1.0);
  mu_gminus_b ~ normal(logit(0.02), 1.0);
  sd_gplus_b  ~ normal(0, 0.5) T[0,];
  sd_gminus_b ~ normal(0, 0.5) T[0,];
  gplus_raw   ~ std_normal();
  gminus_raw  ~ std_normal();

  if (N_pron_anchor > 0) z[pron_anchor_idx] ~ normal(+3, sd_anchor);
  if (N_det_anchor  > 0) z[det_anchor_idx]  ~ normal(-3, sd_anchor);
  for (i in 1:N_items) target += normal_lpdf(z[i] | 0, 1);

  for (i in 1:N_items) {
    real eta_i = inv_logit(z[i]);
    for (j in 1:N_features) {
      int b = block_id[j];
      real lp = alpha_b[b] + beta_b[b] * eta_i + u_feat[j] + v_item[i];
      real p  = inv_logit(lp);
      real q  = (1 - gminus[j]) * p + gplus[j] * (1 - p);
      y_obs[i, j] ~ bernoulli(q);
    }
  }
}
generated quantities {
  matrix[N_items, N_features] q;
  vector[N_items] row_sum_exp;
  vector[N_features] col_sum_exp;
  vector[N_mask] log_lik_mask;
  vector[N_mask] q_mask_mean;
  vector[N_items] eta;

  for (i in 1:N_items) {
    eta[i] = inv_logit(z[i]);
    row_sum_exp[i] = 0;
  }
  for (j in 1:N_features) col_sum_exp[j] = 0;

  for (i in 1:N_items) {
    for (j in 1:N_features) {
      int b = block_id[j];
      real lp = alpha_b[b] + beta_b[b] * eta[i] + u_feat[j] + v_item[i];
      real p  = inv_logit(lp);
      real qij = (1 - gminus[j]) * p + gplus[j] * (1 - p);
      q[i, j] = qij;
      row_sum_exp[i] += qij;
      col_sum_exp[j] += qij;
    }
  }
  for (k in 1:N_mask) {
    int i = idx_i[k];
    int j = idx_j[k];
    log_lik_mask[k] = bernoulli_lpmf(y_obs[i, j] | q[i, j]);
    q_mask_mean[k]  = q[i, j];
  }
}
