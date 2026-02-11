data {
  int<lower=1> N_items;
  int<lower=1> N_features;
  int<lower=1> N_blocks;

  array[N_features] int<lower=1, upper=N_blocks> block_id;

  array[N_items, N_features] int<lower=0, upper=1> y_obs;
  array[N_items, N_features] int<lower=0, upper=1> is_obs;
  array[N_items, N_features] int<lower=0, upper=1> y_holdout;

  int<lower=0> K_pron_anchors;
  array[K_pron_anchors] int<lower=1, upper=N_items> pron_anchor_idx;

  int<lower=0> K_det_anchors;
  array[K_det_anchors] int<lower=1, upper=N_items> det_anchor_idx;

  real<lower=0> eta_anchor_sd;
}
parameters {
  vector<lower=0, upper=1>[N_items] eta;

  vector[N_blocks] alpha;
  vector[N_blocks] beta;

  vector[N_features] u;     // feature difficulties
  vector[N_items]    v_raw; // item heterogeneity

  real<lower=0> tau_f;
  real<lower=0> tau_i;

  vector<lower=0, upper=1>[N_blocks] gamma_pos; // FP
  vector<lower=0, upper=1>[N_blocks] gamma_neg; // FN
}
transformed parameters {
  vector[N_items] v = tau_i * v_raw;
}
model {
  // Priors
  alpha ~ normal(0, 2);
  beta  ~ normal(0, 2);
  tau_f ~ normal(0, 1);
  tau_i ~ normal(0, 1);
  u     ~ normal(0, tau_f);
  v_raw ~ normal(0, 1);
  eta   ~ beta(1, 1);

  // Anchor priors
  for (k in 1:K_pron_anchors)
    target += normal_lpdf(eta[pron_anchor_idx[k]] | 1, eta_anchor_sd);
  for (k in 1:K_det_anchors)
    target += normal_lpdf(eta[det_anchor_idx[k]] | 0, eta_anchor_sd);

  gamma_pos ~ beta(2, 20);
  gamma_neg ~ beta(2, 20);

  // Likelihood with measurement error
  for (i in 1:N_items) {
    for (j in 1:N_features) {
      if (is_obs[i, j] == 1) {
        int b = block_id[j];
        real lin = alpha[b] + beta[b] * eta[i] + u[j] + v[i];
        real p   = inv_logit(lin);                          // P(x=1)
        real q   = (1 - gamma_neg[b]) * p + gamma_pos[b] * (1 - p); // P(y=1)
        y_obs[i, j] ~ bernoulli(q);
      }
    }
  }
}
generated quantities {
  array[N_items, N_features] int y_rep;
  real lpd_mask = 0;

  for (i in 1:N_items) {
    for (j in 1:N_features) {
      int b = block_id[j];
      real lin = alpha[b] + beta[b] * eta[i] + u[j] + v[i];
      real p   = inv_logit(lin);
      real q   = (1 - gamma_neg[b]) * p + gamma_pos[b] * (1 - p);

      y_rep[i, j] = bernoulli_rng(q);

      if (is_obs[i, j] == 0)
        lpd_mask += bernoulli_lpmf(y_holdout[i, j] | q);
    }
  }
}
