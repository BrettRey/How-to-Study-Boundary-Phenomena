# reciprocals2 — measurement-first analysis & robustness checks

This repo holds the analysis pipeline for English reciprocals. It is measurement-first and transparency-first: one binary item×feature matrix is interrogated via ordination, supervised calibration, permutation with preserved margins, specification curves, posterior predictive checks (PPCs), matched-subset robustness, and a small fake-data recovery/calibration study.

## Repo layout

- `data/` — input matrices (e.g. `matrix_clean.csv`, `feature_blocks_by_order.csv`)
- `code/` — R scripts
- `plots/` — generated figures (PNG)
- Root `.txt` outputs — small, human-readable artefacts that back figures and appendix statements:
  - `matched_subset_manifest.txt` — canonical 6+6 item lists (predeclared)
  - `matched_subset_robustness.txt` — one-line appendix summary of canonical p and rotation distribution
  - `ppc_summary.txt` — quantiles + Pr(anchor overlap > reciprocal overlap)
  - `sessionInfo.txt` — R environment snapshot
  - `stan_fit_summary.txt` — posterior table
  - `stan_version.txt` — CmdStan version
  - `weight_calibration.txt` — simulator calibration summaries

## Environment

- R 4.5.x (see `sessionInfo.txt` for exact snapshot)
- Packages: `tidyverse`, `Matrix`, `proxy`, `vegan`, `glmnet`, `cmdstanr`, `posterior`
- CmdStan version: `2.36.0`a

Install CmdStan if needed:
```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
writeLines(cmdstanr::cmdstan_version(), "stan_version.txt")
