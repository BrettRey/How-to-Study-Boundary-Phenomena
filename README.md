# reciprocals2 — measurement-first analysis & robustness checks

This repo holds the analysis pipeline for English reciprocals. It is measurement-first and transparency-first: one binary item×feature matrix is interrogated via ordination, supervised calibration, permutation with preserved margins, specification curves, posterior predictive checks (PPCs), matched-subset robustness, and a small fake-data recovery/calibration study.

## Repo layout

- `data/` — input matrices (e.g. `matrix_clean.csv`, `feature_blocks_by_order.csv`)
  - `inventory_annotations.csv` — override layer for follow-on pronoun/determinative work; keeps seed class, revised class, personhood profile, and anchor notes separate from the original matrix
  - `anchor_sets.csv` — derived first-pass anchor buckets (`clear_pronoun_anchor`, `clear_determinative_anchor`, `review_item`, `excluded_item`)
  - `clear_pronoun_anchors.txt`, `clear_determinative_anchors.txt` — lemma manifests for conservative first-pass anchor sets
  - `anchor_review_items.txt`, `anchor_excluded_items.txt` — items held out from clear anchors
- `code/` — R scripts
  - `build_inventory_annotations.py` — generates `data/inventory_annotations.csv` from `matrix_clean.csv` plus paper-based overrides
  - `build_anchor_sets.py` — derives concrete anchor files from `inventory_annotations.csv`
  - `06_pron_det_anchor_baseline.R` — anchor-only pronoun vs determinative baseline; fits a ridge classifier to clear anchors, scores review items, and writes a Jaccard MDS plot
- `plots/` — generated figures (PNG)
  - `pron_det_anchor_clustermap.png` — clustered Jaccard map for the current pronoun vs determinative baseline
  - `pron_det_anchor_heatmap.png` — anchor-plus-review Jaccard heat map for the current pronoun vs determinative baseline
  - `pron_det_anchor_mds.png` — anchor/review ordination for the first pronoun vs determinative pass
  - `pron_det_review_focus_clustermap.png` — labeled cluster map for review items plus their nearest pronoun and determinative anchors
- Root `.txt` outputs — small, human-readable artefacts that back figures and appendix statements:
  - `matched_subset_manifest.txt` — canonical 6+6 item lists (predeclared)
  - `matched_subset_robustness.txt` — one-line appendix summary of canonical p and rotation distribution
  - `ppc_summary.txt` — quantiles + Pr(anchor overlap > reciprocal overlap)
  - `pron_det_anchor_summary.txt` — anchor counts, cross-validated baseline performance, and review-item scores
  - `sessionInfo.txt` — R environment snapshot
  - `stan_fit_summary.txt` — posterior table
  - `stan_version.txt` — CmdStan version
  - `weight_calibration.txt` — simulator calibration summaries
  - `data/pron_det_anchor_scores.csv` — per-item baseline scores, nearest anchors, and MDS coordinates

## Environment

- R 4.5.x (see `sessionInfo.txt` for exact snapshot)
- Packages: `tidyverse`, `Matrix`, `proxy`, `vegan`, `glmnet`, `cmdstanr`, `posterior`
- CmdStan version: `2.36.0`a

Install CmdStan if needed:
```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
writeLines(cmdstanr::cmdstan_version(), "stan_version.txt")
