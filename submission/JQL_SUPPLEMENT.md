# JQL Supplement Package Draft

This submission should be accompanied by a neutral supplement, separate from the CJL provenance package.

## Contents To Include

- `data/matrix_clean.csv` - item-by-feature matrix.
- `data/feature_blocks_by_order.csv` - feature block mapping.
- `data/inventory_annotations.csv` - inventory annotation layer.
- `data/anchor_sets.csv` - derived anchor buckets.
- `data/clear_pronoun_anchors.txt` and `data/clear_determinative_anchors.txt` - anchor manifests.
- `data/anchor_review_items.txt` and `data/anchor_excluded_items.txt` - held-out and excluded item manifests.
- `matched_subset_manifest.txt` - prespecified comparator items.
- `matched_subset_robustness.txt`, `ppc_summary.txt`, `weight_calibration.txt`, `stan_fit_summary.txt`, and `sessionInfo.txt` - human-readable analysis summaries.
- `code/` - R, Python, and Stan scripts needed to reproduce the analyses.
- `plots/` - generated figure files used in the manuscript.

## Build Notes

The JQL manuscript wrapper is `main-jql.tex`. The historical CJL wrapper and anonymous CJL source package should remain untouched.

Compile with:

```bash
pdflatex main-jql.tex
biber main-jql
pdflatex main-jql.tex
pdflatex main-jql.tex
```

## Submission Note

If the Taylor & Francis system asks for data availability wording, use:

> The data matrix, feature-block mapping, anchor manifests, analysis scripts, generated figures, and derived robustness summaries are supplied as supplementary material.

