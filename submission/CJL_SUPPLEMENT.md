# CJL Supplement Package

This submission is accompanied by an anonymized supplement containing the files below.

## Core Measurement Files

- `data/matrix_clean.csv` — full item-by-feature matrix used in the manuscript
- `data/feature_blocks_by_order.csv` — feature-to-block map used for blockwise analyses
- `data/inventory_annotations.csv` — item-level inventory annotations and anchor notes
- `data/anchor_sets.csv` — derived anchor-bucket assignments
- `data/clear_pronoun_anchors.txt` — conservative pronoun anchor manifest
- `data/clear_determinative_anchors.txt` — conservative compound-determinative anchor manifest
- `data/anchor_review_items.txt` — items held out from clear anchors pending review
- `data/anchor_excluded_items.txt` — items excluded from clear anchors

## Analysis Outputs

- `matched_subset_manifest.txt` — canonical 6+6 comparison-set manifest
- `data/matched_subset_robustness.csv` — observed `Δ` summaries for the canonical set and 100 rotations
- `matched_subset_robustness.txt` — human-readable summary of the rotation distribution and quasiswap tail-area summary
- `data/quasiswap_reference_delta.csv` — canonical quasiswap reference draws for `Δ`
- `data/weight_calibration_curve.csv` — predictive-fit curve across mixture weights
- `weight_calibration.txt` — predictive-blend minima for the two reciprocals

## Key Analysis Scripts

- `code/03_matched_set_robustness.R`
- `code/04_weight_calibration.R`
- `code/06_pron_det_anchor_baseline.R`

## Notes

- The supplement is anonymized for review and omits repository history plus direct author identifiers.
- The quasiswap reference draws were generated with `B = 5000`, `burn-in = 1500`, and `set.seed(2025)` in `code/03_matched_set_robustness.R`.
