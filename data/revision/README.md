# Revised reciprocal analysis

This directory supports *Calibrating Diagnostic Conflict: English Reciprocals in Grammatical Feature Space*. The main analysis is saved in `results/`, the exploratory equal-block-total comparison in `followup/`, and the later deterministic diagnostics in `posthoc/`. The dated records distinguish these stages; the commands below reproduce them from the supplied matrix.

## Reproduce

From the repository or supplementary-package root:

```bash
python3 code/revision_analysis.py
python3 code/revision_tables.py
python3 code/revision_followup.py
python3 code/revision_review_figures.py
python3 code/revision_posthoc.py
cd code
python3 -m unittest -v test_revision_analysis.py test_revision_followup.py
```

The analysis requires NumPy, pandas, and Matplotlib; the tests additionally require SciPy. Recorded versions are in `environment.json` and `results/run_metadata.json`. The full run uses 20,000 randomized profiles for each item, retained-feature specification, and reference scheme. A smoke test can use `--draws 100 --output /path/to/scratch`; it records the draw override and is not the manuscript analysis. Running the full commands replaces the derived results, never the original inputs.

Optional workspace plot styling is used when available; otherwise the script uses Matplotlib's serif defaults. Fonts and figure bytes can vary across environments. Numerical CSV files, their definitions, and the recorded random seeds are the reproducibility targets. Manuscript distances, contrasts, and resampling ranges are displayed to two decimal places; the CSV files retain full precision. Contrasts are calculated before rounding. Displayed digits do not measure the reliability of the grammatical coding.

## Inputs and design

- `../matrix_clean.csv`: unchanged 138-item, 155-feature archived extract, plus item identifier and inherited category label. Rows are coded forms or readings, not 138 distinct lexemes.
- `../feature_blocks_by_order.csv`: unchanged feature-to-block mapping.
- `comparators.csv`: one row per original item, with its analytical role, family, inclusion/exclusion rationale, source, and locative-compound indicator. The primary groups have 41 core personal-pronoun forms/readings and 16 central compound determinatives. Both reciprocal targets are excluded from all comparator pools.
- `protocol.json`: choices recorded before the corrected run, including input hashes, seed, reference-draw count, filtering, holdouts, measures, and sensitivity design. Earlier exploratory results were known when this protocol was written; it is not an independent preregistered replication.

The current reproduction starts from the archived 155-feature extract; it does not reconstruct the earlier reduction from the 232-feature source matrix. The present 155→154 filter is applied directly to the supplied archive.

Filtering removes the singleton identity feature `each_other`, leaving 154 features: 65 morphological, 3 phonological, 36 semantic, and 50 syntactic. Filtering uses the fixed inventory, not a fitted training-data selection procedure.

The personal-pronoun forms are grouped into seven analysis families, pooling the recorded singular/plural uses of *you* and *they*. Compound families share the first component (*every*, *some*, *any*, *no*). Family holdout removes all forms of that family. It addresses one source of dependence, not every dependency among the remaining items.

## Outputs

| File in `results/` | Contents |
|---|---|
| `all_specifications.csv` | 9,440 item/specification rows: 59 items × 5 feature choices × 2 references × 2 comparator pools × 2 holdouts × 4 measures. Repeated rows for the same item are not independent observations. |
| `primary_item_results.csv` | All 59 items under full-pool, all-feature, item-held-out Jaccard with fixed-total randomization. |
| `control_recovery_by_specification.csv` | Directional control recovery, nearest control-range endpoints, and whether both targets lie between those endpoints. This is a descriptive comparison, not a boundary classifier. |
| `target_reference_draws.csv.gz` | Jaccard distances and contrasts for 20,000 draws per reciprocal under each reference scheme (80,000 rows). Other measures' summaries are saved in the specification file; their complete draws can be regenerated. |
| `coding_sensitivity.csv` | Each target under each single-feature deletion, global polarity reversal, and target-cell flip, for all four measures (3,696 rows). |
| `comparator_resampling.csv` | 2,000 family resamples per target at the original number of families. These describe composition sensitivity, not population uncertainty. |
| `feature_dictionary.csv`, `feature_filter.json` | Retained identifiers, block assignments, positive counts, and the filtering record. |
| `illustrative_matrix.csv` | Actual source cells used in manuscript Table 1. |
| `figures/` | PDF and PNG copies of Figures 1 and 2, plus an earlier rendering of Figure 3. The manuscript uses Figure 3 from `review_corrections/`, described below. |
| `numbers.tex`, `*_table.tex` | Generated manuscript numbers and table bodies. Table bodies include their closing `bottomrule`. |
| `run_metadata.json`, `table_metadata.json` | Input, protocol, script, and table provenance hashes. |

In the specification file, `d_pronoun` and `d_compound` are equal-form-weighted mean distances to the remaining comparator pools; `delta` is their difference. Positive values indicate greater proximity to the compound group. The `*_percentile` fields give empirical midranks against target randomization, not category-membership probabilities. The `*_lower_tail` fields use `(1 + number of reference draws no greater than observed)/(B + 1)`, with absolute numerical tolerance `1e-12`.

`fixed_total` preserves the target's total number of presences; `within_block` preserves that total in each retained feature block. Anchors stay fixed in both. These artificial references destroy feature dependencies and are not models of possible English words. IDF-Jaccard weights are estimated from the remaining anchors separately in every fold and then held fixed for the observed target and its random counterparts.

The nonlocative pool has 12 compound anchors. Excluded locatives remain probes in the output, with `eligible_control=False`; they do not enter that specification's control-recovery denominator.

## Coding notes and current Figure 3

[CODING_NOTES.md](CODING_NOTES.md) explains the feature definitions and their limits, source scope, target identity outside morphology, control-range endpoints, and matrix provenance.

`code/revision_review_figures.py` reads the saved specification results and writes Figure 3 to `review_corrections/`. The directory contains the PDF/PNG, the 40 plotted values, and script/input/output metadata. The legend identifies IDF-weighted Jaccard; the caption specifies the full pools, item holdout, and fixed-total placement. The earlier rendering remains in `results/figures/`. This script performs no numerical analysis or randomization.

## Interpretation

The post-review addendum is documented in [FOLLOWUP.md](FOLLOWUP.md), with its own dated protocol and outputs in `followup/`. It adds one equal-block-total weighting comparison and extracts the complete original grid, control margins, exact Monte Carlo counts, standard errors, and zero-hit bounds from saved results. It generates no new random draws and leaves `results/` unchanged.

With all features, reciprocal contrasts fall between the observed control ranges. Unweighted Jaccard recovers all 57 controls under both holdouts. Each reciprocal is nearer each group than every member of the other group is, but farther from it than its own members, apart from dummy it. Low random-placement tails are common among controls as well: 50 of 57 have lower tails below 0.05 against their non-designated group.

Removing morphology reverses the contrasts. Within-pronoun-range placements depend on the dummy-it endpoint; every control failure across the original grid and recorded block-weight settings is also dummy it. The positive equal-block-total contrast depends on two heavily weighted phonological columns. Balancing only semantics and syntax strengthens the pronoun preference. Coordinated morphological deletions locate the dependence in whole-word codes rather than component identities, while the meaning of `does_not_inflect` remains uncertain.

[POSTHOC.md](POSTHOC.md) gives the new deterministic diagnostics and feature glosses. `posthoc_protocol.json` records their scope retrospectively, and `code/revision_posthoc.py` writes separately to `posthoc/`. No random draws or original outputs are changed. The supported claim remains conditional on feature selection, weighting, comparator design, and inherited coding.
