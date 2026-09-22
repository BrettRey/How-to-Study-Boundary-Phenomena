# Post-review weighting comparison and reporting details

This exploratory addendum supplements the original analysis without changing its inputs, scripts, or saved results. `followup_protocol.json` records the scope and input hashes before the added calculation. The comparison was selected in response to review; it is not an independently preregistered test.

Run `python3 code/revision_followup.py` from the repository or supplementary-package root after the original results are available. It reads the archived matrix, block map, comparator manifest, original protocol, and original results, and writes only to `data/revision/followup/`. It generates no new random draws. From `code/`, run `python3 -m unittest -v test_revision_followup.py` for the focused tests.

## Equal total feature weight per block

For each retained feature in block b, the weight is 1 divided by the number of retained features in b. Each included block consequently has total feature weight one. Pairwise weighted Jaccard is the sum of weighted mismatches divided by the sum of weighted presences in the union. The implementation returns zero for an empty weighted union. This ratio is not an average of the block-specific distances, and blocks need not contribute equally to a particular pair's union denominator.

The all-feature weights are 1/65 for morphology, 1/3 for phonology, 1/36 for semantics, and 1/50 for syntax. Removing morphology retains the other three weights. Giving three phonological features the same total weight as a much larger block is an alternative specification, not an assertion of linguistic optimality. The later diagnostic in [POSTHOC.md](POSTHOC.md) shows that two phonological columns produce the positive contrasts. Equalizing nominal block totals doesn't isolate block size from coded content, dependence, or validity.

The scope is the full 41-personal-pronoun and 16-compound pools, item and family holdouts, and all features or morphology removed. Both targets are excluded from every pool; controls are evaluated after removing themselves or their entire family. The 59 evaluated items across four settings yield 236 new rows. The comparison file adds their 236 matching original unweighted-Jaccard rows.

## Findings

| Retained features | Feature weighting | Each other contrast | One another contrast | Controls: item / family |
|---|---|---:|---:|---:|
| All | Equal feature weights | 0.04 | 0.05 | 57/57 / 57/57 |
| All | Equal block totals | 0.07 | 0.07 | 56/57 / 56/57 |
| No morphology | Equal feature weights | −0.08 | −0.08 | 57/57 / 57/57 |
| No morphology | Equal block totals | 0.02 | 0.02 | 56/57 / 56/57 |

Positive contrasts mean closer average proximity to compounds. Contrasts are calculated before rounding; both absolute distances and full precision are in the CSV files. With all features, the targets fall between the control ranges under either weighting and holdout. With morphology removed and equal block totals, they remain between the ranges under item holdout but fall within the pronoun range under family holdout. The sole control failure in each reweighted setting is dummy `it_dum`. Its contrast also sets the upper pronoun-range endpoint. In the no-morphology comparison, that endpoint passes the unchanged target contrasts when control holdout changes; every other pronoun control remains negative. See [CODING_NOTES.md](CODING_NOTES.md) for the endpoints and related representation clarifications.

The positive block-total contrasts depend on the phonological block. Removing morphology and phonology gives both targets Δ = −0.112700 when semantics and syntax have equal total weight; removing phonology alone gives −0.018615 and −0.011833. Both checks recover all 57 controls under both holdouts. Balancing the large non-morphological blocks strengthens the pronoun preference. These additional diagnostics are explicitly post hoc; their outputs are separate from the recorded weighting comparison.

## Post hoc diagnostics

[POSTHOC.md](POSTHOC.md) documents the phonology exclusions, coordinated morphological deletions, dummy-it subject-cell sensitivity, and calibration of the benchmark against saved control rows. Run `python3 code/revision_posthoc.py`. The retrospective reporting record is `posthoc_protocol.json`, and outputs are in `posthoc/`; these don't overwrite `followup/` or `results/` and generate no random draws.

## Output dictionary

All files below are in `followup/`.

| File | Contents |
|---|---|
| `feature_weights.csv` | 243 feature weights: 154 all-feature and 89 no-morphology rows. |
| `block_weight_item_results.csv` | 236 new item/specification rows, with comparator counts, both distances, contrasts, and control margins. |
| `weighting_comparison.csv` | 472 rows: the new results and their matching original equal-feature results. |
| `block_weight_summary.csv` | Eight feature/holdout/weighting summaries, with recovery, control-range endpoints, and each target's location. |
| `complete_specification_grid.csv` | All 160 original pool/holdout/feature/metric/reference settings, with recovery, target contrasts and percentiles, and control ranges. No new analysis settings. |
| `original_control_margins.csv` | All 4,400 eligible-control rows across the 80 original observed-distance specifications, without duplicating distances for the two reference benchmarks. |
| `original_control_margin_summary.csv` | 160 group-specific count/minimum/median/maximum summaries of those signed margins. |
| `target_monte_carlo_details.csv` | Twelve original primary-Jaccard summaries: two targets × two benchmarks × three quantities, computed from the saved reference draws. |
| `numbers.tex`, `block_weight_table.tex` | Generated manuscript numbers and table; the table fragment includes its closing `bottomrule`. |
| `metadata.json` | Script, protocol, input, and output hashes, package versions, and row/draw counts. |

`d_pronoun` and `d_compound` are mean distances over the remaining forms. `delta` is the first minus the second. `signed_control_margin` is −delta for pronoun controls and delta for compound controls, so positive means correct directional recovery. Target margins are blank. `p_min/p_max` and `c_min/c_max` describe raw control contrasts, not signed margins. Target range flags use inclusive endpoints; being between the ranges uses the open gap between the highest pronoun contrast and lowest compound contrast. The `drop` value `none` means all retained features; `morph` means morphology removed.

## Monte Carlo quantities and uncertainty

For each saved set of B = 20,000 draws, numerical ties use absolute tolerance 1e-12. The percentile is `(strictly_below + tied/2)/B`. Its Monte Carlo standard error is the sample standard deviation (denominator B−1) of scores 1, 0.5, and 0, divided by √B. The lower-tail hit count includes ties, and the smoothed lower tail is `(1 + lower_tail_hits)/(B + 1)`.

The primary contrast percentiles are 0.900600 and 0.928550, with Monte Carlo standard errors about 0.0021 and 0.0018. Both distances for both targets have zero lower-tail hits under both Jaccard benchmarks. For zero hits, the exact one-sided 95% upper bound on the placement probability is `1 − 0.05^(1/B)`, or 0.000149775396223. Zero observed variance in the hit indicator does not establish zero placement probability.

These errors and bounds concern approximation to the specified artificial benchmark. Family resampling describes comparator-composition sensitivity; feature changes describe sensitivity to coding or weighting choices. None estimates independent coding reliability, population variation among lexical types, or probability of grammatical category membership.
