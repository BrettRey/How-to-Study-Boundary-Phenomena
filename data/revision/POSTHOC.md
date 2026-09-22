# Post hoc diagnostics of weighting, morphology, and control calibration

These deterministic checks were identified after independent review of the revised manuscript. They were computed before `posthoc_protocol.json` was recorded; that file is a retrospective reporting record, not a prospective protocol or preregistration. They add no random draws and preserve the original matrix, original scripts and results, and recorded equal-block-total comparison.

Run `python3 code/revision_posthoc.py` from the repository or supplementary-package root. The script reads the unchanged inputs and saved results, verifies their hashes, and writes only to `data/revision/posthoc/`. It implements weighted mismatch/union distances directly, without importing the original analysis helpers. Its four baseline comparisons reproduce 472 saved item rows to floating-point precision.

## Weighting and phonology

Full pools are used throughout. Both targets are excluded from all anchors; each control is evaluated with its item or family held out. Equal block totals assign each feature weight `1 / retained features in its block`, recomputed after each deletion. The distance is the weighted mismatch sum divided by the weighted union sum, with zero for an empty union.

| Features | Weights | Each other Δ | One another Δ | Controls, item / family |
|---|---|---:|---:|---|
| No morphology | Equal block totals | 0.018662 | 0.018662 | 56/57 / 56/57 |
| No morphology or phonology | Equal block totals | −0.112700 | −0.112700 | 57/57 / 57/57 |
| No morphology or particle-stress column | Block totals re-equalized | −0.060271 | −0.060271 | 57/57 / 57/57 |
| All except phonology | Equal block totals | −0.018615 | −0.011833 | 57/57 / 57/57 |
| No morphology | Equal feature weights | −0.076409 | −0.076409 | 57/57 / 57/57 |
| No morphology or phonology | Equal feature weights | −0.091319 | −0.091319 | 57/57 / 57/57 |

`Must_be_stressed_as_object_after_particle` has 27 positive entries among the 41 personal-pronoun controls; `Start_with_th` has 10. Both are zero for every compound and target. `Start_with_hw` is zero for all anchors and targets. With all three phonological columns retained, each has weight 1/3, compared with 1/36 for semantics and 1/50 for syntax. This concentrates weight on mismatches with pronouns. Balancing semantics and syntax strengthens the negative contrast; the positive block-total contrast depends on the phonological block.

Haas's fixed internal-stress argument for *each other* concerns a different construction from the particle-stress column. It doesn't independently validate that column's coding.

## Coordinated morphological deletions

These are full-pool Jaccard comparisons with equal feature weights. Target distances use the same pools under either control holdout.

| Columns deleted | Each other Δ | One another Δ |
|---|---:|---:|
| All 49 component-identity columns | 0.056006 | 0.056006 |
| Other 16 morphological columns | −0.073562 | −0.058472 |
| `compound_word`, `does_not_inflect`, `has_derivationtionally_related_words` | −0.028361 | −0.018792 |
| Remaining 13 morphological columns | 0.002409 | 0.012440 |

The three named columns are positive for both targets and all 16 compounds, compared with 9, 1, and 12 personal pronouns, respectively. Deleting them together changes the sign; deleting component identities strengthens the compound preference. These effects include changes to the Jaccard denominator and aren't additive contributions.

The exact 49-, 16-, 3-, and 13-column sets are in `posthoc_protocol.json`; the script checks their membership. The remaining 13 include `Monmorphemic`, `in_word_formation_e_g_a_many_sided_shape`, and `has_corresponding_word`, alongside case and paradigm properties.

## Influential labels and their limits

These glosses identify the source labels. They don't supply missing operational definitions or replace `results/feature_dictionary.csv`.

| Feature | Gloss and evidence | Definition limit |
|---|---|---|
| `compound_word` | Coded compound wordhood. The reciprocal values are consistent with CGEL Ch. 5 §10.1.2, which treats the two orthographic words as one grammatical word. Haas (2007), p. 42, treats fixed internal stress in *each other* as lexicalization evidence. | The cited evidence supports reciprocal wordhood; operational rules for coding the full inventory remain unspecified. |
| `has_derivationtionally_related_words` | Judgement that derivationally related forms exist. Identifier spelling is preserved. | No archived decision rule identifies which relations count. |
| `inflects_for_case` | Case-inflection judgement. CGEL documents reciprocal genitives. | A general label is not a complete rule for coding every form or paradigm. |
| `does_not_inflect` | Negative inflection judgement. Both targets and 12 compounds also have 1 for case inflection. | The original definition cannot be recovered. Haas p. 48 n. 32 distinguishes absence of person, number, and gender inflection, but adopting that as the historical code definition would be retrospective. |
| `paradigm_has_distinct_acc_form` | The item's paradigm has a distinct accusative form; reflexive *herself* receives 1. | This is a paradigm-level property, not a claim that the row is itself accusative. |
| `appears_in_subject` | Subject-use judgement. | The source assigns 0 to reciprocals and dummy pronouns despite relevant subject constructions; a uniform constructional scope isn't documented. |

## Dummy-it endpoints

Dummy `it_dum` is the highest pronoun contrast in all 80 original observed specifications and all four recorded block-total settings. Every control failure in those settings is dummy `it_dum`; every other eligible control is recovered. This statement covers the declared grid, not arbitrary new specifications.

| Jaccard setting | Targets Δ | Dummy it Δ | Next-highest pronoun Δ |
|---|---:|---:|---:|
| All features, item holdout | 0.039776 / 0.045544 | −0.043 | −0.132 |
| No morphology, item holdout | −0.076409 / −0.076409 | −0.074579 | −0.167250 |
| No morphology, family holdout | −0.076409 / −0.076409 | −0.066928 | −0.159808 |

Without morphology, full-pool Jaccard, IDF-weighted Jaccard, and Hamming place targets inside the range only because it extends to this one control. Dice leaves them above it. With nonlocative compounds, Jaccard, Dice, and Hamming leave them between the ranges under either holdout. Full precision and next-control identities for all 80 settings are in `posthoc/original_grid_endpoints.csv`.

For the separate subject-cell intervention, set dummy `it_dum` to 1 on `appears_in_subject` and `Functions_as_subject_of_interrogative_tag` in a copied matrix and recompute every evaluated item. This updates dummy it as both a control and an anchor. Its item-held-out contrast changes to −0.036376 with all features and −0.065304 without morphology. It remains the highest pronoun contrast. Equal-feature placement descriptions are unchanged. Under equal block totals without morphology, item-held-out placement becomes within-range, while the target contrasts change slightly because one anchor was recoded. The archived cells remain unchanged.

## Calibration of the random-placement benchmark

The following summaries use saved primary Jaccard, fixed-total, item-held-out control rows. No new reference draws are generated.

| Controls | Distance to non-designated group: median lower tail | Maximum | Below 0.05 |
|---|---:|---:|---:|
| Personal pronouns → compounds | 0.001850 | 0.102145 | 35/41 |
| Compounds → personal pronouns | 0.000400 | 0.052997 | 15/16 |

Five controls have zero hits against their non-designated group: `everyone`, `hers`, `itself`, `theirs_sing`, and `yours_plur`. Thus 50 of 57 controls satisfy the same low-tail criterion used for the reciprocals against the group they don't belong to. The benchmark describes departure from random placement; it provides little discrimination between these targets and clear controls.

The observed distances provide a more informative comparison. Both reciprocals are nearer each group than every member of the other group is. They are farther from each group than all of its own controls, except dummy it. `control_benchmark.csv` supplies the extrema: the minimum compound-control distance to pronouns is 0.733661, and the minimum pronoun-control distance to compounds is 0.722245; maxima of own-group distances are 0.607794 for pronouns excluding dummy it and 0.570149 for compounds. Target distances retain their original values.

## Output dictionary

| File in `posthoc/` | Contents |
|---|---|
| `deterministic_item_results.csv` | 1,888 item rows: 12 baseline/deletion comparisons and four subject-cell variants, each with both holdouts and all 59 items. |
| `diagnostic_summary.csv` | 32 setting summaries with recovery, failures, endpoints, next-highest pronoun controls, and target values. |
| `original_grid_endpoints.csv` | Endpoints, next-highest controls, failures, and target contrasts for the 80 original observed settings. |
| `feature_counts.csv` | Positive counts for the three shared morphological codes, three phonological codes, and case inflection. |
| `control_benchmark.csv` | Saved control lower-tail summaries and observed-distance extrema. |
| `checks.json` | Structured full-precision results, exact column sets, baseline comparisons, and input/script hashes. |
| `metadata.json` | Reporting-protocol, script, input and output hashes, package versions, and row/draw counts. |

The CSV summaries contain list-valued fields for targets and failures; `checks.json` provides the corresponding structured values. `altered_dummy_it=False` denotes unchanged coding. A positive Δ means nearer compounds. Control failure means Δ ≥ 0 for pronouns or Δ ≤ 0 for compounds.
