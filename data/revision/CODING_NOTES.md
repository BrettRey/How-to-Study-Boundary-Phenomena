# Coding scope and matrix provenance

These notes describe the source matrix, the scope of its coding, and the consequences for target comparison. All archived cells are retained in the baseline; the interventions reported here are sensitivity checks on separate copies.

## Subject position

The `appears_in_subject` column assigns 0 to both reciprocals. Haas (2007), §6.2, pp. 44–45, and Hurst and Nordlinger's author preprint, §2.1, report restricted embedded-subject uses of *each other*. The latter gives *They know what each other wants*, example (5), printed p. 5 / PDF p. 6.

The Hurst–Nordlinger source used here is the 2007 author preprint of the chapter published in 2011 (DOI 10.1075/tsl.98.04hur). Its printed page number is not the published chapter's page number. Haas is the published 2007 article (DOI 10.1017/S1360674306002103).

Without a more specific constructional definition, the inherited 0 does not establish an unrestricted prohibition. The sources do not establish that all speakers accept every subject use, nor do they supply a uniform new definition for the column. The archive and Table 1 preserve the legacy value and identify this scope limitation explicitly.

The existing `results/coding_sensitivity.csv` includes flipping this target cell under full-pool, full-feature Jaccard. The contrasts remain positive:

| Target | Baseline contrast | Contrast after target-cell flip |
|---|---:|---:|
| each other | 0.039776 | 0.056508 |
| one another | 0.045544 | 0.061321 |

These values summarize the previously computed mechanical perturbations. The positive shift does not make subject use intrinsically compound-like; it reflects this form-based comparison and its current feature profiles. Full precision remains in the source CSV files.

Dummy `it_dum` and `there` also have 0 on `appears_in_subject`; dummy `it_dum` has 0 on `Functions_as_subject_of_interrogative_tag`. These values conflict with ordinary dummy-it subject and tag use and with CGEL Ch. 5 §10.1.1. A copied-matrix intervention setting both dummy-it cells to 1 leaves it the highest pronoun contrast under primary and morphology-removed Jaccard. [POSTHOC.md](POSTHOC.md) reports this intervention and its effects under block weighting.

## Wordhood and inflection

The `compound_word` values are consistent with CGEL's single-grammatical-word analysis. Haas (2007), p. 42, treats the fixed internal stress of *each other* as lexicalization evidence; this concerns a different construction from the particle-stress column.

The original definition of `does_not_inflect` cannot be recovered. Both targets and 12 compounds also have 1 on `inflects_for_case`. Sources distinguish genitives from absent person, number, and gender inflection, but that distinction does not recover the historical coding rule. [POSTHOC.md](POSTHOC.md) gives the influential feature glosses, phonology exclusions, and coordinated morphological deletions.

## Identity outside morphology

After filtering, the targets differ only in `a`, `another`, `each`, and `one`, all morphological component columns. All 89 non-morphological entries are identical. Their identical observed results after morphology removal follow from the coded vectors; separately simulated percentile estimates can differ because the target-specific random draws differ.

## Control-range endpoints

Dummy it sets the upper pronoun endpoint in all 80 observed specifications and all four equal-block-total settings. Every control-recovery failure in those settings is that item. Within-pronoun-range placements depend on this endpoint in both the equal-feature morphology ablation and the block-weighted comparison. [POSTHOC.md](POSTHOC.md) reports the next-highest controls, the Dice and nonlocative exceptions, and the benchmark calibration.

### Equal block totals without morphology

With equal block totals and morphology removed, both target contrasts are 0.018662 under either holdout. The target comparator pools do not change. Dummy `it_dum` is the sole incorrectly recovered pronoun control and sets its group's upper range endpoint:

| Holdout | Dummy-it contrast | Largest other pronoun contrast |
|---|---:|---:|
| Item | 0.016605 | −0.070151 |
| Family | 0.020448 | −0.048641 |

The change from between-range to within-pronoun-range placement results from this endpoint crossing the unchanged target contrasts. It does not show that the targets move into the main part of the pronoun-control distribution. Values come from `followup/block_weight_item_results.csv`.

## Additional representation details

- `where` (morphology) and `locative` (semantics) are identical columns across the 138 rows. Equal block totals do not remove such redundancy.
- The seven personal-pronoun analysis families pool the recorded singular/plural uses of *you* and *they*.
- In `results/comparator_resampling.csv`, 8 of 2,000 *each other* resamples and 2 of 2,000 *one another* resamples have negative contrasts. Both central 95% ranges remain positive. These are descriptive counts under the recorded scheme, not population probabilities.

## Reproducibility boundary

The current reproduction begins with `matrix_clean.csv`, the archived 155-feature extract. It does not reconstruct the earlier 232→155 reduction. The present singleton filter removes `each_other`, yielding 154 features. Source, archive, and current retained block counts differ as follows:

| Block | Published 232-feature source | Archived extract | Retained analysis |
|---|---:|---:|---:|
| Morphology | 139 | 66 | 65 |
| Phonology | 3 | 3 | 3 |
| Semantics | 36 | 36 | 36 |
| Syntax | 54 | 50 | 50 |

The published counts come from Reynolds (2021), §2.2. The latter two columns are directly inspectable in the supplied matrix, block map, and feature dictionary. Reproducibility of the current finite-instrument calculations does not establish the upstream reduction history.

## Figure 3 and its plotted values

Run `python3 code/revision_review_figures.py` from the repository or supplementary-package root. It reads the saved specification CSV and writes the current Figure 3 to `data/revision/review_corrections/`. `plotted_values.csv` records its 40 values, and `metadata.json` records settings and script/input/output hashes. The legend identifies IDF-weighted Jaccard; the caption specifies full pools, item holdout, and fixed-total placement. The original plots and all numerical outputs are retained unchanged. No analysis or randomization is rerun by this reporting script.
