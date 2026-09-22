# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Status

**Title:** Calibrating Diagnostic Conflict: English Reciprocals in Grammatical Feature Space
**Author:** Brett Reynolds
**Status:** JQL major-revision manuscript and response updated 2026-09-22 after the approved source intake, Elicit follow-up, Pro corrections, and Claude-review corrections. Brett accepted the complete whole-submission editorial scar-tissue pass and administrative text as is. All edits are applied and verified: 28-page PDF, aligned response versions, refreshed 64-file analysis supplement, and current cover-letter/title-page/form/image-description text in `submission/revision-2026-09-22/`. The final inclusive PDF count is 6,393; the source inventory is 5,803. Completion evidence is in `reviews/scar-tissue-2026-09-22/completion.md`. Numerical code/results, June submission records, and historical full-review ZIPs are unchanged. Ready for resubmission preparation when instructed (due 2027-03-21); nothing resubmitted.
**Preprint:** [LingBuzz 009294](https://ling.auf.net/lingbuzz/009294) (Sept 2025)
**Presented:** Paris, April 10, 2026

## Build

```bash
# Paper
pdflatex main-jql.tex && biber main-jql && pdflatex main-jql.tex && pdflatex main-jql.tex

# Revised analysis
python3 code/revision_analysis.py
python3 code/revision_tables.py
python3 code/revision_followup.py
python3 code/revision_review_figures.py
python3 code/revision_posthoc.py
# Tests: from code/, python3 -m unittest -v test_revision_analysis.py test_revision_followup.py
```

## Structure

- `main.tex` - working paper source with JQL-facing identified branch and preserved CJL conditional branch
- `main-jql.tex` - JQL build wrapper
- `main-cjl.tex` - anonymous CJL build wrapper
- `submission/` - CJL provenance materials plus JQL submission package and tracking
- `refs.bib` - bibliography (unified style)
- `code/` - R scripts + Stan models
- `data/` - CSV matrices and outputs
- `plots/` - generated figures

## Current Argument and Analysis

The revised manuscript compares English reciprocals with explicit personal-pronoun and compound-determinative groups. The unchanged archived matrix yields 154 retained features after excluding its singleton identity column. The primary anchors are 41 personal-pronoun forms/readings and all 16 central CGEL compound determinatives; 12 nonlocatives form a sensitivity pool.

The primary Jaccard instrument recovers all 57 controls. Each reciprocal is nearer each group than every member of the other group is, but farther from it than its own controls, apart from dummy it. Whole-word morphological codes, rather than component identities, support the modest compound preference. Removing morphology reverses it; equal block totals offset that reversal through two heavily weighted phonological columns. Balancing semantics and syntax gives Δ = −0.11. Dummy it sets every original pronoun-range endpoint and is the only control ever misrecovered. Low random-placement tails also occur against the non-designated group for 50/57 controls. The original `does_not_inflect` definition remains uncertain; archived cells and earlier results are preserved. Both targets have identical non-morphological vectors, and the reproducible chain begins at the 155-feature archive. The conclusion is conditional on coding, feature selection, weighting, and comparator design.

- `code/revision_analysis.py`: current analysis and figures.
- `code/revision_tables.py`: manuscript numbers/tables generated from saved results.
- `data/revision/README.md`: reproducibility instructions and output dictionary.
- `data/revision/protocol.json` and `comparators.csv`: recorded design and item manifest.
- `data/revision/results/`: original corrected results and provenance hashes.
- `code/revision_followup.py`, `data/revision/followup_protocol.json`, and `data/revision/followup/`: separately recorded weighting check and expanded summaries of saved results.
- `data/revision/FOLLOWUP.md`: addendum design, findings, and output dictionary.
- `data/revision/CODING_NOTES.md`: source/coding scope, target identity, range endpoints, and provenance boundary.
- `code/revision_review_figures.py` and `data/revision/review_corrections/`: current Figure 3 with an explicit IDF label, generated from 40 unchanged saved values. The original figure remains in `results/figures/`.
- `code/revision_posthoc.py`, `data/revision/posthoc_protocol.json`, `data/revision/POSTHOC.md`, and `data/revision/posthoc/`: retrospective reporting of the phonology, morphology, endpoint, and control-benchmark diagnostics; no new random draws.

The original R/Stan scripts and prior outputs remain as historical records. They do not generate the revised manuscript's results. The separate 2026-04 pronoun/determinative inventory work also remains unchanged.

## Dependencies

The revised pipeline uses Python, NumPy, pandas, and Matplotlib; tests additionally use SciPy. Exact versions are recorded in `data/revision/environment.json`. Legacy R/Stan dependencies remain documented in the repository README and original submission package.

## Multi-Agent Dispatch (MANDATORY)

Before dispatching multiple agents, ALWAYS ask Brett:
1. **Which model(s)?** Claude, Codex, Gemini, Copilot
2. **Redundant outputs?** Multiple models on same task for different perspectives?

See portfolio-level `CLAUDE.md` for CLI command patterns and full workflow.
