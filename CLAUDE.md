# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Status

**Title:** Measuring Stable Diagnostic Ambiguity: A Quantitative Workflow for Small-n Grammatical Boundary Phenomena
**Author:** Brett Reynolds
**Status:** Submitted to the Journal of Quantitative Linguistics on 2026-06-09 (submission ID 269804392) after Canadian Journal of Linguistics/Revue canadienne de linguistique no-reviewer rejection on 2026-06-09
**Preprint:** [LingBuzz 009294](https://ling.auf.net/lingbuzz/009294) (Sept 2025)
**Presented:** Paris, April 10, 2026

## Build

```bash
# Paper
pdflatex main-jql.tex && biber main-jql && pdflatex main-jql.tex && pdflatex main-jql.tex

# Analysis (R + Stan)
# See code/ directory - run scripts in numbered order
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

## Core Argument

English reciprocals (*each other*, *one another*) are the proof-of-concept case for a quantitative workflow measuring stable diagnostic ambiguity in small-n grammatical boundary phenomena. Rather than forcing a binary decision, the paper measures **stability of ambiguity** using:

1. 155-feature binary matrix (morphology, syntax, semantics, phonology)
2. Distance measurement with permutation testing
3. Specification curve analysis (multiverse)
4. Comparison group robustness
5. Mixture calibration (Bayesian, Stan)
6. Posterior predictive checks

Result: ~50/50 pronoun-determinative mixture, stable across analytic choices. Supports a property-cluster interpretation of grammatical categories.

## Analysis Pipeline

1. `data_prep_and_exploration.R` - load matrix
2. `01_fake_data_recovery.R` - calibration
3. `02_ppc_cmdstan.R` - posterior predictive checks
4. `03_matched_set_robustness.R` - robustness
5. `04_weight_calibration.R` - mixture weights
6. Stan models: `reciprocals_hpc.stan`, `m0_nomisclass.stan`, etc.

## Dependencies

- R 4.5.x
- tidyverse, Matrix, proxy, vegan, glmnet, cmdstanr, posterior
- CmdStan 2.36.0


## Multi-Agent Dispatch (MANDATORY)

Before dispatching multiple agents, ALWAYS ask Brett:
1. **Which model(s)?** Claude, Codex, Gemini, Copilot
2. **Redundant outputs?** Multiple models on same task for different perspectives?

See portfolio-level `CLAUDE.md` for CLI command patterns and full workflow.
