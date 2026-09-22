# Decisions Log
<!-- SUMMARY: Accepted scar-tissue pass applied across submission text; current manuscript, response, supplement and administrative files verified · status: revising · updated: 2026-09-22 -->

Append-only record of project decisions. Agents: add an entry whenever a non-trivial decision is made during a session (structural changes, venue choices, theoretical commitments, scope changes, reviewer feedback acted on). Keep entries short.

Format: `## YYYY-MM-DD` then bullet points with **bold topic** and brief rationale.

---

## 2026-04-08

- **Ball-in-valley vs spinning-top metaphor needs nuance at book level.** Cross-linguistically convergent categories like NOUN may genuinely be valleys (deep attractors in the space of possible grammars), not just spinning tops. The spinning top fits boundary items and less universal categories. The book could use a richer landscape: valleys for convergent categories, spinning tops for actively maintained ones, boundary zones where valley walls are shallow. Connects to Powell's convergence-under-constraint framework. Not for the slides -- for the book's stability chapter.

- **Mechanism list reframed from formal systems to processes.** Slide 5's list changed from "morphological realization rules, agreement/case systems, entrenched distributional patterns, grammaticalization pathways, community norms" to the HPC book's process-oriented list: acquisition, entrenchment, interactive alignment, iterated transmission, functional pressure. Reason: the formal-systems framing conflates mechanisms with their products. The book (Ch. 4) explicitly treats these as functional roles ("stabilizers"), not formal properties. Q&A point: formal systems like morphology are *products* of these mechanisms, not mechanisms themselves.

- **Projectibility gap fixed on slides 15 and 17.** Review board (projectibility reviewer) flagged projectibility as "decorative" -- introduced on slide 3, illustrated on slide 4, then absent for 13 slides, with the conclusion naming "maintenance" as the HPC answer. Fixed: slide 15 thesis now cashes out the projective consequence ("roughly half a pronoun's behaviour, half a determinative's"), slide 17 adds step 6 ("Cash out the projective consequences") and reframes the concluding thesis around projectibility as the payoff. Reason: Boyd's slogan is "profile, stabilised by mechanisms, projectible for a purpose" -- the third clause was missing.

- **"Back to the spinning top" slide removed.** Mapped specific mechanisms (entrenchment/transmission, functional pressure/acquisition) to specific feature dimensions (morphology, semantics). The paper measures stability of boundary position, not which mechanisms maintain which properties. The mechanism-to-dimension mapping was speculation the data doesn't support. Loop-closing now happens across slides 16-17 without the intermediate speculative slide.

- **Added teaser slide 2 "The puzzle".** Four theory slides before any mention of reciprocals was too much setup. New slide gives the audience something concrete to hold onto while the HPC framework builds. Addresses Harris's delivery note from the February reviews.

- **New data slide (slide 7).** Introduces the 155-feature × 138-item matrix with a mini-table of 6 features × 4 items (each other, one another, they, somebody). Makes the mixed pattern visible before the dimensional breakdown.

## 2026-04-13

- **Separated the inventory override layer from the feature matrix.** New pronoun/determinative follow-on work will use a derived annotation table rather than rewriting `matrix_clean.csv`, so the original reciprocals analysis remains reproducible while newer category and personhood claims can be layered on cleanly.

- **Adopted conservative first-pass anchor buckets.** Clear anchors are now derived mechanically from the override table: `retain` rows become class-specific anchors, `review` rows stay out of the first pass, and known boundary or non-core cases are excluded. This keeps the first pronoun/determinative comparison high-precision while preserving a transparent queue of items for later expansion.

- **Standardized `what` as determinative throughout.** Exclamative `what a N` is now treated as the count-singular realization of the same determinative seen in plural and mass exclamatives, rather than as an adjective-like residue or a separate `what a` item. The legacy `what_pron` row remains only as an excluded compatibility row to avoid double counting the seed inventory.

- **First pronoun/determinative baseline is anchor-only.** The initial comparison now trains only on clear pronoun and determinative anchors, then scores review items out of sample with a ridge classifier plus Jaccard MDS. This keeps borderline rows from contaminating the first separation while still producing an ordered queue for theoretical cleanup.

- **Temporal deictics removed from the pronoun inventory.** `today`, `tomorrow`, `tonight`, and `yesterday` are now treated as nouns, not pronouns, and are excluded from the pronoun/determinative baseline. This removes a spurious temporal subclass from the review queue and keeps the comparison focused on the intended lexical categories.

- **Promoted personal and residual wh-determinatives out of review.** `we_det`, `us_det`, `you_det`, `whatever_det`, `whichever_rel_det`, and `whichever_int_det` are now treated as determinatives throughout and retained as clear anchors. This aligns the operational inventory with the intended theoretical analysis instead of letting pronoun-overlap heuristics keep them artificially in review.

- **Removed `there` from the pronoun/determinative comparison.** The current seed row conflates existential/tag pronoun uses with the broader preposition use, so `there` is now excluded from the inventory until a context-sensitive split is added. This prevents the baseline from treating a mixed row as a meaningful boundary item.

## 2026-04-15

- **Follow-on inventory stabilized for handoff.** The override layer now treats temporal deictics as nouns, promotes the remaining personal and wh-determinatives to clear determinative anchors, and excludes conflated `there` from the comparison. This reduces the operational residue to a single review item (`whatever_pron`) before corpus retagging moves to a separate project folder.

- **This repo is frozen as the theory/inventory baseline, not the corpus-engineering workspace.** The pronoun/determinative baseline, derived inventories, and diagnostic plots were committed and pushed as `8fe9ec2` (`Add pronoun-determinative inventory baseline`). Future EWT/GUM retagging work should proceed in the separate retagging project, using this repo as the inventory/decision source.

## 2026-04-23

- **Journal-version workflow additions are prior predictive checking, a canonical model ladder, and SBC.** The manuscript now treats these as the next serious workflow upgrades: show what the anchor and misclassification priors imply before fitting, move through one explicit sequence from simple classifier to fuller generative model, and validate the final Stan implementation with simulation-based calibration rather than relying on fit diagnostics alone.

- **PSIS-LOO is secondary in this project.** If predictive comparison enters the workflow, it should follow posterior predictive checks and report Pareto-`k` diagnostics. The main question here is boundary location, not model ranking, so cross-validation should not drive the rhetoric of the paper.

## 2026-04-24

- **House style is binding even against reviewer preference.** The CJL draft keeps contractions, removes banned connective/modals (`however`, `yet`, `must`, `cannot`, `therefore`, etc.), and avoids letting external copyedit suggestions override the explicit house rules. Reason: this project's style guide is an intentional prose system, not a loose preference list.

- **Table 1 uses ragged-right text in its first two columns.** The first two `p{}` columns were switched to `>{\raggedright\arraybackslash}p{...}` with `array` loaded in the preamble. Reason: narrow fully justified cells were producing visibly bad spacing in the diagnostic and illustration columns.

- **Statistics prose is shifting from tool-name-first to question-first exposition.** Reader-facing prose now prefers `row- and column-preserving reference distribution` / `randomizations` and explains what each lens is supposed to show before naming the statistic; `quasiswap` is retained mainly where the algorithm itself must be identified. Reason: CJL readers need the motivation for each method before the terminology.

- **Methods sections now define the distance contrast before benchmarking it.** The manuscript introduces `\Delta` before the constrained reference distribution, explains the row- and column-preserving benchmark before naming `quasiswap`, separates numerical specification sensitivity from comparator-set rotation, and retitles the blend section around predictive calibration. Reason: each statistical lens should answer a clear reader-facing question before its label or implementation detail appears.

- **HPC stays as a compact projectibility frame, not the paper's headline.** The introduction now unpacks the homeostatic-property-cluster claim as a claim about when partial diagnostic profiles support projection, and the discussion qualifies mechanism language as an interpretation consistent with the matrix rather than something directly identified by it. Reason: HPC should explain why stable mixed placement matters without turning the CJL article into a philosophy-of-science paper.

- **CJL is the submission venue for the journal version.** The paper was submitted to the Canadian Journal of Linguistics/Revue canadienne de linguistique on 2026-04-24. Reason: the paper now fits CJL as a methodologically explicit English-grammar article rather than as a broader philosophy-of-linguistics piece.

- **Upload source must be a scrubbed anonymous source package, not raw `main.tex`.** The working `main.tex` keeps identified and anonymous branches, so the upload package is `submission/latex_source_anonymous/` plus `submission/cjl_anonymous_latex_source.zip`. Reason: raw source would expose author metadata and non-anonymous branch content.

- **Submission materials are split by review function.** Identifying material lives in `CJL_TITLE_PAGE.*`, reviewer suggestions in `CJL_REVIEWERS.*`, accessibility descriptions in `CJL_ACCESSIBILITY.md`, and source files in the anonymous source package. Reason: CJL/ScholarOne wants anonymity in the manuscript but separate author/contact/declaration materials.

- **Excluded reviewers are framework conflicts rather than intellectual opponents.** Pullum and Huddleston are listed as opposed/potential conflicts because of coauthorship and their central role in the CGEL framework, while Aarts, Denison, Keizer, Payne, and Wallis are recommended as content/method fits. Reason: this gives the editor usable expertise without violating ordinary conflict expectations.

## 2026-05-27

- **Classify dirty state as submission/admin tracking residue.** `submission/CJL_SUBMISSION.md` now records manuscript ID `CJLRCL-2026-0026` and the ScholarOne author-centre URL; `writing-style.md` is a local style-rule symlink now ignored by `.gitignore`. No manuscript source, anonymous source package, PDF, data, or analysis file changed.

## 2026-06-09

- **Treat the CJL/RCL decision as no-reviewer / venue-fit, not substantive rejection.** CJL/RCL rejected `CJLRCL-2026-0026` because the editors could not secure reviewers despite several attempts and suggested a more specialized journal. Reason: there are no reviewer reports to answer, so the next step is venue retargeting rather than manuscript revision driven by feedback.

- **Retarget to the Journal of Quantitative Linguistics.** The JQL pitch should foreground the paper as a quantitative workflow for stable diagnostic ambiguity in small-*n* grammatical boundary phenomena, with English reciprocals as the proof-of-concept case. Reason: this matches JQL's methods-plus-theoretical-understanding scope better than another broad English-grammar venue.

- **JQL refit makes the measurement workflow the headline.** The working manuscript now uses the title `Measuring Stable Diagnostic Ambiguity: A Quantitative Workflow for Small-n Grammatical Boundary Phenomena`, opens with the small-*n* measurement problem, adds an explicit seven-step workflow, and moves HPC to theoretical payoff rather than premise. The CJL branch is preserved under the existing `\CJLSubmission` conditional, with `main-jql.tex` as the new build wrapper. Reason: Brett accepted the measurement-forward title in Roughdraft and asked for a strong JQL fit without sunk-cost thinking.

- **JQL submission is now the active state.** The paper was submitted to the *Journal of Quantitative Linguistics* on 2026-06-09 with submission ID `269804392`. Reason: the Taylor & Francis confirmation establishes the new tracking state; future revisions should start from the JQL package rather than the old CJL package.

2026-08-07 — Source routed in, R&R material only (under review at JQL as 269804392): Goldsmith-Pinkham, Hull & Kolesár, "Leniency Designs: An Operator's Manual," *JEP* 40(3), 2026, 213–240. Hook at `notes/source-hooks/goldsmith-pinkham-2026-leniency-designs.md`. The leniency-design literature opens on this paper's exact phenomenon — expert decision-makers who agree on clear cases and diverge "systematically" on close calls — but **inverts** its treatment: economics uses that divergence as an instrument for identifying the effect of the decision on later outcomes, never asking what the experts are disagreeing about, whereas this paper makes the divergence the measurand. Two possible uses in an R&R: as unusually strong outside corroboration that stable expert disagreement is real, measurable, and systematic enough to have carried a decade of applied work in a field with no stake in the linguistic question; and as a contrast that states the contribution more sharply (rater variation as object of study rather than as nuisance parameter). Caution recorded in the hook: the designs are not commensurable — leniency designs need as-good-as-random case assignment and an exclusion restriction, and there is no instrument in a grammaticality-judgement study — so this is a conceptual parallel and outside evidence for the premise, not a method to import.

## 2026-09-22

- **JQL decision: reject with invitation to resubmit after major revisions (269804392).** Editor Emmerich Kelih; "some merit," would reconsider a major-revision resubmission, re-reviewed, due 2027-03-21 (extension on request). Three reviewers, 2 positive / 1 major-but-constructive: **R1** recommends acceptance as it stands (only reservation: subjectivity of the feature matrix, mitigation deferred to future work); **R3** positive with minor asks (be more cautious on the epistemic-vs-ontological framing; add analysis of coding-bias impact; reconcile the Delta-magnitude instability in Section 5 with the "stable" claim; integrate acceptability/frequency literature on reciprocals); **R2** a long, constructive methodological critique who expects it can become publishable after revision. **Make-or-break (R2, the stated credibility test):** run the pipeline on clear anchors (EVERYONE, SHE) and show they are *not* classified as boundary cases; if clear anchors also come out "boundary," the method is undermined. Other load-bearing R2 asks: the Section 3 null/reference distribution is the wrong null (keep anchors fixed, randomize only the target; report distance-to-each-set, not just Delta); the Section 6 mixture model is likely misspecified (justify or drop); normalize Delta across metrics via reference-distribution quantiles; define set-distance aggregation and justify the 0/1 coding given Jaccard's asymmetry; add an illustrative data table (a few words x a few features) and clarify the 155 features / 138 words and the lemma-vs-wordform inconsistency; the Delta metric ignores absolute distances (2,2 and 50,50 both give Delta=0). Consistency/proofing: a flat contradiction about the second category (Fig 1 caption / Table 1 vs p.5); terminology "determinative-like" to be replaced with "anybody-ward"/"she-ward". Owner: Brett, decision pending (undertake the major revision or not). Assisting system: Claude Code (Opus 4.8), triage only. Record: `reviews/jql-decision-2026-09-22.md` (gitignored, public repo).

- **Undertake the JQL revision.** Brett's instruction on 2026-09-22 resolves the preceding decision to proceed. Codex prepared `reviews/jql-revision-plan-2026-09-22.md` and opened it in Roughdraft. Comparator definitions, replacement validation, and removal of the blend from the main argument remain proposals pending that review; no new analytical result or revised manuscript conclusion is claimed.

- **Revision proposal approved.** Brett confirmed completion of the Roughdraft review and approved all proposals on 2026-09-22. The source check fixes the primary pools at 41 core personal-pronoun forms/readings and all 16 central compound determinatives in CGEL Ch. 5 §9.6, with the 12 nonlocatives as a sensitivity check. Original inputs remain unchanged; the documented singleton filter yields 154 retained features.

- **Corrected results narrow the claim.** The new full-feature comparison uses 154 retained features, 41 core personal-pronoun forms/readings, and 16 central compound determinatives. Both targets fall between observed control contrast ranges with all features; deleting morphology moves them into the pronoun range under several measures. The working title is now *Calibrating Diagnostic Conflict: English Reciprocals in Grammatical Feature Space*. A stable category-boundary diagnosis, the small-target-count obstacle to classification, and the mixture-weight argument have been withdrawn.

- **Validation is conditional on the coded instrument.** Direct item/family holdout replaces the generative-geometry argument; fixed-anchor target randomization replaces the old benchmark. A single GPT-5.6 Terra audit selected by Brett reproduced the principal results and found no material numerical or method defect. Generated tables, final script hashes, and a self-contained supplementary archive provide traceability. These checks do not independently validate the inherited coding or establish population accuracy.

## 2026-09-22 — Author feedback on revision

Preserve the small-n limitation for characterizing a proposed subset from two lexical types while distinguishing it from comparing those types against larger reference groups. Display distances, contrasts, and resampling ranges to two decimal places, retain full-precision saved results, and compute contrasts before rounding. Both requested external studies are now local and checked; Hurst and Nordlinger is the 2007 author preprint of the 2011 chapter. These bounded edits leave the audited numerical analysis unchanged.

## 2026-09-22 — Full source intake approved

Brett approved the full Haas / Hurst–Nordlinger intake without changes. Applied its two discussion paragraphs and Reviewer 3 response clarification. Acknowledge Haas's actual corpus frequency comparison while preserving the distinction between evidence about use and validation of our distance contrasts. The manuscript now specifies Hurst–Nordlinger's nine-speaker elicitation and task limits. The approved interpretation and all audited numerical results are unchanged.

## 2026-09-22 — Approved post-review weighting check

- **Separate the morphology result from the weighting rule.** Equal block totals preserve a small compound preference after morphology is removed (about 0.02), whereas equal feature weights reverse it. The target–control comparison also depends on holdout: the reweighted no-morphology targets fall within the pronoun range only under family holdout. The conclusion is conditional on feature selection, weighting, and comparator design.
- **Keep the added calculation bounded and exploratory.** The approved dated protocol specifies one alternative weighting over the full pools, two holdouts, and two feature sets: 236 new rows, no new random draws. All original corrected outputs and scripts are preserved byte-for-byte. The same approved GPT-5.6 Terra auditor independently verified the addendum with no defect found.
- **Distinguish error sources in reporting.** Monte Carlo standard errors and the zero-hit binomial bound concern the artificial placement benchmark. Family resampling and coding/weighting perturbations describe sensitivity; they do not estimate independent coding reliability or population uncertainty. The supplement exposes the full original grid and every eligible-control margin.

## 2026-09-22 — Independent-review clarification pass

- **Treat subject-position coding as contestable legacy scope.** Restricted embedded-subject uses of `each other` are documented in the cited sources. Preserve the archived 0, identify its constructional ambiguity, and point to the existing target-cell perturbation; do not invent an original main-clause-only definition or silently recode.
- **Explain identical targets and the dummy-it endpoint.** The only retained target differences are four morphological component columns. The 89 non-morphological cells are identical. In the equal-block no-morphology comparison, target contrasts do not change across holdouts; dummy `it` moves the upper pronoun-range endpoint past them. This qualifies the range interpretation without altering any result.
- **Finish provenance and display clarification separately from analysis.** The current reproducible chain begins at the archived 155-feature extract, not a reconstructed 232→155 reduction. A separate reporting script renders Figure 3 with the explicit IDF label and settings from unchanged saved values. Original numerical scripts, outputs, and the private review packet remain frozen.
- **Include the approved small transparency additions.** Report 8/2,000 and 2/2,000 negative family resamples; use seven analysis families rather than an ambiguous paradigm count; identify `where`/`locative` as a concrete duplicate. Brett approved the complete wording plan without comments or edits. No new study or further agent audit was needed.

## 2026-09-22 — Claude-review manuscript corrections approved

- **Prioritize the manuscript and apply the approved wording.** Brett approved all 27 proposed replacements without edits or comments. The review copy and pre-review snapshot match. The response letter is a subsequent alignment task.
- **Explain the source of the weighting result.** The positive equal-block-total contrast depends on two phonological columns; balancing semantics and syntax produces Δ = −0.11. Coordinated morphological deletions show that whole-word codes, rather than component identities, support the primary compound preference.
- **Calibrate range and benchmark statements with controls.** Dummy it sets every original pronoun endpoint and is the sole recovery failure across the declared grid and recorded weighting settings. Fifty of 57 controls have low random-placement tails against their non-designated group. The manuscript now foregrounds observed two-distance comparisons.
- **Retain documented uncertainty.** Brett confirmed that the original `does_not_inflect` definition is uncertain. Source facts about genitives and missing person/number/gender inflection aren't substituted for a recovered coding rule. The archive remains unchanged. The 13-column deletion is described as mixed morphological properties, including word formation and monomorphemicity.
- **Keep the new computations retrospective and separate.** The post hoc script and reporting record preserve earlier numerical scripts/results and generate no new analytical random draws. The author approved the findings and wording after local verification of the independent Claude report. Full provenance is in `data/revision/POSTHOC.md` and `reviews/claude-manuscript-implementation-2026-09-22.json`.
- **Disclose actual AI roles.** Local session records identify GPT-6 Astra for editing, Brett identifies Pro as Sol 6, Elicit is dated without an inferred model, Claude's report identifies Opus 5.5, and the earlier separate Terra numerical audit remains attributed accurately. No commit, push, or resubmission.

## 2026-09-22 — Response-letter alignment completed

- **Align the letter to the approved manuscript.** Brett's instruction to proceed authorized this follow-through. Corrected both letter bodies at 15 locations, including the phonological source of the block-weight result, dummy-it range endpoints, observed-control benchmark calibration, source motivation, and uncertain coding definitions.
- **Preserve the review history and keep the stage bounded.** The original four comments and YAML replies are unchanged; the clean and annotated bodies agree. No manuscript, analysis, supplement, or historical full-review packet changed. The clean response is ready for Roughdraft review; no new approval or submission is recorded. Verification: `reviews/response-alignment-verification-2026-09-22.json`.

## 2026-09-22 — Response approved; scar-tissue pass requested

- **Record the approval and review the whole submission text.** Brett approved the aligned response without comments or edits; the saved file matches its review hash. He then requested an editorial scar-tissue pass on all text to be submitted and explicitly included updated administrative text.
- **Repair prose without reconstructing the argument.** The audit proposes grouping related coding notes, leading with the observed control comparison, and removing repeated cautions and references to superseded wording. The abstract, conclusion, substantive limitations, numerical results, and source set remain intact. Protocol chronology is retained as methodological provenance.
- **Preserve historical submission material.** The new cover letter, title page, form text, and image descriptions are separate revision drafts. June files and historical review packages stay frozen. Proposed edits are staged for one Roughdraft checkpoint before application; no final pass completion, build, or submission is claimed yet.

## 2026-09-22 — Accepted scar-tissue pass completed

- **Apply the accepted whole-submission edits.** Brett accepted the complete review and new administrative text without changes. Applied all 50 replacements; retained the argument, abstract, conclusion, numerical findings, source set, coding uncertainty, and reviewer-response history.
- **Finish the current submission copies.** The new cover letter, title page, form sheet, and image descriptions are in `submission/revision-2026-09-22/`. The 28-page PDF builds cleanly and has been visually checked. Its inclusive word count is 6,393; the comparable source count falls from 6,045 to 5,803. Layout-only changes leave accepted wording intact.
- **Preserve the analytical and historical records.** The refreshed 64-file supplement changes narrative Markdown and checksums only; 58 entries are byte-identical. Both historical full-review ZIPs and June administrative files remain unchanged. The completed pass is recorded against the final source. Evidence: `reviews/scar-tissue-2026-09-22/completion.md`. No commit, push, or resubmission.
