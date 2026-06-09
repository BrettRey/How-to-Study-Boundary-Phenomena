# Decisions Log

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
