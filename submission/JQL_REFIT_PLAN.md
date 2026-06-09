# JQL Refit Plan
Date: 2026-06-09
## Starting Point
CJL/RCL rejected the manuscript without external review because the editors could not secure reviewers. Treat this as a venue-fit/no-reviewer outcome, not substantive feedback.

The existing paper is already close to JQL in methods, but its CJL-facing version still reads partly as an English grammar article with quantitative support. For Journal of Quantitative Linguistics, the paper should read first as a contribution to linguistic measurement: a portable workflow for small-n grammatical boundary phenomena, with English reciprocals as the proof-of-concept.
## Venue Constraint
The official JQL page says the journal wants work that systematically applies or develops mathematical and/or statistical concepts and methods for theoretical understanding of language phenomena. It also warns that ad hoc computational implementation is not enough, and that simple application of existing theory to language data without substantive theoretical advancement is unlikely to fit.

Implication: the refit should make the paper's theoretical-quantitative contribution explicit. The argument cannot be only "I used several quantitative tools on English reciprocals." It needs to be "stable diagnostic ambiguity is a measurable object, and this workflow shows how to quantify it under small-n constraints."
## Load-Bearing Assumptions
1. **JQL will consider grammar/category measurement within scope.**
  
  - Why plausible: the scope includes morphology, syntax, semantics, theoretical linguistics, methodological problems of linguistic measurement, and mathematical/statistical models.
    
  - Failure condition: if editors interpret JQL as primarily quantitative laws or text/corpus regularities, this paper may still look peripheral.
    
  - Design response: foreground measurement theory, not just reciprocals.
    
2. **The existing analysis can be reframed without new empirical work.**
  
  - Why plausible: the current manuscript already has distance metrics, constrained randomization, specification curves, comparison-set rotations, predictive calibration, and anchor checks.
    
  - Failure condition: if JQL expects a more formal mathematical contribution or broader empirical validation.
    
  - Design response: define the workflow more formally and identify scope conditions and disconfirmation criteria.
    
3. **The hand-coded small-n dataset is acceptable if presented as a measurement problem.**
  
  - Why plausible: the whole contribution is about what can be responsibly inferred when the target inventory is intrinsically tiny.
    
  - Failure condition: if reviewers treat small target n and single-coder features as fatal.
    
  - Design response: make single-coder, theory-dependent coding an explicit limitation, and emphasize calibration against anchors and robustness surfaces.
    
## Recommended Refit
### 1. Title
Use a measurement-forward title. I recommend:

> {==Measuring Stable Diagnostic Ambiguity: A Quantitative Workflow for Small-n Grammatical Boundary Phenomena==}{>>this one's good<<}{id="c1" by="user" at="2026-06-09T16:22:28.437Z"}

Subtitle option if we want the case visible:

> English Reciprocals between Pronouns and Determinatives

This is stronger for JQL than the current title because it names the quantitative object and workflow before the case.
### 2. Abstract
Replace the CJL abstract with a JQL-facing abstract that opens on the general measurement problem. Keep the concrete result, but make the portable workflow the first claim.

Core wording to adapt from `submission/JQL_PITCH.md`:

> This paper develops a quantitative workflow for small-n grammatical boundary cases. Rather than treating inconsistent diagnostics as failed classification, it asks whether diagnostic ambiguity is stable across independently varied measurement choices.
### 3. Introduction
Revise the introduction in three moves:

1. Start with the methodological problem: grammatical categories often have tiny boundary inventories, so ordinary large-sample classification and single-diagnostic argument both fail.
  
2. Define "stable diagnostic ambiguity" as a measurement target.
  
3. Introduce English reciprocals only after the method question is established.
  

Keep HPC/projectibility, but move it from headline frame to theoretical payoff. JQL readers need the math/statistics problem before the philosophy-of-categorization frame.
### 4. Methods Structure
The current section sequence is defensible, but the section titles should become more method-explicit:

- "The Challenge of Measuring Grammatical Similarity" -> "A Feature-Space Measurement Model for Grammatical Categories"
  
- "Benchmarking the Distance Contrast..." can stay, but tighten the transition to explain why constrained randomization is a model check.
  
- "The Garden of Forking Paths" -> "Specification-Curve Analysis of Metric and Feature-Block Dependence"
  
- "The Importance of Comparison Groups" -> "Comparator-Set Sensitivity"
  
- "Predictive Calibration of Boundary Position" can stay.
  
- "Checking the Instruments" -> "Anchor Calibration and Predictive Checks"
  

The aim is to make each section read as a component of one quantitative workflow, not as a list of robustness checks.
### 5. Formal Workflow Box
Add a short numbered workflow near the end of the introduction or start of methods:

1. Encode a theory-informed item-by-property matrix.
  
2. Locate targets relative to clear anchors in full-dimensional feature space.
  
3. Define an interpretable distance contrast.
  
4. Compare the contrast to a constrained reference distribution.
  
5. Vary metrics, feature blocks, and comparators.
  
6. Calibrate boundary location with anchor recovery and predictive blending.
  
7. Report stable ambiguity only when uncertainty localizes to the boundary items and clear anchors remain recoverable.
  

This is the main JQL adaptation. It converts the current pipeline into an explicit contribution.
### 6. Results and Discussion
Retitle or restructure the results/discussion so the conclusion is not mainly "reciprocals sit at the interface." For JQL, the stronger conclusion is:

> The workflow distinguishes arbitrary diagnostic conflict from stable boundary placement under a specified measurement model.

Keep the reciprocal result as the demonstration.

Discussion should add one paragraph on how this relates to quantitative linguistics:

- linguistic measurement under sparse target inventories;
  
- distance metrics for categorical feature spaces;
  
- sensitivity analysis as part of theory construction;
  
- stable ambiguity as a measurable property of category systems.
  
### 7. Submission Materials
Create JQL-specific files without overwriting the CJL package:

- `main-jql.tex`: identified JQL build wrapper or direct working version.
  
- `submission/JQL_TITLE_PAGE.md`: title page, declarations, preprint note, AI-use disclosure.
  
- `submission/JQL_COVER_LETTER.md`: adapted from `submission/JQL_PITCH.md`.
  
- `submission/JQL_SUPPLEMENT.md`: neutral supplement checklist, not CJL-branded.
  

Do not modify or delete the CJL provenance package.
## Edits I Would Make First
1. Add `main-jql.tex` as a new wrapper rather than altering `main-cjl.tex`.
  
2. Add a `\JQLSubmission` conditional only if needed; otherwise use the identified `main.tex` branch and keep JQL-specific wording in the main source.
  
3. Rewrite title, abstract, first two introduction paragraphs, and the paragraph that previews the analyses.
  
4. Add the formal workflow box.
  
5. Retitle methods sections and smooth transitions.
  
6. Revise the final conclusion to foreground the measurement contribution.
  
7. Draft cover letter and title page after the manuscript language is stable.
  
## Risk Notes
- The paper should not oversell independence: all analyses reuse one hand-coded matrix. Keep "internal robustness" language.
  
- Do not bury single-coder feature coding. JQL reviewers may be more method-sensitive than CJL reviewers.
  
- Avoid making HPC sound like the required theoretical buy-in. The measurement contribution should stand without accepting HPC.
  
- Avoid describing the methods as a generic multiverse. The specific object is stability of diagnostic ambiguity under constrained measurement choices.
  
## Sources Checked
- JQL official Taylor & Francis page: [https://www.tandfonline.com/journals/njql20/about-this-journal](https://www.tandfonline.com/journals/njql20/about-this-journal)
  
- JQL Taylor & Francis journal overview: [https://www.tandfonline.com/journals/njql20](https://www.tandfonline.com/journals/njql20)
  
- JQL current/latest article listings from Taylor & Francis, including recent measurement/method articles: [https://www.tandfonline.com/toc/njql20/current](https://www.tandfonline.com/toc/njql20/current)
  
- Existing local pitch: `submission/JQL_PITCH.md`
  
- Existing manuscript: `main.tex`
  
- CJL submission notes: `submission/CJL_SUBMISSION.md`

---
comments:
  c2:
    body: Let's make it a strong fit. No sunk-cost thinking.
    by: user
    at: 2026-06-09T16:23:41.391Z
