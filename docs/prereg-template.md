# Preregistration template

Use this for every calibration, validation or confirmation study before anything runs. Copy it to
`docs/notes/<study>/prereg.md` and fill each section; delete none. A section that does not apply says
why. Established 2026-10-12 (user); the practices it collects are recorded in the development record.

## 0. Status line

`DRAFT, NOT FROZEN, NOT YET EXECUTED` until the user freezes it; then `FROZEN <date> by the user`,
committed **before** the run. Any later change is labelled "changed after freezing".

## 1. Question and label

What the study can establish, and what it cannot (the user's scope sentence). Whether it is an
extension of a frozen study, and if so, labelled as such. Anything chosen after seeing data (a cut,
a width, a design) is labelled **chosen after measuring** wherever it appears.

## 2. Inputs and unit of replication

Data or seeds; the unit of replication (fits vs distinct datasets: report both counts); any prior
information (pilot, pre-freeze results), recorded as prior information, never as evidence.

## 3. Procedure and acceptance rule

Every constant fixed: tolerances, bounds and their conventions (exact, one- or two-sided, family
correction), replication counts, seeds and RNG kinds, B. Outcome classification per replication,
denominators, and how unfavourable outcomes count. Nothing excluded silently.

## 3a. Comparator, acceptance and diagnostics (numerical components; user, 2026-10-25)

Required when the study validates a numerical component (a probability routine, a derivative, a
bound, an estimator's numerics). Otherwise, say why it does not apply.

1. **The independent comparator.** Name it, and say on which axes it differs from what is tested:
   **representation** (a different mathematical form), **algorithm** (a different quadrature,
   series or library) and **code** (no helpers shared with the thing tested). List every piece it
   does share (a library's `pnorm`, a constant). Agreement between two variants that share a
   representation and an algorithm is a diagnostic, not certification.
2. **Acceptance versus diagnostics.** Name the one criterion that decides. Everything else —
   convergence checks, self-checks, secondary comparators, envelope regions — is reported as a
   diagnostic and never fails the study. Before freezing, ask of each criterion whether it can
   fail for a reason the measurement never enters (a convergence check stricter than its rule can
   meet, a floor below the reference's own accuracy). If it can, it is a floor (§4), not a
   criterion.
3. **Tolerances, failure criteria and domain coverage**, fixed before the decisive run. The domain
   is enumerated (the cells, points or regions), with the validity conditions of any bound used
   (e.g. first order: δP ≪ P, with its number). A point outside the domain is a domain failure,
   not an excluded point.
4. **Provenance**: precision, package and version, method agreement, and machine, recorded with
   the result. A reference that becomes a committed value goes through `fixtures/generate/*.R`
   and its manifest.
5. **What this study does not certify.** Name the neighbouring layers (the reference, the
   estimator, the inference, rotation) and the conditions it leaves open, so that a pass here is
   not read as a pass there.

*Why:* the 2026-10-25 Φ₂ tail reference's two methods shared their integrand and quadrature, and
its half-step convergence check, written as a pass/fail criterion, failed 569 corners by design
(E-063, section 8 (e)).

## 4. Design-based failure floors (REQUIRED before any model fitting or selection; user, 2026-10-12)

Computed from the design alone, with **no observed statistic of interest** (no cut statistics, no
estimates). For a procedure that must pass on several cells or sides:

1. **Per cell:** the probability that one cell fails under representative event counts,
   **including zero observed events** (e.g. the predictive tail under the study's own prior at 0/R),
   and under true rates fixed in advance as fractions of the bound (0, ¼, ½, 1).
2. **Per side:** aggregate over the complete set of cells that must pass.
3. **Feasibility of the decision rule:** check whether the side-level floors permit the intended
   comparison (e.g. two risk curves crossing, a pass probability above a target) at all. If they
   cannot, revise the design before freezing.
4. Include every cut-independent component (operational failures, non-convergence) as its own floor,
   reported separately and never read as evidence about the quantity studied.

Record the numbers and the script that produced them (`design_*.R` and its log) in the study's
folder. Only then fit tail models or run the selection computation.

*Why:* the 2026-10-12 predictive-risk selection (E-093) was frozen
without its harmful-side floor; at 0 misses in 100 the frozen prior gave 0.947 across 24 cells, so
the equal-risk rule could not select whatever the data. A general design check, not a correction to
that result.

## 5. Decision rule and fallback

The outcome for each case, including "inconclusive", and its consequence (release scope, refusal).
Sensitivity analyses named in advance, with the threshold for "materially sensitive".

## 6. Script and smoke run

The analysis script, smoke-run on synthetic data only (the real statistic never read), its SHA-256
recorded here at freezing. The run uses the committed file, verified against the hash.

## 7. Reporting

Report the outcome whatever it is. A preregistered failure is preserved, never re-run with other
seeds, constants or models. Explanations computed after the result are labelled descriptive.
