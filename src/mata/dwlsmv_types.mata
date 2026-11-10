*! dwlsmv -- the locked type boundary (spec Section 5.1, Gate 6)
*!
*! Executed under StataMP 18 on 2026-09-18 and it loads clean: tests/unit/test_statmap.do
*! builds these structs and exercises them. Nothing here computes -- they are
*! type declarations -- but the declarations themselves are now known to parse.
*!
*! This file defines the objects that cross module boundaries, and nothing
*! else. It contains no function bodies deliberately: a header full of stubs
*! that "work" is how a module gets mistaken for implemented. Function
*! signatures are specified in spec Section 5.1 and are binding; the bodies
*! land in one file per function, as the spec requires.
*!
*! Changing a field here changes every module that touches it. If you are
*! tempted to add a field, check first whether it is derivable from what is
*! already present -- a derived quantity stored in a struct is a second copy
*! that can disagree with the first.
*!
*! This file does NOT issue `mata clear`, and must not. A definition file that
*! clears can only ever be loaded first: whatever was loaded before it is
*! silently wiped. Clearing belongs to the caller, and every test here opens
*! with `clear all`, which clears Mata too -- verified, not assumed. See
*! E-072, 2026-09-18.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// ---------------------------------------------------------------------------
// The sample-statistic index map (Gate 1, spec Section 4.2).
//
// Built from C alone -- the per-variable threshold counts. No data, no
// estimates, no model. That is what lets Module 0 construct it before Stage 1
// runs, and lets every later module take it as an input rather than rebuild
// it. Two reconstructions in two places is how the copies drift apart.
//
// Ordering is lavaan's, verified against lavInspect(fit, "gamma"): all
// thresholds first, grouped by variable and ascending within variable, then
// all pairs in column-major lower-triangle order.
// ---------------------------------------------------------------------------
struct dwlsmv_statmap {
    real scalar      p        // number of observed variables
    real colvector   K        // p x 1  categories per variable
    real colvector   C        // p x 1  thresholds per variable, K :- 1
    real colvector   Coff     // p x 1  prefix sum of C, Coff[j] = sum_{m<j} C[m]
    real scalar      T        // total thresholds, sum(C)
    real scalar      P        // number of pairs, p(p-1)/2
    real scalar      D        // length of s, T + P

    // one row per element of s, all D x 1
    real colvector   stat_type   // 1 = threshold, 2 = polychoric
    real colvector   var_i       // first variable
    real colvector   var_j       // second variable; 0 for a threshold
    real colvector   thr_k       // threshold index; 0 for a polychoric
    string colvector label       // "y1|t1" / "y1~~y2" -- the harness's
                                 // alignment key. Section 4.2 requires
                                 // comparison by name, never by position.
}

// ---------------------------------------------------------------------------
// One pair's Stage 2 output (Gate 3, spec Section 4.5.2).
//
// Module 2 fills this; Module 4 consumes it. The split is deliberate: Module 2
// ESTIMATES, Module 4 ASSEMBLES, and Module 4 is the only place the threshold
// influence functions and these quantities are combined.
//
// No numerical second derivatives appear anywhere in producing this. The
// information quantities are closed-form in the cell probabilities and their
// first derivatives. Reaching for a numerical Hessian of the bivariate
// log-likelihood means something upstream went wrong -- re-read Section 4.5.2.
//
// The three informations are EMPIRICAL (outer-product / BHHH): weighted by the
// observed cell proportions, not by the model-implied ones. That is lavaan's
// construction and it is what makes Gamma match lavInspect(fit, "gamma");
// weighting by pi instead gives the expected information the Gate 3 result is
// written with, which differs at O(N^-1/2) -- 6.56 relative on case_D's
// correlation block. dwlsmv_pairinfo() takes the weight matrix as an argument
// so both stay reachable, and the two coincide exactly for a saturated table,
// which is every binary pair. No field records which was used: the shipped
// path is one call site, and a flag here would be a second copy of a fact the
// caller already knows.
// ---------------------------------------------------------------------------
struct dwlsmv_pairfit {
    real scalar    i, j       // which variables, i < j
    real scalar    rho        // the polychoric estimate
    real matrix    u          // K_i x K_j, the per-cell score contribution
    real scalar    I_rr       // information, rho on rho
    real colvector I_ra       // C_i x 1, rho on variable i's thresholds
    real colvector I_rb       // C_j x 1, rho on variable j's thresholds
    real scalar    status     // 1..8, spec Module 2's status table; 1 = root
    real scalar    converged  // 0/1, set from status (status == 1), never independently
    real scalar    iterations // evaluations, the endpoint check's included
}

// ---------------------------------------------------------------------------
// The free-parameter map (Gate 2, spec Section 4.8).
//
// `fixed` is carried, not dropped. Sorbom modification indices (v2, Module 12)
// need derivatives with respect to currently-FIXED parameters, and retrofitting
// that into a table shaped around free parameters only is the expensive kind of
// change. Carrying it now costs one column.
//
// Theta does NOT appear. Under delta, theta_jj = 1 - (Lambda Phi Lambda')_jj is
// derived, and a free residual variance is unidentified from correlations
// alone. A request to free one is an error, not an option.
//
// Neither does Phi_eta, and ptype 2 is PSI for that reason. This comment read
// "2 phi" until 2026-09-21, which was stale rather than a second convention:
// spec Section 4.8 derives Phi_eta = (I-B)^-1 Psi (I-B)^-T for every factor and
// lists total factor variances as DERIVED, never free. The free (co)variances
// are Psi's -- disturbance variances for endogenous factors, and the same
// number as Phi only when a factor is exogenous. Writing Module 3 against a
// code labelled `phi` is how "unit variance" gets read as a statement about
// Phi again, which is the error that made case E necessary. Nothing had yet
// been written against the old label: ptype had zero uses outside this file.
// ---------------------------------------------------------------------------
struct dwlsmv_parmap {
    real scalar      npar     // number of FREE parameters
    real scalar      nall     // free + fixed

    string colvector name     // nall x 1, human-readable
    real colvector   ptype    // 1 lambda, 2 psi, 3 beta, 4 tau,
                              // 5 nu (y, 1), 6 omega (y, 1), 7 alpha (F, 1)
    real colvector   row      // position in its own matrix
    real colvector   col
    real colvector   free     // 1 free, 0 fixed
    real colvector   value    // start value if free, fixed value if not
    real colvector   fidx     // index into theta if free, 0 if fixed;
                              // tied members share one fidx (spec 4.8 G2)
    real colvector   grp      // nall x 1, the row's group 1..G (spec 4.8 G6:
                              // the key is (ptype, grp, row, col)). Module 3
                              // fills it; Modules 5/6 read it only when more
                              // than one map is passed (G = 1 never reads it)
}

// ---------------------------------------------------------------------------
// Stage 1's complete output. Model-free by construction.
//
// Nothing model-implied may enter this struct. That boundary is what makes the
// paired validation design of Section 7.1 meaningful: Stage 1 can be checked
// against lavaan without any structural model being specified at all.
// ---------------------------------------------------------------------------
struct dwlsmv_stage1 {
    struct dwlsmv_statmap scalar map
    real colvector  s         // D x 1, the sample statistics
    real matrix     Gamma     // D x D, N times the asymptotic covariance
    real scalar     N
    real colvector  psibar    // D x 1, the centring diagnostic of Section 4.5.5.
                              // Nonzero by construction at the rho entry of a
                              // one-empty-cell pair (O6; spec 4.5.5).
                              // Carried, not asserted on: the naive influence
                              // function is ALSO exactly centred, so a small
                              // psibar is not evidence the propagation term is
                              // present. It catches derivative errors only.
    real colvector  bound     // P x 1, pair order: 1 if a status-2 pair was
                              // included at +-B under keep_bound, else 0.
                              // Stored, not derived: a status-1 root may be
                              // accepted at +-B itself (spec 5.1, 2026-10-03).
}

// ---------------------------------------------------------------------------
// Fitted model. Modules 7-9.
// ---------------------------------------------------------------------------
struct dwlsmv_fit {
    struct dwlsmv_stage1 vector stage1   // one per group (spec 4.8 O14); G = 1: a 1-vector
    struct dwlsmv_parmap scalar par

    real colvector  theta     // npar x 1
    real matrix     acov      // npar x npar, the sandwich
    real matrix     Delta     // D x npar
    real scalar     fmin      // the DWLS objective at theta-hat
    real scalar     converged
    real scalar     iterations

    // Section 4.9. Three statistics share the machinery; T3 is the v1 target.
    // T2 and T3 target DIFFERENT chi-square distributions -- report the df
    // actually used, never D in both cases.
    real scalar     T_raw     // the unadjusted quadratic form, (N-G) F-hat, N the
                              // total; G = 1: (N-1) F-hat, lavaan's scaling, not
                              // Mplus's N F-hat (spec 4.9; MG3, decided 2026-10-09)
    real scalar     T3        // a*T_raw + b, the reported statistic
    real scalar     T3_df     // equals D_stat for T3; still stored explicitly
    real scalar     shift_a
    real scalar     shift_b
    real scalar     D_stat    // degrees of freedom, D_s - npar
}

end
