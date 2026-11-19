*! Executed 2026-11-15 by the main session (StataMP 18) after the re-base for conditional D11 recovery: passes on its first run. Re-base drafted by a read-only agent 2026-11-15.
*! Executed 2026-11-14 by the main session (Stata 18 MP), before the re-base; 7 of 7 mutants run are killed. Formerly: NOT YET EXECUTED. Drafted by a read-only agent 2026-11-14.
*! Conditional mode under the marker metric: admissibility is checked and reported, and a fit that does not converge is
*! qualified, with no standard errors shown (the owner's requirements of 2026-11-14; the suppression of T3 and p confirmed by the owner the same day).
*!
*! Needs edits.md of the same date (dwlsmv_cx_admissnote.mata, the edits to dwlsmv_m11_estimate.mata, dwlsmv_m11_fit.mata,
*! dwlsmv.ado). Simulated data only, fixed seeds, no fixture files. The designs and seeds are those of
*! E-097: configuration D (N 500, loadings .7 .7 -.7 -.7, the reverse-keyed trap
*! design) did not converge under the marker metric on 4 of 4 seeds in probe_inadm.log (iterations 301, T3 in the hundreds,
*! factor variance below 0 or a loading in the hundreds); configuration A seed 3 did not converge either but is admissible.
*! A premise that fails to reproduce is reported as a failure of THIS test (the data design no longer provokes the
*! failure), never skipped. (Re-based 2026-11-15: on the trap design the marker fit now converges through D11's recovery,
*! which is checked as a positive property; the nonconverged-display checks use a spec recovery cannot repair, see (i).)
*!
*! Part A, Mata level (no fit): dwlsmv_cx_admissnote() on hand-built conditional tables, as test_condx_admissible.do builds
*! them: a Heywood loading, a negative disturbance variance, the closed forms of the note.
*! Part B, through the command:
*!   (i)    on the four trap seeds, the marker fit with the original spec now RECOVERS (e(recovery) recovered,
*!          e(converged) 1). The display checks need a fit recovery leaves nonconverged: the same data with a fixed
*!          nonzero gamma, (F <- x@.3), which is ineligible for recovery (test_condx_recovery.do C3). On each seed
*!          where that fit does not converge: e(converged) = 0, e(recovery) ineligible; the prominent line sits directly
*!          under the title line, in the output and on replay and after estimates restore; e(V), e(T3), e(p_T3) are not
*!          posted; no "Std. err." column and no "Prob > chi2" line; the header says no recovery refit was attempted
*!          because the model is not eligible; the admissibility note, where it appears, goes with the last-iterate
*!          paragraph. A seed whose fixed-gamma fit converges is skipped and counted; at least one seed must be
*!          exercised. (i-b): configuration A seed 3 (nonconverged first fit, admissible) now triggers recovery on
*!          the nonconvergence: recovered (standard errors shown) or failed (qualified all the same).
*!   (ii)   a converged, admissible conditional fit: no admissibility note, standard errors, T3 and p, and the "not computed"
*!          sentence names only the aliasing check (the diagnostics of the weight and derivative matrices are computed since 2026-11-15).
*!   (iii)  Part A, plus Part A2: a deterministic converged INADMISSIBLE conditional fit through dwlsmv_m11_estimate() (the
*!          driver the command calls), built from an exact-fit Heywood population (s = sigma(theta0) with theta0 inadmissible,
*!          so dwls_fit converges to theta0 with F = 0), under metric unitvariance (lambda_1 = 1.2) and metric marker
*!          (lambda_1 = 1, psi = 1.44). Since 2026-11-15 the driver runs D11 on a conditional marker fit, which reads the
*!          external dwlsmv_m11_st, so the test creates it first: under marker the flagged converged fit is recovered
*!          (the population is its own unit-variance solution; admissibility stays flagged, spec D11), under
*!          unitvariance rec is none; plus (iii-b) a converged INADMISSIBLE fit through the command, searched for under
*!          metric(unitvariance) (a Heywood loading); if none is found in 60 seeds the test says so and the command-level case
*!          is NOT exercised (stated in the closing line); Part A2 does not depend on that search.
*!   (iv)   the same nonconverging data under metric(unitvariance): converged, no qualification, standard errors present.
*!   (v)    an ordinary (no covariate) fit with e(converged) set to 0 by hand (the only way to force it, spec D8: no
*!          iteration limit option exists): since 2026-11-18 displayed as a conditional one, without its covariate lines.
*!
*! WHAT THIS CERTIFIES AND DOES NOT. It certifies the output path (what is displayed and posted). It says nothing about the
*! estimates, about whether metric(unitvariance) cures any particular failure beyond these data, or about Module 10's
*! diagnostics on a conditional fit.
*!
*! MUTANTS for the main session (edit, run, revert; each must turn this file red). Files are under src/.
*!   M1  mata/dwlsmv_m11_estimate.mata: `A     = dwlsmv_admissible(fit.theta, par, st.map)` -> `A     = J(6, 4, 0)`.
*!       Killed by Part A2 (A[1, 1] would be 0 and no recovery would be triggered under marker) and, where a fixed-gamma
*!       fit's last iterate is inadmissible, by (i). (Before 2026-11-15 it was killed by (i) D seed 1 alone.)
*!   M2  mata/dwlsmv_m11_estimate.mata: `A     = dwlsmv_admissible(fit.theta, par, st.map)` -> `A     = (fit.converged ? J(6, 4, 0) : dwlsmv_admissible(fit.theta, par, st.map))`
*!       (a converged fit bypasses the check). Killed by Part A2 (converged, inadmissible: A[1, 1] would be 0 and the note
*!       empty), under both metrics, whether or not the (iii-b) search finds a fit.
*!   M3  RETIRED 2026-11-15 (was dwlsmv_cx_admissnote.mata: `if (metric == "marker") {` -> `if (metric == "none") {`, killed
*!       by the recovery sentence). The sentence no longer exists; D11 recovery runs for a conditional marker fit.
*!   M4  mata/dwlsmv_cx_admissnote.mata: delete the `if (converged != 1) {` block. Killed by Part A D (the closing paragraph
*!       about standard errors survives for a nonconverged fit) and, where a note appears, by (i).
*!   M5  dwlsmv.ado, header, new prominent block: `e(converged) == 0 {` (the block that prints "convergence not achieved:
*!       the optimizer stopped after") -> `e(converged) == 2 {`. Killed by (i), the replay and the restore.
*!   M6  dwlsmv.ado, posting: delete the `if `converged' == 0 {` branch (since 2026-11-18 every fit; was covariate fits only)
*!       so that e(V) is posted. Killed by (i): rows(e(V)) != 0 and a "Std. err." column.
*!   M7  dwlsmv.ado, scalars: `if `converged' != 0 {` -> `if 1 {` (was the covariate condition). Killed by (i): e(T3) is
*!       posted.
*!   M8  RETIRED 2026-11-18: dropping the covariate condition from the header block is now the rule (owner: non-convergence
*!       is an inference boundary for every fit); (v) checks the uniform display instead.
*!   M9  mata/dwlsmv_m11_fit.mata: restore the old two-line "Not computed ... the admissibility" text. Killed by (ii).
*!   M10 (added 2026-11-15) this file, tm_hit: delete the `if "`phrase'" == "Std. err." {` block (17+ spelling only).
*!       Killed by the self-check before Part B (a "Std. Err." header gives hit 0).
*!   M11 (new 2026-11-15, NOT YET RUN) mata/dwlsmv_m11_estimate.mata: make the recovery call conditional on `!cx`
*!       (a conditional marker fit takes no recovery). Killed by the positive check in (i) (four trap seeds not
*!       recovered), by Part A2 marker (rec none) and by (i-b).
*!   M12 (new 2026-11-15, NOT YET RUN) dwlsmv.ado, header: swap the texts of the `ineligible` and `failed` branches.
*!       Killed by (i) (the "not eligible" header lines) and, if (i-b) fails to recover, by (i-b).
*!   M13 (new 2026-11-15, NOT YET RUN) dwlsmv.ado, header: restore the old two lines "No recovery refit is attempted /
*!       for a model with a covariate." Killed by (i) (hasnt of the old sentence; has of the new lines).
*!
*! CHANGED 2026-11-15 after the Stata 16.1 milestone (d5d1549): 16.1 prints "Std. Err." where 17+ prints "Std. err.",
*! so the six "has" checks failed there and the "hasnt" checks could not see a 16.1 column. "Std. err." now matches
*! either spelling (tm_hit), with a self-check on hand-written headers. A test-text change only; no estimator code.
*!
*! CHANGED 2026-11-15 (recovery): D11 marker recovery now runs for a conditional (one-covariate) fit (the owner's
*! decision of 2026-11-15; working tree, uncommitted at drafting). This file asserted the old behaviour (no recovery
*! refit for a covariate) and failed 11 checks. Re-based, with every property kept: Part A, the sentence "No recovery
*! refit was attempted for a model with a covariate" is wanted ABSENT (0), mutant M3 retired; Part A2, the external
*! dwlsmv_m11_st is created around the direct call to dwlsmv_m11_estimate() and rec is wanted "recovered" under
*! marker, "none" under unitvariance; Part B (i), the premise (the trap seeds end nonconverged under marker) is the
*! very failure recovery repairs, so it becomes a positive check (recovered, converged) and the nonconverged-display
*! checks move to the same data with a fixed nonzero gamma (F <- x@.3), which is ineligible, with the premise checked
*! and a seed whose fit converges skipped (at least one required); (i-b) and (iii-b) re-based on e(recovery). The
*! fixed-gamma spec and the seeds were chosen for this test after test_condx_recovery.do C3 and before any run of this
*! file; whether those fits stay nonconverged is measured here, not assumed.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

global TM_NFAIL 0

* =====================================================================================================
* Part A: dwlsmv_cx_admissnote() at the Mata level. Sources are loaded, types first (as test_condx_admissible.do).
* =====================================================================================================
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_condx_marker: `f' failed to load"
        exit 9
    }
}
capture mata: (void) dwlsmv_cx_admissnote(J(6, 4, 0), "y" \ "z", "F" \ "G", "marker", 1)
if _rc {
    display as error "test_condx_marker: dwlsmv_cx_admissnote() is not defined (src/mata/dwlsmv_cx_admissnote.mata missing?)"
    exit 9
}

mata:

real scalar tm_fail(string scalar what)
{
    printf("{err}test_condx_marker: %s\n", what)
    return(1)
}

// the note msg contains (want = 1) or does not contain (want = 0) the phrase; returns the number of failures
real scalar tm_ck(string scalar msg, string scalar ph, real scalar want, string scalar tag)
{
    real scalar hit

    hit = (strpos(msg, ph) > 0)
    if (hit != want) {
        if (want) {
            return(tm_fail(tag + ": the note lacks: " + ph))
        }
        return(tm_fail(tag + ": the note contains: " + ph))
    }
    return(0)
}

real scalar tm_near(real scalar a, real scalar b, real scalar tol)
{
    if (missing(a) | missing(b)) {
        return(0)
    }
    return(abs(a - b) <= tol)
}

// A hand-built parameter table (as test_condx_admissible.do's ca_mk): free rows get theta in table order
struct dwlsmv_parmap scalar tm_mk(real colvector pt, real colvector rw, real colvector cl,
                                  real colvector fr, real colvector vl)
{
    struct dwlsmv_parmap scalar par
    real scalar i, n, k

    n = rows(pt)
    par.nall  = n
    par.npar  = sum(fr)
    par.name  = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        par.name[i] = "row" + strofreal(i)
    }
    par.ptype = pt
    par.row   = rw
    par.col   = cl
    par.free  = fr
    par.value = vl
    par.fidx  = J(n, 1, 0)
    k = 0
    for (i = 1; i <= n; i++) {
        if (fr[i]) {
            k = k + 1
            par.fidx[i] = k
        }
    }
    par.grp = J(n, 1, 1)
    return(par)
}

// One factor, marker metric, K = (3,3,2): tau (5) | lambda_1 fixed 1, lambda_2, lambda_3 free | psi free | gamma F x |
// kappa y3 x.   theta = (tau, lambda_2, lambda_3, psi, gamma, kappa), 10 free parameters.
struct dwlsmv_parmap scalar tm_par1m()
{
    real colvector pt, rw, cl, fr, vl

    pt = (J(5, 1, 4) \ J(3, 1, 1) \ 2 \ 8 \ 9)
    rw = ((1\1\2\2\3) \ (1\2\3) \ 1 \ 1 \ 3)
    cl = ((1\2\1\2\1) \ J(3, 1, 1) \ 1 \ 1 \ 1)
    fr = (J(5, 1, 1) \ (0\1\1) \ 1 \ 1 \ 1)
    vl = (J(5, 1, 0) \ (1\0\0) \ 0 \ 0 \ 0)
    return(tm_mk(pt, rw, cl, fr, vl))
}

real scalar tm_a()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap scalar map
    real matrix A
    real colvector tau
    string colvector ind, fac
    string scalar m, m2
    real scalar nf

    nf  = 0
    par = tm_par1m()
    map = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
    tau = (-0.5 \ 0.5 \ -0.3 \ 0.6 \ 0.2)
    ind = ("y1" \ "y2" \ "y3")
    fac = ("F")

    // ---- an admissible point: no note, under either metric, converged or not
    A = dwlsmv_admissible((tau \ 0.8 \ 0.7 \ 0.5 \ 0.4 \ 0.25), par, map)
    if (any(A[., 1] :== 1)) {
        nf = nf + tm_fail("A: the admissible point is flagged: " + invtokens(strofreal(A[., 1])'))
    }
    if (dwlsmv_cx_admissnote(A, ind, fac, "marker", 1) != "" | dwlsmv_cx_admissnote(A, ind, fac, "marker", 0) != "" |
        dwlsmv_cx_admissnote(A, ind, fac, "unitvariance", 1) != "") {
        nf = nf + tm_fail("A: an admissible point has a note")
    }

    // ---- B: Heywood at the marker, psi = 1.44, so Theta_1 = 1 - 1.44 = -0.44; converged, marker metric
    A = dwlsmv_admissible((tau \ 0.8 \ 0.7 \ 1.44 \ 0.4 \ 0.25), par, map)
    if (A[1, 1] != 1 | A[1, 3] != 1 | !tm_near(A[1, 2], -0.44, 1e-12)) {
        nf = nf + tm_fail("B: the table's Heywood case is not flagged at item 1 with -0.44")
    }
    m = dwlsmv_cx_admissnote(A, ind, fac, "marker", 1)
    nf = nf + tm_ck(m, "The fitted solution is not admissible", 1, "B")
    nf = nf + tm_ck(m, "y1: residual variance -0.4400, below 0", 1, "B")
    nf = nf + tm_ck(m, "Standard errors and the WLSMV test", 1, "B")
    nf = nf + tm_ck(m, "No recovery refit was attempted for a model with a covariate", 0, "B")
    nf = nf + tm_ck(m, "metric(unitvariance)", 0, "B")   // the sentence that named it is gone (2026-11-15)
    nf = nf + tm_ck(m, "last iterate", 0, "B")

    // ---- C: the same point under metric(unitvariance): the note, and no recovery sentence (none under either metric since 2026-11-15)
    m2 = dwlsmv_cx_admissnote(A, ind, fac, "unitvariance", 1)
    nf = nf + tm_ck(m2, "y1: residual variance -0.4400, below 0", 1, "C")
    nf = nf + tm_ck(m2, "Standard errors and the WLSMV test", 1, "C")
    nf = nf + tm_ck(m2, "No recovery refit", 0, "C")

    // ---- D: the same point, not converged: standard-error paragraph replaced, no recovery sentence (2026-11-15)
    m2 = dwlsmv_cx_admissnote(A, ind, fac, "marker", 0)
    nf = nf + tm_ck(m2, "y1: residual variance -0.4400, below 0", 1, "D")
    nf = nf + tm_ck(m2, "last iterate of a fit that did not converge", 1, "D")
    nf = nf + tm_ck(m2, "Standard errors and the WLSMV test", 0, "D")
    nf = nf + tm_ck(m2, "Point estimates are reported", 0, "D")
    nf = nf + tm_ck(m2, "No recovery refit was attempted for a model with a covariate", 0, "D")

    // ---- E: a negative disturbance variance, psi = -0.3: Phi = psi, minimum eigenvalue -0.3 (row 3)
    A = dwlsmv_admissible((tau \ 0.8 \ 0.7 \ -0.3 \ 0.4 \ 0.25), par, map)
    if (A[3, 1] != 1 | !tm_near(A[3, 2], -0.3, 1e-14)) {
        nf = nf + tm_fail("E: psi = -0.3 is not flagged in row 3 with -0.3; got flag " + strofreal(A[3, 1]) + ", value " + strofreal(A[3, 2]))
    }
    m = dwlsmv_cx_admissnote(A, ind, fac, "marker", 1)
    nf = nf + tm_ck(m, "smallest eigenvalue -0.3000", 1, "E")
    nf = nf + tm_ck(m, "No recovery refit was attempted for a model with a covariate", 0, "E")

    // ---- F: the slopes do not change the note (gamma, kappa are not conditions)
    A  = dwlsmv_admissible((tau \ 0.8 \ 0.7 \ 1.44 \ 5 \ -3), par, map)
    m2 = dwlsmv_cx_admissnote(A, ind, fac, "marker", 1)
    A  = dwlsmv_admissible((tau \ 0.8 \ 0.7 \ 1.44 \ 0.4 \ 0.25), par, map)
    if (m2 != dwlsmv_cx_admissnote(A, ind, fac, "marker", 1)) {
        nf = nf + tm_fail("F: the note changed with (gamma, kappa)")
    }

    // ---- G: the indicator named is the right one with gamma and kappa rows in the table: lambda_3 = 3, psi = 0.5
    A = dwlsmv_admissible((tau \ 0.8 \ 3 \ 0.5 \ 0.4 \ 0.25), par, map)
    m = dwlsmv_cx_admissnote(A, ind, fac, "marker", 1)
    nf = nf + tm_ck(m, "y3: residual variance -3.5000, below 0", 1, "G")
    nf = nf + tm_ck(m, "y1: residual variance", 0, "G")
    return(nf)
}

// Part A2, one case. A conditional Stage 1 object whose statistics are s = sigma(theta0) EXACTLY (Gamma the identity, as
// test_condx_fit.do's cf_st; N 1000) for an INADMISSIBLE theta0, then the command's driver, dwlsmv_m11_estimate(), with the
// canonical text as the command passes it (";"-separated). Exact fit: dwls_fit converges to theta0 with F = 0 (the bound
// 1e-16 is test_condx_fit.do's check 9 bound for the same kind of population), so the solution is converged AND
// inadmissible. Nothing in this path needs an admissible theta0: the implied correlations are lambda_j lambda_k psi (all
// below 1 in absolute value here) and Theta_j enters no statistic.
// Since 2026-11-15 the driver runs D11's recovery on a conditional marker fit, and the recovery re-runs Stage 1's output
// from the external dwlsmv_m11_st (not Stage 1 itself), so the external is created here from st and removed afterwards
// (the command's own cleanup, dwlsmv_m11_clear(), is not run on this path). Under marker the converged inadmissible fit
// fires the trigger (an admissibility flag) and is recovered: the population is its own unit-variance solution, and the
// recovered fit has the same Heywood loading transformed (admissibility stays flagged, spec D11), so A and the note below
// are unchanged. Under unitvariance no recovery exists, so rec is none.
real scalar tm_bcase(string rowvector lines, string scalar metric, real colvector ths, string scalar tag)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, par0
    struct dwlsmv_fit scalar fit
    struct dwlsmv_statmap vector mv
    real matrix A, Bd
    real colvector K
    string scalar canon, rec, rnote, m
    real scalar nf, i, D
    pointer() scalar pst

    nf    = 0
    K     = (3 \ 3 \ 2)
    canon = ""
    for (i = 1; i <= cols(lines); i++) {
        if (i > 1) {
            canon = canon + ";"
        }
        canon = canon + lines[i]
    }
    st.map = dwlsmv_statmap_build(K, 1)
    D      = st.map.D
    par0   = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    if (par0.npar != rows(ths)) {
        return(tm_fail(tag + ": npar is " + strofreal(par0.npar) + ", the constructed theta has " + strofreal(rows(ths))))
    }
    mv = dwlsmv_statmap(1, 1)
    mv[1] = st.map
    st.s      = dwls_model(ths, par0, mv)
    st.Gamma  = I(D)
    st.N      = 1000
    st.psibar = J(D, 1, 0)
    st.bound  = J(st.map.P, 1, 0)

    // the external D11's recovery reads (dwlsmv_m11_recrun), as dwlsmv_m11_stage1() leaves it
    rmexternal("dwlsmv_m11_st")
    rmexternal("dwlsmv_m11_rec")
    pst = crexternal("dwlsmv_m11_st")
    *pst = st
    dwlsmv_m11_estimate(canon, metric, st, par, fit, A, Bd, rec, rnote)
    rmexternal("dwlsmv_m11_st")
    rmexternal("dwlsmv_m11_rec")

    if (fit.converged != 1) {
        nf = nf + tm_fail(tag + ": the exact-fit population did not converge")
    }
    if (!(fit.fmin <= 1e-16)) {
        nf = nf + tm_fail(tag + ": F is " + strofreal(fit.fmin) + " on an exact-fit population, not 0")
    }
    if (rows(A) != 6 | cols(A) != 4) {
        return(nf + tm_fail(tag + ": A is not 6 x 4 (the conditional branch gave " + strofreal(rows(A)) + " x " + strofreal(cols(A)) + ")"))
    }
    if (A[1, 1] != 1 | A[1, 3] != 1 | !tm_near(A[1, 2], -0.44, 1e-9)) {
        nf = nf + tm_fail(tag + ": the converged Heywood fit is not flagged at item 1 with -0.44 (flag " + strofreal(A[1, 1]) +
                          ", item " + strofreal(A[1, 3]) + ", value " + strofreal(A[1, 2]) + "): the check was bypassed")
    }
    if (metric == "marker") {
        if (rec != "recovered") {
            nf = nf + tm_fail(tag + ": rec is " + rec + ", wanted recovered (the exact-fit Heywood population is its own unit-variance solution); note: " + rnote)
        }
        nf = nf + tm_ck(rnote, "restarted from the unit-variance solution", 1, tag + " (recovery note)")
    }
    else {
        if (rec != "none") {
            nf = nf + tm_fail(tag + ": rec is " + rec + ", wanted none under metric unitvariance")
        }
        if (rnote != "") {
            nf = nf + tm_fail(tag + ": a recovery note under metric unitvariance: " + rnote)
        }
    }
    m = dwlsmv_cx_admissnote(A, ("y1" \ "y2" \ "y3"), ("F"), metric, fit.converged)
    nf = nf + tm_ck(m, "The fitted solution is not admissible", 1, tag)
    nf = nf + tm_ck(m, "y1: residual variance -0.4400, below 0", 1, tag)
    nf = nf + tm_ck(m, "Standard errors and the WLSMV test", 1, tag)
    nf = nf + tm_ck(m, "last iterate", 0, tag)
    nf = nf + tm_ck(m, "No recovery refit was attempted for a model with a covariate", 0, tag)
    return(nf)
}

// Part A2: both metrics. Unitvariance: lambda = (1.2, .7, .6), gamma .4, so Theta_1 = 1 - 1.44. Marker: lambda_1 = 1 fixed,
// lambda_2 = .5, lambda_3 = .4, psi = 1.44, gamma .4, so Theta_1 = 1 - 1.44 as well (implied correlations 0.72, 0.576, 0.288).
real scalar tm_a2()
{
    real colvector tau
    string rowvector lu, lm
    real scalar nf

    nf  = 0
    tau = (-0.5 \ 0.5 \ -0.3 \ 0.6 \ 0.2)
    lu  = ("factor F", "ind y1 y2 y3", "lambda F y1 free", "lambda F y2 free", "lambda F y3 free",
           "metric unitvariance F", "covariate x", "gamma F x free")
    lm  = ("factor F", "ind y1 y2 y3", "lambda F y2 free", "lambda F y3 free", "psi F F free",
           "metric marker F y1", "covariate x", "gamma F x free")
    nf = nf + tm_bcase(lu, "unitvariance", (tau \ 1.2 \ 0.7 \ 0.6 \ 0.4), "A2 unitvariance")
    nf = nf + tm_bcase(lm, "marker", (tau \ 0.5 \ 0.4 \ 0.4 \ 1.44), "A2 marker")   // theta order tau, lambda, gamma, psi (spec 4.8); the draft had gamma and psi swapped
    return(nf)
}

end

capture noisily mata: st_local("nA", strofreal(tm_a()))
if _rc {
    display as error "test_condx_marker: Part A aborted, rc " _rc
    global TM_NFAIL = ${TM_NFAIL} + 1
}
else {
    display as text "  Part A (dwlsmv_cx_admissnote at the Mata level): `nA' failure(s)"
    global TM_NFAIL = ${TM_NFAIL} + `nA'
}

* =====================================================================================================
* Part B: through the command
* =====================================================================================================
capture noisily mata: st_local("nA2", strofreal(tm_a2()))
if _rc {
    display as error "test_condx_marker: Part A2 aborted, rc " _rc
    global TM_NFAIL = ${TM_NFAIL} + 1
}
else {
    display as text "  Part A2 (a converged inadmissible conditional fit through dwlsmv_m11_estimate): `nA2' failure(s)"
    global TM_NFAIL = ${TM_NFAIL} + `nA2'
}

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

mata:

// the whole text of a log, lines joined by one space
string scalar tm_txt(string scalar fn)
{
    return(invtokens(cat(fn)', " "))
}

real scalar tm_has(string scalar fn, string scalar ph)
{
    return(strpos(tm_txt(fn), ph) > 0)
}

// does some line equal ph exactly (after trimming)?
real scalar tm_hasline(string scalar fn, string scalar ph)
{
    string colvector L
    real scalar i

    L = cat(fn)
    for (i = 1; i <= rows(L); i++) {
        if (strtrim(L[i]) == ph) {
            return(1)
        }
    }
    return(0)
}

// the number of the first line containing ph, 0 if none
real scalar tm_idx(string scalar fn, string scalar ph)
{
    string colvector L
    real scalar i

    L = cat(fn)
    for (i = 1; i <= rows(L); i++) {
        if (strpos(L[i], ph) > 0) {
            return(i)
        }
    }
    return(0)
}

// is the line right after the first line containing p1 a line containing p2?
real scalar tm_next(string scalar fn, string scalar p1, string scalar p2)
{
    string colvector L
    real scalar i

    L = cat(fn)
    i = tm_idx(fn, p1)
    if (i == 0 | i == rows(L)) {
        return(0)
    }
    return(strpos(L[i + 1], p2) > 0)
}

end

* run a command with its output written to the log file fn; r(rc) is its return code
capture program drop tm_run
program tm_run, rclass
    gettoken fn 0 : 0
    quietly log using "`fn'", text replace name(tmlog)
    capture noisily `0'
    local rc = _rc
    quietly log close tmlog
    return scalar rc = `rc'
end

* the log fn has (kind has) or lacks (kind hasnt) a phrase.
* The coefficient table's header reads "Std. Err." up to Stata 16 and "Std. err." from Stata 17 (seen 2026-11-15:
* the 16.1 milestone printed "Std. Err." and failed the six "has" checks below, while every "hasnt" check passed
* without being able to see a 16.1 column). So the phrase "Std. err." stands for EITHER spelling: has = either is
* present, hasnt = neither is. The self-check after this program exercises the 16.1 spelling on any Stata.
capture program drop tm_hit
program tm_hit, rclass
    args fn phrase
    mata: st_local("hit", strofreal(tm_has(st_local("fn"), st_local("phrase"))))
    if "`phrase'" == "Std. err." {
        mata: st_local("hit2", strofreal(tm_has(st_local("fn"), "Std. Err.")))
        local hit = `hit' | `hit2'
    }
    return scalar hit = `hit'
end

capture program drop tm_expect
program tm_expect
    args fn kind phrase tag
    tm_hit "`fn'" "`phrase'"
    local hit = r(hit)
    if ("`kind'" == "has" & `hit' == 0) | ("`kind'" == "hasnt" & `hit' == 1) {
        display as error "test_condx_marker: `tag': the output " cond("`kind'" == "has", "lacks", "contains") ": `phrase'"
        global TM_NFAIL = ${TM_NFAIL} + 1
    }
end

* one condition, with the failure counted
capture program drop tm_cond
program tm_cond
    args ok tag
    if `ok' != 1 {
        display as error "test_condx_marker: `tag'"
        global TM_NFAIL = ${TM_NFAIL} + 1
    }
end

* the data of probe_inadm.do (pi_data): N, four loadings, seed; factor shifted .5 by x, direct effect +.3 of x on y4*,
* three categories (thresholds -.5, .5)
capture program drop tm_data
program tm_data
    args n l1 l2 l3 l4 seed
    clear
    set seed `seed'
    quietly set obs `n'
    quietly generate byte x = runiform() < .5
    quietly generate double eta = .5 * x + rnormal()
    local j 0
    foreach l in `l1' `l2' `l3' `l4' {
        local ++j
        local sd = sqrt(1 - (`l')^2)
        quietly generate double z`j' = `l' * eta + `sd' * rnormal() + (`j' == 4) * .3 * x
        quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
    }
    drop eta z1 z2 z3 z4
end

* the data of test_condx_output.do: N 1500, loading .7, residual variance .51, factor shifted .4 by x, direct effect +.3
capture program drop tm_data2
program tm_data2
    clear
    set seed 20261114
    quietly set obs 1500
    quietly generate byte x = runiform() < .5
    quietly generate double eta = .4 * x + rnormal()
    forvalues j = 1/4 {
        quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
        quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
    }
    drop eta z1 z2 z3 z4
end

* sets e(converged) to 0 on the posted results, as test_dwlsmv_gaps.do G10 does
capture program drop tm_unconv
program tm_unconv, eclass
    ereturn scalar converged = 0
end

tempfile LG LG2 LG3

* self-check of the matching tm_expect uses, on hand-written header lines (no fit): the 16.1 spelling and the 17+
* spelling both count as "Std. err.", a different column does not
tempname fh
foreach sp in "Std. Err." "Std. err." "Std. dev." {
    quietly file open `fh' using "`LG'", write text replace
    file write `fh' "             |      Coef.   `sp'      z    P>|z|" _n
    file close `fh'
    tm_hit "`LG'" "Std. err."
    local want = ("`sp'" != "Std. dev.")
    tm_cond `=(r(hit) == `want')' "self-check: a header with `sp' gives hit `=r(hit)', wanted `want'"
}

local mB "(F -> y1 y2 y3 y4) (F <- x) (y4 <- x)"
local mBg "(F -> y1 y2 y3 y4) (F <- x@.3) (y4 <- x)"
local mA "(F -> y1 y2 y3 y4) (F <- x)"
local title "DWLS estimation, WLSMV test"
local first 1
local nexi 0
local nskip 0
local nnote 0

* ---- (i) the trap design, four seeds (probe_inadm.log): the marker fit did not converge before D11's recovery ran for a
* conditional fit (2026-11-15). Now (i-pos) it is recovered; the nonconverged-display checks (i-nc) use the same data with a
* fixed nonzero gamma, which recovery refuses (ineligible); then (iv) unitvariance ----
foreach seed in 20262415 20262516 20262617 20262718 {
    display as text _n "  (i) and (iv): configuration D, seed `seed'"
    tm_data 500 .70 .70 -.70 -.70 `seed'

    * (i-pos) the original spec: the first fit does not converge, recovery repairs it
    tm_run "`LG'" dwlsmv `mB'
    local rc = r(rc)
    tm_cond `=`rc' == 0' "(i-pos) seed `seed': the marker fit gave rc `rc', wanted 0"
    if !`rc' {
        tm_cond `=(e(converged) == 1 & "`e(recovery)'" == "recovered")' "(i-pos) seed `seed': e(converged) is `=e(converged)' and e(recovery) is `e(recovery)', wanted 1 and recovered"
        tm_expect "`LG'" hasnt "convergence not achieved" "(i-pos) seed `seed'"
        tm_expect "`LG'" hasnt "No recovery refit" "(i-pos) seed `seed'"
        tm_expect "`LG'" has "Std. err." "(i-pos) seed `seed'"
    }

    * (i-nc) a fixed nonzero gamma: ineligible for recovery; the first fit must still not converge (premise, checked here)
    tm_run "`LG'" dwlsmv `mBg'
    local rc = r(rc)
    tm_cond `=`rc' == 0' "(i-nc) seed `seed': the fixed-gamma marker fit gave rc `rc', wanted 0"
    local nc 0
    if !`rc' {
        if e(converged) == 0 {
            local nc 1
        }
        else {
            display as text "      seed `seed': the fixed-gamma fit converged (e(recovery) `e(recovery)'): the nonconverged-display checks are not exercised on this seed"
        }
    }
    if !`nc' {
        local ++nskip
    }
    if `nc' {
        local ++nexi

        * what is posted
        mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
        mata: st_local("nb", strofreal(cols(st_matrix("e(b)"))))
        tm_cond `=`nv' == 0' "(i) seed `seed': e(V) has `nv' rows; none is wanted"
        tm_cond `=`nb' == 15' "(i) seed `seed': e(b) has `nb' columns, wanted 15"
        tm_cond `=("`e(T3)'" == "")' "(i) seed `seed': e(T3) is posted"
        tm_cond `=("`e(p_T3)'" == "")' "(i) seed `seed': e(p_T3) is posted"
        tm_cond `=(e(df_T3) == 5)' "(i) seed `seed': e(df_T3) = `=e(df_T3)', wanted 5 (the fixed gamma frees one degree of freedom)"
        tm_cond `=(e(iterations) > 0 & "`e(recovery)'" == "ineligible" & "`e(covariates)'" == "x")' "(i) seed `seed': iterations, recovery (`e(recovery)', wanted ineligible) or covariates not as expected"
        tm_cond `=("`e(admissible)'" == "")' "(i) seed `seed': an admissibility indicator is posted; ordinary fits post none"

        * what is printed: the prominent block directly under the title line
        mata: st_local("nx", strofreal(tm_next(st_local("LG"), st_local("title"), "convergence not achieved: the optimizer stopped after")))
        tm_cond `nx' "(i) seed `seed': the line right under the title line is not the convergence line"
        tm_expect "`LG'" has "The estimates below are the last iterate, not a solution." "(i) seed `seed'"
        tm_expect "`LG'" has "No standard errors, T3 or p-value are reported for them." "(i) seed `seed'"
        tm_expect "`LG'" has "Refit with metric(unitvariance). No recovery refit was attempted:" "(i) seed `seed'"
        tm_expect "`LG'" has "the model is not eligible for it (see the notes)." "(i) seed `seed'"
        tm_expect "`LG'" hasnt "No recovery refit is attempted" "(i) seed `seed': the old header text"
        tm_expect "`LG'" has "gamma[F,x] is fixed at" "(i) seed `seed': the recovery note names the row"
        tm_expect "`LG'" hasnt "restarted from the unit-variance solution" "(i) seed `seed'"
        tm_expect "`LG'" has "Scaled-and-shifted test (T3): not reported (not a solution)" "(i) seed `seed'"
        tm_expect "`LG'" hasnt "Prob > chi2" "(i) seed `seed'"
        tm_expect "`LG'" hasnt "Std. err." "(i) seed `seed'"
        tm_expect "`LG'" hasnt "P>|z|" "(i) seed `seed'"
        tm_expect "`LG'" has "Last iterate, not a solution; no standard errors" "(i) seed `seed'"
        mata: st_local("nl", strofreal(tm_hasline(st_local("LG"), "convergence not achieved")))
        tm_cond `=(`nl' == 0)' "(i) seed `seed': the old closing line 'convergence not achieved' is printed as well (a duplicate)"

        * the admissibility note is computed at the last iterate. Whether a fixed-gamma last iterate is inadmissible is not
        * known in advance (the premise was measured for the free-gamma fits), so the note's companions are checked where
        * the note appears, and the number of seeds where it appears is reported
        mata: st_local("hasnote", strofreal(tm_has(st_local("LG"), "The fitted solution is not admissible")))
        if `hasnote' {
            local ++nnote
            tm_expect "`LG'" has "last iterate of a fit that did not converge" "(i) seed `seed'"
            tm_expect "`LG'" hasnt "Standard errors and the WLSMV test" "(i) seed `seed'"
            tm_expect "`LG'" hasnt "No recovery refit was attempted for a model with a covariate" "(i) seed `seed'"
        }
        tm_expect "`LG'" has "Not computed for a model with a covariate in this version: the check for thresholds" "(i) seed `seed'"
        tm_expect "`LG'" hasnt "the admissibility" "(i) seed `seed'"

        * replay: the same qualification
        tm_run "`LG2'" dwlsmv
        tm_cond `=(r(rc) == 0)' "(i) seed `seed': replay gave rc `=r(rc)'"
        mata: st_local("nx", strofreal(tm_next(st_local("LG2"), st_local("title"), "convergence not achieved: the optimizer stopped after")))
        tm_cond `nx' "(i) seed `seed': replay: the line under the title line is not the convergence line"
        tm_expect "`LG2'" has "The estimates below are the last iterate, not a solution." "(i) replay, seed `seed'"
        tm_expect "`LG2'" has "Last iterate, not a solution; no standard errors" "(i) replay, seed `seed'"
        tm_expect "`LG2'" has "the model is not eligible for it (see the notes)." "(i) replay, seed `seed'"
        tm_expect "`LG2'" hasnt "Std. err." "(i) replay, seed `seed'"
        tm_expect "`LG2'" hasnt "Prob > chi2" "(i) replay, seed `seed'"

        * estimates store and restore keep the qualification (first exercised seed only)
        if `first' {
            quietly estimates store tmnc
            quietly dwlsmv (F -> y1 y2 y3 y4), metric(unitvariance)
            quietly estimates restore tmnc
            tm_run "`LG3'" dwlsmv
            tm_expect "`LG3'" has "convergence not achieved: the optimizer stopped after" "(i) after estimates restore"
            tm_expect "`LG3'" hasnt "Std. err." "(i) after estimates restore"
            mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
            tm_cond `=(`nv' == 0 & e(converged) == 0)' "(i) after estimates restore: e(V) or e(converged) not as stored"
            capture estimates drop tmnc
            local first 0
        }
    }

    * (iv) the same data, the original spec, metric(unitvariance): converges, no qualification, standard errors present
    tm_data 500 .70 .70 -.70 -.70 `seed'
    tm_run "`LG'" dwlsmv `mB', metric(unitvariance)
    local rc = r(rc)
    tm_cond `=`rc' == 0' "(iv) seed `seed': the unitvariance fit gave rc `rc', wanted 0"
    if !`rc' {
        tm_cond `=(e(converged) == 1)' "(iv) seed `seed': the unitvariance fit did not converge"
        mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
        tm_cond `=(`nv' == 15)' "(iv) seed `seed': e(V) has `nv' rows, wanted 15"
        tm_cond `=("`e(T3)'" != "" & "`e(p_T3)'" != "")' "(iv) seed `seed': e(T3) or e(p_T3) not posted"
        tm_expect "`LG'" hasnt "convergence not achieved" "(iv) seed `seed'"
        tm_expect "`LG'" hasnt "last iterate" "(iv) seed `seed'"
        tm_expect "`LG'" hasnt "not a solution" "(iv) seed `seed'"
        tm_expect "`LG'" hasnt "No recovery refit" "(iv) seed `seed'"
        tm_expect "`LG'" hasnt "The fitted solution is not admissible" "(iv) seed `seed'"
        tm_expect "`LG'" has "Std. err." "(iv) seed `seed'"
        tm_expect "`LG'" has "Prob > chi2" "(iv) seed `seed'"
    }
}
display as text _n "  (i): the nonconverged display was exercised on `nexi' of 4 seeds (`nskip' skipped: the fixed-gamma fit converged or gave an error); the admissibility note appeared on `nnote'"
tm_cond `=(`nexi' >= 1)' "(i): the fixed-gamma fit converged on all four seeds, so the nonconverged display, replay and restore were NOT exercised (choose another spec, or the hand-set route of test_condx_recovery.do C8)"

* ---- (i-b) configuration A seed 3: the first fit does not converge and is admissible (probe_inadm.log: psi 1.6e-5, Theta_j .36
* and above). Since 2026-11-15 the nonconvergence alone triggers D11's recovery: recovered (a converged fit, standard errors shown) or
* failed (the first fit is reported, qualified all the same, whatever the flag says) ----
display as text _n "  (i-b): configuration A, seed 3 (nonconverged first fit, admissible)"
tm_data 150 .20 .80 .80 .80 20261417
tm_run "`LG'" dwlsmv `mB'
tm_cond `=(r(rc) == 0)' "(i-b): rc `=r(rc)'"
if !r(rc) {
    local rec "`e(recovery)'"
    display as text "      e(recovery) = `rec', e(converged) = `=e(converged)'"
    tm_cond `=inlist("`rec'", "recovered", "failed")' "(i-b): e(recovery) is `rec', wanted recovered or failed (none: the premise, a nonconverged first fit, failed)"
    if "`rec'" == "recovered" {
        mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
        tm_cond `=(e(converged) == 1 & `nv' == 15 & "`e(T3)'" != "")' "(i-b): recovered, but e(converged), e(V) or e(T3) not as for a converged fit"
        tm_expect "`LG'" has "The marker-metric fit did not converge" "(i-b)"
        tm_expect "`LG'" hasnt "convergence not achieved" "(i-b)"
        tm_expect "`LG'" has "Std. err." "(i-b)"
    }
    else if "`rec'" == "failed" {
        tm_cond `=(e(converged) == 0)' "(i-b): recovery failed but e(converged) is `=e(converged)'"
        mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
        tm_cond `=(`nv' == 0 & "`e(T3)'" == "")' "(i-b): e(V) or e(T3) posted"
        mata: st_local("nx", strofreal(tm_next(st_local("LG"), st_local("title"), "convergence not achieved: the optimizer stopped after")))
        tm_cond `nx' "(i-b): the line under the title line is not the convergence line"
        tm_expect "`LG'" has "A recovery refit was attempted" "(i-b)"
        tm_expect "`LG'" has "and failed (see the notes)." "(i-b)"
        tm_expect "`LG'" hasnt "Std. err." "(i-b)"
        tm_expect "`LG'" hasnt "The fitted solution is not admissible" "(i-b): the qualification must not depend on a flag"
    }
}

* ---- (ii) a converged, admissible conditional marker fit ----
display as text _n "  (ii): converged, admissible"
tm_data2
tm_run "`LG'" dwlsmv `mA'
local rc = r(rc)
tm_cond `=(`rc' == 0)' "(ii): rc `rc'"
if !`rc' {
    tm_cond `=(e(converged) == 1)' "(ii): the fit did not converge"
    mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
    tm_cond `=(`nv' == 14)' "(ii): e(V) has `nv' rows, wanted 14"
    tm_cond `=(_se[y2:F] > 0 & _se[y2:F] < .)' "(ii): no standard error for y2:F"
    tm_cond `=("`e(T3)'" != "" & "`e(p_T3)'" != "" & e(df_T3) == 5)' "(ii): T3, p or df not as in test_condx_output.do"
    tm_expect "`LG'" hasnt "The fitted solution is not admissible" "(ii)"
    tm_expect "`LG'" hasnt "No recovery refit" "(ii)"
    tm_expect "`LG'" hasnt "convergence not achieved" "(ii)"
    tm_expect "`LG'" hasnt "last iterate" "(ii)"
    tm_expect "`LG'" has "Std. err." "(ii)"
    tm_expect "`LG'" has "Prob > chi2" "(ii)"
    tm_expect "`LG'" has "Not computed for a model with a covariate in this version: the check for thresholds" "(ii)"
    tm_expect "`LG'" has "that move together on every observation." "(ii)"
    tm_expect "`LG'" hasnt "the admissibility" "(ii): the sentence must name only what is still not computed"
    tm_expect "`LG'" hasnt "admissibility" "(ii)"
}

* ---- (iii-b) a converged INADMISSIBLE conditional fit through the command: searched for under metric(unitvariance) ----
* a Heywood case: a loading above 1 in absolute value with the conditional factor variance fixed at 1 (Theta_j < 0)
display as text _n "  (iii-b): searching for a converged inadmissible fit (N 100, loadings .95)"
local hw 0
forvalues s = 1/60 {
    if `hw' continue
    tm_data 100 .95 .95 .95 .95 `=20263000 + `s''
    capture quietly dwlsmv `mB', metric(unitvariance)
    if _rc continue
    if e(converged) != 1 continue
    local bad 0
    forvalues j = 1/4 {
        if abs(_b[y`j':F]) > 1 local bad 1
    }
    if `bad' local hw `s'
}
if `hw' {
    display as text "      found at seed " `=20263000 + `hw''
    tm_data 100 .95 .95 .95 .95 `=20263000 + `hw''
    tm_run "`LG'" dwlsmv `mB', metric(unitvariance)
    tm_cond `=(r(rc) == 0 & e(converged) == 1)' "(iii-b): the search's fit did not reproduce"
    tm_expect "`LG'" has "The fitted solution is not admissible" "(iii-b)"
    tm_expect "`LG'" has "residual variance" "(iii-b)"
    tm_expect "`LG'" has "Heywood" "(iii-b)"
    tm_expect "`LG'" has "Standard errors and the WLSMV test" "(iii-b)"
    tm_expect "`LG'" has "Std. err." "(iii-b): a converged fit keeps its standard errors"
    tm_expect "`LG'" hasnt "convergence not achieved" "(iii-b)"
    tm_expect "`LG'" hasnt "No recovery refit" "(iii-b): no such sentence under unitvariance"
    tm_expect "`LG'" hasnt "the admissibility" "(iii-b)"
    mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
    tm_cond `=(`nv' == 15 & "`e(T3)'" != "")' "(iii-b): e(V) or e(T3) missing on a converged fit"
    * the same data under the marker metric, whatever it does: the admissibility note is printed only for a fit that
    * triggered D11's recovery (e(recovery) not none: recovered, failed or ineligible), and the removed sentence is not printed
    tm_run "`LG'" dwlsmv `mB'
    if !r(rc) {
        mata: st_local("hasnote", strofreal(tm_has(st_local("LG"), "The fitted solution is not admissible")))
        if `hasnote' {
            tm_cond `=("`e(recovery)'" != "none" & "`e(recovery)'" != "")' "(iii-b): under marker the admissibility note appears but e(recovery) is `e(recovery)'"
        }
        tm_expect "`LG'" hasnt "No recovery refit was attempted for a model with a covariate" "(iii-b): marker"
    }
}
else {
    display as text "      none found in 60 seeds: (iii-b) is NOT exercised at the command level (Part A still covers the note)"
}

* ---- (v) an ordinary fit with e(converged) set to 0 by hand: since 2026-11-18 (owner: non-convergence is an
* inference boundary for every fit) it is displayed as a conditional one is, without the covariate's lines.
* Display only: e(V) is still in e() here (set by hand); the posting is tested on a real fit in test_unconv_post.do
display as text _n "  (v): ordinary fit, e(converged) set to 0 by hand"
tm_data2
quietly dwlsmv (F -> y1 y2 y3 y4)
tm_cond `=(e(converged) == 1)' "(v): the ordinary fit did not converge before the hand-set"
tm_unconv
tm_run "`LG'" dwlsmv
tm_cond `=(r(rc) == 0)' "(v): replay gave rc `=r(rc)'"
tm_expect "`LG'" has "convergence not achieved: the optimizer stopped" "(v)"
tm_expect "`LG'" has "last iterate, not a solution" "(v)"
tm_expect "`LG'" has "not reported (not a solution)" "(v)"
tm_expect "`LG'" has "Refit with metric(unitvariance)." "(v): the marker remedy"
tm_expect "`LG'" hasnt "Std. err." "(v): an ordinary nonconverged fit shows the last-iterate table"
tm_expect "`LG'" hasnt "Prob > chi2" "(v)"
tm_expect "`LG'" hasnt "RMSEA (scaled)" "(v)"
tm_expect "`LG'" hasnt "No further remedy is built for a model with a covariate" "(v)"
tm_expect "`LG'" hasnt "No recovery refit" "(v)"
* a real ordinary fit: nothing of the conditional text is printed
quietly dwlsmv (F -> y1 y2 y3 y4)
tm_run "`LG'" dwlsmv (F -> y1 y2 y3 y4)
tm_expect "`LG'" hasnt "Not computed for a model with a covariate" "(v)"
tm_expect "`LG'" hasnt "No recovery refit" "(v)"

adopath - "`here'/src/ado"
adopath - "`here'/build"

if ${TM_NFAIL} {
    display as error _n "test_condx_marker: ${TM_NFAIL} failure(s)"
    global TM_NFAIL
    exit 9
}
global TM_NFAIL
local ex = cond(`hw' > 0, "a converged inadmissible fit is reported with its note", "the converged inadmissible fit through the command was NOT found by the search, so the command-level display is not certified for it (Part A2 certifies the driver's A and the note)")
display as result "test_condx_marker: the trap-design marker fits are recovered (converged); a nonconverged conditional marker fit that recovery cannot repair (fixed gamma, `nexi' seed(s) exercised) is qualified directly under the title line, on replay and after estimates restore, posts no e(V), e(T3) or e(p_T3), prints no standard errors and says why no recovery refit was attempted; a converged admissible fit prints no note; the same data under unitvariance converge without qualification; an ordinary fit is displayed as before; `ex'"
exit 0
