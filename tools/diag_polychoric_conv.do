*! Drafted 2026-09-28 by a read-only agent. Executed unchanged by the main
*! session on StataMP 18, 2026-09-28 (nrep 5). Results and verdicts on the
*! predictions below: E-008, 2026-09-28, Module 2 slow path.
*! Result lines wrap at linesize 255; join "> " continuation lines to read.
*!
*! dwlsmv -- why do near-bound pairs come back converged = 0 from
*! dwls_polychoric()?  A diagnostic, not a test: no pass/fail, and the suite
*! runner ignores tools/.
*!
*!     cd <repo root>   (or via a %TEMP% wrapper, see below)
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/diag_polychoric_conv.do [nrep]
*!     grep "^CONV|" diag_polychoric_conv.log
*!
*! nrep = replicates per condition (default 5).
*!
*! OBSERVATION. On replicates 1-5 of mc_c_998, mc_B_995, mc_B_9975, 21 pairs
*! with rho-hat in 0.988-0.9982 (inside the 0.999 bound) return converged = 0,
*! agree with Kolenikov's polychoric to ~1e-6, and are the slow path. Four
*! pairs of mc_c_999 sit at the bound.
*!
*! HYPOTHESIS UNDER TEST (not assumed). Near rho ~ .996, drho/dz =
*! 0.999 (1 - tanh(z)^2) ~ 0.006, so the d0 optimizer in z cannot meet
*! ptol/vtol 1e-13 / nrtol 1e-11, exhausts maxiter 200, reports not
*! converged, and the Fisher polish is skipped; the returned rho is near, not
*! at, the root of the analytic score.
*!
*! ---- PREDICTIONS, written before any run -------------------------------------
*! Stated per quantity: H = what the hypothesis predicts; R = what refutes it;
*! A = what an alternative predicts. A1 = genuinely flat likelihood /
*! near-non-identification; A2 = multiple roots; A3 = the feasibility guard
*! (polish infeasible); A4 = optimizer error (_optimize rc != 0, e.g. a
*! missing log-likelihood at a trial point).
*!
*! P1 rc / errortext. H: rc = 0, no error text, for all 21 inside-bound pairs
*!    (non-convergence is reported via optimize_result_converged, not rc).
*!    R: rc != 0 on those pairs. A4 predicts rc != 0 with a "missing values"
*!    or "could not calculate" style message.
*! P2 optimize_result_converged. H: 0 on the 21 pairs, 1 on controls.
*!    R: 1 on the flagged pairs (then the cause is bound or polish, see P4).
*! P3 iterations (optimizer only). H: = 200 (maxiter) on the 21 pairs;
*!    controls well below 200. R: flagged pairs stop below 200 with
*!    converged 0. The slowness observation also follows from H: 200 d0
*!    iterations, each several evaluations.
*! P4 cause. H: "optnotconv" for the 21; "bound" for the 4 mc_c_999 pairs.
*!    R: "rc" or "polish". A3 predicts "polish".
*! P5 |rho_returned - root|. H: > 1e-8 for the 21 (optimizer floor in rho is
*!    ~1e-8 when converged; unpolished and unconverged should be worse), but
*!    < ~1e-6 (Kolenikov agreement). Controls (polished) < 1e-12.
*!    R: |rho - root| < 1e-10 on the flagged pairs (then they are at the root
*!    and only the flag is wrong). A1 predicts a large |rho - root| with a
*!    tiny score (flat); A2 predicts a large distance to the root NEAREST
*!    the returned point is absent but distance to some other root is small.
*! P6 score S/N at rho_returned. H: small but clearly nonzero, of order
*!    |rho - root| * I (P9 consistent with P5). Controls ~1e-15.
*!    R: |S/N| ~ 1e-15 with converged 0.
*! P7 expected information I at rho_returned. H: large and finite (the
*!    likelihood in rho is well curved; the difficulty is the z-map).
*!    A1 predicts I small relative to controls.
*! P8 root / sign changes on the grid. H: exactly one sign change, root
*!    inside (-0.999, 0.999) for the 21; for mc_c_999 bound pairs, no root
*!    inside the bound, score positive at +0.999. R/A2: two or more sign
*!    changes on any flagged pair.
*! P9 one-step Fisher distance S/(N*I). H: agrees with rho_root - rho_ret in
*!    sign and to within a few percent in size (locally quadratic).
*!    A1 predicts poor agreement (curvature not captured).
*! P10 drho/dz at the returned point. H: ~0.002-0.02 on the 21 (rho
*!    0.988-0.9982 gives 0.024-0.0026); controls (mc_d, mc_c_98 lower-rho
*!    pairs) ~0.04 or more; the failing pairs cluster at the smallest values.
*!    R: failing and converging pairs overlap fully in drho/dz.
*! P11 wall time. H: flagged pairs markedly slower than controls (200
*!    iterations vs a few tens). R: similar times.
*!
*! Interpretation is kept out of this file: it measures and prints.
*!
*! WHAT IS MEASURED AGAINST. The score is production's estimating equation,
*! built from dwlsmv_cellprob_raw() and dwlsmv_dcellprob_rho():
*!     S(rho)/N = sum over OBSERVED cells (n > 0, pi > 0) of n * dpi / pi, / N
*! and I(rho) = sum over cells with pi > 0 of dpi^2 / pi (the polish's
*! scaler, the rk_rho rule). The root is found here by bisection inside a
*! sign-change bracket on the grid, then safeguarded Newton (Fisher) steps
*! kept inside the bracket until |step| < 1e-15 -- the rk_rho method of
*! derivation @sec-gamma-rank, used as method only. dwls_polychoric() is not
*! modified; its optimizer is REPLAYED here with identical settings to read
*! rc, error text, converged and iterations, and the replay's rho is printed
*! beside dwls_polychoric()'s so any divergence is visible.
*!
*! DESIGNS. Replicates 1..nrep of mc_c_998, mc_c_999, mc_B_995, mc_B_9975,
*! and healthy controls mc_d and mc_c_98 (fixtures/replication/, gitignored;
*! regenerate with fixtures/generate/50_replicates.R). K = (4,4,3,3) as in
*! tools/compare_polychoric.do; distinct values per variable are printed so a
*! mismatch is visible (assumed, not verified, for mc_d).
*!
*! Line formats (all start CONV|):
*!   CONV|cond|rep|pair a b|rc %g|err <text>|optconv %g|iter %g|polycon %g|
*!        cause <c>|rho_dp ..|rho_replay ..|S/N ..|I ..|root ..|
*!        absdiff ..|fisher1 ..|drho/dz ..|optgrad_z ..|time ..|signchg %g|
*!        S(-bound) ..|S(+bound) ..
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
args nrep
if "`nrep'" == "" local nrep 5
clear all
set more off
set linesize 255

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
quietly do src/mata/dwlsmv_types.mata
foreach f of local srcs {
    quietly do "src/mata/`f'"
}

display "CONV|machine|stata `c(stata_version)' `c(flavor)'|processors `c(processors)'|`c(current_date)' `c(current_time)'"

mata:

// Score / N over observed cells with pi > 0 (as the evaluator's sum runs).
real scalar dc_score(real matrix n, real colvector a, real colvector b,
                     real scalar rho)
{
    real matrix P, D, m, Pp
    real scalar N

    N  = sum(n)
    P  = dwlsmv_cellprob_raw(a, b, rho)
    m  = (n :> 0) :* (P :> 0)
    Pp = P :* m :+ (1 :- m)
    D  = dwlsmv_dcellprob_rho(a, b, rho) :* m
    return(sum(n :* (D :/ Pp)) / N)
}

// Expected information over cells with pi > 0 (the polish's scaler).
real scalar dc_info(real colvector a, real colvector b, real scalar rho)
{
    real matrix P, D, m, Pp

    P  = dwlsmv_cellprob_raw(a, b, rho)
    m  = (P :> 0)
    Pp = P :* m :+ (1 :- m)
    D  = dwlsmv_dcellprob_rho(a, b, rho) :* m
    return(sum((D :* D) :/ Pp))
}

// Root inside [lo, hi] (sign change assumed): bisection to width 1e-10,
// then Fisher steps kept inside the bracket until |step| < 1e-15.
real scalar dc_root(real matrix n, real colvector a, real colvector b,
                    real scalar lo0, real scalar hi0)
{
    real scalar lo, hi, slo, mid, sm, r, s, st, k

    lo  = lo0
    hi  = hi0
    slo = dc_score(n, a, b, lo)
    for (k = 1; k <= 200; k++) {
        if (hi - lo < 1e-10) {
            break
        }
        mid = (lo + hi) / 2
        sm  = dc_score(n, a, b, mid)
        if (sm == 0) {
            return(mid)
        }
        if (sign(sm) == sign(slo)) {
            lo  = mid
            slo = sm
        }
        else {
            hi = mid
        }
    }
    r = (lo + hi) / 2
    for (k = 1; k <= 50; k++) {
        s  = dc_score(n, a, b, r)
        st = s / dc_info(a, b, r)
        if (r + st <= lo0 | r + st >= hi0) {
            break
        }
        r = r + st
        if (abs(st) < 1e-15) {
            break
        }
    }
    return(r)
}

// One pair.
void dc_pair(string scalar tag, real scalar ia, real scalar ib)
{
    external real matrix DC_Y
    external pointer(real colvector) colvector DC_TAU
    struct dwlsmv_pairfit scalar f
    transmorphic S
    real matrix n
    real colvector a, b, g, sg
    real scalar K1, K2, N, t, rc, oc, it, rrep, mi, mj, vi, vj, r0, z0
    real scalar sn, inf, rt, fis1, dz, gz, nsc, k, best, bd, lo, hi, sprev
    real scalar slb, sub, nref, x, sx
    string scalar err, cause, rts

    a  = *DC_TAU[ia]
    b  = *DC_TAU[ib]
    K1 = rows(a) + 1
    K2 = rows(b) + 1
    N  = rows(DC_Y)
    n  = dwlsmv_crosstab(DC_Y[., ia], DC_Y[., ib], K1, K2)

    // production call, timed
    timer_clear(93)
    timer_on(93)
    f = dwls_polychoric(DC_Y[., ia], DC_Y[., ib], a, b)
    timer_off(93)
    t = timer_value(93)[1]

    // replay of the optimizer, settings copied from dwls_polychoric()
    mi  = sum(DC_Y[., ia]) / N
    mj  = sum(DC_Y[., ib]) / N
    vi  = sum((DC_Y[., ia] :- mi) :^ 2) / N
    vj  = sum((DC_Y[., ib] :- mj) :^ 2) / N
    r0  = sum((DC_Y[., ia] :- mi) :* (DC_Y[., ib] :- mj)) / N / sqrt(vi * vj)
    if (r0 >  0.8) r0 =  0.8
    if (r0 < -0.8) r0 = -0.8
    z0 = atanh(r0 / 0.999)
    S = optimize_init()
    optimize_init_evaluator(S, &dwlsmv_polychoric_ll())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_which(S, "max")
    optimize_init_params(S, z0)
    optimize_init_argument(S, 1, n)
    optimize_init_argument(S, 2, a)
    optimize_init_argument(S, 3, b)
    optimize_init_conv_maxiter(S, 200)
    optimize_init_conv_ptol(S, 1e-13)
    optimize_init_conv_vtol(S, 1e-13)
    optimize_init_conv_nrtol(S, 1e-11)
    optimize_init_tracelevel(S, "none")
    rc   = _optimize(S)
    err  = (rc != 0 ? optimize_result_errortext(S) : "none")
    oc   = optimize_result_converged(S)
    it   = optimize_result_iterations(S)
    rrep = 0.999 * tanh(optimize_result_params(S)[1])
    g    = optimize_result_gradient(S)
    gz   = (length(g) ? g[1] : .)

    if (rc != 0)                    cause = "rc"
    else if (!oc)                   cause = "optnotconv"
    else if (abs(rrep) >= 0.9989)   cause = "bound"
    else if (f.converged == 0)      cause = "polish"
    else                            cause = "none"

    sn   = dc_score(n, a, b, f.rho)
    inf  = dc_info(a, b, f.rho)
    fis1 = sn / inf
    dz   = 0.999 - f.rho * f.rho / 0.999

    // grid scan, 2000 points over [-0.999, 0.999]
    sg  = J(2000, 1, .)
    for (k = 1; k <= 2000; k++) {
        sg[k] = dc_score(n, a, b, -0.999 + 1.998 * (k - 1) / 1999)
    }
    nsc   = 0
    best  = .
    bd    = .
    sprev = .
    for (k = 2; k <= 2000; k++) {
        if (missing(sg[k - 1]) | missing(sg[k])) {
            continue
        }
        if (sign(sg[k - 1]) != sign(sg[k]) & sg[k - 1] != 0) {
            nsc++
            lo = -0.999 + 1.998 * (k - 2) / 1999
            hi = -0.999 + 1.998 * (k - 1) / 1999
            rt = dc_root(n, a, b, lo, hi)
            if (missing(bd) | abs(rt - f.rho) < bd) {
                bd   = abs(rt - f.rho)
                best = rt
            }
        }
    }
    slb = sg[1]
    sub = sg[2000]
    // refined grid: 2000 points within +-1e-3 of the chosen root (or of
    // rho_dp if none); counts extra sign changes a coarse grid could miss
    rt = (missing(best) ? f.rho : best)
    sprev = .
    nref = 0
    for (k = 1; k <= 2000; k++) {
        x = rt - 1e-3 + 2e-3 * (k - 1) / 1999
        if (x <= -0.999 | x >= 0.999) {
            continue
        }
        sx = dc_score(n, a, b, x)
        if (!missing(sprev) & !missing(sx)) {
            if (sign(sprev) != sign(sx) & sprev != 0) {
                nref++
            }
        }
        sprev = sx
    }
    rts = (missing(best) ? "no root inside bound" : sprintf("%18.15f", best))

    printf("CONV|%s|pair %g %g|rc %g|err %s|optconv %g|iter %g|polycon %g|polyiter %g|cause %s|rho_dp %18.15f|rho_replay %18.15f|S/N %10.3e|I %10.3e|root %s|absdiff %10.3e|fisher1 %10.3e|drho/dz %10.3e|optgrad_z %10.3e|time %8.4f|signchg %g|signchg_refined %g|S(-bound) %10.3e|S(+bound) %10.3e\n",
           tag, ia, ib, rc, err, oc, it, f.converged, f.iterations, cause,
           f.rho, rrep, sn, inf, rts, bd, fis1, dz, gz, t, nsc, nref, slb, sub)
    displayflush()
}

void dc_prep(string scalar tag)
{
    external real matrix DC_Y
    external pointer(real colvector) colvector DC_TAU
    real colvector K, t
    real scalar j

    DC_Y = st_data(., .)
    K    = (4 \ 4 \ 3 \ 3)
    printf("CONV|%s|N %g|p %g|missing %g|distinct", tag, rows(DC_Y), cols(DC_Y), missing(DC_Y))
    for (j = 1; j <= cols(DC_Y); j++) {
        printf(" %g", rows(uniqrows(DC_Y[., j])))
    }
    printf("\n")
    DC_TAU = J(cols(DC_Y), 1, NULL)
    for (j = 1; j <= cols(DC_Y); j++) {
        t = dwlsmv_thresholds(DC_Y[., j], K[j])
        DC_TAU[j] = &(t[.])
    }
    displayflush()
}

end

foreach c in mc_c_998 mc_c_999 mc_B_995 mc_B_9975 mc_d mc_c_98 {
    forvalues r = 1/`nrep' {
        local f : display "fixtures/replication/`c'/`c'_" %04.0f `r' ".dta"
        capture confirm file "`f'"
        if _rc {
            display "CONV|`c'|`r'|absent: `f'"
            continue
        }
        quietly use "`f'", clear
        capture noisily mata: dc_prep("`c'|`r'")
        if _rc {
            display "CONV|`c'|`r'|prep aborted rc " _rc
            continue
        }
        local p = c(k)
        forvalues a = 1/`p' {
            forvalues b = `=`a'+1'/`p' {
                capture noisily mata: dc_pair("`c'|`r'", `a', `b')
                if _rc {
                    display "CONV|`c'|`r'|pair `a' `b'|aborted rc " _rc
                }
            }
        }
    }
}
display "CONV|done"
