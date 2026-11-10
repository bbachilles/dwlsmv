*! Drafted 2026-09-28 by a read-only agent. Executed by the main session on
*! StataMP 18, 2026-09-28 (nrep 5); the only change is one header label
*! (nr_noise is printed on C2O, not C2Q). Results and verdicts on the
*! predictions below: E-008, 2026-09-28, Module 2 slow path.
*! 2026-09-29: the BOUND group printed no C2X because rt/0.999 = 1 made
*! atanh() missing and the replay aborted (rc 3498). Added: a BOUND replay
*! (dc2x_bound) before that point, and a z0-missing guard. Executed by the
*! main session on StataMP 18 on the mc_c_999 group only (nrep 5), with the
*! condition loop narrowed in a %TEMP% copy; the full loop is not re-run.
*!
*! dwlsmv -- second diagnostic: what separates near-bound pairs on which
*! dwls_polychoric()'s d0 optimizer exhausts maxiter (converged 0, rc 0) from
*! pairs at similar rho on which it converges?  Measurement only: no
*! pass/fail, and the suite runner ignores tools/.
*!
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/diag_polychoric_conv2.do [nrep]
*!     grep "^C2" diag_polychoric_conv2.log      (wrapped lines continue with "> ")
*!
*! nrep = replicates per condition (default 5).
*!
*! KNOWN FROM THE FIRST DIAGNOSTIC (tools/diag_polychoric_conv.do, Stata 18):
*! 21 pairs inside the bound on reps 1-5 of mc_c_998, mc_B_995, mc_B_9975 hit
*! maxiter with rc 0; the polish is skipped; returned rho is 8.6e-10..5.2e-8
*! from the analytic score's unique root; drho/dz does NOT separate them from
*! converging pairs; they end with a nonzero numerical z-gradient. mc_c_999's
*! 4 bound pairs fail with rc 5 (missing values in numerical derivatives).
*! 7 healthy pair (3,4) cases show a second coarse-grid score sign change only
*! near rho = -0.999, where |S| ~ 1e-20.
*!
*! ---- PREDICTIONS, written before any run -------------------------------------
*! Groups: FAIL = converged 0 with |rho| < 0.9989 (expected 21); CTRL =
*! converged 1 with rho > 0.98, same reps of mc_c_998, mc_c_999, mc_B_995,
*! mc_B_9975; BOUND = |rho| >= 0.9989 (expected the 4 mc_c_999 pairs).
*!
*! (a) NOISY-LL HYPOTHESIS. FAIL pairs' ll(z) is non-smooth at the d0 step
*!     scale because some OBSERVED cell has a tiny or inaccurately computed
*!     probability (binormal() absolute floor, rectangle cancellation).
*!     Predicts: (1) FAIL tables have a small nonzero count (1-3) in an
*!     off-diagonal cell; (2) at the root that observed cell's pi is tiny
*!     (< ~1e-10) and/or its relative disagreement with the reflected route
*!     (C2P reldiff_obs) is >> 1e-12, while CTRL observed cells agree to
*!     ~1e-13 and have min pi_obs >> 1e-10; (3) FAIL: the curvature estimate
*!     D2(h)/h^2 departs from its large-h value at h where CTRL's is still
*!     flat, and the quadratic-fit residual rms (C2Q) is >> eps*|ll| (~2e-12);
*!     (4) CTRL: residual rms ~ eps*|ll|, D2/h^2 stable down to h ~ 1e-6.
*!     Refuted if FAIL and CTRL have comparable min pi_obs and reldiff_obs and
*!     comparable residual rms / curvature stability.
*! (b) EMPTY-CELL-PATTERN ALTERNATIVE. FAIL pairs differ from CTRL only in the
*!     empty-cell pattern (e.g. an empty interior/edge cell rather than only
*!     corners), with ll smooth. Predicts: (1) the n_empty / location field
*!     separates the groups cleanly; (2) min pi_obs and reldiff_obs similar
*!     across groups; (3)/(4) residual rms ~ eps*|ll| in BOTH groups.
*!     Refuted if patterns overlap between FAIL and CTRL.
*! (c) OPTIMIZER-LOGIC ALTERNATIVE. ll is smooth, but with |ll| ~ 1.1e4 the
*!     vtol 1e-13 relative criterion demands |dll| < ~1.1e-9 per iteration,
*!     and ptol 1e-13 relative in z (|z| ~ 3.5-4.3 here) demands |dz| <
*!     ~4e-13, below what a numerical d0 derivative with rounding noise
*!     eps*|ll| can resolve; nrtol 1e-11 on g H^-1 g' is then the binding
*!     criterion, and whether it is met depends on the gradient noise
*!     eps*|ll|/h_d0 relative to the z-curvature |H_z| = N I (drho/dz)^2.
*!     Predicts: (1),(2) no separation; (3),(4) residual rms ~ eps*|ll| in
*!     both groups; the separator is the scalar noise-to-curvature ratio
*!     R = (eps*|ll|/h_d0)^2 / |H_z| (printed as C2O nr_noise) -- FAIL pairs
*!     have R above ~1e-11 or the final nrtol quantity (C2O nr_final) sits
*!     above 1e-11, CTRL below; the instrumented replay (C2O) shows FAIL
*!     iterates wandering in z within a band where ll differences are at the
*!     rounding level, never a monotone approach. Refuted if nr_final and
*!     nr_noise overlap between FAIL and CTRL.
*! (d) PAIR (3,4) NEAR rho = -0.999. Prediction: rounding noise, not a root.
*!     The terms n*dpi/pi summed in S are individually ~1e-6..1 while S ~
*!     1e-20, i.e. |S| is BELOW eps * sum|terms| (C2M ratio < ~1), and the
*!     sign of S flips irregularly between adjacent grid points. A genuine
*!     root would show ratio >> 1 on both sides, a smooth monotone S through
*!     zero, and persistence under refinement.
*! ---------------------------------------------------------------------------
*!
*! d0 STEP SIZE. Not taken from documentation (Mata's optimize() manual says
*! d0 uses numerical derivatives with an adaptively chosen step, without a
*! formula this file relies on). It is MEASURED: an instrumented COPY of the
*! evaluator (dc2_ll_rec, identical arithmetic, plus logging) is handed to a
*! replay of the optimizer with production settings, and the distinct
*! |z_call - z_iterate| offsets are summarised (C2O hstep).
*!
*! INDEPENDENT ROUTE for cell probabilities: the reflection (X,Y) -> (-X,-Y),
*! i.e. dwlsmv_cellprob_raw on reversed negated thresholds, flipped back. Same
*! binormal() but different corner values, so its rounding is independent.
*! It is not an exact-arithmetic oracle; disagreement bounds the error.
*!
*! Line formats (prefix C2 + letter, fields separated by |):
*!   C2G  group assignment and production result
*!   C2T  table (rows separated by /), n_empty, min nonzero count, empty cells
*!        (c = off-diagonal corner, e = edge, i = interior)
*!   C2P  at the root: min pi observed/empty, count pi<=0, reflected-route
*!        relative differences, corner magnitudes of the min-pi observed cell
*!   C2H  one per offset h: ll(z+h)-ll(z), ll(z-h)-ll(z), D2, D2/h^2, central g
*!   C2Q  dense quadratic fit residuals, eps*|ll|, analytic H_z, nr_noise
*!   C2O  instrumented replay: calls, iterations, conv, nr_final, step sizes,
*!        spread of the last 40 calls
*!   C2X  missing-ll events (rc 5 pairs): z, rho, offending cell and value
*!   C2M  score near rho = -0.999 for pairs with a sign change there
*!   C2W  wall time per pair
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

display "C2|machine|stata `c(stata_version)' `c(flavor)'|processors `c(processors)'|`c(current_date)' `c(current_time)'"

mata:

real scalar dc_score(real matrix n, real colvector a, real colvector b,
                     real scalar rho)
{
    real matrix P, D, m, Pp

    P  = dwlsmv_cellprob_raw(a, b, rho)
    m  = (n :> 0) :* (P :> 0)
    Pp = P :* m :+ (1 :- m)
    D  = dwlsmv_dcellprob_rho(a, b, rho) :* m
    return(sum(n :* (D :/ Pp)) / sum(n))
}

real scalar dc_info(real colvector a, real colvector b, real scalar rho)
{
    real matrix P, D, m, Pp

    P  = dwlsmv_cellprob_raw(a, b, rho)
    m  = (P :> 0)
    Pp = P :* m :+ (1 :- m)
    D  = dwlsmv_dcellprob_rho(a, b, rho) :* m
    return(sum((D :* D) :/ Pp))
}

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

// Reflected route: P(X in cell) computed from the (-X,-Y) corners.
real matrix dc2_cellprob_refl(real colvector a, real colvector b, real scalar rho)
{
    real matrix Q
    real colvector ar, br

    ar = -a[rows(a)::1]
    br = -b[rows(b)::1]
    Q  = dwlsmv_cellprob_raw(ar, br, rho)
    return(Q[rows(Q)::1, cols(Q)::1])
}

// Production evaluator, called as optimize() would call it.
real scalar dc2_ll(real matrix n, real colvector a, real colvector b, real scalar z)
{
    real scalar lnf
    real rowvector g
    real matrix H

    lnf = .
    g   = J(1, 0, .)
    H   = J(0, 0, .)
    dwlsmv_polychoric_ll(0, (z), n, a, b, lnf, g, H)
    return(lnf)
}

// Instrumented COPY of the evaluator: same arithmetic, plus a call log and
// a record of the offending cell whenever it returns missing.
void dc2_ll_rec(real scalar todo, real rowvector z,
                real matrix n, real colvector a, real colvector b,
                real scalar lnf, real rowvector g, real matrix H)
{
    external real matrix DC2_LOG
    external real matrix DC2_MISS
    real matrix P, Pr, W
    real scalar r, s, rmin, smin, vmin, rho

    rho = 0.999 * tanh(z[1])
    Pr  = dwlsmv_cellprob_raw(a, b, rho)
    P   = Pr :* (n :> 0) :+ (n :== 0)
    if (min(P) <= 0) {
        lnf  = .
        vmin = .
        rmin = .
        smin = .
        for (r = 1; r <= rows(P); r++) {
            for (s = 1; s <= cols(P); s++) {
                if (n[r, s] > 0 & (missing(vmin) | Pr[r, s] < vmin)) {
                    vmin = Pr[r, s]
                    rmin = r
                    smin = s
                }
            }
        }
        DC2_MISS = DC2_MISS \ (z[1], rho, rmin, smin, vmin, n[rmin, smin])
    }
    else {
        lnf = sum(n :* ln(P))
    }
    DC2_LOG = DC2_LOG \ (z[1], lnf)
}

string scalar dc2_fmt(real matrix n)
{
    string scalar s
    real scalar r, c

    s = ""
    for (r = 1; r <= rows(n); r++) {
        if (r > 1) {
            s = s + "/"
        }
        for (c = 1; c <= cols(n); c++) {
            s = s + (c > 1 ? "," : "") + strofreal(n[r, c])
        }
    }
    return(s)
}

void dc2_pair(string scalar cond, string scalar tag, real scalar ia, real scalar ib)
{
    external real matrix DC_Y
    external pointer(real colvector) colvector DC_TAU
    external real matrix DC2_LOG
    external real matrix DC2_MISS
    struct dwlsmv_pairfit scalar f
    transmorphic S
    real matrix n, P, Pr, D, Q
    real colvector a, b, hs, llp, llm, gopt, zz, lz, res, ok
    real matrix X, Hopt
    real scalar K1, K2, N, t, t2, rho, z0, l0, k, h, d2, lo, hi, s1, s2, x
    real scalar nemp, minnz, r, c, pobs, pemp, nneg, rdo, rde, ro, co, rt
    real scalar eps, hz, dz, nrn, rc, oc, it, nrf, i0, zr0, rr, nsc, sabs
    real scalar mi, mj, vi, vj, r0, zs, hmin, hmed, hmax, zlo, zhi, llo, lhi
    real scalar inside4, isctrl, isfail, isbnd, j, m, cnr, rtfound
    string scalar grp, emp, pr

    timer_clear(94)
    timer_on(94)
    a  = *DC_TAU[ia]
    b  = *DC_TAU[ib]
    K1 = rows(a) + 1
    K2 = rows(b) + 1
    N  = rows(DC_Y)
    n  = dwlsmv_crosstab(DC_Y[., ia], DC_Y[., ib], K1, K2)
    pr = sprintf("%s|pair %g %g", tag, ia, ib)

    timer_clear(93)
    timer_on(93)
    f = dwls_polychoric(DC_Y[., ia], DC_Y[., ib], a, b)
    timer_off(93)
    t = timer_value(93)[1]

    inside4 = (cond == "mc_c_998" | cond == "mc_c_999" | cond == "mc_B_995" | cond == "mc_B_9975")
    isfail  = inside4 & (f.converged == 0) & (abs(f.rho) < 0.9989)
    isbnd   = inside4 & (abs(f.rho) >= 0.9989)
    isctrl  = inside4 & (f.converged == 1) & (f.rho > 0.98)
    grp = (isfail ? "FAIL" : (isbnd ? "BOUND" : (isctrl ? "CTRL" : "")))

    // (d) any pair: coarse score scan on [-0.999, -0.95]
    nsc  = 0
    s1   = dc_score(n, a, b, -0.999)
    for (k = 1; k <= 49; k++) {
        s2 = dc_score(n, a, b, -0.999 + 0.049 * k / 49)
        if (!missing(s1) & !missing(s2) & s1 != 0 & sign(s1) != sign(s2)) {
            nsc++
        }
        s1 = s2
    }
    if (nsc > 0) {
        printf("C2M|%s|signchanges_in[-0.999,-0.95] %g|group %s|rho_dp %12.9f\n", pr, nsc, (grp == "" ? "other" : grp), f.rho)
        for (k = 0; k <= 20; k++) {
            x  = -0.999 + 0.0005 * k
            P  = dwlsmv_cellprob_raw(a, b, x)
            D  = dwlsmv_dcellprob_rho(a, b, x)
            ok = vec((n :> 0) :* (P :> 0))
            Q  = vec(n) :* vec(D) :/ (vec(P) :* ok :+ (1 :- ok)) :* ok
            sabs = sum(abs(Q))
            s1 = sum(Q) / N
            pobs = min(vec(P) :/ ok)
            printf("C2M|%s|rho %10.7f|S/N %10.3e|sum|terms|/N %10.3e|ratio %10.3e|max|term|/N %10.3e|min_pi_obs %10.3e|nobs_pi<=0 %g\n",
                   pr, x, s1, sabs / N, abs(s1) / (epsilon(1) * sabs / N), max(abs(Q)) / N, pobs,
                   sum(vec(n :> 0) :* vec(P :<= 0)))
        }
    }

    if (grp == "") {
        timer_off(94)
        return
    }

    printf("C2G|%s|group %s|rho_dp %18.15f|converged %g|iter %g\n", pr, grp, f.rho, f.converged, f.iterations)

    // (1) table
    nemp  = sum(n :== 0)
    minnz = min(vec(n) :/ vec(n :> 0))
    emp   = ""
    for (r = 1; r <= K1; r++) {
        for (c = 1; c <= K2; c++) {
            if (n[r, c] == 0) {
                if ((r == 1 & c == K2) | (r == K1 & c == 1)) {
                    emp = emp + sprintf(" (%g,%g)c", r, c)
                }
                else if (r == 1 | r == K1 | c == 1 | c == K2) {
                    emp = emp + sprintf(" (%g,%g)e", r, c)
                }
                else {
                    emp = emp + sprintf(" (%g,%g)i", r, c)
                }
            }
        }
    }
    printf("C2T|%s|n %s|N %g|n_empty %g|min_nonzero %g|empty%s\n", pr, dwlsmv_fmt_safe(dc2_fmt(n)), N, nemp, minnz, emp)

    if (isbnd) {
        dc2x_bound(pr, n, a, b, ia, ib)
    }


    // root: bracket near f.rho, else f.rho itself
    rtfound = 0
    rt = f.rho
    if (!isbnd) {
        lo = max((f.rho - 0.005, -0.99899))
        s1 = dc_score(n, a, b, lo)
        for (k = 1; k <= 200; k++) {
            hi = min((f.rho - 0.005 + 0.01 * k / 200, 0.99899))
            s2 = dc_score(n, a, b, hi)
            if (!missing(s1) & !missing(s2) & s1 != 0 & sign(s1) != sign(s2)) {
                rt = dc_root(n, a, b, lo, hi)
                rtfound = 1
                break
            }
            lo = hi
            s1 = s2
        }
    }

    // (2) probabilities at the root
    P  = dwlsmv_cellprob_raw(a, b, rt)
    Pr = dc2_cellprob_refl(a, b, rt)
    pobs = .
    pemp = .
    nneg = 0
    rdo  = 0
    rde  = 0
    ro   = .
    co   = .
    for (r = 1; r <= K1; r++) {
        for (c = 1; c <= K2; c++) {
            if (P[r, c] <= 0) {
                nneg++
            }
            if (n[r, c] > 0) {
                if (missing(pobs) | P[r, c] < pobs) {
                    pobs = P[r, c]
                    ro = r
                    co = c
                }
                if (P[r, c] > 0) {
                    rdo = max((rdo, abs(P[r, c] - Pr[r, c]) / P[r, c]))
                }
            }
            else {
                if (missing(pemp) | P[r, c] < pemp) {
                    pemp = P[r, c]
                }
                if (P[r, c] > 0) {
                    rde = max((rde, abs(P[r, c] - Pr[r, c]) / P[r, c]))
                }
            }
        }
    }
    // corner magnitudes of the min-pi observed cell (Phi2 at its 4 corners)
    cnr = 0
    if (ro < K1 & co < K2) {
        cnr = binormal(a[ro], b[co], rt)
    }
    printf("C2P|%s|root %18.15f|rootfound %g|min_pi_obs %10.3e at (%g,%g) n=%g|refl %10.3e|UR_corner %10.3e\n",
           pr, rt, rtfound, pobs, ro, co, n[ro, co], Pr[ro, co], cnr)
    printf("C2P|%s|min_pi_emp %10.3e|n_pi<=0 %g|reldiff_obs %10.3e|reldiff_emp %10.3e\n",
           pr, pemp, nneg, rdo, rde)

    // (3)/(4) ll(z) around z_root
    z0  = atanh(rt / 0.999)

    if (missing(z0)) {
        printf("C2X|%s|z0 missing|rt/0.999 %21.18f|atanh undefined; C2H/C2Q/C2O skipped\n", pr, rt / 0.999)
        timer_off(94)
        displayflush()
        return
    }

    l0  = dc2_ll(n, a, b, z0)
    eps = epsilon(1) * abs(l0)
    dz  = 0.999 * (1 - tanh(z0)^2)
    hz  = -N * dc_info(a, b, rt) * dz^2
    hs  = J(13, 1, .)
    for (k = 1; k <= 13; k++) {
        hs[k] = 10^(-9 + 0.5 * (k - 1))
    }
    for (k = 1; k <= 13; k++) {
        h  = hs[k]
        s1 = dc2_ll(n, a, b, z0 + h) - l0
        s2 = dc2_ll(n, a, b, z0 - h) - l0
        d2 = s1 + s2
        printf("C2H|%s|h %8.2e|dll+ %11.4e|dll- %11.4e|D2 %11.4e|D2/h2 %11.4e|g_c %11.4e|D2/h2_over_Hz %10.3e\n",
               pr, h, s1, s2, d2, d2 / h^2, (s1 - s2) / (2 * h), (d2 / h^2) / hz)
    }
    // dense grid, quadratic least-squares fit, at three spacings
    for (j = 1; j <= 3; j++) {
        h  = (j == 1 ? 1e-7 : (j == 2 ? 1e-6 : 1e-5))
        zz = (-20::20) :* h
        lz = J(41, 1, .)
        for (k = 1; k <= 41; k++) {
            lz[k] = dc2_ll(n, a, b, z0 + zz[k]) - l0
        }
        X   = J(41, 1, 1), zz, zz :^ 2
        res = lz - X * (invsym(cross(X, X)) * cross(X, lz))
        printf("C2Q|%s|spacing %8.2e|ll0 %14.6f|eps*|ll| %10.3e|rms_resid %10.3e|max_resid %10.3e|rms/eps %9.2f|Hz_analytic %11.4e|drho/dz %10.3e\n",
               pr, h, l0, eps, sqrt(mean(res :^ 2)), max(abs(res)), sqrt(mean(res :^ 2)) / eps, hz, dz)
    }

    // instrumented replay of the optimizer, production settings
    DC2_LOG  = J(0, 2, .)
    DC2_MISS = J(0, 6, .)
    mi  = sum(DC_Y[., ia]) / N
    mj  = sum(DC_Y[., ib]) / N
    vi  = sum((DC_Y[., ia] :- mi) :^ 2) / N
    vj  = sum((DC_Y[., ib] :- mj) :^ 2) / N
    r0  = sum((DC_Y[., ia] :- mi) :* (DC_Y[., ib] :- mj)) / N / sqrt(vi * vj)
    if (r0 >  0.8) r0 =  0.8
    if (r0 < -0.8) r0 = -0.8
    zs = atanh(r0 / 0.999)
    S = optimize_init()
    optimize_init_evaluator(S, &dc2_ll_rec())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_which(S, "max")
    optimize_init_params(S, zs)
    optimize_init_argument(S, 1, n)
    optimize_init_argument(S, 2, a)
    optimize_init_argument(S, 3, b)
    optimize_init_conv_maxiter(S, 200)
    optimize_init_conv_ptol(S, 1e-13)
    optimize_init_conv_vtol(S, 1e-13)
    optimize_init_conv_nrtol(S, 1e-11)
    optimize_init_tracelevel(S, "none")
    rc   = _optimize(S)
    oc   = optimize_result_converged(S)
    it   = optimize_result_iterations(S)
    gopt = optimize_result_gradient(S)'
    Hopt = optimize_result_Hessian(S)
    nrf  = .
    if (rows(gopt) == 1 & rows(Hopt) == 1) {
        if (Hopt[1, 1] != 0 & !missing(Hopt[1, 1])) {
            nrf = abs(gopt[1]^2 / Hopt[1, 1])
        }
    }
    zr0 = optimize_result_params(S)[1]
    // offsets of calls from the final iterate (step sizes d0 used near the end)
    m   = rows(DC2_LOG)
    i0  = max((1, m - 39))
    Q   = abs(DC2_LOG[i0::m, 1] :- zr0)
    Q   = select(Q, Q :> 0)
    hmin = (rows(Q) ? min(Q) : .)
    hmed = (rows(Q) ? sort(Q, 1)[ceil(rows(Q) / 2)] : .)
    hmax = (rows(Q) ? max(Q) : .)
    zlo = min(DC2_LOG[i0::m, 1]) - z0
    zhi = max(DC2_LOG[i0::m, 1]) - z0
    llo = min(DC2_LOG[i0::m, 2]) - l0
    lhi = max(DC2_LOG[i0::m, 2]) - l0
    printf("C2O|%s|rc %g|conv %g|iter %g|calls %g|zfinal-zroot %10.3e|g_final %10.3e|H_final %11.4e|nr_final %10.3e|nmiss %g\n",
           pr, rc, oc, it, m, zr0 - z0, (rows(gopt) ? gopt[1] : .), (rows(Hopt) ? Hopt[1, 1] : .), nrf, rows(DC2_MISS))
    printf("C2O|%s|last40 hstep min %10.3e med %10.3e max %10.3e|z-zroot [%10.3e, %10.3e]|ll-ll0 [%10.3e, %10.3e]|nr_noise %10.3e\n",
           pr, hmin, hmed, hmax, zlo, zhi, llo, lhi, (hmed < . ? (eps / hmed)^2 / abs(hz) : .))
    for (k = 1; k <= min((rows(DC2_MISS), 8)); k++) {
        printf("C2X|%s|miss %g|z %12.9f|rho %14.12f|cell (%g,%g) n=%g|pi_raw %10.3e\n",
               pr, k, DC2_MISS[k, 1], DC2_MISS[k, 2], DC2_MISS[k, 3], DC2_MISS[k, 4], DC2_MISS[k, 6], DC2_MISS[k, 5])
    }
    if (rows(DC2_MISS)) {
        k = rows(DC2_MISS)
        rr = min(abs(DC2_MISS[., 2]))
        printf("C2X|%s|nmiss %g|smallest |rho| giving missing %14.12f\n", pr, k, rr)
    }
    timer_off(94)
    t2 = timer_value(94)[1]
    printf("C2W|%s|group %s|t_production %8.4f|t_diag %8.4f\n", pr, grp, t, t2)
    displayflush()
}

// guard against a table string that would exceed the line budget
string scalar dwlsmv_fmt_safe(string scalar s)
{
    if (strlen(s) > 150) {
        return(substr(s, 1, 150) + "...")
    }
    return(s)
}

// Threshold bounding a cell edge as text: k = 0 is -inf, k = rows(t)+1 is +inf.
string scalar dc2x_edge(real colvector t, real scalar k)
{
    if (k < 1) {
        return("-inf")
    }
    if (k > rows(t)) {
        return("+inf")
    }
    return(strofreal(t[k], "%12.9f"))
}

// Instrumented copy of dwlsmv_polychoric_ll: same arithmetic and same return
// value, plus one row in DC2X_EV for EVERY observed cell (n > 0) whose raw
// probability is <= 0 or missing, at every call where that happens.
// Columns: call, z, rho, i, j, n_ij, pi_raw, lnf_is_missing.
void dc2x_ll_rec(real scalar todo, real rowvector z,
                 real matrix n, real colvector a, real colvector b,
                 real scalar lnf, real rowvector g, real matrix H)
{
    external real matrix DC2X_EV
    external real scalar DC2X_NCALL
    external real matrix DC2X_LAST
    real matrix Pr, P
    real scalar r, s, rho, bad

    DC2X_NCALL = DC2X_NCALL + 1
    rho = 0.999 * tanh(z[1])
    Pr  = dwlsmv_cellprob_raw(a, b, rho)
    P   = Pr :* (n :> 0) :+ (n :== 0)
    bad = (min(P) <= 0)
    if (bad) {
        lnf = .
    }
    else {
        lnf = sum(n :* ln(P))
    }
    for (r = 1; r <= rows(Pr); r++) {
        for (s = 1; s <= cols(Pr); s++) {
            if (n[r, s] > 0) {
                if (missing(Pr[r, s]) | Pr[r, s] <= 0) {
                    DC2X_EV = DC2X_EV \ (DC2X_NCALL, z[1], rho, r, s, n[r, s], Pr[r, s], bad)
                }
            }
        }
    }
    DC2X_LAST = DC2X_LAST \ (DC2X_NCALL, z[1], rho, lnf)
    if (rows(DC2X_LAST) > 12) {
        DC2X_LAST = DC2X_LAST[(rows(DC2X_LAST) - 11)::rows(DC2X_LAST), .]
    }
}

// BOUND pairs: replay the production optimizer with the recording evaluator
// and report the observed cells found with pi <= 0. Needs no root and no
// atanh of the returned rho, so it cannot hit the z0 problem.
void dc2x_bound(string scalar pr, real matrix n, real colvector a,
                real colvector b, real scalar ia, real scalar ib)
{
    external real matrix DC_Y
    external real matrix DC2X_EV
    external real scalar DC2X_NCALL
    external real matrix DC2X_LAST
    transmorphic S
    real matrix U, E
    real scalar N, mi, mj, vi, vj, r0, zs, rc, k, i, j, u, nu, zf
    string scalar etxt

    DC2X_EV    = J(0, 8, .)
    DC2X_NCALL = 0
    DC2X_LAST  = J(0, 4, .)
    N  = rows(DC_Y)
    mi = sum(DC_Y[., ia]) / N
    mj = sum(DC_Y[., ib]) / N
    vi = sum((DC_Y[., ia] :- mi) :^ 2) / N
    vj = sum((DC_Y[., ib] :- mj) :^ 2) / N
    r0 = sum((DC_Y[., ia] :- mi) :* (DC_Y[., ib] :- mj)) / N / sqrt(vi * vj)
    if (r0 > 0.8) {
        r0 = 0.8
    }
    if (r0 < -0.8) {
        r0 = -0.8
    }
    zs = atanh(r0 / 0.999)

    S = optimize_init()
    optimize_init_evaluator(S, &dc2x_ll_rec())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_which(S, "max")
    optimize_init_params(S, zs)
    optimize_init_argument(S, 1, n)
    optimize_init_argument(S, 2, a)
    optimize_init_argument(S, 3, b)
    optimize_init_conv_maxiter(S, 200)
    optimize_init_conv_ptol(S, 1e-13)
    optimize_init_conv_vtol(S, 1e-13)
    optimize_init_conv_nrtol(S, 1e-11)
    optimize_init_tracelevel(S, "none")
    rc   = _optimize(S)
    zf   = optimize_result_params(S)[1]
    etxt = optimize_result_errortext(S)
    printf("C2X|%s|replay rc %g|errortext %s|iter %g|calls %g|z_final %14.9f|rho_final %16.13f|events %g\n",
           pr, rc, etxt, optimize_result_iterations(S), DC2X_NCALL, zf,
           0.999 * tanh(zf), rows(DC2X_EV))

    // last calls before the optimizer stopped
    for (k = 1; k <= rows(DC2X_LAST); k++) {
        printf("C2X|%s|lastcall %g|z %14.9f|rho %16.13f|lnf %s\n",
               pr, DC2X_LAST[k, 1], DC2X_LAST[k, 2], DC2X_LAST[k, 3],
               (missing(DC2X_LAST[k, 4]) ? "." : strofreal(DC2X_LAST[k, 4], "%18.9f")))
    }

    if (rows(DC2X_EV) == 0) {
        printf("C2X|%s|no observed cell with pi<=0 at any call\n", pr)
        displayflush()
        return
    }

    // first 8 events in full
    for (k = 1; k <= min((rows(DC2X_EV), 8)); k++) {
        i = DC2X_EV[k, 4]
        j = DC2X_EV[k, 5]
        printf("C2X|%s|event %g|call %g|z %14.9f|rho %16.13f|cell (%g,%g)|row z [%s, %s]|col z [%s, %s]|n %g|pi %10.3e|lnf_missing %g\n",
               pr, k, DC2X_EV[k, 1], DC2X_EV[k, 2], DC2X_EV[k, 3], i, j,
               dc2x_edge(a, i - 1), dc2x_edge(a, i), dc2x_edge(b, j - 1), dc2x_edge(b, j),
               DC2X_EV[k, 6], DC2X_EV[k, 7], DC2X_EV[k, 8])
    }

    // one summary line per distinct offending cell
    U  = uniqrows(DC2X_EV[., 4..5])
    nu = rows(U)
    for (u = 1; u <= nu; u++) {
        i = U[u, 1]
        j = U[u, 2]
        E = select(DC2X_EV, (DC2X_EV[., 4] :== i) :& (DC2X_EV[., 5] :== j))
        k = selectindex(abs(E[., 3]) :== min(abs(E[., 3])))[1]
        printf("C2X|%s|cell (%g,%g)|row z [%s, %s]|col z [%s, %s]|n %g|hits %g|smallest |rho| %16.13f|pi there %10.3e|largest |rho| %16.13f\n",
               pr, i, j, dc2x_edge(a, i - 1), dc2x_edge(a, i),
               dc2x_edge(b, j - 1), dc2x_edge(b, j), E[k, 6], rows(E),
               abs(E[k, 3]), E[k, 7], max(abs(E[., 3])))
    }
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
    printf("C2|%s|N %g|p %g|missing %g|distinct", tag, rows(DC_Y), cols(DC_Y), missing(DC_Y))
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
            display "C2|`c'|`r'|absent: `f'"
            continue
        }
        quietly use "`f'", clear
        capture noisily mata: dc_prep("`c'|`r'")
        if _rc {
            display "C2|`c'|`r'|prep aborted rc " _rc
            continue
        }
        local p = c(k)
        forvalues a = 1/`p' {
            forvalues b = `=`a'+1'/`p' {
                capture noisily mata: dc2_pair("`c'", "`c'|`r'", `a', `b')
                if _rc {
                    display "C2|`c'|`r'|pair `a' `b'|aborted rc " _rc
                }
            }
        }
    }
}
display "C2|done"
