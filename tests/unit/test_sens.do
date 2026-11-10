*! test_sens.do -- dwlsmv_sens(): the S1-S5 predictions of @sec-mg-sensitivity
*!
*! Executed under StataMP 18 on 2026-10-14; green at six processors.
*!
*! Mutants (scratch worktree, each diff printed before its run,
*! %TEMP%\sens\mutate.sh and mutate.log): 22 of 24 red.
*!   M1 J = P (s2, s3)                       M13 S_D without the 2 (s7)
*!   M2 J = P + R (s2, s3)                   M14 S_D without 1/(1 - c) (s7)
*!   M3 kappa_w without r (s4)               M15 S_D at h = 1e-4 (s7)
*!   M4 no 1 + |theta| (s2, s3, s6, s7)      M16 S_D without the epsV term (s7)
*!   M5 S_C without K (s6)                   M17 S_D on tau columns too (s7)
*!   M6 S_C floor without 1/(1 - c) (s6)     M18 epsV guard removed (e8)
*!   M7 column 4 with nu, not 2 nu (s3)      M19 iterations guard removed (e6)
*!   M8 c = max eig, not max |eig| (s3)      M20 convergence guard removed (e1)
*!   M9 c >= 1 refusal removed (e7)          M23 S_D group offset dropped (s7)
*!   M10 A without W (s2, s3)                M24 ds/eps length guard removed (e5)
*!   M11 column 5 without r (s2, s4)
*!   M12 nu without the 2 (s5)
*! Survivors, judged equivalent on this design:
*!   M21 R not symmetrized: cholesky() reads the lower triangle, and R's
*!       asymmetry is finite-difference noise.
*!   M22 R's step at 1e-2: sigma is a polynomial per coordinate on a single-
*!       factor design (Lemma MS2.1), so that difference has no truncation
*!       error at any step. A design with beta would separate it; none is
*!       tested here.
*!
*! One fit: design E (fixtures/data/mg_E.csv, G = 2) at the metric level,
*! built as test_mg_difftest.do builds it. E is consistency-only data for the
*! check 6 investigation (E-027, 2026-10-14). Nothing here
*! measures a processor-count residue, and no magnitude of any prediction is
*! asserted: every check is an identity or a derived property. epsV = 2u, the
*! derivation's value (mg-sens-d).
*!
*!   s1  premises: the fit converged in K >= 2 steps; the DWLS run, repeated
*!       through dwlsmv_minimize(), reproduces f.theta bitwise up to the sign
*!       canon; maxiter K - 2 returns theta_{K-1} (optimize() takes m + 1 steps
*!       at maxiter m), and one full Gauss-Newton step from there gives
*!       theta_K (S3(4)'s "full step, no halving")
*!   s2  S1 against central differences of tight refits (own Gauss-Newton, no
*!       halving), a known relative dw (eta 1e-4) and a known ds (1e-5);
*!       column 5 against the same differences in check 6's units. The
*!       refits carry S5's floor (~1e-12), hence the eta
*!   s3  A, c and column 4 against an independent J: the central difference
*!       of dwlsmv_fdwls()'s gradient (the Hessian of F is 2J); a control that
*!       J is not P here
*!   s4  an exact fit (s = sigma(theta-hat)): c, nu, kappa_w, the dw part of
*!       column 5 and column 6 are exactly 0 (Corollary S1.1, S5(b))
*!   s5  nu against dwlsmv_fdwls()'s g H^-1 g'
*!   s6  column 3: linear in |step| (at 2 step, step, 0), its step part
*!       proportional to K, and its floor (step = 0) against S3(4) evaluated
*!       here with the independent c
*!   s7  column 6: affine in epsV (at 2 epsV, epsV, 0); against S5(c)
*!       evaluated here, with the bound built from the parameter table and the
*!       independent c; S5's bound holds on production Delta: moving theta by
*!       one ulp (two patterns) changes Delta by at most twice |E|, plus a
*!       smooth allowance of 1e-13 off the bound's support
*!   e1-e8 refusals, each rc 3498: not converged; step short; step missing;
*!       ds without eps; ds short; iterations 0; c >= 1 (s = sigma + lambda r,
*!       lambda = 2 / mu, mu the eigenvalue of P^-1 R largest in modulus);
*!       epsV negative
*!
*! Tolerances, chosen before this version of the test was run:
*!   - 1e-6 on the identities of s2, s3 and s7, the value the derivation's
*!     chunks use (chosen after measuring there, mg-sens-a);
*!   - 1e-10 on s5;
*!   - 1e-8 of the step on s1's step;
*!   - 1e-12 relative on the linearity checks of s6 and s7;
*!   - 1e-13 for s7's smooth allowance.
*! eta 1e-4 and ds 1e-5 were chosen after a first run at 1e-5 and 1e-6 failed
*! against S5's floor (%TEMP%\sens\diag1.do: 4.7e-8 at eta 1e-4).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index
do tests/unit/sens_unpolished.mata

mata:

real colvector tsn_polish(real colvector s, real colvector w,
                          struct dwlsmv_parmap scalar par,
                          struct dwlsmv_statmap vector maps,
                          real colvector th0, real scalar n)
{
    real colvector th, r
    real matrix D
    real scalar i

    th = th0
    for (i = 1; i <= n; i++) {
        D  = dwls_delta(th, par, maps)
        r  = s - dwls_model(th, par, maps)
        th = th + cholsolve(quadcross(D, w, D), quadcross(D, w :* r))
    }
    return(th)
}

real scalar tsn_rel(real matrix x, real matrix y)
{
    return(max(abs(x - y)) / max(abs(y)))
}

// A fit struct with every group's s replaced by the matching slice of snew
struct dwlsmv_fit scalar tsn_with_s(struct dwlsmv_fit scalar f, real colvector snew)
{
    struct dwlsmv_fit scalar g
    real scalar k, o

    g = f
    o = 0
    for (k = 1; k <= length(g.stage1); k++) {
        g.stage1[k].s = snew[(o + 1)::(o + g.stage1[k].map.D)]
        o = o + g.stage1[k].map.D
    }
    return(g)
}

// The Hessian of F by central differences of dwlsmv_fdwls()'s gradient; J = H / 2
real matrix tsn_jfd(real colvector th, real colvector s, real colvector w,
                    struct dwlsmv_parmap scalar par, struct dwlsmv_statmap vector maps)
{
    real matrix Hf, Hd
    real rowvector gp, gm
    real scalar q, j, hj, v

    q  = rows(th)
    Hf = J(q, q, .)
    for (j = 1; j <= q; j++) {
        hj = 1e-5 * max((1, abs(th[j])))
        dwlsmv_fdwls(1, (th + hj * e(j, q)')', s, w, par, maps, v, gp, Hd)
        dwlsmv_fdwls(1, (th - hj * e(j, q)')', s, w, par, maps, v, gm, Hd)
        Hf[., j] = (gp - gm)' / (2 * hj)
    }
    return((Hf + Hf') / 4)
}

// S5's |E| bound, built here from the parameter table (an independent route)
real matrix tsn_ebound(real colvector th, real scalar epsV, struct dwlsmv_fit scalar f,
                       struct dwlsmv_statmap vector maps, real colvector sg, real matrix D)
{
    real matrix Eb
    real scalar i, j, g, o, hh, k, t

    Eb = J(rows(sg), rows(th), 0)
    for (i = 1; i <= f.par.nall; i++) {
        if (!f.par.free[i]) continue
        if (f.par.ptype[i] != 1 & f.par.ptype[i] != 2 & f.par.ptype[i] != 3) continue
        j  = f.par.fidx[i]
        g  = f.par.grp[i]
        o  = 0
        for (k = 1; k < g; k++) o = o + maps[k].D
        hh = 1e-5 * max((1, abs(th[j])))
        for (t = maps[g].T + 1; t <= maps[g].D; t++) {
            Eb[o + t, j] = Eb[o + t, j] + (epsV * abs(sg[o + t]) + (epsilon(1) / 2) * (abs(th[j]) + hh) * abs(D[o + t, j])) / hh
        }
    }
    return(Eb)
}

void tsn_setup()
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit, f2
    string scalar cc
    real matrix K, Y
    real colvector gv, s, w, start, thu, thK, thK1, ths
    real scalar G, g, N, i, fu, iu, fK, iK, fK1, iK1, cv
    pointer() scalar p

    K   = st_matrix("e(K)")
    cc  = st_global("e(canon_configural)")
    G   = rows(K)
    gv  = uniqrows(st_data(., "g"))
    sts = dwlsmv_stage1(G, 1)
    maps = dwlsmv_statmap(G, 1)
    for (g = 1; g <= G; g++) {
        Y       = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
        sts[g]  = dwls_sample_stats(Y, K[g, .]')
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(cc, K, "metric"), maps, sts)
    fit = dwlsmv_tsu_fit(sts, par)          // unpolished: the diagnostics model optimize()'s stop (sens_unpolished.mata)
    dwlsmv_fit_stack(sts, s, w, N, maps)

    // dwls_fit's sequence again: ULS from the table's starts, then DWLS
    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) start[par.fidx[i]] = par.value[i]
    }
    (void) dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, start, "d2", "nr", 300, thu, fu, iu)
    cv = dwlsmv_minimize(s, w, par, maps, thu, "d2", "nr", 300, thK, fK, iK)
    (void) dwlsmv_minimize(s, w, par, maps, thu, "d2", "nr", fit.iterations - 2, thK1, fK1, iK1)

    ths = tsn_polish(s, w, par, maps, fit.theta, 40)
    f2 = fit
    f2.theta = ths

    p = crexternal("TSN_FIT");  if (p == NULL) p = findexternal("TSN_FIT");  *p = fit
    p = crexternal("TSN_F2");   if (p == NULL) p = findexternal("TSN_F2");   *p = f2
    p = crexternal("TSN_THK");  if (p == NULL) p = findexternal("TSN_THK");  *p = thK
    p = crexternal("TSN_THK1"); if (p == NULL) p = findexternal("TSN_THK1"); *p = thK1
    p = crexternal("TSN_CONV"); if (p == NULL) p = findexternal("TSN_CONV"); *p = (cv, iK, iK1)
}

real scalar tsn_run()
{
    struct dwlsmv_fit scalar fit, f2, f3
    struct dwlsmv_statmap vector maps
    string colvector FAIL
    real colvector s, w, thK, thK1, step, th, r, xi, eps, ds, fw, fs, pw, ps, den, sg, ev, flo, env, t2
    real rowvector cvv, gp
    real matrix out, A, A3, D, P, Jf, Ji, Pi, L, Li, o2, o0, oK, Hd, Eb, D1, D2
    real scalar N, q, j, c, nu, c3, nu3, v, rr, K, cf, u, uq, gq, epsV, viol

    FAIL = J(0, 1, "")
    fit  = *findexternal("TSN_FIT")
    f2   = *findexternal("TSN_F2")
    thK  = *findexternal("TSN_THK")
    thK1 = *findexternal("TSN_THK1")
    cvv  = *findexternal("TSN_CONV")
    dwlsmv_fit_stack(fit.stage1, s, w, N, maps)
    q    = rows(fit.theta)
    K    = fit.iterations
    u    = epsilon(1) / 2
    epsV = 2 * u

    // ---- s1 ------------------------------------------------------------------
    step = thK - thK1
    D  = dwls_delta(thK1, fit.par, maps)
    r  = s - dwls_model(thK1, fit.par, maps)
    rr = max(abs(thK1 + cholsolve(quadcross(D, w, D), quadcross(D, w :* r)) - thK)) / max(abs(step))
    printf("{txt}  s1  K = %f; rerun converged %f in %f, maxiter K-2 gave %f steps; |theta_K| vs f.theta %g; last step %8.2e; full GN step from theta_{K-1}: %8.2e of it\n",
           K, cvv[1], cvv[2], cvv[3], max(abs(abs(thK) - abs(fit.theta))), max(abs(step)), rr)
    if (fit.converged != 1 | K < 2 | cvv[1] != 1 | cvv[2] != K | cvv[3] != K - 1) FAIL = FAIL \ "s1: premises (converged, K >= 2, step counts)"
    if (max(abs(abs(thK) - abs(fit.theta))) != 0) FAIL = FAIL \ "s1: the rerun does not reproduce f.theta up to sign"
    if (!(rr <= 1e-8)) FAIL = FAIL \ "s1: theta_K is not one full GN step from theta_{K-1}"

    // ---- s2 ------------------------------------------------------------------
    th  = f2.theta
    xi  = sign(sin(1::rows(s)))
    eps = 1e-4 :* xi
    ds  = 1e-5 :* xi
    out = dwlsmv_sens(f2, abs(step), epsV, A, c, nu, J(rows(s), 1, 0), eps)
    sg  = dwls_model(th, fit.par, maps)
    r   = s - sg
    den = 1 :+ abs(th)
    fw  = (tsn_polish(s, w :* (1 :+ eps), fit.par, maps, th, 40) - tsn_polish(s, w :* (1 :- eps), fit.par, maps, th, 40)) / 2
    fs  = (tsn_polish(s + ds, w, fit.par, maps, th, 40) - tsn_polish(s - ds, w, fit.par, maps, th, 40)) / 2
    pw  = A * (r :* eps)
    ps  = A * ds
    printf("{txt}  s2  dw: %8.2e   ds: %8.2e   column 5 vs |d theta|/(1+|theta|): %8.2e\n",
           tsn_rel(fw, pw), tsn_rel(fs, ps), tsn_rel(out[., 5], abs(fw) :/ den))
    if (!(tsn_rel(fw, pw) < 1e-6)) FAIL = FAIL \ "s2: S1 dw prediction misses the finite difference"
    if (!(tsn_rel(fs, ps) < 1e-6)) FAIL = FAIL \ "s2: S1 ds prediction misses the finite difference"
    if (!(tsn_rel(out[., 5], abs(fw) :/ den) < 1e-6)) FAIL = FAIL \ "s2: column 5 misses |d theta| / (1 + |theta|)"

    // ---- s3 ------------------------------------------------------------------
    Jf = tsn_jfd(th, s, w, fit.par, maps)
    D  = dwls_delta(th, fit.par, maps)
    P  = quadcross(D, w, D)
    L  = cholesky(P)
    Li = solvelower(L, I(q))
    ev = symeigenvalues(Li * (P - Jf) * Li')'
    cf = max(abs(ev))
    Ji = cholinv(Jf)
    printf("{txt}  s3  A: %8.2e   c %8.5f vs %8.5f   column 4: %8.2e   |J - P| / |P| = %8.2e\n",
           tsn_rel(A, cholsolve(Jf, (D :* w)')), c, cf,
           tsn_rel(out[., 4], sqrt(2 * nu * diagonal(Ji * P * Ji)) :/ den), tsn_rel(Jf, P))
    if (!(tsn_rel(A, cholsolve(Jf, (D :* w)')) < 1e-6)) FAIL = FAIL \ "s3: A misses the independent J"
    if (!(abs(c - cf) < 1e-6)) FAIL = FAIL \ "s3: c misses the independent J's"
    if (!(tsn_rel(out[., 4], sqrt(2 * nu * diagonal(Ji * P * Ji)) :/ den) < 1e-6)) FAIL = FAIL \ "s3: column 4 misses the independent J"
    if (!(tsn_rel(Jf, P) > 1e-3)) FAIL = FAIL \ "s3: control: J is indistinguishable from P on this fit"

    // ---- s4 ------------------------------------------------------------------
    f3 = tsn_with_s(f2, sg)
    o2 = dwlsmv_sens(f3, abs(step), epsV, A3, c3, nu3, J(rows(s), 1, 0), eps)
    printf("{txt}  s4  exact fit: c %g  nu %g  max kappa_w %g  max column 5 %g  max column 6 %g\n",
           c3, nu3, max(o2[., 2]), max(o2[., 5]), max(o2[., 6]))
    if (c3 != 0 | nu3 != 0 | max(o2[., 2]) != 0 | max(o2[., 5]) != 0 | max(o2[., 6]) != 0) {
        FAIL = FAIL \ "s4: an exact fit does not zero the dw channel and S5"
    }

    // ---- s5 ------------------------------------------------------------------
    dwlsmv_fdwls(2, th', s, w, fit.par, maps, v, gp, Hd)
    printf("{txt}  s5  nu %10.4e vs g H^-1 g' %10.4e\n", nu, gp * invsym(Hd) * gp')
    if (!(abs(nu / (gp * invsym(Hd) * gp') - 1) < 1e-10)) FAIL = FAIL \ "s5: nu misses g H^-1 g'"

    // ---- s6 ------------------------------------------------------------------
    o2 = dwlsmv_sens(f2, 2 * abs(step), epsV, A3, c3, nu3)
    o0 = dwlsmv_sens(f2, J(q, 1, 0), epsV, A3, c3, nu3)
    f3 = f2
    f3.iterations = 2 * K
    oK = dwlsmv_sens(f3, abs(step), epsV, A3, c3, nu3)
    uq = epsilon(1)^2
    gq = rows(s) * uq / (1 - rows(s) * uq)
    Pi = abs(cholinv(P))
    flo = (Pi * (2 * gq * quadcross(abs(D), abs(w :* r)))) / (1 - cf) :/ den
    printf("{txt}  s6  linearity %8.2e   K-proportionality %8.2e   floor %8.2e\n",
           tsn_rel(o2[., 3] - out[., 3], out[., 3] - o0[., 3]),
           tsn_rel(oK[., 3] - o0[., 3], 2 * (out[., 3] - o0[., 3])), tsn_rel(o0[., 3], flo))
    if (!(tsn_rel(o2[., 3] - out[., 3], out[., 3] - o0[., 3]) < 1e-12)) FAIL = FAIL \ "s6: S_C is not linear in |step|"
    if (!(tsn_rel(oK[., 3] - o0[., 3], 2 * (out[., 3] - o0[., 3])) < 1e-12)) FAIL = FAIL \ "s6: S_C's step part is not proportional to K"
    if (!(tsn_rel(o0[., 3], flo) < 1e-6)) FAIL = FAIL \ "s6: S_C's floor misses S3(4)"

    // ---- s7 ------------------------------------------------------------------
    o2 = dwlsmv_sens(f2, abs(step), 2 * epsV, A3, c3, nu3)
    o0 = dwlsmv_sens(f2, abs(step), 0, A3, c3, nu3)
    Eb  = tsn_ebound(th, epsV, f2, maps, sg, D)
    env = (2 * (Pi * quadcross(Eb, abs(w :* r))) / (1 - cf)) :/ den
    printf("{txt}  s7  affine in epsV %8.2e   column 6 vs S5(c) %8.2e   max column 6 %8.2e\n",
           tsn_rel(o2[., 6] - out[., 6], out[., 6] - o0[., 6]), tsn_rel(out[., 6], env), max(out[., 6]))
    if (!(tsn_rel(o2[., 6] - out[., 6], out[., 6] - o0[., 6]) < 1e-12)) FAIL = FAIL \ "s7: column 6 is not affine in epsV"
    if (!(tsn_rel(out[., 6], env) < 1e-6)) FAIL = FAIL \ "s7: column 6 misses S5(c)"
    // the bound on production Delta: two 1-ulp patterns
    viol = -1
    rr   = 0
    for (j = 1; j <= 2; j++) {
        t2 = th + (j == 1 ? 1 : -1) * sign(cos(1::q)) :* epsilon(1) :* abs(th)
        D1 = dwls_delta(t2, fit.par, maps)
        D2 = abs(D1 - D) - 2 * Eb - 1e-13 * (Eb :== 0)
        viol = max((viol, max(D2)))
        rr   = max((rr, max(abs(D1 - D) :/ (2 * Eb + (Eb :== 0)) :* (Eb :> 0))))
    }
    printf("{txt}  s7  1-ulp Delta change: at most %6.3f of twice the S5 bound on its support; largest excess %8.2e (must be <= 0)\n", rr, viol)
    if (viol > 0) FAIL = FAIL \ "s7: production Delta moves by more than S5's bound allows"

    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}  FAIL: %s\n", FAIL[j])
    }
    return(rows(FAIL))
}

// The refusals. which: 1 not converged, 2 step short, 3 step missing, 4 ds
// without eps, 5 ds short, 6 iterations 0, 7 c >= 1, 8 epsV negative.
void tsn_err(real scalar which)
{
    struct dwlsmv_fit scalar f, g
    struct dwlsmv_statmap vector maps
    real colvector s, w, th, r, step, sg, ev
    real matrix A, D, P, L, Li, out
    real scalar c, nu, N, q, mu, epsV

    f = *findexternal("TSN_F2")
    q = rows(f.theta)
    dwlsmv_fit_stack(f.stage1, s, w, N, maps)
    step = J(q, 1, 1e-6)
    epsV = epsilon(1)
    if (which == 1) {
        f.converged = 0
        out = dwlsmv_sens(f, step, epsV, A, c, nu)
    }
    else if (which == 2) out = dwlsmv_sens(f, J(q - 1, 1, 1e-6), epsV, A, c, nu)
    else if (which == 3) out = dwlsmv_sens(f, J(q, 1, .), epsV, A, c, nu)
    else if (which == 4) out = dwlsmv_sens(f, step, epsV, A, c, nu, J(rows(s), 1, 0))
    else if (which == 5) out = dwlsmv_sens(f, step, epsV, A, c, nu, J(rows(s) - 1, 1, 0), J(rows(s), 1, 0))
    else if (which == 6) {
        f.iterations = 0
        out = dwlsmv_sens(f, step, epsV, A, c, nu)
    }
    else if (which == 7) {
        th = f.theta
        sg = dwls_model(th, f.par, maps)
        r  = s - sg
        D  = dwls_delta(th, f.par, maps)
        P  = quadcross(D, w, D)
        L  = cholesky(P)
        Li = solvelower(L, I(q))
        ev = symeigenvalues(Li * (P - tsn_jfd(th, s, w, f.par, maps)) * Li')'
        mu = ev[order(-abs(ev), 1)[1]]
        g  = tsn_with_s(f, sg + (2 / mu) * r)
        printf("{txt}  e7  mu = %g, so lambda = %g\n", mu, 2 / mu)
        out = dwlsmv_sens(g, step, epsV, A, c, nu)
    }
    else if (which == 8) out = dwlsmv_sens(f, step, -1, A, c, nu)
}

end

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
mata: tsn_setup()

local nf 0
mata: st_local("nf", strofreal(tsn_run()))

forvalues k = 1/8 {
    capture noisily mata: tsn_err(`k')
    if _rc != 3498 {
        display as error "    e`k' FAILED: rc `=_rc', expected 3498"
        local ++nf
    }
    else {
        display as text "    e`k' refused, rc 3498"
    }
}

if `nf' {
    display as error "test_sens: `nf' failure(s)"
    exit 9
}
display as result "test_sens: dwlsmv_sens() reproduces S1 by finite differences, its J, c and stopping envelope against an independent Hessian, zeroes the dw channel and S5 at an exact fit, its S_C is linear in the step and proportional to K, and its S5 envelope matches the formula and bounds production Delta's 1-ulp movement"
