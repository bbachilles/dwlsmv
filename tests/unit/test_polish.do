*! Module 7: the post-convergence polish, dwlsmv_polish() (src/mata/dwls_fit.mata).
*!
*! Spec Module 7 (adopted 2026-10-23); @sec-polish; validation E-014.
*! Fixture mg_E (two groups), configural and scalar, built as the front end builds them.
*!
*!   1   wiring: dwls_fit()'s theta equals ULS -> dwlsmv_minimize() -> dwlsmv_polish() -> signcanon
*!       replicated by hand, bitwise; iterations and converged are dwlsmv_minimize()'s (P6, P7). ULS
*!       is replicated unpolished, so a polished ULS would break the equality (P2).
*!   2   fixed point: polishing dwls_fit()'s theta stops by the threshold after one step (status 0).
*!   3   contraction: from theta-hat + a planted offset of 1e-7 (relative, check 6's unit), the polish
*!       returns to within 1e-10 of theta-hat (relative). 1e-10 is 10x the validation's S1 threshold,
*!       chosen after measuring (E-014: max 9.7e-12 across counts under 16.1; 1.4e-11 under 18, polish-s18b).
*!   4   abandonment by size: from an offset of 1e-3 the first step exceeds 1e-6; status -2 and
*!       theta returned exactly as passed.
*!   5   abandonment by a missing step: w = 0 makes H zero, cholsolve() returns missing; status -1
*!       and theta returned exactly as passed.
*!   6   the first-step limit: a first step of about 3e-6 is abandoned (-2).
*!   7   abandonment by growth (-3), theta untouched; 8 the cap (status 2), theta kept, configural.
*!       Both replace s by s + m r, and m is chosen after measuring (2026-11-02).
*!   9   the loop against a hand-written one with the same evaluator and solver, bitwise, so it
*!       tests the stopping logic and is not an independent numerical check.
*!   10  G = 1: dwls_fit is polished, as in check 1. 11 a saturated G = 1 fit stops by the threshold.
*!   Not tested here: that no F comparison decides a step (a property of the code, read in review),
*!   and that an unconverged run is left unpolished (no reachable run isolates it yet).
version 16
clear all
set more off
local here `c(pwd)'
quietly do tools/build_mlib.do
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)

mata:
real scalar tp_rel(real colvector x, real colvector y)
{
    return(max(abs(x - y) :/ (1 :+ abs(y))))
}

real scalar tp_level(string scalar lev)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f
    real matrix K, Y
    real colvector gv, s, w, start, thu, th, t0, t1, off, r, d
    real rowvector pr, gg
    real matrix H
    real scalar G, g, N, i, fu, iu, fm, it, conv, nf, rel

    nf = 0
    K = st_matrix("e(K)")
    G = rows(K)
    gv = uniqrows(st_data(., "g"))
    sts = dwlsmv_stage1(G, 1)
    maps = dwlsmv_statmap(G, 1)
    for (g = 1; g <= G; g++) {
        Y = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
        sts[g] = dwls_sample_stats(Y, K[g, .]')
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(st_global("e(canon_configural)"), K, lev), maps, sts)
    f = dwls_fit(sts, par)
    dwlsmv_fit_stack(sts, s, w, N, maps)

    // 1 wiring
    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) start[par.fidx[i]] = par.value[i]
    }
    (void) dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, start, "d2", "nr", 300, thu, fu, iu)
    conv = dwlsmv_minimize(s, w, par, maps, thu, "d2", "nr", 300, th, fm, it)
    pr = dwlsmv_polish(s, w, par, maps, th)
    th = dwlsmv_signcanon(th, par)
    if (th == f.theta & it == f.iterations & conv == f.converged & pr[1] == 0) {
        printf("  1. %s: dwls_fit = ULS, minimize, polish (%f steps), signcanon, bitwise; iterations %f\n", lev, pr[2], it)
    }
    else {
        printf("{err}  1. %s FAILED: replication differs (status %f, iterations %f vs %f)\n", lev, pr[1], it, f.iterations)
        nf++
    }

    // 2 fixed point
    t0 = f.theta
    pr = dwlsmv_polish(s, w, par, maps, t0)
    if (pr[1] == 0 & pr[2] == 1) {
        printf("  2. %s: theta-hat is a fixed point (one step, %9.2e)\n", lev, pr[3])
    }
    else {
        printf("{err}  2. %s FAILED: status %f, steps %f, last %9.2e\n", lev, pr[1], pr[2], pr[3])
        nf++
    }

    // 3 contraction from a planted 1e-7 offset
    off = 1e-7 :* (1 :+ abs(f.theta)) :* (mod(1::rows(f.theta), 2) :* 2 :- 1)
    t1 = f.theta + off
    pr = dwlsmv_polish(s, w, par, maps, t1)
    if (pr[1] == 0 & tp_rel(t1, f.theta) <= 1e-10) {
        printf("  3. %s: 1e-7 offset contracts to %9.2e of theta-hat in %f steps\n", lev, tp_rel(t1, f.theta), pr[2])
    }
    else {
        printf("{err}  3. %s FAILED: status %f, distance %9.2e\n", lev, pr[1], tp_rel(t1, f.theta))
        nf++
    }

    // 4 abandonment: first step too large
    t1 = f.theta + 1e-3 :* (1 :+ abs(f.theta))
    t0 = t1
    pr = dwlsmv_polish(s, w, par, maps, t1)
    if (pr[1] == -2 & t1 == t0) {
        printf("  4. %s: a 1e-3 offset is abandoned (first step %9.2e), theta untouched\n", lev, pr[3])
    }
    else {
        printf("{err}  4. %s FAILED: status %f, theta changed %f\n", lev, pr[1], t1 != t0)
        nf++
    }

    // 5 abandonment: a missing step
    t1 = f.theta
    pr = dwlsmv_polish(s, J(rows(w), 1, 0), par, maps, t1)
    if (pr[1] == -1 & t1 == f.theta) {
        printf("  5. %s: w = 0 gives a missing step; abandoned, theta untouched\n", lev)
    }
    else {
        printf("{err}  5. %s FAILED: status %f\n", lev, pr[1])
        nf++
    }

    // 6 the first-step limit: 3e-6 is above it and contracts if allowed to run
    t1 = f.theta + 3e-6 :* (1 :+ abs(f.theta)) :* (mod(1::rows(f.theta), 2) :* 2 :- 1)
    t0 = t1
    pr = dwlsmv_polish(s, w, par, maps, t1)
    if (pr[1] == -2 & t1 == t0 & pr[3] < 1e-4) {
        printf("  6. %s: a first step of %9.2e is abandoned, theta untouched\n", lev, pr[3])
    }
    else {
        printf("{err}  6. %s FAILED: status %f, first step %9.2e\n", lev, pr[1], pr[3])
        nf++
    }

    // 7 abandonment: a step grows. s + m r keeps theta-hat stationary and scales
    // R by 1 + m, so the contraction rate by 1 + m; m chosen after measuring
    r = s - dwls_model(f.theta, par, maps)
    off = 1e-8 :* (1 :+ abs(f.theta)) :* (mod(1::rows(f.theta), 2) :* 2 :- 1)
    t1 = f.theta + off
    t0 = t1
    pr = dwlsmv_polish(s + 127 :* r, w, par, maps, t1)
    if (pr[1] == -3 & t1 == t0 & pr[3] > 1e-11) {
        printf("  7. %s: a growing step (%9.2e) is abandoned after %f, theta untouched\n", lev, pr[3], pr[2])
    }
    else {
        printf("{err}  7. %s FAILED: status %f, theta changed %f\n", lev, pr[1], t1 != t0)
        nf++
    }

    // 8 the cap, configural only: m = 31 chosen after measuring; at scalar it grows
    if (lev == "configural") {
        t1 = f.theta + off
        t0 = t1
        pr = dwlsmv_polish(s + 31 :* r, w, par, maps, t1)
        th = t0
        for (i = 1; i <= 20; i++) {
            dwlsmv_fdwls(2, th', s + 31 :* r, w, par, maps, fm, gg, H)
            th = th + cholsolve(H, -gg')
        }
        if (pr[1] == 2 & pr[2] == 20 & t1 == th & pr[3] > 1e-11) {
            printf("  8. %s: stops at the cap of 20, theta kept, last step %9.2e\n", lev, pr[3])
        }
        else {
            printf("{err}  8. %s FAILED: status %f, steps %f, theta changed %f\n", lev, pr[1], pr[2], t1 != t0)
            nf++
        }
    }

    // 9 the loop against a hand-written one: same steps, the threshold step applied
    t1 = f.theta + off
    pr = dwlsmv_polish(s + 7 :* r, w, par, maps, t1)
    th = f.theta + off
    for (i = 1; i <= 20; i++) {
        dwlsmv_fdwls(2, th', s + 7 :* r, w, par, maps, fm, gg, H)
        d = cholsolve(H, -gg')
        rel = max(abs(d) :/ (1 :+ abs(th)))
        th = th + d
        if (rel <= 1e-11) break
    }
    if (pr[1] == 0 & pr[2] == i & pr[3] == rel & pr[3] <= 1e-11 & t1 == th) {
        printf("  9. %s: hand-written loop matches bitwise, %f steps\n", lev, i)
    }
    else {
        printf("{err}  9. %s FAILED: status %f, steps %f vs %f, last %9.2e\n", lev, pr[1], pr[2], i, pr[3])
        nf++
    }
    return(nf)
}

// 10 G = 1: dwls_fit's theta is ULS, minimize, polish, signcanon, as in check 1;
// 11 a saturated fit stops by the threshold
real scalar tp_single(string scalar canon, real rowvector K, string rowvector vars, real scalar chk)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f
    real colvector s, w, start, thu, th
    real rowvector pr
    real scalar i, N, fu, iu, fm, it, conv

    sts = dwlsmv_stage1(1, 1)
    maps = dwlsmv_statmap(1, 1)
    sts[1] = dwls_sample_stats(select(st_data(., vars), st_data(., "g") :== 1), K')
    maps[1] = sts[1].map
    par = dwls_identify(subinstr(canon, ";", char(10)), maps, sts)
    f = dwls_fit(sts, par)
    dwlsmv_fit_stack(sts, s, w, N, maps)
    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) start[par.fidx[i]] = par.value[i]
    }
    (void) dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, start, "d2", "nr", 300, thu, fu, iu)
    conv = dwlsmv_minimize(s, w, par, maps, thu, "d2", "nr", 300, th, fm, it)
    pr = dwlsmv_polish(s, w, par, maps, th)
    th = dwlsmv_signcanon(th, par)
    if (th == f.theta & pr[1] == 0 & it == f.iterations & conv == f.converged) {
        printf("  %f. G = 1: dwls_fit is polished (%f steps), bitwise\n", chk, pr[2])
        return(0)
    }
    printf("{err}  %f. G = 1 FAILED: status %f\n", chk, pr[1])
    return(1)
}
end

mata: st_local("nfail", strofreal(tp_level("configural") + tp_level("scalar")))
mata: TPK = st_matrix("e(K)")
quietly dwlsmv (F -> y1 y2 y3 y4 y5) if g == 1, metric(unitvariance)
mata: st_local("nfail", strofreal(`nfail' + tp_single(st_global("e(canon)"), TPK[1, .], ("y1", "y2", "y3", "y4", "y5"), 10)))
quietly dwlsmv (F -> y1 y2 y3) if g == 1, metric(unitvariance)
if e(df_T3) != 0 {
    display as error "  11 precondition FAILED: df " e(df_T3)
    local ++nfail
}
mata: st_local("nfail", strofreal(`nfail' + tp_single(st_global("e(canon)"), TPK[1, 1..3], ("y1", "y2", "y3"), 11)))
adopath - "`here'/src/ado"
adopath - "`here'/build"
if `nfail' {
    display as error "test_polish: `nfail' failure(s)"
    exit 9
}
display as result "test_polish: dwls_fit polishes after convergence, theta-hat is the polish's fixed point, a 1e-7 offset contracts back, and both abandonment paths return theta untouched"
exit 0
