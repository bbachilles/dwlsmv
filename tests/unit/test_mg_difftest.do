*! Executed under StataMP 18 on 2026-10-14: passes (second run; see the redesign note).
*!
*! DIFFTEST statistic at G > 1 (gate step 3b): dwlsmv_difftest() fed the
*! stacked Delta1, Delta0, Omega = dwlsmv_omega(st) and Td = T0 - T1 with each
*! T = (N - G)F, against the tier-2 fixture 99e
*! (fixtures/generate/99e_lavaan_mg_difftest.R; lavaan 0.7.2; designs E, B, D,
*! two groups, N 1200/2800). Unlike test_difftest.do (G = 1, fed lavaan's own
*! matrices), every input here is ours: one Stage 1 per design, shared by its
*! four levels, each level fitted from dwlsmv_levelgen's text inside Mata.
*!
*! Pairs: the 12 fixture rows (configural-thresholds, thresholds-metric,
*! metric-scalar, configural-scalar per design). The d = 0 rows (B
*! thresholds-configural, D metric-thresholds) are lavaan placeholders: only d
*! is asserted on them, and they never reach the statistic.
*!
*! Checks, on every pair with d >= 1:
*!   1  d = q1 - q0 exactly the fixture's (all 12 rows).
*!   2  T_d3 against T_d3_NG (the fixture's, from (N - G)(F0 - F1)) within
*!        tol = a_d (B0 + B1) + 1e-5 a_d |Td| + 1e-6
*!      Bj is fit j's propagated T bound, test_fit_groups.do check 3's:
*!      (N - G)(2 sum w |r| (1e-8 + 1e-8 |s|) + 1e-8 F), s moving by its tier;
*!      1e-5 and 1e-6 are spec 7.2's a and b rows. Set from the architecture
*!      BEFORE measuring (E-027, 2026-10-14).
*!   3  Tr(M1), Tr(M1^2), a_d at 1e-5 relative; b_d at 1e-6 absolute (7.2).
*!   4  negative control: lavaan's own-convention T_d3 (Td from
*!      sum_g (N_g - 1) F_g) lies OUTSIDE tol, so check 2 can tell (N - G)
*!      from C. If it does not, the bound is too wide to serve: stop and return
*!      to the user (E-027, 2026-10-14), never tighten it here.
*!   5  scale: Omega / N gives every output to 1e-12 (relative; b_d absolute)
*!      (identity; the derivation's mg-difftest-omega figure).
*!   6  negative control: blockdiag(Gamma_g), no pi, moves T_d3 by more than
*!      1e-3 relative (the derivation's figure, chosen after measuring there).
*!   7  p against p_NG within tol x fmax, fmax the largest chi2(d) density on
*!      [T_d3 - tol, T_d3 + tol] (ends, and the mode d - 2 when inside): the
*!      mean value theorem carries check 2's bound to p.
*!
*!   REDESIGNED AFTER THE FIRST RUN (2026-10-14), a prediction of mine that was
*!   wrong, kept here: check 7 first used tol / 2 and asserted every tested
*!   d >= 2 (the chi2 density is <= 1/2 only there). The fixture has d = 1 on D
*!   configural-thresholds, so that assertion failed (the p gap itself was
*!   3.5e-9 against tol / 2 = 3.0e-5). The density bound above holds for every
*!   d and is never wider than tol / 2 where that applied. No other check
*!   changed; first-run ratios: check 2 worst 3.3e-2, check 4 narrowest 2.10x
*!   (D configural-thresholds), check 5 worst 8.2e-14, check 6 narrowest 4.8e-3.
*!
*! Mutants on dwlsmv_difftest (2026-10-14, scratch worktree, diff printed
*! before each run, %TEMP%\exgate\mutate3b.sh), all red:
*!   X1  W-weighted h-hat         7 pairs (not the configural-H1 ones, as MG7(2)
*!                                predicts: there every left inverse agrees)
*!   X2  shift b_d dropped        9 pairs (all but D configural-thresholds,
*!                                whose b_d is 2.7e-15)
*!   X3  unit weights             every pair, checks 2, 3, 5, 7
*! The caller-side mutants (Td in lavaan's convention, Omega without pi, Delta0
*! at theta1, the pair swapped, the statistic at d = 0 or status 2) are
*! test_dtpair.do's.
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

* reference fixture: stringcols(_all) + strtoreal(), delimiters(",") always
import delimited using "fixtures/expected/lavaan_mg_difftest.csv", ///
    clear varnames(1) case(lower) delimiters(",") stringcols(_all)
mata: R_KEY = st_sdata(., ("design", "h1", "h0"))
mata: R_V   = strtoreal(st_sdata(., ("d", "td", "td_ng", "tr_m1", "tr_m1sq", "a_d", "b_d", "t_d3", "t_d3_ng", "p_ng")))

mata:

struct tmd_lev {
    real colvector theta
    real matrix    Delta
    real scalar    T, B, npar, conv
}

void tmd_keep(string scalar X)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    struct tmd_lev vector L
    string rowvector lev
    string scalar cc
    real matrix K, Y
    real colvector gv, s, w, r
    real scalar G, g, k, N
    pointer() scalar p

    K   = st_matrix("e(K)")
    cc  = st_global("e(canon_configural)")
    G   = rows(K)
    gv  = uniqrows(st_data(., "g"))
    sts = dwlsmv_stage1(G, 1)
    maps = dwlsmv_statmap(G, 1)
    N = 0
    s = J(0, 1, .)
    w = J(0, 1, .)
    for (g = 1; g <= G; g++) {
        Y       = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
        sts[g]  = dwls_sample_stats(Y, K[g, .]')
        maps[g] = sts[g].map
        N = N + sts[g].N
    }
    for (g = 1; g <= G; g++) {
        s = s \ sts[g].s
        w = w \ ((sts[g].N / N) :/ diagonal(sts[g].Gamma))
    }
    lev = ("configural", "thresholds", "metric", "scalar")
    L   = tmd_lev(4, 1)
    for (k = 1; k <= 4; k++) {
        par = dwls_identify(dwlsmv_levelgen(cc, K, lev[k]), maps, sts)
        fit = dwls_fit(sts, par)
        r   = s - dwls_model(fit.theta, par, maps)
        L[k].theta = fit.theta
        L[k].Delta = fit.Delta
        L[k].T     = fit.T_raw
        L[k].npar  = par.npar
        L[k].conv  = fit.converged
        L[k].B     = (N - G) * (2 * quadcross(w, abs(r) :* (1e-8 :+ 1e-8 :* abs(s))) + 1e-8 * fit.fmin)
    }
    p = crexternal("TMD_L_" + X)
    if (p == NULL) p = findexternal("TMD_L_" + X)
    *p = L
    p = crexternal("TMD_S_" + X)
    if (p == NULL) p = findexternal("TMD_S_" + X)
    *p = sts
}

real scalar tmd_k(string scalar lev)
{
    if (lev == "configural") return(1)
    if (lev == "thresholds") return(2)
    if (lev == "metric")     return(3)
    if (lev == "scalar")     return(4)
    _error("test_mg_difftest: unknown level")
    return(.)
}

real scalar tmd_rel(real scalar a, real scalar b)
{
    return(abs(a - b) / max((abs(b), 1e-300)))
}

// largest chi2(d) density on [T - tol, T + tol]
real scalar tmd_fmax(real scalar d, real scalar T, real scalar tol)
{
    real scalar lo, hi, f

    lo = max((T - tol, 1e-300))
    hi = T + tol
    f  = max((chi2den(d, lo), chi2den(d, hi)))
    if (d > 2) {
        if (lo <= d - 2 & d - 2 <= hi) {
            f = max((f, chi2den(d, d - 2)))
        }
    }
    return(f)
}

real scalar tmd_run()
{
    external string matrix R_KEY
    external real matrix R_V
    struct tmd_lev vector L
    struct dwlsmv_stage1 vector sts
    string colvector FAIL
    string scalar X, tag
    real matrix Om, Onp
    real rowvector r, rN, rp
    real scalar i, k1, k0, d, Td, tol, Tc, ntest, g, o, n, worst2, worst5

    FAIL   = J(0, 1, "")
    ntest  = 0
    worst2 = 0
    worst5 = 0
    for (i = 1; i <= rows(R_KEY); i++) {
        X   = R_KEY[i, 1]
        tag = X + " " + R_KEY[i, 2] + "-" + R_KEY[i, 3]
        L   = *findexternal("TMD_L_" + X)
        sts = *findexternal("TMD_S_" + X)
        k1  = tmd_k(R_KEY[i, 2])
        k0  = tmd_k(R_KEY[i, 3])
        if (L[k1].conv != 1 | L[k0].conv != 1) {
            FAIL = FAIL \ (tag + ": a fit did not converge")
            continue
        }

        // ---- 1 -------------------------------------------------------------
        d = L[k1].npar - L[k0].npar
        printf("{txt}  %-26s 1. d %f (fixture %f)\n", tag, d, R_V[i, 1])
        if (d != R_V[i, 1]) {
            FAIL = FAIL \ (tag + ": d is not the fixture's")
            continue
        }
        if (d == 0) {
            printf("{txt}  %-26s    d = 0: lavaan placeholder row, not sent to the statistic\n", tag)
            continue
        }
        ntest++

        Om = dwlsmv_omega(sts)
        Td = L[k0].T - L[k1].T
        r  = dwlsmv_difftest(L[k1].Delta, L[k0].Delta, Om, Td)

        // ---- 2 -------------------------------------------------------------
        tol = r[4] * (L[k0].B + L[k1].B) + 1e-5 * r[4] * abs(Td) + 1e-6
        printf("{txt}  %-26s 2. T_d3 %14.10f ; fixture NG %14.10f ; gap %8.2e ; tol %8.2e ; ratio %8.2e ; Td %12.8f (NG %12.8f)\n",
               tag, r[1], R_V[i, 9], abs(r[1] - R_V[i, 9]), tol, abs(r[1] - R_V[i, 9]) / tol, Td, R_V[i, 3])
        worst2 = max((worst2, abs(r[1] - R_V[i, 9]) / tol))
        if (abs(r[1] - R_V[i, 9]) > tol) {
            FAIL = FAIL \ (tag + ": T_d3 disagrees with the fixture's (N - G) statistic")
        }

        // ---- 3 -------------------------------------------------------------
        printf("{txt}  %-26s 3. trM1 rel %8.2e ; trM1^2 rel %8.2e ; a rel %8.2e ; b abs %8.2e\n",
               tag, tmd_rel(r[6], R_V[i, 4]), tmd_rel(r[7], R_V[i, 5]),
               tmd_rel(r[4], R_V[i, 6]), abs(r[5] - R_V[i, 7]))
        if (tmd_rel(r[6], R_V[i, 4]) > 1e-5 | tmd_rel(r[7], R_V[i, 5]) > 1e-5 |
            tmd_rel(r[4], R_V[i, 6]) > 1e-5 | abs(r[5] - R_V[i, 7]) > 1e-6) {
            FAIL = FAIL \ (tag + ": traces, a_d or b_d disagree with the fixture")
        }

        // ---- 4 -------------------------------------------------------------
        Tc = R_V[i, 8]
        printf("{txt}  %-26s 4. lavaan-convention T_d3 %14.10f ; gap %8.2e = %8.2f x tol (must exceed 1)\n",
               tag, Tc, abs(r[1] - Tc), abs(r[1] - Tc) / tol)
        if (abs(r[1] - Tc) <= tol) {
            FAIL = FAIL \ (tag + ": the bound cannot tell (N - G) from C; stop and return to the user")
        }

        // ---- 5 -------------------------------------------------------------
        n  = sts[1].N + sts[2].N
        rN = dwlsmv_difftest(L[k1].Delta, L[k0].Delta, Om :/ n, Td)
        o  = max((tmd_rel(rN[1], r[1]), tmd_rel(rN[3], r[3]), tmd_rel(rN[4], r[4]),
                  abs(rN[5] - r[5]), tmd_rel(rN[6], r[6]), tmd_rel(rN[7], r[7])))
        worst5 = max((worst5, o))
        printf("{txt}  %-26s 5. Omega / N: worst %8.2e (1e-12)\n", tag, o)
        if (o > 1e-12 | rN[2] != r[2]) {
            FAIL = FAIL \ (tag + ": Omega / N changes the statistic")
        }

        // ---- 6 -------------------------------------------------------------
        Onp = J(rows(Om), cols(Om), 0)
        o   = 0
        for (g = 1; g <= rows(sts); g++) {
            Onp[(o + 1)::(o + rows(sts[g].Gamma)), (o + 1)::(o + rows(sts[g].Gamma))] = sts[g].Gamma
            o = o + rows(sts[g].Gamma)
        }
        rp = dwlsmv_difftest(L[k1].Delta, L[k0].Delta, Onp, Td)
        printf("{txt}  %-26s 6. no-pi T_d3 rel change %8.2e (must exceed 1e-3)\n", tag, tmd_rel(rp[1], r[1]))
        if (tmd_rel(rp[1], r[1]) <= 1e-3) {
            FAIL = FAIL \ (tag + ": blockdiag(Gamma_g) without pi gives the same statistic")
        }

        // ---- 7 -------------------------------------------------------------
        printf("{txt}  %-26s 7. p %12.10f ; fixture p_NG %12.10f ; gap %8.2e ; bound %8.2e\n",
               tag, r[3], R_V[i, 10], abs(r[3] - R_V[i, 10]), tol * tmd_fmax(d, r[1], tol))
        if (abs(r[3] - R_V[i, 10]) > tol * tmd_fmax(d, r[1], tol)) {
            FAIL = FAIL \ (tag + ": p disagrees with the fixture's")
        }
    }
    printf("{txt}  tested pairs %f ; worst check-2 ratio %8.2e ; worst check-5 %8.2e\n", ntest, worst2, worst5)
    if (ntest != 10) {
        FAIL = FAIL \ ("expected 10 pairs with d >= 1, tested " + strofreal(ntest))
    }
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}  FAIL: %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

foreach X in E B D {
    import delimited using "fixtures/data/mg_`X'.csv", clear varnames(1) delimiters(",")
    quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
    mata: tmd_keep("`X'")
    display as text "    design `X': four levels fitted"
}

mata: st_local("nf", strofreal(tmd_run()))
if `nf' {
    display as error "test_mg_difftest: `nf' failure(s)"
    exit 9
}
display as result "test_mg_difftest: dwlsmv_difftest() fed stacked Delta, Omega = blockdiag(Gamma_g / pi_g) and Td on (N - G) F reproduces 99e's (N - G) statistic on designs E, B, D within the propagated bound, tells it from lavaan's convention, is invariant to Omega's scale, and needs pi"
