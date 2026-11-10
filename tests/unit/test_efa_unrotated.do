*! ESEM groundwork, step 0: an UNROTATED two-factor EFA through v1's pipeline.
*!
*! Executed under StataMP 16 on 2026-10-24 and it passes (the only Stata on
*! this machine). Measured, every ratio of its tier: (N-1)F 2.1e-5, T3 2.2e-4,
*! sigma 2.9e-4, Theta 3.5e-4, loadings 5.1e-4 (on c6-contract-work); gaps ~2e-8. Both
*! anchorings fail the rank check at theta0 and pass at theta2.
*!
*! An unrotated m-factor EFA is a CFA with every loading free except
*! m(m-1)/2 echelon anchor zeros, Phi = I. Nothing here is new mathematics:
*! Stage 1, Gamma, W_D, F, T3 and df are rotation-invariant, so v1 should fit
*! it as it stands. This test checks that claim against lavaan. Rotation and
*! rotated SEs are out of scope. No src/ code is added or changed.
*!
*! Case E, six 3-category items, m = 2: 11 free loadings, df = 4. Written by
*! hand in canonical form (spec 4.8): unitvariance on both factors, no psi
*! f2 f1 line, so Phi = I; one anchor zero lambda f2 <anchor> fixed 0.
*!
*! Reference: fixtures/expected/lavaan_efa_case_E.csv (99g_lavaan_efa.R),
*! lavaan's explicit echelon fit, itself checked there against lavaan's
*! efa(rotation = "none") on T, T3 and sigma.
*!
*! Tolerances are spec 7.2's tiers, chosen before the first run:
*!   F, T3, sigma, Theta, loadings -- theta-hat tier, 1e-5 abs + 1e-4 rel,
*!     except T3 at the test-statistic tier, 1e-6 abs + 1e-5 rel.
*!   df exact.
*!
*! Checks:
*!   1  anchor y1 (lavaan's): dwls_identify() accepts it (rank check passes at
*!      theta0 or theta2; which one is printed), dwls_fit() converges; (N-1)F,
*!      T3, df, sigma(theta-hat) (thresholds and correlations, every entry),
*!      Theta = 1 - diag(Lambda Lambda'), and the loadings up to a column sign
*!      (same anchors as lavaan, so the sign is the only freedom left). The
*!      sign agreement itself is printed, not asserted: see check 3.
*!   2  anchor y4 (F2's own strongest item, a different but equally valid
*!      echelon): F, T3, df, sigma and Theta equal lavaan's at the same tiers.
*!      Rotation invariance of the fit, measured on a second anchoring.
*!   3  sign canonicalization on F2: its anchor under dwlsmv_signcanon() is
*!      its first free loading, here the CROSS-loading on y2 (~0.03 in
*!      lavaan). Printed so the fragility is visible; not asserted.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 120

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_mg_check.mata          ///
             src/mata/dwlsmv_mg_sub.mata            ///
             src/mata/dwlsmv_mg_parts.mata          ///
             src/mata/dwlsmv_mg_sigma.mata          ///
             src/mata/dwlsmv_mg_dcols.mata          ///
             src/mata/dwlsmv_id_parse.mata          ///
             src/mata/dwlsmv_id_rank.mata           ///
             src/mata/dwlsmv_id_tie.mata            ///
             src/mata/dwlsmv_id_groups.mata         ///
             src/mata/dwlsmv_id_opname.mata         ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_catdegen.mata          ///
             src/mata/dwlsmv_cumprop.mata           ///
             src/mata/dwlsmv_thresholds.mata        ///
             src/mata/dwlsmv_psi_tau_table.mata     ///
             src/mata/dwlsmv_cellprob.mata          ///
             src/mata/dwlsmv_dcellprob_rho.mata     ///
             src/mata/dwlsmv_wmat.mata              ///
             src/mata/dwlsmv_crosstab.mata          ///
             src/mata/dwlsmv_pairinfo.mata          ///
             src/mata/dwls_polychoric.mata          ///
             src/mata/dwlsmv_koff.mata              ///
             src/mata/dwlsmv_psi_tau_stack.mata     ///
             src/mata/dwlsmv_psi_rho_table.mata     ///
             src/mata/dwlsmv_psi_matrix.mata        ///
             src/mata/dwlsmv_psibar.mata            ///
             src/mata/dwls_gamma.mata               ///
             src/mata/dwlsmv_emptycell1.mata      ///
             src/mata/dwlsmv_emptycell_repair.mata ///
             src/mata/dwlsmv_emptycell_pairs.mata ///
             src/mata/dwls_sample_stats.mata        ///
             src/mata/dwlsmv_nfac.mata              ///
             src/mata/dwlsmv_parvals.mata           ///
             src/mata/dwlsmv_matrices.mata          ///
             src/mata/dwlsmv_rystar.mata            ///
             src/mata/dwlsmv_corrvec.mata           ///
             src/mata/dwls_model.mata               ///
             src/mata/dwlsmv_delta_cols.mata        ///
             src/mata/dwls_delta.mata               ///
             src/mata/dwlsmv_canon_lines.mata       ///
             src/mata/dwlsmv_canon_find.mata        ///
             src/mata/dwlsmv_canon_status.mata      ///
             src/mata/dwlsmv_numrank.mata           ///
             src/mata/dwls_identify.mata            ///
             src/mata/dwlsmv_fdwls.mata             ///
             src/mata/dwlsmv_minimize.mata          ///
             src/mata/dwlsmv_signcanon.mata         ///
             src/mata/dwls_fit.mata                 ///
             src/mata/dwlsmv_asym.mata              ///
             src/mata/dwlsmv_acov.mata              ///
             src/mata/dwlsmv_traces.mata            ///
             src/mata/dwlsmv_scaleshift.mata        ///
             src/mata/dwls_inference.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_efa_unrotated: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/expected/lavaan_efa_case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: EQ = st_sdata(., "quantity")
mata: EL = st_sdata(., "label")
mata: EV = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

real scalar te_get(string colvector Q, string colvector L, real colvector V,
                   string scalar q, string scalar l)
{
    real colvector k

    k = selectindex((Q :== q) :& (L :== l))
    if (rows(k) != 1 | cols(k) != 1) {
        return(.)
    }
    return(V[k])
}

string scalar te_spec(string scalar anchor)
{
    string scalar s
    real scalar j

    s = "factor f1 f2" + char(10) + "ind y1 y2 y3 y4 y5 y6" + char(10)
    for (j = 1; j <= 6; j++) {
        s = s + "lambda f1 y" + strofreal(j) + " free" + char(10)
        if ("y" + strofreal(j) == anchor) {
            s = s + "lambda f2 y" + strofreal(j) + " fixed 0" + char(10)
        }
        else {
            s = s + "lambda f2 y" + strofreal(j) + " free" + char(10)
        }
    }
    s = s + "metric unitvariance f1" + char(10) + "metric unitvariance f2" + char(10)
    return(s)
}

// Lambda (6 x 2) at theta, from the parameter table.
real matrix te_lambda(real colvector theta, struct dwlsmv_parmap scalar par)
{
    real matrix L
    real scalar i

    L = J(6, 2, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1) {
            L[par.row[i], par.col[i]] = (par.free[i] ? theta[par.fidx[i]] : par.value[i])
        }
    }
    return(L)
}

// worst gap ratio against a abs + r rel; missing if any value is missing
real scalar te_ratio(real colvector x, real colvector y, real scalar a, real scalar r)
{
    if (missing(x) | missing(y) | rows(x) != rows(y)) {
        return(.)
    }
    return(max(abs(x - y) :/ (a :+ r :* abs(y))))
}

real scalar te_run(string colvector Q, string colvector L, real colvector V,
                   real matrix Y)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    string rowvector anchors
    string scalar an
    real colvector sig, lsig, th, lth, lam, llam, sgn
    real matrix Lm
    real scalar a, i, j, nf, N, r, nfree, nst

    nf      = 0
    st      = dwls_sample_stats(Y, J(6, 1, 3))
    N       = st.N
    anchors = ("y1", "y4")

    lsig = J(st.map.D, 1, .)
    for (i = 1; i <= st.map.D; i++) {
        lsig[i] = te_get(Q, L, V, (st.map.stat_type[i] == 1 ? "tau" : "sigma"), st.map.label[i])
    }
    lth = J(6, 1, .)
    for (j = 1; j <= 6; j++) {
        lth[j] = te_get(Q, L, V, "theta", "y" + strofreal(j))
    }
    if (missing(lsig) | missing(lth)) {
        printf("{err}  the fixture did not fill sigma/theta by label\n")
        return(1)
    }

    for (a = 1; a <= cols(anchors); a++) {
        an  = anchors[a]
        par = dwls_identify(te_spec(an), st.map)
        // theta0 passed iff every free loading still starts at 0.7
        nfree = 0
        nst   = 0
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] == 1 & par.free[i]) {
                nfree++
                nst = nst + (par.value[i] == 0.7)
            }
        }
        printf("{txt}  anchor %s: npar %g, rank check passed at %s\n",
               an, par.npar, (nst == nfree ? "theta0" : "theta2 (theta0 deficient)"))

        fit = dwls_fit(st, par)
        if (fit.converged != 1) {
            printf("{err}  anchor %s: dwls_fit() did not converge\n", an)
            nf++
            continue
        }
        dwls_inference(fit)

        r = dwlsmv_gapratio((N - 1) * fit.fmin, te_get(Q, L, V, "test", "T"), 1e-5, 1e-4)
        printf("{txt}    (N-1)F %16.12f vs %16.12f  ratio %8.2e\n",
               (N - 1) * fit.fmin, te_get(Q, L, V, "test", "T"), r)
        if (!(r <= 1)) nf++

        r = dwlsmv_gapratio(fit.T3, te_get(Q, L, V, "test", "T3"), 1e-6, 1e-5)
        printf("{txt}    T3     %16.12f vs %16.12f  ratio %8.2e\n",
               fit.T3, te_get(Q, L, V, "test", "T3"), r)
        if (!(r <= 1)) nf++

        printf("{txt}    df     %g vs %g\n", fit.T3_df, te_get(Q, L, V, "test", "df"))
        if (fit.T3_df != te_get(Q, L, V, "test", "df")) nf++

        sig = dwls_model(fit.theta, par, st.map)
        r   = te_ratio(sig, lsig, 1e-5, 1e-4)
        printf("{txt}    sigma(theta-hat), %g entries, ratio %8.2e (max abs %8.2e)\n",
               rows(sig), r, max(abs(sig - lsig)))
        if (!(r <= 1)) nf++

        Lm = te_lambda(fit.theta, par)
        th = 1 :- rowsum(Lm :^ 2)
        r  = te_ratio(th, lth, 1e-5, 1e-4)
        printf("{txt}    Theta, ratio %8.2e (max abs %8.2e)\n", r, max(abs(th - lth)))
        if (!(r <= 1)) nf++

        printf("{txt}    F2 sign anchor (first free F2 loading): %s = %9.6f\n",
               (an == "y1" ? "y2" : "y1"), (an == "y1" ? Lm[2, 2] : Lm[1, 2]))

        if (an == "y1") {
            llam = J(12, 1, .)
            for (j = 1; j <= 6; j++) {
                llam[j]     = te_get(Q, L, V, "lambda", "f1=~y" + strofreal(j))
                llam[6 + j] = te_get(Q, L, V, "lambda", "f2=~y" + strofreal(j))
            }
            // the column sign is the one freedom the anchors leave
            sgn = J(2, 1, 1)
            for (j = 1; j <= 2; j++) {
                i = (j == 1 ? 1 : 4)
                if (Lm[i, j] * llam[(j - 1) * 6 + i] < 0) sgn[j] = -1
            }
            lam = vec(Lm :* sgn')
            r   = te_ratio(lam, llam, 1e-5, 1e-4)
            printf("{txt}    loadings (column signs %g, %g vs lavaan), ratio %8.2e (max abs %8.2e)\n",
                   sgn[1], sgn[2], r, max(abs(lam - llam)))
            if (!(r <= 1)) nf++
        }
    }
    return(nf)
}

end

mata: st_local("nfail", strofreal(te_run(EQ, EL, EV, YE)))

if `nfail' {
    display as error _n "test_efa_unrotated: `nfail' failure(s)"
    exit 9
}
display as result "test_efa_unrotated: unrotated 2-factor EFA on Case E matches lavaan (F, T3, df, sigma, Theta, loadings) under two echelon anchors"
exit 0
