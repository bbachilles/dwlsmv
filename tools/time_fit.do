*! Module 7 timing. NOT a test -- no pass/fail, and the suite runner ignores it
*! (the runner globs tests/unit for test_*.do only).
*!
*!     cd <repo root>
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/time_fit.do
*!     grep "^TIME|" time_fit.log
*!
*! Two questions, both answered by measurement rather than argument:
*!
*!   1. WHICH OPTIMIZER CONFIGURATION dwls_fit() uses. dwlsmv_fdwls() can serve
*!      optimize() as "d1" (gradient only) or "d2" (gradient plus the
*!      Gauss-Newton Hessian 2 Delta' W Delta). Three configurations are timed
*!      end to end, ULS then DWLS, from the same starting values:
*!
*!          d2 / nr     damped Gauss-Newton                 (dwls_fit's choice)
*!          d1 / nr     Newton, H by differencing OUR gradient
*!          d1 / bfgs   quasi-Newton, no H at all
*!
*!      None of the three lets optimize() difference F itself. d1/nr differences
*!      the gradient, which costs npar extra Delta evaluations per iteration.
*!
*!   2. WHETHER ANALYTIC DELTA IS WORTH BUILDING -- Module 6's open call. Answered
*!      by the share of dwls_fit's time spent inside dwls_delta(), set beside the
*!      Stage 1 time on the same data. Analytic Delta can remove at most that
*!      share, and only of Module 7's part.
*!
*! Two models. Case E, the suite's own two-factor fixture (6 indicators,
*! npar 19, D 27), repeated so a single timer tick does not dominate. And a
*! synthetic model sized like a real questionnaire: 40 five-category
*! indicators on 4 correlated factors, unit-variance metric, N = 5000 -- T 160,
*! P 780, D 940, npar 206. Its data are SIMULATED HERE from a stated population
*! under a fixed seed and put through the real Stage 1, so the residuals the
*! optimizer sees are sampling-sized rather than invented. Nothing is compared
*! against anything; the population is stated only so the fit is a
*! realistic one.
*!
*! Timings are wall-clock on whatever machine runs this and move with load.
*! Read ratios, not absolutes.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

foreach f in src/mata/dwlsmv_types.mata             ///
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
             src/mata/dwls_fit.mata {
    quietly do "`f'"
}

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

string scalar tf_line(string scalar s)
{
    return(s + char(10))
}

// ULS then DWLS under one configuration. Returns total seconds over 'reps'
// repetitions; fills the final DWLS iteration count and objective.
real scalar tf_run(struct dwlsmv_stage1 scalar st,
                   struct dwlsmv_parmap  scalar par,
                   string scalar ev, string scalar te, real scalar reps,
                   real scalar itu, real scalar itd, real scalar fm,
                   real scalar conv, real colvector th)
{
    real colvector start, thu, w
    real scalar i, r, fu

    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            start[par.fidx[i]] = par.value[i]
        }
    }
    w = 1 :/ diagonal(st.Gamma)

    timer_clear(90)
    timer_on(90)
    for (r = 1; r <= reps; r++) {
        (void) dwlsmv_minimize(st.s, J(st.map.D, 1, 1), par, st.map, start,
                               ev, te, 300, thu, fu, itu)
        conv = dwlsmv_minimize(st.s, w, par, st.map, thu,
                               ev, te, 300, th, fm, itd)
    }
    timer_off(90)
    return(timer_value(90)[1] / reps)
}

// Seconds per call of dwls_model() and dwls_delta() at theta, averaged.
void tf_unit(struct dwlsmv_stage1 scalar st, struct dwlsmv_parmap scalar par,
             real colvector th, real scalar reps,
             real scalar tm, real scalar td)
{
    real scalar r
    real colvector sg
    real matrix Dm

    timer_clear(91)
    timer_on(91)
    for (r = 1; r <= reps; r++) {
        sg = dwls_model(th, par, st.map)
    }
    timer_off(91)
    tm = timer_value(91)[1] / reps

    timer_clear(92)
    timer_on(92)
    for (r = 1; r <= reps; r++) {
        Dm = dwls_delta(th, par, st.map)
    }
    timer_off(92)
    td = timer_value(92)[1] / reps
}

void tf_model(string scalar tag, struct dwlsmv_stage1 scalar st,
              struct dwlsmv_parmap scalar par, real scalar reps,
              real scalar t_stage1)
{
    struct dwlsmv_fit scalar fit
    string rowvector ev, te
    real scalar k, t, itu, itd, fm, conv, tm, td, tfit
    real colvector th, th1

    printf("TIME|%s|size|T %g P %g D %g npar %g N %g\n",
           tag, st.map.T, st.map.P, st.map.D, par.npar, st.N)
    printf("TIME|%s|stage1|%9.4f s\n", tag, t_stage1)

    ev = ("d2", "d1", "d1")
    te = ("nr", "nr", "bfgs")
    for (k = 1; k <= 3; k++) {
        t = tf_run(st, par, ev[k], te[k], (k == 2 ? max((1, reps / 10)) : reps),
                   itu, itd, fm, conv, th)
        // Compared after dwlsmv_signcanon(): d1/bfgs lands on the REFLECTED
        // solution of the p40 model (theta 1.6 away, F equal), which is the
        // arbitrary sign dwls_fit() canonicalizes, not a disagreement.
        th = dwlsmv_signcanon(th, par)
        if (k == 1) {
            th1 = th
        }
        printf("TIME|%s|%s/%-4s|%9.4f s  ULS it %3.0f  DWLS it %3.0f  conv %g  F %21.15g  max|canonical theta - d2/nr| %8.2e\n",
               tag, ev[k], te[k], t, itu, itd, conv, fm, max(abs(th - th1)))
    }

    // dwls_fit itself, which is what a user pays for
    timer_clear(93)
    timer_on(93)
    for (k = 1; k <= reps; k++) {
        fit = dwls_fit(st, par)
    }
    timer_off(93)
    tfit = timer_value(93)[1] / reps
    printf("TIME|%s|dwls_fit|%9.4f s  (includes the final Delta at theta-hat)\n", tag, tfit)

    tf_unit(st, par, fit.theta, reps, tm, td)
    printf("TIME|%s|per call|dwls_model %9.6f s  dwls_delta %9.6f s  ratio %6.1f\n",
           tag, tm, td, td / tm)
}

void tf_main(real matrix YE)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    string scalar spec
    real matrix Lam, Phi, Ch, Z, Y, eta
    real colvector tau, lam, K
    real scalar p, nf, N, j, f, c, t1, i

    // ---- case E ----------------------------------------------------------
    timer_clear(94)
    timer_on(94)
    st = dwls_sample_stats(YE, J(6, 1, 3))
    timer_off(94)
    t1 = timer_value(94)[1]

    spec = tf_line("factor f1 f2") + tf_line("ind y1 y2 y3 y4 y5 y6")
    for (j = 1; j <= 6; j++) {
        spec = spec + tf_line("lambda f" + strofreal(j <= 3 ? 1 : 2) +
                              " y" + strofreal(j) + " free")
    }
    spec = spec + tf_line("beta f2 f1 free") +
           tf_line("metric unitvariance f1") + tf_line("metric unitvariance f2")
    par = dwls_identify(spec, st.map)
    tf_model("caseE", st, par, 50, t1)

    // ---- synthetic: 40 indicators, 4 factors, K = 5 ----------------------
    p  = 40
    nf = 4
    N  = 5000
    rseed(20260923)

    Phi = J(nf, nf, 0.3) + 0.7 * I(nf)
    lam = J(p, 1, .)
    Lam = J(p, nf, 0)
    for (j = 1; j <= p; j++) {
        f          = ceil(j / (p / nf))
        lam[j]     = 0.5 + 0.3 * mod(j - 1, 10) / 9
        Lam[j, f]  = lam[j]
    }
    tau = (-1.2 \ -0.4 \ 0.4 \ 1.2)

    Ch  = cholesky(Phi)
    eta = rnormal(N, nf, 0, 1) * Ch'
    Z   = eta * Lam' + rnormal(N, p, 0, 1) :* sqrt(1 :- lam' :^ 2)
    Y   = J(N, p, 0)
    for (c = 1; c <= rows(tau); c++) {
        Y = Y + (Z :> tau[c])
    }
    K = J(p, 1, 5)

    timer_clear(94)
    timer_on(94)
    st = dwls_sample_stats(Y, K)
    timer_off(94)
    t1 = timer_value(94)[1]

    spec = tf_line("factor f1 f2 f3 f4")
    spec = spec + "ind"
    for (j = 1; j <= p; j++) {
        spec = spec + " y" + strofreal(j)
    }
    spec = spec + char(10)
    for (j = 1; j <= p; j++) {
        spec = spec + tf_line("lambda f" + strofreal(ceil(j / (p / nf))) +
                              " y" + strofreal(j) + " free")
    }
    for (f = 1; f <= nf; f++) {
        for (i = f + 1; i <= nf; i++) {
            spec = spec + tf_line("psi f" + strofreal(f) + " f" + strofreal(i) + " free")
        }
        spec = spec + tf_line("metric unitvariance f" + strofreal(f))
    }
    par = dwls_identify(spec, st.map)
    tf_model("p40", st, par, 3, t1)
}

tf_main(YE)

end
