*! Module 7: dwls_fit() -- the DWLS optimizer.
*!
*! Executed under StataMP 18 on 2026-09-23 and it passes.
*!
*! WHAT THIS FILE IS FOR. Spec Section 5, Module 7: minimize
*! F = (s - sigma)' W_D (s - sigma), W_D = diag(Gamma)^-1, from the ULS
*! solution, and validate the point estimates against lavaan. Five fits: case D
*! (single factor, marker) and case E under all four of its lavaan fits --
*! correlated and F2 <- F1, each under unit-variance and marker metrics.
*!
*!   1   the gradient. dwlsmv_fdwls()'s analytic g = -2 r' W Delta against a
*!       central difference of F itself, at the STARTING values, where g is
*!       far from zero. This is the check that would catch a sign or a factor
*!       of two; it cannot catch an inverted W, since both sides use the same
*!       w. Check 5 does that.
*!   2   convergence from ULS starts. ULS (W = I) and DWLS each declare
*!       convergence on all five fits.
*!   3   theta-hat against lavaan's DWLS/WLSMV estimates (the point estimates
*!       are the same under both), every free parameter, thresholds included,
*!       at spec Section 7.2's theta-hat tier: 1e-5 absolute, 1e-4 relative.
*!   4   sigma(theta-hat) against lavaan's wls.est on case D, same tier.
*!   5   THE WEIGHT CONVENTION AND T's SCALING, against lavaan's unscaled
*!       chi-square. lavaan 0.7.2 reports (N - 1) F-hat -- read from its source,
*!       lav_model_objective and lav_model_test -- and T_raw follows it (spec
*!       4.9, decided 2026-09-23); Mplus's N F-hat is 1/N away. The check
*!       asserts BOTH halves: T_raw matches, and N F-hat does not. The
*!       tolerance is not a tier picked from the table but propagated: F-hat
*!       moves by at most 2 sum_d w_d |r_d| |ds_d| + F |dw/w| when s and Gamma
*!       move by their end-to-end tiers (1e-8 abs + 1e-8 rel on s, 1e-8 rel on
*!       Gamma's diagonal), computed from this fit's own w and r. An inverted
*!       W moves theta-hat and F-hat by orders of magnitude and fails at once.
*!   6   thresholds. Delta is block diagonal (Gate 2, Finding 3) and W_D is
*!       diagonal, so F separates into sum w (s_tau - tau)^2 plus a part free
*!       of tau: at the optimum tau-hat = s_tau exactly. A Gauss-Newton step
*!       solves that linear block exactly, so the measured gap is 0.0 (asserted at 1e-15).
*!       Both the residual s_tau - sigma_tau and the parameters tau-hat - s_tau
*!       (the latter added 2026-10-07: a changed sigma_tau keeps the residual
*!       at 0).
*!   7   stationarity: max |g(theta-hat)| is at most 1e-8 of max |g(start)|.
*!   8   NON-CONVERGENCE IS REPORTED, NOT FATAL. dwlsmv_minimize() capped at
*!       one iteration returns 0, keeps the last iterate, and does not abort.
*!   9   the inputs are not written. st and par are compared bitwise before
*!       and after dwls_fit(). The routine-by-routine audit is recorded in
*!       dwlsmv_fdwls()'s header; this is the executable half of it.
*!   10  the fit's own bookkeeping: D_stat, Delta at theta-hat,
*!       T_raw = (N - 1) F, and Modules 8-9's fields left untouched.
*!   12  ADMISSIBILITY (spec Section 5), through dwlsmv_admissible() and
*!       dwlsmv_admissadvice(). 12a: every fitted solution above passes all
*!       six conditions and gets an empty note. 12b-f: five inadmissible
*!       points built by hand, each with an exactly known answer -- a Heywood
*!       case, a factor correlation of 1.2, reversed thresholds, a singular
*!       (I - B), and an R_y* with a correlation of 1.05, which is the only
*!       case that drives row 4 -- and the rendered notes are printed, because
*!       a strpos() test once passed on unreadable text (dwlsmv_celladvice's
*!       header). 12d and 12f were added after the 2026-09-23 audit found row 4
*!       never fired and the row 4 and row 6 notes never rendered.
*!   13  FACTOR SIGNS, canonicalized after convergence by dwlsmv_signcanon().
*!       Reflections spelled out by hand -- a correlated factor, an
*!       endogenous one (beta flips), both at once (beta does not) -- come
*!       back to theta-hat bit for bit; a marker factor, whose sign a fixed
*!       loading identifies, is left alone; and a fit started at negative
*!       loadings, which the unconstrained optimizer carries to the reflected
*!       solution, is returned canonical by dwls_fit(). 13e-h, added after the
*!       audit: a FIXED nonzero psi_12 blocks the reflection while a fixed 0
*!       does not; three factors, with a beta joining two reversed factors
*!       left untouched and sigma checked unchanged; an anchor of exactly 0;
*!       and a one-factor unit-variance model.
*!
*!       NOT REACHED, and why: converged = 0 THROUGH dwls_fit() (its maxiter
*!       is fixed at 300 and case D/E converge in under 10, so only
*!       dwlsmv_minimize() is capped, in check 8); ULS non-convergence
*!       inside dwls_fit() for the same reason; and a non-zero _optimize()
*!       code with optimize_result_converged() = 1, which no input here
*!       produces.
*!   14  THE RANK CHECK'S SECOND POINT IS THE START (spec 4.8, 2026-10-04): a
*!       random-intercept model deficient only at the symmetric starts
*!       starts at theta2, converges admissibly from the default starts, and
*!       agrees with a distinct-start fit; a model full rank at the starts
*!       keeps them. From theta0 the fit did not converge.
*!   11  the error paths -- all seven _error() sites in the Module 7 files,
*!       signcanon's and admissadvice's added after the audit found them
*!       unvisited -- at rc 3498 and not merely nonzero: an errprintf
*!       with a format/argument mismatch aborts at rc 3001 AT THE PRINTF, and
*!       only a runtime visit tells the two apart (a project coding rule).
*!
*! The correlated/structural equivalence on OUR fits is not here: it is spec
*! Section 7.4's structural-path condition and lives in
*! test_structural_metric.do beside the lavaan half it completes.
*!
*! Reference CSVs are read with stringcols(_all) and converted with
*! strtoreal(), and with delimiters(",") specified, for the two reasons
*! The project's coding rules record: float truncation near 4e-08, and '|' auto-detected as
*! the delimiter when labels like y1|t1 are present -- which they are here.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

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
             src/mata/dwlsmv_admissible.mata        ///
             src/mata/dwlsmv_admissadvice.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_fit: `f' failed to load"
        exit 9
    }
}

* ---- tier 2: lavaan's case D estimates and model-implied vector -------------
import delimited using "fixtures/expected/lavaan_case_D_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: DKEY = st_sdata(., "label")
mata: DVAL = strtoreal(st_sdata(., "est"))

import delimited using "fixtures/expected/lavaan_case_D_wls_est.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: WKEY = st_sdata(., "label")
mata: WVAL = strtoreal(st_sdata(., "value"))

* ---- tier 2: lavaan's four case E fits, estimates and unscaled chi-square ---
import delimited using "fixtures/expected/lavaan_case_E_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: PFIT = st_sdata(., "fit")
mata: PKEY = st_sdata(., "label")
mata: PVAL = strtoreal(st_sdata(., "est"))

import delimited using "fixtures/expected/lavaan_case_E_equivalence.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: EFIT = st_sdata(., "fit")
mata: ECHI = strtoreal(st_sdata(., "chisq"))

* ---- the data ----------------------------------------------------------------
import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

// Long-format lookup; missing when absent, which every caller reports.
real scalar tf_get(string colvector fitcol, string colvector keycol,
                   real colvector valcol, string scalar wfit, string scalar wkey)
{
    real scalar i

    for (i = 1; i <= rows(keycol); i++) {
        if ((wfit == "" | fitcol[i] == wfit) & keycol[i] == wkey) {
            return(valcol[i])
        }
    }
    return(.)
}

string scalar tf_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The five fits, as hand-written canonical specifications -- the same
// models lavaan fitted in fixtures/generate/20_lavaan_reference.R. lavaan's
// marker is the first-listed indicator and std.lv is unit variance on every
// factor; both are written out explicitly here, as the canonical form
// requires.
string scalar tf_canon(string scalar fit)
{
    string rowvector head, lam6

    head = ("factor f1 f2", "ind y1 y2 y3 y4 y5 y6")
    lam6 = ("lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
            "lambda f2 y4 free", "lambda f2 y5 free", "lambda f2 y6 free")

    if (fit == "D_marker") {
        return(tf_spec(("factor f", "ind y1 y2 y3 y4",
                        "lambda f y2 free", "lambda f y3 free",
                        "lambda f y4 free", "psi f f free",
                        "metric marker f y1")))
    }
    if (fit == "cor_unitvar") {
        return(tf_spec((head, lam6, "psi f1 f2 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    if (fit == "str_unitvar") {
        return(tf_spec((head, lam6, "beta f2 f1 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    if (fit == "cor_marker") {
        return(tf_spec((head,
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "psi f1 f1 free", "psi f2 f2 free", "psi f1 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    if (fit == "str_marker") {
        return(tf_spec((head,
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "beta f2 f1 free",
                        "psi f1 f1 free", "psi f2 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    _error("tf_canon: unknown fit")
}

// lavaan's label for one row of our table. Lambda is (indicator, factor);
// beta is (to, from); psi is symmetric and lavaan names the lower index first.
string scalar tf_lab(struct dwlsmv_parmap scalar par, real scalar i,
                     string rowvector fac)
{
    real scalar r, c

    r = par.row[i]
    c = par.col[i]
    if (par.ptype[i] == 4) return(par.name[i])
    if (par.ptype[i] == 1) return(fac[c] + "=~y" + strofreal(r))
    if (par.ptype[i] == 3) return(fac[r] + "~" + fac[c])
    return(fac[min((r, c))] + "~~" + fac[max((r, c))])
}

real colvector tf_start(struct dwlsmv_parmap scalar par)
{
    real colvector st0
    real scalar i

    st0 = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            st0[par.fidx[i]] = par.value[i]
        }
    }
    return(st0)
}

real scalar tf_F(real colvector th, real colvector s, real colvector w,
                 struct dwlsmv_parmap scalar par, struct dwlsmv_statmap scalar map)
{
    real colvector r

    r = s - dwls_model(th, par, map)
    return(quadcross(r, w, r))
}

real rowvector tf_grad(real colvector th, real colvector s, real colvector w,
                       struct dwlsmv_parmap scalar par,
                       struct dwlsmv_statmap scalar map)
{
    real scalar v
    real rowvector g
    real matrix H

    dwlsmv_fdwls(1, th', s, w, par, map, v, g, H)
    return(g)
}

// A case E theta built by hand: thresholds (-0.5, 0.5) on every indicator
// except as overridden, loadings lam, and any psi(1,2) or beta row set to the
// given value. Walks the table rather than assuming fidx positions.
real colvector tf_theta(struct dwlsmv_parmap scalar par, real colvector lam,
                        real scalar link, real colvector tau2)
{
    real colvector th
    real scalar i

    th = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        if (par.ptype[i] == 4) {
            th[par.fidx[i]] = (par.col[i] == 1 ? -0.5 : 0.5)
            if (par.row[i] == 2) {
                th[par.fidx[i]] = tau2[par.col[i]]
            }
        }
        else if (par.ptype[i] == 1) {
            th[par.fidx[i]] = lam[par.row[i]]
        }
        else {
            th[par.fidx[i]] = link
        }
    }
    return(th)
}

// ---- 12b. the inadmissible cases, each at a stated theta with a known answer
string colvector tf_admiss(struct dwlsmv_statmap scalar map)
{
    struct dwlsmv_parmap scalar pc, ps
    string colvector FAIL, ind, fac
    string scalar msg
    real matrix A
    real colvector lam, ok2
    real scalar n

    FAIL = J(0, 1, "")
    ind  = ("y1" \ "y2" \ "y3" \ "y4" \ "y5" \ "y6")
    fac  = ("f1" \ "f2")
    pc   = dwls_identify(tf_canon("cor_unitvar"), map)
    ps   = dwls_identify(tf_canon("str_unitvar"), map)
    lam  = J(6, 1, 0.7)
    ok2  = (-0.5 \ 0.5)

    // Heywood: lambda_y4 = 1.1 on a unit-variance factor, so
    // Theta_44 = 1 - 1.21 = -0.21 exactly (to rounding).
    lam[4] = 1.1
    A   = dwlsmv_admissible(tf_theta(pc, lam, 0.3, ok2), pc, map)
    msg = dwlsmv_admissadvice(A, ind, fac)
    printf("{txt} 12b. Heywood on y4             : flag %g, Theta %g at indicator %g; the note reads:\n%s",
           A[1, 1], A[1, 2], A[1, 3], msg)
    if (A[1, 1] != 1 | A[1, 3] != 4 | abs(A[1, 2] + 0.21) > 1e-15 |
        A[5, 1] != 0 | A[6, 1] != 0) {
        FAIL = FAIL \ "admissibility: the Heywood case on y4 is not reported as such"
    }
    if (!strpos(msg, "y4:") | !strpos(msg, "y4: residual variance -0.2100, below 0") | !strpos(msg, "Heywood") |
        strpos(msg, char(92) + "n")) {
        FAIL = FAIL \ "admissibility: the Heywood note does not read 'y4: ... = -0.2100 <', or carries a literal backslash-n"
    }

    // A factor correlation of 1.2: psi(1,2) = 1.2 between two unit-variance
    // exogenous factors. Phi = [1, 1.2; 1.2, 1], eigenvalues 2.2 and -0.2.
    lam[4] = 0.7
    A   = dwlsmv_admissible(tf_theta(pc, lam, 1.2, ok2), pc, map)
    msg = dwlsmv_admissadvice(A, ind, fac)
    printf("{txt} 12c. factor correlation 1.2     : corr flag %g value %g (%g,%g), Phi flag %g min eig %g\n",
           A[2, 1], A[2, 2], A[2, 3], A[2, 4], A[3, 1], A[3, 2])
    if (A[2, 1] != 1 | abs(A[2, 2] - 1.2) > 1e-15 | A[2, 3] != 1 | A[2, 4] != 2 |
        A[3, 1] != 1 | abs(A[3, 2] + 0.2) > 1e-15 | A[1, 1] != 0) {
        FAIL = FAIL \ "admissibility: a factor correlation of 1.2 is not reported as such"
    }
    if (!strpos(msg, "f1 <-> f2") | !strpos(msg, "not positive definite")) {
        FAIL = FAIL \ "admissibility: the factor-correlation note does not name f1 <-> f2"
    }

    // Disordered thresholds on y2 -- reachable at an arbitrary theta, never at
    // a converged one (check 6).
    A   = dwlsmv_admissible(tf_theta(ps, lam, 0.4, (0.5 \ -0.5)), ps, map)
    msg = dwlsmv_admissadvice(A, ind, fac)
    printf("{txt} 12d. y2 thresholds reversed    : flag %g, gap %g at indicator %g; the note reads:\n%s",
           A[6, 1], A[6, 2], A[6, 3], msg)
    if (A[6, 1] != 1 | A[6, 2] != -1 | A[6, 3] != 2 | A[1, 1] != 0) {
        FAIL = FAIL \ "admissibility: reversed thresholds on y2 are not reported"
    }
    if (!strpos(msg, "y2: thresholds are not in increasing order (smallest gap -1.0000)")) {
        FAIL = FAIL \ "admissibility: the reversed-threshold note is not rendered as expected"
    }

    // R_y* not positive definite -- the one row no earlier case drove to 1.
    // lambda_y4 = 1.5 on a unit-variance factor gives R(y4, y5) =
    // 1.5 x 0.7 = 1.05, a correlation above 1, so R_y* cannot be positive
    // definite. As the header says, row 4 arrives WITH row 1 here
    // (Theta_44 = 1 - 2.25 = -1.25): a unit-diagonal R that is not PD needs
    // some communality above 1.
    lam[4] = 1.5
    A   = dwlsmv_admissible(tf_theta(pc, lam, 0.3, ok2), pc, map)
    msg = dwlsmv_admissadvice(A, ind, fac)
    printf("{txt} 12f. R_y* not PD (lambda 1.5)  : R flag %g min eig %g ; Theta flag %g value %g; the note reads:\n%s",
           A[4, 1], A[4, 2], A[1, 1], A[1, 2], msg)
    if (A[4, 1] != 1 | !(A[4, 2] < 0) | A[1, 1] != 1 | abs(A[1, 2] + 1.25) > 1e-15) {
        FAIL = FAIL \ "admissibility: an R_y* with a correlation of 1.05 is not reported as not positive definite"
    }
    if (!strpos(msg, "underlying continuous") | !strpos(msg, "y4:")) {
        FAIL = FAIL \ "admissibility: the R_y* note is not rendered"
    }
    lam[4] = 0.7

    // Singular (I - B): beta(2,1) = 1 free and a hand-appended FIXED
    // beta(1,2) = 1, so I - B = [1, -1; -1, 1]. No canonical specification
    // reaches this -- the rank check refuses the reciprocal pair at its zero
    // starts -- which is why the row is appended by hand, as test_delta.do
    // does. Rows 1-4 must come back MISSING, not passing.
    n           = ps.nall + 1
    ps.nall     = n
    ps.name     = ps.name  \ "beta[f1,f2]"
    ps.ptype    = ps.ptype \ 3
    ps.row      = ps.row   \ 1
    ps.col      = ps.col   \ 2
    ps.free     = ps.free  \ 0
    ps.value    = ps.value \ 1
    ps.fidx     = ps.fidx  \ 0
    A   = dwlsmv_admissible(tf_theta(ps, lam, 1, ok2), ps, map)
    msg = dwlsmv_admissadvice(A, ind, fac)
    printf("{txt} 12e. singular (I - B)          : flag %g, rows 1-4 flags missing %g; the note reads:\n%s",
           A[5, 1], missing(A[1::4, 1]), msg)
    if (A[5, 1] != 1 | missing(A[1::4, 1]) != 4 | A[6, 1] != 0) {
        FAIL = FAIL \ "admissibility: a singular (I - B) is not reported, or rows 1-4 were not left missing"
    }
    if (!strpos(msg, "feedback loop") | strpos(msg, "Heywood")) {
        FAIL = FAIL \ "admissibility: the singular-(I - B) note is wrong"
    }

    return(FAIL)
}

// Negate the one free parameter at (ptype, row, col); psi matches either
// orientation. Reflections in check 13 are spelled out with this, parameter
// by parameter, rather than generated by the rule under test.
real colvector tf_neg(real colvector th, struct dwlsmv_parmap scalar par,
                      real scalar t, real scalar r, real scalar c)
{
    real colvector o
    real scalar i

    o = th
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.ptype[i] == t &
            ((par.row[i] == r & par.col[i] == c) |
             (t == 2 & par.row[i] == c & par.col[i] == r))) {
            o[par.fidx[i]] = -o[par.fidx[i]]
            return(o)
        }
    }
    _error("tf_neg: no such free parameter")
}

// ---- 13. factor signs are canonicalized after convergence -------------------
string colvector tf_sign(struct dwlsmv_stage1 scalar st)
{
    struct dwlsmv_parmap scalar pc, ps, pm, pn
    struct dwlsmv_fit    scalar fc, fs, fm, fn
    string colvector FAIL
    string scalar spec
    real colvector w, thr, thu, thn
    real scalar j, fr, f0, cu, cn, fu, iu, fnn, itn, a4

    FAIL = J(0, 1, "")
    w    = 1 :/ diagonal(st.Gamma)

    // 13a. correlated, unit variance: reverse f2 by hand -- its three
    // loadings and psi(1,2). Same F; canonicalizes back BIT FOR BIT, since
    // negation is exact.
    pc  = dwls_identify(tf_canon("cor_unitvar"), st.map)
    fc  = dwls_fit(st, pc)
    thr = fc.theta
    for (j = 4; j <= 6; j++) thr = tf_neg(thr, pc, 1, j, 2)
    thr = tf_neg(thr, pc, 2, 1, 2)
    fr  = tf_F(thr, st.s, w, pc, st.map)
    f0  = tf_F(fc.theta, st.s, w, pc, st.map)
    printf("{txt} 13a. cor_unitvar, f2 reversed : F rel change %8.2e ; canonicalizes back exactly %g ; idempotent %g\n",
           abs(fr - f0) / f0, dwlsmv_signcanon(thr, pc) == fc.theta,
           dwlsmv_signcanon(fc.theta, pc) == fc.theta)
    if (abs(fr - f0) / f0 > 1e-15) {
        FAIL = FAIL \ "signs: reversing f2 by hand changed F, so the hand reflection is wrong"
    }
    if (dwlsmv_signcanon(thr, pc) != fc.theta | dwlsmv_signcanon(fc.theta, pc) != fc.theta) {
        FAIL = FAIL \ "signs: cor_unitvar with f2 reversed does not canonicalize back to theta-hat"
    }

    // 13b. F2 <- F1, unit variance: reversing f1 negates its loadings AND
    // beta; reversing both factors leaves beta alone (D B D).
    ps  = dwls_identify(tf_canon("str_unitvar"), st.map)
    fs  = dwls_fit(st, ps)
    thr = fs.theta
    for (j = 1; j <= 3; j++) thr = tf_neg(thr, ps, 1, j, 1)
    thr = tf_neg(thr, ps, 3, 2, 1)
    cu  = (dwlsmv_signcanon(thr, ps) == fs.theta)
    thr = fs.theta
    for (j = 1; j <= 6; j++) thr = tf_neg(thr, ps, 1, j, (j <= 3 ? 1 : 2))
    cn  = (dwlsmv_signcanon(thr, ps) == fs.theta)
    printf("{txt} 13b. str_unitvar              : f1 reversed back exactly %g ; both reversed, beta untouched, back exactly %g\n",
           cu, cn)
    if (!cu | !cn) {
        FAIL = FAIL \ "signs: str_unitvar reflections do not canonicalize back to theta-hat"
    }

    // 13c. marker metric: y4's fixed lambda = 1 identifies f2's sign, so
    // negating f2's FREE loadings is a different model, not a reflection, and
    // must be left exactly as given.
    pm  = dwls_identify(tf_canon("cor_marker"), st.map)
    fm  = dwls_fit(st, pm)
    thr = tf_neg(tf_neg(tf_neg(fm.theta, pm, 1, 5, 2), pm, 1, 6, 2), pm, 2, 1, 2)
    printf("{txt} 13c. cor_marker, f2 negated   : left unchanged %g (a fixed loading identifies the sign)\n",
           dwlsmv_signcanon(thr, pm) == thr)
    if (dwlsmv_signcanon(thr, pm) != thr) {
        FAIL = FAIL \ "signs: a marker factor was reflected, though its sign is identified"
    }

    // 13d. END TO END. Start f2's loadings at -0.7: the unconstrained
    // optimizer converges to the REFLECTED solution -- asserted, so this
    // check cannot pass by never exercising the rule -- and dwls_fit()
    // returns the canonical one.
    spec = tf_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                    "lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
                    "lambda f2 y4 free -0.7", "lambda f2 y5 free -0.7",
                    "lambda f2 y6 free -0.7", "psi f1 f2 free",
                    "metric unitvariance f1", "metric unitvariance f2"))
    pn = dwls_identify(spec, st.map)
    (void) dwlsmv_minimize(st.s, J(st.map.D, 1, 1), pn, st.map, tf_start(pn),
                           "d2", "nr", 300, thu, fu, iu)
    (void) dwlsmv_minimize(st.s, w, pn, st.map, thu, "d2", "nr", 300, thn, fnn, itn)
    a4 = thn[pn.fidx[selectindex(pn.ptype :== 1 :& pn.row :== 4)]]
    fn = dwls_fit(st, pn)
    printf("{txt} 13d. f2 started at -0.7       : raw lambda_y4 %9.6f ; dwls_fit lambda_y4 %9.6f ; max |theta - canonical fit| %8.2e\n",
           a4, fn.theta[pn.fidx[selectindex(pn.ptype :== 1 :& pn.row :== 4)]],
           max(abs(fn.theta - fc.theta)))
    if (!(a4 < 0)) {
        FAIL = FAIL \ "signs: the negative start did not reach the reflected solution, so 13d tests nothing"
    }
    if (fn.converged != 1 | dwlsmv_gapratio(fn.theta, fc.theta, 1e-8, 0) > 1) {
        FAIL = FAIL \ "signs: dwls_fit from a reflected start did not return the canonical theta-hat"
    }
    if (fn.fmin != tf_F(fn.theta, st.s, w, pn, st.map)) {
        FAIL = FAIL \ "signs: fmin is not F at the canonicalized theta"
    }

    return(FAIL)
}

// Set the one free parameter at (ptype, row, col) to v.
real colvector tf_set(real colvector th, struct dwlsmv_parmap scalar par,
                      real scalar t, real scalar r, real scalar c, real scalar v)
{
    real colvector o
    real scalar i

    o = th
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.ptype[i] == t & par.row[i] == r & par.col[i] == c) {
            o[par.fidx[i]] = v
            return(o)
        }
    }
    _error("tf_set: no such free parameter")
}

// ---- 13e-h. the sign-canonicalization cases check 13 did not reach ----------
// Added 2026-09-23 after the Module 7 audit. Each is at a stated theta; none
// needs data.
string colvector tf_sign2(struct dwlsmv_statmap scalar mapE)
{
    struct dwlsmv_parmap  scalar pf, pz, p3, pc, p1
    struct dwlsmv_statmap scalar map3, mapD
    string colvector FAIL
    real colvector th, th0, thr, lam, ok2
    real scalar j, ok, dsig

    FAIL = J(0, 1, "")
    ok2  = (-0.5 \ 0.5)

    // 13e. A FIXED NONZERO psi_12 touches both factors, so neither is
    // reflectable: f2 with negative loadings is left exactly as given. The
    // same table with psi_12 fixed at ZERO is the control -- a fixed zero
    // does not identify a sign, so there f2 is reflected.
    lam = (0.7 \ 0.7 \ 0.7 \ -0.7 \ -0.7 \ -0.7)
    pf  = dwls_identify(tf_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                                 "lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
                                 "lambda f2 y4 free", "lambda f2 y5 free", "lambda f2 y6 free",
                                 "psi f1 f2 fixed 0.3",
                                 "metric unitvariance f1", "metric unitvariance f2")), mapE)
    pz  = dwls_identify(tf_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                                 "lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
                                 "lambda f2 y4 free", "lambda f2 y5 free", "lambda f2 y6 free",
                                 "psi f1 f2 fixed 0",
                                 "metric unitvariance f1", "metric unitvariance f2")), mapE)
    th  = tf_theta(pf, lam, 0, ok2)
    thr = tf_theta(pz, lam, 0, ok2)
    ok  = (dwlsmv_signcanon(th, pf) == th) &
          (dwlsmv_signcanon(thr, pz) == tf_theta(pz, abs(lam), 0, ok2))
    printf("{txt} 13e. fixed psi_12 = 0.3 / = 0   : left alone at 0.3, reflected at 0 -- %g\n", ok)
    if (!ok) {
        FAIL = FAIL \ "signs: a fixed nonzero psi_12 did not block the reflection, or a fixed zero did"
    }

    // 13f. THREE FACTORS: f1 <-> f2 correlated (psi_21 free), f3 <- f1 and
    // f3 <- f2. Reverse f2 AND f3 by hand: their six loadings, psi_21 (one
    // end reversed) and beta_31 (one end reversed) change sign; beta_32 joins
    // two reversed factors and must NOT. sigma must be unchanged -- which
    // checks the hand reflection itself -- and canonicalization must return
    // theta bit for bit.
    map3 = dwlsmv_statmap_build(J(9, 1, 3))
    p3   = dwls_identify(tf_spec(("factor f1 f2 f3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
                                  "lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
                                  "lambda f2 y4 free", "lambda f2 y5 free", "lambda f2 y6 free",
                                  "lambda f3 y7 free", "lambda f3 y8 free", "lambda f3 y9 free",
                                  "psi f1 f2 free", "beta f3 f1 free", "beta f3 f2 free",
                                  "metric unitvariance f1", "metric unitvariance f2",
                                  "metric unitvariance f3")), map3)
    th0 = J(p3.npar, 1, .)
    for (j = 1; j <= p3.nall; j++) {
        if (!p3.free[j]) continue
        if (p3.ptype[j] == 4) th0[p3.fidx[j]] = (p3.col[j] == 1 ? -0.5 : 0.5)
        if (p3.ptype[j] == 1) th0[p3.fidx[j]] = 0.6 + 0.02 * p3.row[j]
    }
    th0 = tf_set(th0, p3, 2, 2, 1, 0.3)
    th0 = tf_set(th0, p3, 3, 3, 1, 0.2)
    th0 = tf_set(th0, p3, 3, 3, 2, 0.4)
    thr = th0
    for (j = 4; j <= 9; j++) thr = tf_neg(thr, p3, 1, j, (j <= 6 ? 2 : 3))
    thr  = tf_neg(tf_neg(thr, p3, 2, 2, 1), p3, 3, 3, 1)
    dsig = max(abs(dwls_model(thr, p3, map3) - dwls_model(th0, p3, map3)))
    ok   = (dwlsmv_signcanon(thr, p3) == th0)
    printf("{txt} 13f. three factors, f2+f3 reversed: max |d sigma| %g ; back exactly %g (beta_32 untouched)\n",
           dsig, ok)
    if (dsig > 1e-15) {
        FAIL = FAIL \ "signs: the hand reflection of f2 and f3 changed sigma, so the test's reflection is wrong"
    }
    if (!ok) {
        FAIL = FAIL \ "signs: three-factor reflection of f2 and f3 does not canonicalize back"
    }

    // 13g. An ANCHOR OF EXACTLY ZERO leaves the factor alone, even with its
    // other loadings negative.
    pc  = dwls_identify(tf_canon("cor_unitvar"), mapE)
    th  = tf_theta(pc, (0.7 \ 0.7 \ 0.7 \ 0 \ -0.7 \ -0.7), 0.3, ok2)
    ok  = (dwlsmv_signcanon(th, pc) == th)
    printf("{txt} 13g. f2 anchor exactly 0       : left alone %g\n", ok)
    if (!ok) {
        FAIL = FAIL \ "signs: a factor whose anchor loading is exactly 0 was reflected"
    }

    // 13h. ONE FACTOR, unit variance (case D's map): all loadings negative ->
    // all positive; thresholds untouched.
    mapD = dwlsmv_statmap_build((4 \ 4 \ 3 \ 3))
    p1   = dwls_identify(tf_spec(("factor f", "ind y1 y2 y3 y4",
                                  "lambda f y1 free", "lambda f y2 free",
                                  "lambda f y3 free", "lambda f y4 free",
                                  "metric unitvariance f")), mapD)
    th0 = J(p1.npar, 1, .)
    for (j = 1; j <= p1.nall; j++) {
        if (!p1.free[j]) continue
        if (p1.ptype[j] == 4) th0[p1.fidx[j]] = -1 + 0.5 * p1.col[j]
        if (p1.ptype[j] == 1) th0[p1.fidx[j]] = 0.5 + 0.1 * p1.row[j]
    }
    thr = th0
    for (j = 1; j <= 4; j++) thr = tf_neg(thr, p1, 1, j, 1)
    ok  = (dwlsmv_signcanon(thr, p1) == th0)
    printf("{txt} 13h. one factor, all reversed  : back exactly %g\n", ok)
    if (!ok) {
        FAIL = FAIL \ "signs: a one-factor unit-variance reflection does not canonicalize back"
    }

    return(FAIL)
}

real scalar tf_checks(real matrix YD, real matrix YE,
                      string colvector DKEY, real colvector DVAL,
                      string colvector WKEY, real colvector WVAL,
                      string colvector PFIT, string colvector PKEY,
                      real colvector PVAL,
                      string colvector EFIT, real colvector ECHI)
{
    struct dwlsmv_stage1  scalar stD, stE, st, st_copy
    struct dwlsmv_parmap  scalar par, par_copy
    struct dwlsmv_fit     scalar fit
    string colvector FAIL
    string rowvector fits, fac
    string scalar fitn, lab
    real colvector w, th0, thu, th, r, sig, ga, gn, tp, tm
    real scalar k, i, d, h, v, worst, ratio, cu, fu, iu, cv, fm, it
    real scalar chi, B, Tn1, Tn, gap_tau, g0, g1, N

    FAIL = J(0, 1, "")
    stD  = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    stE  = dwls_sample_stats(YE, J(6, 1, 3))
    fits = ("D_marker", "cor_unitvar", "str_unitvar", "cor_marker", "str_marker")

    // ---- 1. the gradient, at the starting values ---------------------------
    // Budget: the analytic side carries Delta's finite-difference noise,
    // ~1e-10 per entry (dwlsmv_delta_cols' header), multiplied by
    // 2 sum w |r| -- O(10) at the starts; the difference quotient of F carries
    // O(h^2) truncation and O(eps F / h) cancellation, ~1e-10 each at
    // h = 1e-5. So ~1e-9 on a gradient whose entries are O(0.1-1). 1e-7
    // absolute is two orders above that and six below a sign or factor-of-2
    // error, which is what this check exists to see.
    par = dwls_identify(tf_canon("str_unitvar"), stE.map)
    w   = 1 :/ diagonal(stE.Gamma)
    th0 = tf_start(par)
    ga  = tf_grad(th0, stE.s, w, par, stE.map)'
    gn  = J(par.npar, 1, .)
    for (k = 1; k <= par.npar; k++) {
        h     = 1e-5 * max((1, abs(th0[k])))
        tp    = th0
        tm    = th0
        tp[k] = th0[k] + h
        tm[k] = th0[k] - h
        gn[k] = (tf_F(tp, stE.s, w, par, stE.map) -
                 tf_F(tm, stE.s, w, par, stE.map)) / (2 * h)
    }
    ratio = dwlsmv_gapratio(ga, gn, 1e-7, 0)
    printf("{txt}  1. gradient vs dF/dtheta      : worst ratio %8.2e of 1e-7 over %f entries, max |g| %g\n",
           ratio, par.npar, max(abs(ga)))
    if (ratio > 1) {
        FAIL = FAIL \ ("the analytic gradient disagrees with a difference of F, ratio " +
                       strofreal(ratio))
    }

    for (d = 1; d <= cols(fits); d++) {
        fitn = fits[d]
        if (fitn == "D_marker") {
            st  = stD
            fac = ("f")
        }
        else {
            st  = stE
            fac = ("f1", "f2")
        }
        N   = st.N
        par = dwls_identify(tf_canon(fitn), st.map)
        w   = 1 :/ diagonal(st.Gamma)
        th0 = tf_start(par)

        // ---- 2. convergence from ULS starts --------------------------------
        cu = dwlsmv_minimize(st.s, J(st.map.D, 1, 1), par, st.map, th0,
                             "d2", "nr", 300, thu, fu, iu)

        st_copy  = st
        par_copy = par
        fit      = dwls_fit(st, par)
        th       = fit.theta

        printf("{txt}  2. %-11s converges       : ULS %g in %2.0f it, DWLS %g in %2.0f it, F-hat %21.15g\n",
               fitn, cu, iu, fit.converged, fit.iterations, fit.fmin)
        if (cu != 1) {
            FAIL = FAIL \ (fitn + ": ULS did not converge from the table's start values")
        }
        if (fit.converged != 1) {
            FAIL = FAIL \ (fitn + ": DWLS did not converge from the ULS solution")
        }

        // ---- 3. theta-hat against lavaan -----------------------------------
        worst = 0
        for (i = 1; i <= par.nall; i++) {
            if (!par.free[i]) continue
            lab = tf_lab(par, i, fac)
            if (fitn == "D_marker") {
                v = tf_get(DKEY, DKEY, DVAL, "", lab)
            }
            else {
                v = tf_get(PFIT, PKEY, PVAL, fitn, lab)
            }
            if (missing(v)) {
                FAIL = FAIL \ (fitn + ": lavaan has no estimate labelled " + lab)
                continue
            }
            worst = max((worst, dwlsmv_gapratio(th[par.fidx[i]], v, 1e-5, 1e-4)))
        }
        printf("{txt}  3. %-11s theta vs lavaan : worst ratio %8.2e of the theta-hat budget over %f parameters\n",
               fitn, worst, par.npar)
        if (worst > 1) {
            FAIL = FAIL \ (fitn + ": theta-hat disagrees with lavaan, ratio " + strofreal(worst))
        }

        sig = dwls_model(th, par, st.map)
        r   = st.s - sig

        // ---- 4. sigma(theta-hat) against lavaan's wls.est, case D ----------
        if (fitn == "D_marker") {
            worst = 0
            for (i = 1; i <= st.map.D; i++) {
                v = tf_get(WKEY, WKEY, WVAL, "", st.map.label[i])
                if (missing(v)) {
                    FAIL = FAIL \ ("case D: lavaan's wls.est has no entry " + st.map.label[i])
                    continue
                }
                worst = max((worst, dwlsmv_gapratio(sig[i], v, 1e-5, 1e-4)))
            }
            printf("{txt}  4. %-11s sigma vs wls.est: worst ratio %8.2e of the theta-hat budget over %f entries\n",
                   fitn, worst, st.map.D)
            if (worst > 1) {
                FAIL = FAIL \ ("case D: sigma(theta-hat) disagrees with lavaan's wls.est")
            }
        }

        // ---- 5. the weight convention, against lavaan's chi-square ---------
        if (fitn != "D_marker") {
            chi = tf_get(EFIT, EFIT, ECHI, "", fitn)
            B   = 2 * quadcross(w, abs(r) :* (1e-8 :+ 1e-8 :* abs(st.s))) +
                  1e-8 * fit.fmin
            Tn1 = fit.T_raw
            Tn  = N * fit.fmin
            printf("{txt}  5. %-11s lavaan chisq    : %16.12f ; T_raw %16.12f ratio %8.2e ; Mplus-style N F %16.12f ratio %8.2e\n",
                   fitn, chi, Tn1, abs(Tn1 - chi) / ((N - 1) * B),
                   Tn, abs(Tn - chi) / ((N - 1) * B))
            if (missing(chi)) {
                FAIL = FAIL \ (fitn + ": the equivalence fixture has no chisq")
            }
            else {
                if (abs(Tn1 - chi) > (N - 1) * B) {
                    FAIL = FAIL \ (fitn + ": T_raw is not lavaan's chisq -- W_D or the (N-1) scaling moved")
                }
                if (abs(Tn - chi) <= (N - 1) * B) {
                    FAIL = FAIL \ (fitn + ": N F-hat also matches lavaan's chisq, so this check cannot tell N from N-1")
                }
            }
        }

        // ---- 6. tau-hat = s_tau --------------------------------------------
        // Measured 0.0 on every fit, and exact in the mathematics: a
        // Gauss-Newton step solves the linear threshold block exactly. Asserted
        // at 1e-15 (~4 ulps at |tau| < 2), leaving room for the last addition's
        // rounding rather than demanding bitwise equality.
        // The residual s_tau - sigma_tau(theta-hat) AND the parameters
        // tau-hat - s_tau: a change to sigma_tau = tau keeps the residual at 0
        // while moving tau-hat (added 2026-10-07). Free tau rows are table
        // rows 1..T in s order (Section 4.2).
        gap_tau = max(abs(r[1::st.map.T]))
        for (i = 1; i <= st.map.T; i++) {
            if (par.ptype[i] != 4 | par.free[i] != 1) {
                FAIL = FAIL \ (fitn + ": table row " + strofreal(i) + " is not a free threshold")
                continue
            }
            gap_tau = max((gap_tau, abs(th[par.fidx[i]] - st.s[i])))
        }
        if (d == 1 | gap_tau != 0) {
            printf("{txt}  6. %-11s tau-hat = s_tau : max |s_tau - tau-hat| = %g\n", fitn, gap_tau)
        }
        if (gap_tau > 1e-15) {
            FAIL = FAIL \ (fitn + ": tau-hat is not s_tau, gap " + strofreal(gap_tau))
        }

        // ---- 7. stationarity -----------------------------------------------
        g0 = max(abs(tf_grad(th0, st.s, w, par, st.map)))
        g1 = max(abs(tf_grad(th,  st.s, w, par, st.map)))
        printf("{txt}  7. %-11s stationarity    : max|g(theta-hat)| %8.2e = %8.2e of max|g(start)|\n",
               fitn, g1, g1 / g0)
        if (g1 / g0 > 1e-8) {
            FAIL = FAIL \ (fitn + ": the gradient at theta-hat is not small, relative " +
                           strofreal(g1 / g0))
        }

        // ---- 9. the inputs are not written ---------------------------------
        if (st.s != st_copy.s | st.Gamma != st_copy.Gamma | st.N != st_copy.N |
            st.map.label != st_copy.map.label |
            par.value != par_copy.value | par.free != par_copy.free |
            par.fidx != par_copy.fidx | par.ptype != par_copy.ptype |
            par.row != par_copy.row | par.col != par_copy.col |
            par.npar != par_copy.npar) {
            FAIL = FAIL \ (fitn + ": dwls_fit() wrote to st or par")
        }

        // ---- 10. bookkeeping -----------------------------------------------
        if (fit.D_stat != st.map.D - par.npar) {
            FAIL = FAIL \ (fitn + ": D_stat is not D - npar")
        }
        if (rows(fit.Delta) != st.map.D | cols(fit.Delta) != par.npar) {
            FAIL = FAIL \ (fitn + ": Delta is not D x npar")
        }
        else if (fit.Delta != dwls_delta(th, par, st.map)) {
            FAIL = FAIL \ (fitn + ": Delta is not dwls_delta() at theta-hat")
        }
        if (fit.T_raw != (N - 1) * fit.fmin) {
            FAIL = FAIL \ (fitn + ": T_raw is not (N - 1) * fmin")
        }

        // ---- 12a. every fitted solution here is admissible -----------------
        if (dwlsmv_admissible(th, par, st.map)[., 1] != (0 \ 0 \ 0 \ 0 \ 0 \ 0) |
            dwlsmv_admissadvice(dwlsmv_admissible(th, par, st.map),
                                J(6, 1, "y"), J(2, 1, "f")) != "") {
            FAIL = FAIL \ (fitn + ": the fitted solution is flagged inadmissible")
        }
        if (fit.fmin != tf_F(th, st.s, w, par, st.map)) {
            FAIL = FAIL \ (fitn + ": fmin is not F_DWLS at theta-hat")
        }
        if (rows(fit.acov) != 0 | !missing(fit.T3) | !missing(fit.T3_df) |
            !missing(fit.shift_a) | !missing(fit.shift_b)) {
            FAIL = FAIL \ (fitn + ": a Module 8-9 field was filled by Module 7")
        }
    }
    printf("{txt}  9. inputs unchanged           : st and par bitwise identical after all five fits\n")
    printf("{txt} 10. bookkeeping                : D_stat, Delta, T_raw, fmin, and Modules 8-9 left alone\n")
    printf("{txt} 12a. admissibility, fitted     : all six conditions pass on all five fits, and the note is empty\n")
    FAIL = FAIL \ tf_admiss(stE.map)
    FAIL = FAIL \ tf_sign(stE)
    FAIL = FAIL \ tf_sign2(stE.map)

    // ---- 8. non-convergence is a return value ------------------------------
    par = dwls_identify(tf_canon("str_unitvar"), stE.map)
    w   = 1 :/ diagonal(stE.Gamma)
    th0 = tf_start(par)
    // No assertion on the iteration COUNT: capped at 1, optimize() reports 2,
    // and [M-5] optimize() says only "the number of iterations used". What
    // matters is that the result is the last iterate and not the start: its F
    // is below the starting F, and it is F at the theta returned beside it.
    cv  = dwlsmv_minimize(stE.s, w, par, stE.map, th0, "d2", "nr", 1, th, fm, it)
    g0  = tf_F(th0, stE.s, w, par, stE.map)
    printf("{txt}  8. capped at one iteration    : converged %g, iterations reported %g, F %g from %g\n",
           cv, it, fm, g0)
    if (cv != 0) {
        FAIL = FAIL \ "one iteration reported convergence"
    }
    if (missing(th) | rows(th) != par.npar | missing(fm)) {
        FAIL = FAIL \ "the unconverged result carries no theta or no F"
    }
    else if (!(fm < g0) | fm != tf_F(th, stE.s, w, par, stE.map)) {
        FAIL = FAIL \ "the unconverged result is not the last iterate and its own F"
    }

    printf("{txt}test_fit: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// ---- 11. the error paths ----------------------------------------------------
// Each builds its own case E inputs, so no struct has to survive at the
// interactive level.
void tf_bad(real matrix YE, real scalar which)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    real colvector th
    real scalar fm, it, cv

    st  = dwls_sample_stats(YE, J(6, 1, 3))
    par = dwls_identify(tf_canon("str_unitvar"), st.map)

    if (which == 1) {
        st.s = st.s[1::(st.map.D - 1)]
        fit  = dwls_fit(st, par)
    }
    else if (which == 2) {
        st.Gamma[3, 3] = 0
        fit = dwls_fit(st, par)
    }
    else if (which == 3) {
        par.npar = par.npar + 1
        fit = dwls_fit(st, par)
    }
    else if (which == 4) {
        cv = dwlsmv_minimize(st.s, J(st.map.D, 1, 1), par, st.map,
                             J(par.npar - 1, 1, 0.5), "d2", "nr", 300, th, fm, it)
    }
    else if (which == 5) {
        cv = dwlsmv_minimize(st.s, J(st.map.D - 1, 1, 1), par, st.map,
                             tf_start(par), "d2", "nr", 300, th, fm, it)
    }
    else if (which == 6) {
        th = dwlsmv_signcanon(J(par.npar + 1, 1, 0.5), par)
    }
    else {
        (void) dwlsmv_admissadvice(J(5, 4, 0), J(6, 1, "y"), J(2, 1, "f"))
    }
}

// ---- 14. the rank check's second point is the fit's start (2026-10-04) ------
// A random intercept G (loadings fixed 1, psi_GG free, orthogonal) with
// psi_21 free loses rank at the symmetric starts only (spec 4.8, derivation
// Finding 5). dwls_identify() passes it at theta2 and makes theta2 the start;
// from theta0 the fit did not converge (psi_21 -> -25.5). Asserted: the start
// moved, a model full rank at theta0 keeps theta0, the default-start fit
// converges admissibly, and it agrees with a fit from distinct user starts.
// The 1e-5 agreement tolerance was chosen after the probe showed the two fits
// equal to the four decimals it printed.
string scalar tf_ri_spec(string scalar lamF)
{
    string scalar s, nl
    nl = char(10)
    s = "factor F1 F2 G" + nl + "ind y1 y2 y3 y4 y5 y6" + nl
    s = s + "lambda F1 y1 free" + lamF + nl + "lambda F1 y2 free" + nl + "lambda F1 y3 free" + nl
    s = s + "lambda F2 y4 free" + lamF + nl + "lambda F2 y5 free" + nl + "lambda F2 y6 free" + nl
    s = s + "metric marker G y1" + nl
    s = s + "lambda G y2 fixed 1" + nl + "lambda G y3 fixed 1" + nl + "lambda G y4 fixed 1" + nl
    s = s + "lambda G y5 fixed 1" + nl + "lambda G y6 fixed 1" + nl
    s = s + "metric unitvariance F1" + nl + "metric unitvariance F2" + nl
    s = s + "psi F2 F1 free" + nl + "psi G G free" + nl
    return(s)
}

real scalar tf_second(real matrix YE)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, pard, pc
    struct dwlsmv_fit    scalar fit, fitd
    real matrix A
    real scalar i, p21, dev
    string colvector FAIL

    FAIL = J(0, 1, "")
    st   = dwls_sample_stats(YE, J(6, 1, 3))
    par  = dwls_identify(tf_ri_spec(""), st.map)
    pard = dwls_identify(tf_ri_spec(" .9"), st.map)
    pc   = dwls_identify(tf_canon("str_unitvar"), st.map)

    p21 = select(par.value, par.name :== "psi[F2,F1]")
    printf("{txt} 14. second point: psi[F2,F1] start %9.4f; control loadings start at 0.7: %f\n",
           p21, all(select(pc.value, pc.ptype :== 1 :& pc.free) :== 0.7))
    if (p21 == 0) {
        FAIL = FAIL \ "14: the random-intercept model kept psi[F2,F1]'s start of 0"
    }
    if (!all(select(pc.value, pc.ptype :== 1 :& pc.free) :== 0.7)) {
        FAIL = FAIL \ "14: a model of full rank at the starts did not keep them"
    }

    fit  = dwls_fit(st, par)
    fitd = dwls_fit(st, pard)
    A    = dwlsmv_admissible(fit.theta, fit.par, fit.stage1[1].map)
    dev  = max(abs(fit.theta - fitd.theta))
    printf("{txt} 14. default starts: converged %f, inadmissible flags %f; max |theta - distinct-start theta| %9.2e\n",
           fit.converged, sum(A[., 1]), dev)
    if (fit.converged != 1) {
        FAIL = FAIL \ "14: the random-intercept fit from the default starts did not converge"
    }
    if (sum(A[., 1]) != 0) {
        FAIL = FAIL \ "14: the random-intercept fit from the default starts is inadmissible"
    }
    if (fitd.converged != 1 | dev > 1e-5) {
        FAIL = FAIL \ "14: default-start and distinct-start fits disagree"
    }

    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(tf_checks(YD, YE, DKEY, DVAL, WKEY, WVAL, ///
                                            PFIT, PKEY, PVAL, EFIT, ECHI)))
mata: st_local("nfail", strofreal(`nfail' + tf_second(YE)))

local why1 "dwls_fit: s shorter than the map's D"
local why2 "dwls_fit: a zero on Gamma's diagonal"
local why3 "dwls_fit: npar larger than the free rows fidx numbers"
local why4 "dwlsmv_minimize: start not npar rows long"
local why5 "dwlsmv_minimize: w not D rows long"
local why6 "dwlsmv_signcanon: theta not npar rows long"
local why7 "dwlsmv_admissadvice: A not 6 x 4 or 7 x 4"
local nbadfail 0
forvalues i = 1/7 {
    capture noisily mata: tf_bad(YE, `i')
    local rc = _rc
    if `rc' == 3498 {
        display as text " 11. reject `i' (rc 3498): `why`i''"
    }
    else {
        display as error " 11. reject `i' FAILED with rc `rc': `why`i''"
        local nbadfail = `nbadfail' + 1
    }
}

if `nfail' | `nbadfail' {
    display as error _n "test_fit: `nfail' check failure(s), `nbadfail' error-path failure(s)"
    exit 9
}

display as result "test_fit: dwls_fit converges from ULS starts on cases D and E, matches lavaan's theta-hat, reproduces lavaan's chisq as (N-1) F-hat, puts tau-hat on s_tau, reports non-convergence without aborting, and leaves its inputs untouched"
exit 0
