*! Module 6: dwls_delta() -- the Jacobian d sigma(theta) / d theta'.
*!
*! Executed under StataMP 18 on 2026-09-22 and it passes.
*!
*! WHAT THIS FILE IS FOR. Delta is numerical (central differences against
*! Module 5), which is what spec Section 5 prescribes and what any later
*! analytic version must be validated against. "Numerical" must not mean
*! "unchecked", so several of the eight checks below use oracles that do not
*! involve differentiating anything:
*!
*!   1   THE STRUCTURAL BLOCKS, exactly: 0.0, not a tolerance. Gate 2's
*!       Finding 3.
*!
*!       BE HONEST ABOUT WHAT THIS ONE CANNOT SEE. Both halves are partly
*!       tautological against *this* implementation. The threshold column is
*!       compared to a unit vector positioned by `dwlsmv_pos_tau()` -- the same
*!       call, on the same arguments, that `dwlsmv_delta_cols()` used to place
*!       the 1 -- so a `dwls_identify()` that emitted the tau rows in scrambled
*!       order would produce a permuted leading block and this would still say
*!       0.0. And `dwlsmv_delta_cols()` only ever writes rows (T+1)..D for a
*!       structural column, so the "zero in every threshold row" half cannot
*!       fire here at all.
*!
*!       It is not worthless: it would catch a Delta that differenced the whole
*!       of sigma rather than only the correlation block, which is the obvious
*!       way to write a numerical Jacobian and would put noise in both blocks.
*!       But THE PERMUTATION DETECTOR IS ELSEWHERE -- test_identify.do check 2
*!       asserts row-by-row that the leading T rows are the thresholds in
*!       index-map order with fidx 1..T, and check 4 below compares against a
*!       LABELLED oracle, which no reordering survives. Do not read check 1 as
*!       covering either.
*!   2   A CLOSED FORM. On a single-factor unit-variance model R_jk =
*!       lambda_j lambda_k, so d R_jk / d lambda_m is lambda_k or lambda_j or
*!       zero. The whole lambda block of Delta has an algebraic oracle.
*!   3   ANOTHER ONE, on the two-factor structural model, for the beta column:
*!       Phi = [[1, beta], [beta, beta^2 + 1]], so d R_jk / d beta is 0 within
*!       F1, lambda_j lambda_k across, and 2 beta lambda_j lambda_k within F2.
*!       This is the column that would be wrong if (I-B)^-1 Psi (I-B)^-T were
*!       mis-assembled, and it is the one case A-D could not pose.
*!   4   AGAINST LAVAAN'S ANALYTIC JACOBIAN, on four configurations.
*!       lavInspect(fit, "delta") is the Jacobian -- verified by dispatch and
*!       against a numerical Jacobian when the oracle was generated, not
*!       assumed; the belief carried into that session was that it was the
*!       delta-parameterization SCALING, and that belief was wrong. The
*!       scaling is lavInspect(fit, "mm.delta"), a p x 1 vector of ones here.
*!   5   THE RANK CHECK'S RAW MATERIAL. rank(Delta) equals the free-parameter
*!       count on all four identified configurations, and the singular values
*!       are reported alongside the oracle's.
*!   6   THE OVERPARAMETERIZED CASE, at a beta Module 3 would never pick.
*!       Section 4.8's psi_12-alongside-beta error, evaluated at the oracle's
*!       stated beta = 0.4 rather than at Module 3's default start of 0. The
*!       point is to separate the two candidate knife edges: it is psi_12 = 0
*!       that makes the collinearity numerically exact, NOT beta = 0. Measured
*!       here, the singular-value ratio is 2.1e-17 at beta = 0.4, against
*!       6.0e-17 at the defaults -- both machine zero. Had beta been the
*!       sensitive one, this check would have shown the floor rise. Check 7
*!       runs the same table at Module 3's own default starts so that both
*!       halves of that comparison are measured here rather than one of them
*!       being quoted from a log.
*!   8   THE SAME TABLE AT A NONZERO psi_12 (0.1 and 0.3), added 2026-09-23.
*!       Still deficiency 1 -- the model is unidentified everywhere, and the
*!       rank check sees it everywhere -- but the null direction rotates off
*!       (beta, psi_12) into a rescaling of F2's loadings, asserted against
*!       the analytic vector of the derivation's Gate 2 Finding 4. This
*!       replaced an untested claim, carried in the spec and the project's notes until
*!       then, that a nonzero psi_12 gives full rank.
*!
*! WHY CHECK 6 BUILDS ITS PARAMETER TABLE BY HAND. dwls_identify() refuses to
*! return a rank-deficient table -- that is Module 3's whole job, and
*! test_identify.do asserts the refusal. So this file appends the offending
*! psi row to a valid table directly. The row is appended rather than inserted
*! in (row, col) order, which is fine because nothing here reads table order:
*! Delta's columns are keyed by fidx and the null direction is read back
*! through the same index.
*!
*! TOLERANCES. Checks 1-3 compare against exact algebra, but through a central
*! difference, so they get the finite-difference budget rather than a machine
*! -precision one: 1e-8 absolute and relative, which is spec Section 7.2's
*! bivariate-correlation tier and is roughly two orders above the ~1e-10 noise
*! floor of the step size dwlsmv_delta_cols() documents. Check 4 compares a
*! numerical Jacobian against an analytic one and gets the same tier, for the
*! reason Section 7.2 permits: the two implementations use genuinely different
*! numerical approximations at that step, stated here rather than loosened
*! silently. Every measured value is printed.
*!
*! Two tolerances here are NOT from Section 7.2, and both are named rather than
*! buried. The null-direction checks in 6 use 1e-3, where the generator asserts
*! 1e-8 on its own analytic null vector: a singular vector's accuracy degrades
*! like the noise divided by the GAP to the neighbouring singular value, so a
*! finite-difference Delta cannot be held to the analytic figure, and no
*! Section 7.2 tier describes an eigenvector. The measured components are
*! printed, and the largest off-direction entry is printed too, so the real
*! headroom is visible rather than implied. The singular-value RATIO is asserted
*! separately at 1e-12, which is the quantity the rank check actually depends
*! on and where the margin matters.
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
             src/mata/dwls_identify.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_delta: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/expected/lavaan_modelimplied_input.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: IFIT = st_sdata(., "fit")
mata: IKEY = st_sdata(., "element")
mata: IVAL = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/lavaan_modelimplied_jacobian.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: JFIT = st_sdata(., "fit")
mata: JROW = st_sdata(., "row")
mata: JCOL = st_sdata(., "col")
mata: JVAL = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/lavaan_modelimplied_rank.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: RFIT = st_sdata(., "fit")
mata: RRNK = strtoreal(st_sdata(., "rank"))
mata: RNPR = strtoreal(st_sdata(., "n_free"))

import delimited using "fixtures/expected/lavaan_modelimplied_nullvec.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: NFIT = st_sdata(., "fit")
mata: NKEY = st_sdata(., "param")
mata: NVAL = strtoreal(st_sdata(., "value"))

mata:

real scalar td_get(string colvector fitcol, string colvector keycol,
                   real colvector valcol,
                   string scalar wfit, string scalar wkey)
{
    real scalar i

    for (i = 1; i <= rows(fitcol); i++) {
        if (fitcol[i] == wfit) {
            if (keycol[i] == wkey) return(valcol[i])
        }
    }
    return(.)
}

string scalar td_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The oracle's key for one row of our parameter table. Thresholds keep their
// Section 4.2 label; structural parameters are named by matrix element.
string scalar td_key(struct dwlsmv_parmap scalar par, real scalar i)
{
    string scalar rc

    if (par.ptype[i] == 4) {
        return(par.name[i])
    }
    rc = strofreal(par.row[i]) + "_" + strofreal(par.col[i])
    if (par.ptype[i] == 1) return("lambda_" + rc)
    if (par.ptype[i] == 2) return("psi_"    + rc)
    return("beta_" + rc)
}

// One column key per free parameter, in theta order.
string colvector td_colkeys(struct dwlsmv_parmap scalar par)
{
    string colvector ck
    real scalar i

    ck = J(par.npar, 1, "")
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) ck[par.fidx[i]] = td_key(par, i)
    }
    return(ck)
}

string scalar td_canon(string scalar fit)
{
    if (fit == "D_marker") {
        return(td_spec(("factor f", "ind y1 y2 y3 y4",
                        "lambda f y2 free", "lambda f y3 free",
                        "lambda f y4 free", "psi f f free",
                        "metric marker f y1")))
    }
    if (fit == "D_unitvar") {
        return(td_spec(("factor f", "ind y1 y2 y3 y4",
                        "lambda f y1 free", "lambda f y2 free",
                        "lambda f y3 free", "lambda f y4 free",
                        "metric unitvariance f")))
    }
    if (fit == "E_marker") {
        return(td_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "beta f2 f1 free",
                        "psi f1 f1 free", "psi f2 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    return(td_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                    "lambda f1 y1 free", "lambda f1 y2 free",
                    "lambda f1 y3 free", "lambda f2 y4 free",
                    "lambda f2 y5 free", "lambda f2 y6 free",
                    "beta f2 f1 free",
                    "metric unitvariance f1", "metric unitvariance f2")))
}

// theta at the oracle's stated values.
real colvector td_theta(struct dwlsmv_parmap scalar par, string scalar fit,
                        string colvector IFIT, string colvector IKEY,
                        real colvector IVAL)
{
    real colvector th
    real scalar i

    th = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            th[par.fidx[i]] = td_get(IFIT, IKEY, IVAL, fit, td_key(par, i))
        }
    }
    return(th)
}

// The oracle's Jacobian, scattered into our (row, column) layout by label.
// One pass over the fixture, so nothing here is quadratic in its length.
real matrix td_fixjac(string colvector JFIT, string colvector JROW,
                      string colvector JCOL, real colvector JVAL,
                      string scalar fit,
                      struct dwlsmv_statmap scalar map,
                      struct dwlsmv_parmap scalar par)
{
    real matrix G
    string colvector ck
    real scalar i, r, c

    G  = J(map.D, par.npar, .)
    ck = td_colkeys(par)

    for (i = 1; i <= rows(JFIT); i++) {
        if (JFIT[i] != fit) continue
        r = dwlsmv_canon_find(map.label, JROW[i])
        c = dwlsmv_canon_find(ck, JCOL[i])
        if (r > 0 & c > 0) G[r, c] = JVAL[i]
    }
    return(G)
}

end

mata:

real scalar td_checks(string colvector IFIT, string colvector IKEY, real colvector IVAL,
                      string colvector JFIT, string colvector JROW,
                      string colvector JCOL, real colvector JVAL,
                      string colvector RFIT, real colvector RRNK, real colvector RNPR,
                      string colvector NFIT, string colvector NKEY, real colvector NVAL)
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar map, mapD, mapE
    string colvector FAIL, ck
    string rowvector fits
    string scalar fit
    real matrix Dm, G, U, V, Lam, Psi, B
    real colvector theta, sv, nullv, vals, th6, ana, p12s, dbp
    real matrix Jm
    real scalar f, i, j, k, m, worst, gap, oracle, lam_j, lam_k, bet
    real scalar rank, tol, nfree, kb, kp, want, q, p12, f2max

    FAIL = J(0, 1, "")
    mapD = dwlsmv_statmap_build((4 \ 4 \ 3 \ 3))
    mapE = dwlsmv_statmap_build((3 \ 3 \ 3 \ 3 \ 3 \ 3))
    fits = ("D_marker", "D_unitvar", "E_marker", "E_unitvar")

    for (f = 1; f <= cols(fits); f++) {
        fit = fits[f]
        if (substr(fit, 1, 1) == "D") {
            map = mapD
        }
        else {
            map = mapE
        }

        par   = dwls_identify(td_canon(fit), map)
        theta = td_theta(par, fit, IFIT, IKEY, IVAL)
        if (missing(theta)) {
            FAIL = FAIL \ (fit + ": theta was not fully filled from the oracle")
            continue
        }

        Dm = dwls_delta(theta, par, map)
        if (rows(Dm) != map.D | cols(Dm) != par.npar) {
            FAIL = FAIL \ (fit + ": Delta is not D x npar")
            continue
        }

        // ---- 1. the structural blocks, exactly --------------------------
        worst = 0
        gap   = 0
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] == 4) {
                // this threshold's column must be a unit vector
                for (j = 1; j <= map.D; j++) {
                    want  = (j == dwlsmv_pos_tau(map, par.row[i], par.col[i]))
                    worst = max((worst, abs(Dm[j, par.fidx[i]] - want)))
                }
            }
            else if (par.free[i]) {
                // a structural column must be exactly zero in every
                // threshold row
                for (j = 1; j <= map.T; j++) {
                    gap = max((gap, abs(Dm[j, par.fidx[i]])))
                }
            }
        }
        printf("{txt}  1. %-10s blocks          : max |tau block - I| = %g, max |d(tau)/d(structural)| = %g\n",
               fit, worst, gap)
        if (worst != 0) {
            FAIL = FAIL \ (fit + ": the threshold block of Delta is not a literal identity")
        }
        if (gap != 0) {
            FAIL = FAIL \ (fit + ": d(tau)/d(structural) is not exactly zero")
        }

        // ---- 4. against lavaan's analytic Jacobian -----------------------
        G     = td_fixjac(JFIT, JROW, JCOL, JVAL, fit, map, par)
        worst = 0
        if (missing(G)) {
            FAIL = FAIL \ (fit + ": the oracle Jacobian has entries we could not key")
        }
        else {
            worst = dwlsmv_gapratio(Dm, G, 1e-8, 1e-8)
        }
        printf("{txt}  4. %-10s vs lavaan       : worst ratio %8.2e of the 1e-8 budget over %f x %f entries\n",
               fit, worst, map.D, par.npar)
        if (worst > 1) {
            FAIL = FAIL \ (fit + ": Delta disagrees with lavaan, ratio " + strofreal(worst))
        }

        // ---- 5. rank ------------------------------------------------------
        U  = J(0, 0, 0)
        V  = J(0, 0, 0)
        sv = J(0, 1, 0)
        svd(Dm, U, sv, V)
        tol   = sv[1] * 1e-8
        rank  = sum(sv :> tol)
        nfree = td_get(RFIT, RFIT, RNPR, fit, fit)
        printf("{txt}  5. %-10s rank            : %f of %f free, singular values %g down to %g\n",
               fit, rank, par.npar, sv[1], sv[par.npar])
        if (rank != par.npar) {
            FAIL = FAIL \ (fit + ": Delta is rank deficient on an identified model")
        }
        if (nfree != par.npar) {
            FAIL = FAIL \ (fit + ": our npar " + strofreal(par.npar) +
                           " differs from the oracle's " + strofreal(nfree))
        }

        // ---- 2. the single-factor closed form ----------------------------
        // R_jk = lambda_j lambda_k under a unit-variance single factor, so
        // d R_jk / d lambda_m = lambda_k if m == j, lambda_j if m == k, else 0.
        if (fit == "D_unitvar") {
            dwlsmv_matrices(dwlsmv_parvals(theta, par), par, map.p,
                            dwlsmv_nfac(par), Lam, Psi, B)
            worst = 0
            for (i = 1; i <= par.nall; i++) {
                // the free guard is not decoration: a FIXED row has fidx == 0,
                // and par.fidx[i] would then abort with a bare subscript error
                // instead of a diagnosis. D_unitvar has no fixed loading today,
                // which is exactly why this would go unnoticed.
                if (par.ptype[i] != 1 | !par.free[i]) continue
                m = par.row[i]
                for (j = 1; j <= map.p - 1; j++) {
                    for (k = j + 1; k <= map.p; k++) {
                        lam_j  = Lam[j, 1]
                        lam_k  = Lam[k, 1]
                        oracle = 0
                        if (m == j) oracle = lam_k
                        if (m == k) oracle = lam_j
                        worst = max((worst,
                                     dwlsmv_gapratio(Dm[dwlsmv_pos_rho(map, j, k),
                                                        par.fidx[i]],
                                                     oracle, 1e-8, 1e-8)))
                    }
                }
            }
            printf("{txt}  2. closed form, lambda      : worst ratio %8.2e of the 1e-8 budget (R_jk = lambda_j lambda_k)\n",
                   worst)
            if (worst > 1) {
                FAIL = FAIL \ "the lambda block of Delta disagrees with its closed form"
            }
        }

        // ---- 3. the beta column, in closed form --------------------------
        // Phi = [[1, beta], [beta, beta^2 + 1]] under unit variance on both
        // factors, so d R_jk / d beta is 0 within F1, lambda_j lambda_k
        // across, and 2 beta lambda_j lambda_k within F2.
        if (fit == "E_unitvar") {
            vals = dwlsmv_parvals(theta, par)
            dwlsmv_matrices(vals, par, map.p, dwlsmv_nfac(par), Lam, Psi, B)
            bet   = B[2, 1]
            kb    = 0
            for (i = 1; i <= par.nall; i++) {
                if (par.ptype[i] == 3 & par.free[i]) kb = par.fidx[i]
            }
            if (kb == 0) {
                FAIL = FAIL \ (fit + ": no free beta row to take a Delta column for")
                continue
            }
            worst = 0
            for (j = 1; j <= map.p - 1; j++) {
                for (k = j + 1; k <= map.p; k++) {
                    if (j <= 3 & k <= 3) {
                        oracle = 0
                    }
                    else if (j <= 3) {
                        oracle = Lam[j, 1] * Lam[k, 2]
                    }
                    else {
                        oracle = 2 * bet * Lam[j, 2] * Lam[k, 2]
                    }
                    worst = max((worst,
                                 dwlsmv_gapratio(Dm[dwlsmv_pos_rho(map, j, k), kb],
                                                 oracle, 1e-8, 1e-8)))
                }
            }
            printf("{txt}  3. closed form, beta        : worst ratio %8.2e of the 1e-8 budget at beta = %g\n",
                   worst, bet)
            if (worst > 1) {
                FAIL = FAIL \ "the beta column of Delta disagrees with its closed form"
            }
        }
    }

    // =====================================================================
    // 6. The overparameterized case, at the oracle's beta rather than at
    //    Module 3's default start. dwls_identify() will not return this
    //    table, so it is built by appending the offending row to a valid one.
    // =====================================================================
    par   = dwls_identify(td_canon("E_unitvar"), mapE)
    theta = td_theta(par, "E_unitvar", IFIT, IKEY, IVAL)

    par.nall  = par.nall + 1
    par.name  = par.name  \ "psi[f2,f1]"
    par.ptype = par.ptype \ 2
    par.row   = par.row   \ 2
    par.col   = par.col   \ 1
    par.free  = par.free  \ 1
    par.value = par.value \ 0
    par.npar  = par.npar + 1
    par.fidx  = par.fidx  \ par.npar
    theta     = theta \ 0
    th6       = theta

    Dm = dwls_delta(theta, par, mapE)
    U  = J(0, 0, 0)
    V  = J(0, 0, 0)
    sv = J(0, 1, 0)
    svd(Dm, U, sv, V)
    tol  = sv[1] * 1e-8
    rank = sum(sv :> tol)

    printf("{txt}  6. overparameterized E      : rank %f of %f free; singular values %g down to %g, ratio %8.2e\n",
           rank, par.npar, sv[1], sv[par.npar], sv[par.npar] / sv[1])
    if (rank != par.npar - 1) {
        FAIL = FAIL \ ("the overparameterized case E has deficiency " +
                       strofreal(par.npar - rank) + ", expected 1")
    }
    if (td_get(RFIT, RFIT, RRNK, "E_over", "E_over") != rank) {
        FAIL = FAIL \ "our rank for E_over differs from the oracle's"
    }

    // The ratio is ASSERTED and not merely printed. Without this the rank
    // check's margin could decay by eight orders of magnitude and the
    // deficiency test above would still pass, because it only asks which side
    // of 1e-8 the value falls on. 1e-12 is far looser than the ~2e-17
    // observed and far tighter than the threshold it protects.
    if (sv[par.npar] / sv[1] > 1e-12) {
        FAIL = FAIL \ ("the E_over singular-value ratio is " +
                       strofreal(sv[par.npar] / sv[1]) + ", no longer machine zero")
    }

    // The whole Delta against the oracle, not just its rank. The fixture
    // carries all 27 x 20 entries for E_over, and the machinery to key them
    // is already here, so there is no reason to compare only a summary.
    G = td_fixjac(JFIT, JROW, JCOL, JVAL, "E_over", mapE, par)
    if (missing(G)) {
        FAIL = FAIL \ "the oracle Jacobian for E_over has entries we could not key"
    }
    else {
        worst = dwlsmv_gapratio(Dm, G, 1e-8, 1e-8)
        printf("{txt}     whole Delta vs lavaan    : worst ratio %8.2e of the 1e-8 budget over %f x %f entries\n",
               worst, mapE.D, par.npar)
        if (worst > 1) {
            FAIL = FAIL \ ("E_over: Delta disagrees with lavaan, ratio " + strofreal(worst))
        }
    }

    // The null direction must name beta and psi_21 and nothing else. At
    // beta != 0 this holds only because psi_21 itself is 0 -- see the header.
    ck    = td_colkeys(par)
    nullv = V[par.npar, .]'
    kb    = dwlsmv_canon_find(ck, "beta_2_1")
    kp    = dwlsmv_canon_find(ck, "psi_2_1")
    worst = 0
    for (i = 1; i <= par.npar; i++) {
        if (i != kb & i != kp) worst = max((worst, abs(nullv[i])))
    }
    printf("{txt}     null direction           : beta %9.6f, psi_21 %9.6f, largest other |entry| %g\n",
           nullv[kb], nullv[kp], worst)

    if (abs(abs(nullv[kb]) - sqrt(0.5)) > 1e-3) {
        FAIL = FAIL \ "the null direction does not load +/- 1/sqrt(2) on beta"
    }
    if (abs(abs(nullv[kp]) - sqrt(0.5)) > 1e-3) {
        FAIL = FAIL \ "the null direction does not load +/- 1/sqrt(2) on psi_21"
    }
    if (nullv[kb] * nullv[kp] > 0) {
        FAIL = FAIL \ "the null direction has the same sign on beta and psi_21"
    }
    if (worst > 1e-3) {
        FAIL = FAIL \ ("the null direction loads " + strofreal(worst) +
                       " on a parameter that is not beta or psi_21")
    }

    // and it must agree with the oracle's, up to sign
    gap = abs(abs(nullv[kb]) - abs(td_get(NFIT, NKEY, NVAL, "E_over", "beta_2_1")))
    printf("{txt}     vs the oracle's           : |beta component| differs by %g\n", gap)
    if (gap > 1e-3) {
        FAIL = FAIL \ "our null direction disagrees with the oracle's on beta"
    }

    // ---- 7. THE SAME TABLE AT MODULE 3'S OWN DEFAULT STARTS ---------------
    // Everything above sits at the oracle's beta = 0.4. dwls_identify() runs
    // its rank check at the DEFAULT starts instead -- lambda 0.7, beta 0,
    // psi off-diagonal 0 -- and that is the configuration the development record
    // quotes a ratio for. Measuring it here turns that quotation into an
    // assertion: reject 1 of test_identify.do only establishes rc 3498, and
    // the number it prints reaches the log through an errprintf inside
    // `capture noisily`, where nothing is checking it.
    //
    // This is also the comparison that isolates WHICH parameter is the knife
    // edge. Both configurations hold psi_21 = 0 and differ only in beta; if
    // beta were the sensitive one, the two ratios would differ by orders.
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) theta[par.fidx[i]] = par.value[i]
    }

    Dm = dwls_delta(theta, par, mapE)
    U  = J(0, 0, 0)
    V  = J(0, 0, 0)
    sv = J(0, 1, 0)
    svd(Dm, U, sv, V)
    rank = sum(sv :> sv[1] * 1e-8)

    printf("{txt}  7. same table at the default starts: rank %f of %f, ratio %8.2e (beta = 0 here, beta = 0.4 above)\n",
           rank, par.npar, sv[par.npar] / sv[1])
    if (rank != par.npar - 1) {
        FAIL = FAIL \ "the overparameterized table is not rank deficient at the default starts"
    }
    if (sv[par.npar] / sv[1] > 1e-12) {
        FAIL = FAIL \ ("at the default starts the ratio is " +
                       strofreal(sv[par.npar] / sv[1]) + ", no longer machine zero")
    }

    // ---- 8. THE SAME TABLE AT A NONZERO psi_21 -----------------------------
    // Added 2026-09-23 with the Gate 2 rewrite (derivation, Finding 4). Until
    // then the spec and the project's notes said a Delta evaluated at a nonzero psi_21
    // comes back FULL RANK on a model that is still unidentified, and nothing
    // measured it. It does not. With psi_11 = psi_22 = 1 fixed,
    // (beta, psi_21) -> (Phi_21, Phi_22) has Jacobian
    //     J = [1, 1; 2(beta + psi_21), 2 beta],   det J = -2 psi_21,
    // so at psi_21 != 0 the pair can move Phi_22, which reopens F2's scale
    // against its loadings: deficiency still exactly one, but the null
    // direction now rescales F2's loadings, d lambda_k = e lambda_k, with
    // (d beta, d psi_21) = -e J^-1 (Phi_21, 2 Phi_22)'. psi_21 = 0 is only where
    // that direction collapses onto (beta, psi_21) alone.
    //
    // Asserted: deficiency 1 at psi_21 = 0.1 and 0.3; the null direction
    // against that analytic vector, up to sign, at check 6's 1e-3 singular-
    // vector tolerance; and a ratio below 1e-10. That bound is not machine
    // zero, and deliberately: here beta and psi_21 are no longer the same
    // arithmetic, so the floor is the central difference's rounding at
    // h = 1e-5 (sigma is at most quadratic in any one parameter, so there is
    // no truncation term). Measured 2026-09-23: 4.8e-13 at 0.1 and 8.9e-13 at
    // 0.3, next singular value 0.32 and 0.27 of the largest, null vector
    // within 4.3e-12 of the analytic one. 1e-10 is two orders above that
    // floor and two under the 1e-8 threshold the rank check uses -- so the
    // check keeps a four-order margin at a nonzero psi_21, not the sixteen it
    // has at psi_21 = 0.
    p12s = (0.1 \ 0.3)
    for (q = 1; q <= rows(p12s); q++) {
        p12       = p12s[q]
        theta     = th6
        theta[kp] = p12
        bet       = theta[kb]

        Dm = dwls_delta(theta, par, mapE)
        U  = J(0, 0, 0)
        V  = J(0, 0, 0)
        sv = J(0, 1, 0)
        svd(Dm, U, sv, V)
        rank  = sum(sv :> sv[1] * 1e-8)
        nullv = V[par.npar, .]'

        // the analytic null direction, in theta order
        ana = J(par.npar, 1, 0)
        for (i = 1; i <= par.nall; i++) {
            if (par.free[i] & par.ptype[i] == 1 & par.col[i] == 2) {
                ana[par.fidx[i]] = theta[par.fidx[i]]
            }
        }
        Jm  = (1, 1 \ 2 * (bet + p12), 2 * bet)
        dbp = -luinv(Jm) * ((bet + p12) \ 2 * (1 + bet * bet + 2 * bet * p12))
        ana[kb] = dbp[1]
        ana[kp] = dbp[2]
        ana     = ana / sqrt(quadcross(ana, ana))
        if (nullv[kb] * ana[kb] < 0) {
            ana = -ana
        }
        gap   = max(abs(nullv - ana))
        f2max = 0
        for (i = 1; i <= par.nall; i++) {
            if (par.free[i] & par.ptype[i] == 1 & par.col[i] == 2) {
                f2max = max((f2max, abs(nullv[par.fidx[i]])))
            }
        }

        printf("{txt}  8. same table at psi_21 = %3.1f : rank %f of %f, ratio %8.2e, next %8.2e ; null vs analytic %8.2e ; beta %9.6f psi_21 %9.6f, largest F2 loading entry %8.6f\n",
               p12, rank, par.npar, sv[par.npar] / sv[1], sv[par.npar - 1] / sv[1],
               gap, nullv[kb], nullv[kp], f2max)
        if (rank != par.npar - 1) {
            FAIL = FAIL \ ("at psi_21 = " + strofreal(p12) + " the overparameterized table has rank " +
                           strofreal(rank) + " of " + strofreal(par.npar) + ", not deficiency 1")
        }
        if (sv[par.npar] / sv[1] > 1e-10) {
            FAIL = FAIL \ ("at psi_21 = " + strofreal(p12) + " the ratio is " +
                           strofreal(sv[par.npar] / sv[1]) + ", above the 1e-10 floor bound")
        }
        if (gap > 1e-3) {
            FAIL = FAIL \ ("at psi_21 = " + strofreal(p12) +
                           " the null direction is not the loading rescaling, gap " + strofreal(gap))
        }
    }

    printf("{txt}test_delta: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(td_checks(IFIT, IKEY, IVAL,       ///
                                            JFIT, JROW, JCOL, JVAL, ///
                                            RFIT, RRNK, RNPR,       ///
                                            NFIT, NKEY, NVAL)))

if `nfail' {
    display as error _n "test_delta: `nfail' failure(s)"
    exit 9
}

display as result "test_delta: the threshold block of Delta is a literal identity and both off-blocks are exactly zero, the lambda and beta columns match their closed forms, Delta matches lavaan's analytic Jacobian on four configurations, and the 4.8 overparameterization shows deficiency 1 with a null direction of +/- 1/sqrt(2) on (beta, psi_21) at psi_21 = 0, and still deficiency 1 at a nonzero psi_21, where the null direction rescales F2's loadings"
exit 0
