*! Module 5: dwls_model() -- the model-implied moment vector sigma(theta).
*!
*! Executed under StataMP 18 on 2026-09-22 and it passes.
*!
*! WHAT THIS FILE IS FOR. Spec Section 4.6 maps theta to Sigma_y*, to R_y*, to
*! the model-implied polychorics, and states that under the delta
*! parameterization Theta is DERIVED and the standardization is a no-op. Each
*! of those claims is checked here, and the whole mapping is checked against
*! lavaan.
*!
*!   1   sigma(theta) against lavaan's own model-implied vector, on four
*!       configurations -- case D and case E, each under both factor metrics.
*!       The oracle is fixtures/expected/lavaan_modelimplied_*.csv, generated
*!       at a STATED (Lambda, Psi, B) with `do.fit = FALSE', so nothing in the
*!       comparison depends on an optimizer having converged and nothing
*!       depends on Module 3's parameter ordering: the fixture is keyed on
*!       matrix elements and Section 4.2 labels, and this file does the
*!       translation.
*!   2   the correspondence itself -- every parameter our table carries is
*!       present in the fixture with the same value, and every fixture element
*!       our table does NOT carry is zero. Without this, check 1 could pass
*!       while comparing a different model that happens to imply the same
*!       correlations.
*!   3   the threshold entries of sigma are literal copies of theta's, exactly.
*!       Note this one cannot fail against the current implementation:
*!       dwls_model() writes sigma[pos_tau(...)] = vals[i] using the same call
*!       this check positions with, so it is a statement about the code's
*!       shape rather than a test of it. It would catch a dwls_model() that
*!       recomputed thresholds instead of copying them, which is the mistake
*!       worth guarding, and it costs nothing. The check that the ORDER is
*!       right is test_identify.do check 2, and the check that the values are
*!       right is check 1 here, against a labelled oracle.
*!   4   R_y* has a unit diagonal, and D^-1/2 Sigma D^-1/2 equals Sigma, both
*!       to machine zero. Gate 2's Finding 1.
*!   5   the other half of Finding 1: scaling lambda_j by c and theta_jj by
*!       c^2 leaves R_y* untouched, which is WHY a free Theta is unidentified
*!       from correlations alone and not merely unconventional.
*!   6   derived Theta against lavaan's fitted residual variances on case E.
*!       An independent fixture, and the one that would catch us deriving
*!       Theta correctly from the wrong Phi.
*!   7   Var(F2) = beta^2 + psi_22 for the ENDOGENOUS factor, against lavaan's
*!       implied variance. The 2026-09-21 metric decision, arriving here as an
*!       arithmetic consequence of Phi = (I-B)^-1 Psi (I-B)^-T.
*!   8   the equivalent-models identity at the level of sigma: a correlated
*!       two-factor model and an F2 <- F1 model, at lavaan's own estimates for
*!       each, imply the same R_y*.
*!
*! TOLERANCES. Check 1 compares two implementations of the same closed-form
*! mapping evaluated at identical inputs, so it is a machine-precision
*! comparison and gets spec Section 7.2's exact-identity row: 1e-12 absolute
*! with no relative term. Checks 6-8 compare across two separately converged
*! lavaan fits, so they get the parameter-estimate tier, 1e-5 absolute and
*! 1e-4 relative. Neither number is tuned to what was observed; both are
*! printed so a near-miss is visible.
*!
*! Reference CSVs are read with stringcols(_all) and converted with
*! strtoreal(), and with delimiters(",") specified. `import delimited' types a
*! numeric column as FLOAT without warning, capping agreement near 4e-08, and
*! it auto-detects `|' as the delimiter when labels like y1|t1 make pipes as
*! frequent as commas -- which they do in every file read here.
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
        display as error "test_model: `f' failed to load"
        exit 9
    }
}

* ---- tier 2: the stated (Lambda, Psi, B, tau) behind each configuration ----
import delimited using "fixtures/expected/lavaan_modelimplied_input.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: IFIT = st_sdata(., "fit")
mata: IMAT = st_sdata(., "matrix")
mata: IKEY = st_sdata(., "element")
mata: IVAL = strtoreal(st_sdata(., "value"))

* ---- tier 2: lavaan's model-implied sigma for each configuration -----------
import delimited using "fixtures/expected/lavaan_modelimplied_sigma.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: SFIT = st_sdata(., "fit")
mata: SKEY = st_sdata(., "label")
mata: SVAL = strtoreal(st_sdata(., "value"))

* ---- tier 2: lavaan's FITTED case E parameters, for checks 6-8 -------------
import delimited using "fixtures/expected/lavaan_case_E_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: PFIT = st_sdata(., "fit")
mata: PKEY = st_sdata(., "label")
mata: PVAL = strtoreal(st_sdata(., "est"))

mata:

// ---------------------------------------------------------------------------
// Long-format lookup. Returns missing when the (fit, key) pair is absent,
// which every caller reports rather than treating as zero.
// ---------------------------------------------------------------------------
real scalar tm_get(string colvector fitcol, string colvector keycol,
                   real colvector valcol,
                   string scalar wfit, string scalar wkey)
{
    real scalar i

    for (i = 1; i <= rows(fitcol); i++) {
        if (fitcol[i] == wfit & keycol[i] == wkey) {
            return(valcol[i])
        }
    }
    return(.)
}

string scalar tm_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The fixture's key for one row of our parameter table. Thresholds keep their
// Section 4.2 label; structural parameters are named by matrix element, which
// is the interface the oracle was deliberately built on so that it does not
// depend on theta's ordering.
string scalar tm_key(struct dwlsmv_parmap scalar par, real scalar i)
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

// The canonical specification for each configuration of the oracle. Written
// out rather than generated from the fixture, so that check 2 has something
// independent to verify the correspondence against.
string scalar tm_canon(string scalar fit)
{
    if (fit == "D_marker") {
        return(tm_spec(("factor f", "ind y1 y2 y3 y4",
                        "lambda f y2 free", "lambda f y3 free",
                        "lambda f y4 free", "psi f f free",
                        "metric marker f y1")))
    }
    if (fit == "D_unitvar") {
        return(tm_spec(("factor f", "ind y1 y2 y3 y4",
                        "lambda f y1 free", "lambda f y2 free",
                        "lambda f y3 free", "lambda f y4 free",
                        "metric unitvariance f")))
    }
    if (fit == "E_marker") {
        return(tm_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "beta f2 f1 free",
                        "psi f1 f1 free", "psi f2 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    if (fit == "E_unitvar") {
        return(tm_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                        "lambda f1 y1 free", "lambda f1 y2 free",
                        "lambda f1 y3 free", "lambda f2 y4 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "beta f2 f1 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    // The two case E parameterizations of checks 6-8, at lavaan's FITTED
    // estimates rather than at the oracle's stated values.
    if (fit == "cor_unitvar") {
        return(tm_spec(("factor f1 f2", "ind y1 y2 y3 y4 y5 y6",
                        "lambda f1 y1 free", "lambda f1 y2 free",
                        "lambda f1 y3 free", "lambda f2 y4 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "psi f1 f2 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    return(tm_canon("E_unitvar"))
}

end

mata:

real scalar tm_checks(string colvector IFIT, string colvector IMAT,
                      string colvector IKEY, real colvector IVAL,
                      string colvector SFIT, string colvector SKEY,
                      real colvector SVAL,
                      string colvector PFIT, string colvector PKEY,
                      real colvector PVAL)
{
    struct dwlsmv_parmap  scalar par, parc
    struct dwlsmv_statmap scalar map, mapD, mapE
    string colvector FAIL
    string rowvector fits
    string scalar fit, key, lab
    real colvector theta, sigma, th2, Th, Thc
    real matrix Lam, Psi, B, Lam2, R, R2, Rc, Sig, A, Phi, Phic
    real scalar f, i, j, k, v, worst, ratio, nT, cc, gap, b2, v2, lv2

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

        par   = dwls_identify(tm_canon(fit), map)
        nT    = map.T
        theta = J(par.npar, 1, .)

        // ---- 2. correspondence, in both directions ----------------------
        // Forward: every row we carry is in the fixture with the same value.
        for (i = 1; i <= par.nall; i++) {
            key = tm_key(par, i)
            v   = tm_get(IFIT, IKEY, IVAL, fit, key)
            if (missing(v)) {
                FAIL = FAIL \ (fit + ": the oracle has no element " + key)
                continue
            }
            if (par.free[i]) {
                theta[par.fidx[i]] = v
            }
            else if (dwlsmv_gapratio(par.value[i], v, 1e-12, 0) > 1) {
                FAIL = FAIL \ (fit + ": fixed " + key + " is " +
                               strofreal(par.value[i]) + " here, " + strofreal(v) +
                               " in the oracle")
            }
        }
        if (missing(theta)) {
            FAIL = FAIL \ (fit + ": theta was not fully filled from the oracle")
            continue
        }

        // Reverse: every structural element the oracle states and we do NOT
        // carry a row for must be zero. This is the half that catches a
        // specification quietly missing a loading.
        for (i = 1; i <= rows(IFIT); i++) {
            if (IFIT[i] != fit | IMAT[i] == "tau") {
                continue
            }
            k = 0
            for (j = 1; j <= par.nall; j++) {
                if (tm_key(par, j) == IKEY[i]) k = j
            }
            if (k == 0 & IVAL[i] != 0) {
                FAIL = FAIL \ (fit + ": the oracle states " + IKEY[i] + " = " +
                               strofreal(IVAL[i]) + " but our table has no such row")
            }
        }

        // ---- 1. sigma(theta) against lavaan ------------------------------
        sigma = dwls_model(theta, par, map)
        worst = 0
        for (i = 1; i <= map.D; i++) {
            v = tm_get(SFIT, SKEY, SVAL, fit, map.label[i])
            if (missing(v)) {
                FAIL = FAIL \ (fit + ": the oracle has no sigma entry " + map.label[i])
                continue
            }
            ratio = dwlsmv_gapratio(sigma[i], v, 1e-12, 0)
            worst = max((worst, ratio))
        }
        printf("{txt}  1. %-10s sigma vs lavaan : worst ratio %8.2e of the 1e-12 budget over %f entries\n",
               fit, worst, map.D)
        if (worst > 1) {
            FAIL = FAIL \ (fit + ": sigma disagrees with lavaan, ratio " + strofreal(worst))
        }

        // ---- 3. the threshold entries are literal copies ------------------
        // Exact, not to a tolerance. Section 4.8 orders theta with thresholds
        // first precisely so that no arithmetic happens here.
        worst = 0
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] != 4) continue
            worst = max((worst,
                         abs(sigma[dwlsmv_pos_tau(map, par.row[i], par.col[i])] -
                             theta[par.fidx[i]])))
        }
        if (f == 1) {
            printf("{txt}  3. thresholds are copies    : max |sigma_tau - theta_tau| = %g\n", worst)
        }
        if (worst != 0) {
            FAIL = FAIL \ (fit + ": a threshold entry of sigma is not theta's own value")
        }

        // ---- 4. Finding 1: the standardization is a no-op ----------------
        dwlsmv_matrices(dwlsmv_parvals(theta, par), par, map.p, dwlsmv_nfac(par),
                        Lam, Psi, B)
        A   = luinv(I(rows(Psi)) - B)
        Phi = A * Psi * A'
        R   = dwlsmv_rystar(Lam, Psi, B)
        Sig = Lam * Phi * Lam'
        Th  = 1 :- diagonal(Sig)
        Sig = Sig + diag(Th)

        worst = max(abs(diagonal(R) :- 1))
        gap   = 0
        for (i = 1; i <= map.p; i++) {
            for (j = 1; j <= map.p; j++) {
                gap = max((gap, abs(R[i, j] -
                                    Sig[i, j] / sqrt(Sig[i, i] * Sig[j, j]))))
            }
        }
        printf("{txt}  4. %-10s Finding 1       : max |diag(R) - 1| = %g, max |R - D^-1/2 Sigma D^-1/2| = %g\n",
               fit, worst, gap)
        if (worst != 0) {
            FAIL = FAIL \ (fit + ": R_y* does not have an exactly unit diagonal")
        }
        if (gap > 1e-15) {
            FAIL = FAIL \ (fit + ": standardizing Sigma is not a no-op, gap " + strofreal(gap))
        }

        // ---- 5. and why a free Theta would be unidentified ---------------
        // Scale lambda_1 by c and theta_11 by c^2. R is untouched, so
        // correlations alone cannot separate the two.
        cc         = 2.5
        Lam2       = Lam
        Lam2[1, .] = Lam[1, .] :* cc
        Th         = 1 :- diagonal(Lam * Phi * Lam')
        Th[1]      = Th[1] * cc * cc
        Sig        = Lam2 * Phi * Lam2' + diag(Th)

        // Standardize from an untouched copy: dividing R2 in place would
        // corrupt the very diagonal the later divisions read.
        R2 = Sig
        for (i = 1; i <= map.p; i++) {
            for (j = 1; j <= map.p; j++) {
                R2[i, j] = Sig[i, j] / sqrt(Sig[i, i] * Sig[j, j])
            }
        }
        worst = max(abs(R2 - R))
        if (f == 1) {
            printf("{txt}  5. free Theta unidentified  : scaling lambda_1 by %g and theta_11 by %g moves R by %g\n",
                   cc, cc * cc, worst)
        }
        if (worst > 1e-14) {
            FAIL = FAIL \ (fit + ": the lambda/theta rescaling moved R, so Theta is separable after all")
        }
    }

    // =====================================================================
    // Checks 6-8, at lavaan's FITTED case E estimates.
    // =====================================================================
    par   = dwls_identify(tm_canon("str_unitvar"), mapE)
    theta = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 4) {
            lab = par.name[i]
        }
        else if (par.ptype[i] == 1) {
            lab = "f" + strofreal(par.col[i]) + "=~y" + strofreal(par.row[i])
        }
        else if (par.ptype[i] == 3) {
            lab = "f" + strofreal(par.row[i]) + "~f" + strofreal(par.col[i])
        }
        else {
            lab = "f" + strofreal(par.row[i]) + "~~f" + strofreal(par.col[i])
        }
        v = tm_get(PFIT, PKEY, PVAL, "str_unitvar", lab)
        if (par.free[i]) {
            if (missing(v)) {
                FAIL = FAIL \ ("str_unitvar: no lavaan estimate labelled " + lab)
            }
            else {
                theta[par.fidx[i]] = v
            }
        }
    }

    if (!missing(theta)) {
        dwlsmv_matrices(dwlsmv_parvals(theta, par), par, mapE.p, dwlsmv_nfac(par),
                        Lam, Psi, B)
        A   = luinv(I(2) - B)
        Phi = A * Psi * A'
        R   = dwlsmv_rystar(Lam, Psi, B)
        Th  = 1 :- diagonal(Lam * Phi * Lam')

        // ---- 6. derived Theta against lavaan's residual variances --------
        worst = 0
        for (j = 1; j <= mapE.p; j++) {
            v = tm_get(PFIT, PKEY, PVAL, "str_unitvar",
                       "y" + strofreal(j) + "~~y" + strofreal(j))
            worst = max((worst, dwlsmv_gapratio(Th[j], v, 1e-5, 1e-4)))
        }
        printf("{txt}  6. derived Theta vs lavaan  : worst ratio %8.2e of the theta-hat budget over %f indicators\n",
               worst, mapE.p)
        if (worst > 1) {
            FAIL = FAIL \ ("derived Theta disagrees with lavaan's residual variances")
        }

        // ---- 7. Var(F2) = beta^2 + psi_22, the endogenous metric ---------
        b2  = tm_get(PFIT, PKEY, PVAL, "str_unitvar", "f2~f1")
        v2  = Phi[2, 2]
        lv2 = tm_get(PFIT, PKEY, PVAL, "str_unitvar", "f2~~f2")
        printf("{txt}  7. endogenous Var(F2)       : Phi_22 = %18.15f, beta^2 + psi_22 = %18.15f, psi_22 = %g\n",
               v2, b2 * b2 + 1, lv2)
        if (dwlsmv_gapratio(v2, b2 * b2 + 1, 1e-12, 0) > 1) {
            FAIL = FAIL \ "Phi_22 is not beta^2 + psi_22, so (I-B)^-1 Psi (I-B)^-T is wrong"
        }
        if (abs(v2 - 1) <= 0.1) {
            FAIL = FAIL \ "implied Var(F2) is 1, which is the TOTAL-variance reading 4.8 rejects"
        }

        // ---- 8. the equivalent-models identity, at the sigma level -------
        parc = dwls_identify(tm_canon("cor_unitvar"), mapE)
        th2  = J(parc.npar, 1, .)
        for (i = 1; i <= parc.nall; i++) {
            if (parc.ptype[i] == 4) {
                lab = parc.name[i]
            }
            else if (parc.ptype[i] == 1) {
                lab = "f" + strofreal(parc.col[i]) + "=~y" + strofreal(parc.row[i])
            }
            else {
                lab = "f" + strofreal(parc.col[i]) + "~~f" + strofreal(parc.row[i])
            }
            v = tm_get(PFIT, PKEY, PVAL, "cor_unitvar", lab)
            if (parc.free[i]) {
                if (missing(v)) {
                    FAIL = FAIL \ ("cor_unitvar: no lavaan estimate labelled " + lab)
                }
                else {
                    th2[parc.fidx[i]] = v
                }
            }
        }

        if (!missing(th2)) {
            dwlsmv_matrices(dwlsmv_parvals(th2, parc), parc, mapE.p,
                            dwlsmv_nfac(parc), Lam2, Psi, B)
            Rc    = dwlsmv_rystar(Lam2, Psi, B)
            worst = 0
            for (i = 1; i <= mapE.p; i++) {
                for (j = 1; j <= mapE.p; j++) {
                    worst = max((worst, dwlsmv_gapratio(Rc[i, j], R[i, j], 1e-5, 1e-4)))
                }
            }
            printf("{txt}  8. equivalent models        : worst ratio %8.2e of the theta-hat budget over R's %f entries\n",
                   worst, mapE.p * mapE.p)
            if (worst > 1) {
                FAIL = FAIL \ "the correlated and structural parameterizations imply different R_y*"
            }
        }
        else {
            FAIL = FAIL \ "cor_unitvar theta was not fully filled from lavaan's estimates"
        }
    }
    else {
        FAIL = FAIL \ "str_unitvar theta was not fully filled from lavaan's estimates"
    }

    printf("{txt}test_model: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(tm_checks(IFIT, IMAT, IKEY, IVAL, ///
                                            SFIT, SKEY, SVAL,       ///
                                            PFIT, PKEY, PVAL)))

if `nfail' {
    display as error _n "test_model: `nfail' failure(s)"
    exit 9
}

display as result "test_model: sigma(theta) matches lavaan's model-implied vector on four configurations, thresholds pass through exactly, the delta standardization is a no-op, derived Theta matches lavaan's residual variances, Var(F2) = beta^2 + psi_22, and the two case E parameterizations imply the same R_y*"
exit 0
