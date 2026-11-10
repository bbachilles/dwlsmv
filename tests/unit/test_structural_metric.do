*! Case E: the endogenous-factor metric, and the equivalent-models identity.
*!
*! Executed under StataMP 18 on 2026-09-21 and it passes. Part D, checks 12-15,
*! executed 2026-09-23 with Module 7 and passes; check 15's T3 assertions
*! executed 2026-09-24 with Module 9 and pass, as do checks 7-8 at the
*! tightened test-statistic tier.
*!
*! WHAT THIS FILE IS FOR. Two questions that need a factor to be PREDICTED by
*! another factor, which no case in the suite before E could supply.
*!
*!   A. THE METRIC QUESTION (spec 4.8). Under unit-variance identification,
*!      what is fixed to 1 for an ENDOGENOUS factor -- its TOTAL variance, or
*!      its DISTURBANCE variance? For an exogenous factor the two coincide, so
*!      a single-factor or correlated-factor design cannot tell them apart;
*!      once F2 <- F1 is in the model they differ by beta^2 and every loading
*!      on F2 moves with the answer.
*!
*!      Settled empirically against lavaan 0.7.2 rather than reasoned about:
*!      lavaan fixes the DISTURBANCE variance. Checks 1-5 below read that back
*!      out of the committed fixture, so a lavaan release that changes the
*!      convention fails here instead of silently invalidating spec 4.8.
*!
*!   B. THE EQUIVALENT-MODELS IDENTITY. A correlated two-factor model and an
*!      F2 <- F1 model are reparameterizations of each other: both spend one
*!      parameter on the link between the factors, and each maps onto the
*!      other without changing sigma(theta). Fitted to the same data they must
*!      give the same npar, the same df, the same fitted statistics and the
*!      same T3. Checks 6-8.
*!
*!      Note the two are INDEPENDENT. B holds whichever way A comes out, which
*!      is what makes B usable as an internal consistency check rather than a
*!      second way of asking A.
*!
*!   C. STAGE 1 ON CASE E. Case E is new to the suite and is the first design
*!      with six indicators, and the first whose population correlation matrix
*!      is Lambda Phi Lambda' rather than a single lambda^2. Nothing above is
*!      worth trusting if our own Stage 1 cannot process the data those fits
*!      were made on, so checks 9-11 put case_E.csv through
*!      dwls_sample_stats() before anything is concluded from the fixtures.
*!
*!   D. B ON OUR OWN FITS (spec 7.4's structural-path condition), added
*!      2026-09-23 with Module 7. Both parameterizations fitted with
*!      dwls_fit() under both metrics, and compared WITH EACH OTHER -- no
*!      external reference at all. Checks 12-15: both converge; same npar and
*!      df as each other and as lavaan; the same sigma(theta-hat); the same
*!      T_raw; and (check 15, with Module 9 on 2026-09-24) the same T3, which
*!      also matches lavaan's chisq_scaled at spec 7.2's test-statistic tier.
*!      Until Module 9 existed check 15 printed "NOT ASSERTED"; spec 4.8's
*!      four quantities -- npar, df, sigma(theta-hat), T3 -- are now all
*!      asserted on our own fits. This test never legitimately skips.
*!
*! Reference CSVs are read with stringcols(_all) and converted with
*! strtoreal(). `import delimited` types a numeric column as FLOAT without
*! warning, which caps agreement near 4e-08 -- looser than check 2 asserts.
*! delimiters(",") is specified for the reason the project's coding rules record: the labels in
*! these fixtures contain no pipes, but the rule is not per-file.
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
             src/mata/dwlsmv_asym.mata              ///
             src/mata/dwlsmv_acov.mata              ///
             src/mata/dwlsmv_traces.mata            ///
             src/mata/dwlsmv_scaleshift.mata        ///
             src/mata/dwls_inference.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_structural_metric: `f' failed to load"
        exit 9
    }
}

* ---- tier 2: lavaan's parameter-table readings, long format ----------------
import delimited using "fixtures/expected/lavaan_case_E_metric.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: MFIT = st_sdata(., "fit")
mata: MQTY = st_sdata(., "quantity")
mata: MVAL = strtoreal(st_sdata(., "value"))

* ---- tier 2: the four fits' scalar fit results -----------------------------
import delimited using "fixtures/expected/lavaan_case_E_equivalence.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: EFIT = st_sdata(., "fit")
mata: ECHI = strtoreal(st_sdata(., "chisq"))
mata: EDF  = strtoreal(st_sdata(., "df"))
mata: ET3  = strtoreal(st_sdata(., "chisq_scaled"))
mata: EDF3 = strtoreal(st_sdata(., "df_scaled"))
mata: ENP  = strtoreal(st_sdata(., "npar"))

* ---- tier 1: case E's POPULATION threshold block, for the bias screen ------
import delimited using "fixtures/expected/case_E_gamma_tau.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local gcols `r(varlist)'
mata: G_POP   = strtoreal(st_sdata(., "`gcols'"))
mata: LAB_POP = st_sdata(., "label")

import delimited using "fixtures/expected/case_E_thresholds.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TAU_POP = strtoreal(st_sdata(., "tau"))
mata: Q_POP   = strtoreal(st_sdata(., "q"))

* ---- the data --------------------------------------------------------------
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

// Lookup into the long-format metric fixture. Returns missing if the
// (fit, quantity) pair is absent, which the caller reports rather than
// silently treating as zero.
real scalar dwlsmv_m_get(string colvector fitcol, string colvector qtycol,
                         real colvector valcol,
                         string scalar wfit, string scalar wqty)
{
    real scalar i

    for (i = 1; i <= rows(fitcol); i++) {
        if (fitcol[i] == wfit & qtycol[i] == wqty) {
            return(valcol[i])
        }
    }
    return(.)
}

real scalar dwlsmv_m_row(string colvector fitcol, string scalar wfit)
{
    real scalar i

    for (i = 1; i <= rows(fitcol); i++) {
        if (fitcol[i] == wfit) {
            return(i)
        }
    }
    return(0)
}

real scalar dwlsmv_m_find(string colvector pool, string scalar want)
{
    real scalar i

    for (i = 1; i <= rows(pool); i++) {
        if (pool[i] == want) {
            return(i)
        }
    }
    return(0)
}

string scalar dwlsmv_eq_spec(string scalar link, string scalar idm)
{
    string scalar s

    s = "factor f1 f2" + char(10) + "ind y1 y2 y3 y4 y5 y6" + char(10)
    if (idm == "unitvar") {
        s = s + "lambda f1 y1 free" + char(10) + "lambda f2 y4 free" + char(10)
    }
    s = s + "lambda f1 y2 free" + char(10) + "lambda f1 y3 free" + char(10) +
            "lambda f2 y5 free" + char(10) + "lambda f2 y6 free" + char(10)
    if (link == "cor") {
        s = s + "psi f1 f2 free" + char(10)
    }
    else {
        s = s + "beta f2 f1 free" + char(10)
    }
    if (idm == "unitvar") {
        s = s + "metric unitvariance f1" + char(10) + "metric unitvariance f2" + char(10)
    }
    else {
        s = s + "psi f1 f1 free" + char(10) + "psi f2 f2 free" + char(10) +
                "metric marker f1 y1" + char(10) + "metric marker f2 y4" + char(10)
    }
    return(s)
}

// Checks 12-15. Fits the correlated and the F2 <- F1 parameterization of case
// E under each metric with dwls_fit() and compares the two with each other.
//
// TOLERANCES. npar and df are integers and compared exactly, with each other
// and with lavaan's. sigma(theta-hat) is compared at spec 7.2's theta-hat tier
// (1e-5 abs, 1e-4 rel), since its entries live on the parameters' scale. T_raw
// is (N - 1) F-hat, so it takes the DWLS-objective tier, 1e-8 relative. Spec 4.8
// quotes the test-statistic tier for this identity because what caps it
// is two separately converged optimizations; dwlsmv_minimize()'s tolerances
// put ours far inside both tiers, and the measured gaps are printed so that
// stays visible.
//
// T3 (check 15, added 2026-09-24 with Module 9). Between our own two fits,
// 1e-8 relative with a 1e-8 absolute term: T_raw agrees at the objective tier
// above, and a and b are functions of Gamma and of Delta's column space, which
// the two parameterizations share up to Delta's central-difference error
// (~1e-10 per entry, dwlsmv_delta_cols' header). Against lavaan's
// chisq_scaled, spec 7.2's test-statistic tier: 1e-6 absolute + 1e-5
// relative, as test_inference.do check 3; T3_df against df_scaled exactly.
string colvector dwlsmv_equiv_ours(struct dwlsmv_stage1 scalar st,
                                   string colvector EFIT, real colvector EDF,
                                   real colvector ENP, real colvector ET3,
                                   real colvector EDF3)
{
    struct dwlsmv_parmap scalar pa, pb
    struct dwlsmv_fit    scalar fa, fb
    string colvector FAIL
    string rowvector idm
    real colvector sa, sb
    real scalar i, ia, ib, rs, rt, r3, la, lb

    FAIL = J(0, 1, "")
    idm  = ("unitvar", "marker")

    for (i = 1; i <= cols(idm); i++) {
        pa = dwls_identify(dwlsmv_eq_spec("cor", idm[i]), st.map)
        pb = dwls_identify(dwlsmv_eq_spec("str", idm[i]), st.map)
        fa = dwls_fit(st, pa)
        fb = dwls_fit(st, pb)
        ia = dwlsmv_m_row(EFIT, "cor_" + idm[i])

        // ---- 12. both converged ------------------------------------------
        if (fa.converged != 1 | fb.converged != 1) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": a fit did not converge")
            continue
        }

        // ---- 13. same npar and df, and lavaan's ---------------------------
        printf("{txt} 13. %-7s ours npar/df      : cor %g/%g, str %g/%g, lavaan %g/%g\n",
               idm[i], pa.npar, fa.D_stat, pb.npar, fb.D_stat,
               (ia ? ENP[ia] : .), (ia ? EDF[ia] : .))
        if (pa.npar != pb.npar | fa.D_stat != fb.D_stat) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": npar or df differ")
        }
        if (ia == 0) {
            FAIL = FAIL \ ("equivalence fixture is missing cor_" + idm[i])
        }
        else if (pa.npar != ENP[ia] | fa.D_stat != EDF[ia]) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": npar or df differ from lavaan's")
        }

        // ---- 14. same sigma(theta-hat) and same T_raw ----------------------
        sa = dwls_model(fa.theta, pa, st.map)
        sb = dwls_model(fb.theta, pb, st.map)
        rs = dwlsmv_gapratio(sa, sb, 1e-5, 1e-4)
        rt = abs(fa.T_raw - fb.T_raw) / abs(fb.T_raw)
        printf("{txt} 14. %-7s ours sigma, T_raw : sigma ratio %8.2e of the theta-hat budget ; T_raw %18.12f vs %18.12f, rel %8.2e\n",
               idm[i], rs, fa.T_raw, fb.T_raw, rt)
        if (rs > 1) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": sigma(theta-hat) differs, ratio " +
                           strofreal(rs))
        }
        if (rt > 1e-8) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": T_raw differs, rel " + strofreal(rt))
        }

        // ---- 15. same T3, and lavaan's ---------------------------------------
        dwls_inference(fa)
        dwls_inference(fb)
        ib = dwlsmv_m_row(EFIT, "str_" + idm[i])
        r3 = dwlsmv_gapratio(fa.T3, fb.T3, 1e-8, 1e-8)
        la = (ia ? dwlsmv_gapratio(fa.T3, ET3[ia], 1e-6, 1e-5) : .)
        lb = (ib ? dwlsmv_gapratio(fb.T3, ET3[ib], 1e-6, 1e-5) : .)
        printf("{txt} 15. %-7s ours T3           : cor %16.12f, str %16.12f, ratio %8.2e of 1e-8 ; vs lavaan %8.2e, %8.2e of the T3 tier ; df %g/%g\n",
               idm[i], fa.T3, fb.T3, r3, la, lb, fa.T3_df, fb.T3_df)
        if (r3 > 1 | missing(r3)) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": T3 differs, ratio " + strofreal(r3))
        }
        if (la > 1 | lb > 1 | missing((la, lb))) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": T3 outside the test-statistic tier of lavaan's")
        }
        else if (fa.T3_df != EDF3[ia] | fb.T3_df != EDF3[ib]) {
            FAIL = FAIL \ ("our equivalence, " + idm[i] + ": T3_df is not lavaan's df_scaled")
        }
    }

    return(FAIL)
}

real scalar dwlsmv_structural_metric(
    real matrix Y, real colvector K,
    string colvector MFIT, string colvector MQTY, real colvector MVAL,
    string colvector EFIT, real colvector ECHI, real colvector EDF,
    real colvector ET3, real colvector EDF3, real colvector ENP,
    real matrix G_POP, string colvector LAB_POP,
    real colvector TAU_POP, real colvector Q_POP)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_statmap scalar m
    string colvector FAIL
    string rowvector idm
    string scalar a_tag, b_tag
    real scalar psi_free, psi_est, beta, var2, var1, var2_cor
    real scalar i, ia, ib, rel, j, c, d, s, t, r_within, nw
    real scalar sym, mineig, maxz, z, se, sa, sb, sab, qa, qb, hit, hjt, N
    real colvector n, q, qh, phi
    real scalar p

    FAIL = J(0, 1, "")
    N    = rows(Y)
    p    = rows(K)

    // =====================================================================
    // A. the metric question
    // =====================================================================

    psi_free = dwlsmv_m_get(MFIT, MQTY, MVAL, "str_unitvar", "psi_f2_f2_is_free")
    psi_est  = dwlsmv_m_get(MFIT, MQTY, MVAL, "str_unitvar", "psi_f2_f2_est")
    beta     = dwlsmv_m_get(MFIT, MQTY, MVAL, "str_unitvar", "beta_f2_f1")
    var2     = dwlsmv_m_get(MFIT, MQTY, MVAL, "str_unitvar", "var_f2_implied")
    var1     = dwlsmv_m_get(MFIT, MQTY, MVAL, "str_unitvar", "var_f1_implied")
    var2_cor = dwlsmv_m_get(MFIT, MQTY, MVAL, "cor_unitvar", "var_f2_implied")

    if (missing((psi_free, psi_est, beta, var2, var1, var2_cor))) {
        FAIL = FAIL \ "metric fixture is missing a (fit, quantity) row"
        printf("{txt}test_structural_metric: %f failure(s)\n", rows(FAIL))
        printf("{err}    %s\n", FAIL[1])
        return(rows(FAIL))
    }

    // ---- 1. something IS fixed for the endogenous factor -------------------
    printf("{txt}  1. psi(f2,f2) is fixed      : free = %g, estimate = %g\n",
           psi_free, psi_est)
    if (psi_free != 0) {
        FAIL = FAIL \ "psi(f2,f2) is free under std.lv -- the metric is not what 4.8 documents"
    }
    if (dwlsmv_gapratio(psi_est, 1, 1e-12, 0) > 1) {
        FAIL = FAIL \ ("psi(f2,f2) is fixed to " + strofreal(psi_est) + ", not to 1")
    }

    // ---- 2. what is fixed is the DISTURBANCE: Var(F2) = beta^2 + psi -------
    // Exact by the 4.6 mapping Phi = (I-B)^-1 Psi (I-B)^-T, so 1e-12 with no
    // relative term. A failure here is arithmetic, not convention.
    printf("{txt}  2. Var(F2) = beta^2 + psi   : %20.15g vs %20.15g\n",
           var2, beta * beta + psi_est)
    if (dwlsmv_gapratio(var2, beta * beta + psi_est, 1e-12, 0) > 1) {
        FAIL = FAIL \ "implied Var(F2) is not beta^2 + psi(f2,f2)"
    }

    // ---- 3. THE DISCRIMINATING CHECK: the TOTAL variance is not 1 ----------
    // Checks 1 and 2 are consistent with either candidate answer. This is the
    // one that separates them, and it is deliberately a coarse threshold: the
    // two answers differ by beta^2 = 0.28 here, so there is nothing marginal
    // about it and a tight tolerance would only invite a tolerance argument.
    printf("{txt}  3. total Var(F2) is NOT 1   : %20.15g (beta = %20.15g)\n", var2, beta)
    if (abs(var2 - 1) <= 0.1) {
        FAIL = FAIL \ ("total Var(F2) is " + strofreal(var2) +
                       " -- lavaan now fixes the TOTAL variance; spec 4.8 is stale")
    }
    if (abs(beta) <= 0.1) {
        FAIL = FAIL \ "beta is near zero, so this fixture cannot separate the two answers"
    }

    // ---- 4/5. why no earlier case could have settled this ------------------
    // For the EXOGENOUS factor of the same fit, and for both factors of the
    // correlated fit, the disturbance and the total variance coincide. Every
    // case before E is of that kind.
    printf("{txt}  4. exogenous Var(F1) = 1    : %20.15g\n", var1)
    if (dwlsmv_gapratio(var1, 1, 1e-12, 0) > 1) {
        FAIL = FAIL \ "exogenous Var(F1) is not 1 under std.lv"
    }
    printf("{txt}  5. correlated Var(F2) = 1   : %20.15g\n", var2_cor)
    if (dwlsmv_gapratio(var2_cor, 1, 1e-12, 0) > 1) {
        FAIL = FAIL \ "Var(F2) is not 1 in the correlated fit, where f2 is exogenous"
    }

    // =====================================================================
    // B. the equivalent-models identity, under BOTH identification methods
    // =====================================================================

    idm = ("unitvar", "marker")
    for (i = 1; i <= cols(idm); i++) {
        a_tag = "cor_" + idm[i]
        b_tag = "str_" + idm[i]
        ia    = dwlsmv_m_row(EFIT, a_tag)
        ib    = dwlsmv_m_row(EFIT, b_tag)

        if (ia == 0 | ib == 0) {
            FAIL = FAIL \ ("equivalence fixture is missing " + a_tag + " or " + b_tag)
            continue
        }

        // ---- 6. same free-parameter count and same df ----------------------
        // The reparameterization spends the same budget. If these differ, the
        // two models are not equivalent and nothing below is meaningful.
        if (ENP[ia] != ENP[ib] | EDF[ia] != EDF[ib] | EDF3[ia] != EDF3[ib]) {
            FAIL = FAIL \ ("equivalence, " + idm[i] +
                           ": npar or df differ -- the models are not reparameterizations")
            continue
        }

        // ---- 7. same unscaled chi-square -----------------------------------
        rel = abs(ECHI[ia] - ECHI[ib]) / abs(ECHI[ia])
        printf("{txt}  7. %-7s chisq            : %18.12f vs %18.12f, rel %10.2e, df %g, npar %g\n",
               idm[i], ECHI[ia], ECHI[ib], rel, EDF[ia], ENP[ia])
        if (dwlsmv_gapratio(ECHI[ia], ECHI[ib], 1e-6, 1e-5) > 1) {
            FAIL = FAIL \ ("equivalence, " + idm[i] + ": chisq relative gap " +
                           strofreal(rel))
        }

        // ---- 8. same T3, at spec 7.2's test-statistic tier -----------------
        // 1e-6 absolute + 1e-5 relative since 2026-09-24 (it was 1e-3
        // relative). What caps the agreement is that these are two separately
        // converged optimizations; the reparameterization itself is exact.
        // The measured gap is ~1e-09 relative, four orders inside the tier,
        // so the tier is not doing any work here -- which is the point of
        // quoting it rather than a number tuned to what was observed.
        rel = abs(ET3[ia] - ET3[ib]) / abs(ET3[ia])
        printf("{txt}  8. %-7s T3               : %18.12f vs %18.12f, rel %10.2e, df %g\n",
               idm[i], ET3[ia], ET3[ib], rel, EDF3[ia])
        if (dwlsmv_gapratio(ET3[ia], ET3[ib], 1e-6, 1e-5) > 1) {
            FAIL = FAIL \ ("equivalence, " + idm[i] + ": T3 relative gap " +
                           strofreal(rel))
        }
    }

    // =====================================================================
    // C. Stage 1 on case E, before anything above is trusted
    // =====================================================================

    st = dwls_sample_stats(Y, K)
    m  = st.map

    // ---- 9. structural checks that need no reference ----------------------
    sym    = dwlsmv_symmetry(st.Gamma)
    mineig = dwlsmv_mineig(st.Gamma)
    printf("{txt}  9. Gamma structure          : symmetry %g, min eigenvalue %g\n",
           sym, mineig)
    if (sym != 0) {
        FAIL = FAIL \ ("Gamma is not symmetric to 0.0 on case E: " + strofreal(sym))
    }
    if (mineig <= 0) {
        FAIL = FAIL \ ("Gamma is not positive definite on case E: min eigenvalue " +
                       strofreal(mineig))
    }

    // ---- 10. the within-variable threshold identity -----------------------
    // Exact arithmetic at any N (spec 4.5.1), so 1e-12 absolute with no
    // relative term. First time it is exercised on a six-indicator design.
    qh  = J(m.T, 1, .)
    phi = J(m.T, 1, .)
    for (j = 1; j <= p; j++) {
        n = dwlsmv_catcounts(Y[., j], K[j])
        q = dwlsmv_cumprop(n)
        for (c = 1; c <= m.C[j]; c++) {
            s      = dwlsmv_pos_tau(m, j, c)
            qh[s]  = q[c]
            phi[s] = normalden(invnormal(q[c]))
        }
    }

    r_within = 0
    nw       = 0
    for (j = 1; j <= p; j++) {
        for (c = 1; c <= m.C[j]; c++) {
            for (d = 1; d <= m.C[j]; d++) {
                s        = dwlsmv_pos_tau(m, j, c)
                t        = dwlsmv_pos_tau(m, j, d)
                r_within = max((r_within,
                                dwlsmv_gapratio(st.Gamma[s, t],
                                                (min((qh[s], qh[t])) - qh[s] * qh[t]) /
                                                (phi[s] * phi[t]), 1e-12, 0)))
                nw = nw + 1
            }
        }
    }
    printf("{txt} 10. within-variable identity : ratio %g over %f entries (1e-12 abs)\n",
           r_within, nw)
    if (r_within > 1) {
        FAIL = FAIL \ ("case E within-variable identity: ratio " + strofreal(r_within))
    }

    // ---- 11. bias against the tier-1 population block ----------------------
    // NOT an equality check. case_E_gamma_tau.csv holds the POPULATION
    // threshold block, built from the declared category probabilities and
    // Lambda Phi Lambda'; Gamma-hat differs from it by sampling error by
    // construction. Judged as a z statistic against each entry's own sampling
    // error, exactly as test_gamma_tau.do check 3 does for case D.
    maxz = 0
    for (s = 1; s <= m.T; s++) {
        hit = dwlsmv_m_find(LAB_POP, m.label[s])
        if (hit == 0) {
            FAIL = FAIL \ ("bias: tier 1 has no entry labelled " + m.label[s])
            continue
        }
        for (t = 1; t <= m.T; t++) {
            hjt = dwlsmv_m_find(LAB_POP, m.label[t])
            if (hjt == 0) {
                continue
            }
            qa  = Q_POP[hit]
            qb  = Q_POP[hjt]
            sa  = qa * (1 - qa)
            sb  = qb * (1 - qb)
            sab = G_POP[hit, hjt] * normalden(TAU_POP[hit]) * normalden(TAU_POP[hjt])
            se  = sqrt((sa * sb + sab * sab) / N) /
                  (normalden(TAU_POP[hit]) * normalden(TAU_POP[hjt]))
            z    = (st.Gamma[s, t] - G_POP[hit, hjt]) / se
            maxz = max((maxz, abs(z)))
        }
    }
    printf("{txt} 11. tier-1 bias              : worst |z| %g against each entry's own SE\n",
           maxz)
    if (maxz > 5) {
        FAIL = FAIL \ ("case E bias: worst |z| is " + strofreal(maxz))
    }

    // =====================================================================
    // D. the equivalent-models identity on OUR fits
    // =====================================================================
    FAIL = FAIL \ dwlsmv_equiv_ours(st, EFIT, EDF, ENP, ET3, EDF3)

    printf("{txt}test_structural_metric: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_structural_metric(                 ///
          YE, (3 \ 3 \ 3 \ 3 \ 3 \ 3),                                      ///
          MFIT, MQTY, MVAL, EFIT, ECHI, EDF, ET3, EDF3, ENP,                ///
          G_POP, LAB_POP, TAU_POP, Q_POP)))

if `nfail' {
    display as error _n "test_structural_metric: `nfail' failure(s)"
    exit 9
}

display as result "test_structural_metric: endogenous metric is the DISTURBANCE variance; correlated and structural parameterizations agree on npar, df, chisq and T3 in lavaan, and on npar, df, sigma(theta-hat), T_raw and T3 in our own fits; Stage 1 verified on case E"
exit 0
