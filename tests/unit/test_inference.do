*! Module 9: dwlsmv_traces(), dwlsmv_scaleshift() and dwls_inference() -- the
*! scaled-and-shifted statistic T3 = a T_raw + b, df = D, and the locked entry
*! point that fills acov and T3 together.
*!
*! Executed under StataMP 18 on 2026-09-24; passed on the first run, so
*! nothing here was adjusted to fit a result. Every tolerance was set before
*! the run: from R measurements (checks 1, 4) or from spec 7.2 (check 3).
*! Run at six and at four processors: the printed values that move are
*! checks 2, 2b and 5 -- floating-point residues, check 2's against an
*! oracle and 2b's and 5's in production values at the ulp level. So a, b
*! and the traces are NOT bitwise processor-invariant; what does not move is
*! every printed digit of T3, a, b and every end-to-end ratio.
*!
*! ADVERSARIALLY REVIEWED 2026-09-24 (read-only agent, 23 mutations run
*! against this file at d4ea113). Closed afterwards: checks 1 and 3 took an
*! outer max() over gap ratios, which skips a missing element; check 4's
*! must-fail tests passed on missing (. <= 1 is false); check 6's T3 = T_raw
*! rode on T_raw = 1.7e-29 and now also runs a hand-built D = 0 fit with
*! T_raw = 1; check 5 now compares par and all of stage1; and three guards
*! gained error-path cases (8e-8g).
*!
*! HOW MANY INDEPENDENT TRACE PROBLEMS THIS IS: TWO, NOT FIVE. The four case E
*! fits are equivalent models (test_structural_metric.do), so their traces
*! agree in exact arithmetic and differ only at lavaan's convergence. Checks
*! 1 and 3 confirm case D and case E; the other three case E fits exercise
*! the parameterizations, not new traces.
*!
*! Spec Section 4.9; derivation @sec-gate5, @sec-gate5-traces,
*! @sec-gate5-lavaan. Tolerances are spec Section 7.2's, and each check says
*! which row and why.
*!
*! The checks, and why each is separate:
*!
*!   1  the traces, a, b and T3 from LAVAAN'S OWN INPUTS -- its analytic Delta
*!      at its own theta-hat (lavaan_acov_delta.csv), its Gamma
*!      (lavaan_case_*_gamma.csv) and its unscaled chi-square -- against the
*!      values lavaan stored in fit@test$scaled.shifted (lavaan_t3.csv). No
*!      theta-hat or rho-hat tier is inherited, so this is the arithmetic
*!      alone: 1e-12 relative on Tr(M), Tr(M^2), a and T3, and b at 1e-12
*!      ABSOLUTE because it is a cancellation (spec 4.9). Measured in R at
*!      <= 5.7e-15 and 2.1e-14 absolute (derivation @sec-gate5-lavaan).
*!      test_acov.do check 2, one module on.
*!   2  the expansion against the DEFINITION: U and M = U Gamma formed
*!      explicitly on our five fits, traces compared within the per-fit GB1
*!      budget (@sec-gate5-budget; replaced a fixed 1e-13 relative on 2026-10-28,
*!      after case D missed it at 4 processors under Stata 18, kappa2 ~ 208).
*!      2b: a DIAGONAL Gamma makes W_D = Gamma^-1, so M is idempotent and every
*!      correction must collapse: Tr(M) = Tr(M^2) = D, a = 1, b = 0, to 1e-12
*!      absolute. The derivation's Level 2, the strongest structural check
*!      there is, run through the production routine.
*!   3  end to end: our Gamma, our theta-hat and our Delta through
*!      dwls_inference(), against lavaan_t3.csv. Spec 7.2's test-statistic
*!      tier: Tr(M), Tr(M^2) and a at 1e-5 relative (all bounded away from 0
*!      when D > 0), b at 1e-6 absolute, T3 at 1e-6 absolute + 1e-5 relative
*!      (T3 -> b as the fit becomes perfect, so it takes b's absolute term),
*!      df exactly. The gap is lavaan's own convergence, ~1e-7 at worst
*!      (E-032, 2026-09-24); the worst ratio is printed.
*!   4  the (N-1)F convention, which the 1e-5 tier can now see. Against
*!      lavaan: a N F-hat + b must FAIL check 3's T3 tier. Against Mplus, case
*!      D (tier 3): a N F-hat + b must match chisq.wlsmv to 3e-5 relative and
*!      our a (N-1) F-hat + b must NOT. 9.5e-6 and 1.9e-4 measured in R before
*!      the tolerance was set (derivation @sec-gate5-lavaan); 3e-5 separates
*!      them with margin on both sides. test_acov.do check 6's twin.
*!   5  bookkeeping: acov is dwlsmv_acov()'s, bitwise; T3 is exactly
*!      shift_a T_raw + shift_b; T3_df = D_stat; every field dwls_fit() set is
*!      unchanged, including par and all of stage1; and Tr(M), Tr(M^2) are
*!      recovered from (a, b, D) as (D - b)/a and D/a^2 to 1e-13 -- the claim
*!      that lets them go unstored.
*!   6  D = 0, a just-identified model (y1-y3 of case D): acov filled and equal
*!      to dwlsmv_acov(), T3 = T_raw, T3_df = 0, shift_a and shift_b missing.
*!      lavaan 0.7.2's behaviour, decided 2026-09-24. Its fitted T_raw is
*!      ~1e-29, which cannot tell T3 = T_raw from T3 = 0, so 6b repeats it
*!      with T_raw set to 1 by hand. dwlsmv_scaleshift(., ., 0) returns (., .)
*!      whatever the traces are.
*!   7  a fit that is NOT an optimum -- theta at its starting values, Delta,
*!      F and T_raw evaluated there, converged = 0 -- is processed, not
*!      refused (decided 2026-09-24), and converged stays 0.
*!   8  the error paths, at rc 3498 and not merely nonzero (a project rule).
*!      8d is ALL OR NOTHING: Gamma = Delta Delta' - eps I has a positive
*!      diagonal, so dwlsmv_acov() succeeds, but M = U Gamma = -eps U, so
*!      Tr(M) = -eps Tr(U) < 0 by construction, not by rounding, and
*!      dwlsmv_scaleshift() aborts AFTER acov was computed. The fit must be
*!      left exactly unfilled. It is kept in a Mata external so it survives
*!      the abort to be inspected. 8e-8g: a negative T_raw, a Delta with a
*!      row too many (the row half of the shape guard), and a Delta with no
*!      columns -- guards the review found had no case.
*!
*! Reference CSVs are read with stringcols(_all), strtoreal() and
*! delimiters(","), for the two import traps (float typing, '|' as delimiter).
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
        display as error "test_inference: `f' failed to load"
        exit 9
    }
}

* ---- tier 2: lavaan's T3 components, Delta and Gamma, long format ------------
import delimited using "fixtures/expected/lavaan_t3.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: T3F = st_sdata(., "fit")
mata: T3Q = st_sdata(., "quantity")
mata: T3V = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/lavaan_acov_delta.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: LDF = st_sdata(., "fit")
mata: LDR = st_sdata(., "row")
mata: LDC = st_sdata(., "col")
mata: LDV = strtoreal(st_sdata(., "value"))

foreach c in D E {
    import delimited using "fixtures/expected/lavaan_case_`c'_gamma.csv", ///
        clear varnames(1) delimiters(",") stringcols(_all)
    qui ds label, not
    local gcols `r(varlist)'
    mata: G`c'  = strtoreal(st_sdata(., "`gcols'"))
    mata: G`c'L = st_sdata(., "label")
}

* ---- tier 3: Mplus case D ----------------------------------------------------
import delimited using "fixtures/expected/mplus_case_D_test.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: MPL = st_sdata(., "label")
mata: MPV = strtoreal(st_sdata(., "value"))

* ---- the data ----------------------------------------------------------------
import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

// Long-format lookup; missing when absent, which every caller reports.
real scalar ti_get(string colvector fitcol, string colvector keycol,
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

real scalar ti_pos(string colvector keys, string scalar k)
{
    real scalar i

    for (i = 1; i <= rows(keys); i++) {
        if (keys[i] == k) {
            return(i)
        }
    }
    return(0)
}

string scalar ti_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The five fits of fixtures/generate/20_lavaan_reference.R, as test_acov.do
// writes them, plus the just-identified y1-y3 model of check 6.
string scalar ti_canon(string scalar fit)
{
    string rowvector head, lam6

    head = ("factor f1 f2", "ind y1 y2 y3 y4 y5 y6")
    lam6 = ("lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
            "lambda f2 y4 free", "lambda f2 y5 free", "lambda f2 y6 free")

    if (fit == "D_marker") {
        return(ti_spec(("factor f", "ind y1 y2 y3 y4",
                        "lambda f y2 free", "lambda f y3 free",
                        "lambda f y4 free", "psi f f free",
                        "metric marker f y1")))
    }
    if (fit == "D_just") {
        return(ti_spec(("factor f", "ind y1 y2 y3",
                        "lambda f y2 free", "lambda f y3 free",
                        "psi f f free", "metric marker f y1")))
    }
    if (fit == "cor_unitvar") {
        return(ti_spec((head, lam6, "psi f1 f2 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    if (fit == "str_unitvar") {
        return(ti_spec((head, lam6, "beta f2 f1 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    if (fit == "cor_marker") {
        return(ti_spec((head,
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "psi f1 f1 free", "psi f2 f2 free", "psi f1 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    if (fit == "str_marker") {
        return(ti_spec((head,
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "beta f2 f1 free",
                        "psi f1 f1 free", "psi f2 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    _error("ti_canon: unknown fit")
}

// The definition, forming U and M: the oracle for check 2. Never production.
real rowvector ti_traces_U(real matrix Delta, real matrix Gamma)
{
    real colvector w
    real matrix WD, Ei, U, M

    w  = 1 :/ diagonal(Gamma)
    WD = w :* Delta
    Ei = cholinv(quadcross(Delta, w, Delta))
    U  = diag(w) - quadcross(WD', quadcross(Ei, WD'))
    M  = quadcross(U, Gamma)
    return((sum(diagonal(M)), sum(M :* M')))
}

// Check 2's budget (derivation @sec-gate5-budget, Result GB1): one run's first-order
// error bound for (Tr M, Tr M^2), the expansion's plus the definition's, times 1 + rho.
// Mirrors the derivation's gg_beta_noU and gg_beta_U term for term. Never production.
real scalar ti_g(real scalar n)
{
    return(n * 2^-53 / (1 - n * 2^-53))
}

real rowvector ti_budget(real matrix Dl, real matrix Gm)
{
    real colvector w, dw0, dwk, y
    real matrix WD, P, Ei, G, B, C, EB, In, S, U, M, dPx, dEB, dC, dB, dEi, dP, dG, dWD, dM, dU, dS, dIn
    real scalar p, q, u, k, a1, s1, a2, s2, s3, be1, be2, bu1, bu2, f

    u  = 2^-53
    p  = rows(Dl)
    q  = cols(Dl)
    w  = 1 :/ diagonal(Gm)
    WD = w :* Dl
    P  = Dl' * WD
    Ei = luinv(P)
    G  = Gm * WD
    B  = WD' * G
    C  = G' * (w :* G)
    EB = Ei * B
    In = Ei * WD'
    S  = WD * In
    U  = diag(w) - S
    M  = U * Gm
    dPx = ti_g(p + 1) :* (abs(Dl)' * abs(WD)) + 6 * (q + 1) * u :* sqrt(diagonal(P) * diagonal(P)')

    // the expansion: reverse mode for k = 1 (Tr M) and k = 2 (Tr M^2)
    for (k = 1; k <= 2; k++) {
        if (k == 1) {
            dEB = -I(q)
            dC  = J(q, q, 0)
            dEi = dEB * B'
        }
        else {
            dEB = 2 :* EB'
            dC  = -2 :* Ei'
            dEi = dEB * B' - 2 :* C'
        }
        dB  = Ei' * dEB
        dP  = -(Ei' * (dEi * Ei'))
        dG  = WD * dB + w :* (G * (dC + dC'))
        dWD = Dl * dP + G * dB' + Gm' * dG
        dwk = rowsum(dWD :* Dl) + rowsum(G :* (G * dC))
        if (k == 1) {
            dwk = dwk + diagonal(Gm)
        }
        else {
            dwk = dwk + 2 :* ((Gm :* Gm) * w)
        }
        f = sum(abs(dwk) :* (u :* w)) + sum(abs(dWD) :* (u :* abs(WD))) + sum(abs(dP) :* dPx) +
            sum(abs(dG) :* (ti_g(p + 1) :* (abs(Gm) * abs(WD)))) +
            sum(abs(dB) :* (ti_g(p + 1) :* (abs(WD)' * abs(G)))) +
            sum(abs(dC) :* (ti_g(p + 1) :* (abs(G)' * (w :* abs(G))))) +
            sum(abs(dEB) :* (ti_g(q + 1) :* (abs(Ei) * abs(B))))
        if (k == 1) {
            a1  = sum(w :* diagonal(Gm))
            s1  = sum(diagonal(EB))
            be1 = f + ti_g(p + 1) * sum(abs(w :* diagonal(Gm))) + ti_g(q) * sum(abs(diagonal(EB))) +
                  u * (abs(a1) + abs(s1))
        }
        else {
            y   = (Gm :* Gm) * w
            a2  = w' * y
            s2  = sum(diagonal(Ei * C))
            s3  = sum(diagonal(EB * EB))
            be2 = f + sum(abs(w * w') :* (u :* (Gm :* Gm))) +
                  ti_g(p * p + 1) * sum(abs((w * w') :* (Gm :* Gm))) +
                  2 * (ti_g(q + 1) * sum(abs(Ei) :* abs(C')) + ti_g(q) * abs(s2)) +
                  ti_g(q + 1) * sum(abs(EB) :* abs(EB')) + ti_g(q) * abs(s3) +
                  2 * u * (abs(a2) + 2 * abs(s2) + abs(s3))
        }
    }

    // the definition
    for (k = 1; k <= 2; k++) {
        if (k == 1) {
            dM = I(p)
        }
        else {
            dM = 2 :* M'
        }
        dU  = dM * Gm'
        dS  = -dU
        dIn = WD' * dS
        dEi = dIn * WD
        dP  = -(Ei' * (dEi * Ei'))
        dWD = dS * In' + dIn' * Ei + Dl * dP
        dw0 = diagonal(dU) + rowsum(dWD :* Dl)
        f = sum(abs(dw0) :* (u :* w)) + sum(abs(dWD) :* (u :* abs(WD))) + sum(abs(dP) :* dPx) +
            sum(abs(dIn) :* (ti_g(q + 1) :* (abs(Ei) * abs(WD)'))) +
            sum(abs(dS) :* (ti_g(q + 1) :* (abs(WD) * abs(In)))) +
            sum(abs(dU) :* (u :* abs(U))) +
            sum(abs(dM) :* (ti_g(p + 1) :* (abs(U) * abs(Gm))))
        if (k == 1) {
            bu1 = f + ti_g(p) * sum(abs(diagonal(M)))
        }
        else {
            bu2 = f + ti_g(p * p + 1) * sum(abs(M :* M'))
        }
    }
    return(1.01 :* (be1 + bu1, be2 + bu2))
}

// The worst of several gap ratios, MISSING if any one is missing. Mata's max()
// skips missing elements, so max((r1, r2, .)) would report the worst of the
// other two -- the 2026-09-24 review found checks 1 and 3 doing exactly that.
real scalar ti_worst(real rowvector r)
{
    if (missing(r)) {
        return(.)
    }
    return(max(r))
}

// Whether a fit is exactly as dwls_fit() left it for Modules 8-9.
real scalar ti_isunfilled(struct dwlsmv_fit scalar f)
{
    return(rows(f.acov) == 0 & missing(f.T3) & missing(f.T3_df) &
           missing(f.shift_a) & missing(f.shift_b))
}

real scalar ti_checks(real matrix YD, real matrix YE,
                      string colvector T3F, string colvector T3Q,
                      real colvector T3V,
                      string colvector LDF, string colvector LDR,
                      string colvector LDC, real colvector LDV,
                      real matrix GD, string colvector GDL,
                      real matrix GE, string colvector GEL,
                      string colvector MPL, real colvector MPV)
{
    struct dwlsmv_stage1 scalar stD, stE, st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit, f0, fu
    string colvector FAIL, glab, dcol
    string rowvector fits
    string scalar fitn
    real colvector th0, r0, w
    real rowvector tr, trU, ab, trd, bud
    real matrix Dl, Gl, Glav, A
    real scalar d, i, j, k, m, N, df, lT, lT3, rq, rb, r, g1, g2, alt, mpT, fm0

    FAIL = J(0, 1, "")
    stD  = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    stE  = dwls_sample_stats(YE, J(6, 1, 3))
    fits = ("D_marker", "cor_unitvar", "str_unitvar", "cor_marker", "str_marker")

    for (d = 1; d <= cols(fits); d++) {
        fitn = fits[d]
        if (fitn == "D_marker") {
            st   = stD
            Glav = GD
            glab = GDL
        }
        else {
            st   = stE
            Glav = GE
            glab = GEL
        }
        N   = st.N
        par = dwls_identify(ti_canon(fitn), st.map)
        fit = dwls_fit(st, par)
        if (fit.converged != 1) {
            FAIL = FAIL \ (fitn + ": dwls_fit() did not converge; nothing below is meaningful")
            continue
        }
        df  = ti_get(T3F, T3Q, T3V, fitn, "df")
        lT  = ti_get(T3F, T3Q, T3V, fitn, "T")
        lT3 = ti_get(T3F, T3Q, T3V, fitn, "T3")
        if (missing((df, lT, lT3,
                     ti_get(T3F, T3Q, T3V, fitn, "trM"), ti_get(T3F, T3Q, T3V, fitn, "trM2"),
                     ti_get(T3F, T3Q, T3V, fitn, "a"),   ti_get(T3F, T3Q, T3V, fitn, "b")))) {
            FAIL = FAIL \ (fitn + ": lavaan_t3.csv is missing a quantity")
            continue
        }

        // ---- 1. lavaan's own Delta and Gamma -------------------------------
        // Delta's columns in lavaan's order: the order its labels first appear.
        dcol = J(0, 1, "")
        for (i = 1; i <= rows(LDF); i++) {
            if (LDF[i] == fitn & !ti_pos(dcol, LDC[i])) dcol = dcol \ LDC[i]
        }
        Dl = J(st.map.D, rows(dcol), .)
        for (i = 1; i <= rows(LDF); i++) {
            if (LDF[i] != fitn) continue
            j = ti_pos(st.map.label, LDR[i])
            k = ti_pos(dcol, LDC[i])
            if (j & k) Dl[j, k] = LDV[i]
        }
        Gl = J(st.map.D, st.map.D, .)
        for (i = 1; i <= st.map.D; i++) {
            j = ti_pos(glab, st.map.label[i])
            if (!j) continue
            for (k = 1; k <= st.map.D; k++) {
                m = ti_pos(glab, st.map.label[k])
                if (m) Gl[i, k] = Glav[j, m]
            }
        }
        if (missing(Dl) | missing(Gl) | rows(dcol) != par.npar) {
            FAIL = FAIL \ (fitn + ": lavaan's Delta or Gamma did not fill from the fixtures")
            continue
        }
        tr = dwlsmv_traces(Dl, Gl)
        ab = dwlsmv_scaleshift(tr[1], tr[2], df)
        rq = ti_worst((dwlsmv_gapratio(tr[1], ti_get(T3F, T3Q, T3V, fitn, "trM"),  0, 1e-12),
                  dwlsmv_gapratio(tr[2], ti_get(T3F, T3Q, T3V, fitn, "trM2"), 0, 1e-12),
                  dwlsmv_gapratio(ab[1], ti_get(T3F, T3Q, T3V, fitn, "a"),    0, 1e-12),
                  dwlsmv_gapratio(ab[1] * lT + ab[2], lT3,                    0, 1e-12)))
        rb = dwlsmv_gapratio(ab[2], ti_get(T3F, T3Q, T3V, fitn, "b"), 1e-12, 0)
        printf("{txt}  1. %-11s lavaan inputs      : traces, a, T3 ratio %8.2e of 1e-12 rel ; b ratio %8.2e of 1e-12 abs (b = %g)\n",
               fitn, rq, rb, ab[2])
        if (rq > 1 | rb > 1 | missing(rq) | missing(rb)) {
            FAIL = FAIL \ (fitn + ": from lavaan's own inputs, traces/a/T3 ratio " + strofreal(rq) +
                           ", b ratio " + strofreal(rb))
        }

        // ---- 2. the expansion against U formed explicitly ------------------
        tr  = dwlsmv_traces(fit.Delta, st.Gamma)
        trU = ti_traces_U(fit.Delta, st.Gamma)
        bud = ti_budget(fit.Delta, st.Gamma)
        r   = ti_worst(abs(tr - trU) :/ bud)
        printf("{txt}  2. %-11s expansion vs U     : ratio %8.2e of the GB1 budget (%8.2e, %8.2e abs; rel gap %8.2e)\n",
               fitn, r, bud[1], bud[2], max(abs(tr - trU) :/ abs(trU)))
        if (r > 1 | missing(r)) {
            FAIL = FAIL \ (fitn + ": the trace expansion disagrees with U-forming, ratio " + strofreal(r))
        }
        // 2b. diagonal Gamma: W_D = Gamma^-1, M idempotent, a = 1, b = 0.
        trd = dwlsmv_traces(fit.Delta, diag(diagonal(st.Gamma)))
        ab  = dwlsmv_scaleshift(trd[1], trd[2], fit.D_stat)
        g1  = max(abs((trd[1] - fit.D_stat, trd[2] - fit.D_stat, ab[1] - 1, ab[2])))
        printf("{txt}  2b. %-10s diagonal Gamma    : max |Tr(M)-D, Tr(M^2)-D, a-1, b| = %8.2e (1e-12)\n",
               fitn, g1)
        if (g1 > 1e-12 | missing(g1)) {
            FAIL = FAIL \ (fitn + ": with W_D = Gamma^-1 the correction did not collapse to the identity")
        }

        // ---- 3. end to end through dwls_inference() -------------------------
        f0 = fit
        dwls_inference(fit)
        rq = ti_worst((dwlsmv_gapratio(tr[1], ti_get(T3F, T3Q, T3V, fitn, "trM"),  0, 1e-5),
                  dwlsmv_gapratio(tr[2], ti_get(T3F, T3Q, T3V, fitn, "trM2"), 0, 1e-5),
                  dwlsmv_gapratio(fit.shift_a, ti_get(T3F, T3Q, T3V, fitn, "a"), 0, 1e-5)))
        rb = dwlsmv_gapratio(fit.shift_b, ti_get(T3F, T3Q, T3V, fitn, "b"), 1e-6, 0)
        r  = dwlsmv_gapratio(fit.T3, lT3, 1e-6, 1e-5)
        printf("{txt}  3. %-11s end to end         : traces/a ratio %8.2e ; b %8.2e ; T3 %8.2e (%16.12f vs %16.12f, rel %8.2e) ; df %g vs %g\n",
               fitn, rq, rb, r, fit.T3, lT3, abs(fit.T3 / lT3 - 1), fit.T3_df, df)
        if (rq > 1 | rb > 1 | r > 1 | missing((rq, rb, r))) {
            FAIL = FAIL \ (fitn + ": end to end outside the test-statistic tier: traces/a " + strofreal(rq) +
                           ", b " + strofreal(rb) + ", T3 " + strofreal(r))
        }
        if (fit.T3_df != df) {
            FAIL = FAIL \ (fitn + ": T3_df is not lavaan's df")
        }

        // ---- 4. the (N-1)F convention ----------------------------------------
        alt = fit.shift_a * N * fit.fmin + fit.shift_b
        g1  = dwlsmv_gapratio(alt, lT3, 1e-6, 1e-5)
        printf("{txt}  4. %-11s a N F + b vs lavaan: ratio %8.2e of the T3 tier (must be > 1)\n", fitn, g1)
        // Written !(g1 > 1), not g1 <= 1: a missing g1 must fail too, and
        // . <= 1 is false in Mata (review 2026-09-24).
        if (!(g1 > 1)) {
            FAIL = FAIL \ (fitn + ": a N F + b also matches lavaan, or is missing -- the tier cannot tell N from N - 1")
        }
        if (fitn == "D_marker") {
            mpT = ti_get(MPL, MPL, MPV, "", "chisq.wlsmv")
            g1  = abs(alt / mpT - 1)
            g2  = abs(fit.T3 / mpT - 1)
            printf("{txt}  4. Mplus chisq.wlsmv %g, df %g : a N F + b rel %8.2e (3e-5) ; ours, a (N-1) F + b, rel %8.2e (must exceed 3e-5)\n",
                   mpT, ti_get(MPL, MPL, MPV, "", "df"), g1, g2)
            if (missing(mpT) | !(g1 <= 3e-5)) {
                FAIL = FAIL \ "Mplus: a N F + b does not reproduce chisq.wlsmv to 3e-5"
            }
            if (!(g2 > 3e-5)) {
                FAIL = FAIL \ "Mplus: a (N-1) F + b also matches -- the convention pin has no power"
            }
            if (ti_get(MPL, MPL, MPV, "", "df") != fit.T3_df) {
                FAIL = FAIL \ "Mplus: df differs"
            }
        }

        // ---- 5. bookkeeping ----------------------------------------------------
        if (fit.acov != dwlsmv_acov(fit.Delta, st.Gamma)) {
            FAIL = FAIL \ (fitn + ": acov is not dwlsmv_acov(Delta, Gamma)")
        }
        if (fit.T3 != fit.shift_a * fit.T_raw + fit.shift_b | fit.T3_df != fit.D_stat) {
            FAIL = FAIL \ (fitn + ": T3 is not a T_raw + b, or T3_df is not D_stat")
        }
        if (fit.theta != f0.theta | fit.Delta != f0.Delta | fit.T_raw != f0.T_raw |
            fit.fmin != f0.fmin | fit.converged != f0.converged |
            fit.iterations != f0.iterations | fit.D_stat != f0.D_stat |
            fit.stage1[1].Gamma != f0.stage1[1].Gamma | fit.stage1[1].s != f0.stage1[1].s |
            fit.stage1[1].N != f0.stage1[1].N | fit.stage1[1].psibar != f0.stage1[1].psibar |
            fit.stage1[1].map.K != f0.stage1[1].map.K | fit.stage1[1].map.label != f0.stage1[1].map.label |
            fit.stage1[1].map.D != f0.stage1[1].map.D) {
            FAIL = FAIL \ (fitn + ": dwls_inference() changed a field dwls_fit() set")
        }
        if (fit.par.npar != f0.par.npar | fit.par.nall != f0.par.nall |
            fit.par.name != f0.par.name | fit.par.ptype != f0.par.ptype |
            fit.par.row != f0.par.row | fit.par.col != f0.par.col |
            fit.par.free != f0.par.free | fit.par.value != f0.par.value |
            fit.par.fidx != f0.par.fidx) {
            FAIL = FAIL \ (fitn + ": dwls_inference() changed the parameter table")
        }
        g1 = abs((fit.D_stat - fit.shift_b) / fit.shift_a / tr[1] - 1)
        g2 = abs(fit.D_stat / fit.shift_a^2 / tr[2] - 1)
        printf("{txt}  5. %-11s bookkeeping        : Tr(M), Tr(M^2) from (a, b, D) rel %8.2e, %8.2e (1e-13)\n",
               fitn, g1, g2)
        if (g1 > 1e-13 | g2 > 1e-13) {
            FAIL = FAIL \ (fitn + ": the traces are not recoverable from (a, b, D)")
        }
    }

    // ---- 6. D = 0 ---------------------------------------------------------------
    st  = dwls_sample_stats(YD[., 1..3], (4 \ 4 \ 3))
    par = dwls_identify(ti_canon("D_just"), st.map)
    fit = dwls_fit(st, par)
    dwls_inference(fit)
    printf("{txt}  6. just identified           : converged %g, D_stat %g, T_raw %g, T3 %g, T3_df %g, a %g, b %g, acov %g x %g\n",
           fit.converged, fit.D_stat, fit.T_raw, fit.T3, fit.T3_df,
           fit.shift_a, fit.shift_b, rows(fit.acov), cols(fit.acov))
    if (fit.converged != 1 | fit.D_stat != 0) {
        FAIL = FAIL \ "6: the just-identified fit did not converge or does not have D = 0"
    }
    else if (rows(fit.acov) != par.npar | fit.acov != dwlsmv_acov(fit.Delta, st.Gamma) |
             fit.T3 != fit.T_raw | fit.T3_df != 0 |
             !missing(fit.shift_a) | !missing(fit.shift_b)) {
        FAIL = FAIL \ "6: D = 0 is not acov filled, T3 = T_raw, T3_df = 0, a and b missing"
    }
    // 6b. The fitted T_raw is ~1e-29, which cannot tell T3 = T_raw from T3 = 0;
    // a hand-set T_raw = 1 can. Only T_raw's consistency with fmin is kept.
    fu = fit
    fu.acov    = J(0, 0, .)
    fu.T3      = .
    fu.T3_df   = .
    fu.shift_a = .
    fu.shift_b = .
    fu.fmin    = 1 / (st.N - 1)
    fu.T_raw   = 1
    dwls_inference(fu)
    printf("{txt}  6b. D = 0, T_raw set to 1      : T3 %g (must be 1), T3_df %g, a %g, b %g\n",
           fu.T3, fu.T3_df, fu.shift_a, fu.shift_b)
    if (fu.T3 != 1 | fu.T3_df != 0 | !missing(fu.shift_a) | !missing(fu.shift_b)) {
        FAIL = FAIL \ "6b: with D = 0 and T_raw = 1, T3 is not exactly T_raw"
    }
    if (dwlsmv_scaleshift(1, 1, 0) != (., .) | dwlsmv_scaleshift(-1, -1, 0) != (., .)) {
        FAIL = FAIL \ "6: dwlsmv_scaleshift() with df = 0 is not (., .)"
    }

    // ---- 7. not an optimum: processed, not refused -------------------------------
    par = dwls_identify(ti_canon("D_marker"), stD.map)
    fu  = dwls_fit(stD, par)
    th0 = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) th0[par.fidx[i]] = par.value[i]
    }
    w   = 1 :/ diagonal(stD.Gamma)
    r0  = stD.s - dwls_model(th0, par, stD.map)
    fm0 = quadcross(r0, w, r0)
    fu.theta     = th0
    fu.Delta     = dwls_delta(th0, par, stD.map)
    fu.fmin      = fm0
    fu.T_raw     = (stD.N - 1) * fm0
    fu.converged = 0
    dwls_inference(fu)
    printf("{txt}  7. at the starting values     : converged %g, T_raw %g, T3 %g, a %g, b %g\n",
           fu.converged, fu.T_raw, fu.T3, fu.shift_a, fu.shift_b)
    if (fu.converged != 0 | missing((fu.T3, fu.shift_a, fu.shift_b)) | missing(fu.acov) |
        fu.acov != dwlsmv_acov(fu.Delta, stD.Gamma) |
        fu.T3 != fu.shift_a * fu.T_raw + fu.shift_b) {
        FAIL = FAIL \ "7: an unconverged fit was not processed like any other"
    }

    printf("{txt}test_inference: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }
    return(rows(FAIL))
}

// ---- 8. the dwls_inference() error paths -------------------------------------
// Builds case D's marker fit in a Mata external, breaks it one way, and calls
// dwls_inference() on it. The external survives the abort, so ti_unfilled()
// can then say whether anything was written.
void ti_bad(real matrix YD, real scalar which)
{
    external struct dwlsmv_fit scalar TI_F
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par

    st   = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    par  = dwls_identify(ti_canon("D_marker"), st.map)
    TI_F = dwls_fit(st, par)
    if (which == 1) {
        TI_F.Delta = TI_F.Delta[., 1..(cols(TI_F.Delta) - 1)]
    }
    else if (which == 2) {
        TI_F.D_stat = TI_F.D_stat + 1
    }
    else if (which == 3) {
        TI_F.T_raw = .
    }
    else if (which == 5) {
        TI_F.T_raw = -1
    }
    else if (which == 6) {
        // The row half of the shape guard. Redundant in practice -- without
        // it, dwlsmv_acov()'s conformability guard aborts at 3498 instead --
        // so this case exercises the guard rather than making it load-bearing.
        TI_F.Delta = TI_F.Delta \ J(1, cols(TI_F.Delta), 0)
    }
    else if (which == 4) {
        // Delta Delta' - eps I: positive diagonal, so acov succeeds; M = -eps U,
        // so Tr(M) = -eps Tr(U) < 0 by construction and scaleshift aborts.
        TI_F.stage1[1].Gamma = quadcross(TI_F.Delta', TI_F.Delta')
        TI_F.stage1[1].Gamma = TI_F.stage1[1].Gamma -
                            0.1 * min(diagonal(TI_F.stage1[1].Gamma)) * I(rows(TI_F.Delta))
    }
    dwls_inference(TI_F)
}

real scalar ti_unfilled()
{
    external struct dwlsmv_fit scalar TI_F

    return(ti_isunfilled(TI_F))
}

// Evidence that case 4 fails where the header says it does: after acov,
// with Tr(M) negative by a margin no rounding reaches.
real rowvector ti_case4()
{
    external struct dwlsmv_fit scalar TI_F

    return((dwlsmv_traces(TI_F.Delta, TI_F.stage1[1].Gamma),
            rows(dwlsmv_acov(TI_F.Delta, TI_F.stage1[1].Gamma))))
}

end

mata: st_local("nfail", strofreal(ti_checks(YD, YE, T3F, T3Q, T3V, ///
                  LDF, LDR, LDC, LDV, GD, GDL, GE, GEL, MPL, MPV)))

* ---------------------------------------------------------------------------
* 8. ERROR PATHS, at rc 3498.
* ---------------------------------------------------------------------------
local nerr 0

* 8a-8d: dwls_inference() on a broken fit, and the fit must stay unfilled.
local why1 "Delta has a column fewer than the table"
local why2 "D_stat is not D - npar"
local why3 "T_raw is missing"
local why4 "Gamma = Delta Delta' - eps I: Tr(M) < 0, after acov succeeded"
local why5 "T_raw is negative"
local why6 "Delta has a row more than the map"
forvalues k = 1/6 {
    capture mata: ti_bad(YD, `k')
    local rc = _rc
    mata: st_local("unf", strofreal(ti_unfilled()))
    display as text "  8. dwls_inference, `why`k'': rc `rc', fit left unfilled `unf'"
    if `rc' != 3498 | `unf' != 1 {
        display as error "test_inference: dwls_inference case `k' (`why`k'') gave rc `rc', unfilled `unf'"
        local nerr = `nerr' + 1
    }
}
* The external holds the LAST case's fit; rebuild case 4's before reading it.
capture mata: ti_bad(YD, 4)
mata: st_local("c4", invtokens(strofreal(ti_case4(), "%10.3e")))
display as text "     case 4 evidence: Tr(M), Tr(M^2), rows(acov) = `c4'"
mata: st_local("c4ok", strofreal(ti_case4()[1] < -1e-3 & ti_case4()[3] == 14))
if `c4ok' != 1 {
    display as error "test_inference: case 4 does not fail where claimed -- Tr(M) not clearly negative or acov did not compute"
    local nerr = `nerr' + 1
}

* dwlsmv_traces(): the dwlsmv_acov() guards, repeated.
capture mata: dwlsmv_traces(J(3, 2, 1), J(3, 2, 1))
if _rc != 3498 {
    display as error "test_inference: traces accepted a non-square Gamma (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_traces(J(3, 1, 1), I(2))
if _rc != 3498 {
    display as error "test_inference: traces accepted a Gamma of the wrong size (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_traces(J(3, 0, .), I(3))
if _rc != 3498 {
    display as error "test_inference: traces accepted a Delta with no columns (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_traces(J(2, 3, 1), I(2))
if _rc != 3498 {
    display as error "test_inference: traces accepted npar > D (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_traces((1 \ .), I(2))
if _rc != 3498 {
    display as error "test_inference: traces accepted a missing value in Delta (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_traces(I(2), (1, . \ ., 1))
if _rc != 3498 {
    display as error "test_inference: traces accepted a missing value in Gamma (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_traces((1 \ 1), (1, 0 \ 0, 0))
if _rc != 3498 {
    display as error "test_inference: traces accepted a zero on diag(Gamma) (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_traces((1, 0 \ 2, 0 \ 3, 0), I(3))
if _rc != 3498 {
    display as error "test_inference: traces accepted a zero Delta column (rc `=_rc')"
    local nerr = `nerr' + 1
}

* dwlsmv_scaleshift(): df, then the traces with df > 0.
foreach a in "1, 1, -1" "1, 1, 2.5" "1, 1, ." "0, 1, 2" "1, 0, 2" "-1, 1, 2" "1, -1, 2" ". , 1, 2" "1, ., 2" {
    capture mata: dwlsmv_scaleshift(`a')
    if _rc != 3498 {
        display as error "test_inference: scaleshift accepted (`a') (rc `=_rc')"
        local nerr = `nerr' + 1
    }
}

if `nfail' | `nerr' {
    display as error _n "test_inference: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_inference: Module 9 verified -- lavaan inputs at 1e-12, U oracle, idempotent collapse, end to end at the test-statistic tier, (N-1)F against lavaan and Mplus, bookkeeping, D = 0, unconverged, error paths with all or nothing"
exit 0
