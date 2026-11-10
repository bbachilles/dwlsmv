*! Module 8: dwlsmv_acov() and dwlsmv_ev() -- the robust sandwich and e(V).
*!
*! Executed under StataMP 18 on 2026-09-23; passed on the first run, so
*! nothing here was adjusted to fit a result. The one tolerance set from a
*! measurement is check 6's, and the measurement was taken in R before the
*! test existed (see check 6 below).
*!
*! Spec Section 4.10; derivation @sec-acov. ACOV is the asymptotic covariance
*! of sqrt(N)(theta-hat - theta_0) from Gamma as stored; e(V) = ACOV/(N-1).
*!
*! The checks, and why each is separate:
*!
*!   1  the block identity ACOV_tautau = Gamma_tautau, on all five of OUR
*!      fits, per entry at a floating-point budget of 1e-14 relative (~45
*!      ulps). Exact in production, not approximate: Delta's threshold block
*!      is a literal identity and its off-diagonal blocks literal zeros
*!      (test_delta.do check 1), so only the sandwich's own rounding remains.
*!      1b: doubling Gamma's threshold rows and columns must leave ACOV's
*!      structural block BIT-IDENTICAL -- every product that could carry the
*!      perturbation is an exact zero -- while the control, ACOV_tautau,
*!      moves by the factor of 4 the perturbation put there.
*!   2  the formula with LAVAAN'S OWN INPUTS -- its analytic Delta at its own
*!      theta-hat and its Gamma, lavaan_acov_delta.csv and lavaan_case_*_
*!      gamma.csv -- against its vcov (lavaan_acov_vcov.csv) and the committed
*!      SEs, at 1e-12. No theta-hat, rho-hat or finite-difference tier is
*!      inherited, so this is the test of the arithmetic, the weight, the
*!      bread and the divisor together; test_gamma.do's 5a, one module up.
*!      2b: /N instead of /(N-1) must FAIL it, and its SEs must miss by
*!      exactly sqrt((N-1)/N) - 1, to 1e-12 -- the whole miss is the divisor.
*!      2b is a POWER check on the fixture and does not call dwlsmv_ev(): it
*!      shows the comparison can tell the two divisors apart. What rejects /N
*!      in the code is checks 2, 4 and 6 (measured by mutation, 2026-09-24).
*!   3  end to end: our theta-hat, our Delta, our Gamma against lavaan's SEs
*!      at spec Section 7.2's SE tier, 1e-5 relative + 1e-8 absolute, every
*!      free parameter including thresholds, all five fits; the worst actual
*!      gap is printed so its distance from the tier is visible.
*!      TIGHTENED 2026-09-24 from 1e-4 abs + 1e-3 rel, whose absolute term
*!      made it ~0.5% relative at these SE sizes: on case D it could not see
*!      a unit weight (0.395 of that tier) or observed information (0.64),
*!      nor the 1.0e-4 /N divisor anywhere. The new tier is set from the
*!      gap's source, lavaan's own convergence (worst 1.3e-7, str_marker),
*!      with ~75x headroom, and sees all three: the /N divisor by 9.2-9.7x
*!      (first written "10x or more", corrected by the Module 9 review), a
*!      unit weight and observed information on case D by far more (measured
*!      in R by read-only agents, one-off; E-032, 2026-09-24). The absolute
*!      term never binds here -- the smallest SE is 0.012, so 1e-5 relative
*!      is 1.2e-7 -- and exists so an SE near zero on a future fixture does
*!      not make the criterion meaningless. Check 2 still pins weight,
*!      bread and divisor at 1e-12; this check adds integration.
*!      3b: the off-diagonal entries, which check 3's SEs cannot see: every
*!      covariance whose lavaan correlation exceeds 1e-3 in magnitude must
*!      carry lavaan's sign, so a Delta column or theta ordering with a
*!      flipped sign fails. A sign test needs no new tolerance; the worst
*!      correlation gap is printed, not asserted.
*!   4  the divisor in isolation: dwlsmv_ev(A, N) is exactly A/(N-1), and is
*!      not A/N -- spec Section 6's stray-N risk, unit-tested as its own step.
*!   5  structure on our fits: ACOV exactly symmetric, positive definite.
*!   6  tier 3, Mplus case D. Mplus divides by N, not N - 1: its threshold
*!      SEs are ours times sqrt((N-1)/N). Pinned on THRESHOLDS because by
*!      check 1 their SEs depend on Gamma_tautau alone, so the convention is
*!      not confounded by Mplus's theta-hat, which differs from lavaan's by up
*!      to ~1e-5. Tolerance 2e-6 on (ratio - 1), 2% of the 1.0e-4 convention
*!      gap; measured worst residual 8.0e-7 (y4|t2) BEFORE the tolerance was
*!      set. The 8-significant-digit rounding floor on SEs near 0.018 is
*!      2.3-2.9e-8 (corrected 2026-09-24; this header first said 2.9e-7):
*!      nine thresholds sit at <= 6.2e-8, and y4|t2 at ~30x the floor is an
*!      UNEXPLAINED outlier. So "not confounded by Mplus's theta-hat" is not
*!      fully true for Mplus, whose own threshold estimates differ from the
*!      sample thresholds in the last digits. None of that approaches the
*!      1.0e-4 gap, and N - 1 is rejected by every threshold. What this pins
*!      is the product of Mplus's Gamma normalization and its divisor, the
*!      same caveat spec 4.10 states for lavaan. The structural SEs are
*!      checked at Section 7.2's TIER-3 SE row (1e-4 abs + 1e-3 rel -- the
*!      old tier-2 row, kept for Mplus on 2026-09-24), since their ~1e-5
*!      scatter is Mplus's theta-hat, not the divisor, and would sit at the
*!      edge of the 1e-5 tier-2 row for a reason that is not ours.
*!      Mplus's SEs come from the SAVEDATA RESULTS block, whose layout is
*!      positional and undocumented (30_mplus_reference.R): estimates at
*!      positions 1-14, SEs at 15-28, thresholds first. The layout is
*!      asserted against the printed 3-decimal table before it is used.
*!   7  the error paths, at rc 3498 and not merely nonzero (a project rule). The
*!      rank-deficient case uses a ZERO Delta column, not two identical
*!      ones: cholinv() returns missing only when a floating-point pivot is
*!      <= 0, and identical columns give a bread [s s; s s] whose second
*!      pivot, s - (s/sqrt(s))^2, happens to round to 0 at s = 14 (the value
*!      this test first used) but is +4.4e-16 at s = 2 -- a finite
*!      inverse near 2.3e15 and no abort (measured 2026-09-24; the review's
*!      mutation 13). A zero column gives a zero pivot under any rounding.
*!      The missing-Gamma and wrong-size-Gamma guards each have their own
*!      case, since mutation showed removing either went unnoticed.
*!
*! Reference CSVs are read with stringcols(_all), strtoreal() and
*! delimiters(","), for the two import traps (float typing, '|' as delimiter); every label here holds '|'
*! or '~'.
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
             src/mata/dwlsmv_admissadvice.mata      ///
             src/mata/dwlsmv_asym.mata              ///
             src/mata/dwlsmv_acov.mata              ///
             src/mata/dwlsmv_ev.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_acov: `f' failed to load"
        exit 9
    }
}

* ---- tier 2: lavaan's sandwich inputs and output, long format ---------------
import delimited using "fixtures/expected/lavaan_acov_delta.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: LDF = st_sdata(., "fit")
mata: LDR = st_sdata(., "row")
mata: LDC = st_sdata(., "col")
mata: LDV = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/lavaan_acov_vcov.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: LVF = st_sdata(., "fit")
mata: LVR = st_sdata(., "row")
mata: LVC = st_sdata(., "col")
mata: LVV = strtoreal(st_sdata(., "value"))

foreach c in D E {
    import delimited using "fixtures/expected/lavaan_case_`c'_gamma.csv", ///
        clear varnames(1) delimiters(",") stringcols(_all)
    qui ds label, not
    local gcols `r(varlist)'
    mata: G`c'  = strtoreal(st_sdata(., "`gcols'"))
    mata: G`c'L = st_sdata(., "label")
}

* ---- tier 2: lavaan's SEs --------------------------------------------------
import delimited using "fixtures/expected/lavaan_case_D_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: DKEY = st_sdata(., "label")
mata: DSE  = strtoreal(st_sdata(., "se"))

import delimited using "fixtures/expected/lavaan_case_E_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: PFIT = st_sdata(., "fit")
mata: PKEY = st_sdata(., "label")
mata: PSE  = strtoreal(st_sdata(., "se"))

* ---- tier 3: Mplus case D --------------------------------------------------
import delimited using "fixtures/expected/mplus_case_D_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: MHDR = st_sdata(., "paramheader")
mata: MPAR = st_sdata(., "param")
mata: MEST = strtoreal(st_sdata(., "est"))
mata: MSE  = strtoreal(st_sdata(., "se"))

import delimited using "fixtures/expected/mplus_case_D_results_raw.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: MRAW = strtoreal(st_sdata(., "value"))

* ---- the data ----------------------------------------------------------------
import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

// Long-format lookup; missing when absent, which every caller reports.
real scalar ta_get(string colvector fitcol, string colvector keycol,
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

real scalar ta_pos(string colvector keys, string scalar k)
{
    real scalar i

    for (i = 1; i <= rows(keys); i++) {
        if (keys[i] == k) {
            return(i)
        }
    }
    return(0)
}

string scalar ta_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The five fits of fixtures/generate/20_lavaan_reference.R, as test_fit.do
// writes them.
string scalar ta_canon(string scalar fit)
{
    string rowvector head, lam6

    head = ("factor f1 f2", "ind y1 y2 y3 y4 y5 y6")
    lam6 = ("lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
            "lambda f2 y4 free", "lambda f2 y5 free", "lambda f2 y6 free")

    if (fit == "D_marker") {
        return(ta_spec(("factor f", "ind y1 y2 y3 y4",
                        "lambda f y2 free", "lambda f y3 free",
                        "lambda f y4 free", "psi f f free",
                        "metric marker f y1")))
    }
    if (fit == "cor_unitvar") {
        return(ta_spec((head, lam6, "psi f1 f2 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    if (fit == "str_unitvar") {
        return(ta_spec((head, lam6, "beta f2 f1 free",
                        "metric unitvariance f1", "metric unitvariance f2")))
    }
    if (fit == "cor_marker") {
        return(ta_spec((head,
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "psi f1 f1 free", "psi f2 f2 free", "psi f1 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    if (fit == "str_marker") {
        return(ta_spec((head,
                        "lambda f1 y2 free", "lambda f1 y3 free",
                        "lambda f2 y5 free", "lambda f2 y6 free",
                        "beta f2 f1 free",
                        "psi f1 f1 free", "psi f2 f2 free",
                        "metric marker f1 y1", "metric marker f2 y4")))
    }
    _error("ta_canon: unknown fit")
}

// lavaan's label for one row of our table, as test_fit.do's tf_lab().
string scalar ta_lab(struct dwlsmv_parmap scalar par, real scalar i,
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

// lavaan's labels for our theta, in theta order.
string colvector ta_thlabels(struct dwlsmv_parmap scalar par, string rowvector fac)
{
    string colvector L
    real scalar i

    L = J(par.npar, 1, "")
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            L[par.fidx[i]] = ta_lab(par, i, fac)
        }
    }
    return(L)
}

real scalar ta_checks(real matrix YD, real matrix YE,
                      string colvector LDF, string colvector LDR,
                      string colvector LDC, real colvector LDV,
                      string colvector LVF, string colvector LVR,
                      string colvector LVC, real colvector LVV,
                      real matrix GD, string colvector GDL,
                      real matrix GE, string colvector GEL,
                      string colvector DKEY, real colvector DSE,
                      string colvector PFIT, string colvector PKEY,
                      real colvector PSE,
                      string colvector MHDR, string colvector MPAR,
                      real colvector MEST, real colvector MSE,
                      real colvector MRAW)
{
    struct dwlsmv_stage1 scalar stD, stE, st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    string colvector FAIL, plab, thl, glab, mlab
    string rowvector fits, fac
    string scalar fitn, s
    real colvector se, sef, seN, ratio_m, tt
    real matrix A, A2, G2, Dl, Gl, Vl, V, Glav
    real scalar d, i, j, k, N, T, r, rN, worst, wrel, g1, g2, m, e_, s_
    real scalar nsig, nbad, cl, co
    real scalar tgt, m_ok

    FAIL = J(0, 1, "")
    stD  = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    stE  = dwls_sample_stats(YE, J(6, 1, 3))
    fits = ("D_marker", "cor_unitvar", "str_unitvar", "cor_marker", "str_marker")

    for (d = 1; d <= cols(fits); d++) {
        fitn = fits[d]
        if (fitn == "D_marker") {
            st   = stD
            fac  = ("f")
            Glav = GD
            glab = GDL
        }
        else {
            st   = stE
            fac  = ("f1", "f2")
            Glav = GE
            glab = GEL
        }
        N   = st.N
        T   = st.map.T
        par = dwls_identify(ta_canon(fitn), st.map)
        fit = dwls_fit(st, par)
        if (fit.converged != 1) {
            FAIL = FAIL \ (fitn + ": dwls_fit() did not converge; nothing below is meaningful")
            continue
        }
        A   = dwlsmv_acov(fit.Delta, st.Gamma)
        tt  = (1::T)
        thl = ta_thlabels(par, fac)

        // ---- 1. ACOV_tautau = Gamma_tautau, per entry ----------------------
        r = dwlsmv_gapratio(A[tt, tt], st.Gamma[tt, tt],
                            1e-15 * max(abs(st.Gamma[tt, tt])), 1e-14)
        printf("{txt}  1. %-11s ACOV_tautau = Gamma_tautau : ratio %8.2e of 1e-14 rel over %f entries\n",
               fitn, r, T * T)
        if (r > 1 | missing(r)) {
            FAIL = FAIL \ (fitn + ": ACOV's threshold block is not Gamma_tautau, ratio " + strofreal(r))
        }

        // ---- 1b. the structural block ignores Gamma's threshold rows -------
        if (fitn == "str_unitvar") {
            G2        = st.Gamma
            G2[tt, .] = 2 :* G2[tt, .]
            G2[., tt] = 2 :* G2[., tt]
            A2 = dwlsmv_acov(fit.Delta, G2)
            g1 = max(abs(A2[(T + 1)::par.npar, (T + 1)::par.npar] -
                         A[(T + 1)::par.npar, (T + 1)::par.npar]))
            g2 = max(abs(A2[tt, tt] :/ A[tt, tt] :- 4))
            printf("{txt}  1b. %-10s structural block moved by %g (must be 0); threshold block ratio off 4 by %8.2e (control)\n",
                   fitn, g1, g2)
            if (g1 != 0) {
                FAIL = FAIL \ ("1b: ACOV's structural block moved when only Gamma's threshold rows changed")
            }
            if (g2 > 1e-13) {
                FAIL = FAIL \ ("1b control: the perturbation did not reach ACOV_tautau as 4x")
            }
        }

        // ---- 2. the formula on lavaan's own Delta and Gamma ----------------
        // lavaan's parameter order: the row labels of its first vcov column.
        plab = J(0, 1, "")
        s    = ""
        for (i = 1; i <= rows(LVF); i++) {
            if (LVF[i] != fitn) continue
            if (s == "") s = LVC[i]
            if (LVC[i] == s) plab = plab \ LVR[i]
        }
        Dl = J(st.map.D, rows(plab), .)
        Vl = J(rows(plab), rows(plab), .)
        for (i = 1; i <= rows(LDF); i++) {
            if (LDF[i] != fitn) continue
            j = ta_pos(st.map.label, LDR[i])
            k = ta_pos(plab, LDC[i])
            if (j & k) Dl[j, k] = LDV[i]
        }
        for (i = 1; i <= rows(LVF); i++) {
            if (LVF[i] != fitn) continue
            j = ta_pos(plab, LVR[i])
            k = ta_pos(plab, LVC[i])
            if (j & k) Vl[j, k] = LVV[i]
        }
        Gl = J(st.map.D, st.map.D, .)
        for (i = 1; i <= st.map.D; i++) {
            j = ta_pos(glab, st.map.label[i])
            if (!j) continue
            for (k = 1; k <= st.map.D; k++) {
                m = ta_pos(glab, st.map.label[k])
                if (m) Gl[i, k] = Glav[j, m]
            }
        }
        if (missing(Dl) | missing(Vl) | missing(Gl) | rows(plab) != par.npar) {
            FAIL = FAIL \ (fitn + ": lavaan's Delta, vcov or Gamma did not fill from the fixtures")
            continue
        }
        sef = J(par.npar, 1, .)
        for (i = 1; i <= par.npar; i++) {
            sef[i] = (fitn == "D_marker" ? ta_get(DKEY, DKEY, DSE, "", plab[i]) :
                                          ta_get(PFIT, PKEY, PSE, fitn, plab[i]))
        }
        // Mata's max() skips missing values, so a label that missed here would
        // silently drop its parameter from checks 2, 2b and 3 -- mutation 12
        // of the 2026-09-24 review went unnoticed until this guard existed.
        if (missing(sef)) {
            FAIL = FAIL \ (fitn + ": a fixture SE lookup missed; checks 2, 2b and 3 would skip it")
            continue
        }

        V   = dwlsmv_ev(dwlsmv_acov(Dl, Gl), N)
        r   = dwlsmv_gapratio(V, Vl, 1e-12 * max(abs(Vl)), 1e-12)
        se  = sqrt(diagonal(V))
        rN  = dwlsmv_gapratio(se, sef, 0, 1e-12)
        printf("{txt}  2. %-11s lavaan inputs -> vcov  : ratio %8.2e ; SEs vs fixture %8.2e (1e-12)\n",
               fitn, r, rN)
        if (r > 1 | rN > 1 | missing(r) | missing(rN)) {
            FAIL = FAIL \ (fitn + ": from lavaan's own Delta and Gamma, vcov ratio " + strofreal(r) +
                           ", SE ratio " + strofreal(rN) + " -- arithmetic, weight, bread or divisor")
        }
        // 2b. /N must fail, and by exactly the divisor.
        A2  = dwlsmv_acov(Dl, Gl) :/ N
        r   = dwlsmv_gapratio(A2, Vl, 1e-12 * max(abs(Vl)), 1e-12)
        seN = sqrt(diagonal(A2))
        g1  = max(abs(seN :/ sef :- sqrt((N - 1) / N)))
        printf("{txt}  2b. %-10s ACOV/N instead          : ratio %8.2e (must be > 1) ; SE/fixture off sqrt((N-1)/N) by %8.2e\n",
               fitn, r, g1)
        if (r <= 1) {
            FAIL = FAIL \ (fitn + ": ACOV/N also matches lavaan -- the divisor check has no power")
        }
        if (g1 > 1e-12) {
            FAIL = FAIL \ (fitn + ": the /N miss is not exactly the divisor ratio")
        }

        // ---- 3. end to end, at the SE tier --------------------------------
        V     = dwlsmv_ev(A, N)
        se    = sqrt(diagonal(V))
        worst = 0
        wrel  = 0
        for (i = 1; i <= par.npar; i++) {
            k = ta_pos(plab, thl[i])
            if (!k) {
                FAIL = FAIL \ (fitn + ": lavaan has no parameter labelled " + thl[i])
                continue
            }
            worst = max((worst, dwlsmv_gapratio(se[i], sef[k], 1e-8, 1e-5)))
            wrel  = max((wrel, abs(se[i] / sef[k] - 1)))
        }
        printf("{txt}  3. %-11s SEs end to end         : ratio %8.2e of the SE tier ; worst relative gap %8.2e over %f\n",
               fitn, worst, wrel, par.npar)
        if (worst > 1) {
            FAIL = FAIL \ (fitn + ": end-to-end SEs disagree with lavaan, ratio " + strofreal(worst))
        }

        // ---- 3b. end-to-end covariance signs --------------------------------
        nsig = 0
        nbad = 0
        wrel = 0
        for (i = 1; i <= par.npar; i++) {
            for (j = 1; j < i; j++) {
                k = ta_pos(plab, thl[i])
                m = ta_pos(plab, thl[j])
                if (!k | !m) continue
                cl = Vl[k, m] / sqrt(Vl[k, k] * Vl[m, m])
                co = V[i, j] / sqrt(V[i, i] * V[j, j])
                wrel = max((wrel, abs(co - cl)))
                if (abs(cl) > 1e-3) {
                    nsig = nsig + 1
                    if (sign(co) != sign(cl)) nbad = nbad + 1
                }
            }
        }
        printf("{txt}  3b. %-10s covariance signs       : %f of %f non-negligible entries disagree (must be 0) ; worst correlation gap %8.2e\n",
               fitn, nbad, nsig, wrel)
        if (nbad | !nsig) {
            FAIL = FAIL \ (fitn + ": an end-to-end covariance has the wrong sign, or none was compared")
        }

        // ---- 5. structure --------------------------------------------------
        g1 = dwlsmv_symmetry(A)
        g2 = dwlsmv_mineig(A)
        printf("{txt}  5. %-11s ACOV asymmetry %g (must be 0), min eigenvalue %g\n", fitn, g1, g2)
        if (g1 != 0 | !(g2 > 0)) {
            FAIL = FAIL \ (fitn + ": ACOV is not exactly symmetric and positive definite")
        }

        // ---- 6. tier 3, Mplus, case D --------------------------------------
        if (fitn == "D_marker") {
            // our theta order is thresholds (the RESULTS block's order), then
            // f=~y2, f=~y3, f=~y4, f~~f -- asserted, since the layout is not
            // documented.
            mlab = thl
            if (mlab[T + 1] != "f=~y2" | mlab[T + 2] != "f=~y3" | mlab[T + 3] != "f=~y4" |
                mlab[T + 4] != "f~~f" | par.npar != 14 | rows(MRAW) < 28) {
                FAIL = FAIL \ ("Mplus: case D's theta is not the order the RESULTS block is read in")
            }
            else {
                m_ok    = 1
                ratio_m = J(par.npar, 1, .)
                for (i = 1; i <= par.npar; i++) {
                    // the printed table's row for this parameter
                    k = 0
                    for (j = 1; j <= rows(MPAR); j++) {
                        if (i <= T) {
                            s = strlower(subinstr(MPAR[j], char(36), "|t", 1))
                            if (MHDR[j] == "Thresholds" & s == mlab[i]) k = j
                        }
                        else if (i < par.npar) {
                            if (MHDR[j] == "FAC.BY" & "f=~" + strlower(MPAR[j]) == mlab[i]) k = j
                        }
                        else if (MHDR[j] == "Variances" & MPAR[j] == "FAC") {
                            k = j
                        }
                    }
                    e_ = MRAW[i]
                    s_ = MRAW[par.npar + i]
                    if (!k) {
                        m_ok = 0
                    }
                    else if (abs(e_ - MEST[k]) > 0.0005 + 1e-12 |
                             abs(s_ - MSE[k])  > 0.0005 + 1e-12) {
                        m_ok = 0
                    }
                    ratio_m[i] = s_ / se[i] - 1
                }
                if (!m_ok) {
                    FAIL = FAIL \ ("Mplus: the RESULTS block no longer matches the printed table -- layout changed")
                }
                else {
                    tgt   = sqrt((N - 1) / N) - 1
                    g1    = max(abs(ratio_m[tt] :- tgt))
                    g2    = min(abs(ratio_m[tt]))
                    worst = 0
                    for (i = T + 1; i <= par.npar; i++) {
                        worst = max((worst, dwlsmv_gapratio(se[i] * sqrt((N - 1) / N),
                                                            MRAW[par.npar + i], 1e-4, 1e-3)))
                    }
                    printf("{txt}  6. Mplus thresholds: SE ratio - 1 within %8.2e of sqrt((N-1)/N) - 1 = %9.4e (2e-6) ; smallest |ratio - 1| %8.2e (N-1 needs ~0)\n",
                           g1, tgt, g2)
                    printf("{txt}     Mplus structural: ratio - 1 from %9.4e to %9.4e ; SE tier under /N, ratio %8.2e\n",
                           min(ratio_m[(T + 1)::par.npar]), max(ratio_m[(T + 1)::par.npar]), worst)
                    if (g1 > 2e-6) {
                        FAIL = FAIL \ ("Mplus: threshold SEs are not ours times sqrt((N-1)/N)")
                    }
                    if (g2 < 5e-5) {
                        FAIL = FAIL \ ("Mplus: some threshold SE is consistent with N - 1 -- the convention pin has no power")
                    }
                    if (worst > 1) {
                        FAIL = FAIL \ ("Mplus: structural SEs outside the SE tier under /N")
                    }
                }
            }
        }
    }

    // ---- 4. the divisor in isolation ---------------------------------------
    A = (1, 2 \ 3, 4.5)
    g1 = max(abs(dwlsmv_ev(A, 5000) - A :/ 4999))
    g2 = max(abs(dwlsmv_ev(A, 5000) - A :/ 5000))
    printf("{txt}  4. dwlsmv_ev(A, 5000) - A/4999 = %g (must be 0) ; - A/5000 = %g (must not be)\n", g1, g2)
    if (g1 != 0 | g2 == 0) {
        FAIL = FAIL \ ("4: dwlsmv_ev() is not exactly ACOV/(N-1)")
    }
    if (max(abs(dwlsmv_ev(A, 2) - A)) != 0) {
        FAIL = FAIL \ ("4: dwlsmv_ev(A, 2) is not A/1")
    }

    printf("{txt}test_acov: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(ta_checks(YD, YE, LDF, LDR, LDC, LDV,  ///
                  LVF, LVR, LVC, LVV, GD, GDL, GE, GEL, DKEY, DSE,        ///
                  PFIT, PKEY, PSE, MHDR, MPAR, MEST, MSE, MRAW)))

* ---------------------------------------------------------------------------
* 7. ERROR PATHS, at rc 3498. Each would otherwise return a plausible matrix
* or die somewhere less explicable.
* ---------------------------------------------------------------------------
local nerr 0

* Gamma not square / not conformable with Delta.
capture mata: dwlsmv_acov(J(3, 2, 1), J(3, 2, 1))
if _rc != 3498 {
    display as error "test_acov: acov accepted a non-square Gamma (rc `=_rc')"
    local nerr = `nerr' + 1
}
* More parameters than statistics.
capture mata: dwlsmv_acov(J(2, 3, 1), I(2))
if _rc != 3498 {
    display as error "test_acov: acov accepted npar > D (rc `=_rc')"
    local nerr = `nerr' + 1
}
* A square Gamma of the wrong size. The case above also fails squareness, so
* it never reached the row-count half of that guard.
capture mata: dwlsmv_acov(J(3, 1, 1), I(2))
if _rc != 3498 {
    display as error "test_acov: acov accepted a Gamma whose size does not match Delta (rc `=_rc')"
    local nerr = `nerr' + 1
}
* Missing values, in Delta and -- separately -- in Gamma off the diagonal,
* where the zero-diagonal guard cannot catch it.
capture mata: dwlsmv_acov((1 \ .), I(2))
if _rc != 3498 {
    display as error "test_acov: acov accepted a missing value in Delta (rc `=_rc')"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_acov(I(2), (1, . \ ., 1))
if _rc != 3498 {
    display as error "test_acov: acov accepted a missing value in Gamma (rc `=_rc')"
    local nerr = `nerr' + 1
}
* A zero on diag(Gamma): W_D does not exist.
capture mata: dwlsmv_acov((1 \ 1), (1, 0 \ 0, 0))
if _rc != 3498 {
    display as error "test_acov: acov accepted a zero on diag(Gamma) (rc `=_rc')"
    local nerr = `nerr' + 1
}
* A rank-deficient Delta: a ZERO column, so the bread's pivot is 0 under any
* rounding (see the header, check 7). invsym() would return SEs of 0.
capture mata: dwlsmv_acov((1, 0 \ 2, 0 \ 3, 0), I(3))
if _rc != 3498 {
    display as error "test_acov: acov accepted a rank-deficient Delta (rc `=_rc')"
    local nerr = `nerr' + 1
}
* dwlsmv_ev(): non-square, N < 2, non-integer N, missing N.
capture mata: dwlsmv_ev(J(2, 3, 1), 100)
if _rc != 3498 {
    display as error "test_acov: ev accepted a non-square acov (rc `=_rc')"
    local nerr = `nerr' + 1
}
foreach n in 1 2.5 . {
    capture mata: dwlsmv_ev(I(2), `n')
    if _rc != 3498 {
        display as error "test_acov: ev accepted N = `n' (rc `=_rc')"
        local nerr = `nerr' + 1
    }
}

if `nfail' | `nerr' {
    display as error _n "test_acov: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_acov: Module 8 verified -- block identity, lavaan inputs at 1e-12, /N rejected, end to end, divisor, structure, Mplus convention, error paths"
exit 0
