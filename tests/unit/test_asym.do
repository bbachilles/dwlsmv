*! Module 10: the Gamma symmetry contract -- dwlsmv_asym() and the guards in
*! dwlsmv_traces() and dwlsmv_acov() that call it.
*!
*! Executed under StataMP 18 on 2026-09-25 and it passes.
*!
*! Derivation @sec-gamma-sym; spec Section 5, Module 10. The measure is
*! a(A) = max |A_ij - A_ji| / sqrt(A_ii A_jj) and the budget 4 D eps. Both are
*! defined ONCE, in dwlsmv_asym(). The guards and (when it lands) the
*! diagnostic carry no tolerance of their own, and check 4 is what makes that
*! a tested property rather than a sentence: the same matrices, just inside and
*! just outside the budget, must get the same verdict from every consumer.
*!
*!   1. OUR Gamma is exactly symmetric: a = 0 bitwise, because quadcross()
*!      mirrors one triangle. The budget is 4 D eps with D = 16.
*!   2. LAVAAN's Gamma is not exactly symmetric, and sits inside the budget.
*!      Both committed fixtures, against the derivation's R figures
*!      (5.761809e-16 and 6.834272e-16) at 1e-6 relative. That is the finest
*!      a 7-significant-digit quotation supports. The two languages also form
*!      sqrt(A_ii A_jj) in different orders, which moves a only in its last
*!      few bits.
*!   3. INVARIANCE to rescaling the statistics, S A S for a positive diagonal
*!      S spanning six decades. Asserted in ABSOLUTE terms at 2 eps, as in the
*!      derivation: the rescaling rounds each entry once, which by itself moves
*!      |A_ij - A_ji| by up to eps sqrt(A_ii A_jj).
*!   4. AGREEMENT. One entry of an exactly symmetric Gamma is moved so that
*!      a = 0.5 and 2 times the budget. At 0.5 dwlsmv_asym() reports inside,
*!      and dwlsmv_traces() and dwlsmv_acov() both run. At 2 it reports
*!      outside, and both abort with rc 3498.
*!   5. ERROR PATHS, each rc 3498: non-square, empty, missing, zero diagonal,
*!      negative diagonal.
*!
*! Reference CSVs are read with stringcols(_all) + strtoreal() and an explicit
*! delimiters(","), for the two import traps (float typing, '|' as delimiter).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
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
             src/mata/dwlsmv_asym.mata              ///
             src/mata/dwlsmv_acov.mata              ///
             src/mata/dwlsmv_traces.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_asym: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/expected/lavaan_case_D_gamma.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
mata: G_LAVD = strtoreal(st_sdata(., "`r(varlist)'"))

import delimited using "fixtures/expected/lavaan_case_E_gamma.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
mata: G_LAVE = strtoreal(st_sdata(., "`r(varlist)'"))

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

// One off-diagonal entry of an exactly symmetric A, moved so that the
// measure is f times the budget. A is symmetrized first as (A + A')/2,
// which is exactly symmetric because floating-point addition commutes.
real matrix dwlsmv_t_bump(real matrix A0, real scalar f)
{
    real matrix A
    real scalar k

    A = (A0 + A0') :/ 2
    k = rows(A)
    A[k, k - 1] = A[k, k - 1] +
                  f * 4 * k * epsilon(1) * sqrt(A[k, k] * A[k - 1, k - 1])
    return(A)
}

real scalar dwlsmv_t_asym(real matrix YD, real matrix G_LAVD, real matrix G_LAVE)
{
    struct dwlsmv_stage1 scalar st
    string colvector FAIL
    real rowvector as, aD, aE, s
    real matrix S, A
    real scalar gapD, gapE, inv, eps, i, f, k, q

    FAIL = J(0, 1, "")
    eps  = epsilon(1)

    // ---- 1. our Gamma --------------------------------------------------------
    st = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    as = dwlsmv_asym(st.Gamma)
    printf("{txt}  1. our Gamma      : a = %g, budget %g (4 D eps, D = %f)\n",
           as[1], as[2], rows(st.Gamma))
    if (as[1] != 0) {
        FAIL = FAIL \ ("1: our Gamma is not bitwise symmetric, a = " + strofreal(as[1]))
    }
    if (as[2] != 4 * 16 * eps) {
        FAIL = FAIL \ ("1: the budget is not 4 D eps: " + strofreal(as[2]))
    }

    // ---- 2. lavaan's Gamma ---------------------------------------------------
    // max() skips missing values, so the fixture's completeness is checked
    // first: a failed strtoreal() would otherwise read as a small a.
    if (missing(G_LAVD) | missing(G_LAVE) | rows(G_LAVD) != 16 | rows(G_LAVE) != 27) {
        FAIL = FAIL \ "2: a lavaan Gamma fixture did not import as a complete matrix"
    }
    else {
        aD   = dwlsmv_asym(G_LAVD)
        aE   = dwlsmv_asym(G_LAVE)
        gapD = abs(aD[1] - 5.761809e-16) / 5.761809e-16
        gapE = abs(aE[1] - 6.834272e-16) / 6.834272e-16
        printf("{txt}  2. lavaan Gamma   : case D a = %g (budget %g), case E a = %g (budget %g)\n",
               aD[1], aD[2], aE[1], aE[2])
        printf("{txt}                      against the derivation's R figures: rel gap %g, %g\n",
               gapD, gapE)
        if (!(aD[1] > 0 & aE[1] > 0)) {
            FAIL = FAIL \ "2: lavaan's Gamma read as exactly symmetric -- the check would be vacuous"
        }
        if (!(aD[1] <= aD[2] & aE[1] <= aE[2])) {
            FAIL = FAIL \ "2: a committed lavaan Gamma is outside the symmetry budget"
        }
        // The R figures are quoted to 7 significant digits, so 1e-6 is the
        // finest this comparison can be; last-bit differences sit far below.
        if (!(gapD <= 1e-6 & gapE <= 1e-6)) {
            FAIL = FAIL \ "2: a disagrees with the derivation's R measurement"
        }
    }

    // ---- 3. invariance to rescaling the statistics ---------------------------
    s   = 10 :^ (-3 :+ 6 :* ((0::26)' :/ 26))
    S   = s' * s
    inv = abs(dwlsmv_asym(G_LAVE :* S)[1] - dwlsmv_asym(G_LAVE)[1])
    printf("{txt}  3. invariance     : |a(SAS) - a(A)| = %g (bound 2 eps = %g)\n", inv, 2 * eps)
    if (!(inv <= 2 * eps)) {
        FAIL = FAIL \ ("3: the measure is not invariant to rescaling: " + strofreal(inv))
    }

    // ---- 4. agreement: the constructed matrices ------------------------------
    // The routines' verdicts are taken at Stata level below, with capture,
    // because an abort cannot be caught inside Mata. Here only the
    // construction itself is checked: a must be 0.5x and 2x the budget.
    // The perturbation, about 1e-15, is added to an entry near 0.1, so the sum
    // is quantized to that entry's ulp. The achievable ratio is therefore f
    // to within ulp(A_kl) / (budget sqrt(A_kk A_ll)), about 0.2% here. That
    // bound is computed from the entry itself rather than assumed.
    for (i = 1; i <= 2; i++) {
        f  = (0.5, 2)[i]
        A  = dwlsmv_t_bump(G_LAVD, f)
        as = dwlsmv_asym(A)
        k  = rows(A)
        q  = epsilon(A[k, k - 1]) / (as[2] * sqrt(A[k, k] * A[k - 1, k - 1]))
        printf("{txt}  4. at %3.1fx budget : a / budget = %g (quantization bound %g)\n",
               f, as[1] / as[2], q)
        if (abs(as[1] / as[2] - f) > q) {
            FAIL = FAIL \ ("4: the constructed asymmetry is not " + strofreal(f) + "x the budget")
        }
    }

    printf("{txt}test_asym: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// A fixed, full-rank Delta: its values do not matter, only that the two
// routines get past every other guard and reach the symmetry one.
real matrix dwlsmv_t_delta()
{
    real matrix Dl
    real scalar i

    Dl = J(16, 3, .)
    for (i = 1; i <= 16; i++) {
        Dl[i, .] = (1, i / 16, (i / 16)^2)
    }
    return(Dl)
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_asym(YD, G_LAVD, G_LAVE)))

* ---------------------------------------------------------------------------
* 4. AGREEMENT -- the same matrices, just inside and just outside the budget,
* get the same verdict from dwlsmv_traces() and dwlsmv_acov(). Plain matrices
* at the Mata prompt, which is safe; the interactive-level struct trap is
* about structs.
* ---------------------------------------------------------------------------
mata: DL_T  = dwlsmv_t_delta()
mata: A_IN  = dwlsmv_t_bump(G_LAVD, 0.5)
mata: A_OUT = dwlsmv_t_bump(G_LAVD, 2)

local nagr 0
capture mata: T_TMP = dwlsmv_traces(DL_T, A_IN)
if _rc != 0 {
    display as error "test_asym: traces refused a Gamma at 0.5x the budget, rc " _rc
    local nagr = `nagr' + 1
}
capture mata: T_TMP = dwlsmv_acov(DL_T, A_IN)
if _rc != 0 {
    display as error "test_asym: acov refused a Gamma at 0.5x the budget, rc " _rc
    local nagr = `nagr' + 1
}
capture mata: T_TMP = dwlsmv_traces(DL_T, A_OUT)
if _rc != 3498 {
    display as error "test_asym: traces on a Gamma at 2x the budget gave rc " _rc ", not 3498"
    local nagr = `nagr' + 1
}
capture mata: T_TMP = dwlsmv_acov(DL_T, A_OUT)
if _rc != 3498 {
    display as error "test_asym: acov on a Gamma at 2x the budget gave rc " _rc ", not 3498"
    local nagr = `nagr' + 1
}
display as text "  4. agreement      : traces and acov, inside and outside the budget -- `nagr' disagreement(s)"
local nfail = `nfail' + `nagr'

* ---------------------------------------------------------------------------
* 5. ERROR PATHS -- rc 3498 each, never merely nonzero.
* ---------------------------------------------------------------------------
local nerr 0

capture mata: dwlsmv_asym(J(2, 3, 1))
if _rc != 3498 {
    display as error "test_asym: non-square A gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_asym(J(0, 0, .))
if _rc != 3498 {
    display as error "test_asym: empty A gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_asym((1, .5 \ ., 1))
if _rc != 3498 {
    display as error "test_asym: A with a missing entry gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_asym((1, 0 \ 0, 0))
if _rc != 3498 {
    display as error "test_asym: A with a zero diagonal entry gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_asym((1, 0 \ 0, -1))
if _rc != 3498 {
    display as error "test_asym: A with a negative diagonal entry gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_asym: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_asym: the symmetry contract holds -- one measure, one budget, every consumer agrees"
exit 0
