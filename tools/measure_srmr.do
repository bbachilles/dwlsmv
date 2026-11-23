*! measure_srmr.do -- SRMR and the largest normalized residuals of dwlsmv EFA fits, by m.
*!
*! A descriptive measurement for E-117 (owner gate,
*! 2026-11-21), not a test: it prints, it asserts only its own preconditions. Written by the main
*! session 2026-11-21 to replace the agent draft, which recoded the items itself (garbled), assumed
*! e(Lambda)'s row order, and was never run.
*!
*! Definitions measured (package section 2, decision 1; derivation section to come):
*!   residual   e_ij = r_ij - rho_ij(theta-hat), i > j, r from Stage 1's s, rho the off-diagonal of
*!              Lambda Phi Lambda' (rotation-invariant, so e(Lambda) and e(Phi) of the rotated fit)
*!   SRMR       sqrt( sum e_ij^2 / (P + p) ), P = p(p-1)/2   (lavaan's categorical srmr)
*!   CRMR       sqrt( sum e_ij^2 / P )                        (lavaan's crmr_nomean)
*!   normalized z_ij = e_ij / sqrt(Gamma_kk / N), Gamma = N x the ACOV of s (types: Gamma field)
*! Stage 1 is dwlsmv's own (dwlsmv_m11_stage1: D9 recoding, then dwls_sample_stats, keepbound 0)
*! on e(sample), so r is the vector the fit used. A pair held at a bound would be flagged (st.bound).
*!
*! Datasets: build/dwlsmv_example.dta (y1-y6, N 5000, two factors) and fixtures/data/case_F.csv
*! (y1-y9, three factors, N 5000). case_F was chosen after measuring: efacompare's 2026-11-18
*! record has its m = 4 not converged and m = 5 rotation failure, so it reaches the rows the
*! residual table must handle.
*!
*! Run from the repository root:  do tools/measure_srmr.do   (log: measure_srmr.log)
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 120
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

mata:
// One fit in e(): print m, SRMR, CRMR, the largest |raw| residual and the top 3 by |normalized|.
void msr_row(string scalar items, string scalar touse, real scalar m)
{
    struct dwlsmv_stage1 scalar st
    pointer() scalar pst
    string rowvector it
    string matrix rn
    real matrix L, Ph, S, Rm, E, Z, out
    real colvector e, z, ii, jj, bnd, ord
    real scalar p, d, k, P, srmr, crmr, kr

    it = tokens(items)
    p  = cols(it)
    L  = st_matrix("e(Lambda)")
    rn = st_matrixrowstripe("e(Lambda)")
    if (rows(L) != p | cols(L) != m) {
        _error("e(Lambda) is not p x m")
    }
    // precondition: e(Lambda) rows are the items in the order given
    for (k = 1; k <= p; k++) {
        if (rn[k, 2] != it[k]) {
            errprintf("e(Lambda) row %f is %s, item %f is %s\n", k, rn[k, 2], k, it[k])
            _error("e(Lambda) row order differs from the item order")
        }
    }
    if (m == 1) {
        Ph = I(1)
    }
    else {
        Ph = st_matrix("e(Phi)")
    }
    S = L * Ph * L'

    dwlsmv_m11_stage1(touse, items, 0)
    pst = findexternal("dwlsmv_m11_st")
    st  = *pst
    if (st.N != st_numscalar("e(N)")) {
        _error("Stage 1 N differs from e(N)")
    }

    P = p * (p - 1) / 2
    e = J(P, 1, .); z = J(P, 1, .); ii = J(P, 1, .); jj = J(P, 1, .); bnd = J(P, 1, 0)
    kr = 0
    for (d = 1; d <= st.map.D; d++) {
        if (st.map.stat_type[d] == 2) {
            kr++
            ii[kr] = st.map.var_i[d]
            jj[kr] = st.map.var_j[d]
            e[kr]  = st.s[d] - S[ii[kr], jj[kr]]
            z[kr]  = e[kr] / sqrt(st.Gamma[d, d] / st.N)
            bnd[kr] = st.bound[kr]
        }
    }
    if (kr != P) {
        _error("Stage 1 does not have p(p-1)/2 correlations")
    }
    if (missing(e) | missing(z)) {
        printf("  m = %f: %f missing residual(s), %f missing normalized\n", m, missing(e), missing(z))
    }
    srmr = sqrt(sum(e :^ 2) / (P + p))
    crmr = sqrt(sum(e :^ 2) / P)
    printf("  m = %f  SRMR %8.5f  CRMR %8.5f  pairs at a bound %f\n", m, srmr, crmr, sum(bnd))
    k = order(-abs(e), 1)[1]
    printf("        largest |raw|: %s-%s  e %9.5f  z %8.2f\n", it[ii[k]], it[jj[k]], e[k], z[k])
    ord = order((-abs(z), ii, jj), (1, 2, 3))
    for (k = 1; k <= 3; k++) {
        d = ord[k]
        printf("        rank %f by |z|: %s-%s  e %9.5f  z %8.2f\n", k, it[ii[d]], it[jj[d]], e[d], z[d])
    }
    dwlsmv_m11_clear()
}
end

capture program drop msr_run
program define msr_run
    args items mmax
    forvalues m = 1/`mmax' {
        local fl
        forvalues k = 1/`m' {
            local fl `fl' F`k'
        }
        capture quietly dwlsmv (`fl' -> `items'), efa(`fl')
        local rc = _rc
        if `rc' {
            display "  m = `m': rc `rc' (no residuals)"
            continue
        }
        if e(converged) != 1 {
            display "  m = `m': not converged (no residuals)"
            continue
        }
        display "  m = `m': df " e(df_T3) ", T3 " %9.3f e(T3) ", RMSEA " %6.4f e(rmsea_s) ", CFI " %6.4f e(cfi_s)
        tempvar touse
        quietly generate byte `touse' = e(sample)
        mata: msr_row("`items'", "`touse'", `m')
    }
end

display _n "== build/dwlsmv_example.dta, y1-y6 (m = 1..3)"
use build/dwlsmv_example.dta, clear
msr_run "y1 y2 y3 y4 y5 y6" 3

display _n "== fixtures/data/case_F.csv, y1-y9 (m = 1..5)"
import delimited using "fixtures/data/case_F.csv", clear varnames(1) delimiters(",")
msr_run "y1 y2 y3 y4 y5 y6 y7 y8 y9" 5
