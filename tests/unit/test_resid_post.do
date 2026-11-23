*! Residual diagnostics: what dwlsmv posts (e(srmr), e(Resid), e(Nresid)) and when.
*!
*! Written by the main session 2026-11-22 (spec, posted results: *Residual diagnostics*;
*! derivation @sec-resid-diag SR5; owner decisions O1-O4, O6). The calculation itself is tested
*! against lavaan in test_resid_diag.do; this test is about the glue:
*!   1  G_efa1 (case G, efa(F1)): posted; e(srmr) within the 99t fixture's srmr_tol (case G's
*!      r and gamma_kk meet their tiers, so srmr_tol holds); e(Resid) and e(Nresid) p x p, rows
*!      and columns e(indicators); Resid symmetric with diagonal 0; Nresid diagonal missing and
*!      nothing else missing; e(srmr) equals sqrt(sum lower e^2 / (P + p)) of the posted e(Resid)
*!      at 1e-14 (posting moves nothing).
*!   2  R_cfa1 (emptycell, one factor): Nresid missing at y1-y2 (repaired, O1) and nowhere else
*!      off the diagonal; Resid present there; e(srmr) within srmr_tol (r and rho asserted for
*!      this fit in test_resid_diag.do).
*!   3  E3_df0 (case E y1-y3, df = 0): posted (O4); e(df_T3) = 0.
*!   4  a fit that did not converge (case F, efa(F1-F4), not converged on 2026-11-18 and
*!      2026-11-21): nothing posted, e(srmr) absent, not missing.
*!   5  a conditional fit (one covariate): nothing posted.
*!   6  G = 2 (mg_E y1-y5, group(g), invariance(metric)): nothing posted, fit indices posted.
*!   7  replay: e(srmr), e(Resid) unchanged.
*! Added after the end-of-module review (2026-11-22; drafted by a read-only agent, run by the
*! main session), for mutants that survived:
*!   7b the header's SRMR line: one line, aligned with TLI, value e(srmr) %9.4f; 5b none for a
*!      conditional fit.
*!   8  emptycell with indicators y3 y1 y4 y2: stripes follow e(indicators), not sorted names;
*!      the repaired y1-y2 is missing at [2,4].
*!   9  mg_P one group, efa(F1 F2), converged and inadmissible: posted (O3).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index
local nfail 0

* fixture scalars, read as strings (stringcols(_all), strtoreal)
import delimited using "fixtures/expected/lavaan_srmr.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
foreach f in G_efa1 R_cfa1 {
    foreach q in srmr_hand srmr_tol {
        quietly levelsof value if fit == "`f'" & quantity == "`q'", local(v) clean
        mata: st_local("`f'_`q'", strofreal(strtoreal("`v'"), "%21.0g"))
    }
}

mata:
// checks on the posted matrices; returns a message, "" when all hold
string scalar trp_mats(real scalar expect_missing_pairs)
{
    real matrix E, Z
    string matrix rs, cs
    string rowvector ind
    real scalar p, P, i, j, s, m

    E   = st_matrix("e(Resid)")
    Z   = st_matrix("e(Nresid)")
    ind = tokens(st_global("e(indicators)"))
    p   = cols(ind)
    P   = p * (p - 1) / 2
    if (rows(E) != p | cols(E) != p | rows(Z) != p | cols(Z) != p) return("not p x p")
    rs = st_matrixrowstripe("e(Resid)")
    cs = st_matrixcolstripe("e(Nresid)")
    if (rs[., 2]' != ind | cs[., 2]' != ind) return("stripes are not e(indicators)")
    if (E != E') return("Resid not symmetric")
    if (any(diagonal(E) :!= 0)) return("Resid diagonal not 0")
    if (missing(E)) return("Resid has a missing entry")
    if (missing(diagonal(Z)) != p) return("Nresid diagonal not missing")
    if (Z != Z') return("Nresid not symmetric")
    m = (missing(Z) - p) / 2
    if (m != expect_missing_pairs) return("Nresid missing pairs " + strofreal(m) + ", expected " + strofreal(expect_missing_pairs))
    if (rows(st_numscalar("e(srmr)")) == 0) return("e(srmr) not posted")
    s = 0
    for (j = 1; j < p; j++) {
        for (i = j + 1; i <= p; i++) {
            s = s + E[i, j]^2
        }
    }
    if (abs(sqrt(s / (P + p)) - st_numscalar("e(srmr)")) > 1e-14) return("e(srmr) is not the SRMR of e(Resid)")
    return("")
}

real scalar trp_has(string scalar nm)
{
    return(st_matrix(nm) != J(0, 0, .))
}

// the SRMR line of a logged header. want = 1: exactly one line holding "SRMR", starting in the
// same column as the TLI line and with its "=" in the same column, and the text after "=" is
// e(srmr) as %9.4f; want = 0: no such line. Returns a message, "" when it holds.
string scalar trp_hdr(string scalar fn, real scalar want)
{
    string colvector L
    string scalar val, wantv
    real scalar k, ns, ct, et, cs, es

    L   = cat(fn)
    ns  = 0
    ct  = .
    et  = .
    cs  = .
    es  = .
    val = ""
    for (k = 1; k <= rows(L); k++) {
        if (strpos(L[k], "TLI (scaled)") > 0) {
            ct = strpos(L[k], "TLI")
            et = strpos(L[k], "=")
        }
        if (strpos(L[k], "SRMR") > 0) {
            ns  = ns + 1
            cs  = strpos(L[k], "SRMR")
            es  = strpos(L[k], "=")
            val = strtrim(substr(L[k], es + 1, .))
        }
    }
    if (want == 0) {
        if (ns != 0) return("a header line holding SRMR where none is posted")
        return("")
    }
    if (ns != 1) return("SRMR lines in the header: " + strofreal(ns) + ", expected 1")
    if (missing(ct)) return("no TLI line to align against")
    if (cs != ct | es != et) return("SRMR line columns (" + strofreal(cs) + ", " + strofreal(es) + ") differ from the TLI line's (" + strofreal(ct) + ", " + strofreal(et) + ")")
    wantv = strtrim(strofreal(st_numscalar("e(srmr)"), "%9.4f"))
    if (val != wantv) return("header SRMR [" + val + "] is not e(srmr) as %9.4f [" + wantv + "]")
    return("")
}
end

* ---- 1. G_efa1 ----------------------------------------------------------------------------
import delimited using "fixtures/data/case_G.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F1 -> y1-y6), efa(F1)
mata: st_local("msg", trp_mats(0))
local gap = abs(e(srmr) - `G_efa1_srmr_hand')
display as text " 1. G_efa1: e(srmr) " %12.10f e(srmr) ", fixture gap " %9.2e `gap' " (srmr_tol " %9.2e `G_efa1_srmr_tol' ") `msg'"
if "`msg'" != "" | !(`gap' <= `G_efa1_srmr_tol') {
    display as error " 1. failed"
    local ++nfail
}
* ---- 7. replay ----------------------------------------------------------------------------
tempname R0 s0
matrix `R0' = e(Resid)
scalar `s0' = e(srmr)
quietly dwlsmv
if e(srmr) != scalar(`s0') | mreldif(`R0', e(Resid)) != 0 {
    display as error " 7. replay changed e(srmr) or e(Resid)"
    local ++nfail
}
else display as text " 7. replay: unchanged"

* ---- 7b. the header's SRMR line (review mutant, 2026-11-22) ---------------------------------
tempfile hdrlg
quietly log using "`hdrlg'", text replace name(hdrlog)
dwlsmv
quietly log close hdrlog
mata: st_local("msg", trp_hdr("`hdrlg'", 1))
if "`msg'" != "" {
    display as error " 7b. header SRMR line: `msg'"
    local ++nfail
}
else display as text " 7b. header: the SRMR line sits under TLI and shows e(srmr)"

* ---- 2. R_cfa1 ----------------------------------------------------------------------------
import delimited using "fixtures/data/emptycell.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4)
mata: st_local("msg", trp_mats(1))
mata: st_local("rep", strofreal(missing(st_matrix("e(Nresid)")[2, 1]) & !missing(st_matrix("e(Resid)")[2, 1]) & st_matrix("e(Resid)")[2, 1] != 0))
local gap = abs(e(srmr) - `R_cfa1_srmr_hand')
display as text " 2. R_cfa1: e(srmr) " %12.10f e(srmr) ", fixture gap " %9.2e `gap' ", y1-y2 Nresid missing and Resid present: `rep' `msg'"
if "`msg'" != "" | `rep' != 1 | !(`gap' <= `R_cfa1_srmr_tol') {
    display as error " 2. failed"
    local ++nfail
}

* ---- 3. E3_df0 ----------------------------------------------------------------------------
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3)
mata: st_local("msg", trp_mats(0))
display as text " 3. df = 0: e(df_T3) " e(df_T3) ", e(srmr) " %9.2e e(srmr) " `msg'"
if "`msg'" != "" | e(df_T3) != 0 | "`e(srmr)'" == "" {
    display as error " 3. failed"
    local ++nfail
}

* ---- 4. not converged ---------------------------------------------------------------------
import delimited using "fixtures/data/case_F.csv", clear varnames(1) delimiters(",")
capture quietly dwlsmv (F1 F2 F3 F4 -> y1-y9), efa(F1 F2 F3 F4)
local rc4 = _rc
mata: st_local("hasm", strofreal(trp_has("e(Resid)") | trp_has("e(Nresid)")))
display as text " 4. not converged: rc `rc4', e(converged) " e(converged) ", e(srmr) [`e(srmr)'], matrices `hasm'"
if `rc4' | e(converged) != 0 | "`e(srmr)'" != "" | `hasm' {
    display as error " 4. failed (expected rc 0, a fit that did not converge, nothing posted)"
    local ++nfail
}

* ---- 5. conditional ------------------------------------------------------------------------
clear
set seed 20261114
quietly set obs 1000
quietly generate byte x = runiform() < .5
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal()
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
quietly dwlsmv (F -> y1 y2 y3 y4) (F <- x)
mata: st_local("hasm", strofreal(trp_has("e(Resid)") | trp_has("e(Nresid)")))
display as text " 5. conditional: e(covariates) `e(covariates)', e(srmr) [`e(srmr)'], matrices `hasm'"
if "`e(covariates)'" != "x" | "`e(srmr)'" != "" | `hasm' {
    display as error " 5. failed"
    local ++nfail
}
quietly log using "`hdrlg'", text replace name(hdrlog)
dwlsmv
quietly log close hdrlog
mata: st_local("msg", trp_hdr("`hdrlg'", 0))
if "`msg'" != "" {
    display as error " 5b. conditional header: `msg'"
    local ++nfail
}
else display as text " 5b. conditional header: no SRMR line"

* ---- 6. G = 2 -----------------------------------------------------------------------------
import delimited using "fixtures/data/mg_E.csv", clear
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(metric) metric(unitvariance)
mata: st_local("hasm", strofreal(trp_has("e(Resid)") | trp_has("e(Nresid)")))
display as text " 6. G = 2: e(srmr) [`e(srmr)'], matrices `hasm', e(rmsea_s) " e(rmsea_s)
if "`e(srmr)'" != "" | `hasm' | "`e(rmsea_s)'" == "" {
    display as error " 6. failed"
    local ++nfail
}

* ---- 8. indicators out of sorted order (review mutant, 2026-11-22) -------------------------
import delimited using "fixtures/data/emptycell.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y3 y1 y4 y2)
local ord "`e(indicators)'"
mata: st_local("msg", trp_mats(1))
mata: st_local("pos", strofreal(missing(st_matrix("e(Nresid)")[2, 4]) & missing(st_matrix("e(Nresid)")[4, 2]) & !missing(st_matrix("e(Nresid)")[1, 2]) & !missing(st_matrix("e(Nresid)")[3, 4])))
display as text " 8. out-of-order indicators: e(indicators) [`ord'], y1-y2 Nresid missing at [2,4]: `pos' `msg'"
if "`ord'" != "y3 y1 y4 y2" | "`msg'" != "" | `pos' != 1 {
    display as error " 8. failed (premise: e(indicators) is y3 y1 y4 y2; stripes of Resid and Nresid are e(indicators) in order)"
    local ++nfail
}

* ---- 9. a converged inadmissible efa() fit posts (O3; review mutant, 2026-11-22) -----------
* mg_P fitted as one group, efa(F1 F2): converged and inadmissible (found by a search of the
* committed data, 2026-11-22; efatie_T p = 8, m = 3 is the other such fit)
import delimited using "fixtures/data/mg_P.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F1 F2 -> y1-y5), efa(F1 F2)
local pre = (e(converged) == 1 & "`e(rot_nose)'" == "inadmissible")
mata: st_local("msg", trp_mats(0))
display as text " 9. inadmissible: premise `pre', e(srmr) [`e(srmr)'] `msg'"
if !`pre' | "`e(srmr)'" == "" | "`msg'" != "" {
    display as error " 9. failed (premise: converged and inadmissible; then posted)"
    local ++nfail
}

adopath - "`here'/src/ado"
adopath - "`here'/build"
if `nfail' {
    display as error "test_resid_post: `nfail' failure(s)"
    exit 9
}
display as result "test_resid_post: posted where SR5 says, absent elsewhere"
exit 0
