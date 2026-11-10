*! dwlsmv -- marker-metric start trap: rate study (E-010, 2026-10-06)
*!
*! Executed 2026-10-06, StataMP 18, at c6b8d48 code (results: E-010).
*! Not a suite test: no pass/fail, and the runner ignores it.
*! Fits every replicate of fixtures/generate/95_marker_trap.R three ways:
*!   a  marker metric (att1 and q1 first-listed, so the markers)
*!   b  metric(unitvariance)
*!   c  marker metric with att3 and att4 reverse-coded (4 - x) in condition P;
*!      in C nothing is reversed, so c repeats a
*! and writes one row per replicate to marker_trap_<cond>.csv in the working
*! directory (run from %TEMP%; nothing here is committed).
*!
*! A fit is flagged BAD when e(converged) = 0, or its notes say "not
*! admissible" or "rank deficient" (read back from a per-fit log, the text the
*! user sees). TRAP (fixed before any data): a BAD, b not BAD, c not BAD, and
*! |T3_c - T3_b| / T3_b < 1e-6.
*!
*!     cd %TEMP%\dwlsmv_trap ; StataMP-64.exe /e do <repo>/tools/marker_trap_study.do
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
local repo "`c(pwd)'"
adopath ++ "`repo'/build"
adopath ++ "`repo'/src/ado"
mata: mata mlib index
if "$MT_NREP" == "" global MT_NREP 200
local m "(MathAb -> q1-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb)"

capture program drop mt_fit
program define mt_fit, rclass
    args model opts
    tempfile lg
    quietly log using "`lg'", text replace name(mtfit)
    capture noisily dwlsmv `model', `opts'
    local rc = _rc
    quietly log close mtfit
    local inad 0
    local rdef 0
    tempname fh
    file open `fh' using "`lg'", read text
    file read `fh' line
    while r(eof) == 0 {
        if strpos(`"`line'"', "not admissible") local inad 1
        if strpos(`"`line'"', "rank deficient") local rdef 1
        file read `fh' line
    }
    file close `fh'
    return scalar rc   = `rc'
    return scalar inad = `inad'
    return scalar rdef = `rdef'
    if `rc' == 0 {
        return scalar conv = e(converged)
        return scalar it   = e(iterations)
        return scalar T3   = e(T3)
        return scalar bad  = (e(converged) == 0 | `inad' | `rdef')
    }
    else {
        return scalar conv = .
        return scalar it   = .
        return scalar T3   = .
        return scalar bad  = 1
    }
end

foreach cond in P C {
    tempname ph
    postfile `ph' str2 cond int rep ///
        double(a_rc a_conv a_it a_inad a_rdef a_bad a_T3 a_att3 a_varA ///
               b_rc b_conv b_bad b_T3 c_rc c_conv c_bad c_T3 c_att3 c_varA) ///
        using "marker_trap_`cond'.csv.dta", replace
    forvalues r = 1/$MT_NREP {
        local f : display "`repo'/fixtures/replication/marker_trap_`cond'/marker_trap_`cond'_" %04.0f `r' ".dta"
        use "`f'", clear
        mt_fit "`m'" ""
        local a_rc = r(rc)
        local a_conv = r(conv)
        local a_it = r(it)
        local a_inad = r(inad)
        local a_rdef = r(rdef)
        local a_bad = r(bad)
        local a_T3 = r(T3)
        local a_att3 = .
        local a_varA = .
        if `a_rc' == 0 {
            local a_att3 = _b[att3:MathAtt]
            local a_varA = _b[/var(MathAtt)]
        }
        mt_fit "`m'" "metric(unitvariance)"
        local b_rc = r(rc)
        local b_conv = r(conv)
        local b_bad = r(bad)
        local b_T3 = r(T3)
        if "`cond'" == "P" {
            quietly replace att3 = 4 - att3
            quietly replace att4 = 4 - att4
        }
        mt_fit "`m'" ""
        local c_rc = r(rc)
        local c_conv = r(conv)
        local c_bad = r(bad)
        local c_T3 = r(T3)
        local c_att3 = .
        local c_varA = .
        if `c_rc' == 0 {
            local c_att3 = _b[att3:MathAtt]
            local c_varA = _b[/var(MathAtt)]
        }
        post `ph' ("`cond'") (`r') (`a_rc') (`a_conv') (`a_it') (`a_inad') ///
            (`a_rdef') (`a_bad') (`a_T3') (`a_att3') (`a_varA') (`b_rc') ///
            (`b_conv') (`b_bad') (`b_T3') (`c_rc') (`c_conv') (`c_bad') ///
            (`c_T3') (`c_att3') (`c_varA')
        display "MTRAP|`cond'|`r'|a_bad `a_bad' b_bad `b_bad' c_bad `c_bad'"
    }
    postclose `ph'
    use "marker_trap_`cond'.csv.dta", clear
    generate byte trap = a_bad == 1 & b_bad == 0 & c_bad == 0 & ///
        abs(c_T3 - b_T3) / b_T3 < 1e-6
    export delimited using "marker_trap_`cond'.csv", replace
    quietly count if trap
    local nt = r(N)
    quietly count if b_bad
    local nb = r(N)
    quietly count if a_bad
    local na = r(N)
    cii proportions $MT_NREP `nt', exact
    display "MTRAP|summary|`cond'|trap `nt' of $MT_NREP|a_bad `na'|b_bad `nb'|exact 95% CI " %6.4f r(lb) " " %6.4f r(ub)
    quietly summarize rep if trap
    display "MTRAP|first_trap|`cond'|" r(min)
}
