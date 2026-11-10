*! Module 11: dwlsmv.ado end to end (spec Module 11 D2, D5, D6, D8, D9).
*!
*! Executed under StataMP 18 on 2026-10-04 and it passes.
*!
*! Builds build/ldwlsmv.mlib first and runs against the library and src/ado,
*! as a user would; both adopath entries are removed at the end.
*!
*! Checks:
*!   1  WIRING: on case E (F2 <- F1, marker) e(b)'s free entries equal
*!      dwls_fit()'s theta-hat bit for bit, in theta order, and e(b)'s names
*!      are exactly D6's (thresholds /y:cutk, y:F, To:From, /var(F),
*!      /var(e.F)).
*!   2  D5's MEMBERSHIP CONTRACT (owed by D5): colsof(e(b)) = free + fixed;
*!      e(rank) = free; e(Cns) one row per fixed parameter selecting that
*!      column alone at its e(b) value; e(V) exactly zero in every fixed row
*!      and column and only there. Fails under e(rank) = colsof(e(b)).
*!   3  D2's INVARIANCE (owed by D2): marker and unit variance give the same
*!      T3 and df, and the same standardized loadings for an exogenous (y2)
*!      and an endogenous (y5) factor, Phi_F2F2 written out. The 1e-6
*!      tolerance was chosen after a smoke run showed T3 equal to the three
*!      decimals it printed.
*!   4  D8a/b/f, A STAGE 1 ABORT through the command: a duplicated item
*!      aborts at rc 3498, with the legend (y4 = q7) and the pointer line;
*!      under keepbound it fits, e(bound) marks exactly y1:q7, and the notes
*!      name y1~~q7.
*!   5  THE NEAR-BOUND NOTE (D8e): a pair generated at rho = .995 enters s at
*!      rho-hat >= .99 with bound = 0; the note names it by variables and
*!      says it is not a diagnosis; e(bound) is all zero.
*!   6  a Module 3 abort through the command (D3's unanchored random
*!      intercept) is rethrown at rc 3498 without Stage 1's legend, followed
*!      by the remedy line for each unanchored factor; not for a df < 0 abort;
*!      and the remedy, metric(unitvariance), fits. Check 4 also asserts the
*!      keepbound line for a status-2 duplicate and its absence for status 3
*!      (both lines decided 2026-10-04).
*!   7  refusal at rc 198: an unknown metric(). (standardized was refused
*!      here until D7 step 3, 2026-10-05; test_standardized.do covers it.)
*!   8  D9: CATEGORY RECODING -- y1 recoded 0,1,2 -> 1,2,4 gives identical
*!      e(b) values; PATH ORDER -- (F1 -> y3 y1 y2) puts /y3:cut1 first.
*!   9  replay redisplays without error.
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

local nfail 0

* run a command, capturing what it prints; r(rc), and r(found) for a phrase
capture program drop tdw_run
program tdw_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tdwlog)
    capture noisily `0'
    local rc = _rc
    quietly log close tdwlog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end

import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
local mE "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)"

* ---- 1. wiring and names ---------------------------------------------------
quietly dwlsmv `mE'
local names : colfullnames e(b)
local want "/y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 /y3:cut1 /y3:cut2 /y4:cut1 /y4:cut2 /y5:cut1 /y5:cut2 /y6:cut1 /y6:cut2 y1:F1 y2:F1 y3:F1 y4:F2 y5:F2 y6:F2 F2:F1 /var(F1) /var(e.F2)"
if "`names'" != "`want'" {
    display as error " 1. e(b) names: `names'"
    local ++nfail
}
quietly _dwlsmv_translate , paths(`mE')
local canon `"`r(canon)'"'
mata:
void tdw_wiring(string scalar canon, real scalar nfail)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    real rowvector b
    real matrix C
    real colvector free
    real scalar j

    st  = dwls_sample_stats(st_data(., tokens("y1 y2 y3 y4 y5 y6")), J(6, 1, 3))
    par = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit = dwls_fit(st, par)
    b    = st_matrix("e(b)")
    C    = st_matrix("e(Cns)")
    free = J(cols(b), 1, 1)
    for (j = 1; j <= rows(C); j++) {
        free = free :* (C[j, 1..cols(b)]' :== 0)
    }
    printf("{txt} 1. wiring: %f free entries of e(b), max |e(b) - theta-hat| %g\n",
           sum(free), max(abs(select(b', free) - fit.theta)))
    if (select(b', free) != fit.theta) {
        st_local("nfail", strofreal(nfail + 1))
    }
}
end
mata: tdw_wiring(st_local("canon"), `nfail')

* ---- 2. D5's membership contract -------------------------------------------
quietly dwlsmv `mE'
mata:
void tdw_contract(real scalar nfail)
{
    real rowvector b
    real matrix C, V
    real colvector fixed
    real scalar j, k, bad

    b = st_matrix("e(b)")
    C = st_matrix("e(Cns)")
    V = st_matrix("e(V)")
    bad = 0
    fixed = J(cols(b), 1, 0)
    for (j = 1; j <= rows(C); j++) {
        if (sum(C[j, 1..cols(b)] :!= 0) != 1) bad = 1
        k = selectindex(C[j, 1..cols(b)] :!= 0)
        if (C[j, k] != 1 | C[j, cols(C)] != b[k]) bad = 1
        fixed[k] = 1
    }
    if (st_numscalar("e(rank)") != cols(b) - rows(C)) bad = 1
    if (sum(fixed) != rows(C)) bad = 1
    for (j = 1; j <= cols(b); j++) {
        if (fixed[j]) {
            if (any(V[j, .] :!= 0) | any(V[., j] :!= 0)) bad = 1
        }
        else if (!(V[j, j] > 0)) bad = 1
    }
    printf("{txt} 2. D5 contract: %f columns, rank %f, %f Cns rows: %s\n",
           cols(b), st_numscalar("e(rank)"), rows(C), (bad ? "VIOLATED" : "holds"))
    if (bad) st_local("nfail", strofreal(nfail + 1))
}
end
mata: tdw_contract(`nfail')

* ---- 3. D2's invariance -------------------------------------------------------
quietly dwlsmv `mE'
local T3m = e(T3)
local dfm = e(df_T3)
local s2m = _b[y2:F1] * sqrt(_b[/var(F1)])
local s5m = _b[y5:F2] * sqrt(_b[F2:F1]^2 * _b[/var(F1)] + _b[/var(e.F2)])
quietly dwlsmv `mE', metric(unitvariance)
local T3u = e(T3)
local s2u = _b[y2:F1] * sqrt(_b[/var(F1)])
local s5u = _b[y5:F2] * sqrt(_b[F2:F1]^2 * _b[/var(F1)] + _b[/var(e.F2)])
display as text " 3. D2: T3 " %12.8f `T3m' " vs " %12.8f `T3u' "; std y2 " %9.6f `s2m' " vs " %9.6f `s2u' ///
    "; std y5 " %9.6f `s5m' " vs " %9.6f `s5u'
if abs(`T3m' - `T3u') > 1e-6 | `dfm' != e(df_T3) | abs(`s2m' - `s2u') > 1e-6 | abs(`s5m' - `s5u') > 1e-6 {
    display as error " 3. the metrics disagree"
    local ++nfail
}

* ---- 4. a Stage 1 abort, and keepbound ------------------------------------------
generate q7 = y1
tdw_run "y4 = q7" dwlsmv (F1 -> y1 y2 y3 q7) (F2 -> y4 y5 y6)
local rc4 = r(rc)
local leg = r(found)
tdw_run "Comparison with lavaan and Mplus" dwlsmv (F1 -> y1 y2 y3 q7) (F2 -> y4 y5 y6)
display as text " 4. Stage 1 abort: rc `rc4', legend `leg', pointer " r(found)
if `rc4' != 3498 | !`leg' | !r(found) {
    local ++nfail
}
* the keepbound line: printed where keepbound would pass, not for status 3
tdw_run "is the dwlsmv option keepbound" dwlsmv (F1 -> y1 y2 y3 q7) (F2 -> y4 y5 y6)
local kl1 = r(found)
set seed 99
generate q8 = cond(runiform() < .005, mod(y1 + 1, 3), y1)
tdw_run "is the dwlsmv option keepbound" dwlsmv (F1 -> y1 y2 y3 q8) (F2 -> y4 y5 y6)
display as text " 4. keepbound line: duplicate (status 2) " `kl1' ", near-duplicate (status 3) rc " r(rc) " line " r(found)
if !`kl1' | r(rc) != 3498 | r(found) {
    local ++nfail
}
tdw_run "y1~~q7" dwlsmv (F1 -> y1 y2 y3 q7) (F2 -> y4 y5 y6), keepbound
local rc4 = r(rc)
local foot = r(found)
tempname B
matrix `B' = e(bound)
local nb = 0
forvalues k = 1/`=colsof(`B')' {
    local nb = `nb' + `B'[1, `k']
}
local at = colnumb(`B', "y1:q7")
display as text " 4. keepbound: rc `rc4', footnote `foot', e(bound) sum `nb', y1:q7 = " `B'[1, `at']
if `rc4' != 0 | !`foot' | `nb' != 1 | `B'[1, `at'] != 1 {
    local ++nfail
}

* ---- 6. a Module 3 abort through the command --------------------------------------
tdw_run "Positional labels" dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1), cov(G*F1@0) cov(G*F2@0)
display as text " 6. Module 3 abort: rc " r(rc) ", legend printed " r(found)
if r(rc) != 3498 | r(found) {
    local ++nfail
}
* the remedy line: for each unanchored factor, and not for another Module 3
* abort (more free parameters than statistics, every factor anchored)
tdw_run "Factor F2 has no fixed scale" dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1), cov(G*F1@0) cov(G*F2@0)
local rm1 = r(found)
tdw_run "(F1 -> y1@1 ...), or specify metric(unitvariance)" dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1), cov(G*F1@0) cov(G*F2@0)
local rm2 = r(found)
tdw_run "has no fixed scale" dwlsmv (F -> y1 y2)
display as text " 6. remedy line: unanchored F2 `rm1', F1 example `rm2'; df < 0 abort rc " r(rc) " line " r(found)
if !`rm1' | !`rm2' | r(rc) != 3498 | r(found) {
    local ++nfail
}
* the remedy works: the same model under metric(unitvariance) fits, which
* also takes it through the rank check's second point (spec 4.8)
capture noisily quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1), cov(G*F1@0) cov(G*F2@0) metric(unitvariance)
display as text " 6. remedy applied, metric(unitvariance): rc " _rc ", converged " e(converged)
if _rc | e(converged) != 1 {
    local ++nfail
}

* ---- 7. refusal at rc 198 ------------------------------------------------------------
* standardized is accepted since D7 step (3); test_standardized.do covers it
capture noisily dwlsmv `mE', metric(total)
display as text " 7. metric(total) rc " _rc
if _rc != 198 {
    local ++nfail
}

* ---- 8. D9: recoding and path order ----------------------------------------------------
quietly dwlsmv `mE'
tempname b1
matrix `b1' = e(b)
generate w = cond(y1 == 2, 4, y1 + 1)
quietly dwlsmv (F1 -> w y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)
tempname b2
matrix `b2' = e(b)
display as text " 8. recoding 0,1,2 -> 1,2,4: max |e(b) difference| " mreldif(`b1', `b2')
if mreldif(`b1', `b2') != 0 {
    local ++nfail
}
quietly dwlsmv (F1 -> y3 y1 y2) (F2 -> y4 y5 y6)
local first : word 1 of `: colfullnames e(b)'
display as text " 8. path order (F1 -> y3 y1 y2): first name `first'"
if "`first'" != "/y3:cut1" {
    local ++nfail
}

* ---- 9. replay ------------------------------------------------------------------------
capture noisily dwlsmv
display as text " 9. replay rc " _rc
if _rc {
    local ++nfail
}

* ---- 5. the near-bound note (fresh data) -------------------------------------------------
clear
set seed 1
quietly set obs 5000
generate double f = rnormal()
forvalues j = 1/4 {
    generate double z`j' = .7*f + sqrt(1 - .49)*rnormal()
}
generate double zd = .995*z1 + sqrt(1 - .995^2)*rnormal()
foreach v in z1 z2 z3 z4 zd {
    generate byte c`v' = (`v' > -.3) + (`v' > .65)
}
tdw_run "cz1~~czd" dwlsmv (F -> cz1 cz2 cz3 cz4 czd)
local rc5 = r(rc)
local nm5 = r(found)
tdw_run "this is not a diagnosis" dwlsmv (F -> cz1 cz2 cz3 cz4 czd)
matrix `B' = e(bound)
local nb = 0
forvalues k = 1/`=colsof(`B')' {
    local nb = `nb' + `B'[1, `k']
}
display as text " 5. near-bound note: rc `rc5', pair named `nm5', reading " r(found) ", e(bound) sum `nb'"
if `rc5' != 0 | !`nm5' | !r(found) | `nb' != 0 {
    local ++nfail
}

* ---- 10. the one-empty-cell repair through the command (O6, 2026-10-06) -------------------
mata:
// Stage 1's repaired rho for y1~~y2 of b1-b4, formatted as the note prints it.
void dwlsmv_t_dw_r10()
{
    struct dwlsmv_stage1 scalar st
    st = dwls_sample_stats(st_data(., "b1 b2 b3 b4"), (2 \ 2 \ 2 \ 2))
    st_local("r10", strtrim(strofreal(st.s[st.map.T + 1], "%9.4f")))
}
end
* b1 rare, b2 forced to 1 where b1 is 1: the b1~~b2 table has exactly one empty
* cell. The fit runs (rc 0), the notes block names the pair by variable, and
* e(bound) stays all zero (the repair is not keepbound).
clear
set seed 6
quietly set obs 400
generate double f = rnormal()
forvalues j = 1/4 {
    generate double z`j' = .7*f + sqrt(1 - .49)*rnormal()
}
generate byte b1 = (z1 > 1)
generate byte b2 = (z2 > -.5)
quietly replace b2 = 1 if b1 == 1
generate byte b3 = (z3 > 0)
generate byte b4 = (z4 > .3)
quietly count if b1 == 1 & b2 == 0
local e10 = r(N)
* the note's rho must be Stage 1's repaired value for THIS pair, to the note's
* %9.4f (adversarial review 2026-10-06: a wrong pair offset kept the name)
mata: dwlsmv_t_dw_r10()
tdw_run "    b1~~b2   rho = `r10'" dwlsmv (F -> b1 b2 b3 b4)
local rv10 = r(found)
tdw_run "One-empty-cell repair applied to:" dwlsmv (F -> b1 b2 b3 b4)
local rc10 = r(rc)
local hd10 = r(found)
tdw_run "    b1~~b2   rho = " dwlsmv (F -> b1 b2 b3 b4)
local nm10 = r(found)
tdw_run "    b1~~b3   rho = " dwlsmv (F -> b1 b2 b3 b4)
local nx10 = r(found)
matrix `B' = e(bound)
local nb = 0
forvalues k = 1/`=colsof(`B')' {
    local nb = `nb' + `B'[1, `k']
}
display as text "10. repair note: empty cell count `e10', rc `rc10', heading `hd10', b1~~b2 `nm10', rho `r10' `rv10', b1~~b3 `nx10', e(bound) sum `nb'"
if `e10' != 0 | `rc10' != 0 | !`hd10' | !`nm10' | !`rv10' | `nx10' | `nb' != 0 {
    local ++nfail
}

* ---- 11. a path fixed at 0 keeps its factor endogenous ----------------------
* sem posts /var(e.F2) and no F2:F1 column; the parameter map has no beta row.
* D6 must take the factor's status from sem's stripe, and D2's injection from
* sem's e(lyvars). Same model as cov(F1*F2@0) numerically, so T3 must agree
* (1e-6, the tolerance of check 3).
import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
capture noisily quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1@0)
local rc11 = _rc
local ok11 = 0
if !`rc11' {
    local n11 : colfullnames e(b)
    local T11 = e(T3)
    local df11 = e(df_T3)
    local ok11 = (strpos(" `n11' ", " /var(e.F2) ") > 0) & (strpos(" `n11' ", " F2:F1 ") == 0)
    quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), cov(F1*F2@0)
    local d11 = abs(`T11' - e(T3))
    local ok11 = `ok11' & (`d11' < 1e-6) & (`df11' == e(df_T3))
    capture noisily quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1@0), metric(unitvariance)
    local rcu11 = _rc
    local ok11 = `ok11' & (`rcu11' == 0)
    if !`rcu11' {
        local ok11 = `ok11' & (abs(e(T3) - `T11') < 1e-6)
    }
}
display as text "11. fixed-0 path: rc `rc11', names and T3 match the cov() fit: `ok11'"
if `rc11' | !`ok11' {
    local ++nfail
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_dwlsmv: `nfail' failure(s)"
    exit 9
}
display as result "test_dwlsmv: e(b) is dwls_fit's theta-hat under D6's names, D5's contract holds, the metrics agree on T3 and standardized loadings, Stage 1 aborts carry the legend and pointer, keepbound and the near-bound note report by variable, and recoding and path order follow D9"
exit 0
