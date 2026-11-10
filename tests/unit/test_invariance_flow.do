*! Module 11: the invariance() flow up to the fit (group surface, build step
*! 3, 2026-10-12). Executed under StataMP 18 on 2026-10-12 and it passes.
*!
*! Decisions under test (E-073, 2026-10-12): D9 under groups (b),
*! supports compared past configural, and its wording (a); F4 (a), equal
*! categories before the binary scalar refusal; T0/T2, the binary scalar
*! refusal rc 198 before any fit, message G-A + E1; no configural fit inside
*! invariance(level); the pooled fit (weighting (a)).
*!
*! Checks (data fixtures/data/mg_E.csv, groups g = 1, 2; F -> y1..y5):
*!   1  equal supports go through the existing path unchanged: at each level
*!      e(canon) equals dwlsmv_levelgen on the configural text, and e(T3) and
*!      the free values of e(b) equal the direct pooled fit on hand-built
*!      Stage 1 (test_fit_groups.do's call) exactly: same integers in, same
*!      code path, so no tolerance applies.
*!   2  equal counts, different supports (y1 = 3 recoded 4 in group 2; no
*!      committed fixture has the case, checked 2026-10-12): 2a configural is
*!      not refused and K is 4 in both groups; 2b-2d thresholds, metric and
*!      scalar refuse at rc 198 with the values line and the closing sentence,
*!      never the count line, before Stage 1 and before any fit.
*!   3  F4 order: y6 binary in group 1, three categories in group 2, in an
*!      otherwise all-binary factor under invariance(scalar): the count line,
*!      never the binary refusal.
*!   4  the binary scalar refusal: an all-binary factor under scalar, rc 198,
*!      G-A's first line and E1, before Stage 1; 4b the same factor under
*!      invariance(thresholds) is not refused.
*!   5  M5 (owed by step 2): three groups, one vanishing after markout; the
*!      values are re-read and the text retranslated: two groups fitted.
*!   6  the d = 0 note (text decided 2026-10-12, user): 6a Result I2's case
*!      (F -> y1 y4 y5, one polytomous item) under metric prints the
*!      equivalence wording; e(df_T3) is equal at both levels and e(T3) is the
*!      direct fit's; the direct fits' T_raw agree within check 10's band
*!      (N - G) nrtol; the T3 difference is printed, not asserted; 6b the thresholds step with 3-category items
*!      only prints the equal-df wording, never the equivalence one; 6c the
*!      verified-case detector on hand-written texts (its false branches are
*!      unreachable through the command: every loading there is free); 6d no
*!      note when d > 0.
*!
*! A run reaching the fit posts (rc 0, e(cmd) dwlsmv; step 4, 2026-10-12);
*! tif_run clears e() first. "Before Stage 1" holds by code order (the checks
*! run before dwlsmv_m11_stage1_g); the test sees only that nothing was posted
*! and no Stage 1 state is left, which every refusal path also guarantees by
*! clearing, so the order itself is not observable here.
*!
*! Mutants (2026-10-12, scratch worktree): F4 order swapped (binary check
*! first) red at 3; supports compared by count only red at 2b-2d (and the
*! mis-tied model then fits without a word); M5, group values not re-read
*! after markout, red at 5 (first as an abort, r(3010), in this test's own
*! check; the check now reads Mata state only after the fit is reached); the
*! d = 0 note never printed, red at 6a, 6b; the equivalence wording for every
*! d = 0, red at 6b.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
mata: YALL = st_data(., .)

local nfail 0
local m5 (F -> y1 y2 y3 y4 y5)

* run a dwlsmv command line after clearing Module 11's state; rc into tif_rc,
* the log text into the Mata global TIFLOG
capture program drop tif_run
program tif_run
    args cmd
    capture mata: dwlsmv_m11_clear()
    ereturn clear
    tempfile lg
    capture log close tif
    quietly log using "`lg'", text replace name(tif)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tif
    mata: TIFLOG = invtokens(cat(st_local("lg"))')
    c_local tif_rc `rc'
    c_local tif_reached = (`rc' == 0 & "`e(cmd)'" == "dwlsmv")
end

* has / lacks: 1 if the text is (is not) in the last log
capture program drop tif_has
program tif_has
    args s
    mata: st_local("h", strofreal(strpos(TIFLOG, st_local("s")) > 0))
    c_local hit `h'
end

mata:

string scalar tif_conf()
{
    string scalar s, nl
    real scalar g, j

    nl = char(10)
    s  = "groups 2" + nl
    for (g = 1; g <= 2; g++) {
        s = s + "group " + strofreal(g) + nl + "factor F" + nl + "ind y1 y2 y3 y4 y5" + nl
        for (j = 1; j <= 5; j++) {
            s = s + "lambda F y" + strofreal(j) + " free" + nl
        }
        s = s + (g == 1 ? "metric unitvariance F" : "psi F F fixed 1") + nl
    }
    return(s)
}

// the direct pooled fit (test_fit_groups.do's call): Stage 1 per group on
// YALL's columns cl (g in column 1), K per group KE, the level generated from
// the configural text conf
struct dwlsmv_fit scalar tif_dfit(string scalar conf, real rowvector cl,
                                  real colvector KE, string scalar level)
{
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    real colvector sel
    real scalar g

    sts  = dwlsmv_stage1(2, 1)
    maps = dwlsmv_statmap(2, 1)
    for (g = 1; g <= 2; g++) {
        sel     = selectindex(YALL[., 1] :== g)
        sts[g]  = dwls_sample_stats(YALL[sel, cl], KE)
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(conf, J(2, 1, KE'), level), maps, sts)
    fit = dwls_fit(sts, par)
    dwls_inference(fit)
    return(fit)
}

// 0 if the posted fit equals the direct pooled fit at the level exactly:
// e(canon) is the generated text; e(T3) is the fit's; the free columns of
// e(b) (nonzero e(V) diagonal) are exactly the set of theta-hat's values
real scalar tif_same(string scalar level)
{
    struct dwlsmv_fit scalar fit
    real rowvector b, fb
    real colvector KE
    real scalar bad, i

    KE  = (4 \ 4 \ 4 \ 2 \ 2)
    fit = tif_dfit(tif_conf(), (2..6), KE, level)
    bad = 0
    if (st_global("e(canon)") != dwlsmv_levelgen(tif_conf(), J(2, 1, KE'), level)) {
        printf("{err}    e(canon) is not the generated text at %s\n", level)
        bad = 1
    }
    if (st_numscalar("e(T3)") != fit.T3) {
        printf("{err}    e(T3) %21.0g differs from the direct fit's %21.0g at %s\n",
               st_numscalar("e(T3)"), fit.T3, level)
        bad = 1
    }
    b  = st_matrix("e(b)")
    fb = select(b, diagonal(st_matrix("e(V)"))' :!= 0)
    for (i = 1; i <= cols(fb); i++) {
        if (!anyof(fit.theta, fb[i])) {
            printf("{err}    free e(b) value %21.0g is not in theta-hat at %s\n", fb[i], level)
            bad = 1
        }
    }
    for (i = 1; i <= rows(fit.theta); i++) {
        if (!anyof(fb, fit.theta[i])) {
            printf("{err}    theta-hat value %21.0g is not posted at %s\n", fit.theta[i], level)
            bad = 1
        }
    }
    return(bad)
}

// keep e(df_T3) and e(T3) for 6a
void tif_keep()
{
    external real rowvector TIF0

    TIF0 = (st_numscalar("e(df_T3)"), st_numscalar("e(T3)"))
}

// 6a: the posted df equal at both levels and e(T3) the direct fit's exactly;
// the direct fits' T_raw within check 10's band (N - G) nrtol, nrtol = 1e-12
// (dwlsmv_minimize); the T3 difference printed, not asserted
real scalar tif_band()
{
    external real rowvector TIF0
    struct dwlsmv_fit scalar f0, f1
    string scalar conf, nl
    real scalar N, ok

    nl   = char(10)
    conf = "groups 2" + nl + "group 1" + nl + "factor F" + nl + "ind y1 y4 y5" + nl +
           "lambda F y1 free" + nl + "lambda F y4 free" + nl + "lambda F y5 free" + nl +
           "metric unitvariance F" + nl + "group 2" + nl + "factor F" + nl + "ind y1 y4 y5" + nl +
           "lambda F y1 free" + nl + "lambda F y4 free" + nl + "lambda F y5 free" + nl +
           "psi F F fixed 1" + nl
    f0 = tif_dfit(conf, (2, 5, 6), (4 \ 2 \ 2), "thresholds")
    f1 = tif_dfit(conf, (2, 5, 6), (4 \ 2 \ 2), "metric")
    N  = f1.stage1[1].N + f1.stage1[2].N
    printf("{txt}    6a T_raw difference %g, band %g; T3 difference %g (printed, not asserted)\n",
           f1.T_raw - f0.T_raw, (N - 2) * 1e-12, f1.T3 - f0.T3)
    ok = abs(f1.T_raw - f0.T_raw) <= (N - 2) * 1e-12
    ok = ok & st_numscalar("e(df_T3)") == TIF0[1] & st_numscalar("e(T3)") == f1.T3 & TIF0[2] == f0.T3
    return(ok)
}

// 1 iff dwlsmv_m11_fit_g_i2 gives 1, 0, 0 on: one polytomous item (y1);
// two (y1, y2); one, with a loading fixed at 0.5
real scalar tif_i2()
{
    string scalar a, b, c, nl
    real matrix K

    nl = char(10)
    K  = J(2, 1, (4, 4, 2, 2))
    a  = "groups 2" + nl + "group 1" + nl + "factor F" + nl + "ind y1 y2 y3 y4" + nl
    b  = a + "lambda F y1 free" + nl + "lambda F y2 free" + nl + "lambda F y3 free" + nl +
         "lambda F y4 free" + nl + "group 2" + nl
    c  = a + "lambda F y1 free" + nl + "lambda F y3 fixed 0.5" + nl + "lambda F y4 free" + nl +
         "group 2" + nl
    a  = a + "lambda F y1 free" + nl + "lambda F y3 free" + nl + "lambda F y4 free" + nl +
         "group 2" + nl
    return(dwlsmv_m11_fit_g_i2(a, K) == 1 & dwlsmv_m11_fit_g_i2(b, K) == 0 &
           dwlsmv_m11_fit_g_i2(c, K) == 0)
}

// 1 if Stage 1 and the fit were never reached
real scalar tif_nofit()
{
    // nothing posted (tif_run cleared e()) and no Stage 1 left behind
    return(st_global("e(cmd)") == "" & findexternal("dwlsmv_m11_st") == NULL)
}

end

display as text "1. equal supports, the existing path unchanged"
foreach lev in configural thresholds metric scalar {
    tif_run `"dwlsmv `m5', group(g) invariance(`lev') metric(unitvariance)"'
    local hit `tif_reached'
    local ok = (`hit')
    if `ok' {
        mata: st_local("d", strofreal(tif_same("`lev'")))
        local ok = (`d' == 0)
    }
    if `ok' display as text "  1 `lev': ok"
    else {
        display as error "  1 `lev' FAILED (rc `tif_rc', fit reached `hit')"
        local ++nfail
    }
}

display as text "2. equal counts, different supports"
preserve
quietly replace y1 = 4 if y1 == 3 & g == 2
tif_run `"dwlsmv `m5', group(g) invariance(configural) metric(unitvariance)"'
local hit `tif_reached'
local kok 0
if `hit' {
    mata: st_local("kok", strofreal(st_matrix("e(K)") == J(2, 1, (4, 4, 4, 2, 2))))
}
if `hit' & `kok' display as text "  2a configural fitted, K 4 in both groups: ok"
else {
    display as error "  2a FAILED (rc `tif_rc', fit reached `hit', K/support `kok')"
    local ++nfail
}
local vline "y1 takes the values 0 1 2 3 in group 1 (g = 1) and 0 1 2 4 in group 2 (g = 2)."
local close "across groups, which needs the same categories in every group."
local i 1
foreach lev in thresholds metric scalar {
    local lab = word("b c d", `i')
    local ++i
    tif_run `"dwlsmv `m5', group(g) invariance(`lev') metric(unitvariance)"'
    tif_has "`vline'"
    local h1 `hit'
    tif_has "`close'"
    local h2 `hit'
    tif_has "y1 has"
    local h3 `hit'
    local hit `tif_reached'
    local h4 `hit'
    mata: st_local("nf", strofreal(tif_nofit()))
    if `tif_rc' == 198 & `h1' & `h2' & !`h3' & !`h4' & `nf' {
        display as text "  2`lab' `lev' refused, nothing posted: ok"
    }
    else {
        display as error "  2`lab' `lev' FAILED (rc `tif_rc'; values `h1' closing `h2' count `h3' fit `h4' nofit `nf')"
        local ++nfail
    }
}
restore

display as text "3. F4: equal categories before the binary refusal"
preserve
quietly generate byte y6 = y4
quietly replace y6 = 2 if g == 2 & y4 == 1 & y5 == 1
tif_run `"dwlsmv (F -> y4 y5 y6), group(g) invariance(scalar) metric(unitvariance)"'
tif_has "y6 has 2 categories in group 1 (g = 1) and 3 in group 2 (g = 2)."
local h1 `hit'
tif_has "every indicator has two categories"
local h2 `hit'
mata: st_local("nf", strofreal(tif_nofit()))
if `tif_rc' == 198 & `h1' & !`h2' & `nf' display as text "  3: ok"
else {
    display as error "  3 FAILED (rc `tif_rc'; count line `h1', binary refusal `h2', nofit `nf')"
    local ++nfail
}
restore

display as text "4. the binary scalar refusal"
preserve
quietly generate byte y7 = (y1 >= 2)
tif_run `"dwlsmv (F -> y4 y5 y7), group(g) invariance(scalar) metric(unitvariance)"'
tif_has "Factor(s) F: every indicator has two categories."
local h1 `hit'
tif_has "invariance(scalar) is not available for an all-binary factor in this version"
local h2 `hit'
mata: st_local("nf", strofreal(tif_nofit()))
if `tif_rc' == 198 & `h1' & `h2' & `nf' display as text "  4a: ok"
else {
    display as error "  4a FAILED (rc `tif_rc'; G-A `h1', E1 `h2', nofit `nf')"
    local ++nfail
}
tif_run `"dwlsmv (F -> y4 y5 y7), group(g) invariance(thresholds) metric(unitvariance)"'
local hit `tif_reached'
if `hit' display as text "  4b thresholds not refused: ok"
else {
    display as error "  4b FAILED (rc `tif_rc', fit reached `hit')"
    local ++nfail
}
restore

display as text "5. M5: a group vanishing among three after markout"
preserve
quietly generate byte g3 = mod(_n, 3)
quietly replace y1 = . if g3 == 1
tif_run `"dwlsmv `m5', group(g3) invariance(configural) metric(unitvariance)"'
local hit `tif_reached'
* read only once the fit was posted
local k2 0
if `hit' {
    mata: st_local("k2", strofreal(rows(st_matrix("e(K)")) == 2 & ///
        substr(st_global("e(canon)"), 1, 9) == "groups 2" + char(10)))
}
if `hit' & `k2' display as text "  5: two groups fitted: ok"
else {
    display as error "  5 FAILED (rc `tif_rc', fit reached `hit', two groups `k2')"
    local ++nfail
}
restore

display as text "6. the d = 0 note"
local vtxt "this step adds no testable restriction"
local utxt "Equal degrees of freedom do not by themselves make two"
preserve
quietly generate byte y8  = min(y1, 2)
quietly generate byte y9  = min(y2, 2)
quietly generate byte y10 = min(y3, 2)
* 6a verified (Result I2): one polytomous item; metric against thresholds,
* unscaled T within check 10's band (N - G) nrtol, nrtol = 1e-12
* (dwlsmv_minimize), df exact; the T3 difference printed, not asserted
tif_run `"dwlsmv (F -> y1 y4 y5), group(g) invariance(thresholds) metric(unitvariance)"'
mata: tif_keep()
tif_run `"dwlsmv (F -> y1 y4 y5), group(g) invariance(metric) metric(unitvariance)"'
tif_has "`vtxt'"
local h1 `hit'
tif_has "`utxt'"
local h2 `hit'
mata: st_local("band", strofreal(tif_band()))
if `h1' & !`h2' & `band' display as text "  6a verified note, T and df equal: ok"
else {
    display as error "  6a FAILED (verified note `h1', equal-df note `h2', band and df `band')"
    local ++nfail
}
* 6b not verified: the thresholds step with 3-category items only
tif_run `"dwlsmv (F -> y8 y9 y10), group(g) invariance(thresholds) metric(unitvariance)"'
tif_has "Note: invariance(thresholds) has the same degrees of freedom as invariance(configural)"
local h1 `hit'
tif_has "`utxt'"
local h2 `hit'
tif_has "`vtxt'"
local h3 `hit'
if `h1' & `h2' & !`h3' display as text "  6b equal-df note: ok"
else {
    display as error "  6b FAILED (`h1' `h2' `h3')"
    local ++nfail
}
* 6c dwlsmv_m11_fit_g_i2 directly: through the command every loading is free
* (a loading fixed at a value other than 1 is refused upstream, and a 1 is
* refused under invariance()), so its false branches are reachable only from
* a hand-written text: one polytomous item 1; two polytomous items 0; a fixed
* loading 0
mata: st_local("i2", strofreal(tif_i2()))
if "`i2'" == "1" display as text "  6c the verified-case detector: ok"
else {
    display as error "  6c FAILED: dwlsmv_m11_fit_g_i2"
    local ++nfail
}
* 6d d > 0: no note
tif_run `"dwlsmv `m5', group(g) invariance(metric) metric(unitvariance)"'
tif_has "Note: invariance("
if !`hit' display as text "  6d no note when d > 0: ok"
else {
    display as error "  6d FAILED: a d = 0 note at d > 0"
    local ++nfail
}
restore

capture mata: dwlsmv_m11_clear()
adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_invariance_flow: `nfail' failure(s)"
    exit 9
}
display as result "test_invariance_flow: equal supports fit as before; differing supports and all-binary scalar factors are refused before any fit, in F4's order; the d = 0 note separates equal df from established equivalence"
exit 0
