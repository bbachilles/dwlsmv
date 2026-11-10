*! Module 11: the group surface's translation and refusals (build step 2,
*! 2026-10-12). Executed under StataMP 18 on 2026-10-12 and it passes.
*!
*! Decisions under test (E-073, 2026-10-12 unless noted): group()
*! without invariance() refused as a version limit; invariance() needs group()
*! and two groups in the estimation sample; invariance() requires an explicit
*! metric(unitvariance), and a loading fixed at 1 is refused (two messages);
*! group-specific paths: every group mirrors group 1, checked in Module 11 at
*! rc 198 before the generator; partial() and ginvariant() refused (2026-10-10);
*! O13's group values (nonnegative integers, a string variable rc 109).
*!
*! Checks:
*!   1  translation: the configural text exactly, one factor and two factors
*!      with a cross-loading (probe case 5.1); 1c the same model written group
*!      by group (case 5.2) gives the identical text; 1d three groups named by
*!      the user's values; 1e the texts meet dwlsmv_levelgen's input contract.
*!   2  translation refusals, rc 198 and the decided message: a loading fixed
*!      at 1; a structural path; group-specific paths (cases 5.3, 5.4, 5.6);
*!      a cross-group tie; a constrained factor mean; 2i precedence (a fixed
*!      loading at 1 is reported before a structural path); 2j, 2k message 8,
*!      a status difference (a covariance, a loading fixed in group 2 only).
*!   3  dwlsmv's surface refusals, rc 198 (109 for a string group variable)
*!      and the decided message; 3l one group left only after the indicators
*!      are marked out.
*!   4  state restoration after a group translation and a group refusal.
*!
*! Mutants (2026-10-12, scratch worktree): M1 structural paths reported before
*! a loading at 1, M2 mirroring check skipped, M3 injection counting a loading
*! fixed at 0 as a normalisation (the defect found on the first run), M4 the
*! stripe regex reading the wrong capture group: all red. M5, the group values
*! not re-read after markout, SURVIVES: nothing after that point runs until
*! build step 3, which owes a test where a group vanishes among three.
*!
*! Expected texts follow tools/probe_sem_group_surface.do (StataMP 18).
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

import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
quietly generate byte g = mod(_n, 2)
quietly generate byte g3 = mod(_n, 3) * 4 + 1

local nfail 0

capture program drop ttg_expect
program ttg_expect
    args label model gvals expect
    capture noisily _dwlsmv_translate , paths(`model') metric(unitvariance) ///
        group(g) gvals(`gvals')
    if _rc {
        display as error "  `label': rc " _rc
        c_local tt_bad 1
        exit
    }
    if `"`r(canon)'"' != `"`expect'"' {
        display as error "  `label': canonical text differs"
        display as error "    got:  `r(canon)'"
        display as error "    want: `expect'"
        c_local tt_bad 1
        exit
    }
    display as text "  `label': ok (`r(injected)')"
    c_local tt_bad 0
    c_local tt_canon `"`r(canon)'"'
end

* a refusal: the command (a full command line) must exit with rc and the
* message must appear in the log
capture program drop ttg_refuse
program ttg_refuse
    args label rcwant msg cmd
    tempfile lg
    capture log close ttg
    quietly log using "`lg'", text replace name(ttg)
    capture noisily `cmd'
    local rc = _rc
    quietly log close ttg
    mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("lg"))'), st_local("msg")) > 0))
    if `rc' != `rcwant' | !`hit' {
        display as error "  `label' FAILED: rc `rc' (wanted `rcwant'), message found `hit'"
        c_local tt_bad 1
        exit
    }
    display as text "  `label': ok"
    c_local tt_bad 0
end

display as text "1. translation"
local A1 "groups 2;group 1;factor F;ind y1 y2 y3;lambda F y1 free;lambda F y2 free;lambda F y3 free;metric unitvariance F"
local A2 "group 2;factor F;ind y1 y2 y3;lambda F y1 free;lambda F y2 free;lambda F y3 free;psi F F fixed 1"
ttg_expect "1a one factor" "(F -> y1 y2 y3)" "0 1" "`A1';`A2'"
local nfail = `nfail' + `tt_bad'

local B1 "groups 2;group 1;factor F G;ind y1 y2 y3 y4 y5;lambda F y1 free;lambda F y2 free;lambda F y3 free;lambda G y3 free;lambda G y4 free;lambda G y5 free;metric unitvariance F;metric unitvariance G;psi G F free"
local B2 "group 2;factor F G;ind y1 y2 y3 y4 y5;lambda F y1 free;lambda F y2 free;lambda F y3 free;lambda G y3 free;lambda G y4 free;lambda G y5 free;psi F F fixed 1;psi G F free;psi G G fixed 1"
ttg_expect "1b cross-loading (5.1)" "(F -> y1 y2 y3) (G -> y3 y4 y5)" "0 1" "`B1';`B2'"
local nfail = `nfail' + `tt_bad'
ttg_expect "1c group by group (5.2)" "(0: F -> y1 y2 y3) (1: F -> y1 y2 y3) (G -> y3 y4 y5)" "0 1" "`B1';`B2'"
local nfail = `nfail' + `tt_bad'

* three groups named by the user's values 1, 5, 9 (O13): the text numbers
* groups 1..G; the values reach only the names and messages
capture noisily _dwlsmv_translate , paths((F -> y1 y2 y3)) metric(unitvariance) ///
    group(g3) gvals(1 5 9)
local want "groups 3;group 1;factor F;ind y1 y2 y3;lambda F y1 free;lambda F y2 free;lambda F y3 free;metric unitvariance F"
foreach k in 2 3 {
    local want "`want';group `k';factor F;ind y1 y2 y3;lambda F y1 free;lambda F y2 free;lambda F y3 free;psi F F fixed 1"
}
if _rc | `"`r(canon)'"' != `"`want'"' {
    display as error "  1d three groups FAILED: rc " _rc
    display as error "    got:  `r(canon)'"
    local ++nfail
}
else display as text "  1d three groups: ok"

* the configural texts meet the generator's input contract (L1 (b))
foreach m in A B {
    local t = subinstr("``m'1';``m'2'", ";", char(10), .) + char(10)
    local p = cond("`m'" == "A", 3, 5)
    capture noisily mata: st_local("o", dwlsmv_levelgen(st_local("t"), J(2, `p', 4), "configural"))
    if _rc {
        display as error "  1e `m' FAILED: dwlsmv_levelgen rc " _rc
        local ++nfail
    }
    else display as text "  1e `m': dwlsmv_levelgen accepts the text"
}

display as text "2. translation refusals"
local tr _dwlsmv_translate , metric(unitvariance) group(g) gvals(0 1) paths
ttg_refuse "2a loading at 1" 198 "Factor F has a loading fixed at 1 (F -> y1@1)" ///
    `"`tr'((F -> y1@1 y2 y3))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2b structural path" 198 "Structural paths between factors (H <- F) are not available" ///
    `"`tr'((F -> y1 y2 y3) (H -> y4 y5 y6) (H <- F))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2c group 2 only (5.4)" 198 "Path F -> y4 is in group 2 (g = 1) but not in group 1 (g = 0)." ///
    `"`tr'((F -> y1 y2 y3) (G -> y3 y4 y5) (1: F -> y4))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2d marker moved (5.3)" 198 "Path G -> y3 is in group 2 (g = 1) but not in group 1 (g = 0)." ///
    `"`tr'((F -> y1 y2 y3) (G -> y4 y5) (1: G -> y3))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2e fixed in group 2 (5.6)" 198 "Path F -> y4 is in group 2 (g = 1) but not in group 1 (g = 0)." ///
    `"`tr'((F -> y1 y2 y3) (G -> y3 y4 y5) (1: F -> y4@0.3))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2f cross-group tie" 198 "equal across groups. No model was fitted." ///
    `"`tr'((0: F -> y1 y2@a y3) (1: F -> y1 y2@a y3))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2g factor mean" 198 "where they are estimated in groups after the first. No model was fitted." ///
    `"_dwlsmv_translate , metric(unitvariance) group(g) gvals(0 1) paths((F -> y1 y2 y3)) semopts(mean(F@0.5))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2i precedence" 198 "Factor F has a loading fixed at 1 (F -> y1@1)" ///
    `"`tr'((F -> y1@1 y2 y3) (H -> y4 y5 y6) (H <- F))"'
local nfail = `nfail' + `tt_bad'
* message 8, a status difference (owed by step 2; sem's syntax probed
* 2026-10-12 with %TEMP%/probe_msg8.do, StataMP 18): a covariance fixed at 0
* and a loading fixed at 0.5, each in group 2 only
ttg_refuse "2j covariance fixed in group 2" 198 "The covariance of G and F is fixed at 0 in group 2 (g = 1) but free in group 1 (g = 0)." ///
    `"_dwlsmv_translate , metric(unitvariance) group(g) gvals(0 1) paths((F -> y1 y2 y3) (G -> y4 y5 y6)) semopts(cov(1: F*G@0))"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "2k loading fixed in group 2" 198 "Path F -> y2 is fixed at .5 in group 2 (g = 1) but free in group 1 (g = 0)." ///
    `"`tr'((F -> y1 y2 y3) (1: F -> y2@0.5))"'
local nfail = `nfail' + `tt_bad'

display as text "3. dwlsmv's surface refusals"
local m (F -> y1 y2 y3)
ttg_refuse "3a partial()" 198 "partial invariance is not available in this version." ///
    `"dwlsmv `m', group(g) invariance(scalar) metric(unitvariance) partial(y1:F)"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "3b ginvariant()" 198 "option ginvariant() is not available with dwlsmv" ///
    `"dwlsmv `m', group(g) ginvariant(mcoef)"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "3c no group()" 198 "invariance() needs group(); it compares groups. No model was fitted." ///
    `"dwlsmv `m', invariance(scalar) metric(unitvariance)"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "3d no invariance()" 198 "group() without invariance() is not available in this version." ///
    `"dwlsmv `m', group(g)"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "3e unknown level" 198 "invariance() must be configural, thresholds, metric or scalar." ///
    `"dwlsmv `m', group(g) invariance(loadings) metric(unitvariance)"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "3f metric absent" 198 "invariance() requires metric(unitvariance) in this version; specify metric(unitvariance)." ///
    `"dwlsmv `m', group(g) invariance(metric)"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "3g metric(marker)" 198 "the marker metric is not yet available with them" ///
    `"dwlsmv `m', group(g) invariance(metric) metric(marker)"'
local nfail = `nfail' + `tt_bad'
quietly generate str1 gs = string(g)
ttg_refuse "3h string group" 109 "group variable gs is a string variable" ///
    `"dwlsmv `m', group(gs) invariance(metric) metric(unitvariance)"'
local nfail = `nfail' + `tt_bad'
quietly generate double gh = g + 0.5
ttg_refuse "3i non-integer group" 198 "group variable gh must take nonnegative integer values" ///
    `"dwlsmv `m', group(gh) invariance(metric) metric(unitvariance)"'
local nfail = `nfail' + `tt_bad'
quietly generate byte gn = g - 1
ttg_refuse "3j negative group" 198 "group variable gn must take nonnegative integer values" ///
    `"dwlsmv `m', group(gn) invariance(metric) metric(unitvariance)"'
local nfail = `nfail' + `tt_bad'
ttg_refuse "3k one group" 198 "group() has only one value in the estimation sample (g = 0)." ///
    `"dwlsmv `m' if g == 0, group(g) invariance(metric) metric(unitvariance)"'
local nfail = `nfail' + `tt_bad'

* one group left only after the indicators are marked out: every g == 1
* observation has y1 missing, so the second count (after markout) refuses
preserve
quietly replace y1 = . if g == 1
ttg_refuse "3l one group after markout" 198 "group() has only one value in the estimation sample (g = 0)." ///
    `"dwlsmv `m', group(g) invariance(metric) metric(unitvariance)"'
local nfail = `nfail' + `tt_bad'
restore

display as text "4. state restoration"
sort y3 y1
quietly _datasignature
local dsig `r(datasignature)'
local rng `c(rngstate)'
quietly _dwlsmv_translate , paths((F -> y1 y2 y3)) metric(unitvariance) group(g) gvals(0 1)
capture _dwlsmv_translate , paths((F -> y1@1 y2 y3)) metric(unitvariance) group(g) gvals(0 1)
quietly _datasignature
if "`r(datasignature)'" != "`dsig'" | "`c(rngstate)'" != "`rng'" | "`: sortedby'" != "y3 y1" {
    display as error "  4 FAILED: data, rngstate or sort order changed"
    local ++nfail
}
else display as text "  4: data, sort order and rngstate restored"

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_translate_groups: `nfail' failure(s)"
    exit 9
}
display as result "test_translate_groups: the group surface translates to the configural text and refuses what 0.3.0 does not fit, each at its decided rc and message"
exit 0
