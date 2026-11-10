*! The EFA surface (0.3.0, Unit 6): dwlsmv (F1 F2 -> y1-y6), efa(F1 F2).
*!
*! Executed under StataMP 18 on 2026-10-29 and it passes. Worst ratios to tier: Lambda/Phi
*! 7.1e-3, T3 0.13 (Case F); fit indices <= 2.8e-4. What broke first: a 1 x 0
*! selectindex() result passed a rows() == 1 test in dwlsmv_m11_efa_rotate()
*! (rc 3301); a fixed covariance went unrefused because sem leaves no psi row
*! for it; df = 0 at m = 3 on Case E is a degenerate over-extraction (check 8).
*!
*! End to end through build/ldwlsmv.mlib and src/ado, as test_dwlsmv.do.
*! Reference: fixtures/expected/lavaan_efa_geomin.csv (99h_lavaan_efa_geomin.R,
*! lavaan 0.7.2): rotated loadings, factor correlations, T3, df and the scaled
*! fit indices on Cases E (m = 2), F (m = 3), G (m = 2, binary).
*!
*! Tolerances, fixed before the first run (2026-10-29):
*!   e(Lambda), e(Phi)          theta-hat tier, 1e-5 abs + 1e-4 rel (as test_geomin.do)
*!   e(T3), e(rmsea_s), e(cfi_s), e(tli_s)   1e-6 abs + 1e-5 rel (spec 7.2's
*!     T3 row; the indices are functions of T3, df and the c1 baseline)
*!   e(df_T3), counts, e(rot_eps)  exact
*! Limitation, recorded with the fixture: every case fits well (T3 < df), so
*! RMSEA and CFI sit at their clamps (0, 1); TLI is the unclamped check.
*!
*! Checks:
*!   1  per case: rc 0; e(Lambda), e(Phi) vs lavaan; T3, df, fit indices vs
*!      lavaan; e(rot_eps) the Mplus rule; e(rot_starts) 30
*!   2  Case E: e(b) holds the rotated loadings, the correlation and the
*!      thresholds, matching e(Lambda) and e(Phi); no e(V) (rc 111); e(efa)
*!   3  epsilon(.001) and rstarts(40) reach e(rot_eps) and e(rot_starts)
*!   4  df = 0 (Case E, one factor on y1-y3) is fitted and the note is printed
*!   5  replay prints the rotated table; replay with standardized is rc 198
*!   6  refusals rc 198, nothing fitted: group(), invariance(), metric(),
*!      standardized, a name not latent, a latent not in efa(), a partial
*!      block, an extra path, a fixed covariance, df < 0, epsilon() without
*!      efa(), epsilon(-1), rstarts(2.5)
*!   7  difftest refuses an EFA fit, rc 3498
*!   8  over-extraction (3 factors on Case E): rc 3498, guidance printed,
*!      nothing posted
*!   9  efa(F2 F1): rotated column 1 is named F2 in e(Lambda) and e(b), as
*!      the help says (review finding, 2026-10-29)
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

capture program drop tes_run
program tes_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(teslog)
    capture noisily `0'
    local rc = _rc
    quietly log close teslog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end

import delimited using "fixtures/expected/lavaan_efa_geomin.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: XC = st_sdata(., "case")
mata: XQ = st_sdata(., "quantity")
mata: XL = st_sdata(., "label")
mata: XV = strtoreal(st_sdata(., "value"))

mata:
real scalar tes_get(string scalar c, string scalar q, string scalar l)
{
    external string colvector XC, XQ, XL
    external real colvector XV
    real colvector k

    k = selectindex((XC :== c) :& (XQ :== q) :& (XL :== l))
    return(rows(k) == 1 ? XV[k] : .)
}

// worst ratio of e(Lambda), e(Phi) to the theta-hat tier against case c
real scalar tes_rot(string scalar c)
{
    real matrix L, P
    real scalar j, k, r, x

    L = st_matrix("e(Lambda)")
    P = st_matrix("e(Phi)")
    r = 0
    for (k = 1; k <= cols(L); k++) {
        for (j = 1; j <= rows(L); j++) {
            x = tes_get(c, "lambda", "f" + strofreal(k) + "=~y" + strofreal(j))
            r = max((r, abs(L[j, k] - x) / (1e-5 + 1e-4 * abs(x))))
        }
        for (j = 1; j < k; j++) {
            x = tes_get(c, "phi", "f" + strofreal(k) + "~~f" + strofreal(j))
            r = max((r, abs(P[k, j] - x) / (1e-5 + 1e-4 * abs(x))))
        }
    }
    return(r)
}

real scalar tes_stat(string scalar c, string scalar e, string scalar q, string scalar l)
{
    real scalar x

    x = tes_get(c, q, l)
    return(abs(st_numscalar(e) - x) / (1e-6 + 1e-5 * abs(x)))
}
end

* 1  per case
local cases  E F G
local files  case_E case_F case_G
local models `""(F1 F2 -> y1-y6)" "(F1 F2 F3 -> y1-y9)" "(F1 F2 -> y1-y6)""'
local efas   `""F1 F2" "F1 F2 F3" "F1 F2""'
local epss   1e-4 1e-3 1e-4
forvalues i = 1/3 {
    local c   : word `i' of `cases'
    local f   : word `i' of `files'
    local mdl : word `i' of `models'
    local ef  : word `i' of `efas'
    local ep  : word `i' of `epss'
    import delimited using "fixtures/data/`f'.csv", clear varnames(1) delimiters(",")
    capture noisily dwlsmv `mdl', efa(`ef')
    if _rc {
        display as error "test_efa_surface 1 (case `c'): rc " _rc
        local ++nfail
        continue
    }
    mata: st_local("r", strofreal(tes_rot("`c'")))
    display as text "  case `c': rotated Lambda/Phi ratio " %8.2e `r'
    if !(`r' <= 1) local ++nfail
    foreach pair in "e(T3) test T3" "e(rmsea_s) fit rmsea.scaled" "e(cfi_s) fit cfi.scaled" "e(tli_s) fit tli.scaled" {
        local e : word 1 of `pair'
        local q : word 2 of `pair'
        local l : word 3 of `pair'
        mata: st_local("r", strofreal(tes_stat("`c'", "`e'", "`q'", "`l'")))
        display as text "  case `c': `e' ratio " %8.2e `r'
        if !(`r' <= 1) local ++nfail
    }
    mata: st_local("df", strofreal(tes_get("`c'", "test", "df")))
    if e(df_T3) != `df' | e(rot_eps) != `ep' | e(rot_starts) != 30 {
        display as error "test_efa_surface 1 (case `c'): df " e(df_T3) " vs `df', eps " e(rot_eps) ", starts " e(rot_starts)
        local ++nfail
    }

    * 2  e(b) and e(V), Case E
    if "`c'" == "E" {
        tempname b L P
        matrix `b' = e(b)
        matrix `L' = e(Lambda)
        matrix `P' = e(Phi)
        local nb = colsof(`b')
        local ok = (`nb' == 12 + 1 + 12)
        local ok = `ok' & reldif(`b'[1, colnumb(`b', "y3:F2")], `L'[3, 2]) == 0
        local ok = `ok' & reldif(`b'[1, colnumb(`b', "/:cov(F1,F2)")], `P'[2, 1]) == 0
        local ok = `ok' & "`e(efa)'" == "F1 F2"
        capture confirm matrix e(V)
        local vrc = _rc
        display as text "  case E: e(b) has `nb' columns; e(V) confirm rc `vrc'"
        if !`ok' | `vrc' != 111 {
            display as error "test_efa_surface 2: e(b)/e(V)/e(efa) not as specified"
            local ++nfail
        }
        estimates store efaE

        * 5  replay
        capture noisily dwlsmv
        if _rc {
            display as error "test_efa_surface 5: replay rc " _rc
            local ++nfail
        }
        capture dwlsmv, standardized
        if _rc != 198 {
            display as error "test_efa_surface 5: replay with standardized rc " _rc ", not 198"
            local ++nfail
        }

        * 3  overrides
        capture noisily dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) epsilon(.001) rstarts(40)
        if _rc | e(rot_eps) != .001 | e(rot_starts) != 40 {
            display as error "test_efa_surface 3: rc " _rc ", eps " e(rot_eps) ", starts " e(rot_starts)
            local ++nfail
        }
        estimates store efaE2

        * 4  df = 0: one factor on three items (m = 1, nothing to rotate).
        * Not m = 3 on Case E's six items, also df = 0: Case E has two true
        * factors, and lavaan's own three-factor fit there is degenerate
        * (optimizer fails, a loading near -252), so no rotation can converge
        * on it; that over-extraction path is a separate question (Unit 6).
        tes_run "With df = 0" dwlsmv (F1 -> y1-y3), efa(F1)
        display as text "  df = 0: rc " r(rc) ", note found " r(found) ", df " e(df_T3)
        if r(rc) != 0 | r(found) != 1 | e(df_T3) != 0 {
            display as error "test_efa_surface 4: df = 0 case not as specified"
            local ++nfail
        }

        * 9  labels follow efa()'s order (help: "The names in efa() are
        * attached to the factors in that final order"; review finding,
        * 2026-10-29): efa(F2 F1) names rotated column 1 F2
        estimates restore efaE
        tempname L1
        matrix `L1' = e(Lambda)
        capture noisily dwlsmv (F1 F2 -> y1-y6), efa(F2 F1)
        local r9 = _rc
        tempname L2 b2
        matrix `L2' = e(Lambda)
        matrix `b2' = e(b)
        local cn : colnames `L2'
        local ok9 = (`r9' == 0) & ("`cn'" == "F2 F1") & ("`e(efa)'" == "F2 F1")
        local ok9 = `ok9' & reldif(`L2'[1, 1], `L1'[1, 1]) == 0
        local ok9 = `ok9' & reldif(`b2'[1, colnumb(`b2', "y1:F2")], `L1'[1, 1]) == 0
        display as text "  efa(F2 F1): rc `r9', Lambda columns `cn'"
        if !`ok9' {
            display as error "test_efa_surface 9: efa() order labelling not as specified"
            local ++nfail
        }
        estimates restore efaE

        * 8  over-extraction: three factors on Case E's six items, whose
        * unrotated fit is degenerate; no start converges, rc 3498 with the
        * guidance (the user's choice, 2026-10-29), and nothing is posted
        ereturn clear
        tes_run "more factors than the data support" dwlsmv (F1 F2 F3 -> y1-y6), efa(F1 F2 F3)
        local r8 = r(rc)
        local f8 = r(found)
        display as text "  over-extraction: rc `r8', guidance found `f8', e(cmd) '`e(cmd)''"
        if `r8' != 3498 | `f8' != 1 | "`e(cmd)'" != "" {
            display as error "test_efa_surface 8: over-extraction not as specified"
            local ++nfail
        }
        estimates restore efaE

        * 7  difftest
        * the EFA refusal itself, not a later check that also fails at 3498
        * (mutation round, 2026-10-29)
        tes_run "is an EFA fit" difftest efaE efaE2
        display as text "  difftest on EFA fits: rc " r(rc) ", EFA refusal found " r(found)
        if r(rc) != 3498 | r(found) != 1 {
            display as error "test_efa_surface 7: difftest rc " r(rc) ", not 3498 with the EFA refusal"
            local ++nfail
        }

        * 6  refusals
        generate g = mod(_n, 2)
        local k 0
        foreach cmd in ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) group(g)"                          ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) invariance(configural)"            ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) metric(unitvariance)"              ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) standardized"                      ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F3)"                                   ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1)"                                      ///
            "dwlsmv (F1 -> y1-y6) (F2 -> y2-y6), efa(F1 F2)"                        ///
            "dwlsmv (F1 F2 -> y1-y6) (F2 <- F1), efa(F1 F2)"                        ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) cov(F1*F2@0)"                      ///
            "dwlsmv (F1 F2 F3 -> y1-y5), efa(F1 F2 F3)"                             ///
            "dwlsmv (F1 F2 -> y1-y6), epsilon(.001)"                                ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) epsilon(-1)"                       ///
            "dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) rstarts(2.5)" {
            local ++k
            capture `cmd'
            display as text "  refusal `k': rc " _rc
            if _rc != 198 {
                display as error "test_efa_surface 6.`k' (`cmd'): rc " _rc ", not 198"
                local ++nfail
            }
        }
    }
}

* ---- check 9, 2026-11-02: an empty efa() is refused; a factor may be named cut* ----
capture dwlsmv (F1 -> y1-y6), efa()
display as text "  9 empty efa(): rc " _rc
if _rc != 198 {
    display as error "test_efa_surface 9 empty efa(): rc " _rc ", not 198"
    local ++nfail
}
capture noisily dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
tempname bn
matrix `bn' = e(b)
capture noisily dwlsmv (cut1 cut2 -> y1-y6), efa(cut1 cut2) nocapslatent latent(cut1 cut2)
display as text "  9 factors named cut*: rc " _rc
if _rc {
    display as error "test_efa_surface 9 factors named cut*: rc " _rc
    local ++nfail
}
else if mreldif(`bn', e(b)) > 1e-10 {
    display as error "test_efa_surface 9 factors named cut*: e(b) differs from F1 F2's"
    local ++nfail
}

if `nfail' {
    display as error _n "test_efa_surface: `nfail' failure(s)"
    exit 9
}
display as result "test_efa_surface: efa() fits E, F, G and matches lavaan's rotated Lambda, Phi, T3 and fit indices; e(b) without e(V); overrides, df = 0, replay, refusals and difftest hold"
exit 0
