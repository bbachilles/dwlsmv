*! Drafted 2026-11-20 by a read-only agent; executed under StataMP 18 on 2026-11-20, all checks pass after
*! three test fixes: the external declarations of YI YH9 YH12 YE in tns_mata(); Case H's K is 4 (codes
*! 0..3), not 3; tns_count() and tns_maxline() consumed the caller's string (Mata passes by reference) and
*! now work on a copy. Check 4's premise holds: Case I under metric(unitvariance) ends e(converged) 0.
*!
*! Release B, B3.8 (spec Module 12 Release B, "Notes (B3.8)"; SO2.3, SO3.3): the two notes of
*! dwlsmv_psi_check.mata and dwlsmv_m11_fit.mata.
*!   (a) an improper fit of G over exactly three factors: "G has exactly three first-order factors (...), so
*!       this model fits exactly as well as one in which those three factors are simply correlated. ..."
*!       after the psi lines.
*!   (b) a fit that did not converge under metric(unitvariance) with a second-order factor: "This model has
*!       a second-order factor. With metric(unitvariance) ..."
*!
*! Builds build/ldwlsmv.mlib first and runs against the library and src/ado. Exit 0 pass, 9 fail.
*!
*! Checks (each names the mutant it kills; N1-N5 are listed at the end):
*!   1  Case I (improper), G over three, metric(marker), at the command: the note appears once, naming
*!      (F1, F2, F3), after the line "F1: disturbance variance ... is negative".            kills N2 (inverse)
*!   1b the same at the Mata level (dwlsmv_psi_note on the converged improper fit): the exact text, every
*!      line at most 78 characters, char(10) endings, inside the note before "Point estimates are reported."
*!   2  Case H, G over three (y1-y9, proper), at the command: no note (a).                     kills N2
*!   2b Case H G over three at the Mata level: no class 2 anywhere -> "". With psi_F1 hand-set to -0.5
*!      (class 2 against Case H's own band) the note appears: the control for 3 and 2.         kills N2
*!   3  Case H, G over four (y1-y12, proper), at the command: no note (a). 3b: no improper G-over-four
*!      fixture exists (READ: Case I is the only improper one, fixtures/data/case_I.csv, three factors), so at
*!      the Mata level psi_F1 of Case H's G-over-four fit is hand-set to -0.5 (class 2): the note is "".
*!      Control: the same edit on G over three gives the note (2b).                           kills N1
*!   4  Case I, G over three, metric(unitvariance), at the command: what it ends in is NOT KNOWN before the
*!      first run. If e(converged) is 0: note (b) appears once, after the admissibility note's lines, with
*!      its exact text. If 1: note (b) is absent. The state is printed. A local switch (strictpremise) makes
*!      "it converged" a failure, as test_unconv_post.do treats a premise that fails.        kills N5
*!   4b Mata level: Case I's converged fit with converged hand-set to 0, metric "unitvariance": the exact
*!      text, lines at most 78 characters.                                                  (control)
*!   5  A first-order fit that does not converge: (command) the marker-trap model of test_unconv_post.do under
*!      metric(marker): e(converged) 0 (premise) and no note (b). (Mata) Case E's first-order fit with
*!      converged set to 0 and metric "unitvariance": "". No known first-order fit fails to converge under
*!      unitvariance at the command, so the second-order condition is killed at the Mata level.  kills N4
*!   6  Case I, G over three, metric(marker), converged (command): no note (b). (Mata) the same fit,
*!      metric "marker", converged set to 0: "".                                            kills N3
*!      (Mata) Case I's fit as it is (converged 1), metric "unitvariance": "".               kills N5
*!
*! Mutants: N1 the 'exactly three' test loosened to 'three or more' (nb >= 3); N2 the negative-class
*! condition dropped (neg = 1); N3 note (b)'s metric condition dropped; N4 its second-order condition
*! dropped; N5 its non-convergence condition dropped.
*! Run 2026-11-20, each in a scratch worktree, all five killed (exit 9): N1 by 3b, N2 by 2b, N3 by 6
*! (marker), N4 by 5, N5 by 6 (converged).
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

* 1 makes "it converged under unitvariance" a failure (check 4); 0 prints it and goes on
local strictpremise 0

foreach c in E H I {
    import delimited using "fixtures/data/case_`c'.csv", clear varnames(1) ///
        delimiters(",") stringcols(_all)
    destring _all, replace
    mata: Y`c' = st_data(., .)
}
mata: YH9  = YH[., 1..9]
mata: YH12 = YH[., 1..12]

mata:

string scalar tns_join(string colvector L)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= rows(L); i++) {
        s = s + L[i] + char(10)
    }
    return(s)
}

// Metric V, nk first-order factors of three items each under G (the same builder as test_psi_check.do's
// tp_canon_v): F1..Fnk, G; marker on each factor's first item; 'metric unitvariance G'.
string scalar tns_canon_v(real scalar nk)
{
    string colvector L
    string scalar fl, il
    real scalar k, j, y

    fl = ""
    for (k = 1; k <= nk; k++) {
        fl = fl + " F" + strofreal(k)
    }
    il = ""
    for (y = 1; y <= 3 * nk; y++) {
        il = il + " y" + strofreal(y)
    }
    L = ("factor" + fl + " G") \ ("ind" + il)
    y = 0
    for (k = 1; k <= nk; k++) {
        for (j = 1; j <= 3; j++) {
            y++
            if (j == 1) {
                L = L \ ("metric marker F" + strofreal(k) + " y" + strofreal(y))
            }
            else {
                L = L \ ("lambda F" + strofreal(k) + " y" + strofreal(y) + " free")
            }
        }
        L = L \ ("psi F" + strofreal(k) + " F" + strofreal(k) + " free")
        L = L \ ("beta F" + strofreal(k) + " G free")
    }
    L = L \ "metric unitvariance G"
    return(tns_join(L))
}

// Case E's first-order model (F2 on F1), marker metric: no second-order factor.
string scalar tns_canon_e()
{
    return(tns_join(("factor F1 F2" \ "ind y1 y2 y3 y4 y5 y6" \
        "metric marker F1 y1" \ "lambda F1 y2 free" \ "lambda F1 y3 free" \
        "metric marker F2 y4" \ "lambda F2 y5 free" \ "lambda F2 y6 free" \
        "psi F1 F1 free" \ "psi F2 F2 free" \ "beta F2 F1 free")))
}

real scalar tns_psik(struct dwlsmv_parmap scalar par, real scalar f)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.row[i] == f & par.col[i] == f & par.free[i]) {
            return(par.fidx[i])
        }
    }
    return(0)
}

real scalar tns_chk(string scalar what, real scalar ok)
{
    if (ok) {
        printf("{txt}  ok    %s\n", what)
        return(0)
    }
    printf("{err}  FAIL  %s\n", what)
    return(1)
}

struct dwlsmv_fit scalar tns_fit(real matrix Y, real colvector K, string scalar canon)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_parmap scalar par

    sts    = dwlsmv_stage1(1, 1)
    sts[1] = dwls_sample_stats(Y, K)
    par    = dwls_identify(canon, sts[1].map)
    return(dwls_fit(sts, par))
}

// how many times needle occurs in s
real scalar tns_count(string scalar s0, string scalar needle)
{
    string scalar s
    real scalar n, p

    s = s0    // Mata passes by reference: consume a copy, never the caller's string

    n = 0
    while ((p = strpos(s, needle)) > 0) {
        n++
        s = substr(s, p + strlen(needle), .)
    }
    return(n)
}

// the longest line of a char(10)-terminated text, and 1 if it ends in char(10)
real scalar tns_maxline(string scalar s0)
{
    string scalar s
    real scalar m, p

    s = s0    // Mata passes by reference: consume a copy, never the caller's string

    m = 0
    while ((p = strpos(s, char(10))) > 0) {
        m = max((m, p - 1))
        s = substr(s, p + 1, .)
    }
    return(m)
}

// the exact sentences, single-spaced
string scalar tns_text_a()
{
    return("G has exactly three first-order factors (F1, F2, F3), so this model fits exactly as well as one in " +
           "which those three factors are simply correlated. The negative disturbance variance means their " +
           "correlations cannot be represented by a proper second-order factor.")
}

string scalar tns_text_b()
{
    return("This model has a second-order factor. With metric(unitvariance) its first-order factors cannot " +
           "have a negative disturbance variance; with metric(marker) they can, so a fit that needs one may " +
           "converge there.")
}

real scalar tns_mata()
{
    struct dwlsmv_fit scalar fI, fH3, fH4, fE, fn
    struct dwlsmv_parmap scalar pI, pH3, pH4
    real colvector th, bandH3, bandH4
    real matrix C
    external real matrix YI, YH9, YH12, YE
    string colvector facI, facH3, facH4, facE
    string scalar note, a, b, flat
    real scalar nf, k3, k4

    nf    = 0
    facI  = ("F1" \ "F2" \ "F3" \ "G")
    facH3 = facI
    facH4 = ("F1" \ "F2" \ "F3" \ "F4" \ "G")
    facE  = ("F1" \ "F2")
    a     = tns_text_a()
    b     = tns_text_b()

    // ---- 1b: Case I, converged improper fit, the note as dwlsmv_psi_note returns it ---------------------
    fI   = tns_fit(YI, J(9, 1, 4), tns_canon_v(3))
    pI   = fI.par
    nf   = nf + tns_chk("1b precondition: Case I converged", fI.converged == 1)
    note = dwlsmv_psi_note(fI, facI, "")
    printf("{txt}%s", note)
    flat = stritrim(subinstr(note, char(10), " "))
    nf   = nf + tns_chk("1b the text of note (a), once, names F1, F2, F3", tns_count(flat, a) == 1)
    nf   = nf + tns_chk("1b after the psi line and before the closing paragraph",
                        strpos(flat, a) > strpos(flat, "F1: disturbance variance") &
                        strpos(flat, a) < strpos(flat, "Point estimates are reported."))
    // line length and ending of note (a) itself (the psi line before it is not B3.8's and runs longer)
    C    = dwlsmv_psi_check(fI.theta, pI, dwlsmv_psi_band(fI.theta, pI, fI.stage1, fI.Delta))
    note = dwlsmv_psi_so3note(C, pI, facI)
    nf   = nf + tns_chk("1b note (a): no line longer than 78 characters", tns_maxline(note) <= 78 & tns_maxline(note) > 0)
    nf   = nf + tns_chk("1b note (a) ends in char(10)", substr(note, strlen(note), 1) == char(10))
    // the same note merged into an existing admissibility note keeps its place
    note = dwlsmv_psi_note(fI, facI, "head" + char(10) + "  Point estimates are reported." + char(10))
    flat = stritrim(subinstr(note, char(10), " "))
    nf   = nf + tns_chk("1b merged into a non-empty note: before 'Point estimates are reported.'",
                        strpos(flat, a) > 0 & strpos(flat, a) < strpos(flat, "Point estimates are reported."))

    // ---- 2b, 3b: Case H's proper fits, and the same fits with psi_F1 hand-set to -0.5 ------------------
    fH3    = tns_fit(YH9, J(9, 1, 4), tns_canon_v(3))
    pH3    = fH3.par
    k3     = tns_psik(pH3, 1)
    bandH3 = dwlsmv_psi_band(fH3.theta, pH3, fH3.stage1, fH3.Delta)
    nf     = nf + tns_chk("2b precondition: Case H G over three converged", fH3.converged == 1 & k3 > 0)
    C      = dwlsmv_psi_check(fH3.theta, pH3, bandH3)
    nf     = nf + tns_chk("2b proper fit: no class 2 and note (a) is empty",
                          !any(C[., 4] :== 2) & dwlsmv_psi_so3note(C, pH3, facH3) == "")
    th     = fH3.theta
    th[k3] = -0.5
    C      = dwlsmv_psi_check(th, pH3, bandH3)
    nf     = nf + tns_chk("2b control: psi_F1 = -0.5 is class 2 and, G over three, note (a) is there",
                          any(C[., 4] :== 2) & tns_count(stritrim(subinstr(dwlsmv_psi_so3note(C, pH3, facH3), char(10), " ")), a) == 1)

    fH4    = tns_fit(YH12, J(12, 1, 4), tns_canon_v(4))
    pH4    = fH4.par
    k4     = tns_psik(pH4, 1)
    bandH4 = dwlsmv_psi_band(fH4.theta, pH4, fH4.stage1, fH4.Delta)
    nf     = nf + tns_chk("3b precondition: Case H G over four converged", fH4.converged == 1 & k4 > 0)
    th     = fH4.theta
    th[k4] = -0.5
    C      = dwlsmv_psi_check(th, pH4, bandH4)
    nf     = nf + tns_chk("3b G over four with psi_F1 = -0.5 (class 2): note (a) is empty",
                          any(C[., 4] :== 2) & dwlsmv_psi_so3note(C, pH4, facH4) == "")

    // ---- 4b: note (b) as the function returns it ----------------------------------------------------------
    fn           = fI
    fn.converged = 0
    note         = dwlsmv_psi_unitnote(fn, "unitvariance")
    printf("{txt}%s", note)
    nf           = nf + tns_chk("4b not converged, unitvariance, second-order: the exact text of note (b)",
                                strtrim(stritrim(subinstr(note, char(10), " "))) == b)
    nf           = nf + tns_chk("4b no output line longer than 78 characters, ends in char(10)",
                                tns_maxline(note) <= 78 & tns_maxline(note) > 0 & substr(note, strlen(note), 1) == char(10))

    // ---- 5, 6: the conditions of note (b), one at a time ------------------------------------------------------
    fE = tns_fit(YE, J(6, 1, 3), tns_canon_e())
    nf = nf + tns_chk("5 precondition: Case E (first-order) converged", fE.converged == 1)
    fE.converged = 0
    nf = nf + tns_chk("5 first-order fit, not converged, unitvariance: no note (b)", dwlsmv_psi_unitnote(fE, "unitvariance") == "")
    nf = nf + tns_chk("6 second-order, not converged, metric marker: no note (b)", dwlsmv_psi_unitnote(fn, "marker") == "")
    nf = nf + tns_chk("6 second-order, CONVERGED, metric unitvariance: no note (b)", dwlsmv_psi_unitnote(fI, "unitvariance") == "")

    return(nf)
}

end

mata: st_numscalar("tns_nmata", tns_mata())
local nfail = scalar(tns_nmata)

* ---- command level -------------------------------------------------------------------------------------------
* run a command, keep its log's text (blanks collapsed) in the Mata global tns_txt, and its rc in r(rc)
capture program drop tns_run
program tns_run, rclass
    tempfile lg
    quietly log using "`lg'", text replace name(tnslog)
    capture noisily `0'
    local rc = _rc
    quietly log close tnslog
    mata: tns_txt = stritrim(invtokens(cat(st_local("lg"))', " "))
    return scalar rc = `rc'
end

* count of a phrase / position of a phrase in tns_txt
capture program drop tns_cnt
program tns_cnt, rclass
    args phrase
    mata: st_numscalar("tns_n", tns_count(tns_txt, st_local("phrase")))
    return scalar n = scalar(tns_n)
end

capture program drop tns_pos
program tns_pos, rclass
    args phrase
    mata: st_numscalar("tns_p", strpos(tns_txt, st_local("phrase")))
    return scalar p = scalar(tns_p)
end

local ta "G has exactly three first-order factors (F1, F2, F3), so this model fits exactly as well as one in which those three factors are simply correlated."
local tb "This model has a second-order factor. With metric(unitvariance) its first-order factors cannot have a negative disturbance variance"
local ta_head "has exactly three first-order factors"
local imodel "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (G -> F1 F2 F3)"
local h3model "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (G -> F1 F2 F3)"
local h4model "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (F4 -> y10 y11 y12) (G -> F1 F2 F3 F4)"

* 1: Case I, marker
import delimited using "fixtures/data/case_I.csv", clear varnames(1) delimiters(",")
tns_run dwlsmv `imodel', metric(marker)
local rc1 = r(rc)
tns_cnt "`ta_head'"
local n1 = r(n)
tns_cnt "`ta'"
local n1x = r(n)
tns_pos "F1: disturbance variance"
local p1psi = r(p)
tns_pos "`ta_head'"
local p1a = r(p)
display as text "  1: rc `rc1', e(converged) " e(converged) ", note (a) x `n1' (full sentence x `n1x'), psi line at `p1psi', note at `p1a'"
if `rc1' | e(converged) != 1 | `n1' != 1 | `n1x' != 1 | `p1psi' == 0 | `p1a' <= `p1psi' {
    display as error "  FAIL 1: Case I, G over three, metric(marker): note (a) once, naming F1 F2 F3, after the psi line"
    local ++nfail
}
tns_cnt "This model has a second-order factor"
if r(n) != 0 {
    display as error "  FAIL 6: Case I, metric(marker), converged: note (b) must be absent"
    local ++nfail
}

* 4: Case I, unitvariance
tns_run dwlsmv `imodel', metric(unitvariance)
local rc4 = r(rc)
local conv4 = e(converged)
tns_cnt "`tb'"
local n4 = r(n)
tns_cnt "This model has a second-order factor"
local n4h = r(n)
display as text "  4: rc `rc4', e(converged) `conv4', e(recovery) `e(recovery)', note (b) x `n4' (head x `n4h')"
if `rc4' {
    display as error "  FAIL 4: dwlsmv returned rc `rc4'"
    local ++nfail
}
else if `conv4' == 0 {
    if `n4' != 1 | `n4h' != 1 {
        display as error "  FAIL 4: not converged under unitvariance, second-order: note (b) must appear once"
        local ++nfail
    }
}
else {
    display as text "  4: PREMISE NOT MET: Case I under metric(unitvariance) converged; the placement of note (b) in the command's output is NOT exercised"
    if `n4h' != 0 {
        display as error "  FAIL 4: converged, yet note (b) appears"
        local ++nfail
    }
    if `strictpremise' {
        display as error "  FAIL 4: the premise (Case I unitvariance does not converge) does not hold"
        local ++nfail
    }
}
tns_cnt "`ta_head'"
if r(n) != 0 & `conv4' == 0 {
    display as text "  4: note (a) appears on a non-converged fit (the psi lines are not printed for it; unexpected)"
    local ++nfail
}

* 2: Case H, G over three (y1-y9), proper
import delimited using "fixtures/data/case_H.csv", clear varnames(1) delimiters(",")
tns_run dwlsmv `h3model', metric(marker)
local rc2 = r(rc)
tns_cnt "`ta_head'"
local n2 = r(n)
tns_cnt "is negative"
local n2neg = r(n)
display as text "  2: rc `rc2', e(converged) " e(converged) ", note (a) x `n2', 'is negative' x `n2neg'"
if `rc2' | e(converged) != 1 | `n2neg' != 0 {
    display as error "  FAIL 2 precondition: Case H, G over three, must fit and be proper"
    local ++nfail
}
if `n2' != 0 {
    display as error "  FAIL 2: a proper fit must carry no note (a)"
    local ++nfail
}

* 3: Case H, G over four (y1-y12), proper
tns_run dwlsmv `h4model', metric(marker)
local rc3 = r(rc)
tns_cnt "`ta_head'"
local n3 = r(n)
tns_cnt "is negative"
local n3neg = r(n)
display as text "  3: rc `rc3', e(converged) " e(converged) ", note (a) x `n3', 'is negative' x `n3neg'"
if `rc3' | e(converged) != 1 | `n3neg' != 0 {
    display as error "  FAIL 3 precondition: Case H, G over four, must fit and be proper"
    local ++nfail
}
if `n3' != 0 {
    display as error "  FAIL 3: G over four must carry no note (a)"
    local ++nfail
}

* 5 (command): a first-order fit that does not converge (marker trap, metric marker): no note (b)
import delimited using "fixtures/data/marker_trap.csv", clear varnames(1) delimiters(",") stringcols(_all)
quietly destring q* att*, replace
local mT "(MathAb -> q1@1 q2@1 q3-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb)"
tns_run dwlsmv `mT'
local rc5 = r(rc)
local conv5 = e(converged)
tns_cnt "This model has a second-order factor"
local n5 = r(n)
display as text "  5: rc `rc5', e(converged) `conv5', note (b) head x `n5'"
if `rc5' | `conv5' != 0 {
    display as error "  FAIL 5 precondition: the marker-trap model must run and not converge (test_unconv_post.do premise)"
    local ++nfail
}
if `n5' != 0 {
    display as error "  FAIL 5: a first-order fit that did not converge must carry no note (b)"
    local ++nfail
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_notes_so: `nfail' failure(s)"
    exit 9
}
display as result "test_notes_so: note (a) only for an improper G over exactly three, note (b) only for a non-converged unit-variance second-order fit"
exit 0
