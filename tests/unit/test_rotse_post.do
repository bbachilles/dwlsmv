*! Release A, unit A5: the posting layer of the rotated-EFA covariance, end to end.
*!
*! Drafted 2026-11-17 by a read-only agent; executed and fixed by the main session under StataMP 18
*! at 6 and 4 processors (2026-11-17). Fixed on first run: the joint test used _b[] names (check 3).
*!
*! WHAT THIS CERTIFIES. dwlsmv (F1 F2 -> y1-y6), efa(F1 F2) and its Case F and G analogues post
*! the rotated e(V), through src/ado/dwlsmv.ado, build/ldwlsmv.mlib and
*! src/mata/dwlsmv_m11_efa_se.mata, which agree with the committed tier-2 fixture
*! (lavaan 0.7-3, pinned settings: std_ov FALSE, information.bread "expected", bb at 1e-12;
*! fixtures/expected/lavaan_efa_se.csv, lavaan_efa_se_vcov.csv):
*!
*! DECIDING CRITERION (one, as the project rule on numerical validation asks; fixed before any run):
*!   every standard error of every rotated e(b) column, and every lower-triangle element of the
*!   rotated covariance, is within its tolerance of lavaan's, on Cases E, F and G.
*!     SE:          |a - b| <= 1e-8 + 1e-5 |b|             spec 7.2, "Standard errors" row
*!     covariance:  |a - b| <= 1e-8 + 1e-5 sqrt(V_L,ii V_L,jj)
*!                  spec 7.2.1's element form, as a tolerance in its own right (7.2, A3.9);
*!                  V_L is lavaan's. 7.2.1's attribution rule is NOT applied (not extended to
*!                  rotation, owner 2026-11-17): a miss is reported as a miss.
*!   A comparator, independent of dwlsmv in code, algorithm and representation: lavaan's bordered
*!   information route (R) against production's J ACOV_u J' (Mata). It shares the unrotated fit's
*!   Gamma only in the sense that both estimate the same population quantity from the same data.
*!
*! PRECONDITION, counted (not the deciding criterion): rotated estimates against the same
*! fixture at the theta-hat tier, 1e-5 abs + 1e-4 rel (spec 7.2). If the point misses, the SE
*! comparison is not interpretable; the file then fails here first. Not chosen after measuring.
*!
*! DIAGNOSTICS (printed; counted only where marked):
*!   - e(rank) = n_u (23, 51, 17), columns n* (25, 57, 19)              [counted: A3.14, A3.5]
*!   - e(V) symmetric, with e(b)'s stripe on both axes                   [counted]
*!   - the numerical rank of e(V) by the project's one 1e-8 cut = n_u    [counted: A3.5]
*!   - e(rot_nose) empty, e(rot_cond) finite and positive                [counted: A3.2]
*!   - lincom, test (single, joint) on ordinary combinations of rotated loadings, against
*!     answers computed from e(V) in Stata matrices                      [counted: A3.5]
*!   - replay, estimates store / restore keep e(V)                       [counted]
*!   - m = 1, df = 0 (one factor on three items): e(V) posted, e(rank) = n* = 9  [counted]
*!   - the no-e(V) path (A3.4): see check 6. Which branch Case G without its (y1 = 0, y2 = 1)
*!     cell reaches is NOT known to the drafter; the check accepts either and PRINTS which.
*!   - the no-converged-start path stays rc 3498 with nothing posted     [counted: A3.4 "unchanged"]
*!
*! WHAT IT DOES NOT CERTIFY. Calibration of the SEs (A6); Cases other than E, F, G; m >= 4;
*! tier 3 (Mplus); a combination in the orbit tangent (A3.5: its outcome follows the sign of a
*! 1e-17 residue and differs by machine, so it is never asserted); wording of any message
*! (A7: placeholders are not asserted except by their presence/absence of the matrix branch).
*!
*! DETERMINISM. Run at the default processor count and at set processors 4 (a project rule: wrapper
*! do-file in %TEMP%). Nothing here compares across runs: the lines "TPDET ..." print trace(V),
*! the largest SE and one off-diagonal element to 15 digits, per case and per processor count,
*! for the main session to diff between the two logs. Production uses ordinary Mata products
*! (J ACOV_u J', Delta*), not quadcross(); a difference at 1e-16 relative is possible and
*! harmless to every tolerance here, but a difference above 1e-12 would be a finding.
*!
*! Checks:
*!   1  per case E, F, G: fit, shapes, e(rank), e(rot_nose), e(rot_cond), symmetry, stripes,
*!      numerical rank
*!   2  per case: precondition (estimates) and the deciding criterion (SEs, covariances)
*!   3  per case: lincom and test against answers from e(V); _se[] equals sqrt of the diagonal
*!   4  Case E: replay, replay with level(90), estimates store/restore; the table branch
*!      (no "Rotated loadings" matrices when e(V) is posted)
*!   5  m = 1, df = 0
*!   6  the A3.4 path
*!   7  no converged start: rc 3498, nothing posted
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

display "TPENV stata_version " c(stata_version) " processors " c(processors)

local nfail 0

capture program drop tpp_run
program tpp_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tpplog)
    capture noisily `0'
    local rc = _rc
    quietly log close tpplog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end

* ---- the tier-2 fixture, as text, converted with strtoreal() (the float-import trap) -----------------
import delimited using "fixtures/expected/lavaan_efa_se.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: XC = st_sdata(., "case")
mata: XQ = st_sdata(., "quantity")
mata: XL = st_sdata(., "label")
mata: XV = strtoreal(st_sdata(., "value"))
import delimited using "fixtures/expected/lavaan_efa_se_vcov.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: WC = st_sdata(., "case")
mata: WR = st_sdata(., "row")
mata: WK = st_sdata(., "col")
mata: WV = strtoreal(st_sdata(., "value"))

mata:

// the Stata column of e(b) for lavaan's label: fK=~yJ -> y J : F K, yJ|tK -> /yJ : cutK,
// fA~~fB -> /: cov(FB,FA) or cov(FA,FB) (whichever the stripe holds). Missing if absent.
real scalar tpp_col(string matrix S, string scalar lab)
{
    string scalar l, r, a, b, f
    real colvector k

    if (strpos(lab, "=~") > 0) {
        l = substr(lab, 1, strpos(lab, "=~") - 1)
        r = substr(lab, strpos(lab, "=~") + 2, .)
        f = strupper(substr(l, 1, 1)) + substr(l, 2, .)
        k = selectindex((S[., 1] :== r) :& (S[., 2] :== f))
    }
    else if (strpos(lab, "|") > 0) {
        l = substr(lab, 1, strpos(lab, "|") - 1)
        r = substr(lab, strpos(lab, "|") + 1, .)
        k = selectindex((S[., 1] :== "/" + l) :& (S[., 2] :== "cut" + substr(r, 2, .)))
    }
    else if (strpos(lab, "~~") > 0) {
        l = substr(lab, 1, strpos(lab, "~~") - 1)
        r = substr(lab, strpos(lab, "~~") + 2, .)
        a = strupper(substr(l, 1, 1)) + substr(l, 2, .)
        b = strupper(substr(r, 1, 1)) + substr(r, 2, .)
        k = selectindex((S[., 1] :== "/") :& ((S[., 2] :== "cov(" + a + "," + b + ")") :|
                                              (S[., 2] :== "cov(" + b + "," + a + ")")))
    }
    else {
        return(.)
    }
    return(length(k) == 1 ? k : .)
}

// lavaan's variance of label lab, case cs, from the vcov fixture's diagonal
real scalar tpp_vdiag(string scalar cs, string scalar lab)
{
    external string colvector WC, WR, WK
    external real colvector WV
    real colvector k

    k = selectindex((WC :== cs) :& (WR :== lab) :& (WK :== lab))
    return(length(k) == 1 ? WV[k] : .)
}

// The deciding criterion and the precondition for one case, against the posted e(b), e(V).
// Fills out = (n_est, worst_est, n_se, worst_se, n_cov, worst_cov, n_unmatched); the "worst"
// values are ratios to the tolerance (<= 1 passes). The worst label is printed.
real rowvector tpp_compare(string scalar cs)
{
    external string colvector XC, XQ, XL, WC, WR, WK
    external real colvector XV, WV
    string matrix S
    real rowvector b, out
    real matrix V
    real colvector ix
    real scalar i, c1, c2, x, y, tol, r, vl
    string scalar wl_est, wl_se, wl_cov

    b = st_matrix("e(b)")
    V = st_matrix("e(V)")
    S = st_matrixcolstripe("e(b)")
    out = (0, 0, 0, 0, 0, 0, 0)
    wl_est = ""
    wl_se  = ""
    wl_cov = ""

    ix = selectindex((XC :== cs) :& (XQ :== "est"))
    for (i = 1; i <= length(ix); i++) {
        c1 = tpp_col(S, XL[ix[i]])
        if (c1 >= .) {
            out[7] = out[7] + 1
            continue
        }
        x = b[c1]
        y = XV[ix[i]]
        r = abs(x - y) / (1e-5 + 1e-4 * abs(y))
        out[1] = out[1] + 1
        if (r > out[2]) {
            out[2] = r
            wl_est = XL[ix[i]]
        }
    }
    ix = selectindex((XC :== cs) :& (XQ :== "se"))
    for (i = 1; i <= length(ix); i++) {
        c1 = tpp_col(S, XL[ix[i]])
        if (c1 >= .) {
            out[7] = out[7] + 1
            continue
        }
        x = sqrt(V[c1, c1])
        y = XV[ix[i]]
        r = abs(x - y) / (1e-8 + 1e-5 * abs(y))
        out[3] = out[3] + 1
        if (r > out[4] | r >= .) {
            out[4] = r
            wl_se  = XL[ix[i]]
        }
    }
    ix = selectindex(WC :== cs)
    for (i = 1; i <= length(ix); i++) {
        c1 = tpp_col(S, WR[ix[i]])
        c2 = tpp_col(S, WK[ix[i]])
        if (c1 >= . | c2 >= .) {
            out[7] = out[7] + 1
            continue
        }
        tol = 1e-8 + 1e-5 * sqrt(tpp_vdiag(cs, WR[ix[i]]) * tpp_vdiag(cs, WK[ix[i]]))
        r   = abs(V[c1, c2] - WV[ix[i]]) / tol
        out[5] = out[5] + 1
        if (r > out[6] | r >= .) {
            out[6] = r
            wl_cov = WR[ix[i]] + " , " + WK[ix[i]]
        }
    }
    printf("  case %s: estimates %g, worst ratio %9.2e (%s); SEs %g, worst %9.2e (%s);" +
           " covariances %g, worst %9.2e (%s); unmatched %g\n", cs, out[1], out[2], wl_est,
           out[3], out[4], wl_se, out[5], out[6], wl_cov, out[7])
    return(out)
}

// rank of e(V) by the project's one criterion, and the smallest diagonal element
real rowvector tpp_vshape()
{
    real matrix V, W
    real colvector sv

    V = st_matrix("e(V)")
    sv = J(0, 1, .)
    W  = J(0, 0, .)
    return((dwlsmv_numrank(V, sv, W), min(diagonal(V)), (max(abs(V - V')) == 0)))
}

// the determinism lines: diffed between the 6- and 4-processor logs by the main session
void tpp_det(string scalar cs)
{
    real matrix V

    V = st_matrix("e(V)")
    printf("TPDET %s trace %22.15e maxdiag %22.15e V21 %22.15e Vlast %22.15e\n", cs, trace(V),
           max(diagonal(V)), V[2, 1], V[rows(V), cols(V) - 1])
}
end

* ---- 1, 2, 3  per case -------------------------------------------------------------------
local cases  E F G
local files  case_E case_F case_G
local models `""(F1 F2 -> y1-y6)" "(F1 F2 F3 -> y1-y9)" "(F1 F2 -> y1-y6)""'
local efas   `""F1 F2" "F1 F2 F3" "F1 F2""'
local nus    23 51 17
local nstars 25 57 19
forvalues i = 1/3 {
    local c   : word `i' of `cases'
    local f   : word `i' of `files'
    local mdl : word `i' of `models'
    local ef  : word `i' of `efas'
    local nu  : word `i' of `nus'
    local ns  : word `i' of `nstars'
    import delimited using "fixtures/data/`f'.csv", clear varnames(1) delimiters(",")
    capture noisily dwlsmv `mdl', efa(`ef')
    if _rc {
        display as error "test_rotse_post 1 (case `c'): rc " _rc
        local ++nfail
        continue
    }

    * 1  shapes, rank, withheld flag, condition number, symmetry, stripes, numerical rank
    tempname b V Vt
    capture confirm matrix e(V)
    if _rc {
        display as error "test_rotse_post 1 (case `c'): no e(V) posted (e(rot_nose) = '`e(rot_nose)'')"
        local ++nfail
        continue
    }
    matrix `b' = e(b)
    matrix `V' = e(V)
    matrix `Vt' = `V''
    local ok = (colsof(`b') == `ns') & (rowsof(`V') == `ns') & (colsof(`V') == `ns')
    local ok = `ok' & (e(rank) == `nu') & ("`e(rot_nose)'" == "")
    local ok = `ok' & (e(rot_cond) < .) & (e(rot_cond) > 0)
    local ok = `ok' & ("`: colfullnames `V''" == "`: colfullnames `b''")
    local ok = `ok' & ("`: rowfullnames `V''" == "`: colfullnames `b''")
    local ok = `ok' & (mreldif(`V', `Vt') == 0)
    mata: st_local("shape", invtokens(strofreal(tpp_vshape(), "%21.0g")))
    local vrank : word 1 of `shape'
    local vmin  : word 2 of `shape'
    local ok = `ok' & (`vrank' == `nu') & (`vmin' > 0)
    display as text "  case `c': n* `ns', e(rank) " e(rank) ", numerical rank of e(V) `vrank', cond " %9.3e e(rot_cond) ", min diag `vmin'"
    if !`ok' {
        display as error "test_rotse_post 1 (case `c'): shape, rank, symmetry, stripes or e(rot_cond) not as specified"
        local ++nfail
    }

    * 2  the precondition and the deciding criterion
    mata: TPR = tpp_compare("`c'")
    mata: st_local("res", invtokens(strofreal(TPR, "%21.0g")))
    local k 0
    foreach t of local res {
        local ++k
        local r`k' `t'
    }
    * r1 n_est, r2 worst_est, r3 n_se, r4 worst_se, r5 n_cov, r6 worst_cov, r7 unmatched

    if `r7' != 0 | `r1' != `ns' | `r3' != `ns' | `r5' != `ns' * (`ns' + 1) / 2 {
        * est and se fixtures hold every rotated column; vcov its lower triangle
        display as error "test_rotse_post 2 (case `c'): fixture rows not all matched (est `r1', se `r3', cov `r5', unmatched `r7')"
        local ++nfail
    }
    if !(`r2' <= 1) {
        display as error "test_rotse_post 2 PRECONDITION (case `c'): estimates miss the theta-hat tier, worst ratio `r2'"
        local ++nfail
    }
    if !(`r4' <= 1) {
        display as error "test_rotse_post 2 DECIDING (case `c'): a standard error misses its tolerance, worst ratio `r4'"
        local ++nfail
    }
    if !(`r6' <= 1) {
        display as error "test_rotse_post 2 DECIDING (case `c'): a covariance misses its tolerance, worst ratio `r6'"
        local ++nfail
    }
    mata: tpp_det("`c'")

    * 3  lincom and test on ordinary combinations, against answers from e(V)
    local c1 = colnumb(`b', "y1:F1")
    local c2 = colnumb(`b', "y2:F1")
    local c3 = colnumb(`b', "y3:F2")
    if `c1' >= . | `c2' >= . | `c3' >= . {
        display as error "test_rotse_post 3 (case `c'): y1:F1, y2:F1 or y3:F2 not in e(b)"
        local ++nfail
        continue
    }
    local b1 = `b'[1, `c1']
    local b2 = `b'[1, `c2']
    local b3 = `b'[1, `c3']
    local v11 = `V'[`c1', `c1']
    local v22 = `V'[`c2', `c2']
    local v12 = `V'[`c1', `c2']
    local v13 = `V'[`c1', `c3']
    local v23 = `V'[`c2', `c3']
    local v33 = `V'[`c3', `c3']

    * _se[] is the square root of the diagonal
    local worst3 = reldif(_se[y1:F1], sqrt(`v11'))
    * a difference of two loadings
    quietly lincom _b[y1:F1] - _b[y2:F1]
    local worst3 = `worst3' + reldif(r(estimate), `b1' - `b2') + reldif(r(se), sqrt(`v11' + `v22' - 2 * `v12'))
    quietly test _b[y1:F1] = _b[y2:F1]
    local worst3 = `worst3' + reldif(r(chi2), (`b1' - `b2')^2 / (`v11' + `v22' - 2 * `v12'))
    local df1 = r(df)
    * a combination across factors with unequal weights
    quietly lincom 2 * _b[y1:F1] - _b[y2:F1] + 3 * _b[y3:F2]
    local ev3 = 4 * `v11' + `v22' + 9 * `v33' - 4 * `v12' + 12 * `v13' - 6 * `v23'
    local worst3 = `worst3' + reldif(r(estimate), 2 * `b1' - `b2' + 3 * `b3') + reldif(r(se), sqrt(`ev3'))
    * a point hypothesis
    quietly test _b[y1:F1] = 0.5
    local worst3 = `worst3' + reldif(r(chi2), (`b1' - 0.5)^2 / `v11')
    * a joint test of two loadings = 0
    quietly test [y1]F1 [y2]F1
    local df2 = r(df)
    tempname bs Vs X
    matrix `bs' = (`b1' \ `b2')
    matrix `Vs' = (`v11', `v12' \ `v12', `v22')
    matrix `X' = `bs'' * invsym(`Vs') * `bs'
    local worst3 = `worst3' + reldif(r(chi2), `X'[1, 1])
    display as text "  case `c': lincom/test worst relative difference " %9.2e `worst3' ", df " `df1' " and " `df2'
    if !(`worst3' < 1e-8) | `df1' != 1 | `df2' != 2 {
        display as error "test_rotse_post 3 (case `c'): lincom or test differs from the answer computed from e(V)"
        local ++nfail
    }

    if "`c'" == "E" {
        * 4  replay, level(), estimates store / restore, and the table branch
        estimates store tppE
        * A3.6 revised (owner 2026-11-17): the default shows the matrices also when e(V) is posted
        * (was: a coefficient table, matrices found 0); se/ci are tested in test_efa_display.do
        tpp_run "is posted in e(V)" dwlsmv
        local rr = r(rc)
        local fm = r(found)
        display as text "  case E: replay rc `rr', e(V) line of the matrices default found `fm' (expected 1)"
        if `rr' != 0 | `fm' != 1 {
            display as error "test_rotse_post 4: replay with e(V) not as specified"
            local ++nfail
        }
        capture noisily dwlsmv, level(90)
        if _rc {
            display as error "test_rotse_post 4: replay with level(90) rc " _rc
            local ++nfail
        }
        capture dwlsmv, standardized
        if _rc != 198 {
            display as error "test_rotse_post 4: replay with standardized rc " _rc ", not 198 (A3.11)"
            local ++nfail
        }
        tempname V0
        matrix `V0' = e(V)
        estimates restore tppE
        if mreldif(e(V), `V0') != 0 {
            display as error "test_rotse_post 4: estimates restore changed e(V)"
            local ++nfail
        }
    }
}

* ---- 5  m = 1, df = 0: nothing to rotate; e(V) posted; e(rank) = n* = 9 -----------------
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
capture noisily dwlsmv (F1 -> y1-y3), efa(F1)
local rc5 = _rc
capture confirm matrix e(V)
local vrc5 = _rc
if `rc5' | `vrc5' {
    display as error "test_rotse_post 5: rc `rc5', e(V) confirm rc `vrc5'"
    local ++nfail
}
else {
    tempname V5 V5t
    matrix `V5' = e(V)
    matrix `V5t' = `V5''
    mata: st_local("sh5", invtokens(strofreal(tpp_vshape(), "%21.0g")))
    local rk5 : word 1 of `sh5'
    display as text "  m = 1: e(V) " rowsof(`V5') " x " colsof(`V5') ", e(rank) " e(rank) ", numerical rank `rk5'"
    if rowsof(`V5') != 9 | e(rank) != 9 | `rk5' != 9 | mreldif(`V5', `V5t') != 0 | "`e(rot_nose)'" != "" {
        display as error "test_rotse_post 5: m = 1 e(V) not as specified"
        local ++nfail
    }
}

* ---- 6  the A3.4 path ---------------------------------------------------------------------
* Case G with its (y1 = 0, y2 = 1) cell deleted (E-109, A4 finding, section 3):
* lavaan's unrotated fit has a Heywood loading (1.007) there. Which branch dwlsmv reaches is not
* known to the drafter; each is PRINTED, and the measured one is pinned below:
*   rc 3498           no converged start (A3.4 "unchanged"): nothing posted
*   rc 0, nose set    A3.4: rotated e(b), no e(V), e(rot_nose) names the reason, header says so
*   rc 0, nose empty  e(V) posted: the A3.4 path was NOT exercised by this data (coverage gap)
import delimited using "fixtures/data/case_G.csv", clear varnames(1) delimiters(",")
drop if y1 == 0 & y2 == 1
capture noisily dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
local rc6 = _rc
if `rc6' == 3498 {
    display as text "  6: BRANCH rc 3498 (no converged start); the A3.4 posting path is not exercised here"
    if "`e(cmd)'" != "" {
        display as error "test_rotse_post 6: rc 3498 but something was posted"
        local ++nfail
    }
}
else if `rc6' == 0 {
    capture confirm matrix e(V)
    local vrc6 = _rc
    if "`e(rot_nose)'" != "" {
        display as text "  6: BRANCH rotated e(b) without e(V), e(rot_nose) = `e(rot_nose)'"
        tpp_run "No standard errors" dwlsmv
        if `vrc6' != 111 | colsof(e(b)) != 19 | e(rank) != 17 | r(found) != 1 {
            display as error "test_rotse_post 6: A3.4 posting not as specified (e(V) rc `vrc6', e(b) " colsof(e(b)) " columns, e(rank) " e(rank) ", header found " r(found) ")"
            local ++nfail
        }
    }
    else {
        display as text "  6: BRANCH e(V) posted; the A3.4 path was NOT exercised by this data"
        if `vrc6' != 0 {
            display as error "test_rotse_post 6: e(rot_nose) empty but no e(V)"
            local ++nfail
        }
    }
}
else {
    display as error "test_rotse_post 6: rc `rc6', neither 0 nor 3498"
    local ++nfail
}

* Pinned, CHOSEN AFTER MEASURING (main session, 2026-11-17, StataMP 18 at 6 and 4 processors):
* this data reaches the A3.4 posting path with e(rot_nose) = unconverged (the unrotated fit does not
* converge). Any other branch now fails, so the no-e(V) path cannot go unexercised silently.
if !(`rc6' == 0 & "`e(rot_nose)'" == "unconverged") {
    display as error "test_rotse_post 6: expected the A3.4 path with e(rot_nose) = unconverged (measured 2026-11-17); got rc `rc6', e(rot_nose) '`e(rot_nose)''"
    local ++nfail
}

* ---- 7  no converged start: unchanged, rc 3498, nothing posted (A3.4, 2026-10-29) --------
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
ereturn clear
capture dwlsmv (F1 F2 F3 -> y1-y6), efa(F1 F2 F3)
local rc7 = _rc
display as text "  7: over-extraction rc `rc7', e(cmd) '`e(cmd)''"
if `rc7' != 3498 | "`e(cmd)'" != "" {
    display as error "test_rotse_post 7: rc `rc7' (expected 3498), e(cmd) '`e(cmd)''"
    local ++nfail
}

if `nfail' {
    display as error _n "test_rotse_post: `nfail' failure(s)"
    exit 9
}
display as result "test_rotse_post: the rotated e(V) of Cases E, F, G matches lavaan 0.7-3's SEs and covariance at the 7.2 / 7.2.1 tolerances; e(rank) = n_u; lincom and test follow e(V); the no-e(V) and rc 3498 paths hold"
exit 0
