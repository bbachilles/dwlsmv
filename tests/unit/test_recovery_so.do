*! Executed under StataMP 18 on 2026-11-20: all checks pass. Drafted by a read-only agent; reviewed, reworded
*! to the owner's failure note (2026-11-20) and run by the main session.
*!
*! Release B, decision B3.3 (spec Module 12, Release B, "Recovery (B3.3)"): D11's marker-metric
*! recovery for a second-order model. A factor scaled by 'metric unitvariance' with no lambda row
*! (G) is eligible, c_G = 1 in D11's transform (derivation @sec-so-recovery SO4.1); when the
*! unit-variance route fails and the converged first fit has a negative variance (SO4.4), the failure
*! note says so.
*!
*! Builds build/ldwlsmv.mlib first and runs against the library and src/ado. Exit 0 pass, 9 fail.
*!
*! Tolerances, all set BEFORE any run: 1e-12 (sigma identity, exact algebra); 1e-5 (the transformed
*! U optimum against the V optimum: two separately converged fits, a diagnostic tier, the maximum is
*! printed); 1e-6 (the recovered fit against the first fit: both converged optima of one problem, the
*! recovered one started from the transform); 5e-4 (the printed -0.3742 of test_psi_check.do check 5).
*!
*! Checks (each names the mutant it kills; R1-R6 are listed at the end):
*!   a1  V table, G over four on Case H's items: eligible, marker[G] == 0, marker[F_k] == the item row
*!       of F_k's marker (1, 4, 7, 10).                                    kills R1
*!   a2  U table (every factor unit variance, all loadings free): ineligible; F1's reason is "has no
*!       loading fixed at 1" (a first-order factor is not exempt), and no reason names G.   kills R1
*!   a3  the V table with psi_GG fixed at 2: ineligible, reason "has no indicator and its variance is
*!       fixed at 2, not 1".                                               kills R6
*!   a4  the V table with psi_GG free: ineligible, reason "variance is free".   kills R6'
*!       (a hand-edited table: Module 3 refuses a metric line next to an explicit psi line, so no
*!       canonical text builds a3 or a4; the edit is to par.value / par.free only, which is all
*!       dwlsmv_m11_receligible() reads.)
*!   a5  an item on G (MD: G has a lambda row, 'metric unitvariance G'): ineligible, ordinary
*!       factor rules, reason "factor G has no loading fixed at 1".        kills R1''
*!   b   the transform on Case H: U fit and V fit by dwls_fit(); dwlsmv_m11_rectransform() returns 0,
*!       c[G] == 1 exactly, c[F_k] is F_k's marker loading in the U fit, sigma(V table, transformed)
*!       equals sigma(U table, U fit) to 1e-12, and the transformed theta is within 1e-5 of the V fit.
*!                                                                         kills R2
*!   c   dwlsmv_m11_recover() on Case H's V fit with converged set to 0 by hand (the trigger's
*!       "did not converge"; test_marker_recovery.do check 12 hand-sets flags the same way):
*!       "recovered", the note's restart sentence, the final fit converged and within 1e-6 of the V
*!       fit.                                                              kills R1, R2
*!   e1  a failed second-order recovery with a PROPER first fit: Case H's V fit converged, the trigger
*!       forced through A[1,1] = 1, a text dwls_identify() refuses at the unit-variance step: "failed",
*!       the note has no "cannot represent".                               kills R3
*!   e2  the same on Case I's improper V fit: "failed", the full two-line sentence for a disturbance
*!       variance, after the failure line, and the note printed for the reader.   kills R3 (inverse)
*!   d   command level on Case I, G over three: metric(marker): e(recovery) "failed", the sentence in
*!       the output, e(b)'s /var(e.F1) negative (-0.3742) and e(converged) 1; metric(unitvariance):
*!       e(recovery) "none", the sentence absent.
*!
*! NOT COVERED (no fixture): the "a negative variance" wording (an exogenous negative psi); a failed
*! first-ORDER recovery with a negative variance (R4); a psi inside its band, negative in sign (R5).
*!
*! Mutants: R1 the G exemption removed from receligible; R2 c[G] = 2 in rectransform; R3 the negative
*! psi condition dropped (the sentence on every failed second-order recovery); R4 the second-order
*! condition dropped; R5 the class replaced by a raw sign; R6 the psi_GG == 1 test loosened to "fixed";
*! R7 the generic note kept beside the new one. Run 2026-11-20 in a scratch worktree: R1 (a1-a4, d1, e2),
*! R2 (b), R3 (e1), R6 (a3; it aborted at 3200 until tsr_has() made the check fail cleanly) and R7
*! (e2) killed; R4 survives (NOT COVERED above; a coverage gap, no code change, survivor rule); R5 not run.
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

local nfail 0

* run a command, capturing what it prints; r(rc), r(found1), r(found2) for two phrases (the log's
* lines joined, blanks collapsed, so a sentence wrapped over two lines is found)
capture program drop tsr_run
program tsr_run, rclass
    gettoken ph1 0 : 0
    gettoken ph2 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tsrlog)
    capture noisily `0'
    local rc = _rc
    quietly log close tsrlog
    mata: st_local("txt", stritrim(invtokens(cat(st_local("lg"))', " ")))
    mata: st_local("f1", strofreal(strpos(st_local("txt"), st_local("ph1")) > 0))
    mata: st_local("f2", strofreal(strpos(st_local("txt"), st_local("ph2")) > 0))
    return scalar rc     = `rc'
    return scalar found1 = `f1'
    return scalar found2 = `f2'
end

mata:

string scalar tsr_join(string colvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= rows(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// G over nk first-order factors of three items each (items y1..y(3nk), then xitems extra), metric
// "U" (every factor unit variance, all loadings free) or "V" (item markers, psi_k free, G unit
// variance). The same builder as test_identify_second_order.do's tso_lines().
string colvector tsr_lines(real scalar nk, string scalar met, real scalar xitems)
{
    string colvector L
    string scalar fl, il, f
    real scalar k, j, y

    fl = "factor"
    for (k = 1; k <= nk; k++) {
        fl = fl + " F" + strofreal(k)
    }
    fl = fl + " G"
    il = "ind"
    for (y = 1; y <= 3 * nk + xitems; y++) {
        il = il + " y" + strofreal(y)
    }
    L = (fl \ il)
    y = 0
    for (k = 1; k <= nk; k++) {
        f = "F" + strofreal(k)
        for (j = 1; j <= 3; j++) {
            y++
            if (met == "V" & j == 1) {
                L = L \ ("metric marker " + f + " y" + strofreal(y))
            }
            else {
                L = L \ ("lambda " + f + " y" + strofreal(y) + " free")
            }
        }
        if (met == "V") {
            L = L \ ("psi " + f + " " + f + " free")
        }
        else {
            L = L \ ("metric unitvariance " + f)
        }
        L = L \ ("beta " + f + " G free")
    }
    L = L \ "metric unitvariance G"
    return(L)
}

// the ';'-separated single line dwlsmv.ado passes (recover reads its first token for the factor names)
string scalar tsr_canon(string colvector L)
{
    return(invtokens(L', ";"))
}

real scalar tsr_find(struct dwlsmv_parmap scalar par, real scalar pt, real scalar r, real scalar c)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == pt & par.row[i] == r & par.col[i] == c) {
            return(i)
        }
    }
    return(0)
}

// one assertion: prints, returns 1 on failure
// 1 if some reason contains txt; 0 for an empty list (strpos() on a 0 x 1 vector aborts at 3200,
// which a mutant that leaves no reasons reached: a clean failure instead of an abort).
real scalar tsr_has(string colvector reasons, string scalar txt)
{
    if (rows(reasons) == 0) {
        return(0)
    }
    return(any(strpos(reasons, txt) :> 0))
}

real scalar tsr_chk(string scalar label, real scalar ok)
{
    if (ok) {
        printf("{txt}  ok   %s\n", label)
        return(0)
    }
    printf("{err}  FAIL %s\n", label)
    return(1)
}

// ---- a: eligibility, from hand-written canonical text -----------------------------------------------
real scalar tsr_a()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar parV, parU, parX, parM
    real colvector marker
    string colvector reasons, fac
    string scalar canonV
    real scalar nf, el, g

    nf = 0
    m  = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(12, 1, 3))
    canonV = tsr_canon(tsr_lines(4, "V", 0))
    fac    = tokens(substr(tokens(canonV, ";")[1], 8, .))'

    // a1
    parV = dwls_identify(tsr_join(tsr_lines(4, "V", 0)), m[1])
    el   = dwlsmv_m11_receligible(parV, fac, marker, reasons)
    nf   = nf + tsr_chk("a1 V table, G over four, is eligible (reasons: " + invtokens(reasons', " | ") + ")",
                        el == 1 & rows(reasons) == 0)
    nf   = nf + tsr_chk("a1 marker[G] == 0 and marker[F_k] == the marker item row (1, 4, 7, 10)",
                        rows(marker) == 5 & marker[5] == 0 & marker[1] == 1 & marker[2] == 4 &
                        marker[3] == 7 & marker[4] == 10)

    // a2: the U table
    parU = dwls_identify(tsr_join(tsr_lines(4, "U", 0)), m[1])
    el   = dwlsmv_m11_receligible(parU, fac, marker, reasons)
    nf   = nf + tsr_chk("a2 U table is ineligible", el == 0)
    nf   = nf + tsr_chk("a2 F1: 'factor F1 has no loading fixed at 1' (a first-order factor is not exempt)",
                        tsr_has(reasons, "factor F1 has no loading fixed at 1"))
    nf   = nf + tsr_chk("a2 no reason names G (G is the exempt factor): " + invtokens(reasons', " | "),
                        !tsr_has(reasons, "factor G "))

    // a3, a4: psi_GG hand-edited
    g   = tsr_find(parV, 2, 5, 5)
    nf  = nf + tsr_chk("a3 setup: V table has a fixed psi_GG row at 1",
                       g > 0 & parV.free[g] == 0 & parV.value[g] == 1)
    parX = parV
    parX.value[g] = 2
    el  = dwlsmv_m11_receligible(parX, fac, marker, reasons)
    nf  = nf + tsr_chk("a3 psi_GG fixed at 2: ineligible with its own reason (" + invtokens(reasons', " | ") + ")",
                       el == 0 & tsr_has(reasons, "factor G has no indicator and its variance is fixed at 2, not 1"))
    parX = parV
    parX.free[g] = 1
    el  = dwlsmv_m11_receligible(parX, fac, marker, reasons)
    nf  = nf + tsr_chk("a4 psi_GG free: ineligible, 'variance is free' (" + invtokens(reasons', " | ") + ")",
                       el == 0 & tsr_has(reasons, "factor G has no indicator and its variance is free"))

    // a5: an item on G, p = 13
    m[1] = dwlsmv_statmap_build(J(13, 1, 3))
    parM = dwls_identify(tsr_join(tsr_lines(4, "V", 1) \ "lambda G y13 free"), m[1])
    el   = dwlsmv_m11_receligible(parM, fac, marker, reasons)
    nf   = nf + tsr_chk("a5 an item on G: G is an ordinary factor, ineligible with 'factor G has no loading fixed at 1' (" +
                        invtokens(reasons', " | ") + ")",
                        el == 0 & tsr_has(reasons, "factor G has no loading fixed at 1"))
    return(nf)
}

// ---- b, c, e1: Case H ---------------------------------------------------------------------------------
real scalar tsr_h()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar parU, parV, parr
    struct dwlsmv_fit scalar fitU, fitV, fitF, fitr
    real colvector marker, thm, c, sigU, sigV, cexp
    real matrix A, Bd
    string colvector reasons, fac
    string scalar canonV, rec, note
    real scalar nf, el, rc, k, i, ku, d1, d2

    nf     = 0
    st     = *findexternal("dwlsmv_m11_st")
    maps   = dwlsmv_statmap(1, 1)
    maps[1] = st.map
    canonV = tsr_canon(tsr_lines(4, "V", 0))
    fac    = tokens(substr(tokens(canonV, ";")[1], 8, .))'
    parV   = dwls_identify(tsr_join(tsr_lines(4, "V", 0)), st.map)
    parU   = dwls_identify(tsr_join(tsr_lines(4, "U", 0)), st.map)
    el     = dwlsmv_m11_receligible(parV, fac, marker, reasons)
    if (!el) {
        printf("{err}  precondition: the V table is not eligible (%s)\n", invtokens(reasons', " | "))
        return(1)
    }
    fitU = dwls_fit(st, parU)
    fitV = dwls_fit(st, parV)
    dwls_inference(fitV)
    nf = nf + tsr_chk("precondition: the U and V fits on Case H converged", fitU.converged == 1 & fitV.converged == 1)
    if (nf) {
        return(nf)
    }

    // ---- b
    rc = dwlsmv_m11_rectransform(parV, parU, fitU.theta, marker, thm, c)
    nf = nf + tsr_chk("b rectransform returns 0", rc == 0)
    nf = nf + tsr_chk("b c[G] == 1 exactly (a factor scaled by its variance)", rows(c) == 5 & c[5] == 1)
    cexp = J(4, 1, .)
    for (k = 1; k <= 4; k++) {
        ku = tsr_find(parU, 1, 3 * (k - 1) + 1, k)
        cexp[k] = fitU.theta[parU.fidx[ku]]
    }
    nf = nf + tsr_chk("b c[F_k] is F_k's marker loading in the U fit", max(abs(c[1::4] - cexp)) == 0)
    sigU = dwls_model(fitU.theta, parU, maps)
    sigV = dwls_model(thm, parV, maps)
    d1   = max(abs(sigU - sigV))
    printf("{txt}  b  max |sigma_U(theta_U) - sigma_V(transformed)| = %g\n", d1)
    nf = nf + tsr_chk("b sigma is unchanged by the transform, 1e-12", d1 <= 1e-12)
    d2 = max(abs(thm - fitV.theta))
    printf("{txt}  b  max |transformed U optimum - V optimum| = %g (diagnostic tier 1e-5, set before running)\n", d2)
    nf = nf + tsr_chk("b the transformed U optimum is the V optimum, 1e-5", d2 <= 1e-5)

    // ---- c
    fitF = fitV
    fitF.converged = 0
    A  = dwlsmv_admissible(fitV.theta, parV, st.map)
    Bd = dwlsmv_breaddiag(fitV.Delta, st.Gamma)
    nf = nf + tsr_chk("c precondition: the V fit converged fires no trigger on its own",
                      dwlsmv_m11_rectrigger(fitV, A, Bd) == "")
    parr.nall = -1
    fitr.converged = -1
    rec = dwlsmv_m11_recover(canonV, parV, fitF, A, Bd, parr, fitr, note)
    printf("{txt}  c  e(recovery) %s; note:\n%s", rec, note)
    nf = nf + tsr_chk("c e(recovery) is 'recovered'", rec == "recovered")
    nf = nf + tsr_chk("c the note says it did not converge and was restarted from the unit-variance solution",
                      strpos(note, "did not converge") > 0 & strpos(note, "it was restarted from the unit-variance solution") > 0)
    nf = nf + tsr_chk("c the final fit converged and is the V table's", fitr.converged == 1 & parr.npar == parV.npar)
    if (fitr.converged == 1 & rows(fitr.theta) == rows(fitV.theta)) {
        d2 = max(abs(fitr.theta - fitV.theta))
        printf("{txt}  c  max |recovered - first fit| = %g (tier 1e-6, set before running)\n", d2)
        nf = nf + tsr_chk("c the recovered fit is the first fit's optimum, 1e-6", d2 <= 1e-6)
    }
    else {
        nf = nf + tsr_chk("c theta has the first fit's length", 0)
    }

    // ---- e1: a failed recovery of a second-order model with a PROPER first fit
    A[1, 1] = 1
    nf = nf + tsr_chk("e1 precondition: the first fit has no negative variance (class 2 or 3)",
                      dwlsmv_m11_recnegpsi(fitV) == 0)
    rec = dwlsmv_m11_recover(canonV + ";lambda NoSuchFactor y1 free", parV, fitV, A, Bd, parr, fitr, note)
    printf("{txt}  e1 e(recovery) %s; note:\n%s", rec, note)
    nf = nf + tsr_chk("e1 'failed' with the failure line", rec == "failed" &
                      strpos(note, "Recovery attempt failed; first fit retained") > 0 &
                      strpos(note, "was not admissible") > 0)
    nf = nf + tsr_chk("e1 no negative-variance wording (the first fit is proper)", strpos(note, "cannot be represented") == 0)
    return(nf)
}

// ---- e2: Case I, the improper first fit, at the Mata level --------------------------------------------
real scalar tsr_i()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, parr
    struct dwlsmv_fit scalar fit, fitr
    real matrix A, Bd
    string scalar canon, rec, note, want, flat
    string rowvector L
    real scalar nf, p, j, wide

    nf    = 0
    st    = *findexternal("dwlsmv_m11_st")
    canon = tsr_canon(tsr_lines(3, "V", 0))
    par   = dwls_identify(tsr_join(tsr_lines(3, "V", 0)), st.map)
    fit   = dwls_fit(st, par)
    dwls_inference(fit)
    A     = dwlsmv_admissible(fit.theta, par, st.map)
    Bd    = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    nf    = nf + tsr_chk("e2 precondition: Case I's V fit converged with a negative DISTURBANCE variance (kind 2)",
                         fit.converged == 1 & dwlsmv_m11_recnegpsi(fit) == 2)
    if (nf) {
        return(nf)
    }
    A[1, 1] = 1
    parr.nall = -1
    fitr.converged = -1
    rec  = dwlsmv_m11_recover(canon + ";lambda NoSuchFactor y1 free", par, fit, A, Bd, parr, fitr, note)
    printf("{txt}  e2 e(recovery) %s; note as printed:\n%s", rec, note)
    // The owner's wording (2026-11-20), naming the reason; compared with the line breaks flattened,
    // since the paragraph is wrapped at 78 columns.
    flat = subinstr(note, char(10), " ")
    L    = tokens(note, char(10))
    wide = 0
    for (j = 1; j <= cols(L); j++) {
        if (L[j] != char(10) & strlen(L[j]) > 78) {
            wide = 1
        }
    }
    nf   = nf + tsr_chk("e2 'failed'", rec == "failed")
    nf   = nf + tsr_chk("e2 opens with the reason, naming F1",
                        substr(flat, 1, 91) == "The marker-metric fit produced an invalid solution: F1 has a negative disturbance variance.")
    nf   = nf + tsr_chk("e2 the retained-fit and parameterization sentences",
                        strpos(flat, "was attempted but failed. The first fit is retained. Its negative disturbance variance cannot be represented by the unit-variance parameterization, although the first fit itself converged.") > 0)
    nf   = nf + tsr_chk("e2 the step that failed is named, last",
                        strpos(flat, "The step that failed: dwls_identify() refused the unit-variance specification") > 0 &
                        strpos(flat, "The step that failed:") > strpos(flat, "itself converged."))
    nf   = nf + tsr_chk("e2 the generic failure wording is replaced", strpos(note, "Recovery attempt failed") == 0)
    nf   = nf + tsr_chk("e2 no line is wider than 78", wide == 0)
    nf   = nf + tsr_chk("e2 parr and fitr were not written", parr.nall == -1 & fitr.converged == -1)
    return(nf)
}

end

* ---- Case H: a, b, c, e1 ------------------------------------------------------------------------------
import delimited using "fixtures/data/case_H.csv", clear varnames(1) delimiters(",")
tempvar tu
generate byte `tu' = 1
mata: dwlsmv_m11_stage1("`tu'", "y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12", 0)
mata: st_local("na", strofreal(tsr_a()))
mata: st_local("nh", strofreal(tsr_h()))
mata: dwlsmv_m11_clear()
local nfail = `nfail' + `na' + `nh'

* ---- Case I: e2 (Mata level), d (command level) ---------------------------------------------------------
import delimited using "fixtures/data/case_I.csv", clear varnames(1) delimiters(",")
tempvar tu2
generate byte `tu2' = 1
mata: dwlsmv_m11_stage1("`tu2'", "y1 y2 y3 y4 y5 y6 y7 y8 y9", 0)
mata: st_local("ni", strofreal(tsr_i()))
mata: dwlsmv_m11_clear()
local nfail = `nfail' + `ni'

local imodel "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (G -> F1 F2 F3)"
local sent "produced an invalid solution: F1 has a negative"

* d1: metric(marker): failed, the sentence printed, the improper first fit reported
tsr_run "`sent'" "convergence not achieved" dwlsmv `imodel', metric(marker)
local rcd = r(rc)
local fnd = r(found1)
display as text "  d1 rc `rcd', e(recovery) `e(recovery)', e(converged) " e(converged) ", sentence found `fnd'"
if `rcd' != 0 {
    display as error "  FAIL d1 dwlsmv returned rc `rcd'"
    local ++nfail
}
else {
    if "`e(recovery)'" != "failed" {
        display as error "  FAIL d1 e(recovery) is `e(recovery)', want failed (the trigger did not fire on Case I?)"
        local ++nfail
    }
    if !`fnd' {
        display as error "  FAIL d1 the output lacks: `sent'"
        local ++nfail
    }
    if e(converged) != 1 {
        display as error "  FAIL d1 e(converged) is " e(converged) ": the first fit is not the converged improper one"
        local ++nfail
    }
    tempname b1
    matrix `b1' = e(b)
    * the column whose NAME is var(e.F1), whatever its equation label (the stripe is read in Mata)
    mata: st_local("j", strofreal(sum(selectindex(st_matrixcolstripe("`b1'")[., 2] :== "var(e.F1)")), "%9.0g"))
    if `j' == 0 {
        display as error "  FAIL d1 no column named var(e.F1) in e(b): " `"`: colfullnames `b1''"'
        local ++nfail
    }
    else {
        display as text "  d1 /var(e.F1) = " %9.6f `b1'[1, `j']
        * -0.3742 is the 4-digit value test_psi_check.do check 5 asserts for dwls_fit's first fit
        if !(`b1'[1, `j'] < 0) | abs(`b1'[1, `j'] + 0.3742) > 5e-4 {
            display as error "  FAIL d1 /var(e.F1) is not the improper first fit's -0.3742 within 5e-4"
            local ++nfail
        }
    }
}

* d2: metric(unitvariance): no recovery, the sentence absent
tsr_run "`sent'" "convergence not achieved" dwlsmv `imodel', metric(unitvariance)
local rcd = r(rc)
local fnd = r(found1)
display as text "  d2 rc `rcd', e(recovery) `e(recovery)', sentence found `fnd'"
if `rcd' != 0 | "`e(recovery)'" != "none" | `fnd' {
    display as error "  FAIL d2 metric(unitvariance): rc `rcd', e(recovery) `e(recovery)', sentence found `fnd'"
    local ++nfail
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_recovery_so: `nfail' failure(s)"
    exit 9
}
display as result "test_recovery_so: a factor scaled by its variance is eligible with c = 1, the transform keeps sigma, Case H recovers, and a failed recovery names the negative disturbance variance (Case I) and says nothing of it on a proper first fit"
exit 0
