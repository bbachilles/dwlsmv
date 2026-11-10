*! Executed under StataMP 18 on 2026-10-07 (suite green, 34 / 0 / 0).
*!
*! D11: marker-metric recovery through unit variance (spec Module 11 D11;
*! derivation @sec-marker-recovery; E-016 "Marker-metric start
*! trap and recovery through unit variance"). Drafted by an agent 2026-10-07
*! and rewritten against the implementation.
*!
*! Builds build/ldwlsmv.mlib first and runs against the library and src/ado.
*!
*! Checks:
*!   1  TRAP FIXTURE (fixtures/data/marker_trap.csv). Precondition: the first
*!      marker fit, run directly through dwls_identify()/dwls_fit(), does not
*!      converge, so the fixture still traps. Through the command:
*!      e(recovery) "recovered", e(converged) 1, e(b) not the first fit, the
*!      recovered note printed, and the first fit's admissibility note not
*!      printed (the notes are the reported fit's). e(iterations) is 2: the
*!      final fit skips ULS (D11, 2026-10-07); 2 was measured, not derived.
*!   2  theta-hat vs lavaan_marker_trap_ref.csv (lavaan's unit-variance fit,
*!      transformed) and lavaan_marker_trap_marker.csv (lavaan's own marker
*!      fit), every e(b) column, at spec 7.2's theta-hat tier (1e-5 + 1e-4).
*!   3  T3 vs lavaan's marker T3 at spec 7.2's T3 tier (1e-6 + 1e-5); df exact.
*!   4  SCOPE: metric(unitvariance) on the trap data posts "none".
*!   5  CASE E: "none"; e(b)'s free entries and e(iterations) equal the
*!      direct first fit's bit for bit (the code before D11 reported exactly
*!      that fit).
*!   6  INELIGIBLE WITH THE TRAP PATTERN: (MathAb -> q1@1 q2@1 ...).
*!      Precondition: its direct first fit does not converge. "ineligible";
*!      e(b)'s free entries are that first fit's bit for bit; e(converged) 0;
*!      the note names the second fixed loading.
*!   7  FAILURE PATH. No end-to-end way to force one exists without a
*!      test-only option, which D10 refuses. So:
*!      7a dwlsmv_m11_recover() itself, on a real trapped first fit, given a
*!         canonical text that dwls_identify() refuses at the unit-variance
*!         step: the abort inside the nested call is caught, "failed", the
*!         note says "Recovery attempt failed; first fit retained" and names
*!         the step, parr/fitr are not written, the external is removed, and
*!         a following command fit is unaffected.
*!      7b dwlsmv_m11_rectransform(): c_f = 0 and c_f missing return 1 (the
*!         map does not exist); c_f = 1e-300 returns 0 (no near-zero cutoff:
*!         that gate was closed 2026-10-07).
*!      7c dwlsmv_m11_recover() on the same trapped fit, given a first-fit
*!         table with one free row (a cross-loading) that the unit-variance
*!         table, built from the text, lacks: the transform aborts, "failed",
*!         and the note names the transform as the step that stopped.
*!      NOT EXERCISED: a unit-variance fit or a final fit that does not
*!      converge, and a theta2 substitution. No input is known that produces
*!      them from a real fit.
*!   8  PORT vs THE R REFERENCE (mr-impl), on committed fixtures only:
*!        trap:   lavaan's unit-variance estimates through
*!                dwlsmv_m11_rectransform() against
*!                lavaan_marker_trap_ref.csv, which R computed from the same
*!                estimates: 1e-12 abs + rel (derivation's identity tier, the
*!                CSV's printed digits the only gap);
*!        case E: str_unitvar transformed against str_marker, two separately
*!                converged lavaan fits: spec 7.2's theta-hat tier.
*!      Eligibility: the trap model and case E eligible; check 6's model
*!      ineligible naming two fixed loadings; the unit-variance table
*!      ineligible.
*!   9  PRECISION: starts needing each of the short form, %21.0g and %21x
*!      read back through dwls_identify() as the identical doubles, and at
*!      least one is written as %21x.
*!  10  THE FINAL FIT'S START (dwlsmv_m11_recstart(), no ULS): on the trap's
*!      Stage 1 and a table whose free values are all shifted by 0.25, every
*!      threshold start is s exactly and every other start is par.value
*!      exactly.
*!  11  TAU-HAT = S_TAU ON THE RECOVERED FIT: dwlsmv_m11_recover() on the
*!      trap, "recovered", then max |s_tau - tau-hat| <= 1e-15, the tier of
*!      test_fit.do check 6, which covers dwls_fit() but not the final fit's
*!      dwlsmv_fit_from() path. Catches a change to Delta's threshold block,
*!      to s_tau, or to the final fit's threshold start.
*!  12  THE TRIGGER, DETERMINISTICALLY: dwlsmv_m11_rectrigger() on case E's
*!      real converged first fit (precondition: it returns nothing), with
*!      hand-set flags: an admissibility flag of 1 in the last row and in the
*!      first row fires "was not admissible"; a missing flag and bread row 2
*!      do not; the two- and three-cause wordings are exact. On real fits the
*!      admissibility branch is masked: every trapped fit also fails to
*!      converge.
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
capture program drop tmr_run
program tmr_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tmrlog)
    capture noisily `0'
    local rc = _rc
    quietly log close tmrlog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end

* ---- reference fixtures: string columns, comma delimiter, strtoreal in Mata ----
import delimited using "fixtures/expected/lavaan_marker_trap_ref.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TMR_ref = st_sdata(., ("label", "est"))
import delimited using "fixtures/expected/lavaan_marker_trap_marker.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TMR_mk = st_sdata(., ("label", "est"))
import delimited using "fixtures/expected/lavaan_marker_trap_uv.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TMR_uv = st_sdata(., ("label", "est"))
import delimited using "fixtures/expected/lavaan_marker_trap_test.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TMR_test = st_sdata(., ("label", "value"))
import delimited using "fixtures/expected/lavaan_case_E_parameters.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TMR_Euv = select(st_sdata(., ("label", "est")), st_sdata(., "fit") :== "str_unitvar")
mata: TMR_Emk = select(st_sdata(., ("label", "est")), st_sdata(., "fit") :== "str_marker")

mata:
// lavaan label -> D6 e(b) name ("" for rows e(b) does not carry)
string scalar tmr_d6(string scalar lab, string scalar endog, string scalar facs)
{
    real scalar k
    string scalar a, b

    if (strpos(lab, "~*~") | substr(lab, -2, 2) == "~1") {
        return("")
    }
    if (strpos(lab, "|t")) {
        k = strpos(lab, "|t")
        return("/" + substr(lab, 1, k - 1) + ":cut" + substr(lab, k + 2, .))
    }
    if (strpos(lab, "=~")) {
        k = strpos(lab, "=~")
        return(substr(lab, k + 2, .) + ":" + substr(lab, 1, k - 1))
    }
    if (strpos(lab, "~~")) {
        k = strpos(lab, "~~")
        a = substr(lab, 1, k - 1)
        b = substr(lab, k + 2, .)
        if (a != b | !anyof(tokens(facs), a)) {
            return("")
        }
        return(a == endog ? "/var(e." + a + ")" : "/var(" + a + ")")
    }
    if (strpos(lab, "~")) {
        k = strpos(lab, "~")
        return(substr(lab, 1, k - 1) + ":" + substr(lab, k + 1, .))
    }
    return("")
}

real scalar tmr_lookup(string matrix L, string scalar lab)
{
    real colvector i

    i = selectindex(L[., 1] :== lab)
    if (rows(i) != 1) {
        return(.)
    }
    return(strtoreal(L[i[1], 2]))
}

// spec 7.2: accept when |a - b| <= abs + rel*|b|; a missing a or b fails
real scalar tmr_tier(real scalar a, real scalar b, real scalar ea, real scalar er)
{
    if (missing(a) | missing(b)) {
        return(0)
    }
    return(abs(a - b) <= ea + er * abs(b))
}

// e(b) against every lavaan row with a D6 name; every e(b) column must match
real scalar tmr_eb_vs(string matrix L, string scalar endog, string scalar facs,
                      real scalar ea, real scalar er, string scalar tag)
{
    real rowvector b
    string matrix nm
    string colvector full
    real colvector j
    real scalar i, nbad, nmatched, worst, a, r
    string scalar d6

    b    = st_matrix("e(b)")
    nm   = st_matrixcolstripe("e(b)")
    full = subinstr(nm[., 1] :+ ":" :+ nm[., 2], "/:", "/")
    nbad = 0
    nmatched = 0
    worst = 0
    for (i = 1; i <= rows(L); i++) {
        d6 = tmr_d6(L[i, 1], endog, facs)
        if (d6 == "") {
            continue
        }
        j = selectindex(full :== d6)
        if (rows(j) != 1) {
            printf("{err}   %s: %s (%s) not in e(b)\n", tag, L[i, 1], d6)
            nbad++
            continue
        }
        nmatched++
        a = b[j[1]]
        r = strtoreal(L[i, 2])
        if (abs(a - r) > worst) {
            worst = abs(a - r)
        }
        if (!tmr_tier(a, r, ea, er)) {
            printf("{err}   %s: %s e(b) %18.12g lavaan %18.12g\n", tag, d6, a, r)
            nbad++
        }
    }
    if (nmatched != cols(b)) {
        printf("{err}   %s: %f of %f e(b) columns matched\n", tag, nmatched, cols(b))
        nbad++
    }
    printf("{txt}   %s: %f rows, max |gap| %g, %f outside the tier\n", tag, nmatched, worst, nbad)
    return(nbad)
}

// the first marker fit, directly: locals fc, fi; returns theta-hat
real colvector tmr_firstfit(string scalar canon, string scalar ind, real colvector K)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit

    st  = dwls_sample_stats(st_data(., tokens(ind)), K)
    par = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit = dwls_fit(st, par)
    st_local("fc", strofreal(fit.converged))
    st_local("fi", strofreal(fit.iterations))
    return(fit.theta)
}

// e(b)'s free entries (columns e(Cns) does not fix), in theta order
real colvector tmr_freeb()
{
    real rowvector b
    real matrix C
    real colvector free
    real scalar j

    b    = st_matrix("e(b)")
    C    = st_matrix("e(Cns)")
    free = J(cols(b), 1, 1)
    for (j = 1; j <= rows(C); j++) {
        free = free :* (C[j, 1..cols(b)]' :== 0)
    }
    return(select(b', free))
}

real scalar tmr_same(real colvector th)
{
    real colvector fb

    fb = tmr_freeb()
    if (rows(fb) != rows(th)) {
        return(0)
    }
    printf("{txt}   max |e(b) free - first-fit theta| %g\n", max(abs(fb - th)))
    return(fb == th)
}

// e(b)'s free-entry names, in theta order
string colvector tmr_freenames()
{
    string matrix nm
    string colvector full
    real matrix C
    real colvector free
    real scalar j

    nm   = st_matrixcolstripe("e(b)")
    full = subinstr(nm[., 1] :+ ":" :+ nm[., 2], "/:", "/")
    C    = st_matrix("e(Cns)")
    free = J(rows(full), 1, 1)
    for (j = 1; j <= rows(C); j++) {
        free = free :* (C[j, 1..rows(full)]' :== 0)
    }
    return(select(full, free))
}

// Port: lavaan's unit-variance estimates, in the unit-variance table's theta
// order (namesU), through dwlsmv_m11_rectransform(), against a marker
// reference by name (namesM, the marker table's theta order).
real scalar tmr_port(string scalar canon, string scalar ind, real colvector K,
                     string colvector namesU, string colvector namesM,
                     string matrix LU, string matrix LM, string scalar endog,
                     string scalar facs, real scalar ea, real scalar er, string scalar tag)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar pm, pu
    real colvector thU, thM, c, marker, j
    string colvector reasons, lu, lm, fac
    real scalar i, nbad, worst, ref, rc

    st  = dwls_sample_stats(st_data(., tokens(ind)), K)
    pm  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fac = tokens(substr(tokens(canon, ";")[1], 8, .))'
    if (!dwlsmv_m11_receligible(pm, fac, marker, reasons)) {
        printf("{err}   %s: an eligible model judged ineligible\n", tag)
        return(1)
    }
    pu = dwls_identify(subinstr(dwlsmv_m11_recuvtext(canon), ";", char(10)), st.map)
    if (pu.npar != rows(namesU) | pm.npar != rows(namesM)) {
        printf("{err}   %s: table sizes %f, %f against names %f, %f\n", tag,
               pu.npar, pm.npar, rows(namesU), rows(namesM))
        return(1)
    }
    lu = J(rows(LU), 1, "")
    for (i = 1; i <= rows(LU); i++) {
        lu[i] = tmr_d6(LU[i, 1], endog, facs)
    }
    lm = J(rows(LM), 1, "")
    for (i = 1; i <= rows(LM); i++) {
        lm[i] = tmr_d6(LM[i, 1], endog, facs)
    }
    thU = J(rows(namesU), 1, .)
    for (i = 1; i <= rows(namesU); i++) {
        j = selectindex(lu :== namesU[i])
        if (rows(j) == 1) {
            thU[i] = strtoreal(LU[j[1], 2])
        }
    }
    if (missing(thU)) {
        printf("{err}   %s: %f unit-variance value(s) unmatched\n", tag, missing(thU))
        return(1)
    }
    rc = dwlsmv_m11_rectransform(pm, pu, thU, marker, thM, c)
    if (rc != 0) {
        printf("{err}   %s: rectransform returned %f\n", tag, rc)
        return(1)
    }
    nbad  = 0
    worst = 0
    for (i = 1; i <= rows(namesM); i++) {
        j   = selectindex(lm :== namesM[i])
        ref = (rows(j) == 1 ? strtoreal(LM[j[1], 2]) : .)
        if (!missing(ref)) {
            if (abs(thM[i] - ref) > worst) {
                worst = abs(thM[i] - ref)
            }
        }
        if (!tmr_tier(thM[i], ref, ea, er)) {
            printf("{err}   %s: %s port %18.12g ref %18.12g\n", tag, namesM[i], thM[i], ref)
            nbad++
        }
    }
    printf("{txt}   %s: c = %s; %f rows, max |gap| %g, %f outside\n", tag,
           invtokens(strofreal(c', "%9.6f")), rows(namesM), worst, nbad)
    return(nbad)
}

// eligibility of a table expected ineligible: returns 1 on a wrong verdict,
// or when no reason contains want
real scalar tmr_inelig(string scalar canon, string scalar ind, real colvector K,
                       string scalar want, string scalar tag)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar pm
    real colvector marker
    string colvector reasons, fac
    real scalar el

    st  = dwls_sample_stats(st_data(., tokens(ind)), K)
    pm  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fac = tokens(substr(tokens(canon, ";")[1], 8, .))'
    el  = dwlsmv_m11_receligible(pm, fac, marker, reasons)
    printf("{txt}   %s: eligible %f; reasons: %s\n", tag, el, invtokens(reasons', " | "))
    if (el != 0 | rows(reasons) == 0) {
        return(1)
    }
    return(!any(strpos(reasons, want)))
}

// 7a: dwlsmv_m11_recover() on the real trapped first fit, with a canonical
// text dwls_identify() refuses at the unit-variance step (an undeclared
// factor on an extra line). The first fit's own table is passed as par, so
// the trigger and eligibility see the real model.
real scalar tmr_abort(string scalar canon, string scalar ind)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, parr
    struct dwlsmv_fit    scalar fit, fitr
    real matrix A, Bd
    string scalar rec, note, bad
    real scalar nbad

    nbad = 0
    st   = *findexternal("dwlsmv_m11_st")
    par  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit  = dwls_fit(st, par)
    dwls_inference(fit)
    A    = dwlsmv_admissible(fit.theta, par, st.map)
    Bd   = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    bad  = canon + ";lambda NoSuchFactor " + tokens(ind)[1] + " free"
    parr.nall = -1
    fitr.converged = -1
    rec  = dwlsmv_m11_recover(bad, par, fit, A, Bd, parr, fitr, note)
    printf("{txt}   7a. e(recovery) %s; note:\n%s", rec, note)
    if (rec != "failed") {
        nbad++
    }
    if (!strpos(note, "Recovery attempt failed; first fit retained")) {
        nbad++
    }
    if (!strpos(note, "refused the unit-variance specification (rc 3498)")) {
        nbad++
    }
    if (parr.nall != -1 | fitr.converged != -1) {
        printf("{err}   7a. parr/fitr were written on a failure\n")
        nbad++
    }
    if (findexternal("dwlsmv_m11_rec") != NULL) {
        printf("{err}   7a. the recovery external was left behind\n")
        nbad++
    }
    return(nbad)
}

// 10: the final fit's start, dwlsmv_m11_recstart(): every threshold is
// Stage 1's sample threshold exactly, every other free row is par.value
// exactly, on the session's Stage 1 output and the given table.
real scalar tmr_recstart(string scalar canon)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    real colvector x
    real scalar i, nt, nbad

    nbad = 0
    nt   = 0
    st   = *findexternal("dwlsmv_m11_st")
    par  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            par.value[i] = par.value[i] + 0.25
        }
    }
    x = dwlsmv_m11_recstart(st, par)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            if (par.ptype[i] == 4) {
                nt++
                if (x[par.fidx[i]] != st.s[i]) {
                    nbad++
                }
            }
            else if (x[par.fidx[i]] != par.value[i]) {
                nbad++
            }
        }
    }
    printf("{txt}  10. recstart: %f threshold(s) of %f at s exactly; %f mismatch(es)\n",
           nt, st.map.T, nbad)
    if (nt != st.map.T | nt == 0) {
        nbad++
    }
    return(nbad)
}

// 11: tau-hat = s_tau on the RECOVERED fit, which no longer goes through
// dwls_fit() (test_fit.do check 6 covers that path). Same tier as check 6:
// 1e-15, about 4 ulps at |tau| < 2; measured 0.0 (2026-10-07).
real scalar tmr_taugap(string scalar canon)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, parr
    struct dwlsmv_fit    scalar fit, fitr
    real matrix A, Bd
    string scalar rec, note
    real scalar gap, i

    st   = *findexternal("dwlsmv_m11_st")
    par  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit  = dwls_fit(st, par)
    dwls_inference(fit)
    A    = dwlsmv_admissible(fit.theta, par, st.map)
    Bd   = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    rec  = dwlsmv_m11_recover(canon, par, fit, A, Bd, parr, fitr, note)
    if (rec != "recovered") {
        printf("{err}  11. precondition: the trap was not recovered (%s)\n", rec)
        return(1)
    }
    // the PARAMETERS against s_tau, not the residual s - sigma: a change to
    // sigma_tau would keep the residual at 0 and move tau-hat. Table rows
    // 1..T are the tau rows in s order (Section 4.2), checked.
    gap = 0
    for (i = 1; i <= st.map.T; i++) {
        if (parr.ptype[i] != 4 | parr.free[i] != 1) {
            printf("{err}  11. table row %f is not a free threshold\n", i)
            return(1)
        }
        gap = max((gap, abs(fitr.theta[parr.fidx[i]] - st.s[i])))
    }
    printf("{txt}  11. recovered fit: max |s_tau - tau-hat| = %g over %f thresholds\n", gap, st.map.T)
    return(gap > 1e-15)
}

real scalar tmr_trig1(string scalar tag, struct dwlsmv_fit scalar fit,
                      real matrix A, real matrix Bd, string scalar want)
{
    string scalar got

    got = dwlsmv_m11_rectrigger(fit, A, Bd)
    printf("{txt}  12. %-28s -> %s\n", tag, (got == "" ? "(none)" : got))
    return(got != want)
}

// 12: the trigger, deterministically, on case E's real converged first fit
// with hand-set diagnostic flags (the admissibility branch is otherwise
// masked: every trapped fit also fails to converge).
real scalar tmr_trigger(string scalar canon, string scalar ind, real colvector K)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit, f0
    real matrix A, Bd, A2, B2
    real scalar nbad

    nbad = 0
    st   = dwls_sample_stats(st_data(., tokens(ind)), K)
    par  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit  = dwls_fit(st, par)
    dwls_inference(fit)
    A    = dwlsmv_admissible(fit.theta, par, st.map)
    Bd   = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    if (fit.converged != 1 | dwlsmv_m11_rectrigger(fit, A, Bd) != "") {
        printf("{err}  12. precondition: case E's first fit is not a clean converged fit\n")
        return(1)
    }
    A2 = A
    A2[rows(A2), 1] = 1
    nbad = nbad + tmr_trig1("last admissibility row", fit, A2, Bd, "was not admissible")
    A2 = A
    A2[1, 1] = 1
    nbad = nbad + tmr_trig1("first admissibility row", fit, A2, Bd, "was not admissible")
    A2 = A
    A2[1, 1] = .
    nbad = nbad + tmr_trig1("missing flag", fit, A2, Bd, "")
    B2 = Bd
    B2[2, 1] = 1
    nbad = nbad + tmr_trig1("bread row 2", fit, A, B2, "")
    f0 = fit
    f0.converged = 0
    A2 = A
    A2[1, 1] = 1
    nbad = nbad + tmr_trig1("not converged, inadmissible", f0, A2, Bd,
                            "did not converge and was not admissible")
    B2 = Bd
    B2[3, 1] = 1
    nbad = nbad + tmr_trig1("all three", f0, A2, B2,
        "did not converge, was not admissible and could not determine all parameters reliably (the information in the fit was insufficient)")
    return(nbad)
}

// does not, so the unit-variance table (built from the text) has no row for
// it and dwlsmv_m11_rectransform() aborts with rc 3498. fit, A and Bd are
// the real trapped first fit's, so the trigger fires as in 7a.
real scalar tmr_abortx(string scalar canon, string scalar ind, string scalar xline)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, parx, parr
    struct dwlsmv_fit    scalar fit, fitr
    real matrix A, Bd
    string scalar rec, note
    real scalar nbad

    nbad = 0
    st   = *findexternal("dwlsmv_m11_st")
    par  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit  = dwls_fit(st, par)
    dwls_inference(fit)
    A    = dwlsmv_admissible(fit.theta, par, st.map)
    Bd   = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    parx = dwls_identify(subinstr(canon + ";" + xline, ";", char(10)), st.map)
    if (parx.npar != par.npar + 1) {
        printf("{err}   7c. the extra line did not add one free row\n")
        nbad++
    }
    parr.nall = -1
    fitr.converged = -1
    rec  = dwlsmv_m11_recover(canon, parx, fit, A, Bd, parr, fitr, note)
    printf("{txt}   7c. e(recovery) %s; note:\n%s", rec, note)
    if (rec != "failed") {
        nbad++
    }
    if (!strpos(note, "Recovery attempt failed; first fit retained")) {
        nbad++
    }
    if (!strpos(note, "the transform to the marker metric stopped with an error (rc 3498)")) {
        nbad++
    }
    if (parr.nall != -1 | fitr.converged != -1) {
        printf("{err}   7c. parr/fitr were written on a failure\n")
        nbad++
    }
    if (findexternal("dwlsmv_m11_rec") != NULL) {
        printf("{err}   7c. the recovery external was left behind\n")
        nbad++
    }
    return(nbad)
}

// 7b: the transform's refusal, on case E's tables at a hand-set point
real scalar tmr_cf(string scalar canon, string scalar ind)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar pm, pu
    real colvector marker, thU, thM, c
    string colvector reasons, fac
    real scalar i, k, nbad, r0, rm, rt

    st  = dwls_sample_stats(st_data(., tokens(ind)), J(6, 1, 3))
    pm  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    pu  = dwls_identify(subinstr(dwlsmv_m11_recuvtext(canon), ";", char(10)), st.map)
    fac = tokens(substr(tokens(canon, ";")[1], 8, .))'
    (void) dwlsmv_m11_receligible(pm, fac, marker, reasons)
    // F1's marker loading in the unit-variance table
    k = 0
    for (i = 1; i <= pu.nall; i++) {
        if (pu.ptype[i] == 1) {
            if (pu.row[i] == marker[1]) {
                if (pu.col[i] == 1) {
                    k = pu.fidx[i]
                }
            }
        }
    }
    thU    = J(pu.npar, 1, 0.5)
    thU[k] = 0
    r0     = dwlsmv_m11_rectransform(pm, pu, thU, marker, thM, c)
    thU[k] = .
    rm     = dwlsmv_m11_rectransform(pm, pu, thU, marker, thM, c)
    thU[k] = 1e-300
    rt     = dwlsmv_m11_rectransform(pm, pu, thU, marker, thM, c)
    printf("{txt}   7b. c_f = 0: %f; c_f missing: %f; c_f = 1e-300: %f (want 1, 1, 0)\n", r0, rm, rt)
    nbad = (r0 != 1) + (rm != 1) + (rt != 0)
    return(nbad)
}

// 9: starts written by dwlsmv_m11_recstarts() read back through
// dwls_identify() as the identical doubles, on values that need each of the
// three formats: short (0.5), %21.0g (1/3) and %21x (1/3 * 1e-4, which
// %21.0g does not carry; measured 2026-10-07)
real scalar tmr_prec(string scalar canon, string scalar ind)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar pm, pk
    real colvector th, v
    real scalar i, k, nbad, nx
    string scalar txt

    st  = dwls_sample_stats(st_data(., tokens(ind)), J(6, 1, 3))
    pm  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    v   = (0.5 \ 1/3 \ (1/3) * 1e-4 \ -(2/3) * 1e-7)
    th  = J(pm.npar, 1, 0)
    k   = 0
    for (i = 1; i <= pm.nall; i++) {
        if (pm.free[i]) {
            if (pm.ptype[i] != 4) {
                k++
                th[pm.fidx[i]] = v[mod(k - 1, rows(v)) + 1]
            }
        }
    }
    txt = dwlsmv_m11_recstarts(canon, pm, th)
    nx  = rows(selectindex(strpos(tokens(txt, ";")', "X") :> 0))
    pk  = dwls_identify(subinstr(txt, ";", char(10)), st.map)
    nbad = 0
    for (i = 1; i <= pk.nall; i++) {
        if (pk.free[i]) {
            if (pk.ptype[i] != 4) {
                if (pk.value[i] != th[pk.fidx[i]]) {
                    printf("{err}   9. %s: start %21x read back as %21x\n", pk.name[i],
                           th[pk.fidx[i]], pk.value[i])
                    nbad++
                }
            }
        }
    }
    printf("{txt}   9. %f start(s) checked, %f written as %%21x, %f not identical\n", k, nx, nbad)
    if (nx == 0) {
        nbad++
    }
    return(nbad)
}
end

* ---- the trap data ------------------------------------------------------------
import delimited using "fixtures/data/marker_trap.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring q* att*, replace
local mT "(MathAb -> q1-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb)"
local indT "q1 q2 q3 q4 q5 q6 q7 q8 att1 att2 att3 att4 att5"
mata: TMR_KT = (J(8, 1, 2) \ J(5, 1, 5))

* ---- 1. the trap: precondition, then recovered ------------------------------------
quietly _dwlsmv_translate , paths(`mT') metric(marker)
local canonT `"`r(canon)'"'
mata: TMR_th1 = tmr_firstfit(st_local("canonT"), "`indT'", TMR_KT)
display as text " 1. first marker fit (direct): converged `fc', iterations `fi'"
if `fc' != 0 {
    display as error " 1. precondition: the trap fixture no longer traps the first fit"
    local ++nfail
}
tmr_run "it was restarted from the unit-variance solution" dwlsmv `mT'
local nt1 = r(found)
display as text " 1. e(recovery) `e(recovery)', e(converged) " e(converged) ///
    ", e(iterations) " e(iterations) " (want 2, measured), note `nt1'"
if "`e(recovery)'" != "recovered" | e(converged) != 1 | e(iterations) != 2 | !`nt1' {
    local ++nfail
}
* the notes are the recovered fit's: the first fit's admissibility warning is gone
tmr_run "The fitted solution is not admissible" dwlsmv `mT'
local na1 = r(found)
display as text " 1. first fit's admissibility note shown after recovery: `na1' (want 0)"
if `na1' {
    local ++nfail
}
mata: st_local("same1", strofreal(tmr_same(TMR_th1)))
if `same1' {
    display as error " 1. e(b) is the trapped first fit"
    local ++nfail
}

* ---- 2. theta-hat vs the transformed reference and lavaan's marker fit -------------
mata: st_local("b2r", strofreal(tmr_eb_vs(TMR_ref, "MathAb", "MathAb MathAtt", 1e-5, 1e-4, " 2. vs ref")))
mata: st_local("b2m", strofreal(tmr_eb_vs(TMR_mk,  "MathAb", "MathAb MathAtt", 1e-5, 1e-4, " 2. vs lavaan marker")))
if `b2r' | `b2m' {
    local ++nfail
}

* ---- 3. T3 and df vs lavaan ---------------------------------------------------------
local T3rec = e(T3)
mata: st_local("T3l", strofreal(tmr_lookup(TMR_test, "marker.chisq.scaled"), "%21.15g"))
mata: st_local("dfl", strofreal(tmr_lookup(TMR_test, "marker.df.scaled"), "%21.15g"))
mata: st_local("ok3", strofreal(tmr_tier(`T3rec', `T3l', 1e-6, 1e-5)))
display as text " 3. T3 " %16.10f `T3rec' " lavaan " %16.10f `T3l' "; df " e(df_T3) " lavaan `dfl'"
if !`ok3' | e(df_T3) != `dfl' {
    local ++nfail
}
mata: TMR_nM = tmr_freenames()

* ---- 4. scope: metric(unitvariance) is never recovered --------------------------------
quietly dwlsmv `mT', metric(unitvariance)
display as text " 4. metric(unitvariance) on the trap data: e(recovery) `e(recovery)'"
if "`e(recovery)'" != "none" {
    local ++nfail
}
mata: TMR_nU = tmr_freenames()

* ---- 6. ineligible with the trap pattern ---------------------------------------------
local mI "(MathAb -> q1@1 q2@1 q3-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb)"
quietly _dwlsmv_translate , paths(`mI') metric(marker)
local canonI `"`r(canon)'"'
mata: TMR_thI = tmr_firstfit(st_local("canonI"), "`indT'", TMR_KT)
display as text " 6. ineligible model, first fit (direct): converged `fc'"
if `fc' != 0 {
    display as error " 6. precondition: the ineligible variant no longer trips the trigger"
    local ++nfail
}
tmr_run "factor MathAb has 2 loadings fixed at nonzero values" dwlsmv `mI'
local nt6 = r(found)
display as text " 6. e(recovery) `e(recovery)', e(converged) " e(converged) ", note `nt6'"
mata: st_local("same6", strofreal(tmr_same(TMR_thI)))
if "`e(recovery)'" != "ineligible" | e(converged) != 0 | !`same6' | !`nt6' {
    local ++nfail
}

* ---- 8a. port and eligibility on the trap model ---------------------------------------
quietly _dwlsmv_translate , paths(`mT') metric(unitvariance)
local canonTU `"`r(canon)'"'
mata: st_local("p8T", strofreal(tmr_port(st_local("canonT"), "`indT'", TMR_KT, TMR_nU, TMR_nM, TMR_uv, TMR_ref, "MathAb", "MathAb MathAtt", 1e-12, 1e-12, " 8. port, trap")))
mata: st_local("e8I", strofreal(tmr_inelig(st_local("canonI"), "`indT'", TMR_KT, "2 loadings fixed", " 8. eligibility, second fixed loading")))
mata: st_local("e8U", strofreal(tmr_inelig(st_local("canonTU"), "`indT'", TMR_KT, "variance fixed", " 8. eligibility, the unit-variance table")))
if `p8T' | `e8I' | `e8U' {
    local ++nfail
}

* ---- 7a. a recovery-time abort is a failed recovery -----------------------------------
tempvar tu
generate byte `tu' = 1
mata: dwlsmv_m11_stage1("`tu'", "`indT'", 0)
mata: st_local("a7", strofreal(tmr_abort(st_local("canonT"), "`indT'")))
mata: st_local("a10", strofreal(tmr_recstart(st_local("canonT"))))
mata: st_local("a11", strofreal(tmr_taugap(st_local("canonT"))))
mata: st_local("a7c", strofreal(tmr_abortx(st_local("canonT"), "`indT'", "lambda MathAb att1 free")))
mata: dwlsmv_m11_clear()
quietly dwlsmv `mT'
display as text " 7a. the next command fit: e(recovery) `e(recovery)'"
if `a7' | `a7c' | `a10' | `a11' | "`e(recovery)'" != "recovered" {
    local ++nfail
}

* ---- 5. case E unchanged; 8b. the port on case E; 7b. the transform's refusal -------------
import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
local mE "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)"
local indE "y1 y2 y3 y4 y5 y6"
quietly _dwlsmv_translate , paths(`mE') metric(marker)
local canonE `"`r(canon)'"'
mata: TMR_thE = tmr_firstfit(st_local("canonE"), "`indE'", J(6, 1, 3))
quietly dwlsmv `mE'
mata: st_local("same5", strofreal(tmr_same(TMR_thE)))
display as text " 5. case E: e(recovery) `e(recovery)', e(iterations) " e(iterations) " vs first fit `fi', bit for bit `same5'"
if "`e(recovery)'" != "none" | !`same5' | e(iterations) != `fi' | e(converged) != 1 {
    local ++nfail
}
* ---- 12. the trigger's branches, deterministically, on case E's converged first fit
mata: st_local("a12", strofreal(tmr_trigger(st_local("canonE"), "`indE'", J(6, 1, 3))))
if `a12' {
    local ++nfail
}
mata: TMR_nM = tmr_freenames()
quietly dwlsmv `mE', metric(unitvariance)
mata: TMR_nU = tmr_freenames()
* lavaan's case E labels are lower case
mata: TMR_Emk[., 1] = subinstr(subinstr(TMR_Emk[., 1], "f1", "F1"), "f2", "F2")
mata: TMR_Euv[., 1] = subinstr(subinstr(TMR_Euv[., 1], "f1", "F1"), "f2", "F2")
mata: st_local("p8E", strofreal(tmr_port(st_local("canonE"), "`indE'", J(6, 1, 3), TMR_nU, TMR_nM, TMR_Euv, TMR_Emk, "F2", "F1 F2", 1e-5, 1e-4, " 8. port, case E")))
mata: st_local("c7", strofreal(tmr_cf(st_local("canonE"), "`indE'")))
mata: st_local("p9", strofreal(tmr_prec(st_local("canonE"), "`indE'")))
if `p8E' | `c7' | `p9' {
    local ++nfail
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_marker_recovery: `nfail' failure(s)"
    exit 9
}
display as result "test_marker_recovery: the trap is recovered to lavaan's marker solution and T3; case E is untouched; an ineligible model and a caught abort report the first fit"
exit 0
