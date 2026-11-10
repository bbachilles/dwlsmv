*! dwlsmv -- measure_uls.do: does the recovery's final fit need ULS?
*!
*! Executed under StataMP 18 on 2026-10-07 (against 794cd4a and against the
*! commit that made the final fit skip ULS). Drafted by a read-only agent and
*! fixed after running. A diagnostic: no pass/fail. Figures and the decision:
*! E-092, "The recovery's final fit skips ULS".
*!
*! Question (spec D11 Procedure step 4): what does skipping ULS cost or save
*! in the final marker fit of D11's recovery, and does it matter whether the
*! thresholds start at Stage 1's sample thresholds or at 0? dwlsmv_minimize()
*! is called directly.
*!
*! For each of the marker-trap fixture (recovery's final fit) and case E:
*!   first   dwls_fit() on the marker table, no starts: the current path, printed
*!           in %21x so it can be compared bit for bit with a later run, and
*!           (case E) with dwlsmv's own e(T3), e(iterations).
*!   (a)     D11's final fit as built: dwls_fit() from the transformed
*!           lambda/beta/psi starts, thresholds at Module 3's 0, i.e. ULS then
*!           DWLS. Baseline from E-092: 7 DWLS iterations, trap
*!           T3 59.9383675. The ULS iterations are shown by repeating the two
*!           dwlsmv_minimize() calls dwls_fit() makes (dwls_fit.mata:120-125).
*!   (b)     DWLS only, from the transformed lambda/beta/psi, thresholds at
*!           Stage 1's sample thresholds (the first T entries of s: the
*!           threshold block of Delta is the identity).
*!   (c)     DWLS only, from the transformed lambda/beta/psi, thresholds at 0.
*! Reported per variant: converged, DWLS iterations, F at the start, max|dtheta|
*! against (a), T3 in %21x and %21.0g, and T3 - T3(a).
*!
*! (b) and (c) rebuild the fit struct's fields the way dwls_fit() does
*! (dwls_fit.mata:132-143) so that dwls_inference() can give T3. That is
*! copied code in a measurement file, not a proposal.
*!
*! Run from the repository root, in StataMP, with the tree at the commit under
*! test:  do tools/measure_uls.do
*! Reads fixtures only; writes nothing. Not a suite test, no pass/fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later

version 16
clear all
set more off
set linesize 120

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

mata:

// starting vector as dwls_fit() builds it (dwls_fit.mata:109-114)
real colvector ms_start(struct dwlsmv_parmap scalar par)
{
    real colvector s0
    real scalar i

    s0 = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            s0[par.fidx[i]] = par.value[i]
        }
    }
    if (missing(s0)) {
        _error("measure_uls: fidx does not number the free rows")
    }
    return(s0)
}

// thresholds at Stage 1's sample thresholds: rows 1..T of the table are the
// tau rows, in the index map's order (dwls_identify.mata:536-548)
void ms_thr(struct dwlsmv_parmap scalar par, struct dwlsmv_stage1 scalar st,
            real colvector start)
{
    real scalar i

    for (i = 1; i <= st.map.T; i++) {
        if (par.ptype[i] != 4 | par.free[i] != 1) {
            _error("measure_uls: a threshold row is not where expected")
        }
        start[par.fidx[i]] = st.s[i]
    }
}

// DWLS only from 'start'; the fit struct is rebuilt as dwls_fit() does it
struct dwlsmv_fit scalar ms_direct(struct dwlsmv_stage1 scalar st,
                                   struct dwlsmv_parmap scalar par,
                                   real colvector w, real colvector start,
                                   struct dwlsmv_fit scalar tmpl)
{
    struct dwlsmv_fit scalar f
    real colvector th
    real rowvector gr
    real matrix H
    real scalar conv, fm, it

    conv = dwlsmv_minimize(st.s, w, par, st.map, start, "d2", "nr", 300, th, fm, it)
    th   = dwlsmv_signcanon(th, par)
    dwlsmv_fdwls(0, th', st.s, w, par, st.map, fm, gr, H)
    f            = tmpl
    f.theta      = th
    f.Delta      = dwls_delta(th, par, st.map)
    f.fmin       = fm
    f.converged  = conv
    f.iterations = it
    f.T_raw      = (st.N - 1) * fm
    f.D_stat     = st.map.D - par.npar
    dwls_inference(f)
    return(f)
}

real scalar ms_f0(struct dwlsmv_stage1 scalar st, struct dwlsmv_parmap scalar par,
                  real colvector w, real colvector start)
{
    real rowvector gr
    real matrix H
    real scalar f0

    dwlsmv_fdwls(0, start', st.s, w, par, st.map, f0, gr, H)
    return(f0)
}

void ms_pr(string scalar lab, struct dwlsmv_fit scalar f, real scalar f0,
           struct dwlsmv_fit scalar a)
{
    real scalar d

    d = max(abs(f.theta - a.theta))
    printf("  %s\n", lab)
    printf("      converged %1.0f   DWLS iterations %3.0f   F at start %s\n",
           f.converged, f.iterations, strofreal(f0, "%21.0g"))
    printf("      max|dtheta| vs (a)  %s   (%s)\n", strofreal(d, "%21.0g"), strofreal(d, "%21x"))
    printf("      T3  %s\n", strofreal(f.T3, "%21x"))
    printf("      T3  %s\n", strofreal(f.T3, "%21.0g"))
    printf("      T3 - T3(a)  %s\n", strofreal(f.T3 - a.T3, "%21.0g"))
    printf("      fmin %s   fmin - fmin(a) %s\n", strofreal(f.fmin, "%21.0g"), strofreal(f.fmin - a.fmin, "%21.0g"))
}

void ms_model(string scalar lab, string scalar canon, string scalar ind,
              real colvector K)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, paru, park
    struct dwlsmv_fit    scalar fit0, fitu, fitA, fitB, fitC, fitD
    real colvector w, marker, c, thetam, startA, startB, startC, startD
    real colvector thu, thM
    real scalar rc, i, t2, fu, iu, fM, itM, convM, dM
    string colvector fac, reasons

    printf("\n======================================================================\n")
    printf("%s\n", lab)
    printf("======================================================================\n")

    st = dwls_sample_stats(st_data(., tokens(ind)), K)
    w  = 1 :/ diagonal(st.Gamma)
    printf("N %f  D %f  T %f\n", st.N, st.map.D, st.map.T)

    // ---- the first fit, no starts: the current path ---------------------
    par  = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    fit0 = dwls_fit(st, par)
    dwls_inference(fit0)
    printf("\nFIRST FIT (dwls_fit on the marker table, table starts): npar %f\n", par.npar)
    printf("  converged %1.0f   iterations %3.0f\n", fit0.converged, fit0.iterations)
    printf("  fmin   %s\n", strofreal(fit0.fmin, "%21x"))
    printf("  T_raw  %s\n", strofreal(fit0.T_raw, "%21x"))
    printf("  T3     %s\n", strofreal(fit0.T3, "%21x"))
    printf("  T3     %s\n", strofreal(fit0.T3, "%21.0g"))
    for (i = 1; i <= rows(fit0.theta); i++) {
        printf("  theta[%3.0f] %s\n", i, strofreal(fit0.theta[i], "%21x"))
    }
    st_numscalar("MS_first_conv", fit0.converged)
    st_numscalar("MS_first_it", fit0.iterations)
    st_numscalar("MS_first_T3", fit0.T3)

    // ---- D11 steps 1-4, as dwlsmv_m11_recrun() does them (recrun.mata:65-109)
    fac = tokens(substr(tokens(canon, ";")[1], 8, .))'
    if (!dwlsmv_m11_receligible(par, fac, marker, reasons)) {
        printf("NOT ELIGIBLE for recovery; variants (a)-(c) not run\n")
        return
    }
    paru = dwls_identify(subinstr(dwlsmv_m11_recuvtext(canon), ";", char(10)), st.map)
    fitu = dwls_fit(st, paru)
    printf("\nunit-variance fit: converged %1.0f   iterations %3.0f\n", fitu.converged, fitu.iterations)
    if (fitu.converged != 1) {
        printf("unit-variance fit did not converge; variants (a)-(c) not run\n")
        return
    }
    rc = dwlsmv_m11_rectransform(par, paru, fitu.theta, marker, thetam, c)
    if (rc != 0) {
        printf("transform returned %f; variants (a)-(c) not run\n", rc)
        return
    }
    park = dwls_identify(subinstr(dwlsmv_m11_recstarts(canon, par, thetam), ";", char(10)), st.map)
    t2 = 0
    for (i = 1; i <= park.nall; i++) {
        if (park.free[i] & park.ptype[i] != 4) {
            if (park.value[i] != thetam[park.fidx[i]]) {
                t2 = 1
            }
        }
    }
    printf("c_f: %s\n", invtokens(strofreal(c', "%12.6f")))
    printf("theta2 replaced the transformed starts: %1.0f\n", t2)

    // ---- (a) ULS then DWLS: dwls_fit() itself ---------------------------
    startA = ms_start(park)
    fitA   = dwls_fit(st, park)
    dwls_inference(fitA)
    printf("\nVARIANTS, final marker fit from the transformed starts (npar %f)\n", park.npar)
    printf("  thresholds at table start: max|tau start| %s\n", strofreal(max(abs(startA[|1 \ st.map.T|])), "%21.0g"))
    printf("  (a) ULS then DWLS = dwls_fit(st, park)\n")
    printf("      converged %1.0f   DWLS iterations %3.0f   F at start %s\n",
           fitA.converged, fitA.iterations, strofreal(ms_f0(st, park, w, startA), "%21.0g"))
    printf("      T3  %s\n", strofreal(fitA.T3, "%21x"))
    printf("      T3  %s\n", strofreal(fitA.T3, "%21.0g"))

    // the same two calls dwls_fit() makes, to show the ULS iterations
    (void) dwlsmv_minimize(st.s, J(st.map.D, 1, 1), park, st.map, startA,
                           "d2", "nr", 300, thu, fu, iu)
    convM = dwlsmv_minimize(st.s, w, park, st.map, thu, "d2", "nr", 300, thM, fM, itM)
    thM   = dwlsmv_signcanon(thM, park)
    dM    = max(abs(thM - fitA.theta))
    printf("      by hand: ULS iterations %3.0f, F_ULS %s; DWLS iterations %3.0f, converged %1.0f;\n",
           iu, strofreal(fu, "%21.0g"), itM, convM)
    printf("      by hand vs dwls_fit: max|dtheta| %s   iterations equal %1.0f\n",
           strofreal(dM, "%21.0g"), itM == fitA.iterations)
    printf("      (a) fmin %s  %s\n", strofreal(fitA.fmin, "%21.0g"), strofreal(fitA.fmin, "%21x"))
    printf("      F at the ULS point under W_D: %s\n", strofreal(ms_f0(st, park, w, thu), "%21.0g"))
    st_numscalar("MS_A_it", fitA.iterations)
    st_numscalar("MS_A_T3", fitA.T3)

    // ---- (b) DWLS only, thresholds at the sample thresholds -------------
    startB = startA
    ms_thr(park, st, startB)
    fitB = ms_direct(st, park, w, startB, fitA)
    ms_pr("(b) DWLS only: transformed lambda/beta/psi, thresholds at s[1..T]",
          fitB, ms_f0(st, park, w, startB), fitA)

    // ---- (c) DWLS only, thresholds at 0 ---------------------------------
    startC = startA
    fitC = ms_direct(st, park, w, startC, fitA)
    ms_pr("(c) DWLS only: transformed lambda/beta/psi, thresholds at table start (0)",
          fitC, ms_f0(st, park, w, startC), fitA)

    // ---- (d) DWLS only, thresholds at the unit-variance fit's -----------
    // thetam carries them unchanged (the transform passes tau rows through)
    startD = startA
    for (i = 1; i <= park.nall; i++) {
        if (park.free[i]) {
            if (park.ptype[i] == 4) {
                startD[park.fidx[i]] = thetam[park.fidx[i]]
            }
        }
    }
    printf("\n  max|tau_uv - s_tau|  %s   max|tau_uv - tau(b)| %s\n",
           strofreal(max(abs(startD[|1 \ st.map.T|] - st.s[|1 \ st.map.T|])), "%21.0g"),
           strofreal(max(abs(startD - startB)), "%21.0g"))
    printf("  max|tau_hat(a) - s_tau|  %s\n",
           strofreal(max(abs(fitA.theta[|1 \ st.map.T|] - st.s[|1 \ st.map.T|])), "%21.0g"))
    fitD = ms_direct(st, park, w, startD, fitA)
    ms_pr("(d) DWLS only: transformed lambda/beta/psi, thresholds at the unit-variance fit's",
          fitD, ms_f0(st, park, w, startD), fitA)
    printf("  (b) vs (d): max|dtheta| %s\n", strofreal(max(abs(fitB.theta - fitD.theta)), "%21x"))
}

end

* ---- the marker-trap fixture: the recovery's final fit ---------------------------
import delimited using "fixtures/data/marker_trap.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring q* att*, replace
local mT "(MathAb -> q1-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb)"
local indT "q1 q2 q3 q4 q5 q6 q7 q8 att1 att2 att3 att4 att5"
quietly _dwlsmv_translate , paths(`mT') metric(marker)
local canonT `"`r(canon)'"'
mata: ms_model("MARKER TRAP (fixtures/data/marker_trap.csv)", st_local("canonT"), "`indT'", (J(8, 1, 2) \ J(5, 1, 5)))

* the command: the first fit does not converge, so it recovers; its final fit is (a)
quietly dwlsmv `mT'
display as text _n "dwlsmv on the trap: e(recovery) `e(recovery)', e(converged) " e(converged) ///
    ", e(iterations) " e(iterations) ", e(T3) " %21.15g e(T3)
display as text "  e(iterations) equals (a)'s: " (e(iterations) == scalar(MS_A_it)) ///
    "   e(T3) equals (a)'s bit for bit: " (e(T3) == scalar(MS_A_T3))

* ---- case E ---------------------------------------------------------------------
import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
local mE "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)"
local indE "y1 y2 y3 y4 y5 y6"
quietly _dwlsmv_translate , paths(`mE') metric(marker)
local canonE `"`r(canon)'"'
mata: ms_model("CASE E", st_local("canonE"), "`indE'", J(6, 1, 3))

* the command on case E: its marker fit converges, so e(recovery) is none and
* the reported fit is the first fit above
quietly dwlsmv `mE'
display as text _n "dwlsmv on case E: e(recovery) `e(recovery)', e(converged) " e(converged) ///
    ", e(iterations) " e(iterations) ", e(T3) " %21.15g e(T3)
display as text "  e(iterations) equals the first fit's: " (e(iterations) == scalar(MS_first_it)) ///
    "   e(T3) equals the first fit's bit for bit: " (e(T3) == scalar(MS_first_T3))

adopath - "`here'/src/ado"
adopath - "`here'/build"
display as result _n "measure_uls: done (a measurement, no pass/fail)"
