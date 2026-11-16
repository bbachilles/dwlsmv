*! Executed 2026-11-14 by the main session (Stata 18 MP); one grid cell is a recorded exception (see check 1). Formerly: NOT YET EXECUTED. Drafted by a read-only agent 2026-11-14 (owner's decision of that date: the p-value leaves the condx_U
*! numerical reference comparison and its TRANSFORMATION is tested here instead). Target path tests/unit/test_pvalue_T3.do.
*! Remove this banner only after the file has been run and what broke has been fixed.
*!
*! WHAT IS TESTED. The p-value of T3 is posted in src/ado/dwlsmv.ado at line 455 (single-group and conditional
*! path, the one this file reaches) and at line 331 (the multi-group posting path, NOT reached here):
*!     ereturn scalar p_T3 = cond(df > 0, chi2tail(df, T3), .)
*! and printed by _dwlsmv_header as "Prob > chi2" with %9.4f. So the transformation is Stata's chi2tail(df, T3), and
*! at df = 0 production posts a MISSING value (not 1, not 0), which the header prints as ".".
*! difftest's p (src/mata/dwlsmv_difftest.mata:102) calls the same function from a different line; it is not the same
*! code path as p_T3 and is not covered here.
*!
*! WHY SEPARATE. A reference-fixture comparison of p cannot say whether a disagreement comes from T3 or from the
*! transformation, and no tier exists for p. T3 is compared in tests/unit/test_condx_fit_ref.do; here only
*! chi2tail is examined, against values that do not come from chi2tail.
*!
*! Checks
*!   1  chi2tail(df, T) against closed-form chi-square survival probabilities Q(df, T), df 1..12, T in
*!      (0.01, 0.5, 1, 3.84, 10, 33, 100, 400); 96 cells. Even df = 2m: Q = exp(-T/2) sum_{k<m} (T/2)^k / k!.
*!      df 1: Q = 2 normal(-sqrt(T)). Odd df: Q(df + 2) = Q(df) + (T/2)^(df/2) exp(-T/2) / Gamma(df/2 + 1), from df 1,
*!      with Gamma(3/2) = sqrt(pi)/2 and Gamma(x + 1) = x Gamma(x) carried in the term ratio, so no lngamma enters.
*!      The deciding criterion is the relative gap |chi2tail - Q| / Q <= tol(n); the bracket check below is a
*!      diagnostic that fails the file only if violated (it is a mathematical inequality, not a tolerance).
*!   1b odd-df bracket: Q(df - 1, T) < chi2tail(df, T) < Q(df + 1, T) with both bounds from the EVEN closed forms
*!      (no normal() involved; chi-square is stochastically increasing in df); for df 1 the lower bound is 0.
*!   1c six textbook constants, Q(df, T) = alpha at the tabulated critical value (below).
*!   2  command level: after dwlsmv (F -> y1 y2 y3 y4) (df 2 by hand count: D = 8 + 6 = 14, 8 + 3 + 1 = 12 free) and
*!      after the conditional fit (F -> y1 y2 y3 y4) (F <- x) (df 5, the hand count of test_condx_output.do):
*!      e(p_T3) == chi2tail(e(df_T3), e(T3)) EXACTLY (== on doubles: the same function on the same two doubles), 0 < p < 1,
*!      the header's "Prob > chi2" line shows the same value at %9.4f, and so does the replay.
*!   3  df = 0 (one factor on y1 y2 y3, K = 3 each: D = 6 + 3 = 9 and 6 + 2 + 1 = 9 free): e(df_T3) == 0 and
*!      e(p_T3) is missing, and the header prints "." on the Prob line, at estimation and on replay.
*!
*! TOLERANCE, derived from the arithmetic and not from a measurement. The reference is a sum of n positive terms,
*! each formed by at most n multiplications and divisions (one rounding each, unit roundoff u = eps/2, eps =
*! 2.2e-16 = epsilon(1)), times exp(-T/2) (a few u) and, for odd df, sqrt, pi and normal(). Positive terms do not
*! cancel, so the relative error of the reference is at most about (3n + 4) u. chi2tail is allowed the same budget
*! again, and the total is rounded up by a factor 4 to 8 eps (3n + 4) = 16 (3n + 4) u, n = number of terms (ceil(df/2)
*! for even df; (df + 1)/2 for odd df, counting the normal() term). For df 12 that is 8 * 2.2e-16 * 22 = 3.9e-14.
*! Odd df rests also on normal(), whose far-tail relative accuracy this file does not establish; if an odd-df
*! far-tail cell fails while the bracket check passes, the question for the owner is normal() versus chi2tail
*! there, not a loosened tolerance.
*! UNDERFLOW. The smallest reference in the grid is 2 normal(-20) = 5.5e-89 and exp(-200) = 1.4e-87, far above the
*! smallest normal double 2.2e-308, so no cell underflows. The test asserts every reference lies in (1e-300, 1] so
*! that a changed grid that underflows fails loudly instead of comparing 0 with 0.
*! CONSTANTS. Each critical value is given to 16 significant digits, so its own uncertainty is at most 4e-15 in T
*! (observed spread between two independent printings of the same quantile: 4e-15 at df 3, 2e-15 at df 1), and the
*! hazard f(T)/Q at these points is below 1 (0.6 for df 1 at 3.84), so the relative uncertainty of Q is at most 4e-15;
*! tol = 4e-15 + 8 eps (3n + 4).
*!   Q(1, 3.841458820694124) = 0.05, Q(2, 5.991464547107979) = 0.05, Q(3, 7.814727903251179) = 0.05,
*!   Q(4, 9.487729036781154) = 0.05, Q(1, 6.634896601021213) = 0.01, Q(2, 9.210340371976184) = 0.01 (= 2 ln 100,
*!   so Q = exp(-T/2) is exact for df 2).
*!   Source: the chi-square upper-tail critical values tabulated in Abramowitz and Stegun, Handbook of
*!   Mathematical Functions, Table 26.8, and in the NIST/SEMATECH e-Handbook of Statistical Methods, section 1.3.6.7.4;
*!   the 16-digit forms are those printed by R's qchisq(alpha, df, lower.tail = FALSE) (checked by the drafting agent
*!   in a scratch R run, 2026-11-14; qchisq's last digit differs from the quoted constant by up to 4e-15, which is the
*!   uncertainty used above). The printed digits are NOT re-derived here.
*!
*! WHAT THIS CERTIFIES AND DOES NOT. It certifies that (i) Stata's chi2tail agrees with independent closed forms to
*! the stated relative tolerance on 96 cells and six constants, (ii) production posts exactly that function of its own
*! df and T3, prints it, and posts missing at df 0. It certifies nothing about T3 itself, about the multi-group posting
*! line 331 (a separate call site with the same expression), about difftest's p, or about chi2tail for non-integer df
*! or df above 12. A pass here is not evidence for T3's value (tests/unit/test_condx_fit_ref.do) or for its
*! asymptotic validity.
*!
*! Return codes: 0 pass, 9 fail.
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

mata:

real scalar pv_chk(real scalar ok, string scalar what)
{
    if (!ok) {
        printf("{err}test_pvalue_T3: %s\n", what)
        return(1)
    }
    return(0)
}

// Closed-form survival probability Q(df, T) of a chi-square(df), and the number of terms n it summed.
// Returns (Q, n). Nothing here calls chi2tail or chi2.
real rowvector pv_ref(real scalar df, real scalar T)
{
    real scalar h, t, q, s, k, j, n

    h = T / 2
    if (mod(df, 2) == 0) {
        t = 1
        s = 1
        for (k = 1; k <= df / 2 - 1; k++) {
            t = t * h / k
            s = s + t
        }
        return((exp(-h) * s, df / 2))
    }
    // odd df: start from df 1, then Q(j + 2) = Q(j) + t(j), t(j) = h^(j/2) exp(-h) / Gamma(j/2 + 1);
    // t(1) = sqrt(h) exp(-h) / (sqrt(pi)/2), t(j + 2) = t(j) h / (j/2 + 1)
    q = 2 * normal(-sqrt(T))
    t = sqrt(h) * exp(-h) / (sqrt(pi()) / 2)
    n = 1
    for (j = 1; j <= df - 2; j = j + 2) {
        q = q + t
        t = t * h / (j / 2 + 1)
        n = n + 1
    }
    return((q, n))
}

// the relative tolerance of the header: 8 eps (3 n + 4)
real scalar pv_tol(real scalar n)
{
    return(8 * epsilon(1) * (3 * n + 4))
}

// ---- check 1 and 1b ----------------------------------------------------------------------------
real scalar pv_grid()
{
    real rowvector Tg, r, rl, ru
    real scalar nf, i, df, p, tol, rel, worst, lo, up, ncell

    Tg    = (0.01, 0.5, 1, 3.84, 10, 33, 100, 400)
    nf    = 0
    worst = 0
    ncell = 0
    for (df = 1; df <= 12; df++) {
        for (i = 1; i <= cols(Tg); i++) {
            r   = pv_ref(df, Tg[i])
            p   = chi2tail(df, Tg[i])
            tol = pv_tol(r[2])
            ncell = ncell + 1
            nf = nf + pv_chk(!missing(r[1]) & r[1] > 1e-300 & r[1] <= 1,
                             "check 1: the reference at df " + strofreal(df) + ", T " + strofreal(Tg[i]) + " is " + strofreal(r[1], "%21.15g") + ", outside (1e-300, 1]")
            if (!missing(r[1]) & r[1] > 1e-300) {
                rel = abs(p - r[1]) / r[1]
                if (!missing(rel)) {
                    worst = max((worst, rel / tol))
                }
                // KNOWN EXCEPTION, recorded after the first run (main session, 2026-11-14), NOT a loosened tier:
                // Stata 18's chi2tail(1, 400) = 5.50724823721231e-89 differs from the closed form by 2.9e-14
                // relative (2.35 times the derived tolerance). An asymptotic erfc series and R's pchisq and
                // pnorm (scratch R, 2026-11-14) agree with the closed form to 4e-16, so the gap is chi2tail's
                // far-tail accuracy at p = 5.5e-89, not the reference's. That one cell is asserted at 1e-13,
                // a bound CHOSEN AFTER MEASURING; the owner chose this documented exception over a post hoc exclusion
                // rule (2026-11-14); every other cell keeps the derived one.
                if (df == 1 & Tg[i] == 400) {
                    printf("check 1: known exception chi2tail(1, 400): relative gap %9.3e (asserted <= 1e-13)\n", rel)
                    nf = nf + pv_chk(!missing(rel) & rel <= 1e-13, "check 1: chi2tail(1, 400) exceeds the recorded exception bound 1e-13")
                    continue
                }
                nf = nf + pv_chk(!missing(rel) & rel <= tol,
                                 "check 1: chi2tail(" + strofreal(df) + ", " + strofreal(Tg[i]) + ") = " + strofreal(p, "%21.15g") +
                                 " vs closed form " + strofreal(r[1], "%21.15g") + ", relative gap " + strofreal(rel, "%9.3e") +
                                 " > " + strofreal(tol, "%9.3e"))
            }
            // 1b: odd df is bracketed by the even closed forms
            if (mod(df, 2) == 1) {
                if (df == 1) {
                    lo = 0
                }
                else {
                    rl = pv_ref(df - 1, Tg[i])
                    lo = rl[1]
                }
                ru = pv_ref(df + 1, Tg[i])
                up = ru[1]
                nf = nf + pv_chk(!missing(p) & p > lo & p < up,
                                 "check 1b: chi2tail(" + strofreal(df) + ", " + strofreal(Tg[i]) + ") = " + strofreal(p, "%21.15g") +
                                 " is not inside (" + strofreal(lo, "%21.15g") + ", " + strofreal(up, "%21.15g") + ")")
            }
        }
    }
    printf("check 1: %g cells; largest (relative gap / tolerance) = %9.3e (passes at <= 1)\n", ncell, worst)
    return(nf)
}

// ---- check 1c ----------------------------------------------------------------------------------
real scalar pv_consts()
{
    real matrix C
    real rowvector r
    real scalar nf, i, p, tol, rel

    C = (1, 3.841458820694124, 0.05 \
         2, 5.991464547107979, 0.05 \
         3, 7.814727903251179, 0.05 \
         4, 9.487729036781154, 0.05 \
         1, 6.634896601021213, 0.01 \
         2, 9.210340371976184, 0.01)
    nf = 0
    for (i = 1; i <= rows(C); i++) {
        p   = chi2tail(C[i, 1], C[i, 2])
        tol = 4e-15 + pv_tol(ceil(C[i, 1] / 2))
        rel = abs(p - C[i, 3]) / C[i, 3]
        printf("check 1c: Q(%g, %s) = %s, relative gap to %g is %9.3e (tolerance %9.3e)\n",
               C[i, 1], strofreal(C[i, 2], "%21.15g"), strofreal(p, "%21.15g"), C[i, 3], rel, tol)
        nf = nf + pv_chk(!missing(rel) & rel <= tol,
                         "check 1c: chi2tail(" + strofreal(C[i, 1]) + ", " + strofreal(C[i, 2], "%21.15g") + ") is " +
                         strofreal(p, "%21.15g") + ", not " + strofreal(C[i, 3]))
        // the closed form agrees with the constant too (df 2 exactly by construction, the others as a cross-check)
        r   = pv_ref(C[i, 1], C[i, 2])
        nf  = nf + pv_chk(abs(r[1] - C[i, 3]) / C[i, 3] <= tol,
                          "check 1c: the closed form at the constant (df " + strofreal(C[i, 1]) + ") is " + strofreal(r[1], "%21.15g"))
    }
    return(nf)
}

// the text after the equals sign on the first "Prob > chi2" line of a log file; "" if there is none
string scalar tp_pline(string scalar fn)
{
    string colvector L
    real scalar i, k

    L = cat(fn)
    for (i = 1; i <= rows(L); i++) {
        if (strpos(L[i], "Prob > chi2") > 0) {
            k = strpos(L[i], "=")
            if (k > 0) {
                return(strtrim(substr(L[i], k + 1, .)))
            }
        }
    }
    return("")
}

end

* ---- checks 1, 1b, 1c: the transformation against closed forms ------------------------------
capture noisily mata: st_local("n1", strofreal(pv_grid()))
if _rc {
    display as error "test_pvalue_T3: pv_grid() aborted, rc " _rc
    local ++nfail
}
else local nfail = `nfail' + `n1'
capture noisily mata: st_local("n1c", strofreal(pv_consts()))
if _rc {
    display as error "test_pvalue_T3: pv_consts() aborted, rc " _rc
    local ++nfail
}
else local nfail = `nfail' + `n1c'

* ---- run a command, keep its log: r(rc) and r(pline) -------------------------------------------
capture program drop tp_run
program tp_run, rclass
    tempfile lg
    quietly log using "`lg'", text replace name(tplog)
    capture noisily `0'
    local rc = _rc
    quietly log close tplog
    mata: st_local("pline", tp_pline(st_local("lg")))
    return scalar rc = `rc'
    return local pline `"`pline'"'
end

* ---- the data (the design of test_condx_output.do) ------------------------------------------------
set seed 20261114
quietly set obs 1500
quietly generate byte x = runiform() < .5
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
drop eta z1 z2 z3 z4

* ---- check 2: command level, an ordinary fit (df 2) and a conditional fit (df 5) ---------------------
local lab1 "ordinary fit (F -> y1 y2 y3 y4)"
local mod1 "(F -> y1 y2 y3 y4)"
local df1 2
local lab2 "conditional fit (F -> y1 y2 y3 y4) (F <- x)"
local mod2 "(F -> y1 y2 y3 y4) (F <- x)"
local df2 5
forvalues m = 1/2 {
    tp_run dwlsmv `mod`m''
    local rc = r(rc)
    local pl = r(pline)
    if `rc' {
        display as error " 2. `lab`m'': rc `rc'"
        local ++nfail
        continue
    }
    scalar pp  = e(p_T3)
    scalar pw  = chi2tail(e(df_T3), e(T3))
    local want = strtrim(string(e(p_T3), "%9.4f"))
    if e(df_T3) != `df`m'' | missing(e(T3)) | missing(e(p_T3)) | !(e(p_T3) > 0 & e(p_T3) < 1) {
        display as error " 2. `lab`m'': df " e(df_T3) " (want `df`m''), T3 " e(T3) ", p " e(p_T3) " (want 0 < p < 1)"
        local ++nfail
    }
    else if pp != pw {
        display as error " 2. `lab`m'': e(p_T3) " %21.17g pp " is not chi2tail(e(df_T3), e(T3)) = " %21.17g pw " bitwise"
        local ++nfail
    }
    else if "`pl'" != "`want'" {
        display as error " 2. `lab`m'': the header's Prob line shows '`pl'', e(p_T3) at %9.4f is '`want''"
        local ++nfail
    }
    else {
        tp_run dwlsmv
        if r(rc) | r(pline) != "`want'" {
            display as error " 2. `lab`m'': replay rc " r(rc) ", Prob line '" r(pline) "', wanted '`want''"
            local ++nfail
        }
        else display as text " 2. `lab`m'': e(p_T3) = chi2tail(df " e(df_T3) ", T3) exactly, header and replay show `want'"
    }
}

* ---- check 3: df = 0 (one factor on three K = 3 indicators: D = 9, 9 free) -----------------------------
tp_run dwlsmv (F -> y1 y2 y3)
local rc0 = r(rc)
local pl0 = r(pline)
if `rc0' {
    display as error " 3. df 0 fit (F -> y1 y2 y3): rc `rc0'; the df = 0 case is not established by this run"
    local ++nfail
}
else {
    if e(df_T3) != 0 | !missing(e(p_T3)) | "`pl0'" != "." {
        display as error " 3. df 0: df " e(df_T3) ", p " e(p_T3) " (want missing), header Prob line '`pl0'' (want .)"
        local ++nfail
    }
    else {
        tp_run dwlsmv
        if r(rc) | r(pline) != "." {
            display as error " 3. df 0 replay: rc " r(rc) ", Prob line '" r(pline) "' (want .)"
            local ++nfail
        }
        else display as text " 3. df 0: e(p_T3) is missing and the header prints ., at estimation and on replay"
    }
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_pvalue_T3: `nfail' failure(s)"
    exit 9
}
display as result "test_pvalue_T3: chi2tail agrees with closed-form chi-square survival probabilities (df 1 to 12, six textbook constants); the command posts and prints exactly chi2tail(df, T3), and a missing p at df 0"
exit 0

* ======================================================================================================
* MUTANTS (apply one at a time, run this file, expect exit 9; restore after)
*  M1 src/ado/dwlsmv.ado line 455 (the line with ereturn scalar p_T3 and seven spaces before the equals sign)
*       old: chi2tail(`df', `T3'), .)
*       new: chi2(`df', `T3'), .)
*       killed by check 2 (the CDF is posted, not the survival probability: the bitwise comparison with
*       chi2tail(e(df_T3), e(T3)) fails).
*  M2 src/ado/dwlsmv.ado line 455
*       old: chi2tail(`df', `T3'), .)
*       new: chi2tail(`df', `T3'), 0)
*       killed by check 3 (df 0 posts p = 0, not missing, and the header prints 0.0000).
*  M3 src/ado/dwlsmv.ado, _dwlsmv_header, the Prob line
*       old: "Prob > chi2" _col(67) "= " as result %9.4f e(p_T3)
*       new: "Prob > chi2" _col(67) "= " as result %9.3f e(p_T3)
*       killed by check 2 (the header shows 3 decimals where e(p_T3) at %9.4f has 4) and by check 3 only if it
*       also changed the missing print (it does not: both print .).
*  Not killed on purpose: line 331 (the multi-group posting) is a separate call site; a mutant there needs a
*  grouped fit, which this file does not run. Check 1 tests Stata's chi2tail, not dwlsmv code; it is killed by
*  editing a closed form (for example the exponent in pv_ref) or a tabulated constant, not by a production mutant.
