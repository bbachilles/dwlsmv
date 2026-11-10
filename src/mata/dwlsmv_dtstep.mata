*! dwlsmv -- DIFFTEST: one adjacent pair after checks 8 and 9 (Module 12)
*!
*!     real rowvector dwlsmv_dtstep(struct dwlsmv_parmap scalar par0,
*!                                  real colvector theta0,
*!                                  real scalar d, real scalar exempt,
*!                                  real scalar T0, real scalar T1,
*!                                  real scalar c, real scalar nrtol,
*!                                  real scalar df0,
*!                                  string scalar cc, real matrix K,
*!                                  string scalar inv0, string scalar inv1,
*!                                  string scalar name0, string scalar name1
*!                                  | string scalar note)
*!
*! Built 2026-10-13 (exemption gate step 2c). For one pair in dwlsmv_dtchain's
*! order (H0 = name0, the more restricted), with d and exempt from its info
*! row: check 11 on H0 with its caller's consequence sentence, the zero-df
*! handling of check 8, and check 10. It does NOT compute the statistic
*! (step 3); status says whether the caller may.
*!
*! Inputs: par0 and theta0 H0's table and estimate; T0, T1 the two fits'
*! T = c F (T_raw); c = N - G; nrtol the fits' Newton-decrement tolerance;
*! df0 H0's e(df_T3) (equal to H1's when d = 0); cc and K the shared
*! e(canon_configural) and e(K) (read only when exempt); inv0, inv1 the levels.
*!
*! Returns (status, d, Td, eps):
*!   status 1  d >= 1, T0 - T1 > eps: the caller computes the statistic
*!   status 2  d >= 1, |T0 - T1| <= eps: Td = 0, a restriction that does not
*!             bind to the fits' precision (what that means for T_d3 is the
*!             front end's, open)
*!   status 0  d = 0: no test (spec check 8, *Zero-df generated pairs*,
*!             option (c)); Td = 0 when |T0 - T1| <= eps, else T0 - T1
*! and aborts rc 3498 where a check refuses. note (optional output) holds the
*! approved d = 0 text, lines ended by char(10), else "". The routine prints
*! only aborts; the front end prints note verbatim (decided 2026-10-13, user,
*! option (ii)).
*!
*! CHECK 11 (spec, messages approved 2026-10-13): the generated-level backstop
*! sentence for an exempt pair, the user-written comparison sentence otherwise.
*!
*! d = 0. Only an exempt pair can reach it (check 8 refuses d < 1 otherwise; a
*! non-exempt d = 0 is refused here too, a backstop). The verified case is the
*! thresholds-metric pair with dwlsmv_m11_fit_g_i2(cc, K) (Result I2; option
*! (a) for non-adjacent pairs, 2026-10-13): within the band, the verified
*! wording; outside, the I2 abort. Any other d = 0 pair: the equal-df wording,
*! and outside the band the information line, never an abort. eps is check
*! 10's, c nrtol; no new tolerance.
*!
*! d >= 1: dwlsmv_dtband, with (exempt, par0, theta0) for the negative-phi
*! branch.
*!
*! Load order: dwlsmv_types, dwlsmv_dtboundary, dwlsmv_dtband, dwlsmv_canon_lines,
*! dwlsmv_m11_fit_g.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_dtstep(struct dwlsmv_parmap scalar par0,
                             real colvector theta0,
                             real scalar d, real scalar exempt,
                             real scalar T0, real scalar T1,
                             real scalar c, real scalar nrtol,
                             real scalar df0,
                             string scalar cc, real matrix K,
                             string scalar inv0, string scalar inv1,
                             string scalar name0, string scalar name1,
                             | string scalar note)
{
    real rowvector b
    real scalar Td, eps, ver, i2
    string scalar nt, dfs

    if (args() == 16) {
        note = ""
    }

    // check 11 on H0, the caller's consequence by whether the pair is generated
    if (exempt) {
        dwlsmv_dtboundary(par0, name0,
            "The level generator introduced this fixed value: the configural specification passed this check, so the forcing value comes from the generated level, not from the model you wrote. This is an internal error; the comparison is refused.")
    }
    else {
        dwlsmv_dtboundary(par0, name0,
            "This restriction lies on the boundary of the parameter space, where the chi-squared reference of the difference test does not hold, so the comparison is refused.")
    }

    if (d >= 1) {
        if (exempt) {
            b = dwlsmv_dtband(T0, T1, c, nrtol, name0, name1, exempt, par0, theta0)
        }
        else {
            b = dwlsmv_dtband(T0, T1, c, nrtol, name0, name1)
        }
        if (b[1] == 1) {
            return((1, d, b[2], b[3]))
        }
        return((2, d, 0, b[3]))
    }

    if (d != 0 | !exempt) {
        errprintf("dwlsmv (%s, %s): d = %g reached the step without the generated-chain exemption.\n",
                  name0, name1, d)
        _error("dwlsmv: a d < 1 pair was not refused by check 8")
    }
    if (missing(T0) | missing(T1) | missing(c) | missing(nrtol)) {
        errprintf("dwlsmv check 10 (%s, %s): T0, T1, c or nrtol is missing.\n", name0, name1)
        _error("dwlsmv check 10: T0, T1, c or nrtol is missing")
    }

    // d = 0: no test (spec check 8, zero-df generated pairs)
    eps = c * nrtol
    Td  = T0 - T1
    ver = 0
    if ((inv0 == "metric" & inv1 == "thresholds") | (inv0 == "thresholds" & inv1 == "metric")) {
        i2  = dwlsmv_m11_fit_g_i2(cc, K)
        ver = i2
    }
    dfs = strofreal(df0, "%21.0g")
    nt  = "No test for " + name0 + " against " + name1 +
          ": the two models have the same degrees of freedom (d = 0)." + char(10)

    if (abs(Td) <= eps) {
        Td = 0
    }
    else if (ver) {
        errprintf("dwlsmv check 8 (%s, %s): for this model, invariance(metric) is mathematically the\n",
                  name0, name1)
        errprintf("same model as invariance(thresholds) (every factor has at most one indicator with\n")
        errprintf("three or more categories), so both fits should agree to the precision they are\n")
        errprintf("computed to. They differ by %g (precision %g). The likeliest cause is that one of\n",
                  abs(Td), eps)
        errprintf("the two fits stopped short of its best solution, though another numerical problem\n")
        errprintf("cannot be ruled out. Refit both models; no comparison is reported.\n")
        _error("dwlsmv check 8: equivalent models (Result I2) fit differently")
    }

    if (ver) {
        nt = nt + "Note: invariance(metric) is equivalent to invariance(thresholds) for this model, so" + char(10) +
             "this step adds no testable restriction: both have df = " + dfs + ", and every factor" + char(10) +
             "has at most one indicator with three or more categories, the case in which the" + char(10) +
             "equivalence is established." + char(10)
    }
    else {
        nt = nt + "Note: invariance(" + inv0 + ") has the same degrees of freedom as invariance(" + inv1 + ")" + char(10) +
             "for this model (df = " + dfs + "). Equal degrees of freedom do not by themselves make two" + char(10) +
             "models equivalent, and equivalence is not established for this pattern." + char(10)
        if (Td != 0) {
            nt = nt + "The two models have the same degrees of freedom but do not fit the data equally" + char(10) +
                 "well (difference in fit " + strofreal(abs(Td), "%9.0g") + ", larger than the " +
                 strofreal(eps, "%9.0g") + " to which the fits are" + char(10) +
                 "computed), so equivalence of the two models could not be established. With no" + char(10) +
                 "difference in degrees of freedom there is nothing to test." + char(10)
        }
    }
    if (args() == 16) {
        note = nt
    }
    return((0, 0, Td, eps))
}

end
