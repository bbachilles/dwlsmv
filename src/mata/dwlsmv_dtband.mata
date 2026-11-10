*! dwlsmv -- DIFFTEST validation check 10: the precision band on T0 - T1
*! (Module 12)
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9; executed by the main
*! session 2026-10-08: tests/unit/test_dtband.do 20/20, 11 mutants killed.
*! Changed 2026-10-09 (user): takes the multiplier c of F in T instead of N,
*! eps = c * nrtol (derivation @sec-multigroup-scaling, MG7; T = (N - G) F).
*!
*!     real rowvector dwlsmv_dtband(real scalar T0, real scalar T1,
*!                                  real scalar c, real scalar nrtol,
*!                                  string scalar name0, string scalar name1)
*!
*! Spec Module 12, Validation contract, check 10, with eps = c * nrtol. T0 and
*! T1 are the restricted and unrestricted fits' T = c F (spec 4.9), c = N - G
*! (N - 1 for one group), nrtol the Newton-decrement tolerance the fits were run
*! under (dwlsmv_minimize owns the number; it is passed in so there is one copy).
*! The caller computes c; any finite c > 0 is accepted.
*!
*! Returns (status, Td, eps):
*!     status 1  T0 - T1 > eps: proceed; Td = T0 - T1
*!     status 0  |T0 - T1| <= eps: a restriction that does not bind to the
*!               fits' precision; Td = 0 and the caller prints the note. A
*!               statement of precision, not truncation (spec).
*! and aborts rc 3498 when T0 - T1 < -eps: a premise (a fit at its global
*! minimum) has failed. The suspect is H1: H0 is nested in H1, so H0's optimum
*! is feasible in H1 and H1's minimum cannot exceed it (derived; the message
*! also says the nesting could be at fault, which check 8 has already tested).
*!
*! GENERATED PAIRS WITH A NEGATIVE H0 FACTOR VARIANCE (spec check 10, decided
*! 2026-10-13, option (c); built the same day, exemption gate step 2c, branch
*! here so check 10 keeps one abort site). Optional (exempt, par0, theta0):
*! when T0 - T1 < -eps AND exempt AND some free later-group psi_gff of H0 is
*! < 0 at theta0 (B = 0, so psi = phi; dwlsmv_dtb_negphi), the abort carries
*! the approved message (nesting premise failed at H0's estimate; not a
*! convergence finding; H1 not named as suspect). Otherwise, the message below.
*! Tests: tests/unit/test_dtstep.do.
*!
*! WHAT IT DOES NOT DO. No fit, no Delta, no Gamma, no statistic. It does not
*! decide what "reported as 0" means for T_d3 and its p-value (the caller;
*! open). It does not check that T0 and T1 are individually nonnegative. The
*! band holds only when the minimizer was handed weighting (a), whose F is
*! O(1) (MG7); that is the caller's obligation, not checked here.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The free later-group factor variances psi_gff of par below 0 at theta, as
// "name = value" lines (the "gN:" prefix replaced by "in group N"); empty if
// none, or if par has a beta row free or fixed nonzero (the branch needs
// B = 0, so psi is phi).
string colvector dwlsmv_dtb_negphi(struct dwlsmv_parmap scalar par,
                                   real colvector theta)
{
    string colvector out
    string scalar nm, gp
    real scalar i, v, k

    out = J(0, 1, "")
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 3) {
            if (par.free[i] | par.value[i] != 0) {
                return(J(0, 1, ""))
            }
        }
    }
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.free[i]) {
            if (par.row[i] == par.col[i] & rows(par.grp) >= i) {
                if (par.grp[i] > 1 & !missing(par.grp[i])) {
                    v = theta[par.fidx[i]]
                    if (v < 0) {
                        nm = par.name[i]
                        k  = strpos(nm, ":")
                        gp = strofreal(par.grp[i])
                        if (k > 0) {
                            nm = substr(nm, k + 1, .)
                        }
                        out = out \ (nm + " in group " + gp + " is estimated at " +
                                     strofreal(v, "%9.0g"))
                    }
                }
            }
        }
    }
    return(out)
}

real rowvector dwlsmv_dtband(real scalar T0, real scalar T1, real scalar c,
                             real scalar nrtol, string scalar name0,
                             string scalar name1,
                             | real scalar exempt,
                             struct dwlsmv_parmap par0,
                             real colvector theta0)
{
    real scalar eps, Td, k
    string colvector neg

    if (missing(T0) | missing(T1) | missing(c) | missing(nrtol)) {
        errprintf("dwlsmv check 10 (%s, %s): T0, T1, c or nrtol is missing.\n", name0, name1)
        _error("dwlsmv check 10: T0, T1, c or nrtol is missing")
    }
    if (c <= 0) {
        errprintf("dwlsmv check 10: c = %g; need the positive multiplier of F in T (N - G).\n", c)
        _error("dwlsmv check 10: c must be positive")
    }
    if (nrtol <= 0) {
        errprintf("dwlsmv check 10: nrtol = %g; need a positive tolerance.\n", nrtol)
        _error("dwlsmv check 10: nrtol must be positive")
    }

    eps = c * nrtol
    Td  = T0 - T1

    if (args() != 6 & args() != 9) {
        errprintf("dwlsmv check 10 (%s, %s): the generated-pair inputs are incomplete.\n", name0, name1)
        _error("dwlsmv check 10: incomplete generated-pair inputs")
    }

    // Generated pairs with a negative H0 factor variance (spec check 10,
    // decided 2026-10-13, option (c); message approved the same day): only
    // for a pair that passed check 8 through the exemption.
    if (Td < -eps & args() == 9) {
        if (exempt) {
            neg = dwlsmv_dtb_negphi(par0[1, 1], theta0)
            if (rows(neg) > 0) {
                errprintf("dwlsmv check 10: %s fits the data better than %s (by %g, beyond the %g to\n",
                          name0, name1, -Td, eps)
                errprintf("which the fits are computed), which a more restricted model should not do.\n")
                errprintf("%s's estimates include a factor variance below 0:\n", name0)
                for (k = 1; k <= rows(neg); k++) {
                    errprintf("    %s\n", neg[k])
                }
                errprintf("The invariance levels are nested only when every factor variance is positive,\n")
                errprintf("so at this estimate the comparison is not valid, whatever else may contribute to\n")
                errprintf("the difference. This is a property of %s's solution, not a sign that the fit\n", name0)
                errprintf("failed to converge. No statistic is reported.\n")
                exit(3498)
            }
        }
    }

    if (Td < -eps) {
        errprintf("dwlsmv check 10: T0 - T1 = %g for %s (H0) against %s (H1) is below -eps = %g,\n",
                  Td, name0, name1, -eps)
        errprintf("eps = c nrtol with c = %g. Since %s is nested in %s, H1 cannot fit\n",
                  c, name0, name1)
        errprintf("worse than H0 beyond the optimizer's precision: %s (H1) did not reach its\n", name1)
        errprintf("global minimum, or the models are not nested as tabled. No statistic is reported.\n")
        exit(3498)
    }

    if (abs(Td) <= eps) {
        return((0, 0, eps))
    }
    return((1, Td, eps))
}

end
