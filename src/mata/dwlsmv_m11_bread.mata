*! dwlsmv -- Module 11: the bread's rows at the caller, missing-safe (dwlsmv_m11_bread)
*!
*! Executed 2026-11-15 through tests/unit/test_bread_missing.do (StataMP 18). Drafted by a read-only agent.
*!
*!     real matrix dwlsmv_m11_bread(real matrix Delta, real matrix Gamma)
*!
*! The caller-level rule (spec Module 10, the bread's rows): dwlsmv_breaddiag() keeps its strict
*! contract and aborts (rc 3498) on a Delta with a missing value. No fit is known to reach one: on a
*! reading of dwls_fit and dwls_delta, Delta is evaluated at a finite iterate that dwls_model() has
*! already accepted. The guard is defensive, so that such a Delta would be reported, not abort the
*! command. Module 11 asks the bread through this function: a Delta with any missing value returns the
*! "not evaluated" 4 x 4 that dwlsmv_diagadvice() documents for a caller with no bread,
*! (J(4, 2, .), J(4, 2, 0)), whose flags are all missing, so it fires no trigger term, prints no bread
*! note and posts e(bread_ratio) missing; any other Delta returns dwlsmv_breaddiag(Delta, Gamma)
*! unchanged, the same call. The shape is 4 x 4, never 0 x 0: dwlsmv_m11_fit reads Bd[4, 1].
*! Callers: dwlsmv_m11_estimate (both call sites; ordinary and conditional fits) and
*! dwlsmv_m11_post_g (G > 1, Gamma = diag of the group-weighted diagonals).
*!
*! Requires dwlsmv_breaddiag.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_m11_bread(real matrix Delta, real matrix Gamma)
{
    if (missing(Delta) > 0) {
        return((J(4, 2, .), J(4, 2, 0)))
    }
    return(dwlsmv_breaddiag(Delta, Gamma))
}

end
