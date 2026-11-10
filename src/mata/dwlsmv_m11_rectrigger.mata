*! dwlsmv -- Module 11: D11's recovery trigger (dwlsmv_m11_rectrigger)
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*!
*!     string scalar dwlsmv_m11_rectrigger(struct dwlsmv_fit scalar fit,
*!                                         real matrix A, real matrix Bd)
*!
*! Spec D11, Trigger. A is dwlsmv_admissible() at fit.theta, Bd is
*! dwlsmv_breaddiag(fit.Delta, Gamma): the caller's own values, computed once
*! and reused for the notes. Returns "" when the first marker fit is not
*! recovered, otherwise what fired, as a phrase completing "The marker-metric
*! fit ...", the causes joined as a list ("a", "a and b", "a, b and c"):
*!   converged != 1          did not converge
*!   any A[., 1] == 1        was not admissible
*!   Bd[3, 1] == 1           could not determine all parameters reliably (the
*!                           information in the fit was insufficient): the
*!                           bread's numerical rank (wording 2026-10-07, user)
*! Only a flag of exactly 1 fires: a missing admissibility flag (a condition
*! that could not be evaluated) does not, and only bread row 3 (numerical
*! rank) counts, not row 2's zero column (both decided by the user
*! 2026-10-07).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_m11_rectrigger(struct dwlsmv_fit scalar fit, real matrix A, real matrix Bd)
{
    string colvector why
    real scalar i, adm

    why = J(0, 1, "")
    if (fit.converged != 1) {
        why = why \ "did not converge"
    }
    adm = 0
    for (i = 1; i <= rows(A); i++) {
        if (A[i, 1] == 1) {
            adm = 1
        }
    }
    if (adm) {
        why = why \ "was not admissible"
    }
    if (Bd[3, 1] == 1) {
        why = why \ "could not determine all parameters reliably (the information in the fit was insufficient)"
    }
    if (rows(why) == 0) {
        return("")
    }
    if (rows(why) == 1) {
        return(why[1])
    }
    if (rows(why) == 2) {
        return(why[1] + " and " + why[2])
    }
    return(why[1] + ", " + why[2] + " and " + why[3])
}

end
