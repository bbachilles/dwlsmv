*! dwlsmv -- the admissibility note of a conditional fit (one covariate), Module 11
*!
*! Executed 2026-11-14 by the main session through tests/unit/test_condx_marker.do; formerly marked not yet executed. Drafted by a read-only agent 2026-11-14. Exercised by tests/unit/test_condx_marker.do
*! (Part A at the Mata level, Part B through the command). Remove this banner only after the file has been
*! run and what broke has been fixed.
*!
*!     string scalar dwlsmv_cx_admissnote(real matrix A,
*!                                        string colvector ind,
*!                                        string colvector fac,
*!                                        string scalar metric,
*!                                        real scalar converged)
*!
*! A is dwlsmv_admissible()'s 6 x 4 result at the REPORTED theta, converged or not; ind and fac are the
*! indicator and factor names in declaration order (the covariate is in neither). Returns "" when no
*! condition holds. Otherwise it is dwlsmv_admissadvice()'s note, unchanged in its lines, with one edit that
*! only a conditional fit needs:
*!
*!   1. converged != 1: admissadvice ends with "Point estimates are reported. Standard errors and the WLSMV
*!      test ... should not be relied on". A conditional fit that did not converge reports no standard
*!      errors and no test, so that closing paragraph is cut and replaced by one sentence saying that the
*!      estimates are the last iterate.
*!   (2026-11-15: the former second edit, a paragraph saying that no recovery refit was attempted for a model
*!   with a covariate, is removed. D11's recovery now runs for a conditional marker fit and its own note
*!   (dwlsmv_m11_recover(), printed before this one) says what was tried or why not. The metric argument is
*!   kept so that the signature is unchanged; it is not read.)
*!
*! WHY ADMISSADVICE IS REUSED AND NOT REWRITTEN. It reads only A's cells, ind and fac: row 1's cell 3 indexes
*! ind (one entry per indicator, the rows of Lambda), row 2's cells 3 and 4 index fac, rows 3 and 4 name no
*! one. dwlsmv_admissible() builds Lambda, Psi and B from ptypes 1 to 4 and ignores gamma and kappa rows
*! (ptype 8, 9), so nothing in A is indexed by a covariate row.
*!
*! The covariate's gamma and kappa never enter a condition (derivation @sec-condx-full C6). No claim is made
*! here that an inadmissible conditional point would be admissible under the other metric.
*!
*! Newlines are char(10). No backtick-word-apostrophe sequence appears in any string.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_cx_admissnote(real matrix A,
                                   string colvector ind,
                                   string colvector fac,
                                   string scalar metric,
                                   real scalar converged)
{
    string scalar msg, nl
    real scalar k

    msg = dwlsmv_admissadvice(A, ind, fac)
    if (msg == "") {
        return("")
    }
    nl = char(10)

    if (converged != 1) {
        k = strpos(msg, "  Point estimates are reported.")
        if (k > 0) {
            msg = substr(msg, 1, k - 1)
        }
        msg = msg + "  The estimates are the last iterate of a fit that did not converge; no standard" + nl +
              "  errors and no test are reported for them." + nl
    }
    return(msg)
}

end
