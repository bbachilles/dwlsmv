*! dwlsmv -- Module 11: D2's injection for a multiple-group configural text
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_translate_groups.do.
*!
*!   dwlsmv_inject_g(P)
*!
*! P is dwlsmv_translate_g(0, ...)'s text from the forcenoanchor pass. D2's
*! rule (dwlsmv_inject) is applied to group 1's section: every factor the user
*! did not normalise gets var(F@1), or var(e.F@1) if endogenous. Written
*! without a group prefix, sem applies it in every group. A loading the user
*! fixed at 1 counts as a normalisation here, as at G = 1, and is then refused
*! on the final pass (dwlsmv_tr_gcheck, decided 2026-10-12).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_inject_g(string colvector P)
{
    string colvector S
    string rowvector t
    real scalar a, b, i

    a = 0
    b = rows(P)
    for (i = 1; i <= rows(P); i++) {
        if (P[i] == "group 1") {
            a = i + 1
        }
        else if (P[i] == "group 2") {
            b = i - 1
        }
    }
    if (a == 0) {
        _error("dwlsmv_inject_g: no group 1 section")
    }
    // a loading fixed at 0 is how e(Cns) records a path written for another
    // group only (probe case 5.4); it is an absence, not a normalisation
    // (found 2026-10-12: counted as one, it left F to sem's own marker and
    // the fixed-at-1 refusal named a loading the user never wrote)
    S = P[a :: a + 1]
    for (i = a + 2; i <= b; i++) {
        t = tokens(P[i])
        if (t[1] == "lambda" & cols(t) == 5) {
            if (t[4] == "fixed" & strtoreal(t[5]) == 0) {
                continue
            }
        }
        S = S \ P[i]
    }
    return(dwlsmv_inject(S))
}

end
