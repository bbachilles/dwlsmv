*! dwlsmv -- Module 11: factors the canonical text leaves without a scale
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*!
*! Returns "F y" pairs separated by ";": each factor with no metric line, and
*! the indicator of its first lambda line (for the example in the remedy).
*! dwlsmv.ado calls it after a Module 3 abort, whose message is rethrown
*! unchanged (D8a), to add the remedy in the command's own terms: fix a
*! loading at 1 or use metric(unitvariance) (decided 2026-10-04). Empty when
*! every factor has a metric line, so other Module 3 aborts get no such line.
*! A factor with no lambda line never reaches here: Module 11 refuses it
*! (D3, second-order factors).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_m11_unanchored(string scalar canon)
{
    string rowvector L, fac, t
    string scalar out, first
    real scalar i, k, has

    L   = tokens(canon, ";")
    fac = tokens(substr(L[1], 8, .))
    out = ""
    for (k = 1; k <= cols(fac); k++) {
        has   = 0
        first = ""
        for (i = 1; i <= cols(L); i++) {
            t = tokens(L[i])
            if (cols(t) < 3) {
                continue
            }
            if (t[1] == "metric") {
                if (t[3] == fac[k]) {
                    has = 1
                }
            }
            else if (t[1] == "lambda") {
                if (t[2] == fac[k] & first == "") {
                    first = t[3]
                }
            }
        }
        if (!has) {
            out = out + (out == "" ? "" : ";") + fac[k] + " " + first
        }
    }
    return(out)
}

end
