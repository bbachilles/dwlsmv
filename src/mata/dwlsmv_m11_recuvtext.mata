*! dwlsmv -- Module 11: the unit-variance canonical text of a marker model
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*!
*!     string scalar dwlsmv_m11_recuvtext(string scalar canon)
*!
*! Spec D11 step 1. 'canon' is Module 11's marker canonical text, lines
*! separated by ";" (r(canon) of _dwlsmv_translate). Returns the same model
*! under unit variance, ";"-separated, for dwls_identify().
*!
*! Three things per
*! marker factor F with marker y:
*!   1. 'metric marker F y'  ->  'metric unitvariance F';
*!   2. add 'lambda F y free': the marker loading was created by the metric
*!      line alone; without this line y's loading on F would not exist in
*!      the unit-variance table (dwls_identify: no directive, no element);
*!   3. drop 'psi F F free': Module 3 rejects a metric directive colliding
*!      with an explicit one (dwls_identify.mata, metric pass), and under the
*!      marker metric Module 11 always writes that line.
*! This is mr_uv_table() (derivation chunk mr-impl) in text form: the marker
*! loading becomes free, psi_ff becomes fixed at 1. Every other line is
*! carried unchanged. Called only after dwlsmv_m11_receligible() passed, so
*! psi_ff is free on every marker factor; a fixed one is an internal error.
*!
*! & does not short-circuit in Mata, so token tests are nested on cols().
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_m11_recuvtext(string scalar canon)
{
    string rowvector L, t, mk
    string colvector out
    real scalar i, keep

    L  = tokens(canon, ";")
    mk = J(1, 0, "")
    for (i = 1; i <= cols(L); i++) {
        if (L[i] == ";") {
            continue
        }
        t = tokens(L[i])
        if (cols(t) == 4) {
            if (t[1] == "metric") {
                if (t[2] == "marker") {
                    mk = mk, t[3]
                }
            }
        }
    }

    out = J(0, 1, "")
    for (i = 1; i <= cols(L); i++) {
        if (L[i] == ";") {
            continue
        }
        t = tokens(L[i])
        if (cols(t) == 0) {
            continue
        }
        keep = 1
        if (t[1] == "metric") {
            if (cols(t) == 4) {
                if (t[2] == "marker") {
                    out  = out \ ("metric unitvariance " + t[3])
                    out  = out \ ("lambda " + t[3] + " " + t[4] + " free")
                    keep = 0
                }
            }
        }
        else if (t[1] == "psi") {
            if (cols(t) >= 4) {
                if (t[2] == t[3]) {
                    if (anyof(mk, t[2])) {
                        if (t[4] != "free") {
                            _error("dwlsmv: internal error, a marker factor has a fixed variance")
                        }
                        keep = 0
                    }
                }
            }
        }
        if (keep) {
            out = out \ strtrim(L[i])
        }
    }
    return(invtokens(out', ";"))
}

end
