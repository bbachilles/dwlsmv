*! dwlsmv -- Module 11: write D11's transformed starts into the marker text
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*!
*!     string scalar dwlsmv_m11_recstarts(string scalar canon,
*!                                        struct dwlsmv_parmap scalar parm,
*!                                        real colvector thetam)
*!
*! Spec D11 step 4. 'canon' is the marker canonical text, ";"-separated;
*! parm its parameter table (the first fit's); thetam the transformed point
*! in parm's theta order (dwlsmv_m11_rectransform()). Every 'lambda', 'beta'
*! or 'psi' line whose status is 'free' becomes '<d> <a> <b> free <start>'
*! (spec 4.8's 'free [<start>]'; an existing start is replaced). Fixed lines,
*! 'factor', 'ind' and 'metric' lines are carried unchanged, so the result is
*! the same model and dwls_identify() builds the same table, starts aside.
*!
*! The line's element is found by its key (ptype, row, col), from the
*! 'factor' and 'ind' declarations, exactly as dwls_identify() reads it: a psi
*! line in either spelling is the lower-triangle element.
*!
*! THRESHOLDS ARE NOT WRITTEN: the grammar has no tau directive (dwls_identify
*! refuses one), and Module 3 starts every tau at 0. Their transformed values
*! (unchanged by the map) are therefore not carried to the final fit.
*!
*! PRECISION (spec D11, decided by the user 2026-10-07): the first of
*! strofreal(x), %21.0g and %21x that reads back as the identical double.
*! Measured 2026-10-07 on StataMP 18, 40,000 random doubles per decade:
*! %21.0g fails for 22% near 1e-3 and 92% at 1e-4 and below (e-notation,
*! about 15 significant digits); %21x read back on every value. The read-back
*! is still asserted, so a start that would not read back identically can
*! never reach dwls_fit().
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_m11_recstarts(string scalar canon,
                                   struct dwlsmv_parmap scalar parm,
                                   real colvector thetam)
{
    string rowvector L, t
    string colvector fac, ind, out
    string scalar s
    real scalar i, k, j, ptc, r, cc, a, b, x

    L   = tokens(canon, ";")
    fac = J(0, 1, "")
    ind = J(0, 1, "")
    for (i = 1; i <= cols(L); i++) {
        if (L[i] == ";") {
            continue
        }
        t = tokens(L[i])
        if (cols(t) >= 2) {
            if (t[1] == "factor") {
                fac = t[|2 \ cols(t)|]'
            }
            else if (t[1] == "ind") {
                ind = t[|2 \ cols(t)|]'
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
        ptc = 0
        if (cols(t) >= 4) {
            if (t[4] == "free") {
                if (t[1] == "lambda") {
                    ptc = 1
                    r   = dwlsmv_canon_find(ind, t[3])
                    cc  = dwlsmv_canon_find(fac, t[2])
                }
                else if (t[1] == "beta") {
                    ptc = 3
                    r   = dwlsmv_canon_find(fac, t[2])
                    cc  = dwlsmv_canon_find(fac, t[3])
                }
                else if (t[1] == "psi") {
                    ptc = 2
                    a   = dwlsmv_canon_find(fac, t[2])
                    b   = dwlsmv_canon_find(fac, t[3])
                    r   = max((a, b))
                    cc  = min((a, b))
                }
            }
        }
        if (ptc == 0) {
            out = out \ strtrim(L[i])
            continue
        }
        k = 0
        for (j = 1; j <= parm.nall; j++) {
            if (parm.ptype[j] == ptc) {
                if (parm.row[j] == r) {
                    if (parm.col[j] == cc) {
                        k = j
                    }
                }
            }
        }
        if (k == 0) {
            errprintf("dwlsmv_m11_recstarts: no table row for the line %s\n", L[i])
            _error("dwlsmv: internal error, a free line has no row in the marker table")
        }
        if (!parm.free[k]) {
            errprintf("dwlsmv_m11_recstarts: %s is fixed in the table, free in the text\n", parm.name[k])
            _error("dwlsmv: internal error, the marker text and table disagree")
        }
        x = thetam[parm.fidx[k]]
        s = strtrim(strofreal(x))
        if (strtoreal(s) != x) {
            s = strtrim(strofreal(x, "%21.0g"))
        }
        if (strtoreal(s) != x) {
            s = strtrim(strofreal(x, "%21x"))
        }
        if (missing(x)) {
            _error("dwlsmv: internal error, a transformed start is missing")
        }
        if (strtoreal(s) != x) {
            errprintf("dwlsmv_m11_recstarts: start %s for %s does not read back exactly\n",
                      s, parm.name[k])
            _error("dwlsmv: internal error, a start does not round-trip through text")
        }
        out = out \ (t[1] + " " + t[2] + " " + t[3] + " free " + s)
    }
    return(invtokens(out', ";"))
}

end
