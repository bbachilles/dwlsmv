*! dwlsmv -- Module 3: the operands of one 'tie' line (spec 4.8 G2)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_identify_groups.do.
*!
*! Tokenizes 'tie <elem> <elem> ...' into an n x 4 matrix (ptype, group, row,
*! col), one row per operand, in the element's directive order:
*!
*!     lambda <g> <F> <y>   ptype 1, row y, col F    (as a lambda row)
*!     psi    <g> <F> <F>   ptype 2, lower triangle  (both spellings one element)
*!     beta   <g> <to> <fr> ptype 3, row to, col from
*!     tau    <g> <y> <k>   ptype 4, row y, col k
*!     nu     <g> <y>       ptype 5, row y, col 1
*!     omega  <g> <y>       ptype 6, row y, col 1
*!     alpha  <g> <F>       ptype 7, row F, col 1
*!
*! All seven spellings are parsed so that a known but unadmitted type gets
*! its own refusal (O4) rather than "unrecognized". Refused here: an unknown
*! type, an incomplete operand, fewer than two operands. A group token that is
*! not an integer is "incomplete"; an integer group outside 1..G and a name
*! that is not declared (row or col 0) are left to the caller, which judges
*! every tie in the precedence order decided 2026-10-10.
*!
*! Load order: dwlsmv_types, dwlsmv_canon_find.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_id_tie(string rowvector w,
                          string colvector fac,
                          string colvector ind)
{
    real matrix ops
    string scalar t
    real scalar k, need, ptc, g, r, c, a, b

    ops = J(0, 4, .)
    k   = 2
    while (k <= cols(w)) {
        t = w[k]
        ptc = (t == "lambda") + 2 * (t == "psi") + 3 * (t == "beta") + 4 * (t == "tau") +
              5 * (t == "nu") + 6 * (t == "omega") + 7 * (t == "alpha")
        need = (ptc <= 4 ? 3 : 2)
        if (ptc == 0) {
            errprintf("dwls_identify: a 'tie' line has '%s' where an element type is expected.\n", t)
            errprintf("Element types: lambda, beta, psi, tau, nu, omega, alpha.\n")
            _error("dwls_identify: unrecognized element type in a tie")
        }
        if (k + need > cols(w)) {
            errprintf("dwls_identify: the tie operand starting '%s' has too few tokens.\n", t)
            _error("dwls_identify: an element in a tie is incomplete")
        }
        g = strtoreal(w[k + 1])
        if (missing(g) | g != floor(g)) {
            errprintf("dwls_identify: the tie operand '%s %s ...' has no integer group number.\n",
                      t, w[k + 1])
            _error("dwls_identify: an element in a tie is incomplete")
        }
        if (ptc == 1) {
            r = dwlsmv_canon_find(ind, w[k + 3])
            c = dwlsmv_canon_find(fac, w[k + 2])
        }
        else if (ptc == 2) {
            a = dwlsmv_canon_find(fac, w[k + 2])
            b = dwlsmv_canon_find(fac, w[k + 3])
            r = max((a, b))
            c = (min((a, b)) == 0 ? 0 : min((a, b)))
            if (a == 0 | b == 0) {
                r = 0
            }
        }
        else if (ptc == 3) {
            r = dwlsmv_canon_find(fac, w[k + 2])
            c = dwlsmv_canon_find(fac, w[k + 3])
        }
        else if (ptc == 4) {
            r = dwlsmv_canon_find(ind, w[k + 2])
            c = strtoreal(w[k + 3])
            if (missing(c) | c != floor(c)) {
                errprintf("dwls_identify: the tie operand 'tau %s %s %s' has no integer threshold number.\n",
                          w[k + 1], w[k + 2], w[k + 3])
                _error("dwls_identify: an element in a tie is incomplete")
            }
        }
        else if (ptc == 7) {
            r = dwlsmv_canon_find(fac, w[k + 2])
            c = 1
        }
        else {
            r = dwlsmv_canon_find(ind, w[k + 2])
            c = 1
        }
        ops = ops \ (ptc, g, r, c)
        k   = k + need + 1
    }

    if (rows(ops) < 2) {
        errprintf("dwls_identify: a 'tie' line names %f element(s); a tie joins at least two.\n",
                  rows(ops))
        _error("dwls_identify: a tie needs at least two elements")
    }
    return(ops)
}

end
