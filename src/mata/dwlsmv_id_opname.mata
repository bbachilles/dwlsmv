*! dwlsmv -- Module 3: a tie operand written back as text, for messages
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_identify_groups.do.
*!
*! op is one row of dwlsmv_id_tie()'s (ptype, group, row, col). Returns the
*! operand in the canonical spelling, e.g. "lambda 2 F1 y1" or "tau 1 y3 2",
*! with "?" for a name that was not declared. Spec 4.8 G7 asks the O3 message
*! to name both operands, every key position included (found 2026-10-10: the
*! first version named only the indicators, so a factor mismatch went unnamed).
*! Indices are range-checked with plain if statements: whether Mata's ?:
*! evaluates both branches is not relied on.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_id_opname(real rowvector op,
                               string colvector fac,
                               string colvector ind)
{
    string scalar g, fr, fc, ir
    real scalar pt, r, c

    pt = op[1]
    g  = strofreal(op[2])
    r  = op[3]
    c  = op[4]
    fr = "?"
    fc = "?"
    ir = "?"
    if (r >= 1 & r <= rows(fac)) {
        fr = fac[r]
    }
    if (c >= 1 & c <= rows(fac)) {
        fc = fac[c]
    }
    if (r >= 1 & r <= rows(ind)) {
        ir = ind[r]
    }
    if (pt == 1) {
        return("lambda " + g + " " + fc + " " + ir)
    }
    if (pt == 4) {
        return("tau " + g + " " + ir + " " + strofreal(c))
    }
    if (pt == 2) {
        return("psi " + g + " " + fr + " " + fc)
    }
    if (pt == 3) {
        return("beta " + g + " " + fr + " " + fc)
    }
    if (pt == 7) {
        return("alpha " + g + " " + fr)
    }
    if (pt == 5) {
        return("nu " + g + " " + ir)
    }
    return("omega " + g + " " + ir)
}

end
