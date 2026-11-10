*! dwlsmv -- how many factors a parameter table describes (Modules 3, 5, 6)
*!
*! Executed under StataMP 18 on 2026-09-22. Exercised by tests/unit/test_identify.do, test_model.do and test_delta.do, which build Lambda, Psi and B through it on one- and two-factor designs.
*!
*! Lambda is p x nfac and both Psi and B are nfac x nfac, so Modules 5 and 6
*! need the factor count before they can allocate anything. It is deliberately
*! a function and NOT a field on struct dwlsmv_parmap, for the same reason
*! dwlsmv_koff() is a function: a derived quantity stored in a struct is a
*! second copy that can disagree with the first. Same asymmetry as there, too
*! -- map.p IS a field, because the map is built from K and every module needs
*! it, whereas nfac is read by three routines that all already hold 'par'.
*!
*! WHY DERIVING IT IS SAFE, which is not obvious and is the whole argument for
*! not storing it. The factor count comes from the canonical form's 'factor'
*! declaration, and a declared factor that appeared in no parameter row at all
*! would be invisible here. It cannot happen: dwls_identify() requires every
*! declared factor to carry exactly one metric constraint, and BOTH forms of
*! that constraint emit a row indexed by the factor -- 'marker F y' a fixed
*! lambda row in column F, 'unitvariance F' a fixed psi row at (F,F). So every
*! factor that survives identification is named by at least one row, and the
*! maximum below equals the declared count. If that invariant is ever relaxed,
*! this function is the thing that breaks, silently, by returning a Lambda one
*! column too narrow -- so relax it by storing the count, not by patching here.
*!
*! ptype 1 is lambda, whose COLUMN is the factor (Lambda is p x nfac, so its
*! row is the indicator). ptype 2 (psi) and 3 (beta) are factor-by-factor, so
*! both indices count. ptype 4 (tau) names a variable and a threshold and is
*! skipped -- reading its row as a factor index is the obvious bug here.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_nfac(struct dwlsmv_parmap scalar par)
{
    real scalar i, n

    n = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1) {
            n = max((n, par.col[i]))
        }
        else if (par.ptype[i] == 2 | par.ptype[i] == 3) {
            n = max((n, par.row[i], par.col[i]))
        }
    }

    return(n)
}

end
