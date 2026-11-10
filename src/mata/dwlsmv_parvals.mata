*! dwlsmv -- the current value of every parameter, free or fixed (Modules 5, 6)
*!
*! Executed under StataMP 18 on 2026-09-22. Exercised by every Module 5 and 6 call, and by test_delta.do's hand-built table, whose appended row has no theta position until fidx is read.
*!
*! Returns an nall x 1 column: theta's entry where the row is free, the table's
*! own 'value' where it is fixed. Every routine below the Stage 2 boundary
*! works from this vector rather than from theta directly, and that is what
*! makes Module 6's design note cheap to honour -- spec Section 5.1 requires
*! dwls_delta() to be able to produce a column for a parameter that is
*! currently FIXED (Sorbom modification indices, Module 12), and a fixed
*! parameter has no theta position to perturb. Perturbing an entry of this
*! vector works for either kind.
*!
*! It also localizes the one place theta's length is checked. A theta of the
*! wrong length would otherwise surface as an out-of-range subscript somewhere
*! inside the matrix fill, far from the caller that got it wrong.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_parvals(real colvector theta,
                              struct dwlsmv_parmap scalar par)
{
    real colvector v
    real scalar i

    if (rows(theta) != par.npar) {
        errprintf("dwlsmv_parvals: theta has %f row(s); the parameter table has %f free parameter(s)\n",
                  rows(theta), par.npar)
        _error("dwlsmv_parvals: theta is not npar rows long")
    }

    v = par.value
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            v[i] = theta[par.fidx[i]]
        }
    }

    return(v)
}

end
