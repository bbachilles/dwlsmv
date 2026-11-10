*! dwlsmv -- fill Lambda, Psi and B from a parameter table (Modules 5, 6)
*!
*! Executed under StataMP 18 on 2026-09-22. The index conventions below are checked against the tier-2 oracle's own matrix-element keys in test_model.do and test_delta.do, which is what would catch a transposed lambda.
*!
*! Mata passes by reference, so the three matrices are arguments and the
*! routine returns void -- the same shape spec Section 5.1 uses for every
*! struct-filling routine.
*!
*! THE INDEX CONVENTIONS, which are the only content here and the only thing
*! that can be wrong:
*!
*!     ptype 1  lambda   row = INDICATOR, col = FACTOR    (Lambda is p x nfac)
*!     ptype 2  psi      row, col = FACTOR, FACTOR        (symmetric, filled
*!                                                         on both sides)
*!     ptype 3  beta     row = TO, col = FROM             (eta = B eta + zeta)
*!     ptype 4  tau      not a structural matrix at all; skipped
*!
*! The canonical form of spec Section 4.8 writes 'lambda <F> <y>' -- factor
*! first -- while Lambda is indicator by factor. dwls_identify() does that
*! transposition once, when it builds the table, so every reader here can take
*! row and col at face value. Reading a lambda row as (factor, indicator)
*! produces a non-conformable multiply on a square design and SILENTLY WRONG
*! numbers on any design where p happens to equal nfac, which is why the
*! convention is stated rather than implied.
*!
*! Anything the table does not name is zero. That is the canonical form's rule,
*! not a convenience: it is a literal serialization of the parameter table, so
*! an element with no directive has no row and is not a parameter. The one
*! exception is Psi's diagonal, which dwls_identify() requires every factor to
*! carry -- a zero disturbance variance is never what the user meant.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_matrices(real colvector vals,
                     struct dwlsmv_parmap scalar par,
                     real scalar p,
                     real scalar nfac,
                     real matrix Lam,
                     real matrix Psi,
                     real matrix B)
{
    real scalar i, r, c

    Lam = J(p,    nfac, 0)
    Psi = J(nfac, nfac, 0)
    B   = J(nfac, nfac, 0)

    for (i = 1; i <= par.nall; i++) {
        r = par.row[i]
        c = par.col[i]

        if (par.ptype[i] == 1) {
            Lam[r, c] = vals[i]
        }
        else if (par.ptype[i] == 2) {
            Psi[r, c] = vals[i]
            Psi[c, r] = vals[i]
        }
        else if (par.ptype[i] == 3) {
            B[r, c] = vals[i]
        }
    }
}

end
