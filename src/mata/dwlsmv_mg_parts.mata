*! dwlsmv -- one group's matrices, nu, omega, alpha and kappa (Modules 5, 6)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_mg_model.do.
*!
*! From group g's rows (dwlsmv_mg_sub) and their current values, fills
*! Lambda, Psi, B (dwlsmv_matrices, unchanged), A = (I - B)^-1, and the p x 1
*! nu and omega and nfac x 1 alpha. An absent row reads its type's reference
*! value, one rule for every consumer (spec 4.8 G4, decided 2026-10-09):
*! omega 1, nu 0, alpha 0. kappa = A alpha is the factor mean
*! (@sec-mg-sigma-delta notation).
*!
*! Load order: dwlsmv_types, dwlsmv_matrices.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_mg_parts(real colvector vals,
                     struct dwlsmv_parmap scalar sub,
                     real scalar p,
                     real scalar nfac,
                     real matrix Lam,
                     real matrix Psi,
                     real matrix B,
                     real matrix A,
                     real colvector nu,
                     real colvector om,
                     real colvector al,
                     real colvector kap)
{
    real scalar i

    dwlsmv_matrices(vals, sub, p, nfac, Lam, Psi, B)

    nu = J(p, 1, 0)
    om = J(p, 1, 1)
    al = J(nfac, 1, 0)
    for (i = 1; i <= sub.nall; i++) {
        if (sub.ptype[i] == 5) {
            nu[sub.row[i]] = vals[i]
        }
        else if (sub.ptype[i] == 6) {
            om[sub.row[i]] = vals[i]
        }
        else if (sub.ptype[i] == 7) {
            al[sub.row[i]] = vals[i]
        }
    }

    A = luinv(I(nfac) - B)
    if (missing(A)) {
        errprintf("dwlsmv_mg_parts: (I - B) is singular in one group, so (I-B)^-1 does not exist.\n")
        _error("dwlsmv_mg_parts: (I - B) is singular")
    }
    kap = A * al
}

end
