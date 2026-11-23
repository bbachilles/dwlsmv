*! dwlsmv -- efacompare's largest residuals of one fit (Module 12, efacompare)
*!
*! Written by the main session 2026-11-22 (owner decision 5 of 2026-11-21; derivation
*! @sec-resid-diag; E-117, package 3c).
*!
*!     void dwlsmv_efacompare_top(real scalar m, string scalar topname)
*!
*! Reads the active e(Resid), e(Nresid) and e(indicators) of a fitted efa() row and writes the
*! Stata matrix topname, up to 3 x 6: (m, rank, i, j, residual, normalized), from
*! dwlsmv_resid_rank(e(Nresid), 3) (|z| descending, ties by smaller i then smaller j), with the
*! signed residual of e(Resid). Sets the caller's local efc_excl to the pairs whose normalized
*! residual is missing (one-empty-cell repair or held at the bound, O1/O2), "a-b" each, item
*! order. A missing or mis-sized e(Resid) or e(Nresid) is an internal error, rc 3498.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_efacompare_top(real scalar m, string scalar topname)
{
    real matrix E, Z, Rk, out
    string rowvector ind
    string scalar ex
    real scalar p, k, i, j

    E   = st_matrix("e(Resid)")
    Z   = st_matrix("e(Nresid)")
    ind = tokens(st_global("e(indicators)"))
    p   = cols(ind)
    if (rows(E) != p | cols(E) != p | rows(Z) != p | cols(Z) != p) {
        errprintf("efacompare: internal error: e(Resid) or e(Nresid) is not %g x %g.\n", p, p)
        _error("efacompare: internal error, residual matrices")
    }
    Rk  = dwlsmv_resid_rank(Z, 3)
    out = J(rows(Rk), 6, .)
    for (k = 1; k <= rows(Rk); k++) {
        i = Rk[k, 1]
        j = Rk[k, 2]
        out[k, .] = (m, k, i, j, E[j, i], Z[j, i])
    }
    if (rows(out)) {
        st_matrix(topname, out)
    }
    ex = ""
    for (i = 1; i < p; i++) {
        for (j = i + 1; j <= p; j++) {
            if (missing(Z[j, i])) {
                ex = ex + (ex == "" ? "" : " ") + ind[i] + "-" + ind[j]
            }
        }
    }
    st_local("efc_excl", ex)
}

end
