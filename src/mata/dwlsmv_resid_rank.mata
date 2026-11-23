*! dwlsmv -- residual diagnostics: the largest normalized residuals
*!
*! Drafted 2026-11-22 by a read-only agent; executed by the main session (StataMP 18, 2026-11-22), see dwlsmv_resid.mata.
*!
*!     real matrix dwlsmv_resid_rank(real matrix Z, real scalar k)
*!
*! Z: the p x p normalized-residual matrix of dwlsmv_resid(). Returns up to k rows
*! (i, j), i < j, of the pairs with a non-missing Z, ordered by |z| descending.
*! Ties are broken as fixtures/generate/99t_lavaan_srmr.R does, order(-abs(z), i, j):
*! the smaller i first, then the smaller j. Fewer than k rows if fewer pairs are
*! eligible; a 0 x 2 matrix if none. The diagonal is ignored. Missing Z (a repaired
*! pair, a pair at the bound) is never ranked (owner O1, O2).
*!
*! Aborts rc 3498: Z not square, k missing or < 1, Z missing at (i,j) but not at
*! (j,i) or the reverse (not a matrix dwlsmv_resid() returns).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_resid_rank(real matrix Z, real scalar k)
{
    real scalar p, i, j, n, m
    real colvector ii, jj, az, ord
    real matrix out

    if (rows(Z) != cols(Z)) {
        errprintf("dwlsmv_resid_rank: Z is not square (%g x %g).\n", rows(Z), cols(Z))
        _error("dwlsmv_resid_rank: Z is not square")
    }
    if (missing(k)) {
        errprintf("dwlsmv_resid_rank: k must be at least 1, it is missing.\n")
        _error("dwlsmv_resid_rank: k must be at least 1")
    }
    if (k < 1) {
        errprintf("dwlsmv_resid_rank: k must be at least 1, it is %g.\n", k)
        _error("dwlsmv_resid_rank: k must be at least 1")
    }
    p  = rows(Z)
    ii = J(0, 1, .)
    jj = J(0, 1, .)
    az = J(0, 1, .)
    for (i = 1; i <= p - 1; i++) {
        for (j = i + 1; j <= p; j++) {
            if (missing(Z[j, i]) != missing(Z[i, j])) {
                errprintf("dwlsmv_resid_rank: Z is not symmetric in its missing pattern at (%g, %g).\n", i, j)
                _error("dwlsmv_resid_rank: Z is not symmetric")
            }
            if (missing(Z[j, i]) == 0) {
                ii = ii \ i
                jj = jj \ j
                az = az \ abs(Z[j, i])
            }
        }
    }
    n = rows(ii)
    if (n == 0) {
        return(J(0, 2, .))
    }
    // ascending in -|z|, then i, then j: R's order(-abs(z), i, j)
    ord = order((-az, ii, jj), (1, 2, 3))
    m   = min((floor(k), n))
    out = (ii, jj)
    return(out[ord[|1 \ m|], .])
}

end
