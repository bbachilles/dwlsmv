*! dwlsmv -- Module 11: rotate the fitted EFA block and build its e() (0.3.0, Unit 6)
*!
*! Executed under StataMP 18 on 2026-10-29: tests/unit/test_efa_surface.do.
*!
*!     void dwlsmv_m11_efa_rotate(string scalar bn, string scalar canon,
*!                                string scalar epsilon, string scalar rstarts,
*!                                string scalar newbn, string scalar lamn,
*!                                string scalar phin, string scalar efalist)
*!
*! bn is the Stata matrix dwlsmv_m11_fit() filled from the echelon fit
*! (dwlsmv_m11_efa_canon()'s text, canon). Reads A (p x m, Phi = I) from bn's
*! loadings, rotates it with dwlsmv_geomin(), and fills:
*!
*!     newbn  e(b): rotated loadings y:F (bn's order), the factor correlations
*!            /cov(F,G) (bn's names), then bn's thresholds /y:cutk unchanged.
*!            No variances (fixed at 1) and no e(V): rotated SEs are not in
*!            0.3.0 (spec section 5 Module 12, the user's choice 2026-10-29).
*!     lamn   e(Lambda), p x m, rows the indicators, columns the factors
*!     phin   e(Phi), m x m
*!
*! and the locals rot_starts, rot_converged, rot_minima, rot_atbest, rot_Q,
*! rot_gpnorm, rot_eps, rot_gap (dwlsmv_geomin()'s info, in order). Rotated
*! column k carries the k-th factor name of the factor line: the rotation
*! defines columns only up to sign and order, which dwlsmv_geomin() fixes.
*!
*! epsilon "" means Mplus's rule (.0001 / .001 / .01 at m = 2 / 3 / >= 4);
*! rstarts "" means 30. dwlsmv.ado has already checked both.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_efa_rotate(string scalar bn, string scalar canon,
                           string scalar epsilon, string scalar rstarts,
                           string scalar newbn, string scalar lamn,
                           string scalar phin, string scalar efalist)
{
    string rowvector lines, fac, ind, lab
    string matrix stp, nst
    real rowvector b, info, nb
    real matrix A, L, Phi
    real scalar i, j, k, p, m, eps, ns, n

    lines = select(tokens(canon, ";"), tokens(canon, ";") :!= ";")
    fac   = tokens(substr(lines[1], 8, .))
    ind   = tokens(substr(lines[2], 5, .))
    lab   = tokens(efalist)
    p     = cols(ind)
    m     = cols(fac)
    b     = st_matrix(bn)
    stp   = st_matrixcolstripe(bn)
    // which stage an abort came from, for dwlsmv.ado: only an abort inside
    // dwlsmv_geomin() (its arguments are already valid) is non-convergence
    st_local("rot_stage", "read")

    A = J(p, m, .)
    for (i = 1; i <= cols(b); i++) {
        j = selectindex(ind :== stp[i, 1])
        k = selectindex(fac :== stp[i, 2])
        if (length(j) == 1 & length(k) == 1) {
            A[j, k] = b[i]
        }
    }
    if (missing(A)) {
        _error("dwlsmv: internal error, a loading of the EFA block is not in e(b)")
    }

    // set before rotating, so dwlsmv.ado can say it if the rotation fails:
    // an indicator whose unrotated communality exceeds 1 (negative Theta)
    st_local("rot_inadm", strofreal(any(rowsum(A:^2) :> 1)))

    eps = (epsilon == "" ? (m == 2 ? 1e-4 : (m == 3 ? 1e-3 : 1e-2)) : strtoreal(epsilon))
    ns  = (rstarts == "" ? 30 : strtoreal(rstarts))
    st_local("rot_stage", "rotate")
    L   = dwlsmv_geomin(A, eps, ns, Phi, info)
    st_local("rot_stage", "post")

    n   = 0
    nb  = J(1, cols(b), .)
    nst = J(cols(b), 2, "")
    for (i = 1; i <= cols(b); i++) {
        j = selectindex(ind :== stp[i, 1])
        k = selectindex(lab :== stp[i, 2])
        if (length(j) == 1 & length(k) == 1) {
            n++
            nb[n]     = L[j, k]
            nst[n, .] = stp[i, .]
        }
    }
    for (i = 1; i <= cols(b); i++) {
        if (substr(stp[i, 2], 1, 4) == "cov(") {
            j = selectindex(lab :== substr(stp[i, 2], 5, strpos(stp[i, 2], ",") - 5))
            k = selectindex(lab :== substr(stp[i, 2], strpos(stp[i, 2], ",") + 1,
                                           strlen(stp[i, 2]) - strpos(stp[i, 2], ",") - 1))
            n++
            nb[n]     = Phi[j, k]
            nst[n, .] = stp[i, .]
        }
    }
    for (i = 1; i <= cols(b); i++) {
        if (substr(stp[i, 1], 1, 1) == "/" & substr(stp[i, 2], 1, 3) == "cut") {
            n++
            nb[n]     = b[i]
            nst[n, .] = stp[i, .]
        }
    }
    st_matrix(newbn, nb[1, 1..n])
    st_matrixcolstripe(newbn, nst[1..n, .])
    st_matrix(lamn, L)
    st_matrixrowstripe(lamn, (J(p, 1, ""), ind'))
    st_matrixcolstripe(lamn, (J(m, 1, ""), lab'))
    st_matrix(phin, Phi)
    st_matrixrowstripe(phin, (J(m, 1, ""), lab'))
    st_matrixcolstripe(phin, (J(m, 1, ""), lab'))

    st_local("rot_starts",    strofreal(info[1], "%21.0g"))
    st_local("rot_converged", strofreal(info[2], "%21.0g"))
    st_local("rot_minima",    strofreal(info[3], "%21.0g"))
    st_local("rot_atbest",    strofreal(info[4], "%21.0g"))
    st_local("rot_Q",         strofreal(info[5], "%21.0g"))
    st_local("rot_gpnorm",    strofreal(info[6], "%21.0g"))
    st_local("rot_eps",       strofreal(info[7], "%21.0g"))
    st_local("rot_gap",       strofreal(info[8], "%21.0g"))
}

end
