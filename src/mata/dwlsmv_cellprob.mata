*! dwlsmv -- bivariate cell probabilities pi_rs (Module 2, spec Section 4.5.2)
*!
*! Executed under StataMP 18 on 2026-09-19. Verified against central
*! differences and against two structural identities that hold EXACTLY: the
*! cells sum to 1 and their rho-derivative sums to 0, both to 0.0 rather than
*! to a tolerance, because the rectangle differencing telescopes. Its
*! non-positive-cell abort was re-executed 2026-09-25 and now reaches the user
*! as rc 3498 with its message (tests/unit/test_polychoric.do, error paths).
*!
*!   pi_rs = Phi2(a_{r+1}, b_{s+1}) - Phi2(a_r, b_{s+1})
*!         - Phi2(a_{r+1}, b_s)     + Phi2(a_r, b_s)
*!
*! Rectangle differencing of the bivariate normal CDF over the cell's four
*! corners, with a_0 = b_0 = -infinity and a_{K1} = b_{K2} = +infinity.
*! Returns K1 x K2, where K1 = rows(a) + 1.
*!
*! MATA HAS NO INFINITY, which is why the corner grid is built by index rather
*! than by augmenting the threshold vector with +/- Inf the way the derivation
*! file's R does. In Mata `1/0` is missing, not Inf, and binormal() handed a
*! missing returns missing -- which would propagate into every cell in the
*! first and last row and column without raising anything. The boundary
*! corners are therefore branched explicitly:
*!
*!     Phi2(-inf, y) = 0        Phi2(x, -inf) = 0
*!     Phi2(+inf, y) = Phi(y)   Phi2(x, +inf) = Phi(x)      Phi2(inf,inf) = 1
*!
*! EACH CORNER IS COMPUTED ONCE. The naive implementation evaluates Phi2 four
*! times per cell, which is four times the bivariate-normal work and, worse,
*! computes the same shared corner by two different routes -- the cancellation
*! in the rectangle difference is then between numbers that are not bitwise
*! identical, and small cells lose accuracy for no reason. Build the (K1+1) x
*! (K2+1) corner grid, then difference it.
*!
*! Every pi_rs is strictly positive for finite thresholds and |rho| < 1 IN
*! EXACT ARITHMETIC ONLY. binormal() has an absolute, not a relative, floor:
*! on a banded table near rho = 1 the off-diagonal corner cells fall below
*! 1e-16 and come back as 0 or +/- an ulp of the corners, whichever route
*! computes them (derivation @sec-emptycells). dwlsmv_cellprob() aborts on a
*! non-positive cell, which is right for a caller that needs EVERY cell.
*! Callers that need only the observed cells -- the log-likelihood, its
*! polish, and dwlsmv_pairinfo() at the empirical weights -- use
*! dwlsmv_cellprob_raw(), the same arithmetic without the abort, and apply
*! positivity to the cells they actually use.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_cellprob(real colvector a, real colvector b, real scalar rho)
{
    real matrix P

    P = dwlsmv_cellprob_raw(a, b, rho)

    // The _error() text must stay under Stata's 100-character limit: beyond it
    // _error() raises rc 3300 and DISCARDS the message. This one was 105
    // characters until 2026-09-24, so every caller saw a bare 3300 (found
    // driving near-duplicate indicators through Stage 1). The detail goes
    // through errprintf() instead.
    if (min(P) <= 0) {
        errprintf("dwlsmv_cellprob: the model probability of a cell is %s at rho = %s.\n",
                  strofreal(min(P), "%10.3e"), strofreal(rho, "%12.9f"))
        errprintf("  Every cell is strictly positive for increasing finite thresholds\n")
        errprintf("  and |rho| < 1, so either two thresholds coincide, or rho is so\n")
        errprintf("  close to +/-1 that a cell underflows in floating point.\n")
        _error("dwlsmv_cellprob: a model cell probability is not positive")
    }

    return(P)
}

// The cells with no positivity check: an entry may be 0 or slightly negative
// where the true probability is below binormal()'s resolution. The caller
// decides which cells must be positive.
real matrix dwlsmv_cellprob_raw(real colvector a, real colvector b, real scalar rho)
{
    real matrix C, P
    real scalar K1, K2, r, s

    if (cols(a) != 1 | cols(b) != 1) {
        _error("dwlsmv_cellprob: thresholds must be column vectors")
    }
    if (rows(a) < 1 | rows(b) < 1) {
        _error("dwlsmv_cellprob: each variable needs at least one threshold")
    }
    if (rho <= -1 | rho >= 1 | rho >= .) {
        _error("dwlsmv_cellprob: rho must lie strictly inside (-1, 1)")
    }

    K1 = rows(a) + 1
    K2 = rows(b) + 1

    // Corner grid: row 1 is -inf, row K1+1 is +inf, row r otherwise a[r-1].
    C = J(K1 + 1, K2 + 1, 0)
    for (r = 1; r <= K1 + 1; r++) {
        for (s = 1; s <= K2 + 1; s++) {
            if (r == 1 | s == 1) {
                C[r, s] = 0
            }
            else if (r == K1 + 1 & s == K2 + 1) {
                C[r, s] = 1
            }
            else if (r == K1 + 1) {
                C[r, s] = normal(b[s - 1])
            }
            else if (s == K2 + 1) {
                C[r, s] = normal(a[r - 1])
            }
            else {
                C[r, s] = binormal(a[r - 1], b[s - 1], rho)
            }
        }
    }

    P = J(K1, K2, 0)
    for (r = 1; r <= K1; r++) {
        for (s = 1; s <= K2; s++) {
            P[r, s] = C[r + 1, s + 1] - C[r, s + 1] - C[r + 1, s] + C[r, s]
        }
    }

    return(P)
}

end
