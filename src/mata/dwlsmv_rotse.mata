*! dwlsmv -- the covariance of the rotated EFA solution, numerical core (release A, unit A5)
*!
*! Drafted 2026-11-17 by a read-only agent. Executed under StataMP 18 by the main session the
*! same day: test_rotse.do passes, and the R sign mutation is caught by checks 3, 4, 6, 7b and 10.
*! The chain entry point was renamed from dwlsmv_rotse_acov to dwlsmv_rotse (test_loading contract).
*!
*! Derivation: spec/derivations/dwlsmv-derivations.qmd @sec-rot-se (RS1 to RS7); the
*! R implementations rs_r, rs_R, rs_Gu, rs_GT, rs_J of chunk rs-impl are the source of
*! every formula below. Decisions: spec Module 12, Release A, A3.13 (production
*! computes J ACOV_u J', RS7 route 3, at the GPA stall point, A3.15), A3.2 (the rank
*! cut). Internal helpers below the locked boundary: no locked signature changes and
*! no struct field added; nothing here is reached from the public API.
*!
*! Everything takes explicit matrices. Nothing reads a fit, a struct, or Stata state.
*!
*! ORDERS (one place, so that no caller guesses):
*!   - theta_u  = (tau, a): thresholds, then the free echelon loadings a = A[i,j],
*!     i >= j, column-major (j outer, i inner). nth = number of thresholds.
*!   - theta*   = (tau, vec Lambda column-major, phi), phi = strict lower triangle of
*!     Phi, column-major: (2,1), (3,1), ..., (m,1), (3,2), ...
*!   - r, and the rows of R: the off-diagonal (j,k) of an m x m matrix, column by
*!     column: (2,1), (3,1), ..., (1,2), (3,2), ... that is rows j outer-to-inner
*!     within each column k, skipping j = k (dwlsmv_rotse_off).
*!   - Lambda block of R, J rows, G_u rows: vec Lambda, then phi; thresholds are NOT
*!     in R, G_u, G_T or J0 (R's threshold columns are 0, J's threshold rows are
*!     [I, 0], RS5(a)); dwlsmv_rotse_cov() adds them.
*!   The caller maps the fit's theta (spec 4.2: thresholds, lambda, beta, psi) to
*!   (tau, a) and the posted e(b) to theta*. That mapping is the posting layer's
*!   (unit A5, later), not this file's.
*!
*!   real colvector dwlsmv_rotse_r(L, Phi, eps)                  RS1  stationarity
*!   real matrix    dwlsmv_rotse_R(L, Phi, eps)                  RS4  dr/dtheta*'
*!   real matrix    dwlsmv_rotse_Gu(L, Minv)                     RS2  dg/dtheta_u'
*!   real matrix    dwlsmv_rotse_GT(L, Phi)                      RS2  dg/dY
*!   real matrix    dwlsmv_rotse_minv(A, L)                      (T')^-1 from A, Lambda
*!   real matrix    dwlsmv_rotse_J(L, Phi, eps, Minv)            RS3(ii) J0
*!   real matrix    dwlsmv_rotse_cov(J0, nth, ACOVu)             RS5(a), RS3 route 3
*!   real matrix    dwlsmv_rotse(A, L, Phi, eps, nth, ACOVu, J0)   the above, chained
*!   real scalar    dwlsmv_rotse_rank(Ds, w, R, ratio, cnd)      A3.2 rank and condition
*!
*! L, Phi are the production geomin output (dwlsmv_geomin(): Lambda = A (T')^-1,
*! Phi = T'T, after the sign and order rule). dwlsmv_geomin() does not return T:
*! (T')^-1 = (A'A)^-1 A'Lambda because A has full column rank (RS7), which
*! dwlsmv_rotse_minv() computes. m = 1: R has no rows, J0 = G_u, which is +-I.
*!
*! A singular R G_T (a degenerate critical point, RS6(i)) makes luinv() return a
*! matrix of missing values, so dwlsmv_rotse_J() returns missing. Callers test
*! missing(J0) and take the A3.4 route (rotated e(b), no e(V)). That path is not
*! exercised by test_rotse.do (no fixture has a degenerate critical point).
*!
*! Load order: dwlsmv_geomin (dwlsmv_geomin_dq), dwlsmv_numrank (rank function only).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The off-diagonal index pairs (row, col) of an m x m matrix, column by column.
real matrix dwlsmv_rotse_off(real scalar m)
{
    real matrix o
    real scalar j, k, q

    o = J(m * (m - 1), 2, .)
    q = 0
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= m; j++) {
            if (j != k) {
                q++
                o[q, 1] = j
                o[q, 2] = k
            }
        }
    }
    return(o)
}

// Position of phi_{ab} in phi, the strict lower triangle of an m x m matrix,
// column-major: m = 3 gives (2,1) 1, (3,1) 2, (3,2) 3. a and b in either order.
real scalar dwlsmv_rotse_phidx(real scalar a, real scalar b, real scalar m)
{
    real scalar lo, hi, c, pos

    lo  = min((a, b))
    hi  = max((a, b))
    pos = 0
    for (c = 1; c < lo; c++) {
        pos = pos + (m - c)
    }
    return(pos + (hi - lo))
}

// (T')^-1 = (A'A)^-1 A' Lambda (RS7). Missing if A'A is singular.
real matrix dwlsmv_rotse_minv(real matrix A, real matrix L)
{
    return(luinv(A' * A) * (A' * L))
}

// Result RS1: r = the off-diagonal of E = M - diag(M) Phi, M = Lambda' G_Q.
real colvector dwlsmv_rotse_r(real matrix L, real matrix Phi, real scalar eps)
{
    real matrix M, E, o
    real colvector out
    real scalar q

    M   = L' * dwlsmv_geomin_dq(L, eps)
    E   = M - diag(diagonal(M)) * Phi
    o   = dwlsmv_rotse_off(cols(L))
    out = J(rows(o), 1, .)
    for (q = 1; q <= rows(o); q++) {
        out[q] = E[o[q, 1], o[q, 2]]
    }
    return(out)
}

// Result RS4: the closed-form Jacobian R = dr/d(vec Lambda, phi)', m(m-1) x (p m + m(m-1)/2).
// With l_ij = lambda_ij^2 + eps, h = lambda / l, h' = (eps - lambda^2) / l^2,
// (G_Q)_ij = (2/m) h_ij Q_i, H_i[k,l] = (2/m) Q_i (delta_kl h'_ik + (2/m) h_ik h_il):
//   dr_jk/dlambda_il = delta_lj (G_Q)_ik + lambda_ij H_i[k,l]
//                      - Phi_jk (delta_lj (G_Q)_ij + lambda_ij H_i[j,l])
//   dr_jk/dphi_{jk}  = -M_jj
real matrix dwlsmv_rotse_R(real matrix L, real matrix Phi, real scalar eps)
{
    real matrix l2, Gq, h, hp, M, o, R
    real colvector Qi, Hk, Hj, dMjk, dMjj
    real scalar p, m, q, j, k, l

    p  = rows(L)
    m  = cols(L)
    l2 = L:^2 :+ eps
    Qi = exp(rowsum(ln(l2)) :/ m)
    Gq = (2 / m) :* L :/ l2 :* Qi
    h  = L :/ l2
    hp = (eps :- L:^2) :/ l2:^2
    M  = L' * Gq
    o  = dwlsmv_rotse_off(m)
    R  = J(rows(o), p * m + m * (m - 1) / 2, 0)
    for (q = 1; q <= rows(o); q++) {
        j = o[q, 1]
        k = o[q, 2]
        for (l = 1; l <= m; l++) {
            Hk   = (2 / m) :* Qi :* ((k == l) :* hp[., k] :+ (2 / m) :* h[., k] :* h[., l])
            Hj   = (2 / m) :* Qi :* ((j == l) :* hp[., j] :+ (2 / m) :* h[., j] :* h[., l])
            dMjk = (l == j) :* Gq[., k] :+ L[., j] :* Hk
            dMjj = (l == j) :* Gq[., j] :+ L[., j] :* Hj
            R[q, ((l - 1) * p + 1)..(l * p)] = (dMjk :- Phi[j, k] :* dMjj)'
        }
        R[q, p * m + dwlsmv_rotse_phidx(j, k, m)] = -M[j, j]
    }
    return(R)
}

// G_u = dg/d(a)' at fixed T, rows (vec Lambda, phi): dLambda = dA (T')^-1, so a
// free entry A[i,j] moves row i of Lambda by row j of Minv = (T')^-1. Phi does
// not move. One column per free echelon loading, column-major (j outer, i >= j).
real matrix dwlsmv_rotse_Gu(real matrix L, real matrix Minv)
{
    real matrix G, dL
    real scalar p, m, c, i, j

    p = rows(L)
    m = cols(L)
    G = J(p * m + m * (m - 1) / 2, p * m - m * (m - 1) / 2, 0)
    c = 0
    for (j = 1; j <= m; j++) {
        for (i = j; i <= p; i++) {
            c++
            dL = J(p, m, 0)
            dL[i, .] = Minv[j, .]
            G[|1, c \ p * m, c|] = vec(dL)
        }
    }
    return(G)
}

// G_T = dg/dY, Y the off-diagonal coordinates of T'dT (diag Y = 0), one column
// per off-diagonal (a,b) in dwlsmv_rotse_off() order: dLambda = -Lambda Y' Phi^-1,
// dPhi = Y + Y' (the single entry phi_{ab} moves by 1). Rows (vec Lambda, phi).
real matrix dwlsmv_rotse_GT(real matrix L, real matrix Phi)
{
    real matrix G, Y, Pi, o
    real scalar p, m, q

    p  = rows(L)
    m  = cols(L)
    o  = dwlsmv_rotse_off(m)
    Pi = luinv(Phi)
    G  = J(p * m + m * (m - 1) / 2, rows(o), 0)
    for (q = 1; q <= rows(o); q++) {
        Y = J(m, m, 0)
        Y[o[q, 1], o[q, 2]] = 1
        G[|1, q \ p * m, q|] = vec(-L * Y' * Pi)
        G[p * m + dwlsmv_rotse_phidx(o[q, 1], o[q, 2], m), q] = 1
    }
    return(G)
}

// Result RS3(ii): J0 = (I - G_T (R G_T)^-1 R) G_u, the Jacobian of the rotation
// map with respect to the free echelon loadings, rows (vec Lambda, phi). Minv =
// (T')^-1 at the same point. Missing if R G_T is singular (RS6(i)). m = 1: G_u.
real matrix dwlsmv_rotse_J(real matrix L, real matrix Phi, real scalar eps,
                           real matrix Minv)
{
    real matrix Gu, GT, R

    Gu = dwlsmv_rotse_Gu(L, Minv)
    if (cols(L) == 1) {
        return(Gu)
    }
    GT = dwlsmv_rotse_GT(L, Phi)
    R  = dwlsmv_rotse_R(L, Phi, eps)
    return(Gu - GT * luinv(R * GT) * (R * Gu))
}

// Result RS5(a): J = blockdiag(I_nth, J0), threshold rows [I, 0].
real matrix dwlsmv_rotse_jfull(real matrix J0, real scalar nth)
{
    real matrix Jf

    Jf = J(nth + rows(J0), nth + cols(J0), 0)
    if (nth > 0) {
        Jf[|1, 1 \ nth, nth|] = I(nth)
    }
    Jf[|nth + 1, nth + 1 \ ., .|] = J0
    return(Jf)
}

// Result RS3 route 3 / RS5(a): ACOV* = J ACOV_u J' over (tau, vec Lambda, phi),
// n* x n*, singular by construction (rank n_u). ACOV_u is in (tau, a) order and
// on whatever scale the caller holds (ACOV or e(V)): the product is linear in it.
// Symmetrized, (V + V')/2, which changes nothing beyond rounding for symmetric input.
real matrix dwlsmv_rotse_cov(real matrix J0, real scalar nth, real matrix ACOVu)
{
    real matrix Jf, V

    if (nth < 0 | nth != trunc(nth) | rows(ACOVu) != nth + cols(J0) |
        cols(ACOVu) != rows(ACOVu)) {
        errprintf("dwlsmv_rotse_cov(): ACOV_u is %g x %g, J0 has %g columns, nth %g\n",
                  rows(ACOVu), cols(ACOVu), cols(J0), nth)
        _error(3498, "dwlsmv_rotse_cov(): ACOV_u does not conform to J")
    }
    if (missing(J0) | missing(ACOVu)) {
        _error(3498, "dwlsmv_rotse_cov(): missing value in J0 or ACOV_u")
    }
    Jf = dwlsmv_rotse_jfull(J0, nth)
    V  = Jf * ACOVu * Jf'
    return((V + V') / 2)
}

// The chain a caller runs: (A, Lambda, Phi) from dwlsmv_geomin() and ACOV_u from
// the echelon fit -> ACOV* (and J0, for diagnostics). Returns missing if R G_T
// is singular; the caller then posts no e(V) (A3.4).
real matrix dwlsmv_rotse(real matrix A, real matrix L, real matrix Phi,
                              real scalar eps, real scalar nth, real matrix ACOVu,
                              real matrix J0)
{
    real matrix Minv

    if (rows(A) != rows(L) | cols(A) != cols(L) | rows(Phi) != cols(L) |
        cols(Phi) != cols(L) | rows(A) < cols(A)) {
        errprintf("dwlsmv_rotse(): A is %g x %g, Lambda %g x %g, Phi %g x %g\n",
                  rows(A), cols(A), rows(L), cols(L), rows(Phi), cols(Phi))
        _error(3498, "dwlsmv_rotse(): A, Lambda, Phi do not conform")
    }
    if (missing(A) | missing(L) | missing(Phi) | missing(eps) | eps <= 0) {
        _error(3498, "dwlsmv_rotse(): missing value or eps <= 0")
    }
    Minv = dwlsmv_rotse_minv(A, L)
    J0   = dwlsmv_rotse_J(L, Phi, eps, Minv)
    if (missing(J0)) {
        return(J(nth + rows(J0), nth + rows(J0), .))
    }
    return(dwlsmv_rotse_cov(J0, nth, ACOVu))
}

// Spec A3.2 (RS6(i)): are the SEs available, and how conditioned is the system.
//   Ds     D x n* Delta* at the rotated point, thresholds included
//   w      D x 1 the diagonal of W (positive)
//   R      k x n* the constraint Jacobian WITH its n_th zero threshold columns,
//          (J(k, nth, 0), dwlsmv_rotse_R(...)); k = 0 allowed
//   ratio  filled: sigma_min / sigma_max of W^1/2 Delta* Z, Z a basis of null(R)
//   cnd    filled: the condition number of the bordered matrix
//          [Delta*'W Delta*, R'; R, 0] (stored as an e() scalar, no cutoff)
// Returns 1 if the covariance exists: R has full row rank AND W^1/2 Delta* Z has
// full column rank, both by dwlsmv_numrank(), the ONE 1e-8 relative criterion
// (spec 4.8). The singular values of W^1/2 Delta* Z are those of W^1/2 Delta* P,
// P = I - R^+ R, the n* - k leading ones, so no null basis is formed. Numerically
// the rank is of the unsquared W^1/2 Delta* Z, not of its cross-product (RS6(i)).
real scalar dwlsmv_rotse_rank(real matrix Ds, real colvector w, real matrix R,
                              real scalar ratio, real scalar cnd)
{
    real matrix X, P, V, B, U, Vt
    real colvector sv, sr, sb
    real scalar n, k, nu, rx, rr

    n = cols(Ds)
    k = rows(R)
    if (rows(w) != rows(Ds) | (k > 0 & cols(R) != n) | n < 1) {
        errprintf("dwlsmv_rotse_rank(): Delta* %g x %g, w %g x 1, R %g x %g\n",
                  rows(Ds), cols(Ds), rows(w), rows(R), cols(R))
        _error(3498, "dwlsmv_rotse_rank(): arguments do not conform")
    }
    if (missing(Ds) | missing(w) | missing(R) | any(w :<= 0)) {
        _error(3498, "dwlsmv_rotse_rank(): missing value or nonpositive w")
    }
    X     = sqrt(w) :* Ds
    nu    = n - k
    ratio = 0
    cnd   = .
    if (nu < 1) {
        return(0)
    }
    if (k > 0) {
        P  = I(n) - pinv(R) * R
        rr = dwlsmv_numrank(R', sr, V)
        B  = (X' * X, R' \ R, J(k, k, 0))
    }
    else {
        P  = I(n)
        rr = 0
        B  = X' * X
    }
    // numrank() is called on a matrix with at least as many rows as columns
    rx = dwlsmv_numrank((rows(X) < n ? (X * P)' : X * P), sv, V)
    if (nu <= length(sv)) {
        ratio = (sv[1] > 0 ? sv[nu] / sv[1] : 0)
    }
    U  = J(0, 0, 0)
    Vt = J(0, 0, 0)
    sb = J(0, 1, 0)
    svd(B, U, sb, Vt)
    if (sb[length(sb)] > 0) {
        cnd = sb[1] / sb[length(sb)]
    }
    return(rr == k & rx == nu)
}

end
