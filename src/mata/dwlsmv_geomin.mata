*! dwlsmv -- oblique geomin rotation of an unrotated EFA solution (Module 12, 0.3.0)
*!
*! Executed under StataMP 18 on 2026-10-29. tests/unit/test_geomin.do: rotated
*! Lambda, Phi and Q match lavaan 0.7.2 geomin on Cases E, F, G at the theta-hat
*! tier (worst ratio 7.1e-3). Nothing in this file broke on the first run.
*!
*!     real matrix dwlsmv_geomin(real matrix A, real scalar eps,
*!                               real scalar nstart, real matrix Phi,
*!                               real rowvector info [, real scalar maxit])
*!
*! Returns the rotated loadings Lambda = A (T')^-1 and sets Phi = T'T, for
*! the T (unit-length columns) that minimizes the geomin criterion
*!
*!     Q(Lambda) = sum_i ( prod_j (lambda_ij^2 + eps) )^(1/m)
*!
*! A is p x m, the unrotated loadings with Phi = I (an echelon EFA fit); the
*! rotation leaves Lambda Phi Lambda' = A A', so sigma, F and T3 are those of
*! the unrotated fit (derivation @sec-geomin, Result GM1). Point estimates
*! only: nothing here touches theta or its ACOV.
*!
*! Conventions, all the user's (spec section 5, Module 12, 2026-10-29):
*!   eps      the caller's; Mplus's rule (.0001 / .001 / .01 at m = 2 / 3 / >= 4)
*!            is applied by the caller, so an override is just another eps
*!   rows     no normalization
*!   starts   the pinned GM4 grid (@sec-geomin): start 0 = I, start k a
*!            product of Givens rotations at angles pi*{k(q+1)c},
*!            c = (sqrt(5)-1)/2. No random-number state is read or changed.
*!   GPA      oblique gradient projection with Armijo step halving, run to
*!            the stall (no halving passes); a start is converged iff its
*!            final ||Gp||_F < 1e-5. The best converged start is selected;
*!            none converged is an error.
*!   output   columns reflected to positive sums (a sum of exactly 0 is not
*!            reflected), then put in the complete column order of spec
*!            Module 12 (dwlsmv_geomin_colorder(), owner 2026-11-18, A3.7
*!            revisited): key 1 ascending mean row index of the loadings with
*!            |lambda| >= 0.8 max (lavaan's "index"); key 2 descending sum of
*!            squared loadings; key 3 the loading columns compared item by
*!            item, larger first; key 4, only for identical loading columns,
*!            the arrangement whose Phi (strict lower triangle, row-major) is
*!            lexicographically largest. The order depends on the solution
*!            alone, never on which start reached it. Phi follows.
*!
*!     real rowvector dwlsmv_geomin_colorder(real matrix L, real matrix Phi)
*!         the permutation: L[., ord], Phi[ord, ord] is the posted order.
*!         Idempotent: on an already ordered solution it returns 1..m.
*!     real rowvector dwlsmv_geomin_tiekeys(real matrix L, real matrix Phi)
*!         for an ordered solution, 1 x (m-1): entry j is the key (1-4) that
*!         decides columns j and j+1, 0 when nothing does (identical L and
*!         Phi either way). Used by the output's tie note.
*!
*! info (1 x 8), the e() diagnostics of spec section 5 Module 12:
*!     1 starts attempted   2 converged starts   3 distinct minima among the
*!     converged   4 converged starts at the selected minimum   5 best Q
*!     6 its final ||Gp||_F   7 eps   8 Q gap from the best to the second
*!     distinct minimum (missing when there is only one)
*! Two converged starts are the same minimum iff their Q differ by at most
*! 1e-8 + 1e-6 Q_best (the distinct-minimum tolerance, spec Module 12).
*!
*! maxit (default 10000) caps GPA iterations per start; tests lower it to
*! reach the none-converged error.
*!
*! m = 1 has nothing to rotate: Lambda = A with the sign rule, Phi = 1.
*!
*! src/mata/dwlsmv_types.mata is not needed; this file stands alone.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// ---- the complete column order (spec Module 12; owner 2026-11-18) -------------------------
// Key 1, lavaan's "index": a mean of integer row indices, so its ties are exact.
real rowvector dwlsmv_geomin_index(real matrix L)
{
    real scalar j
    real rowvector v
    real colvector hi

    v = J(1, cols(L), .)
    for (j = 1; j <= cols(L); j++) {
        hi   = selectindex(abs(L[., j]) :>= 0.8 * max(abs(L[., j])))
        v[j] = mean(hi)
    }
    return(v)
}

// Which of keys 1-3 places column a before column b: -1 a first, +1 b first, 0 identical
// loading columns. key returns the deciding key (1, 2, 3) or 0. Exact comparisons throughout.
real scalar dwlsmv_geomin_cmp(real matrix L, real scalar a, real scalar b,
                              real rowvector idx, real rowvector ss, real scalar key)
{
    real scalar i

    key = 1
    if (idx[a] < idx[b]) {
        return(-1)
    }
    if (idx[a] > idx[b]) {
        return(1)
    }
    key = 2
    if (ss[a] > ss[b]) {
        return(-1)
    }
    if (ss[a] < ss[b]) {
        return(1)
    }
    key = 3
    for (i = 1; i <= rows(L); i++) {
        if (L[i, a] > L[i, b]) {
            return(-1)
        }
        if (L[i, a] < L[i, b]) {
            return(1)
        }
    }
    key = 0
    return(0)
}

// Phi's strict lower triangle, row-major
real rowvector dwlsmv_geomin_philow(real matrix P)
{
    real scalar i, j, k
    real rowvector v

    v = J(1, rows(P) * (rows(P) - 1) / 2, .)
    k = 0
    for (i = 2; i <= rows(P); i++) {
        for (j = 1; j < i; j++) {
            k++
            v[k] = P[i, j]
        }
    }
    return(v)
}

// -1 if u is lexicographically larger than v, +1 if smaller, 0 if equal
real scalar dwlsmv_geomin_lexcmp(real rowvector u, real rowvector v)
{
    real scalar k

    for (k = 1; k <= cols(u); k++) {
        if (u[k] > v[k]) {
            return(-1)
        }
        if (u[k] < v[k]) {
            return(1)
        }
    }
    return(0)
}

real rowvector dwlsmv_geomin_colorder(real matrix L, real matrix Phi)
{
    real scalar m, i, j, key, ng, g, c, more
    real rowvector idx, ss, ord, grp, cand, bestv, bestord, v
    real colvector pc
    real matrix gs, cnt, nper
    transmorphic info
    pointer(real matrix) rowvector perms

    m = cols(L)
    if (m < 2) {
        return(1..m)
    }
    idx = dwlsmv_geomin_index(L)
    ss  = colsum(L :^ 2)
    // keys 1-3: insertion sort with the exact comparator (m is small)
    ord = 1..m
    for (i = 2; i <= m; i++) {
        j = i
        while (j > 1) {
            if (dwlsmv_geomin_cmp(L, ord[j - 1], ord[j], idx, ss, key) <= 0) {
                break
            }
            c = ord[j - 1]
            ord[j - 1] = ord[j]
            ord[j] = c
            j--
        }
    }
    // key 4: runs of identical loading columns (comparator 0 between neighbours)
    grp = J(1, m, 0)
    ng  = 0
    for (i = 1; i < m; i++) {
        if (dwlsmv_geomin_cmp(L, ord[i], ord[i + 1], idx, ss, key) == 0) {
            if (grp[i] == 0) {
                ng++
                grp[i] = ng
            }
            grp[i + 1] = grp[i]
        }
    }
    if (ng == 0) {
        return(ord)
    }
    // enumerate every arrangement within the runs; keep the lexicographically largest Phi
    gs    = J(ng, 2, .)
    perms = J(1, ng, NULL)
    for (g = 1; g <= ng; g++) {
        v = selectindex(grp :== g)
        gs[g, .] = (v[1], cols(v))
        perms[g] = &(J(0, cols(v), .))
        info = cvpermutesetup((1::cols(v)))
        while ((pc = cvpermute(info)) != J(0, 1, .)) {
            *perms[g] = *perms[g] \ pc'
        }
    }
    nper  = J(1, ng, .)
    for (g = 1; g <= ng; g++) {
        nper[g] = rows(*perms[g])
    }
    cnt   = J(1, ng, 1)
    bestv = J(1, 0, .)
    more  = 1
    while (more) {
        cand = ord
        for (g = 1; g <= ng; g++) {
            for (c = 1; c <= gs[g, 2]; c++) {
                cand[gs[g, 1] + c - 1] = ord[gs[g, 1] + (*perms[g])[cnt[g], c] - 1]
            }
        }
        v = dwlsmv_geomin_philow(Phi[cand, cand])
        if (cols(bestv) == 0) {
            bestv = v
            bestord = cand
        }
        else if (dwlsmv_geomin_lexcmp(v, bestv) < 0) {
            bestv   = v
            bestord = cand
        }
        // odometer over the runs
        g = 1
        more = 0
        while (g <= ng) {
            if (cnt[g] < nper[g]) {
                cnt[g] = cnt[g] + 1
                more = 1
                break
            }
            cnt[g] = 1
            g++
        }
    }
    return(bestord)
}

real rowvector dwlsmv_geomin_tiekeys(real matrix L, real matrix Phi)
{
    real scalar m, j, key, r
    real rowvector idx, ss, k, u, v, sw

    m = cols(L)
    if (m < 2) {
        return(J(1, 0, .))
    }
    idx = dwlsmv_geomin_index(L)
    ss  = colsum(L :^ 2)
    k   = J(1, m - 1, .)
    for (j = 1; j < m; j++) {
        r = dwlsmv_geomin_cmp(L, j, j + 1, idx, ss, key)
        if (key == 0) {
            sw = 1..m
            sw[j] = j + 1
            sw[j + 1] = j
            u = dwlsmv_geomin_philow(Phi)
            v = dwlsmv_geomin_philow(Phi[sw, sw])
            key = (dwlsmv_geomin_lexcmp(u, v) != 0 ? 4 : 0)
        }
        k[j] = key
    }
    return(k)
}

real scalar dwlsmv_geomin_q(real matrix L, real scalar eps)
{
    return(sum(exp(rowsum(ln(L:^2 :+ eps)) :/ cols(L))))
}

// dQ/dLambda, elementwise (@sec-geomin GM2)
real matrix dwlsmv_geomin_dq(real matrix L, real scalar eps)
{
    real matrix l2
    real colvector r

    l2 = L:^2 :+ eps
    r  = exp(rowsum(ln(l2)) :/ cols(L))
    return((2 / cols(L)) :* L :/ l2 :* r)
}

// Lambda = A (T')^-1
real matrix dwlsmv_geomin_rot(real matrix A, real matrix T)
{
    return(A * luinv(T'))
}

// d Q(A (T')^-1) / dT = -(Lambda' G_Q T^-1)'
real matrix dwlsmv_geomin_grad(real matrix A, real matrix T, real scalar eps)
{
    real matrix L

    L = dwlsmv_geomin_rot(A, T)
    return(-(L' * dwlsmv_geomin_dq(L, eps) * luinv(T))')
}

real matrix dwlsmv_geomin_unit(real matrix X)
{
    return(X :/ sqrt(colsum(X:^2)))
}

// One start: returns T at the stall; f and s (= ||Gp||_F there) by reference.
real matrix dwlsmv_geomin_gpa(real matrix A, real matrix T0, real scalar eps,
                              real scalar maxit, real scalar f, real scalar s)
{
    real matrix T, Tn, G, Gp
    real scalar al, fn, it, h, ok

    T  = T0
    al = 1
    f  = dwlsmv_geomin_q(dwlsmv_geomin_rot(A, T), eps)
    G  = dwlsmv_geomin_grad(A, T, eps)
    for (it = 0; it <= maxit; it++) {
        Gp = G - T * diag(colsum(T :* G))
        s  = sqrt(sum(Gp:^2))
        al = 2 * al
        ok = 0
        for (h = 1; h <= 60; h++) {
            Tn = dwlsmv_geomin_unit(T - al * Gp)
            fn = dwlsmv_geomin_q(dwlsmv_geomin_rot(A, Tn), eps)
            if (fn < f - 0.5 * s^2 * al) {
                ok = 1
                break
            }
            al = al / 2
        }
        if (!ok) {
            break
        }
        T = Tn
        f = fn
        G = dwlsmv_geomin_grad(A, T, eps)
    }
    // s at the returned T: after a stall T did not move, but a loop that ends
    // at maxit has stepped once past the last s (review finding, 2026-10-29)
    s = sqrt(sum((G - T * diag(colsum(T :* G))):^2))
    return(T)
}

// Start k of the GM4 grid (k = 0 is the identity).
real matrix dwlsmv_geomin_start(real scalar m, real scalar k)
{
    real matrix S, R
    real scalar c, q, i, j, a

    S = I(m)
    if (k == 0) {
        return(S)
    }
    c = (sqrt(5) - 1) / 2
    q = 0
    for (i = 1; i < m; i++) {
        for (j = i + 1; j <= m; j++) {
            q++
            a = pi() * mod(k * (q + 1) * c, 1)
            R = I(m)
            R[i, i] = cos(a)
            R[j, j] = cos(a)
            R[i, j] = -sin(a)
            R[j, i] = sin(a)
            S = S * R
        }
    }
    return(S)
}

real matrix dwlsmv_geomin(real matrix A, real scalar eps, real scalar nstart,
                          real matrix Phi, real rowvector info,
                          | real scalar maxit)
{
    real matrix T, Tb, L
    real colvector f, s, conv, fc, sg
    real rowvector ord
    real scalar m, k, b, nconv, tol, ndist, nbest, gap, fk, sk, fprev, gtol

    if (args() < 6) {
        maxit = 10000
    }
    m = cols(A)
    if (rows(A) < 1 | m < 1 | missing(A)) {
        errprintf("dwlsmv_geomin(): A must be a nonempty p x m matrix with no missing value\n")
        _error(3498, "dwlsmv_geomin(): invalid A")
    }
    if (m > rows(A)) {
        errprintf("dwlsmv_geomin(): %g factors exceed %g indicators\n", m, rows(A))
        _error(3498, "dwlsmv_geomin(): more factors than indicators")
    }
    if (missing(eps) | eps <= 0) {
        errprintf("dwlsmv_geomin(): eps must be positive, got %g\n", eps)
        _error(3498, "dwlsmv_geomin(): eps must be positive")
    }
    if (missing(nstart) | nstart < 1 | nstart != trunc(nstart)) {
        errprintf("dwlsmv_geomin(): nstart must be a positive integer, got %g\n", nstart)
        _error(3498, "dwlsmv_geomin(): nstart must be a positive integer")
    }

    if (m == 1) {
        L   = (sum(A) < 0 ? -A : A)
        Phi = 1
        info = (1, 1, 1, 1, dwlsmv_geomin_q(L, eps), 0, eps, .)
        return(L)
    }

    // the convergence rule, stated once: it both counts converged starts and
    // decides which may be selected (mutation round, 2026-10-29: two copies
    // of 1e-5 could drift apart unseen)
    gtol = 1e-5
    f    = J(nstart, 1, .)
    s    = J(nstart, 1, .)
    b    = 0
    for (k = 1; k <= nstart; k++) {
        T = dwlsmv_geomin_gpa(A, dwlsmv_geomin_start(m, k - 1), eps, maxit, fk, sk)
        f[k] = fk
        s[k] = sk
        if (sk < gtol) {
            if (b == 0) {
                b  = k
                Tb = T
            }
            else if (fk < f[b]) {
                b  = k
                Tb = T
            }
        }
    }
    conv  = (s :< gtol)
    nconv = sum(conv)
    if (nconv == 0) {
        errprintf("dwlsmv_geomin(): none of the %g rotation starts converged (||Gp|| < 1e-5)\n", nstart)
        _error(3498, "dwlsmv_geomin(): no rotation start converged")
    }

    // distinct minima among the converged starts, ascending in Q
    tol   = 1e-8 + 1e-6 * abs(f[b])
    fc    = sort(select(f, conv), 1)
    ndist = 1
    nbest = 0
    gap   = .
    fprev = fc[1]
    for (k = 1; k <= rows(fc); k++) {
        if (fc[k] - f[b] <= tol) {
            nbest++
        }
        if (fc[k] - fprev > tol) {
            ndist++
            if (ndist == 2) {
                gap = fc[k] - f[b]
            }
            fprev = fc[k]
        }
    }

    // sign: positive column sums (a sum of exactly 0 is not reflected); order: the complete
    // rule of spec Module 12 (lavaan's "index", then dwlsmv_geomin_colorder()'s keys 2-4)
    L  = dwlsmv_geomin_rot(A, Tb)
    sg = (colsum(L) :< 0)'
    sg = 1 :- 2 :* sg
    L  = L :* sg'
    Tb = Tb :* sg'
    ord = dwlsmv_geomin_colorder(L, Tb' * Tb)
    L   = L[., ord]
    Tb  = Tb[., ord]

    Phi  = Tb' * Tb
    info = (nstart, nconv, ndist, nbest, f[b], s[b], eps, gap)
    return(L)
}

end
