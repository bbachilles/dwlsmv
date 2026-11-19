*! Rotated-EFA standard errors, numerical core (release A, unit A5): R, J and the
*! rotated covariance against the tier-1 fixture.
*!
*! Drafted 2026-11-17 by a read-only agent. Executed under StataMP 18 by the main session the
*! same day. Fixed on first run: tr_ar1() (Mata has no column-by-row colon outer product) and
*! check 12's expectation at s = -1 (D ACOV_u D, not ACOV_u).
*!
*! Reference: fixtures/expected/rotse_tier1_jacobians.csv (99q_tier1_rot_jacobians.R,
*! base R): R by complex step (h 1e-30), J by complex-step re-rotation, at the Newton
*! root and the GPA stall point, Cases E (m 2), F (m 3), G (binary, m 2). Mata's
*! dwlsmv_rotse_*() functions are the closed forms of @sec-rot-se RS4 and RS3(ii); the
*! fixture values do not use them (different representation and algorithm; both share r,
*! i.e. GM2's gradient). No fit runs here: A, Lambda, Phi, T come from the fixture, so a
*! failure is the rotation-SE code and not an estimator. nth is the fixture's, and the
*! threshold rows of J are [I, 0] (RS5(a)).
*!
*! Tolerances, fixed in this file BEFORE the first run (2026-11-17):
*!   R vs fixture, both points          1e-12 relative to max|R| (spec @sec-rot-se
*!       verification table, layer R; the fixture's own closed-vs-complex-step gap is
*!       1.4e-16 to 4.0e-16, so the Mata port has four decades of room)
*!   J at the root vs fixture           1e-10 relative to max|J| (table, layer J; the
*!       fixture's own gap is 2.0e-16 to 9.3e-16)
*!   (T')^-1 recovered, r, R J = 0, RS5(b) invariants, covariance assembly
*!                                      1e-12 (arithmetic only; no stochastic content)
*!   J at the FIXTURE stall point vs root   1.0e-7  } CHOSEN AFTER MEASURING, for the
*!   J at MATA's stall point vs root        2.5e-7  } owner. Measured (99q diag
*!       J_stall_vs_root, relative to max|J|): E 2.3e-9, F 5.3e-8, G 2.0e-9; the covariance
*!       built from the same J differs from the root covariance by 3.0e-9, 4.3e-8, 1.5e-9
*!       (scratch measurement, AR(1) ACOV_u below). 1.0e-7 is 1.9x the worst J gap and
*!       2.3x the worst covariance gap; 2.5e-7 is 4.7x, because Mata's own stall point need
*!       not equal R's (both stop at ||Gp|| < 1e-5; test_geomin measured <= 1.2e-8). These are
*!       not loosened to pass: a Mata value above 2.5e-7 is a finding. They say the stall
*!       point is within the A3.15 diagnostic (SEs differ by at most 1.1e-7 relative).
*!   rank/condition function            known answers, 1e-12 relative (ratio 1e-6 where
*!       the answer is a power of ten read through a singular value)
*!
*! Checks (per case E, F, G unless stated):
*!   1  (T')^-1 = (A'A)^-1 A'Lambda against solve(T') from the fixture T, both points
*!   2  r against the fixture's residual vector, both points (port check, abs 1e-12)
*!   3  R against fixture R, both points                                    DECIDING for R
*!   4  J at the root point against the fixture's re-rotation J             DECIDING for J
*!   5  R J = 0. Holds for ANY R by construction (RS3(i)); it does not certify R
*!   6  J at the fixture's stall (Lambda, Phi) vs the fixture root J
*!   7  Mata's own dwlsmv_geomin() on the fixture A: same Lambda/Phi as the fixture's
*!      root within 1e-5 (labelling guard), then J there vs the root J
*!   8  RS5(b): grad h J = grad h_u for the communalities h_i = (Lambda Phi Lambda')_ii
*!      (valid at any point, for any R; certifies G_u and G_T), and grad h G_T = 0
*!   9  covariance assembly: synthetic ACOV_u = AR(1), 0.5^|i-j|, J0 = fixture root J;
*!      nine inline summaries from base R (ref.R, same fixture); threshold block equals
*!      ACOV_u,tt (RS5(a)); the tau x loading block equals ACOV_u[tau,a] J0'
*!   10 end to end, dwlsmv_rotse(): from the fixture stall point and from Mata's own
*!      geomin, against the covariance of check 9
*!   11 mutation controls (the comparator is not vacuous): a flipped phi term of R moves
*!      R and J by > 1e-6 against the fixture while R J stays 0 (so check 5 is blind to it)
*!   12 m = 1: J0 = +-I, covariance = ACOV_u
*!   13 dwlsmv_rotse_rank(): known answers, the 1e-8 edge from both sides via the shared
*!      dwlsmv_numrank(), unsquared weights, rank-deficient R, k = 0, bordered condition
*!   14 error paths, rc 3498 with the message
*!
*! WHAT THIS CERTIFIES: the closed forms R, G_u, G_T, J and the product J ACOV J' as coded
*! in Mata reproduce the independent tier-1 values on E, F, G to the tolerances above, at
*! the root and, within the stated gap, at the stall point. WHAT IT DOES NOT: the bordered
*! covariance or Delta* (route 1, needs a fit), agreement with lavaan's SEs (tier 2), the
*! theta to (tau, a) mapping from a fit, calibration (A6), m >= 4, several minima, a
*! degenerate critical point (the missing-J path), that Delta* J = Delta_u (needs Delta*).
*!
*! MUTATION HOOK (to run by hand, not in the suite). In src/mata/dwlsmv_rotse.mata change
*!     R[q, p * m + dwlsmv_rotse_phidx(j, k, m)] = -M[j, j]
*! to
*!     R[q, p * m + dwlsmv_rotse_phidx(j, k, m)] = M[j, j]
*! Expected: exit 9 at check 3 (R off by 2.2e-3, 9.0e-3, 4.3e-3 relative on E, F, G),
*! checks 4, 6, 7, 10 (J off by 3.3e-3, 2.6e-2, 1.1e-2); checks 5 and 8 STILL PASS (R J = 0
*! and RS5(b) hold for any R). Two more, measured in base R at the root: dropping the
*! Phi_jk term of the loading block moves R by 8.1e-3, 3.6e-3, 6.0e-3; replacing
*! (k == l) by (j == l) in the h' term of Hk, 0.92 to 1.0 (all caught by check 3).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 120

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_numrank.mata           ///
             src/mata/dwlsmv_geomin.mata            ///
             src/mata/dwlsmv_rotse.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_rotse: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/expected/rotse_tier1_jacobians.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: FS = st_sdata(., ("case", "point", "quantity"))
mata: FN = strtoreal(st_sdata(., ("row", "col", "value")))

mata:

// one fixture matrix, from the long form (case, point, quantity, row, col, value)
real matrix tr_get(string matrix FS, real matrix FN, string scalar cs,
                   string scalar pt, string scalar qn)
{
    real colvector ix
    real matrix M
    real scalar k

    ix = selectindex((FS[., 1] :== cs) :& (FS[., 2] :== pt) :& (FS[., 3] :== qn))
    if (length(ix) == 0) {
        errprintf("tr_get(): no fixture rows for %s %s %s\n", cs, pt, qn)
        return(J(0, 0, .))
    }
    M = J(max(FN[ix, 1]), max(FN[ix, 2]), .)
    for (k = 1; k <= length(ix); k++) {
        M[FN[ix[k], 1], FN[ix[k], 2]] = FN[ix[k], 3]
    }
    return(M)
}

// max |X - Y| / max |Y|; missing if the shapes differ or either has a missing value
real scalar tr_rel(real matrix X, real matrix Y)
{
    if (rows(X) != rows(Y) | cols(X) != cols(Y) | missing(X) | missing(Y)) {
        return(.)
    }
    return(max(abs(X - Y)) / max(abs(Y)))
}

// print one check and return 1 if it fails (a missing value fails)
real scalar tr_chk(string scalar label, real scalar v, real scalar tol)
{
    printf("{txt}    %-58s %10.3e  (tol %8.1e)%s\n", label, v, tol,
           (v <= tol ? "" : "  FAIL"))
    return(!(v <= tol))
}

// the nine inline summaries of the covariance of check 9, from base R on the same fixture
real colvector tr_expect(string scalar cs)
{
    if (cs == "E") {
        return((24.668142406348228 \ 53.003209355336928 \ 6.7835252152769048 \ 1 \
                0.00011936799479574086 \ 1.2338188575454041 \ 0.22891935640767871 \
                2.0108825794948362 \ 0.58477983665624833))
    }
    if (cs == "F") {
        return((51.475970736614372 \ 113.76124835797513 \ 9.5941347123897192 \ 1 \
                3.656124072143142e-09 \ 1.0673927721664827 \ 0.23305381402922515 \
                2.1295445888001563 \ 0.058871183163861134))
    }
    return((20.302928009938945 \ 34.758538364300769 \ 6.9046520204540913 \ 1 \
            0.0077034290935196438 \ 1.246020047515485 \ 0.27898171738909439 \
            2.9194248050819582 \ 0.68446136163151594))
}

real colvector tr_summ(real matrix V, real scalar nth)
{
    real scalar n

    n = rows(V)
    return((sum(diagonal(V)) \ sum(V) \ sqrt(sum(V:^2)) \ V[1, 1] \ V[1, nth + 2] \
            V[nth + 1, nth + 1] \ V[nth + 1, nth + 3] \ V[n, n] \ V[n, nth + 1]))
}

real matrix tr_ar1(real scalar n)
{
    // Mata's colon operators do not form an outer difference from a column and
    // a row vector, so both are expanded to n x n first.
    return(0.5 :^ abs((1::n) * J(1, n, 1) :- J(n, 1, 1) * (1..n)))
}

// gradient of the communalities h_i = (Lambda Phi Lambda')_ii with respect to
// (vec Lambda, phi), and of h_i = sum_j A_ij^2 with respect to the free a
void tr_gradh(real matrix A, real matrix L, real matrix Phi, real matrix gh,
              real matrix ghu)
{
    real matrix LP
    real scalar p, m, i, j, a, b, c

    p   = rows(L)
    m   = cols(L)
    LP  = L * Phi
    gh  = J(p, p * m + m * (m - 1) / 2, 0)
    ghu = J(p, p * m - m * (m - 1) / 2, 0)
    for (i = 1; i <= p; i++) {
        for (j = 1; j <= m; j++) {
            gh[i, (j - 1) * p + i] = 2 * LP[i, j]
        }
        for (b = 1; b < m; b++) {
            for (a = b + 1; a <= m; a++) {
                gh[i, p * m + dwlsmv_rotse_phidx(a, b, m)] = 2 * L[i, a] * L[i, b]
            }
        }
    }
    c = 0
    for (j = 1; j <= m; j++) {
        for (i = j; i <= p; i++) {
            c++
            ghu[i, c] = 2 * A[i, j]
        }
    }
}

real scalar tr_case(string matrix FS, real matrix FN, string scalar cs)
{
    real matrix dm, A, Lf, Pf, Tf, Rf, Jroot, Lr, Pr, Tr, Rm, Jm, J0, Gu, GT, ACOV
    real matrix Ls, Ps, Ms, Js, Lo, Po, Jo, gh, ghu, V, Vroot, Jfull, Rr, Rn, o
    real rowvector info
    real colvector ex, got, rf
    real scalar nf, p, m, eps, nth, na, nu, nstar, k, r, q, phix, rd
    string rowvector pts

    nf    = 0
    dm    = tr_get(FS, FN, cs, "both", "dims")
    p     = dm[1]
    m     = dm[2]
    eps   = dm[3]
    nth   = dm[4]
    na    = dm[5]
    nu    = nth + na
    nstar = nth + p * m + m * (m - 1) / 2
    printf("{txt}  case %s (p = %g, m = %g, eps = %g, nth = %g, n_u = %g, n* = %g)\n",
           cs, p, m, eps, nth, nu, nstar)
    if (dm[6] != m * (m - 1) | dm[7] != p * m + m * (m - 1) / 2 | na != p * m - m * (m - 1) / 2) {
        printf("{err}    fixture dims are not the expected shapes\n")
        return(1)
    }
    A     = tr_get(FS, FN, cs, "root", "A")
    Jroot = tr_get(FS, FN, cs, "root", "J")
    Jfull = dwlsmv_rotse_jfull(Jroot, nth)
    ACOV  = tr_ar1(nu)

    pts = ("stall", "root")
    for (k = 1; k <= 2; k++) {
        Lf = tr_get(FS, FN, cs, pts[k], "Lambda")
        Pf = tr_get(FS, FN, cs, pts[k], "Phi")
        Tf = tr_get(FS, FN, cs, pts[k], "T")
        Rf = tr_get(FS, FN, cs, pts[k], "R")
        rf = tr_get(FS, FN, cs, pts[k], "r")
        printf("{txt}   point %s\n", pts[k])
        // 1  (T')^-1 from A and Lambda
        nf = nf + tr_chk("1 (T')^-1 recovered vs solve(T') of the fixture T",
                         tr_rel(dwlsmv_rotse_minv(A, Lf), luinv(Tf')), 1e-12)
        // 2  r (port check)
        r  = max(abs(dwlsmv_rotse_r(Lf, Pf, eps) - rf))
        nf = nf + tr_chk("2 r vs fixture residual vector (max abs)", r, 1e-12)
        // 3  R
        nf = nf + tr_chk("3 R vs fixture R, relative to max|R|",
                         tr_rel(dwlsmv_rotse_R(Lf, Pf, eps), Rf), 1e-12)
    }

    // 4  J at the root against the fixture's re-rotation
    Lr = tr_get(FS, FN, cs, "root", "Lambda")
    Pr = tr_get(FS, FN, cs, "root", "Phi")
    Tr = tr_get(FS, FN, cs, "root", "T")
    Rr = tr_get(FS, FN, cs, "root", "R")
    J0 = dwlsmv_rotse_J(Lr, Pr, eps, dwlsmv_rotse_minv(A, Lr))
    nf = nf + tr_chk("4 J at the root vs re-rotation J, relative to max|J|",
                     tr_rel(J0, Jroot), 1e-10)

    // 5  R J = 0 (holds for any R)
    Rn = dwlsmv_rotse_R(Lr, Pr, eps)
    r  = max(abs(Rn * J0)) / (max(abs(Rn)) * max(abs(J0)))
    nf = nf + tr_chk("5 R J = 0, relative to max|R| max|J| (any R passes)", r, 1e-12)

    // 6  J at the fixture's stall point
    Ls = tr_get(FS, FN, cs, "stall", "Lambda")
    Ps = tr_get(FS, FN, cs, "stall", "Phi")
    Js = dwlsmv_rotse_J(Ls, Ps, eps, dwlsmv_rotse_minv(A, Ls))
    nf = nf + tr_chk("6 J at the fixture stall point vs root J", tr_rel(Js, Jroot), 1.0e-7)

    // 7  Mata's own stall point
    Lo = dwlsmv_geomin(A, eps, 30, Po, info)
    r  = max((max(abs(Lo - Lr)), max(abs(Po - Pr))))
    nf = nf + tr_chk("7a Mata geomin Lambda, Phi vs fixture root (max abs)", r, 1e-5)
    Ms = dwlsmv_rotse_minv(A, Lo)
    Jo = dwlsmv_rotse_J(Lo, Po, eps, Ms)
    nf = nf + tr_chk("7b J at Mata's stall point vs root J", tr_rel(Jo, Jroot), 2.5e-7)

    // 8  RS5(b): grad h J = grad h_u at the stall point (any R), and grad h G_T = 0
    tr_gradh(A, Ls, Ps, gh, ghu)
    nf = nf + tr_chk("8a grad h J - grad h_u, relative to max|grad h_u|",
                     max(abs(gh * Js - ghu)) / max(abs(ghu)), 1e-12)
    nf = nf + tr_chk("8b grad h G_T (max abs)",
                     max(abs(gh * dwlsmv_rotse_GT(Ls, Ps))), 1e-12)

    // 9  covariance assembly from the fixture's root J
    V = dwlsmv_rotse_cov(Jroot, nth, ACOV)
    ex  = tr_expect(cs)
    got = tr_summ(V, nth)
    r   = 0
    for (q = 1; q <= 9; q++) {
        r = max((r, dwlsmv_gapratio(got[q], ex[q], 1e-13, 1e-12)))
    }
    nf = nf + tr_chk("9a nine summaries vs base R (gap / (1e-13 + 1e-12 |ref|))", r, 1)
    nf = nf + tr_chk("9b threshold block = ACOV_u,tt (max abs)",
                     max(abs(V[|1, 1 \ nth, nth|] - ACOV[|1, 1 \ nth, nth|])), 1e-14)
    nf = nf + tr_chk("9c tau x rotated block = ACOV_u[tau,a] J0' (relative)",
                     tr_rel(V[|1, nth + 1 \ nth, nstar|],
                            ACOV[|1, nth + 1 \ nth, nu|] * Jroot'), 1e-12)
    nf = nf + tr_chk("9d symmetric (max abs)", max(abs(V - V')), 1e-15)
    Vroot = V

    // 10  end to end
    V = dwlsmv_rotse(A, Ls, Ps, eps, nth, ACOV, J0)
    nf = nf + tr_chk("10a acov() from the fixture stall point vs root cov",
                     tr_rel(V, Vroot), 1.0e-7)
    V = dwlsmv_rotse(A, Lo, Po, eps, nth, ACOV, J0)
    nf = nf + tr_chk("10b acov() from Mata's own geomin vs root cov",
                     tr_rel(V, Vroot), 2.5e-7)

    // 11  mutation controls: flip the phi term of R in the fixture's own R
    Rm = Rr
    o  = dwlsmv_rotse_off(m)
    for (q = 1; q <= m * (m - 1); q++) {
        phix = p * m + dwlsmv_rotse_phidx(o[q, 1], o[q, 2], m)
        Rm[q, phix] = -Rm[q, phix]
    }
    Gu = dwlsmv_rotse_Gu(Lr, dwlsmv_rotse_minv(A, Lr))
    GT = dwlsmv_rotse_GT(Lr, Pr)
    Jm = Gu - GT * luinv(Rm * GT) * (Rm * Gu)
    rd = tr_rel(Rm, Rr)
    nf = nf + tr_chk("11a control: flipped-phi R differs from fixture R by more than 1e-6",
                     (rd < . & rd > 1e-6 ? 0 : 1), 0)
    rd = tr_rel(Jm, Jroot)
    nf = nf + tr_chk("11b control: its J differs from fixture J by more than 1e-6",
                     (rd < . & rd > 1e-6 ? 0 : 1), 0)
    nf = nf + tr_chk("11c control: R J = 0 still holds for it (check 5 is blind)",
                     max(abs(Rm * Jm)) / (max(abs(Rm)) * max(abs(Jm))), 1e-12)
    return(nf)
}

// 12  m = 1: nothing to rotate; J0 = +-I and the covariance is ACOV_u
real scalar tr_m1()
{
    real matrix A, L, Phi, J0, V, ACOV
    real rowvector info
    real scalar nf, s

    nf   = 0
    ACOV = tr_ar1(6)
    for (s = -1; s <= 1; s = s + 2) {
        A  = s :* (0.7 \ 0.6 \ 0.5 \ 0.1)
        L  = dwlsmv_geomin(A, 1e-2, 30, Phi, info)
        V  = dwlsmv_rotse(A, L, Phi, 1e-2, 2, ACOV, J0)
        nf = nf + tr_chk("12 m = 1: J0 = s I, s the sign geomin applied (max abs)",
                         max(abs(J0 - s :* I(4))), 1e-15)
        // A reflection (s = -1) flips the sign of the threshold x loading block,
        // so the expected covariance is D ACOV_u D, D = diag(I_2, s I_4).
        // Corrected by the main session 2026-11-17: the draft expected ACOV_u
        // itself, which fails at s = -1 by exactly 2 x 0.5 (AR(1) neighbours).
        nf = nf + tr_chk("12 m = 1: covariance = D ACOV_u D (max abs)",
                         max(abs(V - diag((1, 1, s, s, s, s)) * ACOV * diag((1, 1, s, s, s, s)))), 1e-15)
    }
    return(nf)
}

// 13  the rank and condition function on matrices with known answers
real scalar tr_rank()
{
    real matrix Ds, R
    real colvector w
    real scalar nf, ok, ratio, cnd, g

    nf = 0
    w  = J(4, 1, 1)
    R  = (0, 0, 0, 1)
    g  = (3 + sqrt(5)) / 2
    ok = dwlsmv_rotse_rank(I(4), w, R, ratio, cnd)
    nf = nf + tr_chk("13a identity, R = e4': returns 1", 1 - ok, 0)
    nf = nf + tr_chk("13a ratio = 1", abs(ratio - 1), 1e-12)
    nf = nf + tr_chk("13a bordered condition = (3 + sqrt 5)/2, relative", abs(cnd - g) / g, 1e-12)

    // the shared 1e-8 relative cut, from both sides
    ok = dwlsmv_rotse_rank(diag((1, 1, 2e-8, 1)), w, R, ratio, cnd)
    nf = nf + tr_chk("13b sigma ratio 2e-8: kept (returns 1)", 1 - ok, 0)
    ok = dwlsmv_rotse_rank(diag((1, 1, 5e-9, 1)), w, R, ratio, cnd)
    nf = nf + tr_chk("13b sigma ratio 5e-9: dropped (returns 0)", ok, 0)
    nf = nf + tr_chk("13b ratio reported, relative to 5e-9", abs(ratio - 5e-9) / 5e-9, 1e-6)

    // weights enter as sqrt(w), unsquared: w = 1e6 gives ratio 1e-3, not 1e-6
    ok = dwlsmv_rotse_rank(I(4), (1e6 \ 1 \ 1 \ 1), R, ratio, cnd)
    nf = nf + tr_chk("13c w = 1e6: ratio 1e-3 (not 1e-6), relative", abs(ratio - 1e-3) / 1e-3, 1e-6)
    ok = dwlsmv_rotse_rank(I(4), (1e18 \ 1 \ 1 \ 1), R, ratio, cnd)
    nf = nf + tr_chk("13c w = 1e18: ratio 1e-9, below the cut (returns 0)", ok, 0)

    // the singular direction lies in range(R'): the covariance still exists
    Ds = diag((1, 1, 1, 0))
    ok = dwlsmv_rotse_rank(Ds, w, R, ratio, cnd)
    nf = nf + tr_chk("13d Delta* null along R' only: returns 1", 1 - ok, 0)
    Ds = diag((1, 1, 0, 1))
    ok = dwlsmv_rotse_rank(Ds, w, R, ratio, cnd)
    nf = nf + tr_chk("13d Delta* null inside null(R): returns 0", ok, 0)

    // R without full row rank, and no constraints at all
    ok = dwlsmv_rotse_rank(I(4), w, (0, 0, 0, 1 \ 0, 0, 0, 2), ratio, cnd)
    nf = nf + tr_chk("13e R rank 1 of 2 rows: returns 0", ok, 0)
    ok = dwlsmv_rotse_rank(I(4), w, J(0, 4, 0), ratio, cnd)
    nf = nf + tr_chk("13f k = 0: returns 1", 1 - ok, 0)
    nf = nf + tr_chk("13f k = 0: condition number 1", abs(cnd - 1), 1e-12)
    return(nf)
}

end

local nfail = 0
foreach cs in E F G {
    mata: st_local("n1", strofreal(tr_case(FS, FN, "`cs'")))
    local nfail = `nfail' + `n1'
}
mata: st_local("n1", strofreal(tr_m1()))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tr_rank()))
local nfail = `nfail' + `n1'

* 14  error paths: rc 3498 AND the message (as test_geomin.do check 9)
mata: TA = (0.7, 0.1 \ 0.6, 0.2 \ 0.1, 0.7 \ 0.2, 0.6 \ 0.5, 0.5)
mata: TL = TA
mata: TP = I(2)
mata: TJ = .
mata: TV = I(11)
mata: TD = I(4)
mata: TW = J(4, 1, 1)
mata: TR = (0, 0, 0, 1)
mata: TQ = .
mata: TC = .
* negative controls: the same calls with valid arguments succeed
capture mata: TO = dwlsmv_rotse(TA, TL, TP, 1e-4, 2, TV, TJ)
if _rc {
    display as error "  error-path control (acov): the valid call failed, rc " _rc
    local ++nfail
}
capture mata: TO = dwlsmv_rotse_rank(TD, TW, TR, TQ, TC)
if _rc {
    display as error "  error-path control (rank): the valid call failed, rc " _rc
    local ++nfail
}
capture program drop tr_run
program tr_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(trlog)
    capture noisily `0'
    local rc = _rc
    quietly log close trlog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end
local call1 "dwlsmv_rotse(TA, TL[., 1], TP, 1e-4, 2, TV, TJ)"
local why1  "do not conform"
local call2 "dwlsmv_rotse(TA, TL, TP, 0, 2, TV, TJ)"
local why2  "eps <= 0"
local call3 "dwlsmv_rotse(TA, TL, TP, 1e-4, 2, TV[|1, 1 \ 9, 9|], TJ)"
local why3  "does not conform to J"
local call4 "dwlsmv_rotse_cov(J(3, 2, 1), 1, I(4))"
local why4  "does not conform to J"
local call5 "dwlsmv_rotse_rank(TD, J(4, 1, 0), TR, TQ, TC)"
local why5  "nonpositive w"
local call6 "dwlsmv_rotse_rank(TD, J(3, 1, 1), TR, TQ, TC)"
local why6  "do not conform"
forvalues k = 1/6 {
    tr_run "`why`k''" mata: TO = `call`k''
    display as text "  error path `k': rc " r(rc) ", reason found " r(found)
    if r(rc) != 3498 | r(found) != 1 {
        display as error "  error path `k' (`call`k''): expected rc 3498 with '`why`k''', got rc " r(rc)
        local ++nfail
    }
}

if `nfail' {
    display as error _n "test_rotse: `nfail' failure(s)"
    exit 9
}
display as result "test_rotse: R, J and the rotated covariance match the tier-1 fixture on Cases E, F, G (root and stall point); assembly, m = 1, rank function, mutation controls and error paths hold"
exit 0
