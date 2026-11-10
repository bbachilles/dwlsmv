*! dwlsmv -- one pair's polychoric correlation and its information quantities
*! (Module 2, Gate 3; spec Sections 4.5.2, 5 Module 2 and 5.1)
*!
*! First executed under StataMP 18 on 2026-09-19 with optimize() and a
*! Fisher-scoring polish. Since 2026-09-30 the root is found by option B, the
*! bracketed root-finder of derivation @sec-rootfind: the estimator is
*! unchanged, only the algorithm that finds its root.
*!
*! Two-stage, following Muthen (1984): Module 1's thresholds arrive already
*! estimated and are HELD FIXED here while rho is estimated from the pairwise
*! likelihood. That the thresholds are fixed at a random quantity is exactly
*! what creates the propagation term of spec 4.5.2, and it is why this routine
*! returns the information quantities alongside rho rather than rho alone.
*!
*! WHAT THIS RETURNS, and what it deliberately does not. Gate 3 settled that
*! Module 2 returns rho-hat, u_rs, I_rho_rho, I_rho_a and I_rho_b per pair, and
*! that Module 4 combines them with the Stage 1 staircase to form psi_rho. The
*! influence function is NOT built here. Since 2026-09-30 the pair also carries
*! `status` (spec Module 2's table, 1 = root), and `converged` is set from it.
*!
*! NO NUMERICAL SECOND DERIVATIVES, anywhere. The information quantities are
*! closed-form in the cell probabilities and their FIRST derivatives, computed
*! by dwlsmv_pairinfo() at the OBSERVED cell proportions (the empirical
*! information lavaan uses; derivation @sec-normaliser).
*!
*! THE ROOT-FINDER (derivation @sec-rootfind; spec Module 2). rho-hat is the
*! root of the observed-cells score S(rho) = sum_O n_rs u_rs. Only USABLE
*! points are used for a decision: every observed cell's computed pi above
*! delta_P, the absolute error of a computed cell probability. From rho = 0 the
*! walk halves the distance to +-1 in the direction of S(0), capped at +-B, until
*! the score changes sign; a safeguarded Fisher loop, scaled by the masked
*! expected information, shrinks that sign-certified bracket until BOTH
*! |S|/N <= tau_S and width <= tau_rho; the endpoint check then requires the
*! right strict sign at the last usable point on each side. Every other exit is
*! a distinct status. It replaces optimize(), whose d0 step collapsed on 21
*! high-rho pairs (E-008, 2026-09-28).
*!
*! The constants: B = 0.999 (a declared hard search limit); delta_Phi, the
*! per-corner error budget of binormal(), CHOSEN AFTER MEASURING on StataMP 18
*! and 16.1 (tools/measure_binormal.do); the root search's cap of 40, CHOSEN
*! AFTER MEASURING; the endpoint check's cap of 54 per side, DERIVED from its
*! four-ulp frontier width.
*!
*! Requires dwlsmv_types, dwlsmv_cellprob, dwlsmv_dcellprob_rho, dwlsmv_wmat,
*! dwlsmv_pairinfo, dwlsmv_crosstab and dwlsmv_catcounts to be loaded.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The spacing of doubles at x (x = 0 maps to the smallest normal's). Mata has
// no frexp: take the exponent from logs and correct it by exact powers of 2.
real scalar dwlsmv_rf_ulp(real scalar x)
{
    real scalar ax, e

    ax = abs(x)
    if (ax < 2^-1022) ax = 2^-1022
    e = floor(ln(ax) / ln(2))
    if (2^e > ax) e = e - 1
    if (2^(e + 1) <= ax) e = e + 1
    return(2^(e - 52))
}

// delta_P = 4 delta_Phi + 3 eps: four corners each within delta_Phi, and three
// roundings of numbers at most 1 (derivation @sec-rootfind, the stopping rule).
real scalar dwlsmv_rf_dP()
{
    return(4 * 9.3908313968469e-16 + 3 * epsilon(1))
}

// One evaluation at x: (x, usable, S/N, c, tau_S). Unusable returns (x, 0, ., ., .).
// Usable: every OBSERVED cell's computed pi exceeds delta_P. S/N sums over the
// observed cells, with the mask applied before the division; c is the expected
// information over cells whose computed pi is positive; tau_S is the score's
// own first-order evaluation error. The derivation's rf_eval().
real rowvector dwlsmv_rf_eval(real scalar x, real matrix n,
                              real colvector a, real colvector b)
{
    return(dwlsmv_rf_evalp(x, n, a, b, &dwlsmv_dcellprob_rho()))
}

// The same with the derivative routine as an argument, as the derivation's
// rf_eval(dS =): production passes dwlsmv_dcellprob_rho(); only the tests'
// negative controls pass another.
real rowvector dwlsmv_rf_evalp(real scalar x, real matrix n,
                               real colvector a, real colvector b,
                               pointer(real matrix function) scalar dS)
{
    real matrix   P, D, O, Po, pos, Pp, u, w, W
    real colvector A, Bt
    real scalar   dP, eps, K1, K2, r, s, p, q, z, tauS

    dP  = dwlsmv_rf_dP()
    eps = epsilon(1)
    P   = dwlsmv_cellprob_raw(a, b, x)
    O   = (n :> 0)
    Po  = P :* O :+ (1 :- O)
    if (missing(Po) | min(Po) <= dP) {
        return((x, 0, ., ., .))
    }
    D   = (*dS)(a, b, x)
    w   = n :/ sum(n)
    u   = (D :/ Po) :* O
    pos = (P :> 0)
    Pp  = P :* pos :+ (1 :- pos)

    // Per cell: sum over its four corners of |phi2| (z + 7), z the density's
    // exponent; corners at infinity weigh 0 (the derivation's rf_absD).
    K1 = rows(a) + 1
    K2 = rows(b) + 1
    A  = (. \ a \ .)
    Bt = (. \ b \ .)
    W  = J(K1 + 1, K2 + 1, 0)
    for (r = 2; r <= K1; r++) {
        for (s = 2; s <= K2; s++) {
            p = A[r]
            q = Bt[s]
            z = (p * p - 2 * x * p * q + q * q) / (2 * (1 - x * x))
            W[r, s] = exp(-z) / (2 * pi() * sqrt(1 - x * x)) * (z + 7)
        }
    }
    W = W[|2, 2 \ K1 + 1, K2 + 1|] + W[|1, 2 \ K1, K2 + 1|] +
        W[|2, 1 \ K1 + 1, K2|] + W[|1, 1 \ K1, K2|]

    tauS = sum(w :* (eps :* W :+ abs(u) :* dP) :/ Po) +
           (sum(O) + 2) * eps * sum(w :* abs(u))
    return((x, 1, sum(w :* u), sum(((D :* D) :/ Pp) :* pos), tauS))
}

// tau_rho = 16 tau_S / c + 8 ulp(x), at an evaluation e.
real scalar dwlsmv_rf_taur(real rowvector e)
{
    return(16 * e[5] / e[4] + 8 * dwlsmv_rf_ulp(e[1]))
}

// The walk (and, if a trial is unusable, the frontier bisection) on side d.
// Returns the inner point in ein and, on a bracket, the outer point in eout;
// k counts evaluations. The value is 0 for a bracket, else the status:
// 2 no sign change up to the bound, 3 up to the usable end, 4 gap.
real scalar dwlsmv_rf_side(pointer(real rowvector function) scalar fe,
                           real matrix n, real colvector a, real colvector b,
                           real rowvector ein, real scalar d, real scalar B,
                           real rowvector eout, real scalar k)
{
    real rowvector e, eB
    real scalar    cand, out, m

    k = 0
    while (1) {
        cand = d * min((B, 1 - (1 - d * ein[1]) / 2))
        if (cand == ein[1]) {
            return(2)
        }
        e = (*fe)(cand, n, a, b)
        k = k + 1
        if (e[2] & d * e[3] < 0) {
            eout = e
            return(0)
        }
        if (e[2]) {
            ein = e
            continue
        }
        if (cand != d * B) {
            eB = (*fe)(d * B, n, a, b)
            k  = k + 1
            if (eB[2]) {
                return(4)
            }
        }
        out = cand
        while (abs(out - ein[1]) > 4 * dwlsmv_rf_ulp(out)) {
            m = (ein[1] + out) / 2
            e = (*fe)(m, n, a, b)
            k = k + 1
            if (!e[2]) {
                out = m
                continue
            }
            if (d * e[3] < 0) {
                eout = e
                return(0)
            }
            ein = e
        }
        return(3)
    }
}

// Usability alone at x: exactly dwlsmv_rf_evalp()'s test (every OBSERVED
// cell's computed pi exceeds delta_P), without the score, information or tau_S.
real scalar dwlsmv_rf_usable(real scalar x, real matrix n,
                             real colvector a, real colvector b)
{
    real matrix P, O, Po

    P  = dwlsmv_cellprob_raw(a, b, x)
    O  = (n :> 0)
    Po = P :* O :+ (1 :- O)
    return(!(missing(Po) | min(Po) <= dwlsmv_rf_dP()))
}

// The last usable evaluation on side d: d * B if usable, else the frontier
// bisection from ein, to four ulps. capc is the endpoint check's own cap per
// side (54, derived). Returns 1 if the cap was reached; the end is left in
// eend and the evaluation count in k.
//
// The bisection's decisions use only usability, so with the production
// evaluator each step tests usability alone (dwlsmv_rf_usable(), the same
// test) and one full evaluation is made at the final inner point. The points
// visited, the decisions, eend and k are unchanged; the repeat evaluation at
// an already-visited point is not counted (2026-10-05, a speed change only:
// E-089). Any other fe, such as a
// test's negative control, is called at every step, as before.
real scalar dwlsmv_rf_end(pointer(real rowvector function) scalar fe,
                          real matrix n, real colvector a, real colvector b,
                          real scalar d, real scalar B, real rowvector ein,
                          real scalar capc, real rowvector eend, real scalar k)
{
    pointer(real rowvector function) scalar prod
    real rowvector e, ei
    real scalar    out, m, xi, u, cheap, moved

    e = (*fe)(d * B, n, a, b)
    k = 1
    if (e[2]) {
        eend = e
        return(0)
    }
    prod  = &dwlsmv_rf_eval()
    cheap = (fe == prod)
    ei    = ein
    xi    = ein[1]
    moved = 0
    out   = d * B
    while (abs(out - xi) > 4 * dwlsmv_rf_ulp(out)) {
        if (k >= capc) {
            eend = ((cheap & moved) ? (*fe)(xi, n, a, b) : ei)
            return(1)
        }
        m = (xi + out) / 2
        if (cheap) {
            u = dwlsmv_rf_usable(m, n, a, b)
        }
        else {
            e = (*fe)(m, n, a, b)
            u = e[2]
        }
        k = k + 1
        if (u) {
            xi    = m
            moved = 1
            if (!cheap) {
                ei = e
            }
        }
        else {
            out = m
        }
    }
    eend = ((cheap & moved) ? (*fe)(xi, n, a, b) : ei)
    return(0)
}

// The endpoint check around a candidate root with bracket ends L and R.
// Returns the status (1 root, 5 more than one sign change, 6 ambiguous,
// 7 cap) and adds its evaluations to k.
real scalar dwlsmv_rf_check(pointer(real rowvector function) scalar fe,
                            real matrix n, real colvector a, real colvector b,
                            real scalar B, real rowvector L, real rowvector R,
                            real scalar capc, real scalar k)
{
    real rowvector lo, hi
    real scalar    clo, chi, klo, khi

    clo = dwlsmv_rf_end(fe, n, a, b, -1, B, L, capc, lo, klo)
    chi = dwlsmv_rf_end(fe, n, a, b,  1, B, R, capc, hi, khi)
    k   = k + klo + khi
    if (clo | chi) return(7)
    if (abs(lo[3]) <= lo[5] | abs(hi[3]) <= hi[5]) return(6)
    if (lo[3] > lo[5] & hi[3] < -hi[5]) return(1)
    return(5)
}

// The internal root finder (spec 5.1). fe is the score evaluator: production
// passes &dwlsmv_rf_eval(); tests may pass a modified one. bracket is (.,.) in
// production; a supplied bracket skips the walk (internal testing only).
// Returns (rho-hat, status, evaluations), the endpoint check's included.
// Broken invariants abort (rc 3498) and are not statuses.
real rowvector dwlsmv_rootfind(pointer(real rowvector function) scalar fe,
                               real matrix n, real colvector a, real colvector b,
                               real rowvector bracket, real scalar B,
                               real scalar cap, real scalar capc)
{
    real rowvector e0, pm, pp, L, R, x, ein, eout, e
    real scalar    k, ks, d, st, tr, xn, dxold, stp

    k = 0
    if (missing(bracket)) {
        e0 = (*fe)(0, n, a, b)
        k  = 1
        if (!e0[2]) {
            errprintf("dwls_polychoric: rho = 0 is not usable (an observed cell probability ")
            errprintf("is within delta_P of 0), so the pair has no starting point.\n")
            _error("dwls_polychoric: rho = 0 is not usable")
        }
        if (e0[3] == 0) {
            // An exact zero at 0: probe both sides at +-tau_rho(0)/2.
            tr = dwlsmv_rf_taur(e0)
            pm = (*fe)(-tr / 2, n, a, b)
            pp = (*fe)( tr / 2, n, a, b)
            k  = 3
            if (pm[2] & pm[3] > pm[5] & pp[2] & pp[3] < -pp[5]) {
                st = dwlsmv_rf_check(fe, n, a, b, B, pm, pp, capc, k)
                return((0, st, k))
            }
            if (pm[2] & pm[3] > pm[5] & pp[2] & pp[3] > pp[5]) {
                e0 = pp
            }
            else if (pm[2] & pm[3] < -pm[5] & pp[2] & pp[3] < -pp[5]) {
                e0 = pm
            }
            else {
                return((0, 6, k))
            }
        }
        d   = sign(e0[3])
        ein = e0
        st  = dwlsmv_rf_side(fe, n, a, b, ein, d, B, eout, ks)
        k   = k + ks
        if (st != 0) {
            return((ein[1], st, k))
        }
        if (d > 0) {
            L = ein
            R = eout
        }
        else {
            L = eout
            R = ein
        }
    }
    else {
        L = (*fe)(bracket[1], n, a, b)
        R = (*fe)(bracket[2], n, a, b)
        k = 2
        if (!L[2] | !R[2]) {
            _error("dwls_polychoric: unusable bracket end")
        }
        if (!(L[3] > 0 & R[3] < 0)) {
            _error("dwls_polychoric: bracket not maximum-oriented")
        }
    }

    x     = (abs(L[3]) < abs(R[3]) ? L : R)
    dxold = R[1] - L[1]
    while (1) {
        tr = dwlsmv_rf_taur(x)
        if (abs(x[3]) <= x[5] & R[1] - L[1] <= tr) {
            st = dwlsmv_rf_check(fe, n, a, b, B, L, R, capc, k)
            return((x[1], st, k))
        }
        if (k >= cap) {
            return((x[1], 7, k))
        }
        if (abs(x[3]) <= x[5]) {
            // Probe the far side.
            xn = (x == L ? x[1] + tr / 2 : x[1] - tr / 2)
            if (!(xn > L[1] & xn < R[1])) xn = (L[1] + R[1]) / 2
        }
        else {
            stp = x[3] / x[4]
            xn  = x[1] + stp
            if (!(xn > L[1] & xn < R[1]) | abs(stp) > abs(dxold) / 2) {
                xn  = (L[1] + R[1]) / 2
                stp = (R[1] - L[1]) / 2
            }
            dxold = stp
        }
        e = (*fe)(xn, n, a, b)
        k = k + 1
        if (!e[2]) {
            return((x[1], 8, k))
        }
        if (e[3] > 0) {
            L = e
        }
        else {
            R = e
        }
        x = (abs(L[3]) < abs(R[3]) ? L : R)
    }
}

struct dwlsmv_pairfit scalar dwls_polychoric(real colvector yi,
                                             real colvector yj,
                                             real colvector tau_i,
                                             real colvector tau_j)
{
    struct dwlsmv_pairfit scalar f
    real matrix    n
    real rowvector r
    real scalar    K1, K2, N

    if (rows(yi) != rows(yj)) {
        _error("dwls_polychoric: the two variables have different lengths")
    }
    if (cols(tau_i) != 1 | cols(tau_j) != 1) {
        _error("dwls_polychoric: thresholds must be column vectors")
    }
    if (rows(tau_i) < 1 | rows(tau_j) < 1) {
        _error("dwls_polychoric: each variable needs at least one threshold")
    }

    K1 = rows(tau_i) + 1
    K2 = rows(tau_j) + 1
    N  = rows(yi)

    n = dwlsmv_crosstab(yi, yj, K1, K2)

    // B = 0.999, the cap 40 (chosen after measuring), the endpoint cap 54 per
    // side (derived): spec Module 2.
    r = dwlsmv_rootfind(&dwlsmv_rf_eval(), n, tau_i, tau_j, (., .), 0.999, 40, 54)

    f.rho        = r[1]
    f.status     = r[2]
    f.converged  = (f.status == 1)
    f.iterations = r[3]

    // ---- information quantities, all at rho-hat -------------------------
    // Weighted by the OBSERVED cell proportions, n:/N: the empirical
    // (outer-product) information lavaan uses (dwlsmv_pairinfo(); derivation
    // @sec-normaliser). rho-hat is a usable point whatever the status, so
    // every observed cell is positive there and pairinfo cannot abort on one.
    dwlsmv_pairinfo(f, tau_i, tau_j, f.rho, n :/ N)

    // The caller knows which variables these are; this routine does not, and
    // the locked signature does not pass them. Module 4 fills them in.
    f.i = 0
    f.j = 0

    return(f)
}

end
