*! dwlsmv -- conditional mode, Module 2: one pair's residual polychoric given
*! a 0/1 covariate, and its Hessian-form blocks U_rho and U_theta
*! (@sec-condx Item 2; @sec-condx-hess "The pair blocks"; @sec-condx-policy
*! Result P2; spec 4.8; prereg.md section 1, the phi2 clause)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*! tests/unit/test_condx_stage1.do compares rho-hat with the
*! reference, prints the blocks against it, and asserts the phi2 copy guard
*! and the exact zeros of the x = 0 slope columns.
*!
*! THE ESTIMATOR (A2). One rho_ij for both groups; the item fits HELD FIXED.
*! The pair log-likelihood pools the two groups' K_i x K_j tables, each with
*! its own shifted bounds a_x = tau_i - b_i x and b_x = tau_j - b_j x. The
*! Module 2 routines take one threshold vector per variable and no x
*! (prereg.md: "The conditional caller passes bounds already shifted").
*!
*! THE ROOT. Module 2's derived root-finder, dwlsmv_rootfind(), unchanged,
*! with B = 0.999 and the endpoint cap 54 as in v1, so the status table and
*! the v1 entry rule apply (decided 2026-10-25, P2: "Detection reuses Module
*! 2's status codes"). It is handed dwlsmv_cx_rf_eval(), which evaluates each
*! group with v1's dwlsmv_rf_eval() and pools them, weighted by N_x / N:
*!     S/N = sum_x (N_x/N) (S_x/N_x),  c = sum_x (N_x/N) c_x,
*!     tau_S = sum_x (N_x/N) tau_S,x + 3 eps sum_x (N_x/N) |S_x/N_x|,
*! usable iff both groups are usable. Derived in @sec-condx-rootfind: the
*! usable rule is @sec-rootfind's applied to the stacked table (Result C2);
*! tau_S bounds the pooled score's error given each group's tau_S, to first
*! order, and given the tau_x to within a factor (1-u)^-4 (Result C1; the
*! 3 eps charges 2u per rounding, as the single-table tau_S does). Premises
*! measured, not proved, for the pooled score: H/c >= 1/2, the usable set
*! beyond [-1/2, 1/2] (@sec-condx-rootfind, Result C3). The uniqueness of
*! the pooled root is not established.
*! THE EVALUATION BUDGET is 200 in this mode (the user, 2026-11-13); v1 keeps
*! 40. It is a margin chosen after measuring over a conditional bound of 180
*! evaluations (@sec-condx-rootcap), not a derived bound: at 40, nine
*! constructed pairs with a root stopped at status 7.
*! No O6 repair: lavaan applies none when x is present (P2).
*!
*! THE BLOCKS, per cell (c_i, c_j, x), closed form (R: cxh_pair_d):
*!     u = P_rho / P,  v = du/drho = P_rhorho / P - u^2,
*!     w = du/dtheta' = dP_rho / P - u dP / P     (theta = (theta_i, theta_j))
*!     U_rho = sum_obs (n/N) v,  U_theta = sum_obs (n/N) w
*! THE FORM IS BINDING (prereg.md, the phi2 clause; tau_U, cxh_cell_err):
*!   * phi2 is the plain density, exponent N = a^2 - 2 rho a b + b^2, written
*!     exactly as dwlsmv_dcellprob_rho() writes it; Phi2 is binormal(), via
*!     dwlsmv_cellprob_raw();
*!   * each of P, P_rho, P_rhorho, dP, dP_rho is ONE four-term sum of
*!     per-corner terms in the production order
*!         C(a2,b2) - C(a1,b2) - C(a2,b1) + C(a1,b1);
*!     P and P_rho are dwlsmv_cellprob_raw() and dwlsmv_dcellprob_rho()
*!     themselves; the other three are formed below in the same order;
*!   * the per-corner terms are the products tau_U counts:
*!         P_rhorho:  f * g,          g = (rho (1 - Q) + a b) / s2
*!         dP_rho:    (f * h) * A,    h_i = -(a - rho b) / s2
*!         dP:        t * A,          t_i = phi(a) Phi((b - rho a) / s)
*!     with s2 = 1 - rho rho, s = sqrt(s2), Q = N / s2, and A, B the bound
*!     gradients, entries 0 and +-1 (x is 0/1), so * A is exact.
*! A different factoring is outside the bound and reopens gate step 3.
*!
*! `ctl` on dwlsmv_cx_pair_blocks() is for internal testing only; production
*! passes 0. 1 zeroes U_theta's slope columns (R cx_psi, prop = "tau"); 2
*! zeroes U_theta (prop = "none"); 4 drops the quotient term -u dP / P
*! (cxh_pair_d, drop = "quot").
*!
*! Requires dwlsmv_cellprob, dwlsmv_dcellprob_rho, dwls_polychoric (for
*! dwlsmv_rootfind and dwlsmv_rf_eval), dwlsmv_crosstab, dwlsmv_catcounts and
*! dwlsmv_cx_item (dwlsmv_cx_grad).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The phi2 corner grid of one group, (K1+1) x (K2+1): row 1 and row K1+1 are
// the infinite bounds, where the density is 0; likewise the columns. The
// expression is dwlsmv_dcellprob_rho()'s, character for character but for
// the names, so that its four-corner difference equals that routine's output
// bit for bit (asserted in test_condx_stage1.do, check 2).
real matrix dwlsmv_cx_phi2grid(real colvector a, real colvector b, real scalar rho)
{
    real matrix F
    real scalar K1, K2, r, s, sig, twopi, p, q

    K1    = rows(a) + 1
    K2    = rows(b) + 1
    sig   = sqrt(1 - rho * rho)
    twopi = 2 * pi()
    F = J(K1 + 1, K2 + 1, 0)
    for (r = 2; r <= K1; r++) {
        for (s = 2; s <= K2; s++) {
            p       = a[r - 1]
            q       = b[s - 1]
            F[r, s] = exp(-(p * p - 2 * rho * p * q + q * q) /
                          (2 * (1 - rho * rho))) / (twopi * sig)
        }
    }
    return(F)
}

// The pooled score evaluator, in the shape dwlsmv_rootfind() calls:
// (x, n, a, b) -> (x, usable, S/N, c, tau_S). The two groups arrive stacked:
// n = (n_0 \ n_1), 2 K1 x K2; a = (a_0 \ a_1), 2 C1 x 1; b = (b_0 \ b_1).
real rowvector dwlsmv_cx_rf_eval(real scalar x, real matrix n,
                                 real colvector a, real colvector b)
{
    real matrix    n0, n1
    real colvector a0, a1, b0, b1
    real rowvector e0, e1
    real scalar    K1, C1, C2, N0, N1, w0, w1, eps

    K1 = rows(n) / 2
    C1 = rows(a) / 2
    C2 = rows(b) / 2
    n0 = n[|1, 1 \ K1, cols(n)|]
    n1 = n[|K1 + 1, 1 \ 2 * K1, cols(n)|]
    a0 = a[|1 \ C1|]
    a1 = a[|C1 + 1 \ 2 * C1|]
    b0 = b[|1 \ C2|]
    b1 = b[|C2 + 1 \ 2 * C2|]
    N0 = sum(n0)
    N1 = sum(n1)
    w0 = N0 / (N0 + N1)
    w1 = N1 / (N0 + N1)
    eps = epsilon(1)

    e0 = dwlsmv_rf_eval(x, n0, a0, b0)
    e1 = dwlsmv_rf_eval(x, n1, a1, b1)
    if (!e0[2] | !e1[2]) {
        return((x, 0, ., ., .))
    }
    return((x, 1, w0 * e0[3] + w1 * e1[3], w0 * e0[4] + w1 * e1[4],
            w0 * e0[5] + w1 * e1[5] + 3 * eps * (w0 * abs(e0[3]) + w1 * abs(e1[3]))))
}

// The pair's 2 K_i K_j cell counts, x-major then row-major:
// n[x K_i K_j + c_i K_j + c_j + 1]. Inputs assumed validated.
real colvector dwlsmv_cx_pair_counts(real colvector yi, real colvector yj,
                                     real colvector x, real scalar Ki,
                                     real scalar Kj)
{
    real colvector n
    real scalar i, r

    n = J(2 * Ki * Kj, 1, 0)
    for (i = 1; i <= rows(yi); i++) {
        r    = x[i] * Ki * Kj + yi[i] * Kj + yj[i] + 1
        n[r] = n[r] + 1
    }
    return(n)
}

// rho-hat for the pair (thi, thj = theta_i, theta_j: thresholds then slope).
// Returns (rho-hat, status, evaluations), as dwlsmv_rootfind().
real rowvector dwlsmv_cx_pair_fit(real colvector yi, real colvector yj,
                                  real colvector x, real colvector thi,
                                  real colvector thj, | real rowvector ik)
{
    real matrix    n0, n1
    real colvector ti, tj, ai, bj
    real scalar    Ki, Kj
    string scalar  who

    Ki = rows(thi)
    Kj = rows(thj)
    ti = thi[|1 \ Ki - 1|]
    tj = thj[|1 \ Kj - 1|]
    n0 = dwlsmv_crosstab(select(yi, x :== 0), select(yj, x :== 0), Ki, Kj)
    n1 = dwlsmv_crosstab(select(yi, x :== 1), select(yj, x :== 1), Ki, Kj)

    // Bounds shifted per group, b x: at x = 0 the shift is 0.
    ai = (ti \ (ti :- thi[Ki]))
    bj = (tj \ (tj :- thj[Kj]))

    // ik = (i, k), optional, for the message only. dwlsmv_rootfind() evaluates
    // dwlsmv_cx_rf_eval() at rho = 0 first and refuses if that point is unusable;
    // it is unusable iff either group is, and a group is unusable iff
    // dwlsmv_rf_usable() says so (the test dwlsmv_rf_eval() applies). Each line
    // below is printed only then, so only when the call that follows refuses.
    who = "a pair"
    if (args() >= 6) {
        if (cols(ik) == 2) {
            who = "y" + strofreal(ik[1]) + "~~y" + strofreal(ik[2])
        }
    }
    if (!dwlsmv_rf_usable(0, n0, ti, tj)) {
        errprintf("dwls_sample_stats_cx: the pair %s has no starting point: rho = 0 is not usable in the group x = 0.\n",
                  who)
    }
    if (!dwlsmv_rf_usable(0, n1, ti :- thi[Ki], tj :- thj[Kj])) {
        errprintf("dwls_sample_stats_cx: the pair %s has no starting point: rho = 0 is not usable in the group x = 1.\n",
                  who)
    }

    return(dwlsmv_rootfind(&dwlsmv_cx_rf_eval(), (n0 \ n1), ai, bj,
                           (., .), 0.999, 200, 54))
}

// Per-cell table at (theta_i, theta_j, rho), 2 K_i K_j rows in the counts'
// order, columns (P, u, v, w_1 .. w_{K_i + K_j}). Cells whose P is not
// positive come back with whatever the arithmetic gives; callers use only
// observed cells and check them.
real matrix dwlsmv_cx_pair_tab(real colvector thi, real colvector thj,
                               real scalar rho, real scalar ctl)
{
    real matrix    T, Pm, Prm, F, Gg, Hi, Hj, Ti, Tj
    real colvector a, b, A, B, A1, A2, B1, B2, dPi, dPj, dPri, dPrj, w
    real scalar    Ki, Kj, xg, r, s, row, s2, sd, p, q, f, Nn, Q
    real scalar    P, Pr, Prr, u, v, afin, bfin

    Ki = rows(thi)
    Kj = rows(thj)
    s2 = 1 - rho * rho
    sd = sqrt(1 - rho * rho)
    T  = J(2 * Ki * Kj, 3 + Ki + Kj, .)

    for (xg = 0; xg <= 1; xg++) {
        a   = thi[|1 \ Ki - 1|] :- thi[Ki] * xg
        b   = thj[|1 \ Kj - 1|] :- thj[Kj] * xg
        Pm  = dwlsmv_cellprob_raw(a, b, rho)
        Prm = dwlsmv_dcellprob_rho(a, b, rho)
        F   = dwlsmv_cx_phi2grid(a, b, rho)

        // Corner grids of the per-corner terms; 0 wherever the term is 0
        // (an infinite bound kills phi2, and phi(a) needs a finite a).
        A  = (. \ a \ .)
        B  = (. \ b \ .)
        Gg = J(Ki + 1, Kj + 1, 0)
        Hi = J(Ki + 1, Kj + 1, 0)
        Hj = J(Ki + 1, Kj + 1, 0)
        Ti = J(Ki + 1, Kj + 1, 0)
        Tj = J(Ki + 1, Kj + 1, 0)
        for (r = 1; r <= Ki + 1; r++) {
            for (s = 1; s <= Kj + 1; s++) {
                afin = (r >= 2 & r <= Ki)
                bfin = (s >= 2 & s <= Kj)
                if (afin & bfin) {
                    p  = A[r]
                    q  = B[s]
                    f  = F[r, s]
                    Nn = p * p - 2 * rho * p * q + q * q
                    Q  = Nn / s2
                    Gg[r, s] = f * ((rho * (1 - Q) + p * q) / s2)
                    Hi[r, s] = f * (-(p - rho * q) / s2)
                    Hj[r, s] = f * (-(q - rho * p) / s2)
                }
                // t_i = phi(a) Phi((b - rho a)/s); Phi is 0 at b = -inf and
                // 1 at b = +inf.
                if (afin) {
                    p = A[r]
                    if (s == Kj + 1) {
                        Ti[r, s] = normalden(p)
                    }
                    else if (s >= 2) {
                        Ti[r, s] = normalden(p) * normal((B[s] - rho * p) / sd)
                    }
                }
                if (bfin) {
                    q = B[s]
                    if (r == Ki + 1) {
                        Tj[r, s] = normalden(q)
                    }
                    else if (r >= 2) {
                        Tj[r, s] = normalden(q) * normal((A[r] - rho * q) / sd)
                    }
                }
            }
        }

        // Cells. Corner (k, l): rows r (a1) and r+1 (a2), columns s (b1) and
        // s+1 (b2). Every sum is C(a2,b2) - C(a1,b2) - C(a2,b1) + C(a1,b1).
        for (r = 1; r <= Ki; r++) {
            A1 = dwlsmv_cx_grad(r - 1, Ki, xg)
            A2 = dwlsmv_cx_grad(r, Ki, xg)
            for (s = 1; s <= Kj; s++) {
                B1  = dwlsmv_cx_grad(s - 1, Kj, xg)
                B2  = dwlsmv_cx_grad(s, Kj, xg)
                row = xg * Ki * Kj + (r - 1) * Kj + s

                P   = Pm[r, s]
                Pr  = Prm[r, s]
                Prr = Gg[r + 1, s + 1] - Gg[r, s + 1] - Gg[r + 1, s] + Gg[r, s]
                dPi  = Ti[r + 1, s + 1] :* A2 :- Ti[r, s + 1] :* A1 :-
                       Ti[r + 1, s] :* A2 :+ Ti[r, s] :* A1
                dPj  = Tj[r + 1, s + 1] :* B2 :- Tj[r, s + 1] :* B2 :-
                       Tj[r + 1, s] :* B1 :+ Tj[r, s] :* B1
                dPri = Hi[r + 1, s + 1] :* A2 :- Hi[r, s + 1] :* A1 :-
                       Hi[r + 1, s] :* A2 :+ Hi[r, s] :* A1
                dPrj = Hj[r + 1, s + 1] :* B2 :- Hj[r, s + 1] :* B2 :-
                       Hj[r + 1, s] :* B1 :+ Hj[r, s] :* B1

                u = Pr / P
                v = Prr / P - u * u
                if (ctl == 4) {
                    w = (dPri \ dPrj) :/ P
                }
                else {
                    w = (dPri \ dPrj) :/ P :- u :* (dPi \ dPj) :/ P
                }
                T[row, .] = (P, u, v, w')
            }
        }
    }
    return(T)
}

// The pair's blocks from its observed cells. Outputs (overwritten): Ur
// (U_rho), Ut (U_theta, 1 x (K_i + K_j)), tab (dwlsmv_cx_pair_tab()) and n
// (the cell counts). Aborts if an observed cell has no positive P or a
// missing value: at an accepted root every observed cell is usable.
// Names the pair and the observed cell of row r (the counts' order: x-major, then
// row-major, r = x Ki Kj + c_i Kj + c_j + 1, categories 0-based) for a refusal of
// dwlsmv_cx_pair_blocks(). ik = (i, k) or empty. Called only on the line before
// an _error().
void dwlsmv_cx_badcell(real scalar r, real scalar Ki, real scalar Kj,
                       real rowvector ik, string scalar what)
{
    real scalar   xg, ci, cj, rem
    string scalar who, ni, nk

    rem = r - 1
    xg  = floor(rem / (Ki * Kj))
    rem = rem - xg * Ki * Kj
    ci  = floor(rem / Kj)
    cj  = rem - ci * Kj
    who = "a pair"
    ni  = "item i"
    nk  = "item k"
    if (cols(ik) == 2) {
        who = "y" + strofreal(ik[1]) + "~~y" + strofreal(ik[2])
        ni  = "y" + strofreal(ik[1])
        nk  = "y" + strofreal(ik[2])
    }
    errprintf("dwlsmv_cx_pair_blocks: %s: %s at the observed cell (x = %s, %s = %s, %s = %s).\n",
              who, what, strofreal(xg), ni, strofreal(ci), nk, strofreal(cj))
}

void dwlsmv_cx_pair_blocks(real colvector yi, real colvector yj,
                           real colvector x, real colvector thi,
                           real colvector thj, real scalar rho,
                           real scalar ctl, real scalar Ur,
                           real rowvector Ut, real matrix tab,
                           real colvector n, | real rowvector ik)
{
    real scalar Ki, Kj, N, r, wt
    real rowvector id

    id = J(1, 0, .)
    if (args() >= 12) {
        id = ik
    }
    Ki  = rows(thi)
    Kj  = rows(thj)
    N   = rows(yi)
    n   = dwlsmv_cx_pair_counts(yi, yj, x, Ki, Kj)
    tab = dwlsmv_cx_pair_tab(thi, thj, rho, ctl)
    Ur  = 0
    Ut  = J(1, Ki + Kj, 0)
    for (r = 1; r <= rows(n); r++) {
        if (n[r] == 0) {
            continue
        }
        if (missing(tab[r, .])) {
            dwlsmv_cx_badcell(r, Ki, Kj, id, "a missing value")
            _error("dwlsmv_cx_pair_blocks: an observed cell has a missing value")
        }
        if (tab[r, 1] <= 0) {
            dwlsmv_cx_badcell(r, Ki, Kj, id, "no positive probability")
            _error("dwlsmv_cx_pair_blocks: an observed cell has no positive probability")
        }
        wt = n[r] / N
        Ur = Ur + wt * tab[r, 3]
        Ut = Ut :+ wt :* tab[|r, 4 \ r, 3 + Ki + Kj|]
    }
    if (ctl == 1) {
        Ut[Ki]      = 0
        Ut[Ki + Kj] = 0
    }
    if (ctl == 2) {
        Ut = J(1, Ki + Kj, 0)
    }
}

end
