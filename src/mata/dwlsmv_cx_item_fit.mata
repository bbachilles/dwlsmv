*! dwlsmv -- conditional mode, Module 1: one item's ordered probit on a 0/1
*! covariate, its score, its Hessian block H_j, and its fit
*! (@sec-condx Item 1; @sec-condx-hess "The per-item block"; spec 4.8)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*! tests/unit/test_condx_stage1.do compares the fitted thresholds and
*! slopes, and H_j, with fixtures/expected/condx_ref_{A,D}_*.csv.
*!
*! THE MODEL (A1, a parallel shift): y*_j | x ~ N(b_j x, 1), thresholds tau_j
*! the same in both groups, so the bounds of category c are shifted by -b_j x:
*!     l = tau_{j,c} - b_j x,   u = tau_{j,c+1} - b_j x,
*! with tau_{j,0} = -inf and tau_{j,K} = +inf. theta_j = (tau_j1 .. tau_j,K-1,
*! b_j), K_j entries. The thresholds are those at x = 0 (spec 4.2).
*!
*! THE CELL PROBABILITY uses the upper-tail form when l > 0 (adopted
*! 2026-10-25; prereg.md, the phi2 clause, "Univariate P"):
*!     P = Phi(-l) - Phi(-u)   if l > 0,   else   P = Phi(u) - Phi(l).
*! Mata has no infinity, so the infinite bounds are branched explicitly
*! (Phi(+inf) = 1, Phi(-inf) = 0, phi(+-inf) = phi'(+-inf) = 0), never formed
*! as a missing value and never as Inf * 0.
*! The R reference's Stage 1 solver uses Phi(u) - Phi(l) throughout; its H_j
*! uses the branch (prereg [V4], accepted; far below the s tier there).
*!
*! THE SCORE AND ITS DERIVATIVE, per cell (c, x), closed form
*! (@sec-condx-hess; R: cxh_uni_dg):
*!     a_l = e_c - x e_K,  a_u = e_{c+1} - x e_K  (e_0, e_K inside a for an
*!            infinite bound are never used: phi there is 0)
*!     g   = (phi(u) a_u - phi(l) a_l) / P
*!     dg  = (phi'(u) a_u a_u' - phi'(l) a_l a_l') / P - g g'
*!     H_j = - sum_cells (n_cell / N) dg       (positive definite at the MLE)
*! The data enter only through the 2 K_j cell counts n_{c,x}; the ordering of
*! the cells is x-major: cell index x K_j + c + 1.
*!
*! THE FIT. Newton on the analytic Hessian, theta <- theta + H^-1 gbar, from
*! the pooled thresholds and b = 0 (the R reference's cx_uni_fit start), with
*! step halving only when a trial step leaves the thresholds non-increasing,
*! makes an observed cell's P non-positive, or lowers the log-likelihood by
*! more than 1e-12 (1 + |ll|). The log-likelihood is concave in theta_j, so the
*! path rule cannot change the point found. Stop on the Newton step, not the
*! score (@sec-condx-itemstop): inside the loop, when max |step| <= d1, where
*! step = H^-1 gbar at the current iterate; at iterate 50, accept if
*! max |step| <= d2; else refuse at rc 3498. A step or score with a missing
*! entry never passes (dwlsmv_cx_item_stop). d1 = 1e-12 and d2 = 1e-10: the form
*! of the certificate, (1 + h)(d + beta gamma), is derived; the two decades are
*! CHOSEN (1e-4 and 1e-2 of the 1e-8 tier), not derived. The point returned is
*! the iterate tested, theta_k, not theta_k + step.
*!
*! LIMITATION (kept open): the rule is measured in Mata only on items with
*! ||H^-1|| near 9.5 (condx_A, condx_D). For weakly informative items the
*! reachable range of d1 (||H^-1|| up to about 600) is a model fitted to R's
*! arithmetic, not a Mata result; L, gamma and xi are measured premises.
*!
*! `tolc` (internal testing only; production omits it): (d1, d2) in place of
*! the constants, so a test can reach the iterate-50 branch and the refusal.
*! Approved by the user 2026-11-11 as a test hook only: no production caller
*! passes it, and it is refused (rc 3498) unless it is 1 x 2 with no missing.
*!
*! `ctl` (internal testing only; production passes 0): 3 drops H_j's
*! (tau, b) cross entries, the derivation's drop_h = "slope" negative control.
*!
*! Requires dwlsmv_catcounts, dwlsmv_catdegen, dwlsmv_cumprop and
*! dwlsmv_thresholds.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// P of one cell. lofin / upfin: whether the lower / upper bound is finite.
// l and u are read only when finite.
real scalar dwlsmv_cx_uniP(real scalar l, real scalar u,
                           real scalar lofin, real scalar upfin)
{
    real scalar hi, lo

    if (lofin) {
        if (l > 0) {
            if (upfin) {
                return(normal(-l) - normal(-u))
            }
            return(normal(-l))
        }
    }
    hi = 1
    if (upfin) {
        hi = normal(u)
    }
    lo = 0
    if (lofin) {
        lo = normal(l)
    }
    return(hi - lo)
}

// The gradient of a bound in theta_j (K x 1): e_cat when 1 <= cat <= K-1,
// and -x in the slope position. R: cxh_grad(cat, K, x). The entries are 0
// and +-1 because x is 0/1, which the tau_U bound assumes.
real colvector dwlsmv_cx_grad(real scalar cat, real scalar K, real scalar xg)
{
    real colvector a

    a = J(K, 1, 0)
    if (cat >= 1) {
        if (cat <= K - 1) {
            a[cat] = 1
        }
    }
    a[K] = -xg
    return(a)
}

// The 2K cell counts of one item, x-major: n[x K + c + 1]. y and x are
// assumed validated (dwlsmv_cx_check).
real colvector dwlsmv_cx_item_counts(real colvector y, real colvector x,
                                     real scalar K)
{
    real colvector n
    real scalar i

    n = J(2 * K, 1, 0)
    for (i = 1; i <= rows(y); i++) {
        n[x[i] * K + y[i] + 1] = n[x[i] * K + y[i] + 1] + 1
    }
    return(n)
}

// Everything about one item at theta_j, from its 2K cell counts n.
// Outputs (overwritten): P (2K x 1), G (2K x K, row r the score g' of cell
// r), H (K x K), gbar (K x 1, the mean score), ll (the log-likelihood; missing
// if an observed cell's P is not positive). Rows of P and G for unobserved
// cells are computed but enter nothing.
void dwlsmv_cx_item_eval(real colvector th, real scalar K, real colvector n,
                         real scalar ctl, real colvector P, real matrix G,
                         real matrix H, real colvector gbar, real scalar ll)
{
    real colvector al, au, g
    real matrix D
    real scalar N, xg, c, r, lofin, upfin, l, u, fl, fu, dl, du, w

    N    = sum(n)
    P    = J(2 * K, 1, .)
    G    = J(2 * K, K, .)
    H    = J(K, K, 0)
    gbar = J(K, 1, 0)
    ll   = 0

    for (xg = 0; xg <= 1; xg++) {
        for (c = 0; c <= K - 1; c++) {
            r     = xg * K + c + 1
            lofin = (c >= 1)
            upfin = (c <= K - 2)
            l  = 0
            u  = 0
            fl = 0
            fu = 0
            dl = 0
            du = 0
            if (lofin) {
                l  = th[c] - th[K] * xg
                fl = normalden(l)
                dl = -l * fl
            }
            if (upfin) {
                u  = th[c + 1] - th[K] * xg
                fu = normalden(u)
                du = -u * fu
            }
            P[r] = dwlsmv_cx_uniP(l, u, lofin, upfin)
            al   = dwlsmv_cx_grad(c, K, xg)
            au   = dwlsmv_cx_grad(c + 1, K, xg)
            g    = (fu :* au :- fl :* al) :/ P[r]
            G[r, .] = g'

            if (n[r] > 0) {
                if (missing(P[r])) {
                    ll = .
                }
                else if (P[r] <= 0) {
                    ll = .
                }
                else {
                    ll = ll + n[r] * ln(P[r])
                }
                D = (du :* (au * au') :- dl :* (al * al')) :/ P[r] :- g * g'
                if (ctl == 3) {
                    D[|K, 1 \ K, K - 1|] = J(1, K - 1, 0)
                    D[|1, K \ K - 1, K|] = J(K - 1, 1, 0)
                }
                w    = n[r] / N
                H    = H :- w :* D
                gbar = gbar :+ w :* g
            }
        }
    }
}

// Thresholds strictly increasing and nothing missing.
real scalar dwlsmv_cx_item_ok(real colvector th, real scalar K)
{
    real scalar c

    if (missing(th)) {
        return(0)
    }
    for (c = 2; c <= K - 1; c++) {
        if (!(th[c] > th[c - 1])) {
            return(0)
        }
    }
    return(1)
}

// Convergence test on the Newton step: 1 iff neither gbar nor step has a
// missing entry and max |step| <= tol, else 0. max() skips missing values
// and sum() counts them as 0, so the missing checks are what make a missing
// entry refuse convergence; max(abs(step)) alone cannot see one.
real scalar dwlsmv_cx_item_stop(real colvector gbar, real colvector step,
                                real scalar tol)
{
    if (missing(gbar)) {
        return(0)
    }
    if (missing(step)) {
        return(0)
    }
    if (max(abs(step)) <= tol) {
        return(1)
    }
    return(0)
}

// The fit of item j (j is for messages only). Returns theta_j (K x 1);
// it is the number of Newton iterations taken. Aborts at rc 3498 if the
// information is not positive definite or the fit does not converge.
real colvector dwlsmv_cx_item_fit(real colvector y, real colvector x,
                                  real scalar K, real scalar j,
                                  real scalar it, | real matrix tolc)
{
    real colvector n, th, thn, P, Pn, gbar, gn, step
    real matrix G, Gn, H, Hn
    real scalar ll, lln, lam, h, acc, gmax, smax, d1, d2

    // @sec-condx-itemstop: the decades are chosen, the certificate derived.
    d1 = 1e-12
    d2 = 1e-10
    if (args() >= 6) {
        if (rows(tolc) != 1) {
            _error("dwlsmv_cx_item_fit: tolc must be a 1 x 2 row vector")
        }
        if (cols(tolc) != 2) {
            _error("dwlsmv_cx_item_fit: tolc must be a 1 x 2 row vector")
        }
        if (missing(tolc)) {
            _error("dwlsmv_cx_item_fit: tolc has a missing value")
        }
        d1 = tolc[1]
        d2 = tolc[2]
    }

    n  = dwlsmv_cx_item_counts(y, x, K)
    th = dwlsmv_thresholds(y, K) \ 0
    dwlsmv_cx_item_eval(th, K, n, 0, P, G, H, gbar, ll)
    if (missing(ll) | missing(gbar) | missing(H)) {
        errprintf("dwls_sample_stats_cx: item y%s cannot be evaluated at its starting values.\n",
                  strofreal(j))
        _error("dwls_sample_stats_cx: an item fails at its starting values")
    }

    for (it = 0; it < 50; it++) {
        step = cholsolve(H, gbar)
        if (missing(step)) {
            errprintf("dwls_sample_stats_cx: item y%s: the information is not positive definite\n",
                      strofreal(j))
            errprintf("  at iteration %s of its fit.\n", strofreal(it))
            _error("dwls_sample_stats_cx: an item's information is not positive definite")
        }
        if (dwlsmv_cx_item_stop(gbar, step, d1)) {
            break
        }
        lam = 1
        acc = 0
        for (h = 0; h < 30; h++) {
            thn = th :+ lam :* step
            if (dwlsmv_cx_item_ok(thn, K)) {
                dwlsmv_cx_item_eval(thn, K, n, 0, Pn, Gn, Hn, gn, lln)
                if (!missing(lln)) {
                    if (!missing(gn) & !missing(Hn)) {
                        if (lln >= ll - 1e-12 * (1 + abs(ll))) {
                            acc = 1
                            break
                        }
                    }
                }
            }
            lam = lam / 2
        }
        if (!acc) {
            errprintf("dwls_sample_stats_cx: item y%s: no acceptable step at iteration %s.\n",
                      strofreal(j), strofreal(it))
            _error("dwls_sample_stats_cx: an item's fit cannot move")
        }
        th   = thn
        P    = Pn
        G    = Gn
        H    = Hn
        gbar = gn
        ll   = lln
    }

    // The loop ran to 50 without a break (a break leaves it < 50): test the
    // step of the last iterate, from its own H and gbar, against d2.
    if (it == 50) {
        step = cholsolve(H, gbar)
        if (missing(step)) {
            errprintf("dwls_sample_stats_cx: item y%s: the information is not positive definite\n",
                      strofreal(j))
            errprintf("  at its last iterate.\n")
            _error("dwls_sample_stats_cx: an item's information is not positive definite")
        }
        if (!dwlsmv_cx_item_stop(gbar, step, d2)) {
            gmax = max(abs(gbar))
            smax = max(abs(step))
            errprintf("dwls_sample_stats_cx: item y%s did not converge: max |step| %s, max |mean score| %s after %s iterations.\n",
                      strofreal(j), strofreal(smax, "%9.2e"),
                      strofreal(gmax, "%9.2e"), strofreal(it))
            _error("dwls_sample_stats_cx: an item's fit did not converge")
        }
    }
    return(th)
}

end
