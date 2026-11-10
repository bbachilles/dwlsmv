*! Oblique geomin rotation (0.3.0, Unit 5): dwlsmv_geomin() against lavaan.
*!
*! Executed under StataMP 18 on 2026-10-29 and it passes. Worst ratios to tier:
*! Lambda 7.1e-3, Phi 2.2e-3, Q 3.5e-4, T3 0.13 (all Case F); sigma invariance
*! <= 2.8e-16; ||Gp|| <= 1.2e-8; gradients vs central differences <= 9.8e-10.
*! What broke first: linear subscripts on a matrix in tg_fd() (rc 3201), and
*! undeclared output arguments, which made every error path rc 3499 -- hence
*! the declarations and the valid-call control before check 9.
*!
*! Each case is fitted as an echelon EFA through v1 (as test_efa_unrotated.do:
*! m(m-1)/2 anchor zeros, unitvariance on every factor, no psi covariance, so
*! Phi = I), A is read off theta-hat, and dwlsmv_geomin() rotates it with the
*! Mplus-rule eps. The fixture is lavaan's geomin under the same conventions,
*! and lavaan's sign and "index" order ARE dwlsmv's convention, so the rotated
*! columns are compared directly, with no alignment.
*!
*! Reference: fixtures/expected/lavaan_efa_geomin.csv (99h_lavaan_efa_geomin.R,
*! lavaan 0.7.2). Cases E (m = 2), F (m = 3, data 91_efa_data.R), G (m = 2,
*! all binary).
*!
*! Tolerances, fixed before the first run (2026-10-29):
*!   rotated Lambda, Phi, Q  theta-hat tier, 1e-5 abs + 1e-4 rel (spec 7.2):
*!     they are functions of theta-hat, whose own gap to lavaan is that tier;
*!     the rotation adds lavaan's gpa_tol 1e-6 stopping error, measured at
*!     <= 3.5e-7 in @sec-geomin from a common A, inside the tier.
*!   T3 1e-6 abs + 1e-5 rel; df exact (the fit is v1 unchanged).
*!   sigma invariance under the rotation  1e-12 abs (algebra, GM1).
*!   gradient at the solution  ||Gp||_F < 1e-5, the convergence rule itself.
*!   analytic vs central-difference gradients  1e-6 relative to the largest
*!     entry (h = 1e-6, loadings bounded away from 0).
*!   determinism, nesting: exact equality; nested Q within the
*!     distinct-minimum tolerance 1e-8 + 1e-6 Q.
*!
*! Checks:
*!   1  per case: fit converges; T3, df vs lavaan; rotated Lambda, Phi and Q
*!      vs lavaan; eps in info is the Mplus rule and the fixture's
*!   2  per case: Lambda Phi Lambda' = A A'
*!   3  per case: ||Gp|| at the returned solution, recomputed here from
*!      Lambda and A (T' = (L'L)^-1 L'A), below 1e-5; info[6] below 1e-5
*!   4  per case: sign and order convention, recomputed here
*!   5  per case: info well-formed (counts, best <= converged <= starts)
*!   6  dQ/dLambda and dQ/dT against central differences, m = 2 and 3
*!   7  determinism (two calls identical) and nesting (40 starts never worse
*!      than 30, and here equal within tolerance), Case F
*!   8  m = 1: Phi = 1, positive column sum, A unchanged up to sign
*!   9  error paths, rc 3498: eps 0 and missing, nstart 0 and 2.5, A with a
*!      missing value, m > p, and no converged start (maxit 0)
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 120

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_mg_check.mata          ///
             src/mata/dwlsmv_mg_sub.mata            ///
             src/mata/dwlsmv_mg_parts.mata          ///
             src/mata/dwlsmv_mg_sigma.mata          ///
             src/mata/dwlsmv_mg_dcols.mata          ///
             src/mata/dwlsmv_id_parse.mata          ///
             src/mata/dwlsmv_id_rank.mata           ///
             src/mata/dwlsmv_id_tie.mata            ///
             src/mata/dwlsmv_id_groups.mata         ///
             src/mata/dwlsmv_id_opname.mata         ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_catdegen.mata          ///
             src/mata/dwlsmv_cumprop.mata           ///
             src/mata/dwlsmv_thresholds.mata        ///
             src/mata/dwlsmv_psi_tau_table.mata     ///
             src/mata/dwlsmv_cellprob.mata          ///
             src/mata/dwlsmv_dcellprob_rho.mata     ///
             src/mata/dwlsmv_wmat.mata              ///
             src/mata/dwlsmv_crosstab.mata          ///
             src/mata/dwlsmv_pairinfo.mata          ///
             src/mata/dwls_polychoric.mata          ///
             src/mata/dwlsmv_koff.mata              ///
             src/mata/dwlsmv_psi_tau_stack.mata     ///
             src/mata/dwlsmv_psi_rho_table.mata     ///
             src/mata/dwlsmv_psi_matrix.mata        ///
             src/mata/dwlsmv_psibar.mata            ///
             src/mata/dwls_gamma.mata               ///
             src/mata/dwlsmv_emptycell1.mata        ///
             src/mata/dwlsmv_emptycell_repair.mata  ///
             src/mata/dwlsmv_emptycell_pairs.mata   ///
             src/mata/dwls_sample_stats.mata        ///
             src/mata/dwlsmv_nfac.mata              ///
             src/mata/dwlsmv_parvals.mata           ///
             src/mata/dwlsmv_matrices.mata          ///
             src/mata/dwlsmv_rystar.mata            ///
             src/mata/dwlsmv_corrvec.mata           ///
             src/mata/dwls_model.mata               ///
             src/mata/dwlsmv_delta_cols.mata        ///
             src/mata/dwls_delta.mata               ///
             src/mata/dwlsmv_canon_lines.mata       ///
             src/mata/dwlsmv_canon_find.mata        ///
             src/mata/dwlsmv_canon_status.mata      ///
             src/mata/dwlsmv_numrank.mata           ///
             src/mata/dwls_identify.mata            ///
             src/mata/dwlsmv_fdwls.mata             ///
             src/mata/dwlsmv_minimize.mata          ///
             src/mata/dwlsmv_signcanon.mata         ///
             src/mata/dwls_fit.mata                 ///
             src/mata/dwlsmv_asym.mata              ///
             src/mata/dwlsmv_acov.mata              ///
             src/mata/dwlsmv_traces.mata            ///
             src/mata/dwlsmv_scaleshift.mata        ///
             src/mata/dwls_inference.mata           ///
             src/mata/dwlsmv_geomin.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_geomin: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/expected/lavaan_efa_geomin.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: GC = st_sdata(., "case")
mata: GQ = st_sdata(., "quantity")
mata: GL = st_sdata(., "label")
mata: GV = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)
import delimited using "fixtures/data/case_F.csv", clear varnames(1) delimiters(",")
mata: YF = st_data(., .)
import delimited using "fixtures/data/case_G.csv", clear varnames(1) delimiters(",")
mata: YG = st_data(., .)

mata:

real scalar tg_get(string colvector C, string colvector Q, string colvector L,
                   real colvector V, string scalar c, string scalar q,
                   string scalar l)
{
    real colvector k

    k = selectindex((C :== c) :& (Q :== q) :& (L :== l))
    if (rows(k) != 1 | cols(k) != 1) {
        return(.)
    }
    return(V[k])
}

real scalar tg_eps(real scalar m)
{
    return(m == 2 ? 1e-4 : (m == 3 ? 1e-3 : 1e-2))
}

// echelon canonical spec: factor k is anchored at zero on items 1..k-1
string scalar tg_spec(real scalar p, real scalar m)
{
    string scalar s
    real scalar j, k

    s = "factor"
    for (k = 1; k <= m; k++) {
        s = s + " f" + strofreal(k)
    }
    s = s + char(10) + "ind"
    for (j = 1; j <= p; j++) {
        s = s + " y" + strofreal(j)
    }
    s = s + char(10)
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            s = s + "lambda f" + strofreal(k) + " y" + strofreal(j) +
                (j < k ? " fixed 0" : " free") + char(10)
        }
    }
    for (k = 1; k <= m; k++) {
        s = s + "metric unitvariance f" + strofreal(k) + char(10)
    }
    return(s)
}

real scalar tg_ratio(real colvector x, real colvector y, real scalar a, real scalar r)
{
    if (missing(x) | missing(y) | rows(x) != rows(y)) {
        return(.)
    }
    return(max(abs(x - y) :/ (a :+ r :* abs(y))))
}

// ||Gp||_F at (A, Lambda), with T' recovered by least squares from L T' = A
real scalar tg_gpnorm(real matrix A, real matrix L, real scalar eps)
{
    real matrix Tt, T, G

    Tt = luinv(L' * L) * L' * A
    T  = Tt'
    G  = dwlsmv_geomin_grad(A, T, eps)
    return(sqrt(sum((G - T * diag(colsum(T :* G))):^2)))
}

real scalar tg_case(string scalar cs, real matrix Y, real scalar m,
                    string colvector C, string colvector Q, string colvector L,
                    real colvector V)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    real matrix A, Psi, B, Lr, Phi
    real rowvector info
    real colvector K, lam, llam, phv, lphv, avg
    real scalar p, j, k, r, nf, eps, gp

    nf  = 0
    p   = cols(Y)
    K   = J(p, 1, .)
    for (j = 1; j <= p; j++) {
        K[j] = max(Y[., j]) + 1
    }
    st  = dwls_sample_stats(Y, K)
    par = dwls_identify(tg_spec(p, m), st.map)
    fit = dwls_fit(st, par)
    if (fit.converged != 1) {
        printf("{err}  case %s: dwls_fit() did not converge\n", cs)
        return(1)
    }
    dwls_inference(fit)
    printf("{txt}  case %s (p = %g, m = %g)\n", cs, p, m)

    r = dwlsmv_gapratio(fit.T3, tg_get(C, Q, L, V, cs, "test", "T3"), 1e-6, 1e-5)
    printf("{txt}    T3 %12.8f vs %12.8f  ratio %8.2e;  df %g vs %g\n", fit.T3,
           tg_get(C, Q, L, V, cs, "test", "T3"), r, fit.T3_df,
           tg_get(C, Q, L, V, cs, "test", "df"))
    if (!(r <= 1)) nf++
    if (fit.T3_df != tg_get(C, Q, L, V, cs, "test", "df")) nf++

    dwlsmv_matrices(dwlsmv_parvals(fit.theta, par), par, p, m, A, Psi, B)
    eps = tg_eps(m)
    if (eps != tg_get(C, Q, L, V, cs, "rot", "eps")) {
        printf("{err}    eps %g is not the fixture's %g\n", eps,
               tg_get(C, Q, L, V, cs, "rot", "eps"))
        nf++
    }
    Lr = dwlsmv_geomin(A, eps, 30, Phi, info)

    // 1  against lavaan, columns compared directly
    lam  = vec(Lr)
    llam = J(p * m, 1, .)
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            llam[(k - 1) * p + j] = tg_get(C, Q, L, V, cs, "lambda",
                                           "f" + strofreal(k) + "=~y" + strofreal(j))
        }
    }
    r = tg_ratio(lam, llam, 1e-5, 1e-4)
    printf("{txt}    rotated Lambda ratio %8.2e (max abs %8.2e)\n", r, max(abs(lam - llam)))
    if (!(r <= 1)) nf++

    phv  = J(0, 1, .)
    lphv = J(0, 1, .)
    for (k = 2; k <= m; k++) {
        for (j = 1; j < k; j++) {
            phv  = phv \ Phi[k, j]
            lphv = lphv \ tg_get(C, Q, L, V, cs, "phi",
                                 "f" + strofreal(k) + "~~f" + strofreal(j))
        }
    }
    r = tg_ratio(phv, lphv, 1e-5, 1e-4)
    printf("{txt}    Phi ratio %8.2e (max abs %8.2e)\n", r, max(abs(phv - lphv)))
    if (!(r <= 1)) nf++

    r = dwlsmv_gapratio(info[5], tg_get(C, Q, L, V, cs, "rot", "Q"), 1e-5, 1e-4)
    printf("{txt}    Q %14.10f vs %14.10f  ratio %8.2e\n", info[5],
           tg_get(C, Q, L, V, cs, "rot", "Q"), r)
    if (!(r <= 1)) nf++

    // 2  sigma unchanged
    r = max(abs(Lr * Phi * Lr' - A * A'))
    printf("{txt}    Lambda Phi Lambda' - A A'  max abs %8.2e\n", r)
    if (!(r <= 1e-12)) nf++

    // 3  stationarity at the returned solution
    gp = tg_gpnorm(A, Lr, eps)
    printf("{txt}    ||Gp|| recomputed %8.2e, reported %8.2e\n", gp, info[6])
    if (!(gp < 1e-5) | !(info[6] < 1e-5)) nf++

    // 4  sign and order, recomputed
    if (any(colsum(Lr) :<= 0)) {
        printf("{err}    a column sum is not positive\n")
        nf++
    }
    avg = J(m, 1, .)
    for (k = 1; k <= m; k++) {
        avg[k] = mean(selectindex(abs(Lr[., k]) :>= 0.8 * max(abs(Lr[., k]))))
    }
    for (k = 2; k <= m; k++) {
        if (avg[k] < avg[k - 1]) {
            printf("{err}    index order violated at column %g\n", k)
            nf++
        }
    }

    // 5  info well-formed
    printf("{txt}    info: starts %g, converged %g, distinct %g, at best %g, gap %g\n",
           info[1], info[2], info[3], info[4], info[8])
    if (info[1] != 30 | info[2] < 1 | info[2] > 30 | info[4] < 1 |
        info[4] > info[2] | info[3] < 1 | info[7] != eps) nf++
    if ((info[3] == 1) != missing(info[8])) nf++
    return(nf)
}

real scalar tg_fd(real scalar seed)
{
    real matrix A, T, D, Ap, Am, Tp, Tm
    real scalar m, i, j, eps, h, r, nf

    nf  = 0
    h   = 1e-6
    eps = 1e-3
    rseed(seed)
    for (m = 2; m <= 3; m++) {
        A = (0.2 :+ 0.7 :* runiform(9, m)) :* (2 :* (runiform(9, m) :> 0.3) :- 1)
        T = dwlsmv_geomin_unit(rnormal(m, m, 0, 1))
        D = J(9, m, .)
        for (i = 1; i <= 9; i++) {
            for (j = 1; j <= m; j++) {
                Ap = A
                Am = A
                Ap[i, j] = Ap[i, j] + h
                Am[i, j] = Am[i, j] - h
                D[i, j] = (dwlsmv_geomin_q(Ap, eps) - dwlsmv_geomin_q(Am, eps)) / (2 * h)
            }
        }
        r = max(abs(D - dwlsmv_geomin_dq(A, eps))) / max(abs(dwlsmv_geomin_dq(A, eps)))
        printf("{txt}  m = %g: dQ/dLambda vs central difference, rel %8.2e\n", m, r)
        if (!(r <= 1e-6)) nf++
        D = J(m, m, .)
        for (i = 1; i <= m; i++) {
            for (j = 1; j <= m; j++) {
                Tp = T
                Tm = T
                Tp[i, j] = Tp[i, j] + h
                Tm[i, j] = Tm[i, j] - h
                D[i, j] = (dwlsmv_geomin_q(dwlsmv_geomin_rot(A, Tp), eps) -
                           dwlsmv_geomin_q(dwlsmv_geomin_rot(A, Tm), eps)) / (2 * h)
            }
        }
        r = max(abs(D - dwlsmv_geomin_grad(A, T, eps))) /
            max(abs(dwlsmv_geomin_grad(A, T, eps)))
        printf("{txt}  m = %g: dQ/dT vs central difference, rel %8.2e\n", m, r)
        if (!(r <= 1e-6)) nf++
    }
    return(nf)
}

real scalar tg_det(real matrix Y)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    real matrix A, Psi, B, L1, L2, L3, P1, P2, P3
    real rowvector i1, i2, i3
    real scalar nf, tol

    nf  = 0
    st  = dwls_sample_stats(Y, J(9, 1, 4))
    par = dwls_identify(tg_spec(9, 3), st.map)
    fit = dwls_fit(st, par)
    dwlsmv_matrices(dwlsmv_parvals(fit.theta, par), par, 9, 3, A, Psi, B)
    L1 = dwlsmv_geomin(A, 1e-3, 30, P1, i1)
    L2 = dwlsmv_geomin(A, 1e-3, 30, P2, i2)
    L3 = dwlsmv_geomin(A, 1e-3, 40, P3, i3)
    if (L1 != L2 | P1 != P2 | i1 != i2) {
        printf("{err}  two identical calls differ\n")
        nf++
    }
    tol = 1e-8 + 1e-6 * i1[5]
    printf("{txt}  determinism: identical %s; nesting: Q(40) - Q(30) = %9.2e (tol %8.2e)\n",
           (nf ? "NO" : "yes"), i3[5] - i1[5], tol)
    if (!(i3[5] <= i1[5] + tol) | !(abs(i3[5] - i1[5]) <= tol)) nf++
    return(nf)
}

real scalar tg_m1()
{
    real matrix A, L, Phi
    real rowvector info

    A = (-0.7 \ -0.6 \ -0.5 \ 0.1)
    L = dwlsmv_geomin(A, 1e-2, 30, Phi, info)
    printf("{txt}  m = 1: Phi %g, column sum %g\n", Phi, sum(L))
    return((Phi != 1) + (sum(L) <= 0) + (L != -A))
}

// 10  the convergence rule and the reported ||Gp|| (mutation round,
// 2026-10-29: every fixture start stalls near 1e-8, so neither the 1e-5
// threshold nor the post-loop recompute was exercised). For maxit 1..40:
// a call either aborts (no start converged) or returns a solution whose
// reported ||Gp|| is below 1e-5 and equals the norm recomputed from L.
real scalar tg_maxit(real matrix A)
{
    real matrix L, Phi
    real rowvector info
    real scalar k, nf, nret, gp, rc

    nf   = 0
    nret = 0
    for (k = 1; k <= 40; k++) {
        rc = _stata("mata: TL = dwlsmv_geomin(TA, 1e-4, 30, TP, TI, " + strofreal(k) + ")", 1)
        if (rc) {
            // the only allowed abort is no converged start
            if (rc != 3498) {
                printf("{err}  maxit %g: aborted at rc %g, not 3498\n", k, rc)
                nf++
            }
            continue
        }
        nret++
        L    = st_matrix_tmp_L()
        info = st_matrix_tmp_I()
        gp   = tg_gpnorm(A, L, 1e-4)
        if (!(info[6] < 1e-5) | !(abs(info[6] - gp) <= 1e-6 * gp + 1e-12)) {
            printf("{err}  maxit %g: reported ||Gp|| %8.2e, recomputed %8.2e\n", k, info[6], gp)
            nf++
        }
    }
    printf("{txt}  maxit sweep 1..40: %g calls returned, %g violated the rule\n", nret, nf)
    return(nf + (nret == 0) + (nret == 40))
}

real matrix st_matrix_tmp_L()
{
    external real matrix TL
    return(TL)
}

real rowvector st_matrix_tmp_I()
{
    external real rowvector TI
    return(TI)
}

// 11  the GM4 grid as pinned (@sec-geomin GM4), recomputed independently:
// angles pi*frac(k(q+1)c), pairs in combn order, product left to right
real scalar tg_grid()
{
    real matrix S, R, P
    real scalar k, q, a, c, i, j, nf, d

    nf = 0
    c  = (sqrt(5) - 1) / 2
    P  = (1, 2 \ 1, 3 \ 2, 3)
    for (k = 1; k <= 7; k++) {
        S = I(3)
        for (q = 1; q <= 3; q++) {
            a = pi() * (k * (q + 1) * c - floor(k * (q + 1) * c))
            i = P[q, 1]
            j = P[q, 2]
            R = I(3)
            R[i, i] = cos(a)
            R[j, j] = cos(a)
            R[i, j] = -sin(a)
            R[j, i] = sin(a)
            S = S * R
        }
        d = max(abs(S - dwlsmv_geomin_start(3, k)))
        if (!(d <= 1e-14)) {
            printf("{err}  GM4 start %g differs from the pinned grid by %8.2e\n", k, d)
            nf++
        }
    }
    printf("{txt}  GM4 starts 1..7 (m = 3) match the pinned definition: %s\n", (nf ? "NO" : "yes"))
    return(nf)
}

// 12  the "index" order where its 0.8 cutoff decides it: column 1 has its
// largest loading on y6 and moderate ones on y1, y2; column 2 is large on
// y3, y4. At 0.8, column 2 comes first; a 0.2 cutoff would put column 1
// first. Checked on the returned L, with the cutoff's sensitivity printed.
real scalar tg_order()
{
    real matrix A, L, Phi
    real rowvector info
    real colvector a8, a2
    real scalar k, nf

    A = (0.35, 0 \ 0.35, 0 \ 0, 0.75 \ 0, 0.75 \ 0, 0.05 \ 0.9, 0)
    L = dwlsmv_geomin(A, 1e-4, 30, Phi, info)
    a8 = J(2, 1, .)
    a2 = J(2, 1, .)
    for (k = 1; k <= 2; k++) {
        a8[k] = mean(selectindex(abs(L[., k]) :>= 0.8 * max(abs(L[., k]))))
        a2[k] = mean(selectindex(abs(L[., k]) :>= 0.2 * max(abs(L[., k]))))
    }
    printf("{txt}  index order: mean index at 0.8 (%g, %g), at 0.2 (%g, %g)\n",
           a8[1], a8[2], a2[1], a2[2])
    nf = (a8[1] > a8[2])
    // the case must discriminate, or the check proves nothing
    if (!(a2[1] > a2[2])) {
        printf("{err}  the synthetic case does not separate the 0.8 and 0.2 cutoffs\n")
        nf++
    }
    return(nf)
}

end

local nfail = 0
mata: st_local("n1", strofreal(tg_case("E", YE, 2, GC, GQ, GL, GV)))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tg_case("F", YF, 3, GC, GQ, GL, GV)))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tg_case("G", YG, 2, GC, GQ, GL, GV)))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tg_fd(20261029)))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tg_det(YF)))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tg_m1()))
local nfail = `nfail' + `n1'

* 9  error paths
mata: TA = (0.7, 0.1 \ 0.6, 0.2 \ 0.1, 0.7 \ 0.2, 0.6 \ 0.5, 0.5)
* the output arguments must exist, or every call fails at rc 3499 instead
mata: TP = .
mata: TI = .
* negative control: the same call with valid arguments succeeds
capture mata: TL = dwlsmv_geomin(TA, 1e-4, 30, TP, TI)
if _rc {
    display as error "  error-path control: the valid call failed, rc " _rc
    local ++nfail
}
* each error path must abort at rc 3498 AND for its own reason: its message
* is checked too (mutation round, 2026-10-29: eps = 0 still reached rc 3498
* by another route when the eps check was weakened)
capture program drop tg_run
program tg_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tglog)
    capture noisily `0'
    local rc = _rc
    quietly log close tglog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end
local call1 "dwlsmv_geomin(TA, 0, 30, TP, TI)"
local why1  "eps must be positive"
local call2 "dwlsmv_geomin(TA, ., 30, TP, TI)"
local why2  "eps must be positive"
local call3 "dwlsmv_geomin(TA, 1e-4, 0, TP, TI)"
local why3  "nstart must be a positive integer"
local call4 "dwlsmv_geomin(TA, 1e-4, 2.5, TP, TI)"
local why4  "nstart must be a positive integer"
local call5 "dwlsmv_geomin(TA :* (1 \ 1 \ . \ 1 \ 1), 1e-4, 30, TP, TI)"
local why5  "no missing value"
local call6 "dwlsmv_geomin(TA[1, .], 1e-4, 30, TP, TI)"
local why6  "exceed"
local call7 "dwlsmv_geomin(TA, 1e-4, 30, TP, TI, 0)"
local why7  "none of the"
forvalues k = 1/7 {
    tg_run "`why`k''" mata: TL = `call`k''
    display as text "  error path `k': rc " r(rc) ", reason found " r(found)
    if r(rc) != 3498 | r(found) != 1 {
        display as error "  error path `k' (`call`k''): expected rc 3498 with '`why`k''', got rc " r(rc)
        local ++nfail
    }
}

* 10-12  convergence rule and reported ||Gp||, the pinned grid, the index cutoff
mata: st_local("n1", strofreal(tg_maxit(TA)))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tg_grid()))
local nfail = `nfail' + `n1'
mata: st_local("n1", strofreal(tg_order()))
local nfail = `nfail' + `n1'

if `nfail' {
    display as error _n "test_geomin: `nfail' failure(s)"
    exit 9
}
display as result "test_geomin: rotated Lambda, Phi and Q match lavaan geomin on Cases E, F, G; sigma unchanged; stationary; convention, gradients, determinism, nesting, m = 1 and error paths hold"
exit 0
