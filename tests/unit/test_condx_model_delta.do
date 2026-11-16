*! Conditional mode (one 0/1 covariate, one group): dwls_model() and dwls_delta().
*!
*! NOT YET EXECUTED. Drafted 2026-11-13 by a read-only agent; the main session
*! runs it, fixes what breaks, and only then removes this banner.
*!
*! Needs slopes_changes.md applied (dwlsmv_cx1_check.mata, dwls_model, dwls_delta,
*! dwlsmv_corrvec, dwlsmv_delta_cols). It does NOT depend on the dwls_identify
*! draft: dwls_identify() does not parse covariates yet, so the parameter table
*! is built by hand (cx_par()) with the identify draft's ptype conventions:
*! gamma ptype 8 (row = factor, col = covariate 1), kappa ptype 9 (row =
*! indicator, col = covariate 1), tau ptype 4 (row = variable, col = threshold),
*! lambda ptype 1 (row = indicator, col = factor), beta ptype 3 (row = receiving
*! factor, col = sending factor), psi ptype 2. Table order is theta order: tau,
*! lambda, beta, gamma, kappa, psi; fixed rows carry fidx 0.
*!
*!   (a) sigma(theta) of a 3-item one-factor model with gamma and one kappa, and
*!       of a 6-item two-factor model, against closed forms computed here with
*!       explicit indices (not the position helpers): tau, then slopes
*!       b = Lambda (I-B)^-1 Gamma + kappa, then rho = offdiag(Lambda Phi Lambda').
*!   (b) Delta against central finite differences of dwls_model() at a generic
*!       theta. Tolerance 1e-8 absolute on every entry. The analytic slope rows
*!       and the existing finite-difference correlation rows (step 1e-5*max(1,|v|),
*!       error floor about 1e-10, per dwlsmv_delta_cols.mata) are compared with a
*!       central difference at step 1e-6*max(1,|v|) (truncation ~1e-12, rounding
*!       ~1e-10). The two floors add to about 3e-10, so 1e-8 leaves a factor of
*!       30; it was fixed before running, and is not to be widened to make a
*!       failure pass. The slope rows are ALSO compared with explicit closed-form
*!       derivatives at 1e-13 (analytic against analytic, rounding only).
*!   (c) the threshold block of Delta is the identity, exactly, and the cross
*!       blocks are exactly zero.
*!   (d) two factors with beta F2 <- F1 and x predicting F1: the slopes of F2's
*!       items pass through (I-B)^-1 (b_j = lambda_j * beta * gamma).
*!   (e) CANNOT be tested in one tree: "a model without covariates gives the
*!       same sigma and Delta as before the change, bit for bit". What is
*!       asserted instead: (i) at q = 0 the correlation offset D - P equals T
*!       exactly and the unconditional model matches its closed form; (ii) the
*!       same model at q = 1 with gamma = 0 and kappa = 0 has threshold and
*!       correlation rows of sigma, and every threshold and correlation entry of
*!       Delta's structural columns, EXACTLY (difference 0) equal to the q = 0
*!       model's, and slope rows exactly 0; (iii) gamma and kappa columns are
*!       exactly 0 in every correlation row. The real bit-for-bit check is the
*!       main session's: run the existing test_model.do, test_delta.do and the
*!       suite before and after the edits and diff their logs.
*!   (f) refusals, each asserted at rc 3498 exactly: two groups with a covariate
*!       (model, delta); q = 2; gamma rows on a q = 0 map; a nu row at q = 1;
*!       gamma with covariate 2; gamma on a factor that does not exist; kappa on
*!       an indicator that does not exist.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

* Load every source, types first (as test_s2_guard.do does).
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_condx_model_delta: `f' failed to load"
        exit 9
    }
}

capture mata: (void) dwlsmv_cx1_check(dwlsmv_parmap(), dwlsmv_statmap(0), "load check")
if _rc {
    display as error "test_condx_model_delta: dwlsmv_cx1_check() is not defined (src/mata/dwlsmv_cx1_check.mata missing?)"
    exit 9
}

mata:

// ---- helpers -----------------------------------------------------------------

real scalar cx_fail(string scalar what)
{
    printf("{err}test_condx_model_delta: %s\n", what)
    return(1)
}

// A hand-built parameter table, all rows in theta order; fidx numbers the free rows.
struct dwlsmv_parmap scalar cx_par(real colvector pt, real colvector rw, real colvector cl,
                                   real colvector fr, real colvector vl)
{
    struct dwlsmv_parmap scalar par
    real scalar i, n, k

    n = rows(pt)
    par.nall  = n
    par.npar  = sum(fr)
    par.name  = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        par.name[i] = "row" + strofreal(i)
    }
    par.ptype = pt
    par.row   = rw
    par.col   = cl
    par.free  = fr
    par.value = vl
    par.fidx  = J(n, 1, 0)
    k = 0
    for (i = 1; i <= n; i++) {
        if (fr[i]) {
            k = k + 1
            par.fidx[i] = k
        }
    }
    par.grp = J(n, 1, 1)
    return(par)
}

// Model 1: K = (3,3,2), one factor, unit disturbance variance (psi fixed 1).
// theta = (tau: 5 | lambda: 3 | gamma | kappa_3), 10 free. g, k set gamma and kappa.
struct dwlsmv_parmap scalar cx_par1(real scalar g, real scalar k)
{
    return(cx_par((4 \ 4 \ 4 \ 4 \ 4 \ 1 \ 1 \ 1 \ 8 \ 9 \ 2),
                  (1 \ 1 \ 2 \ 2 \ 3 \ 1 \ 2 \ 3 \ 1 \ 3 \ 1),
                  (1 \ 2 \ 1 \ 2 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1),
                  (1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 0),
                  (-0.5 \ 0.6 \ -0.4 \ 0.5 \ 0.2 \ 0.8 \ 0.7 \ 0.6 \ g \ k \ 1)))
}

// Model 0: the same without x: theta = (tau: 5 | lambda: 3), 8 free.
struct dwlsmv_parmap scalar cx_par0()
{
    return(cx_par((4 \ 4 \ 4 \ 4 \ 4 \ 1 \ 1 \ 1 \ 2),
                  (1 \ 1 \ 2 \ 2 \ 3 \ 1 \ 2 \ 3 \ 1),
                  (1 \ 2 \ 1 \ 2 \ 1 \ 1 \ 1 \ 1 \ 1),
                  (1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 0),
                  (-0.5 \ 0.6 \ -0.4 \ 0.5 \ 0.2 \ 0.8 \ 0.7 \ 0.6 \ 1)))
}

// Model 2: K = (3,2,2,3,2,2), F1 -> y1..y3, F2 -> y4..y6, beta F2 <- F1,
// gamma F1 <- x, kappa on y2 and y5, psi_11 fixed 1, psi_22 free.
// theta = (tau: 8 | lambda: 6 | beta | gamma | kappa_2, kappa_5 | psi_22), 19 free.
struct dwlsmv_parmap scalar cx_par2()
{
    return(cx_par((4\4\4\4\4\4\4\4 \ 1\1\1\1\1\1 \ 3 \ 8 \ 9\9 \ 2 \ 2),
                  (1\1\2\3\4\4\5\6 \ 1\2\3\4\5\6 \ 2 \ 1 \ 2\5 \ 1 \ 2),
                  (1\2\1\1\1\2\1\1 \ 1\1\1\2\2\2 \ 1 \ 1 \ 1\1 \ 1 \ 2),
                  (J(18, 1, 1) \ 0 \ 1),
                  (-0.5\0.6\-0.4\0.3\-0.2\0.7\0.1\-0.3 \ 0.8\0.7\0.6\0.75\0.65\0.55 \ 0.4 \ 0.5 \ -0.3\0.25 \ 1 \ 0.8)))
}

// table order is theta order for the free rows, so theta is the free values
real colvector cx_theta(struct dwlsmv_parmap scalar par)
{
    return(select(par.value, par.free))
}

struct dwlsmv_statmap vector cx_maps(struct dwlsmv_statmap scalar m)
{
    struct dwlsmv_statmap vector mv

    mv = dwlsmv_statmap(1, 1)
    mv[1] = m
    return(mv)
}

// central differences of dwls_model(), step 1e-6 * max(1, |theta_k|)
real matrix cx_fd(real colvector th, struct dwlsmv_parmap scalar par,
                  struct dwlsmv_statmap scalar m)
{
    struct dwlsmv_statmap vector mv
    real matrix Dn
    real colvector tp, tm
    real scalar k, h

    mv = cx_maps(m)
    Dn = J(m.D, rows(th), .)
    for (k = 1; k <= rows(th); k++) {
        h     = 1e-6 * max((1, abs(th[k])))
        tp    = th
        tm    = th
        tp[k] = th[k] + h
        tm[k] = th[k] - h
        Dn[., k] = (dwls_model(tp, par, mv) - dwls_model(tm, par, mv)) / (2 * h)
    }
    return(Dn)
}

// ---- closed forms with explicit indices -----------------------------------------

// Model 1: rows 1..5 tau, 6..8 slopes, 9..11 rho(1,2), rho(1,3), rho(2,3)
real colvector cx_expect1(real colvector th)
{
    real colvector lam, b
    real scalar gam, kap

    lam = th[6::8]
    gam = th[9]
    kap = th[10]
    b   = lam * gam
    b[3] = b[3] + kap
    return(th[1::5] \ b \ (lam[1] * lam[2] \ lam[1] * lam[3] \ lam[2] * lam[3]))
}

// Model 2: rows 1..8 tau, 9..14 slopes, 15..29 rho in column-major
// lower-triangle order (pairs by first index, then second)
real colvector cx_expect2(real colvector th)
{
    real matrix Lam, B, Psi, A, Phi, R
    real colvector Gam, Kap, b, rho
    real scalar i, j, n

    Lam = J(6, 2, 0)
    Lam[1, 1] = th[9];  Lam[2, 1] = th[10]; Lam[3, 1] = th[11]
    Lam[4, 2] = th[12]; Lam[5, 2] = th[13]; Lam[6, 2] = th[14]
    B   = (0, 0 \ th[15], 0)
    Gam = (th[16] \ 0)
    Kap = (0 \ th[17] \ 0 \ 0 \ th[18] \ 0)
    Psi = (1, 0 \ 0, th[19])
    A   = luinv(I(2) - B)
    b   = Lam * A * Gam + Kap
    Phi = A * Psi * A'
    R   = Lam * Phi * Lam'
    rho = J(15, 1, .)
    n = 0
    for (i = 1; i <= 5; i++) {
        for (j = i + 1; j <= 6; j++) {
            n = n + 1
            rho[n] = R[i, j]
        }
    }
    return(th[1::8] \ b \ rho)
}

// Explicit slope-block derivatives (rows = items, columns = all 19 theta of model 2)
real matrix cx_dslope2(real colvector th)
{
    real matrix S
    real colvector lamv
    real scalar bet, gam, j

    lamv = th[9::14]
    bet  = th[15]
    gam  = th[16]
    S    = J(6, 19, 0)
    // b_j = lambda_j * gamma (j = 1..3), lambda_j * beta * gamma (j = 4..6), + kappa
    for (j = 1; j <= 3; j++) {
        S[j, 8 + j] = gam
        S[j, 16]    = lamv[j]
    }
    for (j = 4; j <= 6; j++) {
        S[j, 8 + j] = bet * gam
        S[j, 15]    = lamv[j] * gam
        S[j, 16]    = lamv[j] * bet
    }
    S[2, 17] = 1
    S[5, 18] = 1
    return(S)
}

// ---- (a) to (d) ----------------------------------------------------------------

real scalar cx_check_a()
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar m
    real colvector th, sg, ex
    real scalar nf, tol

    nf  = 0
    tol = 1e-12     // the exact-identity row of spec 7.2: one closed form against another

    // model 1, the table's own values and a generic theta
    par = cx_par1(0.5, -0.3)
    m   = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
    th  = cx_theta(par)
    sg  = dwls_model(th, par, cx_maps(m))
    if (rows(sg) != 11) {
        nf = nf + cx_fail("(a) model 1: sigma has " + strofreal(rows(sg)) + " rows, wanted 11")
    }
    else if (max(abs(sg - cx_expect1(th))) > tol) {
        nf = nf + cx_fail("(a) model 1: sigma differs from the closed form by " +
                          strofreal(max(abs(sg - cx_expect1(th)))))
    }
    th = (-0.45 \ 0.55 \ -0.35 \ 0.45 \ 0.15 \ 0.75 \ 0.62 \ 0.55 \ 0.35 \ -0.2)
    sg = dwls_model(th, par, cx_maps(m))
    if (max(abs(sg - cx_expect1(th))) > tol) {
        nf = nf + cx_fail("(a) model 1 at a generic theta: sigma differs from the closed form by " +
                          strofreal(max(abs(sg - cx_expect1(th)))))
    }

    // model 2 (d): slopes through (I-B)^-1
    par = cx_par2()
    m   = dwlsmv_statmap_build((3 \ 2 \ 2 \ 3 \ 2 \ 2), 1)
    th  = cx_theta(par)
    sg  = dwls_model(th, par, cx_maps(m))
    if (rows(sg) != 29) {
        nf = nf + cx_fail("(d) model 2: sigma has " + strofreal(rows(sg)) + " rows, wanted 29")
    }
    else {
        ex = cx_expect2(th)
        if (max(abs(sg - ex)) > tol) {
            nf = nf + cx_fail("(d) model 2: sigma differs from the closed form by " +
                              strofreal(max(abs(sg - ex))))
        }
        // the point of the design: an F2 item's slope is lambda * beta * gamma (+ kappa)
        if (abs(sg[9 + 3] - th[12] * th[15] * th[16]) > tol) {
            nf = nf + cx_fail("(d) model 2: slope of y4 is not lambda_4 * beta * gamma")
        }
    }
    return(nf)
}

real scalar cx_check_bcd()
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar m
    real matrix Dm, Dn, S, Z
    real colvector th
    real scalar nf, tolfd, tolcf, pass, T, D, np

    nf    = 0
    tolfd = 1e-8
    tolcf = 1e-13

    for (pass = 1; pass <= 2; pass++) {
        if (pass == 1) {
            par = cx_par1(0.5, -0.3)
            m   = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
            th  = (-0.45 \ 0.55 \ -0.35 \ 0.45 \ 0.15 \ 0.75 \ 0.62 \ 0.55 \ 0.35 \ -0.2)
        }
        else {
            par = cx_par2()
            m   = dwlsmv_statmap_build((3 \ 2 \ 2 \ 3 \ 2 \ 2), 1)
            th  = (-0.45\0.55\-0.35\0.25\-0.15\0.65\0.05\-0.25 \ 0.78\0.66\0.58\0.72\0.61\0.52 \ 0.35 \ 0.45 \ -0.25\0.2 \ 0.85)
        }
        T  = m.T
        D  = m.D
        np = par.npar

        Dm = dwls_delta(th, par, cx_maps(m))
        if (rows(Dm) != D | cols(Dm) != np) {
            nf = nf + cx_fail("(b) model " + strofreal(pass) + ": Delta is " + strofreal(rows(Dm)) +
                              " x " + strofreal(cols(Dm)) + ", wanted " + strofreal(D) + " x " + strofreal(np))
            continue
        }

        // (b) against central differences of dwls_model
        Dn = cx_fd(th, par, m)
        if (max(abs(Dm - Dn)) > tolfd) {
            nf = nf + cx_fail("(b) model " + strofreal(pass) + ": Delta differs from finite differences by " +
                              strofreal(max(abs(Dm - Dn))))
        }

        // (c) the threshold block is the identity; the cross blocks are exactly zero
        if (max(abs(Dm[1::T, 1::T] - I(T))) != 0) {
            nf = nf + cx_fail("(c) model " + strofreal(pass) + ": the threshold block of Delta is not the identity")
        }
        if (max(abs(Dm[1::T, (T + 1)::np])) != 0) {
            nf = nf + cx_fail("(c) model " + strofreal(pass) + ": a structural column has a nonzero threshold row")
        }
        if (max(abs(Dm[(T + 1)::D, 1::T])) != 0) {
            nf = nf + cx_fail("(c) model " + strofreal(pass) + ": a threshold column has a nonzero slope or correlation row")
        }

        // slope block (rows T+1 .. T+p) against explicit closed forms
        if (pass == 1) {
            S = J(3, 10, 0)
            S[1, 6] = th[9]; S[2, 7] = th[9]; S[3, 8] = th[9]     // d b_j / d lambda_j = gamma
            S[1, 9] = th[6]; S[2, 9] = th[7]; S[3, 9] = th[8]     // d b_j / d gamma = lambda_j
            S[3, 10] = 1                                          // d b_3 / d kappa_3
        }
        else {
            S = cx_dslope2(th)
        }
        if (max(abs(Dm[(T + 1)::(T + m.p), .] - S)) > tolcf) {
            nf = nf + cx_fail("(b) model " + strofreal(pass) + ": slope rows differ from the closed-form derivatives by " +
                              strofreal(max(abs(Dm[(T + 1)::(T + m.p), .] - S))))
        }

        // (e iii) gamma and kappa do not enter any correlation: exactly zero there.
        // The gamma and kappa columns are the last ones before psi.
        if (pass == 1) {
            Z = Dm[(T + m.p + 1)::D, 9::10]
        }
        else {
            Z = Dm[(T + m.p + 1)::D, 16::18]
        }
        if (max(abs(Z)) != 0) {
            nf = nf + cx_fail("(e) model " + strofreal(pass) + ": a gamma or kappa column has a nonzero correlation row")
        }
        // psi does not enter any slope (pass 2 has a free psi_22, the last column)
        if (pass == 2) {
            if (max(abs(Dm[(T + 1)::(T + m.p), 19])) != 0) {
                nf = nf + cx_fail("(e) model 2: the psi column has a nonzero slope row")
            }
        }
    }
    return(nf)
}

// ---- (e) structural invariants in place of a bit-for-bit comparison ---------------

real scalar cx_check_e()
{
    struct dwlsmv_parmap  scalar p0, p1
    struct dwlsmv_statmap scalar m0, m1
    real matrix D0, D1
    real colvector th0, th1, s0, s1, ex
    real scalar nf, i, j, k, tol

    nf  = 0
    tol = 1e-12

    m0  = dwlsmv_statmap_build((3 \ 3 \ 2))
    m1  = dwlsmv_statmap_build((3 \ 3 \ 2), 1)

    // (i) q = 0: the correlation offset is T exactly, and sigma is the closed form
    if (m0.D - m0.P != m0.T) {
        nf = nf + cx_fail("(e) q = 0: D - P is not T")
    }
    if (m1.D - m1.P != m1.T + 3) {
        nf = nf + cx_fail("(e) q = 1: D - P is not T + p")
    }
    p0  = cx_par0()
    th0 = cx_theta(p0)
    s0  = dwls_model(th0, p0, cx_maps(m0))
    ex  = th0[1::5] \ (th0[6] * th0[7] \ th0[6] * th0[8] \ th0[7] * th0[8])
    if (rows(s0) != 8 | max(abs(s0 - ex)) > tol) {
        nf = nf + cx_fail("(e) q = 0: sigma differs from the closed form")
    }
    D0 = dwls_delta(th0, p0, cx_maps(m0))
    if (max(abs(D0[1::5, 1::5] - I(5))) != 0 | max(abs(D0[1::5, 6::8])) != 0 |
        max(abs(D0[6::8, 1::5])) != 0) {
        nf = nf + cx_fail("(e) q = 0: the threshold block of Delta is not the identity")
    }

    // (ii) q = 1 with gamma = 0, kappa = 0 against q = 0, exactly
    p1  = cx_par1(0, 0)
    th1 = cx_theta(p1)
    s1  = dwls_model(th1, p1, cx_maps(m1))
    D1  = dwls_delta(th1, p1, cx_maps(m1))
    if (max(abs(s1[1::5] - s0[1::5])) != 0) {
        nf = nf + cx_fail("(e) q = 1, gamma = kappa = 0: threshold rows differ from the q = 0 model")
    }
    if (max(abs(s1[9::11] - s0[6::8])) != 0) {
        nf = nf + cx_fail("(e) q = 1, gamma = kappa = 0: correlation rows differ from the q = 0 model")
    }
    if (max(abs(s1[6::8])) != 0) {
        nf = nf + cx_fail("(e) q = 1, gamma = kappa = 0: the slope rows are not exactly 0")
    }
    if (max(abs(D1[1::5, 1::8] - D0[1::5, 1::8])) != 0) {
        nf = nf + cx_fail("(e) q = 1: threshold rows of Delta differ from the q = 0 model")
    }
    if (max(abs(D1[9::11, 1::8] - D0[6::8, 1::8])) != 0) {
        nf = nf + cx_fail("(e) q = 1: correlation rows of Delta differ from the q = 0 model (form changed)")
    }
    return(nf)
}

// ---- (f) refusals: each scenario builds its inputs and calls one function ----------
// which: 1 = dwls_model, 2 = dwls_delta

void cx_refuse(real scalar which, real scalar scen)
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar m
    struct dwlsmv_statmap vector mv
    real colvector th, out, cc
    real matrix Dm

    par = cx_par1(0.5, -0.3)
    m   = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
    mv  = cx_maps(m)

    if (scen == 1) {                     // two groups with a covariate
        mv = dwlsmv_statmap(2, 1)
        mv[1] = m
        mv[2] = m
    }
    else if (scen == 2) {                // two covariates
        mv[1] = dwlsmv_statmap_build((3 \ 3 \ 2), 2)
    }
    else if (scen == 3) {                // gamma and kappa rows, no covariate in the map
        mv[1] = dwlsmv_statmap_build((3 \ 3 \ 2))
    }
    else if (scen == 4) {                // a nu row (ptype 5) at q = 1
        par = cx_par((4 \ 4 \ 4 \ 4 \ 4 \ 1 \ 1 \ 1 \ 5 \ 2),
                     (1 \ 1 \ 2 \ 2 \ 3 \ 1 \ 2 \ 3 \ 1 \ 1),
                     (1 \ 2 \ 1 \ 2 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1),
                     (1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 1 \ 0),
                     (-0.5 \ 0.6 \ -0.4 \ 0.5 \ 0.2 \ 0.8 \ 0.7 \ 0.6 \ 0 \ 1))
    }
    else if (scen == 5) {                // gamma on covariate 2
        cc = par.col; cc[9] = 2; par.col = cc
    }
    else if (scen == 6) {                // gamma on factor 2 of a one-factor model
        cc = par.row; cc[9] = 2; par.row = cc
    }
    else if (scen == 7) {                // kappa on indicator 4 of three
        cc = par.row; cc[10] = 4; par.row = cc
    }
    else if (scen == 8) {                // gamma and kappa rows, two groups, no covariate in either map
        mv = dwlsmv_statmap(2, 1)
        mv[1] = dwlsmv_statmap_build((3 \ 3 \ 2))
        mv[2] = mv[1]
    }

    th = cx_theta(par)
    if (which == 1) {
        out = dwls_model(th, par, mv)
    }
    else {
        Dm = dwls_delta(th, par, mv)
    }
}

end

* ---- run ---------------------------------------------------------------------

local nfail 0

mata: st_local("n", strofreal(cx_check_a()))
local nfail = `nfail' + `n'
mata: st_local("n", strofreal(cx_check_bcd()))
local nfail = `nfail' + `n'
mata: st_local("n", strofreal(cx_check_e()))
local nfail = `nfail' + `n'

forvalues which = 1/2 {
    forvalues scen = 1/8 {
        capture mata: cx_refuse(`which', `scen')
        if _rc != 3498 {
            display as error "test_condx_model_delta: (f) function `which', scenario `scen': rc is " _rc ", wanted 3498"
            local nfail = `nfail' + 1
        }
    }
}

if `nfail' > 0 {
    display as error "test_condx_model_delta: `nfail' failure(s)"
    exit 9
}
display as text "test_condx_model_delta: pass"
exit 0
