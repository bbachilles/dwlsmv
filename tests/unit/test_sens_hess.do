*! test_sens_hess.do -- Result MS6's exact routines: dwlsmv_sens_hess and the
*! exact variants dwlsmv_sens_ddelta_x, dwlsmv_sens_jmat_x, dwlsmv_sens_x,
*! dwlsmv_sens_t3_x (plan E-031, step 3)
*!
*! Executed under StataMP 16.1 on 2026-10-19: tests/unit/test_sens_hess.do green; suite 62/0/0 at four processors.
*!
*! Design E (fixtures/data/mg_E.csv, G = 2). References: tests/unit/
*! sens_hess_ref.mata, written by tools/export_sens_hess_ref.R from the .qmd's
*! own mg2_* and mg2h_* chunks at 12 FIXED points (g, s and c at each of the
*! four levels), so they do not depend on Stata's version or processor count.
*! Fits for the x3-t2 checks: the metric and scalar fits of
*! tests/unit/sens_t3_setup.mata.
*!
*!   p.k  precondition per point: Stata's table has the reference's npar, and
*!        dwls_model at the reference theta reproduces the reference sigma
*!        (1e-12 relative), so theta is in the same order on both sides
*!   x1.k dwlsmv_sens_hess against the reference H: 1e-12 relative to
*!        max(1, max|H|), and the zero pattern identical
*!   x2   H symmetric in (i, l) and zero on (l, l), bitwise, every point
*!   x3   dwlsmv_sens_ddelta_x(f) is dwlsmv_sens_hess at f.theta, bitwise; and
*!        dwlsmv_sens_ddelta (finite differences) agrees within 1e-6 of max|dD|
*!        (test_sens_t3's h3.b bound; measured noise 3e-8 in R, V3a)
*!   j1   dwlsmv_sens_jmat_x is P - R with R = sum_k (w o r)_k H_k, bitwise, and R
*!        symmetric, bitwise (metric and scalar fits)
*!   j2   dwlsmv_sens_jmat (finite differences) agrees with it within 1e-7 of
*!        max|J| (V3c's threshold)
*!   s1   dwlsmv_sens_x's A is cholsolve(dwlsmv_sens_jmat_x, (Delta :* w)'),
*!        bitwise (test_sens_t3's h2, for the exact J)
*!   s2   dwlsmv_sens_x against dwlsmv_sens, every column: 1e-4 relative
*!   t1   dwlsmv_sens_t3_x against dwlsmv_sens_t3 (metric_full, scalar_full):
*!        columns 1, 2, 4, 5, 6 and 8 bitwise equal (they read no second
*!        derivative); columns 3 and 7 within 1e-4 relative
*!   t2   theta-hat moved by ~1 ulp per element (theta :* (1 + 2 eps cos(k))):
*!        the exact path's beta_2 moves at most 1e-9 relative (added after
*!        mutants M8-M10 survived; threshold set before its first run)
*!   e1-e7 refusals, each rc 3498: dwlsmv_sens_hess with a beta row (c1), with
*!        two factors (c2), with two rows of one group sharing a theta element
*!        (c3); dwlsmv_sens_ddelta_x, dwlsmv_sens_jmat_x, dwlsmv_sens_x and
*!        dwlsmv_sens_t3_x on a fit whose table breaks c3
*!   e8   a hand-built one-factor table with a beta(1, 1) row (dwls_identify
*!        never builds one; dwlsmv_nfac still counts one factor): rc 3498 from
*!        dwlsmv_sens_hdom, with the message naming c1 and not c2 (added
*!        2026-10-19 for mutant M5)
*!   tj   the J inside dwlsmv_sens_t3_x: scaling the step scales Theta's S_C
*!        affinely, so beta_2's second difference over s = 0, L, 2L is
*!        2 a nT v'|J|v exactly; it must match dwlsmv_sens_jmat_x's J to 1e-13
*!        relative, and miss the finite-difference J by at least 1e-11 (so the
*!        check can tell them apart; added 2026-10-19 for mutant M9).
*!        CHOSEN AFTER THE FIRST RUN: that run used 1e-10 and a >= 1e-9 gap
*!        between the two J's (guesses), and failed only the gap: measured 0
*!        and 2.0e-16 against the exact J, 3.2e-10 and 4.3e-10 against the
*!        finite-difference J. The thresholds sit about 3 decades from each
*!        side
*!
*! Mutants (scratch copy %TEMP%\shmut, mutate.sh, each diff printed before its
*! run; mutants.log, mutants2.log, mutants3.log): 11 of 12 red here, M12 red in the suite.
*!   M1  {lambda, alpha} entry +om (x1, x3b, j2, s2, t1b)
*!   M2  {nu, omega} entry +1 (x1, x3b)
*!   M3  psi left out of the correlation monomial (x1, x3b, j2, s2, t1b)
*!   M4  only pairs b > a (x1, j1, j2, s2, t1b)
*!   M6  clause c3 removed (e3-e7)
*!   M7  R from r, not w o r, in dwlsmv_sens_jmat_x (j1, j2, s1)
*!   M8  dwlsmv_sens_t3_x's dD by finite differences (t2)
*!   M10 dwlsmv_sens_t3_x's Theta from dwlsmv_sens (t2)
*!   M11 dwlsmv_sens_x's R scaled by 1 + 1e-6 (s1)
*!   M12 dwlsmv_sens_t3's default path routed to exact = 1: survives here;
*!       test_sens_t3's references catch it (its beta_2 moves ~1e-8 relative)
*!   M5  clause c1 removed (e8; it survived until e8 was added, because every
*!       table dwls_identify builds with beta has a second factor and c2
*!       refused it, so e8 builds a one-factor beta(1, 1) table by hand)
*!   M9  dwlsmv_sens_t3_x's J by finite differences (tj; it survived until tj
*!       was added, because T_J is ~1e-5 of beta_2 and no check read J alone)
*! t2 was added after M8-M10 survived the first round.
*! t2's margin is thin on the metric fit (8.3e-10 against 1e-9): the exact path
*! keeps production Delta_c, whose ~1e-11 noise reaches beta_2 through P, ghat
*! and G_Delta; MS6 replaces only dDelta/dtheta and R.
*! Division of responsibility with test_sens_t3.do (clarified 2026-10-19): this
*! test owns Hessian correctness, against references from Result MS6's closed
*! form. test_sens_t3's exact references carry Stata's own dDelta/dtheta and J,
*! so it checks reproducibility of the assembly only; its survivors N3 (J by
*! finite differences inside dwlsmv_sens_t3_x) and N5 (a wrong {nu, omega}
*! entry) are caught here, by tj and x1. That split is intended.
*!
*! Tolerances, chosen before the first run: 1e-12 on p and x1 (R's arithmetic
*! against Mata's, entries O(1)); bitwise on x2, x3's first part, j1, s1 and
*! t1's equal columns; 1e-6, 1e-7 from the preregistered noise thresholds;
*! 1e-4 on s2 and t1's moving columns (no measurement exists; they only show
*! the exact path stays near the finite-difference one).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index
do tests/unit/sens_unpolished.mata
do tests/unit/sens_t3_setup.mata
do tests/unit/sens_hess_ref.mata

mata:

void sht_chk(string colvector FAIL, string scalar nm, real scalar ok, real scalar val)
{
    if (ok) {
        printf("{txt}  %-6s ok   (%9.2e)\n", nm, val)
    }
    else {
        printf("{err}  %-6s FAIL (%9.2e)\n", nm, val)
        FAIL = FAIL \ nm
    }
}

// The level's table and maps for design E, built as sens_t3_setup builds its fits.
void sht_par(string scalar lev, struct dwlsmv_parmap scalar par, struct dwlsmv_statmap vector maps)
{
    struct dwlsmv_stage1 vector sts
    real matrix K, Y
    real colvector gv
    real scalar G, g

    K    = st_matrix("e(K)")
    G    = rows(K)
    gv   = uniqrows(st_data(., "g"))
    sts  = dwlsmv_stage1(G, 1)
    maps = dwlsmv_statmap(G, 1)
    for (g = 1; g <= G; g++) {
        Y       = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
        sts[g]  = dwls_sample_stats(Y, K[g, .]')
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(st_global("e(canon_configural)"), K, lev), maps, sts)
}

real scalar sht_run()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    struct dwlsmv_fit scalar f, f2
    string colvector FAIL
    string rowvector id, cs
    real colvector th, sig, s, w, wr, step, ds, eps, Th0, vv
    real matrix nz, H, Hr, D, Jx, Jf, dD, A, Ax, sens, sx, dOm
    real rowvector dq, o, ox, o2, b0, b1, b2
    real scalar k, q, val, ok, l, i, N, c, nu, cx, nux, epsV, j, L, sd, a, nT, px, pf

    FAIL = J(0, 1, "")
    epsV = epsilon(1)

    // ---- p, x1, x2: the fixed points ---------------------------------------------
    for (k = 1; k <= sh_nref(); k++) {
        id = sh_ref(k, th, sig, nz, dq)
        sht_par(id[1], par, maps)
        q   = dq[2]
        val = (par.npar == q ? max(abs(dwls_model(th, par, maps) - sig)) / max((1, max(abs(sig)))) : .)
        sht_chk(FAIL, "p." + strofreal(k), val <= 1e-12, val)
        if (par.npar != q) continue

        Hr = J(dq[1], q * q, 0)
        for (i = 1; i <= rows(nz); i++) {
            Hr[nz[i, 1], nz[i, 2]] = nz[i, 3]
        }
        H   = dwlsmv_sens_hess(th, par, maps)
        val = max(abs(H - Hr)) / max((1, max(abs(Hr))))
        ok  = (val <= 1e-12) & ((H :!= 0) == (Hr :!= 0))
        sht_chk(FAIL, "x1." + strofreal(k), ok, val)

        ok = 1
        for (l = 1; l <= q; l++) {
            ok = ok & (H[., (l - 1) * q + l] == J(dq[1], 1, 0))
            for (i = 1; i <= q; i++) {
                ok = ok & (H[., (l - 1) * q + i] == H[., (i - 1) * q + l])
            }
        }
        sht_chk(FAIL, "x2." + strofreal(k), ok, 0)
    }

    // ---- x3, j1, j2, s1, s2, t1: the fits ----------------------------------------
    cs = ("metric", "scalar")
    for (k = 1; k <= 2; k++) {
        f    = *findexternal("ST3_F_" + cs[k])
        step = *findexternal("ST3_STEP_" + cs[k])
        dwlsmv_fit_stack(f.stage1, s, w, N, maps)
        th = f.theta
        q  = rows(th)
        H  = dwlsmv_sens_hess(th, f.par, maps)

        dD = dwlsmv_sens_ddelta_x(f)
        sht_chk(FAIL, "x3a." + strofreal(k), dD == H, max(abs(dD - H)))
        val = max(abs(dwlsmv_sens_ddelta(f) - dD)) / max(abs(dD))
        sht_chk(FAIL, "x3b." + strofreal(k), val <= 1e-6, val)

        Jx = dwlsmv_sens_jmat_x(f)
        D  = dwls_delta(th, f.par, maps)
        wr = w :* (s - dwls_model(th, f.par, maps))
        Hr = colshape(quadcross(H, wr)', q)'
        ok = (Jx == quadcross(D, w, D) - Hr) & (Hr == Hr')
        sht_chk(FAIL, "j1." + strofreal(k), ok, max(abs(Hr - Hr')))
        Jf  = dwlsmv_sens_jmat(f)
        val = max(abs(Jf - Jx)) / max(abs(Jx))
        sht_chk(FAIL, "j2." + strofreal(k), val <= 1e-7, val)

        sx = dwlsmv_sens_x(f, step, epsV, Ax, cx, nux)
        sht_chk(FAIL, "s1." + strofreal(k), Ax == cholsolve(Jx, (D :* w)'),
                max(abs(Ax - cholsolve(Jx, (D :* w)'))))
        sens = dwlsmv_sens(f, step, epsV, A, c, nu)
        val  = 0
        for (j = 1; j <= 4; j++) {
            val = max((val, max(abs(sx[., j] - sens[., j])) / max(abs(sens[., j]))))
        }
        val = max((val, max(abs(sx[., 6] - sens[., 6])) / max(abs(sens[., 6]))))
        sht_chk(FAIL, "s2." + strofreal(k), val <= 1e-4, val)

        st3_pert(f, ds, eps, dOm)
        o  = dwlsmv_sens_t3(f, step, epsV, 0.01, ds, eps, dOm)
        ox = dwlsmv_sens_t3_x(f, step, epsV, 0.01, ds, eps, dOm)
        ok = (ox[(1, 2, 4, 5, 6, 8)] == o[(1, 2, 4, 5, 6, 8)])
        sht_chk(FAIL, "t1a." + strofreal(k), ok, max(abs(ox[(1, 2, 4, 5, 6, 8)] - o[(1, 2, 4, 5, 6, 8)])))
        val = max(abs(ox[(3, 7)] - o[(3, 7)]) :/ abs(o[(3, 7)]))
        printf("{txt}         beta_2 exact %21.17g, finite differences %21.17g\n", ox[3], o[3])
        sht_chk(FAIL, "t1b." + strofreal(k), val <= 1e-4, val)

        // t2: theta-hat moved by about one ulp per element. The exact path's
        // beta_2 may move only by P's and T_g's own first-order change; the
        // finite-difference path redraws its noise (printed, not asserted).
        f2         = f
        f2.theta   = th :* (1 :+ 2 * epsilon(1) * cos(1::q))
        f2.Delta   = dwls_delta(f2.theta, f.par, maps)
        o2         = dwlsmv_sens_t3_x(f2, step, epsV, 0.01, ds, eps, dOm)
        val        = abs(o2[3] - ox[3]) / abs(ox[3])
        o2         = dwlsmv_sens_t3(f2, step, epsV, 0.01, ds, eps, dOm)
        printf("{txt}         beta_2 move at one ulp: exact %9.2e, finite differences %9.2e\n",
               val, abs(o2[3] - o[3]) / abs(o[3]))
        sht_chk(FAIL, "t2." + strofreal(k), val <= 1e-9, val)
    }


    // ---- tj: the J that dwlsmv_sens_t3_x uses -----------------------------------
    // The step enters Theta only through S_C, affinely, so scaling it by s gives
    // Theta(s) = Theta(0) + s v and beta_2 * (1 + |T3|) = lin(Theta) + a nT Theta'|J|Theta.
    // The second difference at s = 0, L, 2L is exactly 2 a nT v'|J|v: J itself.
    for (k = 1; k <= 2; k++) {
        f    = *findexternal("ST3_F_" + cs[k])
        step = *findexternal("ST3_STEP_" + cs[k])
        dwlsmv_fit_stack(f.stage1, s, w, N, maps)
        th   = f.theta
        sx   = dwlsmv_sens_x(f, step, epsV, Ax, cx, nux)
        L    = 1 / max((1 :+ abs(th)) :* sx[., 3])
        sx   = dwlsmv_sens_x(f, 0 * step, epsV, Ax, cx, nux)
        Th0  = (1 :+ abs(th)) :* (sx[., 3] + sx[., 6])
        sx   = dwlsmv_sens_x(f, L * step, epsV, Ax, cx, nux)
        vv   = (1 :+ abs(th)) :* (sx[., 3] + sx[., 6]) - Th0
        b0   = dwlsmv_sens_t3_x(f, 0 * step, epsV, 0.01)
        b1   = dwlsmv_sens_t3_x(f, L * step, epsV, 0.01)
        b2   = dwlsmv_sens_t3_x(f, 2 * L * step, epsV, 0.01)
        sd   = (1 + abs(b0[1])) * (b2[3] - 2 * b1[3] + b0[3])
        a    = f.shift_a
        nT   = N - length(maps)
        px   = 2 * a * nT * quadcross(vv, abs(dwlsmv_sens_jmat_x(f)) * vv)
        pf   = 2 * a * nT * quadcross(vv, abs(dwlsmv_sens_jmat(f)) * vv)
        val  = abs(sd - px) / abs(px)
        printf("{txt}         T_J second difference against exact J %9.2e, against finite-difference J %9.2e\n",
               val, abs(sd - pf) / abs(pf))
        sht_chk(FAIL, "tj." + strofreal(k), val <= 1e-13 & abs(sd - pf) / abs(pf) >= 1e-11, val)
    }

    if (rows(FAIL)) printf("{err}  failed: %s\n", invtokens(FAIL'))
    return(rows(FAIL))
}

// The refusals, each expected at rc 3498.
void sht_err(real scalar which)
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    struct dwlsmv_fit scalar f
    real colvector th, sig, step, ds, eps
    real matrix nz, X, A, dOm
    real rowvector dq, o
    real scalar i, i1, i2, c, nu

    if (which <= 3) {
        (void) sh_ref(4, th, sig, nz, dq)
        sht_par("thresholds", par, maps)
    }
    else {
        f    = *findexternal("ST3_F_metric")
        step = *findexternal("ST3_STEP_metric")
        par  = f.par
    }
    i1 = i2 = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1 & par.grp[i] == 1 & par.row[i] == 1) i1 = i
        if (par.ptype[i] == 1 & par.grp[i] == 1 & par.row[i] == 2) i2 = i
    }
    if (which == 8) {
        // one factor and a beta row, by hand (dwls_identify never builds it):
        // a group-1 loading turned into beta(1, 1), so dwlsmv_nfac is still 1
        // and only clause c1 can refuse it
        par.ptype[i2] = 3
        par.row[i2]   = 1
        par.col[i2]   = 1
        if (dwlsmv_nfac(par) != 1) {
            _error("test_sens_hess e8: the hand-built table is not one-factor")
        }
        dwlsmv_sens_hdom(par, 2, "test_sens_hess e8")
        return
    }
    if (which == 1) {
        par.ptype[i2] = 3
        par.row[i2]   = 2
        par.col[i2]   = 1
        X = dwlsmv_sens_hess(th, par, maps)
    }
    else if (which == 2) {
        par.col[i2] = 2
        X = dwlsmv_sens_hess(th, par, maps)
    }
    else {
        par.fidx[i2] = par.fidx[i1]
        if (which == 3) X = dwlsmv_sens_hess(th, par, maps)
        f.par = par
        if (which == 4) X = dwlsmv_sens_ddelta_x(f)
        if (which == 5) X = dwlsmv_sens_jmat_x(f)
        if (which == 6) X = dwlsmv_sens_x(f, step, epsilon(1), A, c, nu)
        if (which == 7) {
            st3_pert(f, ds, eps, dOm)
            o = dwlsmv_sens_t3_x(f, step, epsilon(1), 0.01, ds, eps, dOm)
        }
    }
}

end

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
mata: st3_setup()

local nf 0
mata: st_local("nf", strofreal(sht_run()))

forvalues k = 1/7 {
    capture noisily mata: sht_err(`k')
    if _rc != 3498 {
        display as error "    e`k' FAILED: rc `=_rc', expected 3498"
        local ++nf
    }
    else {
        display as text "    e`k' refused, rc 3498"
    }
}

* e8: clause c1 alone. rc 3498, and the message names c1 and not c2.
tempfile e8log
quietly log using "`e8log'", text name(e8) replace
capture noisily mata: sht_err(8)
local e8rc = _rc
quietly log close e8
local e8c1 0
local e8c2 0
tempname fh
file open `fh' using "`e8log'", read text
file read `fh' line
while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "(clause c1)") local e8c1 1
    if strpos(`"`macval(line)'"', "(clause c2)") local e8c2 1
    file read `fh' line
}
file close `fh'
if `e8rc' != 3498 | !`e8c1' | `e8c2' {
    display as error "    e8 FAILED: rc `e8rc' (wanted 3498), c1 named `e8c1', c2 named `e8c2'"
    local ++nf
}
else {
    display as text "    e8 refused, rc 3498, clause c1 named"
}

if `nf' {
    display as error "test_sens_hess: `nf' failure(s)"
    exit 9
}
display as result "test_sens_hess: the exact Hessian reproduces Result MS6's R at 12 fixed points, the exact variants are built from it and stay near their finite-difference versions, and every out-of-domain table is refused at rc 3498"
