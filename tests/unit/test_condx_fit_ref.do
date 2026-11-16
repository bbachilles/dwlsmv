*! Executed 2026-11-14 by the main session (Stata 18 MP); mutants M1 to M5 killed (M4 by an rc 3498 abort); M6 not run.
*! REVISED 2026-11-14 on the owner's decisions (the p-value leaves this comparison; T_raw is on the T3 tier). The
*! revision is NOT YET EXECUTED; M1 and M4 below were killed in their earlier form (p included) and must be re-run.
*!
*! Conditional mode (one 0/1 covariate, one group): dwls_fit() and dwls_inference() against a tier-1
*! (base-R) reference on a dataset that DOES NOT FIT EXACTLY.
*!
*! WHY. tests/unit/test_condx_fit.do covers exact-fit populations only (F = 0, zero residual), where the
*! slope-row weights in W_D, T_raw = (N - 1) F, the slope rows of the inference and the effect of the
*! zero starts for gamma and kappa cannot move a number. Here the fitted model omits a true
*! cross-loading, so F > 0 and each of them does (E-105, "Conditional fit: gamma
*! and kappa start at zero": the required validation item).
*!
*! NOT LINK 2. This file is not link 2 of E-098 and does not read
*! condx_A, condx_D, condx_link3 or condx_ref_A_*, condx_ref_D_*. Its data and reference are the new
*! fixtures of fixtures/generate/99o_condx_ref_nonexact.R (MANIFEST_condx_ref_U.json):
*!     fixtures/data/condx_U.csv                    y1 .. y6, x   (N = 600; K = 3,2,4,3,2,3)
*!     fixtures/expected/condx_ref_U_s.csv          label, value
*!     fixtures/expected/condx_ref_U_gamma.csv      row, col, value  (lower triangle with diagonal)
*!     fixtures/expected/condx_ref_U_theta.csv      label, est, se
*!     fixtures/expected/condx_ref_U_scalars.csv    quantity, value  (N, npar, df, F, T_raw, a, b, T3, p_T3, ...)
*! read with import delimited, stringcols(_all), delimiters(",") and strtoreal().
*!
*! MODEL (two factors, gamma on both, one kappa, no cross-loading; 21 free parameters, D = 32, df = 11):
*!     F1 =~ y1 y2 y3, F2 =~ y4 y5 y6, unitvariance F1 F2, psi F1 F2 free, covariate x,
*!     gamma F1 x, gamma F2 x, kappa y3 x. theta order (spec 4.8): tau (11), lambda (6), gamma (2),
*!     kappa (1), psi_12 (1).
*!
*! TOLERANCES. All are the preregistered link-2 table, set before any result and applied here to a
*! different dataset on the owner's decision; none is chosen from a measurement. |a - b| <=
*! e_abs + e_rel |b|:
*!     s (thresholds, slopes, correlations)   1e-8, 1e-8   E-098
*!                                                         (slopes take the threshold tier, :134, [V5])
*!     Gamma-hat, every entry                 1e-8, 1e-8   prereg.md:126
*!     theta-hat                              1e-5, 1e-4   prereg.md:127; used the same way at
*!                                                         tests/unit/test_condx_fit.do:22 and :265
*!     SE = sqrt(diag(acov) / (N - 1))        1e-8, 1e-5   prereg.md:128; test_condx_fit.do:274
*!     T3                                     1e-6, 1e-5   prereg.md:129; test_condx_fit.do:281
*!     a                                      none, 1e-5   prereg.md:130; test_condx_fit.do:283
*!     b                                      1e-6, none   prereg.md:131; test_condx_fit.do:285
*!     df                                     exact        prereg.md:132; test_condx_fit.do:253
*! NO EXISTING TIER, DECIDED BY THE OWNER 2026-11-14:
*!     T_raw and F   the preregistration has no row for T_raw (test_condx_fit.do:24 records that the
*!                   earlier draft borrowed the T3 tier and that this was the agent's, not the
*!                   prereg's). Accepted by the owner 2026-11-14: the T3 tier on T_raw, and on
*!                   (N - 1) F-hat, because T3 = a T_raw + b is on the same scale as T_raw and has the same
*!                   reference interpretation. A T_raw of N F instead of (N - 1) F differs by
*!                   T_raw / (N - 1), a relative 1.7e-3, which is 170 times the relative term.
*!     p-value       NOT COMPARED here (owner, 2026-11-14): no tier exists for it, and a reference
*!                   comparison of p cannot tell a wrong T3 from a wrong transformation. T3 is compared
*!                   above; the transformation chi2tail(df, T3) is tested separately, against closed-form
*!                   chi-square survival probabilities and at the command, in tests/unit/test_pvalue_T3.do.
*!
*! CHECKS (one dataset):
*!   0  structure: npar 21, the table's theta order (ptype, and row and col for lambda, gamma, kappa), the
*!      first 11 labels equal the reference's; every free gamma and kappa row starts at exactly 0
*!      (the zero-start decision) -- the precondition of check 5.
*!   1  Stage 1 s and Gamma against the reference, by label (precondition: a failure here is Stage 1's).
*!   2  dwls_fit() from the table's starts (ULS, then DWLS): converged, D_stat == T3_df == df exactly.
*!   3  theta-hat, every entry; (N - 1) fmin and T_raw (T_raw and F on the T3 tier, accepted 2026-11-14).
*!   4  dwls_inference(): SE, a, b, T3 (no p: tests/unit/test_pvalue_T3.do).
*!   5  zero starts: dwlsmv_fit_from() from tau = s_tau, lambda = 0.5, gamma = kappa = psi_12 = 0
*!      converges and lands on the same theta-hat, fmin and T_raw (same tiers).
*!   6  F is clearly positive on this dataset (the dataset's design floor; not a tolerance on any
*!      comparison): fmin >= 1e-3, so the checks above are not vacuous.
*!
*! WHAT THIS CERTIFIES AND DOES NOT. It certifies that production dwls_fit() and dwls_inference() agree
*! with a base-R Hessian-form reference on one misfitting two-factor conditional dataset with gamma, kappa
*! and psi_12 (an exact-fit population cannot). It does not certify several covariates, groups,
*! admissibility, fit indices, DIFFTEST or |rho| near 1, and one dataset is one dataset. It does not certify the
*! p-value or its transformation from T3 (tests/unit/test_pvalue_T3.do), and a pass here says nothing about it.
*!
*! MUTANTS (edit, run this file, revert; each must turn it red; line numbers as of the drafting commit):
*!   M1  src/mata/dwls_fit.mata:292  f.T_raw      = (N - G) * fm   ->   f.T_raw      = N * fm
*!       (kills: T_raw and T3, check 3 and 4; fmin is unchanged so the (N - 1) fmin line still passes)
*!   M2  src/mata/dwls_fit.mata:129  return(1 :/ g)   ->
*!       return(1 :/ (g :* (st.map.stat_type :!= 3) :+ (st.map.stat_type :== 3)))
*!       (weight 1 on the slope rows: theta-hat, SE, T_raw move on a misfitting dataset: checks 3, 4, 5)
*!   M3  src/mata/dwls_inference.mata:74  return(sts[1].Gamma)   ->
*!       return(sts[1].Gamma :* (I(rows(sts[1].Gamma)) :+ ((sts[1].map.stat_type :!= 3) * (sts[1].map.stat_type :!= 3)')))
*!       (slope-row covariances dropped from the Gamma the inference sees: SE, a, b, T3: check 4;
*!       the fit itself is unchanged, so checks 2, 3 and 5 pass)
*!   M4  src/mata/dwls_fit.mata:293  f.D_stat     = rows(s) - par.npar   ->   rows(s) - par.npar - sts[1].map.p
*!       (the slopes dropped from the degrees of freedom: df, a, b, T3: checks 2 and 4)
*!   M5  src/mata/dwls_inference.mata:134  ab[1] * f.T_raw + ab[2]   ->   ab[1] * f.T_raw
*!       (the shift dropped from T3: check 4, T3 moves by b)
*!   M6  the zero start: in dwls_identify's default for gamma or kappa rows (the place that writes
*!       par.value for ptype 8 and 9) start at 0.5 instead of 0: check 0 fails at once. If the fit still
*!       converges to the same point, the start is not what a reference check can see; the check-0
*!       assertion is then the only guard, which is the intended division of labour.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

foreach f in tests/dwlsmv_assert.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_condx_fit_ref: `f' failed to load"
        exit 9
    }
}

* Load every source, types first (as test_condx_fit.do does).
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
        display as error "test_condx_fit_ref: `f' failed to load"
        exit 9
    }
}

capture mata: (void) dwlsmv_cx1_fit_check(dwlsmv_stage1(0), "load check")
if _rc {
    display as error "test_condx_fit_ref: dwlsmv_cx1_fit_check() is not defined"
    exit 9
}

* ---- fixtures (stringcols(_all) and strtoreal(); delimiters(",") explicit) -------------
foreach f in data/condx_U.csv expected/condx_ref_U_s.csv expected/condx_ref_U_gamma.csv ///
             expected/condx_ref_U_theta.csv expected/condx_ref_U_scalars.csv {
    capture confirm file "fixtures/`f'"
    if _rc {
        display as error "test_condx_fit_ref: fixtures/`f' is missing (run fixtures/generate/99o_condx_ref_nonexact.R)"
        exit 9
    }
}

import delimited using "fixtures/expected/condx_ref_U_s.csv", clear varnames(1) delimiters(",") stringcols(_all)
mata: SL = st_sdata(., "label")
mata: SV = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/condx_ref_U_gamma.csv", clear varnames(1) delimiters(",") stringcols(_all)
mata: GR = st_sdata(., "row")
mata: GC = st_sdata(., "col")
mata: GV = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/expected/condx_ref_U_theta.csv", clear varnames(1) delimiters(",") stringcols(_all)
mata: TL = st_sdata(., "label")
mata: TE = strtoreal(st_sdata(., "est"))
mata: TS = strtoreal(st_sdata(., "se"))

import delimited using "fixtures/expected/condx_ref_U_scalars.csv", clear varnames(1) delimiters(",") stringcols(_all)
mata: QN = st_sdata(., "quantity")
mata: QV = strtoreal(st_sdata(., "value"))

import delimited using "fixtures/data/condx_U.csv", clear varnames(1) delimiters(",") stringcols(_all)
mata: Y = strtoreal(st_sdata(., ("y1", "y2", "y3", "y4", "y5", "y6")))
mata: X = strtoreal(st_sdata(., "x"))

mata:

// ---- helpers -------------------------------------------------------------------

real scalar cr_fail(string scalar what)
{
    printf("{err}test_condx_fit_ref: %s\n", what)
    return(1)
}

real scalar cr_chk(real scalar ok, string scalar what)
{
    if (!ok) {
        return(cr_fail(what))
    }
    return(0)
}

// A tier comparison reported as the gap ratio (<= 1 passes); a missing ratio fails.
real scalar cr_tier(real matrix a, real matrix b, real scalar eabs, real scalar erel, string scalar what)
{
    real scalar r

    r = dwlsmv_gapratio(a, b, eabs, erel)
    printf("%s: gap ratio %9.3e\n", what, r)
    if (missing(r) | r > 1) {
        return(cr_fail(what + " is outside its tier, ratio " + strofreal(r, "%9.3e")))
    }
    return(0)
}

real scalar cr_look(string colvector QN, real colvector QV, string scalar nm)
{
    real scalar i

    for (i = 1; i <= rows(QN); i++) {
        if (QN[i] == nm) {
            return(QV[i])
        }
    }
    return(.)
}

string scalar cr_join(string rowvector a)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(a); i++) {
        s = s + a[i] + char(10)
    }
    return(s)
}

// The fitted model of the header.
string scalar cr_spec()
{
    return(cr_join(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
                    "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                    "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                    "metric unitvariance F1", "metric unitvariance F2",
                    "covariate x", "gamma F1 x free", "gamma F2 x free",
                    "kappa y3 x free", "psi F1 F2 free")))
}

// The reference's slope label "y<j>~x" is the map's "y<j>~x1".
string scalar cr_lab(string scalar L)
{
    if (strlen(L) >= 2) {
        if (substr(L, strlen(L) - 1, 2) == "~x") {
            return(L + "1")
        }
    }
    return(L)
}

// Position of a reference label in the map; 0 if absent, -1 if ambiguous.
real scalar cr_pos(struct dwlsmv_statmap scalar m, string scalar L)
{
    string scalar M
    real scalar k, hit

    M   = cr_lab(L)
    hit = 0
    for (k = 1; k <= m.D; k++) {
        if (m.label[k] == M) {
            if (hit) {
                return(-1)
            }
            hit = k
        }
    }
    return(hit)
}

// ---- check 0: the table's structure and the zero starts -------------------------------
real scalar cr_struct(struct dwlsmv_parmap scalar par, string colvector TL)
{
    real colvector pt, rw, cl, wpt, wrw, wcl
    real scalar nf, i, nz

    nf = 0
    nf = nf + cr_chk(par.npar == 21, "check 0: npar is " + strofreal(par.npar) + ", not 21")
    if (par.npar != 21) {
        return(nf)
    }
    pt = J(21, 1, .)
    rw = J(21, 1, .)
    cl = J(21, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            pt[par.fidx[i]] = par.ptype[i]
            rw[par.fidx[i]] = par.row[i]
            cl[par.fidx[i]] = par.col[i]
        }
    }
    wpt = (J(11, 1, 4) \ J(6, 1, 1) \ 8 \ 8 \ 9 \ 2)
    nf = nf + cr_chk(all(pt :== wpt), "check 0: theta's ptype order is not tau, lambda, gamma, kappa, psi")
    // lambda: row = indicator, col = factor; gamma: row = factor; kappa: row = indicator, col = covariate
    wrw = (1 \ 2 \ 3 \ 4 \ 5 \ 6 \ 1 \ 2 \ 3)
    wcl = (1 \ 1 \ 1 \ 2 \ 2 \ 2 \ 1 \ 1 \ 1)
    nf = nf + cr_chk(all(rw[12::20] :== wrw) & all(cl[12::20] :== wcl),
                     "check 0: the lambda, gamma, kappa rows are not (y1..y6 on F1, F1, F1, F2, F2, F2), gamma F1, F2, kappa y3")
    for (i = 1; i <= 11; i++) {
        nf = nf + cr_chk(par.name[i] == TL[i], "check 0: table label " + par.name[i] + " vs reference " + TL[i])
    }
    // the zero-start decision (2026-11-13): every free gamma and kappa row starts at exactly 0
    nz = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & (par.ptype[i] == 8 | par.ptype[i] == 9)) {
            nz = nz + 1
            nf = nf + cr_chk(par.value[i] == 0,
                             "check 0: a free gamma or kappa row starts at " + strofreal(par.value[i]) + ", not 0")
        }
    }
    nf = nf + cr_chk(nz == 3, "check 0: " + strofreal(nz) + " free gamma and kappa rows, not 3")
    return(nf)
}

// ---- check 1: Stage 1 s and Gamma, by label -----------------------------------------
real scalar cr_stage1(struct dwlsmv_stage1 scalar st, string colvector SL, real colvector SV,
                      string colvector GR, string colvector GC, real colvector GV)
{
    struct dwlsmv_statmap scalar m
    real colvector R
    real matrix Gr
    real scalar nf, k, D, pa, pb, nexp

    nf = 0
    m  = st.map
    D  = m.D
    nf = nf + cr_chk(D == 32 & rows(SL) == 32 & !missing(SV), "check 1: the s fixture is not D = 32 complete rows")
    R = J(D, 1, .)
    for (k = 1; k <= rows(SL); k++) {
        pa = cr_pos(m, SL[k])
        if (pa <= 0) {
            return(nf + cr_fail("check 1: s label not in the map: " + SL[k]))
        }
        R[pa] = SV[k]
    }
    if (missing(R)) {
        return(nf + cr_fail("check 1: the s fixture does not cover every position"))
    }
    nf = nf + cr_tier(st.s, R, 1e-8, 1e-8, "check 1: s")

    nexp = D * (D + 1) / 2
    Gr   = J(D, D, .)
    if (rows(GR) != nexp | missing(GV)) {
        return(nf + cr_fail("check 1: the Gamma fixture is not D(D+1)/2 complete rows"))
    }
    for (k = 1; k <= nexp; k++) {
        pa = cr_pos(m, GR[k])
        pb = cr_pos(m, GC[k])
        if (pa <= 0 | pb <= 0) {
            return(nf + cr_fail("check 1: Gamma label not in the map: " + GR[k] + ", " + GC[k]))
        }
        Gr[pa, pb] = GV[k]
        Gr[pb, pa] = GV[k]
    }
    if (missing(Gr)) {
        return(nf + cr_fail("check 1: the Gamma fixture does not fill the D x D matrix"))
    }
    nf = nf + cr_tier(st.Gamma, Gr, 1e-8, 1e-8, "check 1: Gamma")
    return(nf)
}

// ---- checks 2 to 6 ----------------------------------------------------------------------
real scalar cr_one(real matrix Y, real colvector X, string colvector SL, real colvector SV,
                   string colvector GR, string colvector GC, real colvector GV,
                   string colvector TL, real colvector TE, real colvector TS,
                   string colvector QN, real colvector QV)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f, f2
    real colvector K, se, start
    real scalar nf, N, df, T3r

    nf = 0
    K  = (3 \ 2 \ 4 \ 3 \ 2 \ 3)
    st = dwls_sample_stats_cx(Y, K, X)
    par = dwls_identify(cr_spec(), st.map)
    sv = dwlsmv_stage1(1, 1)
    sv[1] = st

    // 0 and 1
    nf = nf + cr_struct(par, TL)
    if (par.npar != 21) {
        return(nf)
    }
    nf = nf + cr_stage1(st, SL, SV, GR, GC, GV)

    // 2: the fit from the table's starts, then inference
    N  = st.N
    df = cr_look(QN, QV, "df")
    nf = nf + cr_chk(N == cr_look(QN, QV, "N"), "check 2: N is " + strofreal(N) + ", the reference has " + strofreal(cr_look(QN, QV, "N")))
    f = dwls_fit(sv, par)
    dwls_inference(f)
    nf = nf + cr_chk(f.converged == 1, "check 2: the fit from the table's starts did not converge")
    nf = nf + cr_chk(f.D_stat == df & f.T3_df == df & st.map.D - par.npar == df,
                     "check 2: D_stat is " + strofreal(f.D_stat) + ", T3_df " + strofreal(f.T3_df) + ", wanted " + strofreal(df))
    if (rows(f.theta) != 21 | rows(f.acov) != 21) {
        return(nf + cr_fail("check 2: theta or acov has the wrong size; later checks skipped"))
    }

    // 3: theta-hat, (N - 1) fmin, T_raw
    nf = nf + cr_tier(f.theta, TE, 1e-5, 1e-4, "check 3: theta-hat")
    nf = nf + cr_tier((N - 1) * f.fmin, cr_look(QN, QV, "T_raw"), 1e-6, 1e-5, "check 3: (N - 1) fmin vs reference T_raw")
    nf = nf + cr_tier(f.T_raw, cr_look(QN, QV, "T_raw"), 1e-6, 1e-5, "check 3: T_raw")

    // 4: inference
    se = sqrt(diagonal(f.acov) :/ (N - 1))
    nf = nf + cr_tier(se, TS, 1e-8, 1e-5, "check 4: SE")
    nf = nf + cr_tier(f.shift_a, cr_look(QN, QV, "a"), 0, 1e-5, "check 4: a")
    nf = nf + cr_tier(f.shift_b, cr_look(QN, QV, "b"), 1e-6, 0, "check 4: b")
    T3r = cr_look(QN, QV, "T3")
    nf = nf + cr_tier(f.T3, T3r, 1e-6, 1e-5, "check 4: T3")
    // the p-value is not compared here (owner, 2026-11-14); see tests/unit/test_pvalue_T3.do

    // 5: an explicit zero start: tau = s_tau, lambda 0.5, gamma = kappa = psi_12 = 0
    start = (st.s[1::11] \ J(6, 1, 0.5) \ 0 \ 0 \ 0 \ 0)
    f2 = dwlsmv_fit_from(sv, par, start)
    nf = nf + cr_chk(f2.converged == 1, "check 5: the fit from the explicit zero start did not converge")
    nf = nf + cr_tier(f2.theta, TE, 1e-5, 1e-4, "check 5: theta-hat from the zero start")
    nf = nf + cr_tier(f2.T_raw, cr_look(QN, QV, "T_raw"), 1e-6, 1e-5, "check 5: T_raw from the zero start")

    // 6: the dataset's design floor
    nf = nf + cr_chk(f.fmin >= 1e-3, "check 6: fmin is " + strofreal(f.fmin) + ", the data fit (nearly) exactly")

    printf("condx_U: converged %g, iterations %g, fmin %12.8f, T_raw %12.8f, T3 %12.8f, a %10.6f, b %10.6f\n",
           f.converged, f.iterations, f.fmin, f.T_raw, f.T3, f.shift_a, f.shift_b)
    return(nf)
}

end

local nfail 0
capture noisily mata: st_local("nf", strofreal(cr_one(Y, X, SL, SV, GR, GC, GV, TL, TE, TS, QN, QV)))
if _rc {
    display as error "test_condx_fit_ref: cr_one() aborted, rc " _rc
    local nfail = `nfail' + 1
}
else {
    display as text "  cr_one: `nf' failure(s)"
    local nfail = `nfail' + `nf'
}

if `nfail' {
    display as error _n "test_condx_fit_ref: `nfail' failure(s)"
    exit 9
}

display as result "test_condx_fit_ref: dwls_fit and dwls_inference agree with the tier-1 reference on a misfitting two-factor conditional dataset (theta-hat, T_raw, SE, a, b, T3, df; zero starts)"
exit 0
