*! test_fit_groups_betag1.do -- the pooled fit and inference under groups for a
*! model WITH a structural path FREE in group 1 and ABSENT in group 2 (g1only;
*! spec 4.8 O2 as widened 2026-10-11 after mgb-absent-prereg passed; derivation
*! @sec-mg-beta), against lavaan 0.7.2.
*!
*! FIRST RUN (2026-10-11): FAILED check 6 (e(V)) at thresholds, metric, scalar
*! (worst ratios SE 1.0 / 1.2 / 1.8, cov 0.60 / 0.72 / 1.3); configural and every
*! other check passed. Diagnosis after the failure (E-078, not a
*! test): at LAVAAN's theta our e(V) passes check 6 at every level (ratios <=
*! 8.6e-4); our theta is the better minimum (F_l - F_o = 3.5e-11, 5.9e-11,
*! 1.0e-10); max |theta_o - theta_l| about 1e-5. lavaan stops short of the
*! minimum; with the misspecified group (T about 800) that moves e(V) by about
*! 6e-6 relative against a 1e-5 tier. RESOLUTION (user, 2026-10-11, option (b);
*! chosen AFTER the failure): the reference is regenerated with lavaan's
*! optimizer at a stricter relative tolerance (99d_, control rel.tol 1e-14, set
*! before the rerun); this test, its tolerances and its rules are unchanged.
*! Declined: (a) extending check 4's attribution rule to e(V).
*! (b) ATTEMPTED 2026-10-11 and FAILED: at rel.tol 1e-14 lavaan reports that
*! its optimizer did not find a solution (configural, the first level), so the
*! generator's convergence assertion stopped it before writing. The generator
*! and the first-run reference are unchanged.
*! RESOLUTION (user, 2026-10-11, option (a) as a GENERAL POLICY): check 6 adopts
*! spec 7.2.1, the attribution rule for inference quantities (defined and
*! committed before this rerun). Elements of e(V): every SE sqrt(V_ii) and every
*! covariance V_ij (i < j, theta order). An element passes on its plain check-6
*! tier, or under 7.2.1 conditions 1-3; check 6 fails if any element passes
*! neither. Constants and method exactly as 7.2.1; none changed after a run.
*! RERUN 2026-10-11 (after 910c1db): green. Configural passed plainly; 5 e(V)
*! elements missed the plain tier (thresholds 1, metric 1, scalar 3), all 5
*! attributed by 7.2.1: (1) held at every level; (2) worst 8.6e-4 of tau
*! (required <= 1e-2); (3) worst gap 0.46 of its bound. Every other check passed
*! plainly. NEGATIVE CONTROL: ACOV scaled by (1 + 3e-5) in dwls_inference turned
*! check 6 red at all four levels, through condition 2 (1.5-2.3 of tau). Condition
*! 3 alone did not separate it (its bound carries the condition-2 residual), so
*! the rule's protection against a formula error rests on condition 2.
*!
*! A copy of test_fit_groups_betag2.do with the beta moved to group 1. As that file:
*! with these changes and nothing else: the texts write "beta F2 F1 free" in
*! group 1 only; the reference is fixtures/expected/lavaan_mg_betag1_*.csv from
*! fixtures/generate/99d_lavaan_mg_betag1.R (the same data, mg_beta.csv, whose
*! raw beta is .4 in both groups, so group 2 is deliberately misspecified; the
*! sample statistics still come from lavaan_mg_beta_sampstat.csv, model-free);
*! check 1's df is Result B4 with nb = 1 (17, 20, 21, 25); check b asserts exactly
*! one beta row, free, in group 1 (it is a parameter, so checks 2 and 6 compare it
*! with lavaan's). Checks s and 1-6, their tolerances and check 4's attribution
*! rule are the free-beta test's, applied unchanged and fixed before this test
*! first ran. The free-beta test's header follows.
*!
*! Design M of @sec-mg-beta (K = 4, 4, 2, 4, 2, 2; F1: y1-y3, F2: y4-y6, F2 on
*! F1), two groups of UNEQUAL size (1200, 2800), the four levels configural,
*! thresholds, metric, scalar, written as canonical 'groups' texts by hand
*! (Module 11 is not involved). Reference: see above (betag1 files, 99d_).
*! The checks, comparison rules and
*! tolerances are test_fit_groups.do's checks s and 1-6, unchanged; only the
*! sign equivalence is extended to two factors:
*!
*!   s   Stage 1 per group against lavaan's per-group statistics, by label:
*!       1e-8 abs + 1e-8 rel; a correlation outside that passes only under
*!       test_fit_groups.do's Newton-step rule (user, 2026-10-10, option (a)),
*!       applied unchanged.
*!   1   converged; D_stat = lavaan's df = Result B4 with nb = 1 (17, 20, 21, 25), exact.
*!   2   theta-hat against lavaan by (ptype, group, row, col), every free row,
*!       after sign equivalence: theta tier, 1e-5 abs + 1e-4 rel. Signs are
*!       matched PER FACTOR AND GROUP on that factor's first free loading in
*!       that group; a loading and an alpha take their factor's sign, a beta
*!       the product of its two factors' signs, psi (variances only) none.
*!   3   T_raw = (N - G) F-hat against (N - G) sum_g pi_g F_g from lavaan's own
*!       per-group F_g, with test_fit_groups.do's propagated bound; candidate C
*!       (lavaan's chi-square) must NOT match within that bound.
*!   4   a and b against lavaan's: a 1e-5 rel, b 1e-6 abs. An a or b outside its
*!       tier passes ONLY under the attribution rule below (user, 2026-10-10,
*!       option (a); defined here and committed BEFORE the rerun; chosen after
*!       the first run's failure, recorded under FIRST RUN).
*!
*!       ATTRIBUTION RULE (to lavaan stopping before our optimizer). With
*!       theta_o our theta-hat, theta_l lavaan's (sign-matched as in check 2,
*!       mapped to our theta order; a tie is one element), d = theta_o - theta_l,
*!       x one of a, b, x(.) our dwls_inference() value at a given theta (Delta
*!       rebuilt there; Gamma, pi and D_stat are the fit's), x_L lavaan's value:
*!       (i)   better minimum: F(theta_o) <= F(theta_l) (1 + 1e-12), F our DWLS
*!             objective. The 1e-12 is a rounding allowance on two evaluations
*!             of one function (chosen here, before the rerun).
*!       (ii)  exact formula: |x(theta_l) - x_L| <= 1e-8 (b), relative 1e-8 (a).
*!       (iii) predicted change: |x(theta_o) - x_L| <= 2 P + 1e-8 (b; for a,
*!             2 P + 1e-8 |x_L|), with P = sum_k |dx/dtheta_k| |d_k| (a triangle
*!             bound on the first-order change, so terms cannot cancel), each
*!             derivative a central difference at theta_o with step
*!             h_k = 1e-6 max(1, |theta_o,k|). The factor 2 is the Newton-step
*!             rule's mean-value factor (test_fit_groups.do check s).
*!       Edge cases, fixed with it: d = 0 gives P = 0, so (iii) reduces to
*!       1e-8; a missing or non-finite F, x or derivative FAILS the rule (it
*!       never passes by default); the rule is evaluated and printed for every
*!       level and both of a, b, but decides only where the plain tier failed;
*!       (i)-(iii) must all hold.
*!
*!       FIRST RUN (2026-10-10, before the rule existed): FAILED at the scalar
*!       level only, b gap 1.2e-6 > 1e-6 (a rel 1.9e-7, inside). Diagnosis
*!       after the failure (%TEMP% diag_tfb.do, not a test): at lavaan's theta
*!       our b matched lavaan's to 1.5e-9 and a to 3.7e-10; max |d| = 8.0e-6;
*!       F(theta_l) - F(theta_o) = 8.8e-11 > 0. Every other check passed.
*!   5   T3 = a T_raw + b against a_lav T_B,lav + b_lav: 1e-6 abs + 1e-5 rel.
*!   6   e(V) = ACOV / (N - G) against lavaan's vcov() by key (attribution: spec
*!       7.2.1, see RESOLUTION above): SEs 1e-8 abs +
*!       1e-5 rel; covariances within 1e-8 + 1e-5 sqrt(V_ii V_jj), sign-matched
*!       as in 2.
*!   b   beta: exactly one beta row, free, in group 1 (none in group 2); its
*!       estimate and SE are printed beside lavaan's (compared in checks 2 and 6).
*!
*! test_fit_groups.do's checks 7-10 (configural oracle, the no-pi control, the
*! O16 reflection) concern the pooling and the sign symmetry, not beta, and
*! stay there.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata tests/dwlsmv_assert.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_fit_groups_betag1: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/expected/lavaan_mg_betag1_fits.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: FITS = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_betag1_parameters.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: PARS = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_betag1_vcov.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: VCOV = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_beta_sampstat.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: SAMP = st_sdata(., .)
import delimited using "fixtures/data/mg_beta.csv", clear varnames(1) delimiters(",")
mata: YALL = st_data(., .)

mata:

string scalar tfb_join(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The canonical text of design M at a level, beta free in group 1 only (the
// level syntax of fixtures/generate/99d_lavaan_mg_betag1.R, @sec-mg-beta).
string scalar tfb_text(string scalar level)
{
    string rowvector L, base
    real rowvector C, KM
    real scalar j, k

    KM   = (4, 4, 2, 4, 2, 2)
    C    = KM :- 1
    base = ("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free")
    L = ("groups 2", "group 1", base, "beta F2 F1 free", "metric unitvariance F1", "metric unitvariance F2",
         "group 2", base)
    if (level == "configural" | level == "thresholds") {
        L = L, ("psi F1 F1 fixed 1", "psi F2 F2 fixed 1")
    }
    else {
        L = L, ("psi F1 F1 free", "psi F2 F2 free")
    }
    if (level == "thresholds" | level == "metric") {
        for (j = 1; j <= 6; j++) {
            L = L, ("nu y" + strofreal(j) + " free")
        }
    }
    for (j = 1; j <= 6; j++) {
        if (level == "metric" | level == "scalar" | (level == "thresholds" & KM[j] > 2)) {
            L = L, ("omega y" + strofreal(j) + " free")
        }
    }
    if (level == "scalar") {
        L = L, ("alpha F1 free", "alpha F2 free")
    }
    if (level != "configural") {
        for (j = 1; j <= 6; j++) {
            for (k = 1; k <= C[j]; k++) {
                L = L, ("tie tau 1 y" + strofreal(j) + " " + strofreal(k) +
                        " tau 2 y" + strofreal(j) + " " + strofreal(k))
            }
        }
    }
    if (level == "metric" | level == "scalar") {
        for (j = 1; j <= 6; j++) {
            L = L, ("tie lambda 1 F" + strofreal(j <= 3 ? 1 : 2) + " y" + strofreal(j) +
                    " lambda 2 F" + strofreal(j <= 3 ? 1 : 2) + " y" + strofreal(j))
        }
    }
    return(tfb_join(L))
}

string scalar tfb_key(struct dwlsmv_parmap scalar par, real scalar i)
{
    return(strofreal(par.ptype[i]) + ":" + strofreal(par.grp[i]) + ":" +
           strofreal(par.row[i]) + ":" + strofreal(par.col[i]))
}

real scalar tfb_get(string matrix M, string scalar level, string scalar key, real scalar vc)
{
    real scalar i

    for (i = 1; i <= rows(M); i++) {
        if (M[i, 1] == level & M[i, 2] == key) {
            return(strtoreal(M[i, vc]))
        }
    }
    return(.)
}

transmorphic tfb_vmap(string scalar level)
{
    external string matrix VCOV
    transmorphic A
    real scalar k

    A = asarray_create()
    for (k = 1; k <= rows(VCOV); k++) {
        if (VCOV[k, 1] == level) {
            asarray(A, VCOV[k, 2] + "|" + VCOV[k, 3], strtoreal(VCOV[k, 4]))
        }
    }
    return(A)
}

real scalar tfb_vget(transmorphic A, string scalar ki, string scalar kj)
{
    if (asarray_contains(A, ki + "|" + kj)) {
        return(asarray(A, ki + "|" + kj))
    }
    return(.)
}

real scalar tfb_fitval(string scalar level, real scalar col)
{
    external string matrix FITS
    real scalar i

    for (i = 1; i <= rows(FITS); i++) {
        if (FITS[i, 1] == level) {
            return(strtoreal(FITS[i, col]))
        }
    }
    return(.)
}

// the row of factor f's first free loading in group g
real scalar tfb_anchor(struct dwlsmv_parmap scalar par, real scalar g, real scalar f)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1 & par.grp[i] == g & par.col[i] == f & par.free[i]) {
            return(i)
        }
    }
    return(0)
}

// the sign multiplier of row i, from sgn[g, f]
real scalar tfb_mult(struct dwlsmv_parmap scalar par, real scalar i, real matrix sgn)
{
    real scalar g

    g = par.grp[i]
    if (par.ptype[i] == 1) return(sgn[g, par.col[i]])
    if (par.ptype[i] == 7) return(sgn[g, par.row[i]])
    if (par.ptype[i] == 3) return(sgn[g, par.row[i]] * sgn[g, par.col[i]])
    return(1)
}

real scalar tfb_rel(real scalar a, real scalar b)
{
    return(abs(a - b) / max((abs(b), 1e-300)))
}

// (a, b) from dwls_inference() at theta: Delta rebuilt at theta, everything
// else (Gamma, pi, D_stat) the fit's
real rowvector tfb_ab(struct dwlsmv_fit scalar fit, real colvector th,
                      struct dwlsmv_parmap scalar par, struct dwlsmv_statmap vector maps)
{
    struct dwlsmv_fit scalar fl

    fl = fit
    fl.theta = th
    fl.Delta = dwls_delta(th, par, maps)
    dwls_inference(fl)
    return((fl.shift_a, fl.shift_b))
}

// The attribution rule of check 4, exactly as the header defines it. Returns
// (F_o, F_l, (i), (ii) a rel, (ii) b abs, (iii) a gap, a bound, b gap, b bound,
//  a passes, b passes). Anything missing or non-finite fails.
real rowvector tfb_attrib(struct dwlsmv_fit scalar fit, struct dwlsmv_parmap scalar par,
                          struct dwlsmv_statmap vector maps, struct dwlsmv_stage1 vector sts,
                          real matrix sgn, string scalar level, real scalar alav, real scalar blav)
{
    external string matrix PARS
    struct dwlsmv_statmap vector mp
    real colvector tho, thl, d, s, w, tp, tm
    real rowvector abo, abl, abp, abm, gr, P
    real matrix H
    real scalar i, k, h, Fo, Fl, N, c1, ea, eb, ga, gb, Ba, Bb, pa, pb

    tho = fit.theta
    thl = tho
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) thl[par.fidx[i]] = tfb_mult(par, i, sgn) * tfb_get(PARS, level, tfb_key(par, i), 3)
    }
    d = tho - thl
    dwlsmv_fit_stack(sts, s, w, N, mp)
    dwlsmv_fdwls(0, tho', s, w, par, mp, Fo, gr, H)
    dwlsmv_fdwls(0, thl', s, w, par, mp, Fl, gr, H)
    c1  = (Fo <= Fl * (1 + 1e-12))
    abo = tfb_ab(fit, tho, par, maps)
    abl = tfb_ab(fit, thl, par, maps)
    ea  = abs(abl[1] - alav) / abs(alav)
    eb  = abs(abl[2] - blav)
    P   = (0, 0)
    for (k = 1; k <= rows(tho); k++) {
        if (d[k] == 0) continue
        h  = 1e-6 * max((1, abs(tho[k])))
        tp = tho
        tm = tho
        tp[k] = tho[k] + h
        tm[k] = tho[k] - h
        abp = tfb_ab(fit, tp, par, maps)
        abm = tfb_ab(fit, tm, par, maps)
        P = P + abs((abp - abm) / (2 * h)) * abs(d[k])
    }
    ga = abs(abo[1] - alav)
    gb = abs(abo[2] - blav)
    Ba = 2 * P[1] + 1e-8 * abs(alav)
    Bb = 2 * P[2] + 1e-8
    if (missing((Fo, Fl, abo, abl, P, ga, gb)) | hasmissing(d)) {
        return((Fo, Fl, 0, ea, eb, ga, Ba, gb, Bb, 0, 0))
    }
    pa = c1 & ea <= 1e-8 & ga <= Ba
    pb = c1 & eb <= 1e-8 & gb <= Bb
    return((Fo, Fl, c1, ea, eb, ga, Ba, gb, Bb, pa, pb))
}

// e(V) = ACOV / (N - G) from dwls_inference() at theta (Delta rebuilt there;
// Gamma, pi, D_stat the fit's)
real matrix tfb_Vat(struct dwlsmv_fit scalar fit, real colvector th, struct dwlsmv_parmap scalar par,
                    struct dwlsmv_statmap vector maps, real scalar NG)
{
    struct dwlsmv_fit scalar fl

    fl = fit
    fl.theta = th
    fl.Delta = dwls_delta(th, par, maps)
    dwls_inference(fl)
    return(fl.acov :/ NG)
}

// Spec 7.2.1 applied to check 6. V is ours at theta_o, VL lavaan's (theta order,
// sign-matched). Elements: SE_i = sqrt(V_ii); covariance V_ij, i < j. Returns
// (F_o, F_l, cond 1, cond 2 worst ratio |x(theta_l) - x_L| / tau, cond 2 holds,
//  elements missing the plain tier, of those attributed by cond 3, worst cond-3
//  ratio gap / bound over the missing ones, check 6 passes). Anything missing
// or non-finite fails the rule.
real rowvector tfb_vattrib(struct dwlsmv_fit scalar fit, struct dwlsmv_parmap scalar par,
                           struct dwlsmv_statmap vector maps, struct dwlsmv_stage1 vector sts,
                           real matrix sgn, string scalar level, real matrix V, real matrix VL,
                           real scalar NG)
{
    external string matrix PARS
    struct dwlsmv_statmap vector mp
    real colvector tho, thl, d, s, w, tp, tm, xo, xl, xL, tS, PS, dS
    real matrix Vl, Vp, Vm, tC, PC, H
    real rowvector gr
    real scalar i, j, k, h, q, Fo, Fl, N, c1, w2, c2, nm, na, w3, miss, b3, g3, bad

    q   = rows(V)
    tho = fit.theta
    thl = tho
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) thl[par.fidx[i]] = tfb_mult(par, i, sgn) * tfb_get(PARS, level, tfb_key(par, i), 3)
    }
    d = tho - thl
    dwlsmv_fit_stack(sts, s, w, N, mp)
    dwlsmv_fdwls(0, tho', s, w, par, mp, Fo, gr, H)
    dwlsmv_fdwls(0, thl', s, w, par, mp, Fl, gr, H)
    c1 = (Fo <= Fl * (1 + 1e-12))
    Vl = tfb_Vat(fit, thl, par, maps, NG)

    // tiers at the reference (7.2.1): SE 1e-8 + 1e-5 x_L; cov 1e-8 + 1e-5 sqrt(VL_ii VL_jj)
    xo = sqrt(diagonal(V))
    xl = sqrt(diagonal(Vl))
    xL = sqrt(diagonal(VL))
    tS = 1e-8 :+ 1e-5 :* xL
    tC = 1e-8 :+ 1e-5 :* sqrt(diagonal(VL) * diagonal(VL)')

    // condition 2: every element, at theta_l
    w2 = max(abs(xl - xL) :/ tS)
    for (i = 1; i <= q; i++) {
        for (j = i + 1; j <= q; j++) {
            w2 = max((w2, abs(Vl[i, j] - VL[i, j]) / tC[i, j]))
        }
    }
    c2 = (w2 <= 1e-2)

    // P_e: central differences at theta_o, h_k = 1e-6 max(1, |theta_o,k|)
    PS = J(q, 1, 0)
    PC = J(q, q, 0)
    for (k = 1; k <= q; k++) {
        if (d[k] == 0) continue
        h  = 1e-6 * max((1, abs(tho[k])))
        tp = tho
        tm = tho
        tp[k] = tho[k] + h
        tm[k] = tho[k] - h
        Vp = tfb_Vat(fit, tp, par, maps, NG)
        Vm = tfb_Vat(fit, tm, par, maps, NG)
        dS = (sqrt(diagonal(Vp)) - sqrt(diagonal(Vm))) :/ (2 * h)
        PS = PS + abs(dS) :* abs(d[k])
        PC = PC + abs((Vp - Vm) :/ (2 * h)) :* abs(d[k])
    }

    // the plain check-6 tier at theta_o decides first; condition 3 on the misses
    nm = 0
    na = 0
    w3 = 0
    for (i = 1; i <= q; i++) {
        for (j = i; j <= q; j++) {
            if (i == j) {
                miss = dwlsmv_gapratio(sqrt(V[i, i]), sqrt(VL[i, i]), 1e-8, 1e-5) > 1
                g3   = abs(xo[i] - xL[i])
                b3   = 2 * PS[i] + abs(xl[i] - xL[i])
            }
            else {
                miss = abs(V[i, j] - VL[i, j]) / (1e-8 + 1e-5 * sqrt(V[i, i] * V[j, j])) > 1
                g3   = abs(V[i, j] - VL[i, j])
                b3   = 2 * PC[i, j] + abs(Vl[i, j] - VL[i, j])
            }
            if (miss) {
                nm++
                w3 = max((w3, g3 / b3))
                if (g3 <= b3) na++
            }
        }
    }
    bad = missing((Fo, Fl, w2, w3)) + missing(PS) + missing(PC) + missing(Vl) + missing(d)
    if (bad) {
        return((Fo, Fl, 0, w2, 0, nm, 0, w3, 0))
    }
    return((Fo, Fl, c1, w2, c2, nm, na, w3, (nm == 0) | (c1 & c2 & na == nm)))
}

void tfb_level(string scalar level, struct dwlsmv_stage1 vector sts,
               string colvector FAIL)
{
    external string matrix PARS
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    struct dwlsmv_fit scalar fit
    transmorphic VA
    real colvector th, s, w, r
    real matrix V, sgn
    real scalar G, N, i, j, g, f, a, v, worst, wc, wo, df, F1, F2, TB, TC, Bd,
                alav, blav, T3l, nfree, ki, kj, ok_a, ok_b
    real rowvector att, att6
    real matrix VL

    G    = 2
    maps = dwlsmv_statmap(2, 1)
    maps[1] = sts[1].map
    maps[2] = sts[2].map
    N    = sts[1].N + sts[2].N
    par  = dwls_identify(tfb_text(level), maps, sts)
    fit  = dwls_fit(sts, par)
    dwls_inference(fit)
    th   = fit.theta
    V    = fit.acov :/ (N - G)

    // ---- 1 -----------------------------------------------------------------
    df = tfb_fitval(level, 4)
    printf("{txt}  %-10s 1. converged %f, D_stat %f, lavaan df %f\n", level, fit.converged, fit.D_stat, df)
    if (fit.converged != 1 | fit.D_stat != df) {
        FAIL = FAIL \ (level + ": not converged or D_stat is not lavaan's df")
    }

    // sign equivalence per factor and group
    sgn = J(G, 2, 1)
    for (g = 1; g <= G; g++) {
        for (f = 1; f <= 2; f++) {
            i = tfb_anchor(par, g, f)
            v = tfb_get(PARS, level, tfb_key(par, i), 3)
            if (sign(v) != sign(th[par.fidx[i]])) {
                sgn[g, f] = -1
            }
        }
    }
    printf("{txt}  %-10s    lavaan sign per group and factor: g1 %f %f ; g2 %f %f\n", level,
           sgn[1, 1], sgn[1, 2], sgn[2, 1], sgn[2, 2])

    // ---- 2 -----------------------------------------------------------------
    worst = 0
    nfree = 0
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        nfree++
        v = tfb_get(PARS, level, tfb_key(par, i), 3)
        if (missing(v)) {
            FAIL = FAIL \ (level + ": lavaan has no estimate for " + tfb_key(par, i))
            continue
        }
        v = v * tfb_mult(par, i, sgn)
        worst = max((worst, dwlsmv_gapratio(th[par.fidx[i]], v, 1e-5, 1e-4)))
    }
    printf("{txt}  %-10s 2. theta vs lavaan: worst ratio %8.2e of the tier over %f free rows (npar %f)\n",
           level, worst, nfree, par.npar)
    if (worst > 1) {
        FAIL = FAIL \ (level + ": theta-hat disagrees with lavaan")
    }

    // ---- 3 -----------------------------------------------------------------
    F1 = tfb_fitval(level, 7)
    F2 = tfb_fitval(level, 8)
    TB = (N - G) * ((sts[1].N / N) * F1 + (sts[2].N / N) * F2)
    TC = tfb_fitval(level, 6)
    s  = sts[1].s \ sts[2].s
    w  = ((sts[1].N / N) :/ diagonal(sts[1].Gamma)) \ ((sts[2].N / N) :/ diagonal(sts[2].Gamma))
    r  = s - dwls_model(th, par, maps)
    Bd = (N - G) * (2 * quadcross(w, abs(r) :* (1e-8 :+ 1e-8 :* abs(s))) + 1e-8 * fit.fmin)
    printf("{txt}  %-10s 3. T_raw %16.10f ; lavaan (N-G) sum pi F %16.10f ratio %8.2e ; lavaan chisq (C) %16.10f ratio %8.2e\n",
           level, fit.T_raw, TB, abs(fit.T_raw - TB) / Bd, TC, abs(fit.T_raw - TC) / Bd)
    if (abs(fit.T_raw - TB) > Bd) {
        FAIL = FAIL \ (level + ": T_raw is not (N - G) F at lavaan's per-group F")
    }
    if (abs(fit.T_raw - TC) <= Bd) {
        FAIL = FAIL \ (level + ": T_raw also matches candidate C, so check 3 cannot tell B from C")
    }

    // ---- 4, 5 --------------------------------------------------------------
    alav = 1 / tfb_fitval(level, 10)
    blav = tfb_fitval(level, 11)
    T3l  = alav * TB + blav
    printf("{txt}  %-10s 4. a %14.10f (lavaan %14.10f, rel %8.2e) ; b %14.10f (lavaan %14.10f, abs %8.2e)\n",
           level, fit.shift_a, alav, tfb_rel(fit.shift_a, alav), fit.shift_b, blav, abs(fit.shift_b - blav))
    ok_a = tfb_rel(fit.shift_a, alav) <= 1e-5
    ok_b = abs(fit.shift_b - blav) <= 1e-6
    att  = tfb_attrib(fit, par, maps, sts, sgn, level, alav, blav)
    printf("{txt}  %-10s 4. rule: (i) F_o %20.15e F_l %20.15e holds %f ; (ii) a %8.2e b %8.2e ; (iii) a gap %8.2e <= %8.2e, b gap %8.2e <= %8.2e ; a %f b %f\n",
           level, att[1], att[2], att[3], att[4], att[5], att[6], att[7], att[8], att[9], att[10], att[11])
    if (!ok_a & !att[10]) {
        FAIL = FAIL \ (level + ": a disagrees with lavaan beyond its tier and the attribution rule")
    }
    if (!ok_b & !att[11]) {
        FAIL = FAIL \ (level + ": b disagrees with lavaan beyond its tier and the attribution rule")
    }
    if ((!ok_a & att[10]) | (!ok_b & att[11])) {
        printf("{txt}  %-10s 4. outside the plain tier, passed by the attribution rule\n", level)
    }
    printf("{txt}  %-10s 5. T3 %14.10f ; a_lav T_B,lav + b_lav %14.10f ; ratio %8.2e\n",
           level, fit.T3, T3l, dwlsmv_gapratio(fit.T3, T3l, 1e-6, 1e-5))
    if (dwlsmv_gapratio(fit.T3, T3l, 1e-6, 1e-5) > 1) {
        FAIL = FAIL \ (level + ": T3 disagrees with lavaan's scaled-shifted statistic at T_B")
    }

    // ---- 6 -----------------------------------------------------------------
    VA = tfb_vmap(level)
    wo = 0
    wc = 0
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        for (j = 1; j <= par.nall; j++) {
            if (!par.free[j]) continue
            ki = par.fidx[i]
            kj = par.fidx[j]
            v  = tfb_vget(VA, tfb_key(par, i), tfb_key(par, j))
            if (missing(v)) {
                FAIL = FAIL \ (level + ": lavaan vcov has no entry " + tfb_key(par, i) + ", " + tfb_key(par, j))
                continue
            }
            a = tfb_mult(par, i, sgn) * tfb_mult(par, j, sgn)
            v = a * v
            if (i == j) {
                wo = max((wo, dwlsmv_gapratio(sqrt(V[ki, ki]), sqrt(v), 1e-8, 1e-5)))
            }
            else {
                wc = max((wc, abs(V[ki, kj] - v) / (1e-8 + 1e-5 * sqrt(V[ki, ki] * V[kj, kj]))))
            }
        }
    }
    printf("{txt}  %-10s 6. e(V): SEs worst ratio %8.2e ; covariances worst ratio %8.2e\n", level, wo, wc)
    // spec 7.2.1, evaluated and printed at every level; decides only the misses
    VL = J(par.npar, par.npar, .)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        for (j = 1; j <= par.nall; j++) {
            if (!par.free[j]) continue
            VL[par.fidx[i], par.fidx[j]] = tfb_mult(par, i, sgn) * tfb_mult(par, j, sgn) *
                                           tfb_vget(VA, tfb_key(par, i), tfb_key(par, j))
        }
    }
    att6 = tfb_vattrib(fit, par, maps, sts, sgn, level, V, VL, N - G)
    printf("{txt}  %-10s 6. 7.2.1: (1) F_o %20.15e F_l %20.15e holds %f ; (2) worst %8.2e of tau, holds %f ; (3) %f element(s) missed the plain tier, %f attributed, worst gap/bound %8.2e ; check 6 %f\n",
           level, att6[1], att6[2], att6[3], att6[4], att6[5], att6[6], att6[7], att6[8], att6[9])
    if (!att6[9]) {
        FAIL = FAIL \ (level + ": e(V) = ACOV / (N - G) disagrees with lavaan's vcov beyond its tier and spec 7.2.1")
    }
    if ((wo > 1 | wc > 1) & att6[9]) {
        printf("{txt}  %-10s 6. outside the plain tier, passed by spec 7.2.1 (attribution, not plain agreement)\n", level)
    }

    // ---- b: exactly one beta row, free, in group 1 ---------------------------
    nfree = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] != 3) continue
        nfree++
        if (!par.free[i] | par.grp[i] != 1) {
            FAIL = FAIL \ (level + ": a beta row is not free, or not in group 1")
            continue
        }
        ki = par.fidx[i]
        printf("{txt}  %-10s b. beta g%f: %10.6f (lavaan %10.6f) ; SE %10.6f (lavaan %10.6f)\n", level,
               par.grp[i], th[ki], tfb_mult(par, i, sgn) * tfb_get(PARS, level, tfb_key(par, i), 3),
               sqrt(V[ki, ki]), sqrt(tfb_vget(VA, tfb_key(par, i), tfb_key(par, i))))
    }
    if (nfree != 1) {
        FAIL = FAIL \ (level + ": expected exactly one beta row (group 1)")
    }
}

// one pair's polychoric log-likelihood score per observation at rho, our
// thresholds (as test_fit_groups.do)
real scalar tfb_score(real colvector yi, real colvector yj, real scalar Ki,
                      real scalar Kj, real scalar rho)
{
    real colvector ti, tj
    real matrix n

    ti = dwlsmv_thresholds(yi, Ki)
    tj = dwlsmv_thresholds(yj, Kj)
    n  = dwlsmv_crosstab(yi, yj, Ki, Kj)
    return(sum(n :* (dwlsmv_dcellprob_rho(ti, tj, rho) :/ dwlsmv_cellprob(ti, tj, rho))) / rows(yi))
}

real scalar tfb_run()
{
    external string matrix SAMP
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    string colvector FAIL, lev
    real colvector KM, sel, yi, yj
    real scalar g, i, k, v, worst, n, gr, nat, Ki, Kj, So, Sl, Sp

    FAIL = J(0, 1, "")
    KM   = (4 \ 4 \ 2 \ 4 \ 2 \ 2)
    sts  = dwlsmv_stage1(2, 1)
    for (g = 1; g <= 2; g++) {
        sel    = selectindex(YALL[., 1] :== g)
        sts[g] = dwls_sample_stats(YALL[sel, 2..7], KM)
    }

    // ---- s -----------------------------------------------------------------
    worst = 0
    n     = 0
    nat   = 0
    for (g = 1; g <= 2; g++) {
        sel = selectindex(YALL[., 1] :== g)
        for (i = 1; i <= sts[g].map.D; i++) {
            v = .
            for (k = 1; k <= rows(SAMP); k++) {
                if (strtoreal(SAMP[k, 1]) == g & SAMP[k, 2] == sts[g].map.label[i]) {
                    v = strtoreal(SAMP[k, 3])
                    break
                }
            }
            if (missing(v)) {
                FAIL = FAIL \ ("Stage 1: lavaan has no statistic " + sts[g].map.label[i])
                continue
            }
            n++
            gr = dwlsmv_gapratio(sts[g].s[i], v, 1e-8, 1e-8)
            if (gr > 1 & sts[g].map.stat_type[i] == 2) {
                yi = YALL[sel, 1 + sts[g].map.var_i[i]]
                yj = YALL[sel, 1 + sts[g].map.var_j[i]]
                Ki = KM[sts[g].map.var_i[i]]
                Kj = KM[sts[g].map.var_j[i]]
                So = tfb_score(yi, yj, Ki, Kj, sts[g].s[i])
                Sl = tfb_score(yi, yj, Ki, Kj, v)
                Sp = (tfb_score(yi, yj, Ki, Kj, sts[g].s[i] + 1e-5) -
                      tfb_score(yi, yj, Ki, Kj, sts[g].s[i] - 1e-5)) / 2e-5
                printf("{txt}  s. g%f %s: gap %9.2e ; |S|/N ours %9.2e lavaan's %9.2e ; Newton bound 2|S_L|/|S'| %9.2e ; signs %f %f\n",
                       g, sts[g].map.label[i], v - sts[g].s[i], abs(So), abs(Sl), 2 * abs(Sl) / abs(Sp),
                       sign(v - sts[g].s[i]), sign(Sl / Sp))
                if (abs(So) > 1e-12 | abs(v - sts[g].s[i]) > 2 * abs(Sl) / abs(Sp) |
                    sign(v - sts[g].s[i]) != sign(Sl / Sp)) {
                    FAIL = FAIL \ ("Stage 1: " + sts[g].map.label[i] + " in group " + strofreal(g) +
                                    " misses lavaan by more than lavaan's own Newton step explains")
                }
                nat++
                continue
            }
            worst = max((worst, gr))
        }
    }
    printf("{txt}  s. Stage 1 per group vs lavaan: worst ratio %8.2e over %f statistics; N = %f, %f\n",
           worst, n, sts[1].N, sts[2].N)
    printf("{txt}  s. %f correlation(s) outside the tier, each attributed to lavaan's Newton step above\n", nat)
    if (worst > 1) {
        FAIL = FAIL \ "Stage 1 per group disagrees with lavaan's per-group statistics"
    }

    lev = ("configural" \ "thresholds" \ "metric" \ "scalar")
    for (i = 1; i <= rows(lev); i++) {
        tfb_level(lev[i], sts, FAIL)
    }

    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}  FAILED: %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nf", strofreal(tfb_run()))
if `nf' {
    display as error "test_fit_groups_betag1: `nf' failure(s)"
    exit 9
}
display as result "test_fit_groups_betag1: with beta free in group 1 and absent in group 2 (design M), the pooled fit and inference reproduce lavaan's two-group fits at all four levels (theta incl. beta, T_raw via (N-G) F, a, b, T3, e(V))"
exit 0
