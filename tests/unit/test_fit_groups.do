*! test_fit_groups.do -- the pooled fit and inference under groups (spec 4.8
*! O14, O16; derivation @sec-multigroup-scaling MG1-MG6), against lavaan 0.7.2.
*!
*! Executed under StataMP 18 on 2026-10-10; green on the fourth run (the first three runs and what changed are recorded under checks s, 8 and 9).
*!
*! Design E (K = 4, 4, 4, 2, 2), two groups of UNEQUAL size (1200, 2800), the
*! four levels, written as canonical 'groups' texts by hand (Module 11 is not
*! involved). Reference: fixtures/expected/lavaan_mg_E_*.csv from
*! fixtures/generate/98_lavaan_multigroup.R. Comparison rules (user,
*! 2026-10-10): estimates after sign equivalence; T_raw through the stated
*! conversion; e(V) in the same parameterization and order; T3 as the actual
*! scaled-and-shifted statistic, with a and b compared; the configural
*! variance oracle from Stage 1 alone. Tolerances are spec 7.2's rows unless a
*! check says otherwise.
*!
*! Checks:
*!   s   Stage 1 per group (dwls_sample_stats on each group's rows) against
*!       lavaan's per-group sample statistics, by label: thresholds and
*!       correlations tier, 1e-8 abs + 1e-8 rel. A correlation outside that
*!       tier passes only if its whole gap is lavaan's unfinished Newton step
*!       (user, 2026-10-10, option (a); bound justified and fixed before the
*!       rerun): with S the pair's score per observation, (i) |S(our rho)| <=
*!       1e-12; (ii) |rho_L - rho| <= 2 |S(rho_L)| / |S'(rho)|, S' a central
*!       difference of the analytic score, h = 1e-5 (mean value theorem; the
*!       factor 2 allows S' to change by up to half across the gap, which at
*!       a gap near 1e-8 it cannot); (iii) the gap has the sign of lavaan's
*!       step, sign(S(rho_L) / S'(rho)). FIRST RUN (2026-10-10): failed at the
*!       plain tier, g1 y1~~y3 gap 2.4e-8 (1.65x), with |S|/N 4.7e-14 at ours
*!       and 2.9e-8 at lavaan's: lavaan stops earlier, as test_gamma.do
*!       recorded on case E. SECOND RUN: (i) 4.7e-14 and (ii) held (the gap,
*!       2.41e-8, is exactly |S_L| / |S'|), (iii) failed: as first fixed it read
*!       sign(S(rho_L) / -S'), contradicting the mean-value step it came from
*!       (rho_L - rho = S(rho_L) / S'). A slip in the preregistered assertion,
*!       corrected to the form above with the user's approval (2026-10-10).
*!   1  converged; D_stat = lavaan's df = Result I4 (10, 13, 15, 19), exact.
*!   2   theta-hat against lavaan by (ptype, group, row, col), every free table
*!       row (tied members each), after sign equivalence: within a group the
*!       factor's sign is matched on that group's first loading (in a tied table
*!       it is one sign; in a configural one each group has its own). Theta
*!       tier, 1e-5 abs + 1e-4 rel.
*!   3   T_raw = (N - G) F-hat against lavaan's converted: (N - G) sum_g pi_g F_g
*!       from lavaan's own per-group F_g. Tolerance propagated as test_fit.do
*!       check 5 does (s and Gamma moving by their tiers), pooled. And the
*!       unconverted lavaan chi-square (candidate C, sum_g (N_g - 1) F_g) must
*!       NOT match within the same bound: the check can tell B from C.
*!   4   a and b against lavaan's (a = 1 / scaling_factor, b = shift): a
*!       1e-5 rel, b 1e-6 abs.
*!   5   T3 = a T_raw + b against a_lav T_B,lav + b_lav: 1e-6 abs + 1e-5 rel.
*!   6   e(V) = ACOV / (N - G) against lavaan's vcov() by key, every pair of
*!       free table rows: SEs at the SE tier (1e-8 abs, 1e-5 rel); off-diagonal
*!       entries within 1e-8 + 1e-5 sqrt(V_ii V_jj) (the SE tier carried to a
*!       covariance; chosen here before the first run). Sign-matched as in 2.
*!   7   configural only: the threshold block of e(V) against the Stage-1-only
*!       oracle Gamma_g,kk / (pi_g (N - G)) (MG6's block identity): 1e-9
*!       relative (chosen before the first run: exact algebra through one
*!       cholinv), and that oracle against lavaan's at the SE tier.
*!   8   negative control (MG2's Tell), tied levels: traces from
*!       blockdiag(Gamma_g), no pi, give an a that differs from the correct one
*!       by more than its tier. Configural: the same a to 1e-12 (with no tie U
*!       is block diagonal and pi_g cancels block by block).
*!   9   O16 metamorphic (user's check B), loadings and intercepts levels:
*!       reflect the factor in every group at once (lambda, alpha; each theta
*!       entry once): sigma bitwise unchanged; signcanon() returns theta-hat
*!       bitwise; fits started at theta-hat and at its reflection return the
*!       same canonical theta (1e-10 rel) and the same F, T_raw, T3, a, b,
*!       D_stat (1e-10 rel), tied members still sharing their entry.
*!
*!   REDESIGNED AFTER THE FIRST RUN (2026-10-10), both predictions of mine that
*!   were wrong, kept here: check 8 asserted "must differ" at the configural
*!   level too and failed there (rel 1.2e-14), which is the algebra above, not
*!   a defect; check 9 compared the refit from the reflection against the
*!   ORIGINAL fit, which arrived by a different path (ULS, then DWLS), and
*!   failed on theta at 1e-10 (ratios 3.9e2, 1.4e3) while F agreed to 1e-14:
*!   two separately converged runs agree to the optimizer's stopping precision,
*!   not to 1e-10. Both runs now start from the same point up to reflection.
*!  10   negative control for O16: at the intercepts level negate only group
*!       2's alpha (what a one-group reflection can reach once lambda is tied):
*!       sigma moves by more than 1e-3. Not a symmetry.
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
        display as error "test_fit_groups: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/expected/lavaan_mg_E_fits.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: FITS = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_E_parameters.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: PARS = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_E_vcov.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: VCOV = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_E_sampstat.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: SAMP = st_sdata(., .)
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
mata: YALL = st_data(., .)

mata:

string scalar tfg_join(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The canonical text of design E at a level (spec 4.8; the level table of
// @sec-invariance-id, as fixtures/generate/98_lavaan_multigroup.R fits it).
string scalar tfg_text(string scalar level)
{
    string rowvector L
    real rowvector C
    real scalar j, k

    L = ("groups 2", "group 1", "factor F", "ind y1 y2 y3 y4 y5")
    for (j = 1; j <= 5; j++) {
        L = L, ("lambda F y" + strofreal(j) + " free")
    }
    L = L, ("metric unitvariance F", "group 2", "factor F", "ind y1 y2 y3 y4 y5")
    for (j = 1; j <= 5; j++) {
        L = L, ("lambda F y" + strofreal(j) + " free")
    }
    if (level == "configural" | level == "thresholds") {
        L = L, "psi F F fixed 1"
    }
    else {
        L = L, "psi F F free"
    }
    if (level == "thresholds" | level == "loadings") {
        for (j = 1; j <= 5; j++) {
            L = L, ("nu y" + strofreal(j) + " free")
        }
    }
    if (level == "thresholds") {
        for (j = 1; j <= 3; j++) {
            L = L, ("omega y" + strofreal(j) + " free")
        }
    }
    if (level == "loadings" | level == "intercepts") {
        for (j = 1; j <= 5; j++) {
            L = L, ("omega y" + strofreal(j) + " free")
        }
    }
    if (level == "intercepts") {
        L = L, "alpha F free"
    }
    if (level != "configural") {
        C = (3, 3, 3, 1, 1)
        for (j = 1; j <= 5; j++) {
            for (k = 1; k <= C[j]; k++) {
                L = L, ("tie tau 1 y" + strofreal(j) + " " + strofreal(k) +
                        " tau 2 y" + strofreal(j) + " " + strofreal(k))
            }
        }
    }
    if (level == "loadings" | level == "intercepts") {
        for (j = 1; j <= 5; j++) {
            L = L, ("tie lambda 1 F y" + strofreal(j) + " lambda 2 F y" + strofreal(j))
        }
    }
    return(tfg_join(L))
}

string scalar tfg_key(struct dwlsmv_parmap scalar par, real scalar i)
{
    return(strofreal(par.ptype[i]) + ":" + strofreal(par.grp[i]) + ":" +
           strofreal(par.row[i]) + ":" + strofreal(par.col[i]))
}

// lavaan's value for (level, key) in a fixture whose first two columns are
// level and key and whose value column is vc; missing if absent
real scalar tfg_get(string matrix M, string scalar level, string scalar key, real scalar vc)
{
    real scalar i

    for (i = 1; i <= rows(M); i++) {
        if (M[i, 1] == level & M[i, 2] == key) {
            return(strtoreal(M[i, vc]))
        }
    }
    return(.)
}

// lavaan's vcov for one level as an asarray keyed "key_i|key_j"
transmorphic tfg_vmap(string scalar level)
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

real scalar tfg_vget(transmorphic A, string scalar ki, string scalar kj)
{
    if (asarray_contains(A, ki + "|" + kj)) {
        return(asarray(A, ki + "|" + kj))
    }
    return(.)
}

real scalar tfg_fitval(string scalar level, real scalar col)
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

// sign of factor F in group g: the first free loading's row of group g
real scalar tfg_anchor(struct dwlsmv_parmap scalar par, real scalar g)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1 & par.grp[i] == g & par.free[i]) {
            return(i)
        }
    }
    return(0)
}

// the reflection of factor F in every group: lambda and alpha entries, once each
real colvector tfg_reflect(real colvector th, struct dwlsmv_parmap scalar par)
{
    real colvector neg, out
    real scalar i

    neg = J(par.npar, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            if (par.ptype[i] == 1 | par.ptype[i] == 7) {
                neg[par.fidx[i]] = 1
            }
        }
    }
    out = th
    for (i = 1; i <= par.npar; i++) {
        if (neg[i]) {
            out[i] = -out[i]
        }
    }
    return(out)
}

real scalar tfg_rel(real scalar a, real scalar b)
{
    return(abs(a - b) / max((abs(b), 1e-300)))
}

void tfg_level(string scalar level, struct dwlsmv_stage1 vector sts,
               string colvector FAIL)
{
    external string matrix PARS
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    struct dwlsmv_fit scalar fit, fr, f0
    transmorphic VA
    real colvector th, thr, sig, sigr, s, w, r, sgn
    real matrix V, Om0
    real rowvector tr0, ab0
    real scalar G, N, i, j, g, a, b, v, worst, wc, df, F1, F2, TB, TC, Bd,
                alav, blav, T3l, nfree, k, o, orac, wo, wl, d, pig, ki, kj

    G    = 2
    maps = dwlsmv_statmap(2, 1)
    maps[1] = sts[1].map
    maps[2] = sts[2].map
    N    = sts[1].N + sts[2].N
    par  = dwls_identify(tfg_text(level), maps, sts)
    fit  = dwls_fit(sts, par)
    dwls_inference(fit)
    th   = fit.theta
    V    = fit.acov :/ (N - G)

    // ---- 1 -----------------------------------------------------------------
    df = tfg_fitval(level, 4)
    printf("{txt}  %-10s 1. converged %f, D_stat %f, lavaan df %f\n", level, fit.converged, fit.D_stat, df)
    if (fit.converged != 1 | fit.D_stat != df) {
        FAIL = FAIL \ (level + ": not converged or D_stat is not lavaan's df")
    }

    // sign equivalence: per group, from that group's first loading
    sgn = J(G, 1, 1)
    for (g = 1; g <= G; g++) {
        i = tfg_anchor(par, g)
        v = tfg_get(PARS, level, tfg_key(par, i), 3)
        if (sign(v) != sign(th[par.fidx[i]])) {
            sgn[g] = -1
        }
    }
    printf("{txt}  %-10s    lavaan sign per group: %f %f\n", level, sgn[1], sgn[2])

    // ---- 2 -----------------------------------------------------------------
    worst = 0
    nfree = 0
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        nfree++
        v = tfg_get(PARS, level, tfg_key(par, i), 3)
        if (missing(v)) {
            FAIL = FAIL \ (level + ": lavaan has no estimate for " + tfg_key(par, i))
            continue
        }
        if (par.ptype[i] == 1 | par.ptype[i] == 7) {
            v = v * sgn[par.grp[i]]
        }
        worst = max((worst, dwlsmv_gapratio(th[par.fidx[i]], v, 1e-5, 1e-4)))
    }
    printf("{txt}  %-10s 2. theta vs lavaan: worst ratio %8.2e of the tier over %f free rows (npar %f)\n",
           level, worst, nfree, par.npar)
    if (worst > 1) {
        FAIL = FAIL \ (level + ": theta-hat disagrees with lavaan")
    }

    // ---- 3 -----------------------------------------------------------------
    F1 = tfg_fitval(level, 7)
    F2 = tfg_fitval(level, 8)
    TB = (N - G) * ((sts[1].N / N) * F1 + (sts[2].N / N) * F2)
    TC = tfg_fitval(level, 6)
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
    alav = 1 / tfg_fitval(level, 10)
    blav = tfg_fitval(level, 11)
    T3l  = alav * TB + blav
    printf("{txt}  %-10s 4. a %14.10f (lavaan %14.10f, rel %8.2e) ; b %14.10f (lavaan %14.10f, abs %8.2e)\n",
           level, fit.shift_a, alav, tfg_rel(fit.shift_a, alav), fit.shift_b, blav, abs(fit.shift_b - blav))
    if (tfg_rel(fit.shift_a, alav) > 1e-5 | abs(fit.shift_b - blav) > 1e-6) {
        FAIL = FAIL \ (level + ": a or b disagrees with lavaan")
    }
    printf("{txt}  %-10s 5. T3 %14.10f ; a_lav T_B,lav + b_lav %14.10f ; ratio %8.2e\n",
           level, fit.T3, T3l, dwlsmv_gapratio(fit.T3, T3l, 1e-6, 1e-5))
    if (dwlsmv_gapratio(fit.T3, T3l, 1e-6, 1e-5) > 1) {
        FAIL = FAIL \ (level + ": T3 disagrees with lavaan's scaled-shifted statistic at T_B")
    }

    // ---- 6 -----------------------------------------------------------------
    // every free table row's lavaan key; the fixture's vcov is over lavaan's rows
    VA = tfg_vmap(level)
    wo = 0
    wc = 0
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        for (j = 1; j <= par.nall; j++) {
            if (!par.free[j]) continue
            ki = par.fidx[i]
            kj = par.fidx[j]
            v  = tfg_vget(VA, tfg_key(par, i), tfg_key(par, j))
            if (missing(v)) {
                FAIL = FAIL \ (level + ": lavaan vcov has no entry " + tfg_key(par, i) + ", " + tfg_key(par, j))
                continue
            }
            a = 1
            if (par.ptype[i] == 1 | par.ptype[i] == 7) a = a * sgn[par.grp[i]]
            if (par.ptype[j] == 1 | par.ptype[j] == 7) a = a * sgn[par.grp[j]]
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
    if (wo > 1 | wc > 1) {
        FAIL = FAIL \ (level + ": e(V) = ACOV / (N - G) disagrees with lavaan's vcov")
    }

    // ---- 7 -----------------------------------------------------------------
    if (level == "configural") {
        wo = 0
        wl = 0
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] != 4 | !par.free[i]) continue
            g    = par.grp[i]
            pig  = sts[g].N / N
            // the statistic's position in group g's s: threshold (row, col)
            o    = sts[g].map.Coff[par.row[i]] + par.col[i]
            orac = sts[g].Gamma[o, o] / (pig * (N - G))
            ki   = par.fidx[i]
            wo   = max((wo, tfg_rel(V[ki, ki], orac)))
            v    = tfg_vget(VA, tfg_key(par, i), tfg_key(par, i))
            wl   = max((wl, dwlsmv_gapratio(sqrt(orac), sqrt(v), 1e-8, 1e-5)))
        }
        printf("{txt}  %-10s 7. threshold block vs Gamma_kk / (pi_g (N - G)): worst rel %8.2e ; oracle vs lavaan SE ratio %8.2e\n",
               level, wo, wl)
        if (wo > 1e-9 | wl > 1) {
            FAIL = FAIL \ (level + ": the configural threshold block is not MG6's oracle")
        }
    }

    // ---- 8 -----------------------------------------------------------------
    d   = sts[1].map.D
    Om0 = J(2 * d, 2 * d, 0)
    Om0[1::d, 1::d] = sts[1].Gamma
    Om0[(d + 1)::(2 * d), (d + 1)::(2 * d)] = sts[2].Gamma
    tr0 = dwlsmv_traces(fit.Delta, Om0)
    ab0 = dwlsmv_scaleshift(tr0[1], tr0[2], fit.D_stat)
    if (level == "configural") {
        printf("{txt}  %-10s 8. a without pi %14.10f vs %14.10f : rel %8.2e (no tie: pi cancels, must be below 1e-12)\n",
               level, ab0[1], fit.shift_a, tfg_rel(ab0[1], fit.shift_a))
        if (tfg_rel(ab0[1], fit.shift_a) > 1e-12) {
            FAIL = FAIL \ (level + ": with no tie the traces still depend on pi_g")
        }
    }
    else {
        printf("{txt}  %-10s 8. a without pi %14.10f vs %14.10f : rel %8.2e (must exceed 1e-5)\n",
               level, ab0[1], fit.shift_a, tfg_rel(ab0[1], fit.shift_a))
        if (tfg_rel(ab0[1], fit.shift_a) <= 1e-5) {
            FAIL = FAIL \ (level + ": traces without pi_g give the same a; check 4 cannot see the weighting")
        }
    }

    // ---- 9, 10 -------------------------------------------------------------
    if (level == "loadings" | level == "intercepts") {
        thr  = tfg_reflect(th, par)
        sig  = dwls_model(th, par, maps)
        sigr = dwls_model(thr, par, maps)
        f0   = dwlsmv_fit_from(sts, par, th)
        dwls_inference(f0)
        fr   = dwlsmv_fit_from(sts, par, thr)
        dwls_inference(fr)
        printf("{txt}  %-10s 9. reflection: sigma equal %f ; signcanon(reflected) == theta %f ; refit theta ratio %8.2e\n",
               level, sig == sigr, dwlsmv_signcanon(thr, par) == th, dwlsmv_gapratio(fr.theta, f0.theta, 0, 1e-10))
        printf("{txt}  %-10s    refit F %8.2e T_raw %8.2e T3 %8.2e a %8.2e b %8.2e (rel) ; D_stat %f %f\n",
               level, tfg_rel(fr.fmin, f0.fmin), tfg_rel(fr.T_raw, f0.T_raw), tfg_rel(fr.T3, f0.T3),
               tfg_rel(fr.shift_a, f0.shift_a), tfg_rel(fr.shift_b, f0.shift_b), fr.D_stat, f0.D_stat)
        if (sig != sigr) {
            FAIL = FAIL \ (level + ": sigma moved under the all-groups reflection")
        }
        if (dwlsmv_signcanon(thr, par) != th) {
            FAIL = FAIL \ (level + ": signcanon() does not return theta-hat from its reflection")
        }
        if (dwlsmv_gapratio(fr.theta, f0.theta, 0, 1e-10) > 1 | tfg_rel(fr.fmin, f0.fmin) > 1e-10 |
            tfg_rel(fr.T_raw, f0.T_raw) > 1e-10 | tfg_rel(fr.T3, f0.T3) > 1e-10 |
            tfg_rel(fr.shift_a, f0.shift_a) > 1e-10 | tfg_rel(fr.shift_b, f0.shift_b) > 1e-10 |
            fr.D_stat != f0.D_stat) {
            FAIL = FAIL \ (level + ": a fit from the reflected start is not the same fit")
        }
        if (level == "intercepts") {
            thr = th
            for (i = 1; i <= par.nall; i++) {
                if (par.ptype[i] == 7 & par.grp[i] == 2 & par.free[i]) {
                    thr[par.fidx[i]] = -thr[par.fidx[i]]
                }
            }
            sigr = dwls_model(thr, par, maps)
            printf("{txt}  %-10s 10. group 2's alpha alone negated: max |sigma change| %8.2e (must exceed 1e-3)\n",
                   level, max(abs(sigr - sig)))
            if (max(abs(sigr - sig)) <= 1e-3) {
                FAIL = FAIL \ (level + ": negating group 2's alpha alone left sigma unchanged")
            }
        }
    }
}

// one pair's polychoric log-likelihood score per observation at rho, our
// thresholds (the quantity Module 2 drives to zero)
real scalar tfg_score(real colvector yi, real colvector yj, real scalar Ki,
                      real scalar Kj, real scalar rho)
{
    real colvector ti, tj
    real matrix n

    ti = dwlsmv_thresholds(yi, Ki)
    tj = dwlsmv_thresholds(yj, Kj)
    n  = dwlsmv_crosstab(yi, yj, Ki, Kj)
    return(sum(n :* (dwlsmv_dcellprob_rho(ti, tj, rho) :/ dwlsmv_cellprob(ti, tj, rho))) / rows(yi))
}

real scalar tfg_run()
{
    external string matrix SAMP
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    string colvector FAIL, lev
    real colvector KE, sel, yi, yj
    real scalar g, i, k, v, worst, n, gr, nat, Ki, Kj, So, Sl, Sp

    FAIL = J(0, 1, "")
    KE   = (4 \ 4 \ 4 \ 2 \ 2)
    sts  = dwlsmv_stage1(2, 1)
    for (g = 1; g <= 2; g++) {
        sel    = selectindex(YALL[., 1] :== g)
        sts[g] = dwls_sample_stats(YALL[sel, 2..6], KE)
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
                // the user's option (a): attribute the gap to lavaan's Newton step
                yi = YALL[sel, 1 + sts[g].map.var_i[i]]
                yj = YALL[sel, 1 + sts[g].map.var_j[i]]
                Ki = KE[sts[g].map.var_i[i]]
                Kj = KE[sts[g].map.var_j[i]]
                So = tfg_score(yi, yj, Ki, Kj, sts[g].s[i])
                Sl = tfg_score(yi, yj, Ki, Kj, v)
                Sp = (tfg_score(yi, yj, Ki, Kj, sts[g].s[i] + 1e-5) -
                      tfg_score(yi, yj, Ki, Kj, sts[g].s[i] - 1e-5)) / 2e-5
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

    lev = ("configural" \ "thresholds" \ "loadings" \ "intercepts")
    for (i = 1; i <= rows(lev); i++) {
        tfg_level(lev[i], sts, FAIL)
    }

    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}  FAILED: %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nf", strofreal(tfg_run()))
if `nf' {
    display as error "test_fit_groups: `nf' failure(s)"
    exit 9
}
display as result "test_fit_groups: the pooled fit and inference reproduce lavaan's two-group fits of design E at all four levels (theta, T_raw via (N-G) F, a, b, T3, e(V)), the configural threshold block is MG6's oracle, the no-pi traces are told apart, and an all-groups reflection is a symmetry that canonicalizes back while a one-group sign change is not"
exit 0
