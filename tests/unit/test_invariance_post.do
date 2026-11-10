*! Module 11: what an invariance() fit posts (group surface, build step 4,
*! 2026-10-12). Executed under StataMP 18 on 2026-10-12 and it passes.
*!
*! Decisions under test (E-073, 2026-10-12): the multi-group
*! posting convention S4-1 (a) and its order S4-2 (a) (spec D5); what each fit
*! posts for the exemption (e(canon), e(canon_configural), e(invariance));
*! e(K) (option (ii)); O-T10 (a), the pooled e(bread_ratio); S4-3 (c)
*! narrowed, the groupwise Stage 1 notes (standardized: test_standardized_g.do); names
*! (spec 4.8, 2026-10-10: sem's group() form, gsem's thresholds, /scale(),
*! /mean()).
*!
*! Data fixtures/data/mg_E.csv (design E, groups g = 1, 2 of UNEQUAL size,
*! 1200 and 2800), F -> y1..y5, metric(unitvariance), the four levels.
*!
*! Checks, at every level:
*!   1  e(invariance) is the level; e(canon) is dwlsmv_levelgen's output for
*!      the configural text, byte for byte; e(canon_configural) is its
*!      configural output, and equals e(canon) at configural only.
*!   2  e(K) is 2 x 5, (4 4 4 2 2) in both rows, columns named y1..y5, rows
*!      named by group value (1, 2).
*!   3  membership (D5's contract carried to G > 1): 56 columns at every level
*!      (22 tau, 10 lambda, 2 psi, then 10 nu, 10 omega, 2 alpha); e(rank)
*!      the free count of the direct parameter map; e(Cns) one row per fixed
*!      column, selecting it alone at its value; e(V) exactly zero in every
*!      fixed row and column and positive on every free diagonal.
*!   4  names: unique; each parses to (ptype, group, row, col) and renders
*!      back to itself under the rule; the blocks in S4-2's order (tau,
*!      lambda, psi, then nu, omega, alpha group by group).
*!   5  end to end against lavaan 0.7.2 (fixtures/expected/lavaan_mg_E_*.csv,
*!      levels configural, thresholds, loadings = metric, intercepts =
*!      scalar): every free coefficient of e(b), found by its parsed key,
*!      after per-group sign equivalence on y1's loading, at the theta tier
*!      (1e-5 abs + 1e-4 rel); its SE from e(V) against lavaan's vcov at the
*!      SE tier (1e-8 abs + 1e-5 rel); e(df_T3) equals lavaan's df exactly;
*!      every lavaan row of the level is matched.
*!   6  e(bread_ratio) against a direct SVD of W^1/2 Delta, Delta the direct
*!      pooled fit's, W = blockdiag(pi_g / diag(Gamma_g)): 1e-10 relative
*!      (chosen before the first run: two SVD routines on one matrix); and the
*!      ratio without pi_g differs from it by more than that (the check can
*!      tell the weighting).
*! And once:
*!   7  replay shows the invariance header; standardized runs (rc 0)
*!      at estimation and on replay (the refusal was lifted 2026-11-01).
*!   8  groupwise notes (S4-3): keepbound with y6 = y1 in group 2 only prints
*!      the bound pair under "Group 2 (g = 2):" and nothing for group 1; an
*!      empty cell of y4 x y5 in group 2 only prints the repair under group 2.
*!   9  build/dwlsmv_example_mg.dta (tools/build_example.do, the help's
*!      invariance example, P5-0 (b)) holds mg_E.csv's values exactly.
*!
*! Mutants (2026-10-12, scratch worktree; M-a and M-d first failed to apply
*! and were re-run on code checked to be mutated): M-a pi_g dropped from the
*! bread weights, red at 6; M-b e(V) divided by N - 1, red at 5 (every SE);
*! M-c the omega block before nu, red at 4; M-d group 1's absent defaults not
*! posted, red (D6's check refuses, rc 3498); M-e e(canon_configural) set to
*! e(canon), red at 1 past configural.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
quietly do tests/dwlsmv_assert.mata
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

import delimited using "fixtures/expected/lavaan_mg_E_fits.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: FITS = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_E_parameters.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: PARS = st_sdata(., .)
import delimited using "fixtures/expected/lavaan_mg_E_vcov.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: VCOV = st_sdata(., .)
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
mata: YALL = st_data(., .)

local nfail 0
local m5 (F -> y1 y2 y3 y4 y5)

capture program drop tip_run
program tip_run
    args cmd keep
    capture mata: dwlsmv_m11_clear()
    if "`keep'" == "" ereturn clear
    tempfile lg
    capture log close tip
    quietly log using "`lg'", text replace name(tip)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tip
    mata: TIPLOG = invtokens(cat(st_local("lg"))')
    c_local tip_rc `rc'
end

capture program drop tip_has
program tip_has
    args s
    mata: st_local("h", strofreal(strpos(TIPLOG, st_local("s")) > 0))
    c_local hit `h'
end

mata:

string scalar tip_conf()
{
    string scalar s, nl
    real scalar g, j

    nl = char(10)
    s  = "groups 2" + nl
    for (g = 1; g <= 2; g++) {
        s = s + "group " + strofreal(g) + nl + "factor F" + nl + "ind y1 y2 y3 y4 y5" + nl
        for (j = 1; j <= 5; j++) {
            s = s + "lambda F y" + strofreal(j) + " free" + nl
        }
        s = s + (g == 1 ? "metric unitvariance F" : "psi F F fixed 1") + nl
    }
    return(s)
}

struct dwlsmv_fit scalar tip_dfit(string scalar level)
{
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    real colvector sel, KE
    real scalar g

    KE   = (4 \ 4 \ 4 \ 2 \ 2)
    sts  = dwlsmv_stage1(2, 1)
    maps = dwlsmv_statmap(2, 1)
    for (g = 1; g <= 2; g++) {
        sel     = selectindex(YALL[., 1] :== g)
        sts[g]  = dwls_sample_stats(YALL[sel, 2..6], KE)
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(tip_conf(), J(2, 1, KE'), level), maps, sts)
    fit = dwls_fit(sts, par)
    dwls_inference(fit)
    return(fit)
}

// (ptype, group, row, col) of a posted name, or a row of missing values if
// it follows no rule
real rowvector tip_key(string scalar eq, string scalar nm)
{
    string rowvector ind
    string scalar tg, v, body
    real scalar g, j, k, h

    ind = ("y1", "y2", "y3", "y4", "y5")
    if (eq != "/" & substr(eq, 1, 1) == "/") {
        h = strpos(nm, "#c.cut")
        if (!h) return((., ., ., .))
        tg = substr(nm, 1, h - 1)
        k  = strtoreal(substr(nm, h + 6, .))
        j  = 0
        for (g = 1; g <= 5; g++) if (ind[g] == substr(eq, 2, .)) j = g
        g  = tip_grp(tg)
        return((4, g, j, k))
    }
    if (eq == "/") {
        h    = strpos(nm, ")#")
        if (!h) return((., ., ., .))
        body = substr(nm, 1, h)
        g    = tip_grp(substr(nm, h + 2, .))
        if (body == "var(F)")   return((2, g, 1, 1))
        if (body == "mean(F)")  return((7, g, 1, 1))
        for (j = 1; j <= 5; j++) {
            if (body == "scale(" + ind[j] + ")") return((6, g, j, 1))
        }
        return((., ., ., .))
    }
    j = 0
    for (g = 1; g <= 5; g++) if (ind[g] == eq) j = g
    h = strpos(nm, "#c.F")
    if (h) return((1, tip_grp(substr(nm, 1, h - 1)), j, 1))
    return((5, tip_grp(nm), j, 1))
}

// a group tag 1bn.g / 2.g to its position (values 1, 2 here), else missing
real scalar tip_grp(string scalar tg)
{
    if (tg == "1bn.g") return(1)
    if (tg == "2.g")   return(2)
    return(.)
}

// the name the rule gives a key
string scalar tip_name(real rowvector k)
{
    string scalar tg, y

    tg = (k[2] == 1 ? "1bn.g" : "2.g")
    y  = "y" + strofreal(k[3])
    if (k[1] == 4) return("/" + y + ":" + tg + "#c.cut" + strofreal(k[4]))
    if (k[1] == 1) return(y + ":" + tg + "#c.F")
    if (k[1] == 2) return("/var(F)#" + tg)
    if (k[1] == 5) return(y + ":" + tg)
    if (k[1] == 6) return("/scale(" + y + ")#" + tg)
    return("/mean(F)#" + tg)
}

// checks 1-6 at one level; returns the number of failures
real scalar tip_level(string scalar level, string scalar lav)
{
    external string matrix PARS, VCOV, FITS
    struct dwlsmv_fit scalar fit
    string matrix st
    string colvector nm
    real rowvector b, key, sg, prev
    real matrix V, C, KK, Kp, M
    real colvector fr, gd, gd0, sv, sel
    real scalar bad, n, i, j, k, r, lv, lse, ours, ose, blk, nfix, N, rat, rat0, found, nl
    transmorphic A
    string colvector ks
    real colvector sm
    real scalar lik, nmiss, noff

    bad = 0
    KK  = J(2, 1, (4, 4, 4, 2, 2))
    // ---- 1 --------------------------------------------------------------
    if (st_global("e(invariance)") != level) {
        printf("{err}  1 %s: e(invariance) is %s\n", level, st_global("e(invariance)"))
        bad++
    }
    if (st_global("e(canon)") != dwlsmv_levelgen(tip_conf(), KK, level)) {
        printf("{err}  1 %s: e(canon) is not the generator's output\n", level)
        bad++
    }
    if (st_global("e(canon_configural)") != dwlsmv_levelgen(tip_conf(), KK, "configural")) {
        printf("{err}  1 %s: e(canon_configural) is not the generator's configural output\n", level)
        bad++
    }
    if ((st_global("e(canon)") == st_global("e(canon_configural)")) != (level == "configural")) {
        printf("{err}  1 %s: e(canon) == e(canon_configural) at the wrong level\n", level)
        bad++
    }
    // ---- 2 --------------------------------------------------------------
    Kp = st_matrix("e(K)")
    if (Kp != KK | st_matrixcolstripe("e(K)")[., 2]' != ("y1", "y2", "y3", "y4", "y5") |
        st_matrixrowstripe("e(K)")[., 2]' != ("1", "2")) {
        printf("{err}  2 %s: e(K) or its names\n", level)
        bad++
    }
    // ---- 3 --------------------------------------------------------------
    fit = tip_dfit(level)
    b   = st_matrix("e(b)")
    V   = st_matrix("e(V)")
    n   = cols(b)
    fr  = (diagonal(V) :!= 0)
    nfix = n - sum(fr)
    if (n != 56 | st_numscalar("e(rank)") != fit.par.npar) {
        printf("{err}  3 %s: %f columns (56 wanted), e(rank) %f (map %f)\n", level, n,
               st_numscalar("e(rank)"), fit.par.npar)
        bad++
    }
    C = st_matrix("e(Cns)")
    if (rows(C) != nfix) {
        printf("{err}  3 %s: e(Cns) has %f rows, %f fixed columns\n", level, rows(C), nfix)
        bad++
    }
    else {
        for (r = 1; r <= nfix; r++) {
            sel = selectindex(C[r, 1..n] :!= 0)
            if (cols(sel) != 1) {
                printf("{err}  3 %s: e(Cns) row %f selects %f columns\n", level, r, cols(sel))
                bad++
                continue
            }
            if (fr[sel[1]] | C[r, sel[1]] != 1 | C[r, n + 1] != b[sel[1]]) {
                printf("{err}  3 %s: e(Cns) row %f\n", level, r)
                bad++
            }
        }
    }
    for (i = 1; i <= n; i++) {
        if (!fr[i]) {
            if (any(V[i, .] :!= 0) | any(V[., i] :!= 0)) {
                printf("{err}  3 %s: e(V) not zero in fixed column %f\n", level, i)
                bad++
            }
        }
        else if (V[i, i] <= 0) {
            printf("{err}  3 %s: e(V) diagonal %f not positive\n", level, i)
            bad++
        }
    }
    // ---- 4 --------------------------------------------------------------
    st = st_matrixcolstripe("e(b)")
    nm = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        nm[i] = (st[i, 1] == "/" ? "/" + st[i, 2] : st[i, 1] + ":" + st[i, 2])
    }
    if (rows(uniqrows(nm)) != n) {
        printf("{err}  4 %s: names not unique\n", level)
        bad++
    }
    blk  = 1
    prev = (0, 0)
    for (i = 1; i <= n; i++) {
        key = tip_key(st[i, 1], st[i, 2])
        if (missing(key)) {
            printf("{err}  4 %s: name %s follows no rule\n", level, nm[i])
            bad++
            continue
        }
        if (tip_name(key) != nm[i]) {
            printf("{err}  4 %s: %s renders back as %s\n", level, nm[i], tip_name(key))
            bad++
        }
        // fixed defaults (ledger MG-G13): omega 1, alpha 0
        if (!fr[i] & ((key[1] == 6 & b[i] != 1) | (key[1] == 7 & b[i] != 0))) {
            printf("{err}  4 %s: fixed %s is %f\n", level, nm[i], b[i])
            bad++
        }
        // block order: tau 4, lambda 1, psi 2, nu 5, omega 6, alpha 7; within
        // nu, omega, alpha the group is nondecreasing
        // explicit ifs: Mata's ?: groups left to right, so the chained form
        // this replaced (2026-10-31) gave 3 for tau
        if (key[1] == 4)      k = 1
        else if (key[1] == 1) k = 2
        else if (key[1] == 2) k = 3
        else                  k = key[1] - 1
        if (k < prev[1] | (k == prev[1] & k >= 4 & key[2] < prev[2])) {
            printf("{err}  4 %s: %s out of S4-2's order\n", level, nm[i])
            bad++
        }
        prev = (k, key[2])
    }
    // ---- 5 --------------------------------------------------------------
    sg = J(1, 2, 1)
    for (j = 1; j <= 2; j++) {
        lv = tip_lav(PARS, lav, tip_kstr((1, j, 1, 1)))
        ours = b[tip_col(st, (1, j, 1, 1))]
        sg[j] = sign(ours) * sign(lv)
    }
    found = 0
    for (i = 1; i <= n; i++) {
        if (!fr[i]) continue
        key  = tip_key(st[i, 1], st[i, 2])
        lv   = tip_lav(PARS, lav, tip_kstr(key))
        if (missing(lv)) {
            printf("{err}  5 %s: free %s has no lavaan row\n", level, nm[i])
            bad++
            continue
        }
        found++
        ours = b[i] * ((key[1] == 1 | key[1] == 7) ? sg[key[2]] : 1)
        if (dwlsmv_gapratio(ours, lv, 1e-5, 1e-4) > 1) {
            printf("{err}  5 %s: %s %15.10f vs lavaan %15.10f\n", level, nm[i], ours, lv)
            bad++
        }
        lse = sqrt(tip_v(VCOV, lav, tip_kstr(key)))
        ose = sqrt(V[i, i])
        if (dwlsmv_gapratio(ose, lse, 1e-8, 1e-5) > 1) {
            printf("{err}  5 %s: SE %s %15.10f vs lavaan %15.10f\n", level, nm[i], ose, lse)
            bad++
        }
    }
    // every free off-diagonal of e(V) against lavaan's vcov (ledger MG-G20,
    // 2026-10-31), signs as above: |ours - lavaan| <= 2e-5 sqrt(v_ii v_kk),
    // the SE tier's 1e-5 relative carried to a covariance (chosen before the
    // first run)
    A = asarray_create()
    for (r = 1; r <= rows(VCOV); r++) {
        if (VCOV[r, 1] == lav) asarray(A, VCOV[r, 2] + "|" + VCOV[r, 3], strtoreal(VCOV[r, 4]))
    }
    ks = J(n, 1, "")
    sm = J(n, 1, 1)
    for (i = 1; i <= n; i++) {
        key = tip_key(st[i, 1], st[i, 2])
        if (missing(key)) continue
        ks[i] = tip_kstr(key)
        sm[i] = ((key[1] == 1 | key[1] == 7) ? sg[key[2]] : 1)
    }
    nmiss = 0
    noff  = 0
    for (i = 1; i <= n; i++) {
        if (!fr[i] | ks[i] == "") continue
        for (k = i + 1; k <= n; k++) {
            if (!fr[k] | ks[k] == "") continue
            // explicit ifs: a nested ?: read a stored 0 as false here
            lik = .
            if (asarray_contains(A, ks[i] + "|" + ks[k])) {
                lik = asarray(A, ks[i] + "|" + ks[k])
            }
            else if (asarray_contains(A, ks[k] + "|" + ks[i])) {
                lik = asarray(A, ks[k] + "|" + ks[i])
            }
            if (missing(lik)) {
                nmiss++
                continue
            }
            noff++
            ours = V[i, k] * sm[i] * sm[k]
            if (abs(ours - lik) > 2e-5 * sqrt(V[i, i] * V[k, k])) {
                printf("{err}  5 %s: cov(%s, %s) %12.6e vs lavaan %12.6e\n", level, nm[i], nm[k], ours, lik)
                bad++
            }
        }
    }
    if (nmiss | !noff) {
        printf("{err}  5 %s: %f free pairs without a lavaan vcov entry, %f compared\n", level, nmiss, noff)
        bad++
    }
    nl = 0
    for (i = 1; i <= rows(PARS); i++) nl = nl + (PARS[i, 1] == lav)
    if (found != nl) {
        printf("{err}  5 %s: %f free columns matched, lavaan has %f rows\n", level, found, nl)
        bad++
    }
    for (i = 1; i <= rows(FITS); i++) {
        if (FITS[i, 1] == lav) {
            if (strtoreal(FITS[i, 4]) != st_numscalar("e(df_T3)")) {
                printf("{err}  5 %s: e(df_T3) %f vs lavaan df %s\n", level,
                       st_numscalar("e(df_T3)"), FITS[i, 4])
                bad++
            }
        }
    }
    // ---- 6 --------------------------------------------------------------
    N   = fit.stage1[1].N + fit.stage1[2].N
    gd  = J(0, 1, .)
    gd0 = J(0, 1, .)
    for (j = 1; j <= 2; j++) {
        gd  = gd  \ ((fit.stage1[j].N / N) :/ diagonal(fit.stage1[j].Gamma))
        gd0 = gd0 \ (1 :/ diagonal(fit.stage1[j].Gamma))
    }
    sv   = svdsv(sqrt(gd) :* fit.Delta)
    rat  = min(sv) / max(sv)
    sv   = svdsv(sqrt(gd0) :* fit.Delta)
    rat0 = min(sv) / max(sv)
    printf("{txt}  6 %-10s e(bread_ratio) %12.9g direct %12.9g (rel %8.2e); without pi %12.9g\n",
           level, st_numscalar("e(bread_ratio)"), rat,
           abs(st_numscalar("e(bread_ratio)") - rat) / rat, rat0)
    if (abs(st_numscalar("e(bread_ratio)") - rat) > 1e-10 * rat) {
        printf("{err}  6 %s: e(bread_ratio) differs from the direct SVD\n", level)
        bad++
    }
    if (abs(rat0 - rat) <= 1e-10 * rat) {
        printf("{err}  6 %s: the ratio without pi cannot be told apart\n", level)
        bad++
    }
    return(bad)
}

string scalar tip_kstr(real rowvector k)
{
    return(strofreal(k[1]) + ":" + strofreal(k[2]) + ":" + strofreal(k[3]) + ":" + strofreal(k[4]))
}

real scalar tip_col(string matrix st, real rowvector k)
{
    real scalar i

    for (i = 1; i <= rows(st); i++) {
        if (tip_key(st[i, 1], st[i, 2]) == k) return(i)
    }
    return(.)
}

real scalar tip_lav(string matrix M, string scalar lev, string scalar key)
{
    real scalar i

    for (i = 1; i <= rows(M); i++) {
        if (M[i, 1] == lev & M[i, 2] == key) return(strtoreal(M[i, 3]))
    }
    return(.)
}

real scalar tip_v(string matrix M, string scalar lev, string scalar key)
{
    real scalar i

    for (i = 1; i <= rows(M); i++) {
        if (M[i, 1] == lev & M[i, 2] == key & M[i, 3] == key) return(strtoreal(M[i, 4]))
    }
    return(.)
}

end

local levs configural thresholds metric scalar
local lavs configural thresholds loadings intercepts
forvalues i = 1/4 {
    local lev : word `i' of `levs'
    local lav : word `i' of `lavs'
    tip_run `"dwlsmv `m5', group(g) invariance(`lev') metric(unitvariance)"'
    if `tip_rc' {
        display as error "  `lev' FAILED: rc `tip_rc'"
        local ++nfail
        continue
    }
    mata: st_local("b", strofreal(tip_level("`lev'", "`lav'")))
    if `b' {
        display as error "  checks 1-6 at `lev': `b' failure(s)"
        local nfail = `nfail' + `b'
    }
    else display as text "  checks 1-6 at `lev': ok"
}

display as text "7. replay and standardized"
quietly dwlsmv `m5', group(g) invariance(scalar) metric(unitvariance)
tip_run "dwlsmv" keep
tip_has "Invariance: scalar, 2 groups; factor variances fixed at 1 in the first group (the reference group)"
local h1 = `hit' & `tip_rc' == 0
quietly dwlsmv `m5', group(g) invariance(scalar) metric(unitvariance)
capture noisily dwlsmv, standardized
local r2 = _rc
tip_run `"dwlsmv `m5', group(g) invariance(scalar) metric(unitvariance) standardized"'
if `h1' & `r2' == 0 & `tip_rc' == 0 display as text "  7: ok"
else {
    display as error "  7 FAILED (replay `h1', replay standardized rc `r2', estimation rc `tip_rc')"
    local ++nfail
}

display as text "8. groupwise notes"
preserve
quietly generate byte y6 = cond(g == 2, y1, mod(_n, 4))
tip_run `"dwlsmv (F -> y1 y2 y3 y6), group(g) invariance(configural) metric(unitvariance) keepbound"'
tip_has "Group 2 (g = 2): Pairs included at the bound under keepbound:     y1~~y6"
local h1 `hit'
tip_has "Group 1 (g = 1):"
local h2 `hit'
if `tip_rc' == 0 & `h1' & !`h2' display as text "  8a keepbound pair, group 2 only: ok"
else {
    display as error "  8a FAILED (rc `tip_rc', group 2 pair `h1', group 1 note `h2')"
    local ++nfail
}
restore
preserve
quietly replace y5 = 1 if y4 == 1 & g == 2
tip_run `"dwlsmv `m5', group(g) invariance(configural) metric(unitvariance)"'
tip_has "Group 2 (g = 2): One-empty-cell repair applied to:     y4~~y5"
local h1 `hit'
tip_has "Group 1 (g = 1):"
local h2 `hit'
if `tip_rc' == 0 & `h1' & !`h2' display as text "  8b empty-cell repair, group 2 only: ok"
else {
    display as error "  8b FAILED (rc `tip_rc', group 2 repair `h1', group 1 note `h2')"
    local ++nfail
}
restore

display as text "9. the invariance example dataset is the fixture's data"
* P5-0 (b), 2026-10-12: build/dwlsmv_example_mg.dta (tools/build_example.do)
* holds mg_E.csv's values exactly, g renamed grp
quietly do tools/build_example.do
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
tempfile fx
quietly save `fx'
quietly use "build/dwlsmv_example_mg.dta", clear
rename grp g
capture noisily cf g y1 y2 y3 y4 y5 using `fx'
if _rc == 0 & _N == 4000 display as text "  9: ok"
else {
    display as error "  9 FAILED: the example data differ from mg_E.csv (cf rc " _rc ")"
    local ++nfail
}

display as text "10. group values 0 and 7, rows shuffled; missing group values"
* ledger MG-G14: e(K) rows and e(b) names carry the group values, not their
* positions; the estimates do not depend on the values or the row order
* (1e-7 relative, chosen before the first run: one fit, summation order only)
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv `m5', group(g) invariance(metric) metric(unitvariance)
matrix tip_B12 = e(b)
quietly recode g (1 = 0) (2 = 7)
set seed 20261031
quietly generate double tip_u = runiform()
sort tip_u
drop tip_u
tip_run `"dwlsmv `m5', group(g) invariance(metric) metric(unitvariance)"'
local rn : rownames e(K)
local cn : colfullnames e(b)
mata: st_local("d", strofreal(mreldif(st_matrix("e(b)"), st_matrix("tip_B12"))))
if `tip_rc' == 0 & "`rn'" == "0 7" & e(N) == 4000 & strpos("`cn'", "0bn.g") & ///
    strpos("`cn'", "7.g") & !strpos("`cn'", "1bn.g") & `d' <= 1e-7 {
    display as text "  10a group values 0, 7: ok (e(b) reldif `d')"
}
else {
    display as error "  10a FAILED (rc `tip_rc', e(K) rows '`rn'', e(N) " e(N) ", e(b) reldif `d')"
    local ++nfail
}
* ledger MG-G18: rows with a missing group value are outside the estimation
* sample, and e(N) agrees with e(sample)
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly replace g = . in 1/100
tip_run `"dwlsmv `m5', group(g) invariance(configural) metric(unitvariance)"'
quietly count if e(sample)
if `tip_rc' == 0 & e(N) == 3900 & r(N) == 3900 {
    display as text "  10b missing group values: ok"
}
else {
    display as error "  10b FAILED (rc `tip_rc', e(N) " e(N) ", e(sample) " r(N) ")"
    local ++nfail
}

capture mata: dwlsmv_m11_clear()
adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_invariance_post: `nfail' failure(s)"
    exit 9
}
display as result "test_invariance_post: an invariance() fit posts the full group-expanded table at every level, named and ordered by rule, with e(K), e(canon*), e(invariance) and the pooled bread ratio, and it reproduces lavaan's estimates and SEs end to end"
exit 0
