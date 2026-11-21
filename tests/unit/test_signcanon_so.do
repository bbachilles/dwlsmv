*! Executed under StataMP 18 on 2026-11-20: all checks pass on the first run. Drafted by a read-only
*! agent, reviewed and run by the main session. Each mutant below was run in a scratch worktree and
*! killed by the check named (M5 also by test_translate_second_order.do check 8).
*!
*! Module 7, release B (decision B3.10): dwlsmv_signcanon() reflects a higher-order factor G
*! (no free loading, at least one free outgoing path) when the sum of its children's
*! STANDARDIZED paths, taken after the first-order factors are anchored, is below 0.
*!
*! Checks (each assertion-based; the mutant each kills is named):
*!   a  G over 3, U: all 16 sign images (F1..F3 and G, the D-transform) return the canonical
*!      theta bit for bit.                                       kills M5 (the G-only images)
*!   b  standardized, not raw: paths (3,-1,-1): raw sum 1 > 0, standardized sum < 0, G flips.
*!                                                               kills M1
*!   c  ties: paths (1,-1,0) and (-1,1,0) come back unchanged.   kills M3
*!   d  order with pass 1: F1 reflected in pass 1 changes the sign of the score; G not flipped.
*!                                                               kills M2
*!   e  metric V (item-marked F's, psi_k free): paths (-.8,-.6,-.5) flip; F loadings untouched.
*!                                                               kills M5
*!   f  MD under U: G with a free item loading is left to pass 1 (f1: positive loading, negative
*!      paths, unchanged; f2: negative loading, positive paths, reflected by pass 1 only).
*!                                                               kills M6
*!   g  phi guard: a child with Phi[k,k] < 0 is excluded, the others decide.   kills M4
*!   h  two correlated G's: G1 flips, G2 does not, psi_G1G2 changes sign.     kills M5
*!   i  third level (hand-built table, G2 has the LOWER factor index than G1): G1 must be
*!      decided before G2.                                       kills M7
*!   j  sigma(theta) before and after, 1e-14, inside a, b, e, h (dwls_model, one group).
*!   k  several groups (hand-built table, grp has a 2): pass 2 does not run; the same table
*!      with grp all 1 is reflected (control).                   kills M8
*!
*! Mutants: M1 raw b in the score; M2 s_k := 1 (pass-1 flips ignored); M3 score < 0 -> <= 0;
*! M4 phi guard removed; M5 pass 2 removed; M6 pass 2 also applied to a factor with a free
*! loading; M7 G's decided in factor index order without waiting for G-children; M8 group
*! restriction removed.
*!
*! Exit 0 pass, 9 fail.
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
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    quietly do "`f'"
}

mata:

string scalar tsc_join(string colvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= rows(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

string scalar tsc_gname(real scalar g, real scalar nG)
{
    if (nG == 1) {
        return("G")
    }
    return("G" + strofreal(g))
}

// Same builder as tso_lines() in test_identify_second_order.do: first-order F1..FK with
// nper[k] items each, metric met ("U" or "V"), parent[k] the G that F_k loads on (0 none),
// nG higher-order factors (declared AFTER F1..FK, so F_k has index k, G_g index K + g),
// each 'metric unitvariance' when gmet; xitems extra items for the caller's own lines.
string colvector tsc_lines(real rowvector nper, real rowvector parent,
                           real scalar nG, string scalar met,
                           real scalar gmet, real scalar xitems)
{
    string colvector L
    string scalar fl, il, f
    real scalar k, j, y, g, ntot

    fl = "factor"
    for (k = 1; k <= cols(nper); k++) {
        fl = fl + " F" + strofreal(k)
    }
    for (g = 1; g <= nG; g++) {
        fl = fl + " " + tsc_gname(g, nG)
    }
    ntot = sum(nper) + xitems
    il = "ind"
    for (y = 1; y <= ntot; y++) {
        il = il + " y" + strofreal(y)
    }
    L = (fl \ il)

    y = 0
    for (k = 1; k <= cols(nper); k++) {
        f = "F" + strofreal(k)
        for (j = 1; j <= nper[k]; j++) {
            y++
            if (met == "V" & j == 1) {
                L = L \ ("metric marker " + f + " y" + strofreal(y))
            }
            else {
                L = L \ ("lambda " + f + " y" + strofreal(y) + " free")
            }
        }
        if (met == "V") {
            L = L \ ("psi " + f + " " + f + " free")
        }
        else {
            L = L \ ("metric unitvariance " + f)
        }
        if (parent[k] > 0) {
            L = L \ ("beta " + f + " " + tsc_gname(parent[k], nG) + " free")
        }
    }
    if (gmet) {
        for (g = 1; g <= nG; g++) {
            L = L \ ("metric unitvariance " + tsc_gname(g, nG))
        }
    }
    return(L)
}

struct dwlsmv_parmap scalar tsc_par(string colvector L, real scalar p)
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(p, 1, 3))
    par  = dwls_identify(tsc_join(L), m)
    return(par)
}

real colvector tsc_sig(real colvector th, struct dwlsmv_parmap scalar par, real scalar p)
{
    struct dwlsmv_statmap vector m

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(p, 1, 3))
    return(dwls_model(th, par, m))
}

// A theta from values by role. lamv indexed by item (the lambda row), pathv by the
// factor the beta runs OUT of (the beta row), psiv by factor (free diagonal psi, metric
// V), psig the one free off-diagonal psi. Thresholds get distinct fixed numbers.
real colvector tsc_theta(struct dwlsmv_parmap scalar par, real colvector lamv,
                         real colvector pathv, real colvector psiv, real scalar psig)
{
    real colvector th
    real scalar i, t, r, c

    th = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) {
            continue
        }
        t = par.ptype[i]
        r = par.row[i]
        c = par.col[i]
        if (t == 4) {
            th[par.fidx[i]] = -0.5 + 0.1 * c
        }
        else if (t == 1) {
            th[par.fidx[i]] = lamv[r]
        }
        else if (t == 3) {
            th[par.fidx[i]] = pathv[r]
        }
        else if (t == 2 & r == c) {
            th[par.fidx[i]] = psiv[r]
        }
        else if (t == 2) {
            th[par.fidx[i]] = psig
        }
        else {
            _error("tsc_theta: unexpected free parameter type")
        }
    }
    if (missing(th)) {
        _error("tsc_theta: a free parameter was not given a value")
    }
    return(th)
}

// The D-transform written from the derivation: negate lambda in a reversed column, and
// beta / off-diagonal psi with exactly one reversed end. fl is the 0/1 vector over factors.
real colvector tsc_refl(real colvector theta, struct dwlsmv_parmap scalar par,
                        real colvector fl)
{
    real colvector th, neg
    real scalar i, t, r, c

    th  = theta
    neg = J(par.npar, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) {
            continue
        }
        t = par.ptype[i]
        r = par.row[i]
        c = par.col[i]
        if (t == 1) {
            if (fl[c]) {
                neg[par.fidx[i]] = 1
            }
        }
        else if (t == 3 | (t == 2 & r != c)) {
            if (fl[r] != fl[c]) {
                neg[par.fidx[i]] = 1
            }
        }
    }
    for (i = 1; i <= par.npar; i++) {
        if (neg[i]) {
            th[i] = -th[i]
        }
    }
    return(th)
}

// Position in theta of the free parameter (t, r, c); psi matches either order.
real scalar tsc_pos(struct dwlsmv_parmap scalar par, real scalar t, real scalar r,
                    real scalar c)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.ptype[i] == t &
            ((par.row[i] == r & par.col[i] == c) |
             (t == 2 & par.row[i] == c & par.col[i] == r))) {
            return(par.fidx[i])
        }
    }
    _error("tsc_pos: no such free parameter")
}

real scalar tsc_chk(real scalar ok, string scalar msg)
{
    if (!ok) {
        printf("{err}  FAIL: %s\n", msg)
        return(1)
    }
    return(0)
}

real colvector tsc_bits(real scalar mask, real scalar n)
{
    real colvector b
    real scalar j

    b = J(n, 1, 0)
    for (j = 1; j <= n; j++) {
        b[j] = mod(floor(mask / 2^(j - 1)), 2)
    }
    return(b)
}

// ---- a (and j): G over 3, U, all 16 sign images -----------------------------------
// canonical: loadings positive, paths (0.8, 0.6, 0.5). Hand: Phi[G,G] = 1, Phi[k,k] =
// 1 + b_k^2, so the canonical score is > 0 and the G-reflected image's is < 0; after pass
// 1 the F's are positive in every image, so G alone is decided by the paths.
real scalar tsc_a()
{
    struct dwlsmv_parmap scalar par
    real colvector th0, img, out, fl, lamv, s0, s1
    real scalar mask, nfail, differs

    par  = tsc_par(tsc_lines((3,3,3), (1,1,1), 1, "U", 1, 0), 9)
    lamv = (0.5 \ 0.4 \ 0.6 \ 0.5 \ 0.4 \ 0.6 \ 0.5 \ 0.4 \ 0.6)
    th0  = tsc_theta(par, lamv, (0.8 \ 0.6 \ 0.5), J(4, 1, .), .)
    s0   = tsc_sig(th0, par, 9)
    nfail = 0
    for (mask = 0; mask < 16; mask++) {
        fl  = tsc_bits(mask, 4)
        img = tsc_refl(th0, par, fl)
        out = dwlsmv_signcanon(img, par)
        differs = (max(abs(img - th0)) > 0)
        nfail = nfail + tsc_chk(differs == (mask > 0), "a: an image should differ from theta0 exactly when mask > 0")
        nfail = nfail + tsc_chk(!missing(out), "a: signcanon returned a missing value")
        nfail = nfail + tsc_chk(max(abs(out - th0)) == 0, "a: a sign image did not canonicalize back to theta0 bit for bit")
        s1 = tsc_sig(out, par, 9)
        nfail = nfail + tsc_chk(max(abs(tsc_sig(img, par, 9) - s1)) <= 1e-14, "j(a): sigma changed under canonicalization")
        nfail = nfail + tsc_chk(max(abs(s1 - s0)) <= 1e-14, "j(a): sigma at the canonical theta differs from sigma at theta0")
    }
    printf("{txt}  a: 16 sign images of G over 3 (U) canonicalize back exactly; sigma equal to 1e-14\n")
    return(nfail)
}

// ---- b: standardized, not raw ------------------------------------------------------
// paths (3,-1,-1), loadings 0.3 (so sigma exists: 0.09 * Phi[F1,F1] = 0.9 < 1).
// Phi[k,k] = 10, 2, 2; Phi[G,G] = 1. b* = 3/sqrt(10) = 0.9487, -1/sqrt(2) = -0.7071 twice;
// sum = -0.4655 < 0 -> flip. Raw sum = 3 - 1 - 1 = 1 > 0 (M1 would not flip).
real scalar tsc_b()
{
    struct dwlsmv_parmap scalar par
    real colvector fit, want, out, lamv
    real scalar nfail

    par  = tsc_par(tsc_lines((3,3,3), (1,1,1), 1, "U", 1, 0), 9)
    lamv = J(9, 1, 0.3)
    fit  = tsc_theta(par, lamv, (3 \ -1 \ -1), J(4, 1, .), .)
    want = tsc_theta(par, lamv, (-3 \ 1 \ 1), J(4, 1, .), .)
    out  = dwlsmv_signcanon(fit, par)
    nfail = 0
    nfail = nfail + tsc_chk(sum((3, -1, -1)) > 0, "b: premise: raw sum should be positive")
    nfail = nfail + tsc_chk(3/sqrt(10) - 2/sqrt(2) < 0, "b: premise: standardized sum should be negative")
    nfail = nfail + tsc_chk(max(abs(out - want)) == 0, "b: paths (3,-1,-1) should become (-3,1,1), nothing else moved")
    nfail = nfail + tsc_chk(max(abs(tsc_sig(fit, par, 9) - tsc_sig(out, par, 9))) <= 1e-14, "j(b): sigma changed")
    printf("{txt}  b: paths (3,-1,-1): raw sum 1, standardized sum %7.4f, G reflected %g\n",
           3/sqrt(10) - 2/sqrt(2), max(abs(out - want)) == 0)
    return(nfail)
}

// ---- c: ties -----------------------------------------------------------------------
// Phi = diag-ish with Phi[1,1] = Phi[2,2] = 2: the two terms are s*b*sqrt(1/2) with b = +-1
// and the third is 0 * sqrt(1) = 0; I - B is unit triangular with entries in {0, 1, -1} so
// luinv is exact and the sum is exactly 0 (a + -a).
real scalar tsc_c()
{
    struct dwlsmv_parmap scalar par
    real colvector t1, t2, lamv
    real scalar nfail

    par  = tsc_par(tsc_lines((3,3,3), (1,1,1), 1, "U", 1, 0), 9)
    lamv = J(9, 1, 0.5)
    t1   = tsc_theta(par, lamv, (1 \ -1 \ 0), J(4, 1, .), .)
    t2   = tsc_theta(par, lamv, (-1 \ 1 \ 0), J(4, 1, .), .)
    nfail = 0
    nfail = nfail + tsc_chk(max(abs(dwlsmv_signcanon(t1, par) - t1)) == 0, "c: tie (1,-1,0) was changed")
    nfail = nfail + tsc_chk(max(abs(dwlsmv_signcanon(t2, par) - t2)) == 0, "c: tie (-1,1,0) was changed")
    printf("{txt}  c: ties (1,-1,0) and (-1,1,0) left as fitted\n")
    return(nfail)
}

// ---- d: order with pass 1 ----------------------------------------------------------
// canonical paths (0.8, 0.3, 0.3); image: F1 reflected, so F1's loadings and b_1 are negated:
// raw (-0.8, 0.3, 0.3). Phi[k,k] = 1.64, 1.09, 1.09. Standardized, ignoring pass 1:
//   -0.8/1.2806 + 2 * 0.3/1.0440 = -0.6247 + 0.5747 = -0.0500 < 0 (would flip G: wrong).
// After pass 1 (s_1 = -1): +0.6247 + 0.5747 = 1.1994 > 0, G stays. Expected: theta0.
real scalar tsc_d()
{
    struct dwlsmv_parmap scalar par
    real colvector th0, img, out, lamv
    real scalar nfail, naive

    par  = tsc_par(tsc_lines((3,3,3), (1,1,1), 1, "U", 1, 0), 9)
    lamv = J(9, 1, 0.5)
    th0  = tsc_theta(par, lamv, (0.8 \ 0.3 \ 0.3), J(4, 1, .), .)
    img  = tsc_refl(th0, par, (1 \ 0 \ 0 \ 0))
    out  = dwlsmv_signcanon(img, par)
    naive = -0.8/sqrt(1.64) + 2 * 0.3/sqrt(1.09)
    nfail = 0
    nfail = nfail + tsc_chk(naive < 0, "d: premise: the pass-1-blind score should be negative")
    nfail = nfail + tsc_chk(max(abs(img - th0)) > 0, "d: premise: the image should differ from theta0")
    nfail = nfail + tsc_chk(max(abs(out - th0)) == 0, "d: F1 reflected in pass 1 must give theta0 with G not flipped")
    printf("{txt}  d: pass-1-blind score %8.4f (< 0), after pass 1 positive; G not flipped %g\n",
           naive, max(abs(out - th0)) == 0)
    return(nfail)
}

// ---- e (and j): metric V -----------------------------------------------------------
// F's carry a fixed marker loading (not reflectable); psi_k = 0.3, psi_GG fixed 1.
// Phi[k,k] = b_k^2 + 0.3 = 0.94, 0.66, 0.55 (< 1, so the marker's residual is positive).
// paths (-0.8,-0.6,-0.5): every standardized term negative -> flip -> (0.8,0.6,0.5).
real scalar tsc_e()
{
    struct dwlsmv_parmap scalar par
    real colvector fit, want, out, lamv, psiv
    real scalar nfail, i

    par  = tsc_par(tsc_lines((3,3,3), (1,1,1), 1, "V", 1, 0), 9)
    lamv = (1 \ 0.8 \ 0.7 \ 1 \ 0.8 \ 0.7 \ 1 \ 0.8 \ 0.7)
    psiv = (0.3 \ 0.3 \ 0.3 \ 1)
    fit  = tsc_theta(par, lamv, (-0.8 \ -0.6 \ -0.5), psiv, .)
    want = tsc_theta(par, lamv, (0.8 \ 0.6 \ 0.5), psiv, .)
    out  = dwlsmv_signcanon(fit, par)
    nfail = 0
    nfail = nfail + tsc_chk(max(abs(out - want)) == 0, "e: V paths (-.8,-.6,-.5) should flip to positive, F loadings untouched")
    nfail = nfail + tsc_chk(max(abs(tsc_sig(fit, par, 9) - tsc_sig(out, par, 9))) <= 1e-14, "j(e): sigma changed")
    printf("{txt}  e: metric V, paths (-.8,-.6,-.5) reflected, loadings untouched: %g\n",
           max(abs(out - want)) == 0)
    return(nfail)
}

// ---- f: MD under U -----------------------------------------------------------------
// G over 4 plus y13 on G (the accepted a6 table of test_identify_second_order). G has a free
// loading, so pass 1 owns it and pass 2 must skip it.
//   f1: loading +0.5, paths negative: unchanged (a pass-2 on G would reflect it).
//   f2: loading -0.5, paths positive: pass 1 reflects G: loading and paths negated, once.
real scalar tsc_f()
{
    struct dwlsmv_parmap scalar par
    real colvector fit, want, out, lamv, pv
    real scalar nfail

    par  = tsc_par(tsc_lines((3,3,3,3), (1,1,1,1), 1, "U", 1, 1) \ "lambda G y13 free", 13)
    nfail = 0

    lamv = J(13, 1, 0.5)
    fit  = tsc_theta(par, lamv, (-0.8 \ -0.6 \ -0.5 \ -0.4), J(5, 1, .), .)
    out  = dwlsmv_signcanon(fit, par)
    nfail = nfail + tsc_chk(max(abs(out - fit)) == 0, "f1: MD, G loading positive: pass 2 must not touch G")

    lamv = J(13, 1, 0.5)
    lamv[13] = -0.5
    pv   = (0.8 \ 0.6 \ 0.5 \ 0.4)
    fit  = tsc_theta(par, lamv, pv, J(5, 1, .), .)
    lamv[13] = 0.5
    want = tsc_theta(par, lamv, -pv, J(5, 1, .), .)
    out  = dwlsmv_signcanon(fit, par)
    nfail = nfail + tsc_chk(max(abs(out - want)) == 0, "f2: MD, G loading negative: reflected by pass 1 exactly once")
    printf("{txt}  f: MD: pass 1 only (positive loading kept; negative loading reflected once)\n")
    return(nfail)
}

// ---- g: the phi guard --------------------------------------------------------------
// metric V, psi_k = (-10, 0.5, 0.5), paths (3, -0.6, -0.5). Phi[G,G] = 1;
//   Phi[1,1] = 3^2 + (-10) = -1  (not > 0: excluded)
//   Phi[2,2] = 0.36 + 0.5 = 0.86 ; Phi[3,3] = 0.25 + 0.5 = 0.75
// guarded score = -0.6/sqrt(0.86) - 0.5/sqrt(0.75) = -0.6470 - 0.5774 = -1.2244 < 0 -> flip.
// Unguarded: sqrt(1/-1) is missing, the score is missing, missing < 0 is false -> no flip.
real scalar tsc_g()
{
    struct dwlsmv_parmap scalar par
    real colvector fit, want, out, lamv, psiv
    real scalar nfail

    par  = tsc_par(tsc_lines((3,3,3), (1,1,1), 1, "V", 1, 0), 9)
    lamv = (1 \ 0.8 \ 0.7 \ 1 \ 0.8 \ 0.7 \ 1 \ 0.8 \ 0.7)
    psiv = (-10 \ 0.5 \ 0.5 \ 1)
    fit  = tsc_theta(par, lamv, (3 \ -0.6 \ -0.5), psiv, .)
    want = tsc_theta(par, lamv, (-3 \ 0.6 \ 0.5), psiv, .)
    out  = dwlsmv_signcanon(fit, par)
    nfail = 0
    nfail = nfail + tsc_chk(max(abs(out - want)) == 0, "g: the child with Phi[k,k] < 0 must be excluded and G reflected")
    printf("{txt}  g: phi guard: Phi[1,1] = -1 excluded, score %8.4f -> G reflected %g\n",
           -0.6/sqrt(0.86) - 0.5/sqrt(0.75), max(abs(out - want)) == 0)
    return(nfail)
}

// ---- h (and j): two correlated G's -------------------------------------------------
// F1..F6 = 1..6, G1 = 7 (over F1-F3), G2 = 8 (over F4-F6), psi_G2G1 free = 0.3.
// All four reflections of (G1, G2) return the canonical theta; in the G1-only image G1's
// paths and psi are negated, and only G1 is turned back.
real scalar tsc_h()
{
    struct dwlsmv_parmap scalar par
    real colvector th0, img, out, fl, lamv, pv
    real scalar nfail, mask, pg, p1

    par  = tsc_par(tsc_lines((3,3,3,3,3,3), (1,1,1,2,2,2), 2, "U", 1, 0) \ "psi G2 G1 free", 18)
    lamv = (0.5 \ 0.4 \ 0.6)
    lamv = J(6, 1, 1) # lamv
    pv   = (0.8 \ 0.6 \ 0.5 \ 0.7 \ 0.6 \ 0.5)
    th0  = tsc_theta(par, lamv, pv, J(8, 1, .), 0.3)
    pg   = tsc_pos(par, 2, 8, 7)
    nfail = 0
    for (mask = 0; mask < 4; mask++) {
        fl  = J(6, 1, 0) \ tsc_bits(mask, 2)
        img = tsc_refl(th0, par, fl)
        out = dwlsmv_signcanon(img, par)
        nfail = nfail + tsc_chk(max(abs(out - th0)) == 0, "h: a G-reflection image did not return the canonical theta")
        nfail = nfail + tsc_chk(max(abs(tsc_sig(img, par, 18) - tsc_sig(out, par, 18))) <= 1e-14, "j(h): sigma changed")
    }
    // G1 only: the explicit entries.
    img = tsc_refl(th0, par, (J(6, 1, 0) \ 1 \ 0))
    out = dwlsmv_signcanon(img, par)
    nfail = nfail + tsc_chk(img[pg] == -0.3 & out[pg] == 0.3, "h: psi_G1G2 should change sign when only G1 is turned")
    for (p1 = 1; p1 <= 6; p1++) {
        nfail = nfail + tsc_chk(out[tsc_pos(par, 3, p1, (p1 <= 3 ? 7 : 8))] == pv[p1], "h: a path is not at its canonical value")
    }
    printf("{txt}  h: two correlated G's: all four reflections canonicalize; psi_G1G2 %4.1f -> %4.1f\n",
           img[pg], out[pg])
    return(nfail)
}

// ---- i: third level, hand-built table ---------------------------------------------
// Whether dwls_identify() accepts a second-order-over-second-order table was not determined
// by the drafting agent (no run), so the table is built by hand: signcanon reads only
// ptype, row, col, free, value, fidx, nall, npar and grp. Factors F1..F4 = 1..4, G2 = 5,
// G1 = 6 (the HIGHER G has the LOWER index, so index order decides it first, wrongly).
// Betas: F1,F2 <- G1; G1,F3,F4 <- G2. One free loading (0.5) per F; psi diagonals fixed 1.
// Image: b(F1,G1) = -0.8, b(F2,G1) = -0.6, b(G1,G2) = -0.7, b(F3,G2) = b(F4,G2) = 0.2.
// Phi[G2,G2] = 1, Phi[G1,G1] = 1 + 0.49 = 1.49, Phi[F1,F1] = 1.64, Phi[F2,F2] = 1.36,
// Phi[F3,F3] = Phi[F4,F4] = 1.04.
// Right order: G1 first: -0.8*sqrt(1.49/1.64) - 0.6*sqrt(1.49/1.36) = -0.7626 - 0.6280 =
//   -1.3906 < 0 -> flip G1. Then G2 with s_G1 = -1: (-1)(-0.7)sqrt(1/1.49) + 2*0.2*sqrt(1/1.04)
//   ... = 0.5735 + 0.3922 = 0.9657 > 0 -> G2 stays. Result: the canonical (all positive).
// Index order: G2 first with s_G1 = 1: -0.5735 + 0.3922 = -0.1813 < 0 -> flips G2 (wrong),
//   then G1 flips: a different theta.
real matrix tsc_Si(real colvector b)
{
    return((1,1,1,1,0.5 \ 1,2,2,1,0.5 \ 1,3,3,1,0.5 \ 1,4,4,1,0.5 \
            3,1,6,1,b[1] \ 3,2,6,1,b[2] \ 3,6,5,1,b[3] \ 3,3,5,1,b[4] \ 3,4,5,1,b[5] \
            2,1,1,0,1 \ 2,2,2,0,1 \ 2,3,3,0,1 \ 2,4,4,0,1 \ 2,5,5,0,1 \ 2,6,6,0,1))
}

struct dwlsmv_parmap scalar tsc_hand(real matrix S, real colvector grp)
{
    struct dwlsmv_parmap scalar par
    real scalar i, n, k

    n         = rows(S)
    par.nall  = n
    par.name  = J(n, 1, "x")
    par.ptype = S[., 1]
    par.row   = S[., 2]
    par.col   = S[., 3]
    par.free  = S[., 4]
    par.value = S[., 5]
    par.fidx  = J(n, 1, 0)
    k = 0
    for (i = 1; i <= n; i++) {
        if (S[i, 4]) {
            k++
            par.fidx[i] = k
        }
    }
    par.npar = k
    par.grp  = grp
    return(par)
}

real colvector tsc_hand_theta(struct dwlsmv_parmap scalar par)
{
    real colvector th
    real scalar i

    th = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            th[par.fidx[i]] = par.value[i]
        }
    }
    return(th)
}

real scalar tsc_i()
{
    struct dwlsmv_parmap scalar pi_, pc
    real colvector img, want, out
    real scalar nfail, n

    n     = rows(tsc_Si((0 \ 0 \ 0 \ 0 \ 0)))
    pi_   = tsc_hand(tsc_Si((-0.8 \ -0.6 \ -0.7 \ 0.2 \ 0.2)), J(n, 1, 1))
    pc    = tsc_hand(tsc_Si((0.8 \ 0.6 \ 0.7 \ 0.2 \ 0.2)), J(n, 1, 1))
    img   = tsc_hand_theta(pi_)
    want  = tsc_hand_theta(pc)
    out   = dwlsmv_signcanon(img, pi_)
    nfail = 0
    nfail = nfail + tsc_chk(dwlsmv_nfac(pi_) == 6, "i: the hand-built table should have 6 factors")
    nfail = nfail + tsc_chk(max(abs(out - want)) == 0, "i: G1 must be decided before G2; the result should be the all-positive theta")
    printf("{txt}  i: third level, G2 (index 5) over G1 (index 6): G1 decided first, canonical returned %g\n",
           max(abs(out - want)) == 0)
    return(nfail)
}

// ---- k: several groups --------------------------------------------------------------
// Module 3 refuses group() with a second-order factor (B3.5, test_identify_second_order r8),
// so a table with a G and two groups cannot come from dwls_identify(); this one is built
// by hand. G = 4 over F1..F3, paths (-0.8,-0.6,-0.5), one free loading (0.5) per F.
// grp all 1 (control): reflected. grp with a 2 on the last row: unchanged.
real scalar tsc_k()
{
    struct dwlsmv_parmap scalar p1, p2
    real matrix S
    real colvector grp, out, th, want
    real scalar nfail, n

    S = (1,1,1,1,0.5 \ 1,2,2,1,0.5 \ 1,3,3,1,0.5 \
         3,1,4,1,-0.8 \ 3,2,4,1,-0.6 \ 3,3,4,1,-0.5 \
         2,1,1,0,1 \ 2,2,2,0,1 \ 2,3,3,0,1 \ 2,4,4,0,1)
    n   = rows(S)
    p1  = tsc_hand(S, J(n, 1, 1))
    grp = J(n, 1, 1)
    grp[n] = 2
    p2  = tsc_hand(S, grp)
    th  = tsc_hand_theta(p1)
    want = th
    want[4::6] = -th[4::6]
    nfail = 0
    nfail = nfail + tsc_chk(max(abs(dwlsmv_signcanon(th, p1) - want)) == 0, "k: control (one group): the hand-built G should be reflected")
    out = dwlsmv_signcanon(th, p2)
    nfail = nfail + tsc_chk(max(abs(out - th)) == 0, "k: with a second group present pass 2 must not run")
    printf("{txt}  k: one-group control reflected, grp with a 2 unchanged\n")
    return(nfail)
}

end

local nfail 0
foreach j in a b c d e f g h i k {
    capture noisily mata: st_local("nf", strofreal(tsc_`j'()))
    if _rc {
        display as error "  `j': aborted with rc " _rc
        local ++nfail
    }
    else if `nf' {
        local nfail = `nfail' + `nf'
    }
}

if `nfail' {
    display as error "test_signcanon_so: `nfail' failure(s)"
    exit 9
}
display "test_signcanon_so: pass"
