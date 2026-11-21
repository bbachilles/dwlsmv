*! Module 3, release B (unit B5, first step): dwls_identify() accepts a factor
*! measured only by other factors, from hand-written canonical text.
*!
*! Executed under StataMP 18 on 2026-11-20 and it passes. Mutants run the same
*! day in a scratch copy, each killed by its intended case (see the end of this
*! header): lift removed, 10 failures; nonzero condition dropped, r2; mkind
*! condition dropped, r8.
*!
*! Spec Module 12 *Release B*, B3.1 (decided 2026-11-19): a factor with no
*! observed indicator is scaled by its variance -- metric U ('metric
*! unitvariance' on every factor) or metric V (item markers on the first-order
*! factors, 'metric unitvariance G'). Module 3 lifts only its "a factor has no
*! indicators" abort, for a factor with an outgoing beta (free or fixed nonzero)
*! and a 'metric unitvariance' line. Identification is the rank check's.
*! Nothing here touches Module 11 (the seam rule).
*!
*! Accepted (a1-a8): G over 4 and over 3 in U and V; two correlated G's (M33);
*! an item also loading on G (MD, already accepted in 0.5.0); a factor outside
*! the hierarchy covarying with G; a two-item first-order factor. For each, the
*! table's counts are asserted, rank(Delta) is recomputed at the start the
*! check left in par.value and asserted full, and the smallest singular-value
*! ratio is printed (a diagnostic; the criterion is dwlsmv_numrank()'s). For
*! a1-a5, where @sec-second-order SO2 establishes that theta0 is deficient
*! (B = 0), the second point having been taken is asserted (every free beta
*! start moved off 0); for a6-a8 it is printed only.
*!
*! Rejected, rc 3498 AND the reason's text (r1-r8): no loadings and no
*! outgoing path; an outgoing beta fixed at 0 only; no metric line; an
*! endogenous factor with no loadings and no children; G over two with no
*! constraint, a one-item first-order factor and two uncorrelated G's (all
*! three by the rank check, B3.4); and a later group's factor with an outgoing
*! beta and no loadings (the 'metric unitvariance' condition: later groups
*! carry no metric line, and group() is refused with a second-order factor,
*! B3.5).
*!
*! B3.13 (added 2026-11-20): on r5-r7 the message prints theta2's null
*! directions under a diagnostic label; the direction counts (1, 1, 2 at
*! theta2 against 2, 4, 4 at theta0) tell the two points apart.
*!
*! Mutants this is meant to kill: the lift removed (a1-a5, a7, a8 abort); the
*! nonzero condition on the beta dropped (r2 accepted); the mkind condition
*! dropped (r8 reaches a different abort or is accepted).
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

string scalar tso_join(string colvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= rows(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

string scalar tso_slurp(string scalar path)
{
    real scalar fh
    string scalar line, out

    fh  = fopen(path, "r")
    out = ""
    while ((line = fget(fh)) != J(0, 0, "")) {
        out = out + line + char(10)
    }
    fclose(fh)
    return(out)
}

string scalar tso_gname(real scalar g, real scalar nG)
{
    if (nG == 1) {
        return("G")
    }
    return("G" + strofreal(g))
}

// First-order factors F1..FK with nper[k] items each, in metric met ("U" or
// "V"); parent[k] is the higher-order factor F_k loads on (0: none). nG
// higher-order factors, each with 'metric unitvariance' when gmet is 1. xitems
// extra items are declared after the first-order ones, for the caller's own
// lines (an item on G).
string colvector tso_lines(real rowvector nper, real rowvector parent,
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
        fl = fl + " " + tso_gname(g, nG)
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
            L = L \ ("beta " + f + " " + tso_gname(parent[k], nG) + " free")
        }
    }
    if (gmet) {
        for (g = 1; g <= nG; g++) {
            L = L \ ("metric unitvariance " + tso_gname(g, nG))
        }
    }
    return(L)
}

// One accepted case. Returns the number of failed assertions; an abort inside
// dwls_identify() is caught by the caller's capture and counted there.
real scalar tso_accept(string scalar tag, string colvector L, real scalar p,
                       real scalar wantnpar, real scalar wantnall,
                       real scalar mustmove)
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par
    real colvector theta, sv
    real matrix V, Dm
    real scalar i, rk, nfail, nb, nmoved

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(p, 1, 3))
    par  = dwls_identify(tso_join(L), m)

    nfail = 0
    if (par.npar != wantnpar | par.nall != wantnall) {
        printf("{err}  %s: npar %f nall %f, expected %f and %f\n", tag,
               par.npar, par.nall, wantnpar, wantnall)
        nfail++
    }

    theta = J(par.npar, 1, .)
    nb = 0
    nmoved = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            theta[par.fidx[i]] = par.value[i]
            if (par.ptype[i] == 3) {
                nb++
                nmoved = nmoved + (par.value[i] != 0)
            }
        }
    }
    sv = J(0, 1, 0)
    V  = J(0, 0, 0)
    Dm = dwls_delta(theta, par, m)
    rk = dwlsmv_numrank(Dm, sv, V)
    printf("{txt}  %s: npar %f, rank %f at the start left in par.value, smallest ratio %9.2e; ",
           tag, par.npar, rk, sv[par.npar] / sv[1])
    printf("{txt}%f of %f free beta starts moved off 0\n", nmoved, nb)
    if (rk != par.npar) {
        printf("{err}  %s: rank %f < npar %f at the start\n", tag, rk, par.npar)
        nfail++
    }
    if (mustmove & (nb == 0 | nmoved != nb)) {
        printf("{err}  %s: the second point was not taken (theta0 should be deficient, SO2)\n", tag)
        nfail++
    }
    return(nfail)
}

real scalar tso_a1()
{
    // G over 4, U: 24 tau + 12 lambda + 4 beta free; 5 psi fixed.
    return(tso_accept("a1 G over 4, U", tso_lines((3,3,3,3), (1,1,1,1), 1, "U", 1, 0),
                      12, 40, 45, 1))
}

real scalar tso_a2()
{
    // G over 4, V: 24 tau + 8 lambda + 4 beta + 4 psi free; 4 marker lambda
    // and psi_GG fixed.
    return(tso_accept("a2 G over 4, V", tso_lines((3,3,3,3), (1,1,1,1), 1, "V", 1, 0),
                      12, 40, 45, 1))
}

real scalar tso_a3()
{
    return(tso_accept("a3 G over 3, U", tso_lines((3,3,3), (1,1,1), 1, "U", 1, 0),
                      9, 30, 34, 1))
}

real scalar tso_a4()
{
    return(tso_accept("a4 G over 3, V", tso_lines((3,3,3), (1,1,1), 1, "V", 1, 0),
                      9, 30, 34, 1))
}

real scalar tso_a5()
{
    // M33, U: 36 tau + 18 lambda + 6 beta + psi_G2G1 free; 8 psi fixed.
    return(tso_accept("a5 two correlated G, U",
                      tso_lines((3,3,3,3,3,3), (1,1,1,2,2,2), 2, "U", 1, 0) \ "psi G2 G1 free",
                      18, 61, 69, 1))
}

real scalar tso_a6()
{
    // MD, U: G over 4 plus y13 on G. G has a loading, so this is the 0.5.0
    // path; a regression pin. 26 tau + 13 lambda + 4 beta.
    return(tso_accept("a6 item on G, U",
                      tso_lines((3,3,3,3), (1,1,1,1), 1, "U", 1, 1) \ "lambda G y13 free",
                      13, 43, 48, 0))
}

real scalar tso_a7()
{
    // G over 3 and F4 outside the hierarchy covarying with G, U:
    // 24 tau + 12 lambda + 3 beta + psi_F4,G free; 5 psi fixed.
    return(tso_accept("a7 outside factor with G, U",
                      tso_lines((3,3,3,3), (1,1,1,0), 1, "U", 1, 0) \ "psi F4 G free",
                      12, 40, 45, 0))
}

real scalar tso_a8()
{
    // G over 3, F3 with two items, U: 16 tau + 8 lambda + 3 beta.
    return(tso_accept("a8 two-item factor, U", tso_lines((3,3,2), (1,1,1), 1, "U", 1, 0),
                      8, 27, 31, 0))
}

void tso_reject1(string colvector L, real scalar p)
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(p, 1, 3))
    par  = dwls_identify(tso_join(L), m)
}

void tso_r1()
{
    // G declared and scaled, but no loadings and no path out.
    tso_reject1(tso_lines((3,3,3), (0,0,0), 1, "U", 1, 0), 9)
}

void tso_r2()
{
    // G's only path out is fixed at 0, which is no path.
    tso_reject1(tso_lines((3,3,3), (0,0,0), 1, "U", 1, 0) \ "beta F1 G fixed 0", 9)
}

void tso_r3()
{
    // Paths out, but no metric line for G.
    tso_reject1(tso_lines((3,3,3), (1,1,1), 1, "U", 0, 0), 9)
}

void tso_r4()
{
    // G regressed on F1, with no loadings and no children: endogenous,
    // measured by nothing.
    tso_reject1(tso_lines((3,3,3), (0,0,0), 1, "U", 1, 0) \ "beta G F1 free", 9)
}

void tso_r5()
{
    // G over two, unconstrained: deficient at both points (SO2, B3.4).
    tso_reject1(tso_lines((3,3), (1,1), 1, "U", 1, 0), 6)
}

void tso_r6()
{
    // G over 3 with a one-item first-order factor (B3.4).
    tso_reject1(tso_lines((3,3,1), (1,1,1), 1, "U", 1, 0), 7)
}

void tso_r7()
{
    // Two uncorrelated G's over two factors each (B3.4).
    tso_reject1(tso_lines((3,3,3,3), (1,1,2,2), 2, "U", 1, 0), 12)
}

void tso_r8()
{
    // Two groups. Group 1's G is accepted (metric unitvariance, paths out);
    // group 2's G has paths out but, as a later group's factor, no metric
    // line, so the "no indicators" abort stands for it.
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st
    struct dwlsmv_parmap scalar par
    string colvector L, L2
    real scalar g, k

    m  = dwlsmv_statmap(1, 2)
    st = dwlsmv_stage1(1, 2)
    for (g = 1; g <= 2; g++) {
        m[g]      = dwlsmv_statmap_build(J(9, 1, 3))
        st[g].map = m[g]
        st[g].s   = J(m[g].D, 1, 0.3)
    }
    L  = tso_lines((3,3,3), (1,1,1), 1, "U", 1, 0)
    L2 = ("factor F1 F2 F3 G" \ "ind y1 y2 y3 y4 y5 y6 y7 y8 y9")
    for (k = 1; k <= 9; k++) {
        L2 = L2 \ ("lambda F" + strofreal(ceil(k / 3)) + " y" + strofreal(k) + " free")
    }
    for (k = 1; k <= 3; k++) {
        L2 = L2 \ ("psi F" + strofreal(k) + " F" + strofreal(k) + " fixed 1")
        L2 = L2 \ ("beta F" + strofreal(k) + " G free")
    }
    L2 = L2 \ "psi G G fixed 1"
    par = dwls_identify(tso_join(("groups 2" \ "group 1" \ L \ "group 2" \ L2)), m, st)
}

// B6 review (2026-11-20; added after its mutant I3 survived): r2's twin with the
// path fixed at -1. Decision 2327: a fixed NONZERO beta counts as a path, of either
// sign, so the lift applies and the table is accepted (r2, fixed at 0, is refused).
void tso_neg()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(9, 1, 3))
    par  = dwls_identify(tso_join(tso_lines((3,3,3), (0,0,0), 1, "U", 1, 0) \
                                  "beta F1 G fixed -1"), m)
}
end

local nfail 0

* ---------------------------------------------------------------------------
* Accepted.
* ---------------------------------------------------------------------------
forvalues j = 1/8 {
    capture noisily mata: st_local("nf", strofreal(tso_a`j'()))
    if _rc {
        display as error "  a`j': dwls_identify() aborted with rc " _rc
        local ++nfail
    }
    else if `nf' {
        local nfail = `nfail' + `nf'
    }
}

* ---------------------------------------------------------------------------
* Rejected: rc 3498 and the reason's own text.
* ---------------------------------------------------------------------------
local why1 "has no loadings, so nothing observed measures it"
local why2 "has no loadings, so nothing observed measures it"
local why3 "has no metric directive"
local why4 "has no loadings, so nothing observed measures it"
local why5 "rank(Delta) is"
local why6 "rank(Delta) is"
local why7 "rank(Delta) is"
local why8 "has no loadings, so nothing observed measures it"

* B3.13 (owner, 2026-11-20): on r5-r7 the printed null directions are theta2's,
* labelled a diagnostic and not a proof. theta0 is deficient by 2, 4 and 4
* (B = 0 adds K in U), theta2 by 1, 1 and 2, so the counts tell the points
* apart: want`j' must appear and nowant`j' must not.
local diag "not a proof that the model is unidentified"
local want5 "Set 1 of 1:"
local nowant5 "Set 1 of 2:"
local want6 "Set 1 of 1:"
local nowant6 "Set 1 of 4:"
local want7 "Set 2 of 2:"
local nowant7 "Set 1 of 4:"
* B6 review (2026-11-20; added after its mutants K1 and K2 survived, and written after
* seeing the unmutated output): every set is numbered from 1, so "Set 0 of" and "Set -"
* never appear; and the directions are theta2's in content, not only in count. G over
* two trades its two paths against each other (r5 names both), and a one-item factor's
* deficiency is its own loading (r6 names lambda[y7,F3]); theta0's V names a G path.
local name5a "beta[F1,G]"
local name5b "beta[F2,G]"
local name6a "lambda[y7,F3]"
local name6b "lambda[y7,F3]"
local name7a "beta[F3,G2]"
local name7b "beta[F4,G2]"

forvalues j = 1/8 {
    tempfile lg
    quietly log using "`lg'", text replace name(tso)
    capture noisily mata: tso_r`j'()
    local rc = _rc
    quietly log close tso
    mata: st_local("hit", strofreal(strpos(tso_slurp(st_local("lg")), st_local("why`j'")) > 0))
    display "  r`j': rc `rc', reason text found `hit'"
    if `rc' != 3498 | `hit' != 1 {
        display as error "  r`j' FAILED: want rc 3498 and '`why`j''"
        local ++nfail
    }
    if `j' >= 5 & `j' <= 7 {
        mata: st_local("d1", strofreal(strpos(tso_slurp(st_local("lg")), st_local("diag")) > 0))
        mata: st_local("d2", strofreal(strpos(tso_slurp(st_local("lg")), st_local("want`j'")) > 0))
        mata: st_local("d3", strofreal(strpos(tso_slurp(st_local("lg")), st_local("nowant`j'")) > 0))
        display "  r`j' (B3.13): diagnostic label `d1', theta2 count `d2', theta0 count `d3'"
        if `d1' != 1 | `d2' != 1 | `d3' != 0 {
            display as error "  r`j' FAILED B3.13: want the label and '`want`j''; not '`nowant`j''"
            local ++nfail
        }
        mata: st_local("d4", strofreal(strpos(tso_slurp(st_local("lg")), "Set 0 of") + strpos(tso_slurp(st_local("lg")), "Set -") > 0))
        mata: st_local("d5", strofreal(strpos(tso_slurp(st_local("lg")), st_local("name`j'a")) > 0 & strpos(tso_slurp(st_local("lg")), st_local("name`j'b")) > 0))
        display "  r`j' (B6): sets numbered from 1 `=1-`d4'', theta2's parameters named `d5'"
        if `d4' != 0 | `d5' != 1 {
            display as error "  r`j' FAILED B6: a set not numbered from 1, or '`name`j'a'' / '`name`j'b'' missing"
            local ++nfail
        }
    }
}

* B6: a path out of G fixed at -1 is a path (decision 2327), so the lift applies
capture noisily mata: tso_neg()
if _rc != 0 {
    display as error "  neg FAILED: G with one path fixed at -1 must be accepted, rc " _rc
    local ++nfail
}
else display "  neg: G with one path fixed at -1 accepted"

if `nfail' {
    display as error "test_identify_second_order: `nfail' failure(s)"
    exit 9
}
display "test_identify_second_order: pass"
