*! Module 7, release B (unit B5): the start rule for an item-anchored
*! higher-order factor, dwlsmv_hof_start(), and its effect on the fit.
*!
*! WHAT. dwlsmv_fit_uls() starts the ULS run from par.value; for a factor f
*! with a free psi_ff, exactly one lambda row and two or more free outgoing
*! beta rows (MD in metric V: the second-order factor G anchored on one item,
*! its scale carried by the marker, psi_GG free) dwlsmv_hof_start() moves
*! those betas from Module 3's 0 to 0.7. Module 3's rank check never sees it.
*!
*! WHY. With such a G, every implied correlation has a finite limit along
*! c_k = b_k psi_G fixed, |psi_G| -> infinity (a ridge), and the default betas
*! (b = 0) sit on the ray's direction. Measured on Case H: G declared last,
*! the ULS stage runs to its 300 cap and the fit fails at 1 and 4 processors
*! (E-114, the trace log); from 0.7 it takes 8
*! iterations at every processor count (E-114, the prototype's log). Decided by the
*! owner, 2026-11-20 (A-R, E-114, the start package,
*! section 3); the condition is derived in @sec-so-hof (Result SO7).
*!
*! Checks.
*!  1. Table level, no fit (h1..h7): MD in metric V, G declared first and G
*!     declared last: the helper returns 4, exactly the four free betas out of
*!     G become 0.7, every other start is unchanged. Negative controls, each
*!     returns 0 and changes nothing: (a) two items on G; (b) G with no item
*!     (metric unitvariance G); (c) MD in metric U (psi_GG fixed); (d) Case E,
*!     marker metric, F1 with three items; (e) a factor with one lambda row, a
*!     free psi and ONE outgoing beta. For (e) dwls_identify() cannot be used
*!     (that model is deficient and the rank check refuses it: psi_1, b and
*!     psi_2 enter only through two products), so the parameter table is
*!     edited by hand from Case E: F1's lambda rows for y2 and y3 are moved to
*!     column F2. That is table surgery on a copy, used only to read the
*!     helper's condition; nothing is fitted from it.
*!  2. Fit level (the point of the decision): MD on Case H, metric V, both
*!     declaration orders, at processors 1, 4 and c(processors): dwls_fit()
*!     converges; the ULS stage replicated from the helper-adjusted start
*!     (dwlsmv_fit_uls() does not return its iteration count, so the stage is
*!     replicated, as the prototype in E-114 does)
*!     converges in at most 100 iterations; fmin agrees across all six runs to
*!     1e-12 relative; theta-hat agrees across processor counts within each
*!     order to 1e-8 absolute, and across the two orders after sorting (the
*!     parameter order differs with the declaration order; G's sign is fixed by
*!     its item, so no reflection is expected).
*!     The bound 100 was chosen after measuring 8 iterations with the
*!     prototype (E-114, the prototype's log); the cap is
*!     300. It is a margin, not a derived number.
*!
*! Mutant this must kill. Removing the call
*!     (void) dwlsmv_hof_start(par, start)
*! in dwlsmv_fit_uls() (src/mata/dwls_fit.mata): check 2 fails for the G-last
*! order at processors 1 and 4 (baseline: ULS 301 iterations, fit not
*! converged, E-114). The table-level checks are killed by the helper
*! returning anything but 4 / 0, or moving any other start. A weakened
*! condition is killed by the controls: "two or more free betas" dropped to
*! "one" kills (e); "exactly one lambda" dropped kills (a); "psi free" dropped
*! kills (b) and (c).
*!
*! Return codes: pass 0; failure 9 after the failures are printed; no skips.
*! Error paths do not arise here (an abort inside a Mata function reaches the
*! caller's capture noisily, which reports the rc and counts one failure).
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

import delimited using "fixtures/data/case_H.csv", clear varnames(1) delimiters(",")
mata: YH = st_data(., .)

mata:

// Four first-order factors F1..F4 (items y1..y12, three each, loading on G),
// and G, in the stated declaration order. fmet: "V" (item marker on each F_k,
// psi_k free) or "U" (all loadings free, unitvariance). gkind:
//   "V1" one item y13 on G, metric marker G y13, psi_GG free   (MD, metric V)
//   "V2" two items y13, y14 on G, marker on y13, psi_GG free   (control a)
//   "N"  no item on G, metric unitvariance G                   (control b)
//   "U1" one item y13 on G, free, metric unitvariance G         (control c)
string scalar hs_canon(real scalar gfirst, string scalar fmet, string scalar gkind)
{
    string scalar L, fl, f
    real scalar k, j, y, nitems

    nitems = 12
    if (gkind == "V1" | gkind == "U1") {
        nitems = 13
    }
    if (gkind == "V2") {
        nitems = 14
    }
    fl = " F1 F2 F3 F4"
    if (gfirst) {
        fl = " G" + fl
    }
    else {
        fl = fl + " G"
    }
    L = "factor" + fl + char(10) + "ind"
    for (y = 1; y <= nitems; y++) {
        L = L + " y" + strofreal(y)
    }
    L = L + char(10)
    y = 0
    for (k = 1; k <= 4; k++) {
        f = "F" + strofreal(k)
        for (j = 1; j <= 3; j++) {
            y++
            if (fmet == "V" & j == 1) {
                L = L + "metric marker " + f + " y" + strofreal(y) + char(10)
            }
            else {
                L = L + "lambda " + f + " y" + strofreal(y) + " free" + char(10)
            }
        }
        if (fmet == "V") {
            L = L + "psi " + f + " " + f + " free" + char(10)
        }
        else {
            L = L + "metric unitvariance " + f + char(10)
        }
        L = L + "beta " + f + " G free" + char(10)
    }
    if (gkind == "V1") {
        L = L + "metric marker G y13" + char(10) + "psi G G free" + char(10)
    }
    else if (gkind == "V2") {
        L = L + "metric marker G y13" + char(10) + "lambda G y14 free" + char(10)
        L = L + "psi G G free" + char(10)
    }
    else if (gkind == "U1") {
        L = L + "lambda G y13 free" + char(10) + "metric unitvariance G" + char(10)
    }
    else {
        L = L + "metric unitvariance G" + char(10)
    }
    return(L)
}

// Case E (6 items, 3 categories), marker metric: F1 = y1..y3, F2 = y4..y6,
// beta F2 F1 free.
string scalar hs_caseE()
{
    string scalar L

    L = "factor F1 F2" + char(10) + "ind y1 y2 y3 y4 y5 y6" + char(10)
    L = L + "metric marker F1 y1" + char(10) + "lambda F1 y2 free" + char(10)
    L = L + "lambda F1 y3 free" + char(10)
    L = L + "metric marker F2 y4" + char(10) + "lambda F2 y5 free" + char(10)
    L = L + "lambda F2 y6 free" + char(10)
    L = L + "psi F1 F1 free" + char(10) + "psi F2 F2 free" + char(10)
    return(L + "beta F2 F1 free" + char(10))
}

// theta start as dwlsmv_fit_uls() builds it, from par.value.
real colvector hs_start(struct dwlsmv_parmap scalar par)
{
    real colvector start
    real scalar i

    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            start[par.fidx[i]] = par.value[i]
        }
    }
    return(start)
}

// Table-level case. expect = the theta positions that must become 0.7 (empty
// for a negative control); every other entry must equal its pre-call value.
// Returns the number of failed assertions.
real scalar hs_table(string scalar tag, struct dwlsmv_parmap scalar par,
                     real colvector expect)
{
    real colvector start, before, want
    real scalar nfail, ret, i

    start  = hs_start(par)
    before = start
    want   = before
    for (i = 1; i <= rows(expect); i++) {
        want[expect[i]] = 0.7
    }
    ret   = dwlsmv_hof_start(par, start)
    nfail = 0
    if (ret != rows(expect)) {
        printf("{err}  %s: helper returned %f, expected %f\n", tag, ret, rows(expect))
        nfail++
    }
    if (any(start :!= want)) {
        printf("{err}  %s: start differs from the expected start in %f entries\n",
               tag, sum(start :!= want))
        nfail++
    }
    printf("{txt}  %s: returned %f, %f starts changed\n", tag, ret,
           sum(start :!= before))
    return(nfail)
}

// theta positions of the free betas whose col is the factor index gcol.
real colvector hs_betas_out(struct dwlsmv_parmap scalar par, real scalar gcol)
{
    real colvector out
    real scalar i

    out = J(0, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 3 & par.free[i] & par.col[i] == gcol) {
            out = out \ par.fidx[i]
        }
    }
    return(out)
}

real scalar hs_md(string scalar tag, real scalar gfirst)
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par
    real colvector e
    real scalar gcol, nfail

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(13, 1, 4))
    par  = dwls_identify(hs_canon(gfirst, "V", "V1"), m)
    gcol = (gfirst ? 1 : 5)
    e    = hs_betas_out(par, gcol)
    nfail = 0
    if (rows(e) != 4) {
        printf("{err}  %s: %f free betas out of factor %f, expected 4 (table not as built)\n",
               tag, rows(e), gcol)
        nfail++
    }
    return(nfail + hs_table(tag, par, e))
}

real scalar hs_h1()
{
    return(hs_md("h1 MD, V, G first", 1))
}

real scalar hs_h2()
{
    return(hs_md("h2 MD, V, G last", 0))
}

real scalar hs_neg(string scalar tag, real scalar p, string scalar fmet,
                   string scalar gkind, real scalar gfirst)
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(p, 1, 4))
    par  = dwls_identify(hs_canon(gfirst, fmet, gkind), m)
    return(hs_table(tag, par, J(0, 1, .)))
}

real scalar hs_h3()
{
    // (a) two items on G: lambda_a lambda_b psi_G is a fitted correlation.
    return(hs_neg("h3 two items on G", 14, "V", "V2", 0))
}

real scalar hs_h4()
{
    // (b) G with no item, metric unitvariance G: psi_GG fixed.
    return(hs_neg("h4 no item on G", 12, "V", "N", 0))
}

real scalar hs_h5()
{
    // (c) MD in metric U: psi_GG fixed at 1.
    return(hs_neg("h5 MD, metric U", 13, "U", "U1", 0))
}

real scalar hs_h6()
{
    // (d) Case E, marker metric: F1 has three lambda rows and one free beta.
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(6, 1, 3))
    par  = dwls_identify(hs_caseE(), m)
    return(hs_table("h6 Case E marker metric", par, J(0, 1, .)))
}

real scalar hs_h7()
{
    // (e) one lambda row, free psi, ONE outgoing beta: Case E with F1's
    // lambda rows for items 2 and 3 moved to column F2 (table surgery on a
    // copy; the model is deficient and is never fitted or identified).
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar par
    real scalar i, nlam, nbeta, npsi, nfail

    m    = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(6, 1, 3))
    par  = dwls_identify(hs_caseE(), m)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1 & par.col[i] == 1 & (par.row[i] == 2 | par.row[i] == 3)) {
            par.col[i] = 2
        }
    }
    nlam = 0
    nbeta = 0
    npsi = 0
    for (i = 1; i <= par.nall; i++) {
        nlam  = nlam  + (par.ptype[i] == 1 & par.col[i] == 1)
        nbeta = nbeta + (par.ptype[i] == 3 & par.free[i] & par.col[i] == 1)
        npsi  = npsi  + (par.ptype[i] == 2 & par.free[i] & par.row[i] == 1 & par.col[i] == 1)
    }
    nfail = 0
    if (nlam != 1 | nbeta != 1 | npsi != 1) {
        printf("{err}  h7: surgery left %f lambda rows, %f free betas, %f free psi on F1; expected 1, 1, 1\n",
               nlam, nbeta, npsi)
        nfail++
    }
    return(nfail + hs_table("h7 one lambda, one beta (edited table)", par, J(0, 1, .)))
}

// One fit-level run. Appends a column (order \ fmin \ theta) to the external
// HS and returns the number of failed assertions for this run.
real scalar hs_fit(string scalar tag, real scalar gfirst, real scalar p)
{
    external real matrix YH, HS
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f
    real colvector s, w, start, thu, col
    real scalar N, nfail, cu, fu, iu, nch

    nfail = 0
    st  = dwls_sample_stats(YH, J(cols(YH), 1, 4))
    sts = dwlsmv_stage1(1, 1)
    sts[1] = st
    par = dwls_identify(hs_canon(gfirst, "V", "V1"), st.map)
    dwlsmv_fit_stack(sts, s, w, N, maps, 1)
    start = hs_start(par)
    nch   = dwlsmv_hof_start(par, start)
    if (nch != 4) {
        printf("{err}  %s p=%f: helper returned %f on Case H, expected 4\n", tag, p, nch)
        nfail++
    }
    cu = dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, start, "d2", "nr", 300, thu, fu, iu)
    // 100 chosen after measuring 8 iterations with the prototype
    // (E-114, the prototype's log); the cap is 300.
    if (cu != 1 | iu > 100) {
        printf("{err}  %s p=%f: ULS from the helper start: converged %f, iterations %f (want 1 and <= 100)\n",
               tag, p, cu, iu)
        nfail++
    }
    f = dwls_fit(sts, par)
    if (f.converged != 1) {
        printf("{err}  %s p=%f: dwls_fit did not converge (iterations %f)\n", tag, p, f.iterations)
        nfail++
    }
    printf("{txt}  %s p=%f: ULS converged %f in %3.0f iterations; fit converged %f, %3.0f iterations, fmin %15.10g\n",
           tag, p, cu, iu, f.converged, f.iterations, f.fmin)
    col = (gfirst \ f.fmin \ f.theta)
    if (cols(HS) == 0) {
        HS = col
    }
    else {
        HS = HS, col
    }
    return(nfail)
}

// Agreement across the stored runs: fmin to 1e-12 relative over all of them,
// theta to 1e-8 absolute within each order, and across orders after sorting
// (the parameter order depends on the declaration order).
real scalar hs_compare()
{
    external real matrix HS
    real scalar j, nfail, d, ref, a
    real colvector t0, t1, tref

    nfail = 0
    if (cols(HS) < 2) {
        printf("{err}  hs_compare: %f runs stored\n", cols(HS))
        return(1)
    }
    ref = HS[2, 1]
    for (j = 2; j <= cols(HS); j++) {
        d = abs(HS[2, j] - ref) / abs(ref)
        if (d > 1e-12) {
            printf("{err}  fmin run %f differs from run 1 by %9.2e relative\n", j, d)
            nfail++
        }
    }
    tref = sort(HS[|3, 1 \ ., 1|], 1)
    for (j = 2; j <= cols(HS); j++) {
        t1 = HS[|3, j \ ., j|]
        for (a = 1; a < j; a++) {
            if (HS[1, a] == HS[1, j]) {
                t0 = HS[|3, a \ ., a|]
                if (max(abs(t1 - t0)) > 1e-8) {
                    printf("{err}  theta run %f differs from run %f (same order) by %9.2e\n",
                           j, a, max(abs(t1 - t0)))
                    nfail++
                }
                break
            }
        }
        d = max(abs(sort(t1, 1) - tref))
        if (d > 1e-8) {
            printf("{err}  sorted theta run %f differs from run 1 by %9.2e\n", j, d)
            nfail++
        }
    }
    printf("{txt}  %f runs compared; fmin spread %9.2e relative\n", cols(HS),
           (max(HS[2, .]) - min(HS[2, .])) / abs(ref))
    return(nfail)
}

end

local nfail 0

* ---------------------------------------------------------------------------
* 1. Table level, no fit.
* ---------------------------------------------------------------------------
forvalues j = 1/7 {
    capture noisily mata: st_local("nf", strofreal(hs_h`j'()))
    if _rc {
        display as error "  h`j': aborted with rc " _rc
        local ++nfail
    }
    else if `nf' {
        local nfail = `nfail' + `nf'
    }
}

* ---------------------------------------------------------------------------
* 2. Fit level: MD on Case H, metric V, both declaration orders.
* ---------------------------------------------------------------------------
mata: HS = J(0, 0, .)
local c0 = c(processors)
foreach p in 1 4 `c0' {
    set processors `p'
    foreach o in 1 0 {
        local tag = cond(`o', "MD V G-first", "MD V G-last")
        capture noisily mata: st_local("nf", strofreal(hs_fit("`tag'", `o', `p')))
        if _rc {
            display as error "  `tag' p=`p': aborted with rc " _rc
            local ++nfail
        }
        else if `nf' {
            local nfail = `nfail' + `nf'
        }
    }
}
set processors `c0'
capture noisily mata: st_local("nf", strofreal(hs_compare()))
if _rc {
    display as error "  hs_compare aborted with rc " _rc
    local ++nfail
}
else if `nf' {
    local nfail = `nfail' + `nf'
}

if `nfail' {
    display as error "test_fit_hof_start: `nfail' failure(s)"
    exit 9
}
display "test_fit_hof_start: pass"
