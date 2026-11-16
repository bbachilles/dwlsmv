*! Investigation D's 14 inside-D flush-event pairs (cell 180) as a regression
*! fixture: fixtures/expected/phi2_flush_events.csv and phi2_flush_corners.csv,
*! from fixtures/generate/99l_phi2_flush_events.R (design:
*! E-099). Tests do not call R.
*!
*!   A0  Arm 0, production, a CHARACTERIZATION and not a correctness claim
*!       (the production question on flushes is open): on each pair the
*!       production solver returns the fixture's rho-hat (exactly), status,
*!       and its ends xl0/xh0 are usable points of its own walk whose
*!       rf_check classes are cl0/ch0. A future true-corner evaluator is
*!       expected to fail this block; that is the point of it.
*!   R   The repair: the solver input equals the crosstab after the
*!       production empty-cell repair (if it applies).
*!   F   The flush itself: at every listed corner binormal() is exactly 0
*!       while the certified reference (Method C, checked by M4) is
*!       positive. Read with the corner budget's own scale: ref_hi > 0.
*!   S1  Arm 1 is static study evidence, not executed here and NOT certified
*!       (user, 2026-11-09: the corner values are certified, the arm-1 decision
*!       is not; it rests on the M4 replay alone). Checked only as recorded: its
*!       status differs from arm 0's on every row, and the replay matched.
*!   E   Aborts at rc 3498 with their own message, on this fixture's inputs:
*!       no thresholds for the first item; rho at 1.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

clear all
set more off
set linesize 255

foreach f in src/mata/dwlsmv_types.mata                ///
             src/mata/dwlsmv_statmap_build.mata        ///
             src/mata/dwlsmv_pos_rho.mata              ///
             src/mata/dwlsmv_catcounts.mata            ///
             src/mata/dwlsmv_catdegen.mata             ///
             src/mata/dwlsmv_cumprop.mata              ///
             src/mata/dwlsmv_thresholds.mata           ///
             src/mata/dwlsmv_psi_tau_table.mata        ///
             src/mata/dwlsmv_cellprob.mata             ///
             src/mata/dwlsmv_dcellprob_rho.mata        ///
             src/mata/dwlsmv_wmat.mata                 ///
             src/mata/dwlsmv_crosstab.mata             ///
             src/mata/dwlsmv_pairinfo.mata             ///
             src/mata/dwls_polychoric.mata             ///
             src/mata/dwlsmv_emptycell1.mata           ///
             src/mata/dwlsmv_emptycell_repair.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_phi2_flush_events: `f' failed to load"
        exit 9
    }
}

mata:

// Production evaluator, unmodified; each evaluation kept.
real rowvector xf_fe(real scalar x, real matrix n, real colvector a, real colvector b)
{
    external real matrix XFT
    real rowvector e

    e = dwlsmv_rf_eval(x, n, a, b)
    XFT = XFT \ e[1..5]
    return(e)
}

// rf_check's reading of an end (as the sv_s3 of E-100).
real scalar xf_class(real rowvector e)
{
    if (e[3] > e[5]) {
        return(1)
    }
    if (e[3] < -e[5]) {
        return(-1)
    }
    return(0)
}

real scalar xf_col(string rowvector nm, string scalar v)
{
    real colvector i

    i = selectindex(nm :== v)
    if (length(i) != 1) {
        _error("test_phi2_flush_events: fixture column missing")
    }
    return(i[1])
}

real matrix xf_tab(string rowvector r, string rowvector nm, string scalar p)
{
    real matrix M
    real scalar i, j

    M = J(3, 3, .)
    for (i = 1; i <= 3; i++) {
        for (j = 1; j <= 3; j++) {
            M[i, j] = strtoreal(r[xf_col(nm, p + strofreal(i) + strofreal(j))])
        }
    }
    return(M)
}

string rowvector xf_names()
{
    real scalar i
    string rowvector nm

    nm = J(1, st_nvar(), "")
    for (i = 1; i <= st_nvar(); i++) {
        nm[i] = st_varname(i)
    }
    return(nm)
}

real scalar xf_events()
{
    external real matrix XFT
    string matrix E
    string rowvector nm, r
    real colvector a, b
    real matrix n, f, nf
    real rowvector rr
    real scalar e, bad, nbad, s0, s1, il, ih

    E = st_sdata(., .)
    nm = xf_names()
    if (rows(E) != 14) {
        printf("{err}  fixture rows %f, not 14\n", rows(E))
        return(1)
    }
    nbad = 0
    for (e = 1; e <= rows(E); e++) {
        r = E[e, .]
        bad = 0
        a = strtoreal(r[xf_col(nm, "a1")]) \ strtoreal(r[xf_col(nm, "a2")])
        b = strtoreal(r[xf_col(nm, "b1")]) \ strtoreal(r[xf_col(nm, "b2")])
        n = xf_tab(r, nm, "n")
        f = xf_tab(r, nm, "f")
        // R
        nf = (dwlsmv_emptycell1(n) ? dwlsmv_emptycell_repair(n) : n)
        if (nf != f) bad++
        // A0
        XFT = J(0, 5, .)
        rr = dwlsmv_rootfind(&xf_fe(), f, a, b, (., .), 0.999, 40, 54)
        s0 = strtoreal(r[xf_col(nm, "s0")])
        if (rr[1] != strtoreal(r[xf_col(nm, "rho0")])) bad++
        if (rr[2] != s0) bad++
        il = selectindex(XFT[., 1] :== strtoreal(r[xf_col(nm, "xl0")]) :& XFT[., 2])
        ih = selectindex(XFT[., 1] :== strtoreal(r[xf_col(nm, "xh0")]) :& XFT[., 2])
        if (rows(il) == 0 | rows(ih) == 0) {
            bad++
        }
        else {
            if (xf_class(XFT[il[rows(il)], .]) != strtoreal(r[xf_col(nm, "cl0")])) bad++
            if (xf_class(XFT[ih[rows(ih)], .]) != strtoreal(r[xf_col(nm, "ch0")])) bad++
        }
        // S1
        s1 = strtoreal(r[xf_col(nm, "s1")])
        if (s1 == s0 | s1 >= .) bad++
        if (r[xf_col(nm, "replay_same")] != "1") bad++
        if (bad) {
            printf("{err}  run %s pair %s: %f check(s) failed (status %f, rho %21.17g)\n",
                   r[xf_col(nm, "run")], r[xf_col(nm, "k")], bad, rr[2], rr[1])
        }
        else {
            printf("{txt}  run %s pair %s: status %f, ends %s/%s, as fixed\n", r[xf_col(nm, "run")],
                   r[xf_col(nm, "k")], rr[2], r[xf_col(nm, "cl0")], r[xf_col(nm, "ch0")])
        }
        nbad = nbad + bad
    }
    return(nbad)
}

real scalar xf_corners()
{
    string matrix E
    string rowvector nm
    real scalar e, bad, h, k, rho, ref

    E = st_sdata(., .)
    nm = xf_names()
    bad = 0
    if (rows(E) < 14) {
        printf("{err}  corner rows %f, fewer than one per pair\n", rows(E))
        return(1)
    }
    for (e = 1; e <= rows(E); e++) {
        h = strtoreal(E[e, xf_col(nm, "h")])
        k = strtoreal(E[e, xf_col(nm, "kc")])
        rho = strtoreal(E[e, xf_col(nm, "rho")])
        ref = strtoreal(E[e, xf_col(nm, "ref_hi")])
        if (binormal(h, k, rho) != 0 | !(ref > 0) | abs(h) > 2.5 | abs(k) > 2.5) {
            printf("{err}  corner row %f: binormal %21x, reference %21x\n", e,
                   binormal(h, k, rho), ref)
            bad++
        }
    }
    printf("{txt}  %f flushed corners, each exactly 0 in binormal(), each positive in the reference\n",
           rows(E))
    return(bad)
}

// E: case 1 no thresholds for the first item; case 2 rho = 1.
void xf_abort(real scalar c)
{
    string matrix E
    string rowvector nm, r
    real colvector a, b
    real matrix f

    E = st_sdata(., .)
    nm = xf_names()
    r = E[1, .]
    a = strtoreal(r[xf_col(nm, "a1")]) \ strtoreal(r[xf_col(nm, "a2")])
    b = strtoreal(r[xf_col(nm, "b1")]) \ strtoreal(r[xf_col(nm, "b2")])
    f = xf_tab(r, nm, "f")
    if (c == 1) {
        (void) dwlsmv_cellprob_raw(J(0, 1, .), b, -0.5)
    }
    else {
        (void) dwlsmv_rf_eval(1, f, a, b)
    }
}

end

local nfail = 0
local nerr  = 0

import delimited using fixtures/expected/phi2_flush_events.csv, stringcols(_all) delimiters(",") varnames(1) clear
mata: st_local("w", strofreal(xf_events()))
local nfail = `nfail' + `w'

local msg1 "dwlsmv_cellprob: each variable needs at least one threshold"
local msg2 "dwlsmv_cellprob: rho must lie strictly inside (-1, 1)"
tempfile ablog
forvalues i = 1/2 {
    quietly log using "`ablog'", text replace name(xfabort)
    capture noisily mata: xf_abort(`i')
    local rc = _rc
    quietly log close xfabort
    if `rc' != 3498 {
        display as error "test_phi2_flush_events: abort case `i' gave rc `rc', not 3498"
        local nerr = `nerr' + 1
    }
    mata: st_local("hit", strofreal(any(strpos(cat(st_local("ablog")), st_local("msg`i'")))))
    if !`hit' {
        display as error "test_phi2_flush_events: abort case `i' did not report: `msg`i''"
        local nerr = `nerr' + 1
    }
}

import delimited using fixtures/expected/phi2_flush_corners.csv, stringcols(_all) delimiters(",") varnames(1) clear
mata: st_local("w", strofreal(xf_corners()))
local nfail = `nfail' + `w'

if `nfail' | `nerr' {
    display as error _n "test_phi2_flush_events: `nfail' failure(s), `nerr' error-path failure(s)"
    exit 9
}
display as result "test_phi2_flush_events: production reproduces all 14 flush-event pairs; every corner flushed; aborts at rc 3498"
