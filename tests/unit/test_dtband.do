*! Module 12, DIFFTEST validation check 10: dwlsmv_dtband() classifies T0 - T1
*! against eps = c * nrtol, c the multiplier of F in T (N - G).
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9; executed by the main
*! session 2026-10-08: 20/20 green; mutants B1-B9 (B8 as three) each turned it red.
*! Changed 2026-10-09 (user): the routine takes c, not N. Cases 1-12 keep their
*! numbers with c = 499 (N - 1 at N = 500); 13-17 use c = 998 (N - G at N = 1000,
*! G = 2) and a non-integer c; the old N = 1 / N = 2.5 refusals become c = 0 /
*! c < 0. Mutants B10-B12 added.
*!
*! Rule: spec Module 12 check 10: T0 - T1 < -eps aborts (rc 3498);
*! |T0 - T1| <= eps is reported as 0 (status 0); otherwise proceed (status 1)
*! with Td = T0 - T1. c = 499, nrtol = 1e-12 unless a case says otherwise, so
*! eps = 499e-12. No data and no model: constructed numbers.
*!
*! Cases (T0, T1 are chosen so that T0 - T1 is exact: one of them is 0):
*!   1  Td = 6                       proceed
*!   2  Td = +eps exactly            band     kills "<" for "<=" on the + side
*!   3  Td = -eps exactly            band     kills "<=" on the - side
*!   4  Td = eps (1 + 1e-9)          proceed  kills a band that is too wide
*!   5  Td = -eps (1 + 1e-9)         abort
*!   6  Td = 0                       band
*!   7  Td = +4.995e-10              proceed  between c nrtol and (c+1) nrtol:
*!                                            kills eps = (c + 1) * nrtol
*!   8  Td = -4.995e-10              abort    likewise
*!   9  Td = -eps / 2                band     kills "abort when Td < 0"
*!   10 nrtol = 1e-6, Td = 1e-4      band     eps = 4.99e-4; kills a hard-coded
*!                                            1e-12
*!   11 Td = -3                      abort
*!   12 Td = 506.4                   proceed
*!   13 c = 998, Td = +998e-12       band     c != N - 1 shape: kills (c - 1) nrtol
*!   14 c = 998, Td = -998e-12       band
*!   15 c = 998, Td = +9.985e-10     proceed  between c nrtol and (c+1) nrtol
*!   16 c = 998, Td = -9.985e-10     abort
*!   17 c = 2.5, nrtol = 1e-6,       band     non-integer c accepted (eps 2.5e-6)
*!      Td = 2.4e-6
*!   18-25 error paths, each rc 3498: T0 missing, T1 missing, c missing, c = 0,
*!      c = -5, nrtol = 0, nrtol missing, nrtol < 0.
*! For every accepted case the returned (status, Td, eps) is compared exactly:
*! Td is 0 in the band and T0 - T1 outside it, eps is c * nrtol.
*!
*! MUTANTS, each of which must turn this file red. Killed by (measured
*! 2026-10-09, StataMP 18; r[3] = eps is compared, so formula mutants fail
*! every accepted case):
*!   B1 "<" for "<=" at +eps                      2, 3, 13, 14
*!   B2 "<=" for "<" at -eps (abort at -eps)      3, 14
*!   B3 eps = (c + 1) * nrtol                     1-17 except 11
*!   B4 abort threshold 0 instead of -eps         3, 9, 14
*!   B5 band returns Td, not 0                    2, 3, 9, 10, 13, 14, 17
*!   B6 nrtol hard-coded to 1e-12                 10, 17
*!   B7 band on one side only (Td in [0, eps])    3, 9, 14
*!   B8a no missing-input check                   18, 19 (20 and 24 are
*!       rc-equivalent: eps is missing and the band abort fires, 3498 with
*!       the wrong message; only the message distinguishes them)
*!   B8b no c > 0 check                           21, 22
*!   B8c no nrtol > 0 check                       23, 25
*!   B9 status codes swapped                      every accepted case
*!   B10 eps = (c - 1) * nrtol (the old formula   every accepted case
*!       applied to c)
*!   B11 c required to be an integer              17
*!   B12 c > 0 check as c >= 0 (c = 0 accepted)   21
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

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
        display as error "test_dtband: `f' failed to load"
        exit 9
    }
}

mata:

// columns: T0, T1, c, nrtol, wanted rc, wanted status
real matrix tdd_cases()
{
    real scalar c1, nr, eps, c2, eps2

    c1   = 499
    nr   = 1e-12
    eps  = c1 * nr
    c2   = 998
    eps2 = c2 * nr
    return(( 10,                 4,                   c1,  nr,    0,    1) \
           ( eps,                0,                   c1,  nr,    0,    0) \
           ( 0,                  eps,                 c1,  nr,    0,    0) \
           ( eps * (1 + 1e-9),   0,                   c1,  nr,    0,    1) \
           ( 0,                  eps * (1 + 1e-9),    c1,  nr, 3498,    .) \
           ( 1,                  1,                   c1,  nr,    0,    0) \
           ( 4.995e-10,          0,                   c1,  nr,    0,    1) \
           ( 0,                  4.995e-10,           c1,  nr, 3498,    .) \
           ( 0,                  eps / 2,             c1,  nr,    0,    0) \
           ( 1e-4,               0,                   c1, 1e-6,   0,    0) \
           ( 1.2,                4.2,                 c1,  nr, 3498,    .) \
           ( 506.9,              0.5,                 c1,  nr,    0,    1) \
           ( eps2,               0,                  c2,  nr,    0,    0) \
           ( 0,                  eps2,               c2,  nr,    0,    0) \
           ( 9.985e-10,          0,                  c2,  nr,    0,    1) \
           ( 0,                  9.985e-10,          c2,  nr, 3498,    .) \
           ( 2.4e-6,             0,                 2.5, 1e-6,   0,    0) \
           ( .,                  1,                   c1,  nr, 3498,    .) \
           ( 1,                  .,                   c1,  nr, 3498,    .) \
           ( 1,                  0,                   .,  nr, 3498,    .) \
           ( 1,                  0,                   0,  nr, 3498,    .) \
           ( 1,                  0,                  -5,  nr, 3498,    .) \
           ( 1,                  0,                   c1,   0, 3498,    .) \
           ( 1,                  0,                   c1,   ., 3498,    .) \
           ( 1,                  0,                   c1, -1e-12, 3498, .))
}

// Runs one case; sets the local ok to 1 when the returned row is exactly right.
void tdd_run(real rowvector c)
{
    real rowvector r
    real scalar wantu

    r = dwlsmv_dtband(c[1], c[2], c[3], c[4], "H0", "H1")
    if (c[6] == 0) {
        wantu = 0
    }
    else {
        wantu = c[1] - c[2]
    }
    st_local("ok", strofreal(r[1] == c[6] & r[2] == wantu & r[3] == c[3] * c[4]))
}

end

mata: CASE = tdd_cases()
mata: st_local("ncase", strofreal(rows(CASE)))

local nfail 0

* Guard against an empty or truncated case list passing vacuously.
if `ncase' != 25 {
    display as error "test_dtband: expected 25 cases, built `ncase'"
    exit 9
}

forvalues i = 1/`ncase' {
    mata: st_local("want", strofreal(CASE[`i', 5]))
    local ok 0
    capture noisily mata: tdd_run(CASE[`i', .])
    local rc = _rc
    if `rc' == `want' & (`want' != 0 | `ok' == 1) {
        display as text "    case `i' (rc `rc', as required)"
    }
    else {
        display as error "    case `i' FAILED: rc `rc', wanted `want', row correct `ok'"
        local nfail = `nfail' + 1
    }
}

if `nfail' {
    display as error _n "test_dtband: `nfail' failure(s)"
    exit 9
}

display as result "test_dtband: eps = c nrtol; Td < -eps aborts (rc 3498), |Td| <= eps is reported as 0, otherwise Td proceeds; bad inputs rc 3498"
exit 0
