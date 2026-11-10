*! Module 12, check 8 guard: DIFFTEST is refused for a pair with groups and a
*! structural path (user, 2026-10-10, spec 4.8 O2, option (e); the beta DIFFTEST
*! follow-up was NOT MET, E-004).
*!
*! Written and executed 2026-10-10 by the main session.
*!
*! The guard reads only grp, ptype, free and value, so the tables are built by
*! hand (no dwls_identify(), which needs per-group maps and Stage 1 results).
*!
*! Helper dwlsmv_dtn_mgbeta(par), expected:
*!   h1 one group (grp all 1), beta free                         0
*!   h2 two groups, no beta row                                  0
*!   h3 two groups, beta free                                    1
*!   h4 two groups, beta fixed at 0 only (an absent path)        0
*!   h5 two groups, beta fixed at .4                             1
*!   h6 grp empty (a hand-built single-group table), beta free   0  (max() of an
*!      empty grp is missing, and missing > 1 is true: the trap the guard avoids)
*! dwlsmv_dtnest():
*!   d1 both tables two-group with beta: rc 3498 and the guard's message
*!   d2 only H1 two-group with beta: rc 3498 and the guard's message
*!   d3 control: a single-group pair with beta in both passes the guard (rc 0;
*!      a valid nested pair, d = 1)
*!
*! Mutants (each must turn this test red):
*!   G1 guard removed                                    d1, d2
*!   G2 guard tests par0 only                            d2
*!   G3 groups found by max(grp) > 1                     h6
*!   G4 a beta fixed at 0 counts                         h4
*!   G5 fixed nonzero beta does not count                h5
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

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
        display as error "test_dtnest_mgbeta: `f' failed to load"
        exit 9
    }
}

mata:
// A table: rows (ptype, grp, row, col, free, value); grp may be J(0, 1, .).
struct dwlsmv_parmap scalar tdm_par(real matrix R, real scalar emptygrp)
{
    struct dwlsmv_parmap scalar p
    real scalar i, k

    p.nall  = rows(R)
    p.ptype = R[., 1]
    p.grp   = (emptygrp ? J(0, 1, .) : R[., 2])
    p.row   = R[., 3]
    p.col   = R[., 4]
    p.free  = R[., 5]
    p.value = R[., 6]
    p.name  = J(p.nall, 1, "")
    p.fidx  = J(p.nall, 1, 0)
    k = 0
    for (i = 1; i <= p.nall; i++) {
        p.name[i] = "g" + strofreal(R[i, 2]) + ":t" + strofreal(R[i, 1]) + "[" + strofreal(R[i, 3]) +
                    "," + strofreal(R[i, 4]) + "]"
        if (R[i, 5]) {
            k++
            p.fidx[i] = k
        }
    }
    p.npar = k
    return(p)
}

// lambda rows for two factors in each listed group (3 items each), free
real matrix tdm_lam(real rowvector gs)
{
    real matrix L
    real scalar g, j

    L = J(0, 6, .)
    for (g = 1; g <= cols(gs); g++) {
        for (j = 1; j <= 6; j++) {
            L = L \ (1, gs[g], j, (j <= 3 ? 1 : 2), 1, .7)
        }
    }
    return(L)
}

real scalar tdm_helper(string scalar id)
{
    real matrix R
    real scalar e

    e = 0
    if (id == "h1") {
        R = tdm_lam(1) \ (3, 1, 2, 1, 1, 0)
    }
    else if (id == "h2") {
        R = tdm_lam((1, 2))
    }
    else if (id == "h3") {
        R = tdm_lam((1, 2)) \ (3, 1, 2, 1, 1, 0) \ (3, 2, 2, 1, 1, 0)
    }
    else if (id == "h4") {
        R = tdm_lam((1, 2)) \ (3, 1, 2, 1, 0, 0) \ (3, 2, 2, 1, 0, 0)
    }
    else if (id == "h5") {
        R = tdm_lam((1, 2)) \ (3, 1, 2, 1, 0, .4) \ (3, 2, 2, 1, 0, .4)
    }
    else {
        R = tdm_lam(1) \ (3, 1, 2, 1, 1, 0)
        e = 1
    }
    return(dwlsmv_dtn_mgbeta(tdm_par(R, e)))
}

void tdm_dtnest(string scalar id)
{
    struct dwlsmv_parmap scalar p0, p1
    real matrix B1, B2

    B1 = tdm_lam(1) \ (3, 1, 2, 1, 1, 0)
    B2 = tdm_lam((1, 2)) \ (3, 1, 2, 1, 1, 0) \ (3, 2, 2, 1, 1, 0)
    if (id == "d1") {
        p0 = tdm_par(B2 \ (2, 1, 2, 1, 0, 0), 0)
        p1 = tdm_par(B2 \ (2, 1, 2, 1, 1, 0), 0)
        (void) dwlsmv_dtnest(p0, p1, 10, 9, "H0", "H1")
    }
    else if (id == "d2") {
        p0 = tdm_par(B1, 0)
        p1 = tdm_par(B2, 0)
        (void) dwlsmv_dtnest(p0, p1, 10, 9, "H0", "H1")
    }
    else {
        p0 = tdm_par(B1, 0)
        p1 = tdm_par(B1 \ (2, 1, 2, 1, 1, 0), 0)
        st_numscalar("tdm_d", dwlsmv_dtnest(p0, p1, 10, 9, "H0", "H1"))
    }
}
end

local nfail 0
foreach c in h1:0 h2:0 h3:1 h4:0 h5:1 h6:0 {
    gettoken id want : c, parse(":")
    local want = subinstr("`want'", ":", "", 1)
    mata: st_local("got", strofreal(tdm_helper("`id'")))
    if "`got'" != "`want'" {
        display as error "    `id' FAILED: helper returned `got', wanted `want'"
        local ++nfail
    }
}

tempfile lg
foreach id in d1 d2 {
    capture log close tdm
    log using "`lg'", text replace name(tdm)
    capture noisily mata: tdm_dtnest("`id'")
    local rc = _rc
    log close tdm
    mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("lg"))'), "DIFFTEST with groups and a") > 0))
    if `rc' != 3498 | !`hit' {
        display as error "    `id' FAILED: rc `rc' (wanted 3498), guard message found `hit'"
        local ++nfail
    }
}

capture noisily mata: tdm_dtnest("d3")
if _rc {
    display as error "    d3 FAILED: the single-group control was refused, rc " _rc
    local ++nfail
}
else if scalar(tdm_d) != 1 {
    display as error "    d3 FAILED: d = " scalar(tdm_d) ", wanted 1"
    local ++nfail
}

if `nfail' {
    display as error "test_dtnest_mgbeta: `nfail' failure(s)"
    exit 9
}
display as result "test_dtnest_mgbeta: the check 8 guard refuses a pair with groups and a beta (either table), counts groups without max(), ignores a beta fixed at 0, and leaves a single-group beta pair alone"
exit 0
