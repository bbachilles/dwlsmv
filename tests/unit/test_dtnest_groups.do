*! Module 12, check 8 at G > 1 (user, 2026-10-12, option (b); spec Module 12
*! check 8): rows keyed by (ptype, grp, row, col) (spec 4.8 G6), absent rows read
*! by G4 (omega 1, every other type 0), and any pair with a tie refused, rc 3498.
*!
*! Written and executed 2026-10-12 by the main session. Tables are built by hand
*! (as in test_dtnest_mgbeta.do): check 8 reads only the parameter tables.
*! Every table has two groups, two factors, three loadings per factor, no beta.
*!
*! Cases (want):
*!   k1 configural pair: H0 fixes psi21 at 0 in group 2, H1 frees it; H1 also
*!      lists psi21 in group 1 fixed at 0, which H0 omits      accept, d = 1
*!   k2 H0 frees psi21 in group 2 only, H1 in group 1 only      3498 not nested
*!   w1 omega absent from H0, fixed at 1 in H1 (group 2)        accept, d = 1
*!   w2 omega absent from H0, fixed at .5 in H1                 3498 not nested
*!   w3 omega fixed at 1 in H0, absent from H1                  accept, d = 1
*!   t1 H0 ties lambda across groups, H1 untied                 3498 tie message
*!   t2 H0 untied, H1 ties lambda (the tie only in H1)          3498 tie message
*!   t3 H0 ties lambda, H1 ties tau (the non-nested pair)       3498 tie message
*!   e1 the k1 pair with grp empty in both, one group: control  accept, d = 1
*!
*! Mutants (each must turn this test red). Predicted | observed 2026-10-12, all red:
*!   K1 group dropped from the row key           k1         | k1 k2 w1 w2 w3
*!   K2 absent default 0 for omega (old rule)    w1, w3     | w1 w3
*!   K3 absent default 1 for every type          k1         | k1
*!   K4 tie refusal removed                      t1, t2, t3 | t1 t2 t3
*!   K5 tie refusal tests par0 only              t2         | t2
*!   K6 fixed rows (fidx 0) counted as ties      k1 w1 w3 e1| k1 k2 w1 w2 w3 e1
*! (K1 and K6 also turn k2/w2 red: they refuse with another message.)
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
        display as error "test_dtnest_groups: `f' failed to load"
        exit 9
    }
}

mata:
// A table from rows (ptype, grp, row, col, free, value, tie): free rows with the
// same tie key > 0 share one fidx. grp is emptied when emptygrp is set.
struct dwlsmv_parmap scalar tdg_par(real matrix R, real scalar emptygrp)
{
    struct dwlsmv_parmap scalar p
    real scalar i, k, m
    real colvector seen

    p.nall  = rows(R)
    p.ptype = R[., 1]
    p.grp   = (emptygrp ? J(0, 1, .) : R[., 2])
    p.row   = R[., 3]
    p.col   = R[., 4]
    p.free  = R[., 5]
    p.value = R[., 6]
    p.name  = J(p.nall, 1, "")
    p.fidx  = J(p.nall, 1, 0)
    seen    = J(100, 1, 0)
    k = 0
    for (i = 1; i <= p.nall; i++) {
        p.name[i] = "g" + strofreal(R[i, 2]) + ":t" + strofreal(R[i, 1]) + "[" + strofreal(R[i, 3]) +
                    "," + strofreal(R[i, 4]) + "]"
        if (R[i, 5]) {
            m = R[i, 7]
            if (m > 0) {
                if (seen[m] == 0) {
                    k++
                    seen[m] = k
                }
                p.fidx[i] = seen[m]
            }
            else {
                k++
                p.fidx[i] = k
            }
        }
    }
    p.npar = k
    return(p)
}

// Loadings (free, untied unless tie is set: item j tied on key j), variances
// fixed at 1, in every listed group.
real matrix tdg_base(real rowvector gs, real scalar tie)
{
    real matrix L
    real scalar g, j

    L = J(0, 7, .)
    for (g = 1; g <= cols(gs); g++) {
        for (j = 1; j <= 6; j++) {
            L = L \ (1, gs[g], j, (j <= 3 ? 1 : 2), 1, .7, (tie ? j : 0))
        }
        L = L \ (2, gs[g], 1, 1, 0, 1, 0) \ (2, gs[g], 2, 2, 0, 1, 0)
    }
    return(L)
}

// tau rows, item 1 one threshold per group, tied on key 50 when tie is set
real matrix tdg_tau(real scalar tie)
{
    return((4, 1, 1, 1, 1, 0, (tie ? 50 : 0)) \ (4, 2, 1, 1, 1, 0, (tie ? 50 : 0)))
}

real scalar tdg_run(string scalar id)
{
    real matrix B, R0, R1
    real scalar e

    B = tdg_base((1, 2), 0)
    e = 0
    if (id == "k1") {
        R0 = B \ (2, 2, 2, 1, 0, 0, 0)
        R1 = B \ (2, 1, 2, 1, 0, 0, 0) \ (2, 2, 2, 1, 1, 0, 0)
    }
    else if (id == "k2") {
        R0 = B \ (2, 2, 2, 1, 1, 0, 0)
        R1 = B \ (2, 1, 2, 1, 1, 0, 0)
    }
    else if (id == "w1") {
        R0 = B
        R1 = B \ (6, 2, 1, 1, 0, 1, 0) \ (2, 2, 2, 1, 1, 0, 0)
    }
    else if (id == "w2") {
        R0 = B
        R1 = B \ (6, 2, 1, 1, 0, .5, 0) \ (2, 2, 2, 1, 1, 0, 0)
    }
    else if (id == "w3") {
        R0 = B \ (6, 2, 1, 1, 0, 1, 0)
        R1 = B \ (2, 2, 2, 1, 1, 0, 0)
    }
    else if (id == "t1") {
        R0 = tdg_base((1, 2), 1)
        R1 = B
    }
    else if (id == "t2") {
        R0 = B \ (2, 2, 2, 1, 0, 0, 0)
        R1 = tdg_base((1, 2), 1) \ (2, 2, 2, 1, 1, 0, 0)
    }
    else if (id == "t3") {
        R0 = tdg_base((1, 2), 1) \ tdg_tau(0)
        R1 = B \ tdg_tau(1)
    }
    else {
        B  = tdg_base(1, 0)
        R0 = B \ (2, 1, 2, 1, 0, 0, 0)
        R1 = B \ (2, 1, 2, 1, 1, 0, 0)
        e  = 1
    }
    return(dwlsmv_dtnest(tdg_par(R0, e), tdg_par(R1, e), 10, 9, "H0", "H1"))
}
end

local nfail 0
tempfile lg

* accepted, d = 1
foreach id in k1 w1 w3 e1 {
    capture noisily mata: st_numscalar("tdg_d", tdg_run("`id'"))
    if _rc {
        display as error "    `id' FAILED: refused, rc " _rc " (wanted acceptance)"
        local ++nfail
    }
    else if scalar(tdg_d) != 1 {
        display as error "    `id' FAILED: d = " scalar(tdg_d) ", wanted 1"
        local ++nfail
    }
}

* refused, rc 3498, with the expected message
foreach c in "k2|the models are not nested" "w2|the models are not nested" ///
             "t1|equality constraints (tied" "t2|equality constraints (tied" ///
             "t3|equality constraints (tied" {
    gettoken id msg : c, parse("|")
    local msg = substr("`msg'", 2, .)
    capture log close tdg
    log using "`lg'", text replace name(tdg)
    capture noisily mata: (void) tdg_run("`id'")
    local rc = _rc
    log close tdg
    mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("lg"))'), st_local("msg")) > 0))
    if `rc' != 3498 | !`hit' {
        display as error "    `id' FAILED: rc `rc' (wanted 3498), message found `hit'"
        local ++nfail
    }
}

if `nfail' {
    display as error "test_dtnest_groups: `nfail' failure(s)"
    exit 9
}
display as result "test_dtnest_groups: check 8 keys rows by group, reads an absent omega as 1 and other absent rows as 0, and refuses any pair with a tie (either table)"
exit 0
