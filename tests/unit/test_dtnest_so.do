*! Module 12, check 8, release B (B3.9): DIFFTEST between a model with a
*! second-order factor and one without is refused with its own message, rc 3498.
*!
*! Written and executed 2026-11-20 by the main session. Tables are built by
*! hand, as in test_dtnest_mgbeta.do (the guard reads ptype, grp, col, free and
*! value only). Factors: 1, 2, 3 measured by items (3 each); 4 = G.
*!
*! Helper dwlsmv_dtn_so(par), expected:
*!   h1 correlated three-factor model (psi rows only)            0
*!   h2 G over 1-3, free betas, no lambda on G                   1
*!   h3 G's only beta rows fixed at 0                            0
*!   h4 G with a beta fixed at 1 (nonzero fixed)                 1
*!   h5 a structural model F1 -> F2 (F1 has items)               0
*!   h6 G with a lambda row fixed at 0 only                      1
*!   h7 G with a free lambda row (an item-anchored parent)       0
*! dwlsmv_dtnest():
*!   d1 H0 second-order, H1 its correlated parent: rc 3498, B3.9 message
*!   d2 the reverse order: rc 3498, B3.9 message
*!   d3 control: two second-order models, H0 with one beta fixed at .7 and H1
*!      free: rc 0, d = 1
*!
*! Mutants (each must turn this test red):
*!   S1 guard removed                                  d1 (rule 0 refuses with
*!                                                     another message), d2
*!   S2 a beta fixed at 0 counts as a path             h3
*!   S3 a lambda fixed at 0 counts as an indicator     h6
*!   S4 the guard reads par0 only (so(par0) == 1)      d2
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
        display as error "test_dtnest_so: `f' failed to load"
        exit 9
    }
}

mata:
// A table: rows (ptype, grp, row, col, free, value); names carry factor names
// so rule 0 compares like with like.
struct dwlsmv_parmap scalar tds_par(real matrix R)
{
    struct dwlsmv_parmap scalar p
    real scalar i, k

    p.nall  = rows(R)
    p.ptype = R[., 1]
    p.grp   = R[., 2]
    p.row   = R[., 3]
    p.col   = R[., 4]
    p.free  = R[., 5]
    p.value = R[., 6]
    p.name  = J(p.nall, 1, "")
    p.fidx  = J(p.nall, 1, 0)
    k = 0
    for (i = 1; i <= p.nall; i++) {
        p.name[i] = "t" + strofreal(R[i, 1]) + "[" + strofreal(R[i, 3]) + "," + strofreal(R[i, 4]) + "]"
        if (R[i, 5]) {
            k++
            p.fidx[i] = k
        }
    }
    p.npar = k
    return(p)
}

// lambda rows: items 1-9 on factors 1-3, the first of each fixed at 1
real matrix tds_lam()
{
    real matrix L
    real scalar j, f

    L = J(0, 6, .)
    for (j = 1; j <= 9; j++) {
        f = ceil(j / 3)
        L = L \ (1, 1, j, f, (mod(j, 3) != 1), (mod(j, 3) == 1 ? 1 : .7))
    }
    return(L)
}

// the correlated parent: free psi among factors 1-3
real matrix tds_corr()
{
    return(tds_lam() \ (2, 1, 1, 1, 1, .5) \ (2, 1, 2, 1, 1, .3) \ (2, 1, 2, 2, 1, .5) \
           (2, 1, 3, 1, 1, .3) \ (2, 1, 3, 2, 1, .3) \ (2, 1, 3, 3, 1, .5))
}

// the second-order model: G = factor 4 over 1-3 (betas row = to, col = from)
real matrix tds_so(real scalar b1free, real scalar b1val)
{
    return(tds_lam() \ (3, 1, 1, 4, b1free, b1val) \ (3, 1, 2, 4, 1, .7) \ (3, 1, 3, 4, 1, .7) \
           (2, 1, 1, 1, 1, .5) \ (2, 1, 2, 2, 1, .5) \ (2, 1, 3, 3, 1, .5) \ (2, 1, 4, 4, 0, 1))
}

real scalar tds_helper(string scalar id)
{
    real matrix R

    if (id == "h1") {
        R = tds_corr()
    }
    else if (id == "h2") {
        R = tds_so(1, .7)
    }
    else if (id == "h3") {
        R = tds_lam() \ (3, 1, 1, 4, 0, 0) \ (3, 1, 2, 4, 0, 0) \ (2, 1, 4, 4, 0, 1)
    }
    else if (id == "h4") {
        R = tds_so(0, 1)
    }
    else if (id == "h5") {
        R = tds_lam() \ (3, 1, 2, 1, 1, .4)
    }
    else if (id == "h6") {
        R = tds_so(1, .7) \ (1, 1, 1, 4, 0, 0)
    }
    else {
        R = tds_so(1, .7) \ (1, 1, 10, 4, 1, .7)
    }
    return(dwlsmv_dtn_so(tds_par(R)))
}

void tds_dtnest(string scalar id)
{
    struct dwlsmv_parmap scalar ps, pc, p0, p1

    ps = tds_par(tds_so(1, .7))
    pc = tds_par(tds_corr())
    if (id == "d1") {
        (void) dwlsmv_dtnest(ps, pc, 25, 24, "H0", "H1")
    }
    else if (id == "d2") {
        (void) dwlsmv_dtnest(pc, ps, 24, 25, "H0", "H1")
    }
    else {
        p0 = tds_par(tds_so(0, .7))
        p1 = ps
        st_numscalar("tds_d", dwlsmv_dtnest(p0, p1, 26, 25, "H0", "H1"))
    }
}
end

local nfail 0
foreach c in h1:0 h2:1 h3:0 h4:1 h5:0 h6:1 h7:0 {
    gettoken id want : c, parse(":")
    local want = subinstr("`want'", ":", "", 1)
    mata: st_local("got", strofreal(tds_helper("`id'")))
    if "`got'" != "`want'" {
        display as error "    `id' FAILED: helper returned `got', wanted `want'"
        local ++nfail
    }
    else display as text "    `id' ok"
}

tempfile lg
foreach id in d1 d2 {
    capture log close tds
    log using "`lg'", text replace name(tds)
    capture noisily mata: tds_dtnest("`id'")
    local rc = _rc
    log close tds
    mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("lg"))'), "one model has a second-order factor") > 0))
    if `rc' != 3498 | !`hit' {
        display as error "    `id' FAILED: rc `rc' (wanted 3498), B3.9 message found `hit'"
        local ++nfail
    }
    else display as text "    `id' ok"
}

capture noisily mata: tds_dtnest("d3")
if _rc {
    display as error "    d3 FAILED: the pair of second-order models was refused, rc " _rc
    local ++nfail
}
else if scalar(tds_d) != 1 {
    display as error "    d3 FAILED: d = " scalar(tds_d) ", wanted 1"
    local ++nfail
}
else display as text "    d3 ok"

if `nfail' {
    display as error "test_dtnest_so: `nfail' failure(s)"
    exit 9
}
display as result "test_dtnest_so: check 8 refuses a pair with a second-order factor on one side only (either order) with B3.9's message, ignores paths and loadings fixed at 0, and passes a pair of second-order models to the rules"
exit 0
