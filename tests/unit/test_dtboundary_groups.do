*! Module 12, check 11 at G > 1 (derivation @sec-inv-check11 C11-1 to C11-3;
*! spec Module 12 check 11; exemption gate step 1, 2026-10-13) and the
*! configural preflight in Module 11's invariance() flow (call site decided
*! 2026-10-13, user, P1 (i)).
*!
*! Executed under StataMP 18 on 2026-10-13 and it passes (15 checks; written
*! and run by the main session). run_all.do discovers it by name.
*!
*! Part 1, dwlsmv_dtboundary() on hand-built tables (check 11 reads only the
*! table). Base: two groups, two factors, items 1-3 on F1 and 4-6 on F2, every
*! loading free, psi11 and psi22 fixed at 1 in each group, psi21 absent (fixed
*! 0). Cases (want):
*!   c0 the base                                                      accept
*!   a1 group 1 psi21 free, group 2 psi21 fixed 1                     3498, group 2
*!   a2 a1 with group 2 psi11 free (C11-3: free diagonal)             accept
*!   a3 group 2 has no psi22 row (absent = fixed 0)                   3498
*!   b1 group 2 item 1 loading fixed 1, omega absent (reads 1)        3498, rule (b)
*!   b2 group 2 item 1 loading fixed .9, omega fixed .9 (Theta .3439) accept
*!   b3 group 2 item 1 loading fixed .9, omega fixed 1.2 (-.1664)     3498, group 2
*!   b4 b3 with omega free, start 1.2                                 accept
*!   b5 group 1 item 1 loading fixed .9 (Theta .19), group 2 item 1
*!      loading free and omega fixed 1.2 (omega is group-local)       accept
*!   b6 item 1 loading tied across groups, start 1.05                 accept
*!   e1 one group, grp empty, psi21 fixed 1                           3498
*!   e2 one group, grp empty, the base                                accept
*! The G = 1 rule is covered unchanged by test_dtboundary.do (20 cases).
*!
*! Part 2, the preflight through the command (fixtures/data/mg_E.csv,
*! (F1 -> y1 y2 y3) (F2 -> y4 y5), group(g), metric(unitvariance)):
*!   p1 cov(F1*F2@1) at configural and at metric: rc 3498 with check 11's
*!      message naming the configural specification and group 1 (C11-4:
*!      group 1 refuses first), nothing posted
*!   p2 cov(F1*F2@.5) at configural: rc 0, the fit is posted
*!
*! Mutants (each must turn this test red). Predicted | observed 2026-10-13
*! (scratch worktree, mutated line printed before each run), all red:
*!   B1 only group 1 checked                     a1 a3 b1 b3 | a1 a3 b1 b3
*!   B2 absent omega reads 0                     b1          | b1
*!   B3 omega ignored (m = 1)                    b3          | b3
*!   B4 omega linear, not squared                b3          | b3
*!   B5 free omega read at its start             b4          | b4
*!   B6 omega looked up ignoring the group       b5          | b5
*!   B7 grp fallback removed                     e1 e2       | e1 e2 (rc 3301)
*!   B8 group loop stops at G - 1                a1          | a1 a3 b1 b3 e1
*!   P  preflight call removed from fit_g        p1          | p1
*! Step 2a (2026-10-13, check 11 messages): p1 also asserts the preflight's
*! consequence sentence; b1 and b3 assert rule (b) names item y1 (tdbg_par
*! names loadings as dwls_identify does). Mutants, predicted | observed:
*!   D1 the routine drops the caller's sentence  p1          | p1 (both levels)
*!   D3 rule (b) names the row index, not item   b1 b3       | b1 b3
*! Under P the cov(F1*F2@1) model fits and is posted at rc 0 with no word:
*! the preflight is the only thing between that model and a fit.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

local nfail 0

* run a command under a log; rc into tdg_rc, the log text into Mata TDGLOG
capture program drop tdg_log
program tdg_log
    args cmd
    ereturn clear
    tempfile lg
    capture log close tdg
    quietly log using "`lg'", text replace name(tdg)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tdg
    mata: TDGLOG = invtokens(cat(st_local("lg"))')
    c_local tdg_rc `rc'
end

* 1 if the text is in the last log
capture program drop tdg_has
program tdg_has
    args s
    mata: st_local("h", strofreal(strpos(TDGLOG, st_local("s")) > 0))
    c_local hit `h'
end

mata:
// A table from rows (ptype, grp, row, col, free, value, tie): free rows with
// the same tie key > 0 share one fidx. grp is emptied when emptygrp is set.
struct dwlsmv_parmap scalar tdbg_par(real matrix R, real scalar emptygrp)
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
        if (R[i, 1] == 1) {
            // dwls_identify's form, so rule (b)'s message can name the item
            p.name[i] = "g" + strofreal(R[i, 2]) + ":lambda[y" + strofreal(R[i, 3]) +
                        ",F" + strofreal(R[i, 4]) + "]"
        }
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

// One group's base rows: loadings free (item 1 skipped when skip1), psi11
// and psi22 fixed at 1 unless dropped (d11, d22).
real matrix tdbg_grp(real scalar g, real scalar skip1, real scalar d11, real scalar d22)
{
    real matrix L
    real scalar j

    L = J(0, 7, .)
    for (j = 1; j <= 6; j++) {
        if (j == 1 & skip1) {
            continue
        }
        L = L \ (1, g, j, (j <= 3 ? 1 : 2), 1, .7, 0)
    }
    if (!d11) {
        L = L \ (2, g, 1, 1, 0, 1, 0)
    }
    if (!d22) {
        L = L \ (2, g, 2, 2, 0, 1, 0)
    }
    return(L)
}

void tdbg_run(string scalar id)
{
    real matrix R
    real scalar eg

    eg = 0
    if (id == "c0") {
        R = tdbg_grp(1, 0, 0, 0) \ tdbg_grp(2, 0, 0, 0)
    }
    else if (id == "a1") {
        R = tdbg_grp(1, 0, 0, 0) \ (2, 1, 2, 1, 1, 0, 0) \
            tdbg_grp(2, 0, 0, 0) \ (2, 2, 2, 1, 0, 1, 0)
    }
    else if (id == "a2") {
        R = tdbg_grp(1, 0, 0, 0) \ (2, 1, 2, 1, 1, 0, 0) \
            tdbg_grp(2, 0, 1, 0) \ (2, 2, 1, 1, 1, 1, 0) \ (2, 2, 2, 1, 0, 1, 0)
    }
    else if (id == "a3") {
        R = tdbg_grp(1, 0, 0, 0) \ tdbg_grp(2, 0, 0, 1)
    }
    else if (id == "b1") {
        R = tdbg_grp(1, 0, 0, 0) \ tdbg_grp(2, 1, 0, 0) \ (1, 2, 1, 1, 0, 1, 0)
    }
    else if (id == "b2") {
        R = tdbg_grp(1, 0, 0, 0) \ tdbg_grp(2, 1, 0, 0) \ (1, 2, 1, 1, 0, .9, 0) \
            (6, 2, 1, 1, 0, .9, 0)
    }
    else if (id == "b3") {
        R = tdbg_grp(1, 0, 0, 0) \ tdbg_grp(2, 1, 0, 0) \ (1, 2, 1, 1, 0, .9, 0) \
            (6, 2, 1, 1, 0, 1.2, 0)
    }
    else if (id == "b4") {
        R = tdbg_grp(1, 0, 0, 0) \ tdbg_grp(2, 1, 0, 0) \ (1, 2, 1, 1, 0, .9, 0) \
            (6, 2, 1, 1, 1, 1.2, 0)
    }
    else if (id == "b5") {
        R = tdbg_grp(1, 1, 0, 0) \ (1, 1, 1, 1, 0, .9, 0) \
            tdbg_grp(2, 0, 0, 0) \ (6, 2, 1, 1, 0, 1.2, 0)
    }
    else if (id == "b6") {
        R = tdbg_grp(1, 1, 0, 0) \ (1, 1, 1, 1, 1, 1.05, 7) \
            tdbg_grp(2, 1, 0, 0) \ (1, 2, 1, 1, 1, 1.05, 7)
    }
    else if (id == "e1") {
        R  = tdbg_grp(1, 0, 0, 0) \ (2, 1, 2, 1, 0, 1, 0)
        eg = 1
    }
    else if (id == "e2") {
        R  = tdbg_grp(1, 0, 0, 0)
        eg = 1
    }
    else {
        _error("tdbg_run: unknown case")
    }
    dwlsmv_dtboundary(tdbg_par(R, eg), "case " + id, "TDBG consequence.")
}
end

* ---------------------------------------------------------------- part 1
local ids  c0 a1 a2 a3 b1 b2 b3 b4 b5 b6 e1 e2
local want 0 3498 0 3498 3498 0 3498 0 0 0 3498 0
local ncase : word count `ids'
local nwant : word count `want'
if `ncase' != 12 | `nwant' != 12 {
    display as error "test_dtboundary_groups: expected 12 cases, built `ncase'/`nwant'"
    exit 9
}
forvalues i = 1/`ncase' {
    local id : word `i' of `ids'
    local w  : word `i' of `want'
    tdg_log `"mata: tdbg_run("`id'")"'
    if `tdg_rc' != `w' {
        display as error "    `id' FAILED: rc `tdg_rc', wanted `w'"
        local ++nfail
    }
    else {
        display as text "    `id' ok (rc `tdg_rc')"
    }
    * the refusing group and rule, read from the message
    if inlist("`id'", "a1", "b3") {
        tdg_has "case `id', group 2)"
        if !`hit' {
            display as error "    `id' FAILED: message does not name group 2"
            local ++nfail
        }
    }
    if "`id'" == "a1" {
        tdg_has "the fixed values in Psi force a boundary"
        if !`hit' {
            display as error "    a1 FAILED: not refused by rule (a)"
            local ++nfail
        }
    }
    if inlist("`id'", "b1", "b3") {
        tdg_has "force the residual variance of y1"
        if !`hit' {
            display as error "    `id' FAILED: not refused by rule (b) naming item y1"
            local ++nfail
        }
    }
    if inlist("`id'", "e1") {
        tdg_has ", group "
        if `hit' {
            display as error "    e1 FAILED: a one-group message names a group"
            local ++nfail
        }
    }
}

* ---------------------------------------------------------------- part 2
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
local mdl (F1 -> y1 y2 y3) (F2 -> y4 y5), group(g) metric(unitvariance)

foreach lev in configural metric {
    capture mata: dwlsmv_m11_clear()
    tdg_log `"dwlsmv `mdl' invariance(`lev') cov(F1*F2@1)"'
    if `tdg_rc' != 3498 {
        display as error "    p1 `lev' FAILED: rc `tdg_rc', wanted 3498"
        local ++nfail
    }
    tdg_has "check 11 (invariance() preflight on the configural specification, group"
    local h1 `hit'
    tdg_has "g1:psi[F2,F1] = 1"
    if !`h1' | !`hit' {
        display as error "    p1 `lev' FAILED: not check 11's preflight message on group 1"
        local ++nfail
    }
    * the preflight caller's consequence (spec check 11, approved 2026-10-13)
    tdg_has "The configural specification is refused before fitting: no invariance level can be tested from it."
    if !`hit' {
        display as error "    p1 `lev' FAILED: the preflight's consequence sentence is missing"
        local ++nfail
    }
    if "`e(cmd)'" != "" {
        display as error "    p1 `lev' FAILED: something was posted"
        local ++nfail
    }
    display as text "    p1 `lev' done (rc `tdg_rc')"
}

capture mata: dwlsmv_m11_clear()
tdg_log `"dwlsmv `mdl' invariance(configural) cov(F1*F2@.5)"'
if `tdg_rc' != 0 | "`e(cmd)'" != "dwlsmv" {
    display as error "    p2 FAILED: rc `tdg_rc', e(cmd) `e(cmd)'"
    local ++nfail
}
else {
    display as text "    p2 ok"
}

if `nfail' {
    display as error "test_dtboundary_groups: `nfail' failure(s)"
    exit 9
}
display as result "test_dtboundary_groups: all passed"
