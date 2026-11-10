*! Pins spec 4.8 O7 as a policy: theta(2) moves lambda, beta and psi, not nu, omega
*! or alpha. This is NOT a claim about identification. Model F is deficient at every
*! alpha = 0 point tried and full rank at every nonzero alpha tried (2026-11-01,
*! E-074); the refusal below describes the start. If the O7 policy review
*! changes the rule, this test changes with it.
*!
*!   F    alpha free, start 0: rc 3498, rank deficient. Mutant Mth2 accepts it.
*!   F1   the same with alpha free 0.5: accepted at theta(0), npar 7, values as written.
*!   B    two factors, two indicators each, beta free in both groups, no psi21: passes
*!        only because theta(2) moves beta. A mutant that leaves beta at 0 refuses it.
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
struct dwlsmv_parmap scalar top_id(string scalar al)
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st
    string rowvector a
    string scalar s
    real scalar g, i

    m  = dwlsmv_statmap(1, 2)
    st = dwlsmv_stage1(1, 2)
    for (g = 1; g <= 2; g++) {
        m[g]      = dwlsmv_statmap_build(J(2, 1, 3))
        st[g].map = m[g]
        st[g].s   = J(m[g].D, 1, 0.3)
        st[g].s[1..4] = (-0.4 \ 0.4 \ -0.3 \ 0.5)
    }
    a = ("groups 2", "group 1", "factor F", "ind y1 y2",
         "lambda F y1 fixed 0.7", "lambda F y2 fixed 0.7", "metric unitvariance F",
         "group 2", "factor F", "ind y1 y2", "lambda F y1 free", "lambda F y2 free",
         "psi F F fixed 1", "alpha F " + al,
         "tie tau 1 y1 1 tau 2 y1 1", "tie tau 1 y1 2 tau 2 y1 2",
         "tie tau 1 y2 1 tau 2 y2 1", "tie tau 1 y2 2 tau 2 y2 2")
    s = ""
    for (i = 1; i <= cols(a); i++) {
        s = s + a[i] + char(10)
    }
    return(dwls_identify(s, m, st))
}

void top_F()
{
    struct dwlsmv_parmap scalar p

    p = top_id("free")
}

void top_F1()
{
    struct dwlsmv_parmap scalar p
    real scalar i, bad

    p   = top_id("free 0.5")
    bad = (p.npar != 7)
    for (i = 1; i <= p.nall; i++) {
        if (p.ptype[i] == 7 & p.value[i] != 0.5) {
            bad = 1
        }
        if (p.ptype[i] == 1 & p.value[i] != 0.7) {
            bad = 1
        }
    }
    st_local("bad", strofreal(bad))
}

void top_B()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st
    struct dwlsmv_parmap scalar p
    string rowvector a, l
    string scalar s
    real scalar g, i, bad

    m  = dwlsmv_statmap(1, 2)
    st = dwlsmv_stage1(1, 2)
    for (g = 1; g <= 2; g++) {
        m[g]      = dwlsmv_statmap_build(J(4, 1, 3))
        st[g].map = m[g]
        st[g].s   = J(m[g].D, 1, 0.3)
        st[g].s[1..8] = (-0.4 \ 0.4 \ -0.3 \ 0.5 \ -0.5 \ 0.3 \ -0.2 \ 0.6)
    }
    l = ("factor F1 F2", "ind y1 y2 y3 y4", "lambda F1 y1 free", "lambda F1 y2 free",
         "lambda F2 y3 free", "lambda F2 y4 free", "beta F2 F1 free")
    a = ("groups 2", "group 1", l, "metric unitvariance F1", "metric unitvariance F2",
         "group 2", l, "psi F1 F1 fixed 1", "psi F2 F2 fixed 1")
    s = ""
    for (i = 1; i <= cols(a); i++) {
        s = s + a[i] + char(10)
    }
    p   = dwls_identify(s, m, st)
    bad = (p.npar != 26)
    for (i = 1; i <= p.nall; i++) {
        if (p.ptype[i] == 3 & p.value[i] == 0) {
            bad = 1
        }
    }
    st_local("bad", strofreal(bad))
}
end

local nfail 0
capture noisily mata: top_F()
display "F rc = " _rc
if _rc != 3498 {
    local ++nfail
}
local bad 1
capture noisily mata: top_F1()
display "F1 rc = " _rc " bad = `bad'"
if _rc | `bad' {
    local ++nfail
}
local bad 1
capture noisily mata: top_B()
display "B rc = " _rc " bad = `bad'"
if _rc | `bad' {
    local ++nfail
}
if `nfail' {
    display as error "test_identify_o7_policy: `nfail' failure(s)"
    exit 9
}
display "test_identify_o7_policy: pass"
