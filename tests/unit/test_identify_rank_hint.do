*! dwls_identify's rank refusal: the psi/beta sentence prints only when the model
*! frees a psi covariance and a beta path between one pair of factors in one
*! group. Case F (2026-11-01 Module 3 diagnostics) has neither and printed it.
*! The O7 sentence replaces the artefact claim when a free nu, omega or alpha sat at
*! its reference value (F); the psi/beta model keeps the artefact claim.
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
string scalar trh_join(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

string scalar trh_slurp(string scalar path)
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

void trh_psibeta()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar p

    m = dwlsmv_statmap(1, 1)
    m[1] = dwlsmv_statmap_build(J(6, 1, 3))
    p = dwls_identify(trh_join(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "metric unitvariance F1", "metric unitvariance F2",
        "beta F2 F1 free", "psi F2 F1 free")), m)
}

void trh_F()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st
    struct dwlsmv_parmap scalar p
    real scalar g

    m  = dwlsmv_statmap(1, 2)
    st = dwlsmv_stage1(1, 2)
    for (g = 1; g <= 2; g++) {
        m[g]      = dwlsmv_statmap_build(J(2, 1, 3))
        st[g].map = m[g]
        st[g].s   = J(m[g].D, 1, 0.3)
    }
    p = dwls_identify(trh_join(("groups 2", "group 1", "factor F", "ind y1 y2",
        "lambda F y1 fixed 0.7", "lambda F y2 fixed 0.7", "metric unitvariance F",
        "group 2", "factor F", "ind y1 y2", "lambda F y1 free", "lambda F y2 free",
        "psi F F fixed 1", "alpha F free",
        "tie tau 1 y1 1 tau 2 y1 1", "tie tau 1 y1 2 tau 2 y1 2",
        "tie tau 1 y2 1 tau 2 y2 1", "tie tau 1 y2 2 tau 2 y2 2")), m, st)
}
end

local nfail 0
foreach c in psibeta F {
    tempfile lg
    quietly log using "`lg'", text replace name(trh)
    capture noisily mata: trh_`c'()
    local rc = _rc
    quietly log close trh
    mata: st_local("has", strofreal(strpos(trh_slurp(st_local("lg")), "Freeing a psi covariance") > 0))
    mata: st_local("def", strofreal(strpos(trh_slurp(st_local("lg")), "rank deficient") > 0))
    local want = ("`c'" == "psibeta")
    mata: st_local("o7", strofreal(strpos(trh_slurp(st_local("lg")), "stayed at") > 0))
    mata: st_local("art", strofreal(strpos(trh_slurp(st_local("lg")), "an artefact of") > 0))
    display "`c': reference-value text `o7' artefact text `art'"
    if `o7' == `want' | `art' != `want' {
        local ++nfail
    }
    display "`c': rc `rc' rank-deficient text `def' hint `has' (want `want')"
    if `rc' != 3498 | `def' != 1 | `has' != `want' {
        local ++nfail
    }
}
if `nfail' {
    display as error "test_identify_rank_hint: `nfail' failure(s)"
    exit 9
}
display "test_identify_rank_hint: pass"
