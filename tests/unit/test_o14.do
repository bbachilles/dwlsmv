*! test_o14.do -- spec 4.8 O14: the per-group representation at the public API.
*!
*! dwls_identify(), dwls_model() and dwls_delta() take struct dwlsmv_statmap
*! vector maps; dwls_fit() takes struct dwlsmv_stage1 vector; the fit's stage1
*! field is that vector. A single group is a 1-vector. The group code has landed
*! (Modules 3, 5, 6 and, 2026-10-10, 7-9), so o3-o5 now test the MISMATCH
*! refusals: two maps or Stage 1 results against a one-group table or fit.
*!
*! Checks:
*!   o1  a 1-vector and a scalar map give bitwise the same table, sigma, Delta
*!   o2  dwls_fit() with a 1-vector Stage 1: fit.stage1 has length 1 and holds it
*!   o3  two maps: dwls_identify, dwls_model, dwls_delta each _rc == 3498
*!   o4  two Stage 1 results, a one-group table: dwls_fit _rc == 3498
*!   o5  a one-group fit carrying two Stage 1 results: dwls_inference _rc == 3498
*!   o6  zero maps: dwls_identify _rc == 3498
*!
*! Mutants (each must turn the file red): drop any one length check (o3-o6);
*! test length(maps) > 1 instead of != 1 (o6); store only st, not a 1-vector,
*! in fit.stage1 (o2, which also fails to compile).
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
        display as error "test_o14: `f' failed to load"
        exit 9
    }
}

mata:

string scalar to_spec()
{
    return("factor F" + char(10) + "ind y1 y2 y3 y4" + char(10) +
           "lambda F y1 free" + char(10) + "lambda F y2 free" + char(10) +
           "lambda F y3 free" + char(10) + "lambda F y4 free" + char(10) +
           "metric unitvariance F")
}

struct dwlsmv_statmap scalar to_map()
{
    return(dwlsmv_statmap_build((3 \ 3 \ 2 \ 2)))
}

struct dwlsmv_stage1 scalar to_st()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    real colvector th

    st.map = to_map()
    par    = dwls_identify(to_spec(), st.map)
    th     = J(par.npar, 1, 0)
    th[(par.npar - 3)::par.npar] = (0.6 \ 0.7 \ 0.5 \ 0.8)
    st.s      = dwls_model(th, par, st.map) + 0.01 * (mod(1::st.map.D, 3) :- 1)
    st.Gamma  = I(st.map.D)
    st.N      = 500
    st.psibar = J(st.map.D, 1, 0)
    st.bound  = J(st.map.P, 1, 0)
    return(st)
}

// o1, o2. Returns the number of failures.
real scalar to_ok()
{
    struct dwlsmv_statmap scalar m
    struct dwlsmv_statmap vector mv
    struct dwlsmv_parmap scalar p1, pv
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_fit scalar f
    real colvector th
    real scalar nf

    nf = 0
    m  = to_map()
    mv = dwlsmv_statmap(1)
    mv[1] = m

    p1 = dwls_identify(to_spec(), m)
    pv = dwls_identify(to_spec(), mv)
    if (p1.npar != pv.npar | p1.nall != pv.nall | p1.ptype != pv.ptype |
        p1.row != pv.row | p1.col != pv.col | p1.free != pv.free |
        p1.value != pv.value | p1.fidx != pv.fidx | p1.name != pv.name) {
        printf("{err}o1: the table differs between a scalar map and a 1-vector\n")
        nf++
    }
    th = J(pv.npar, 1, 0.3)
    if (dwls_model(th, p1, m) != dwls_model(th, pv, mv) |
        dwls_delta(th, p1, m) != dwls_delta(th, pv, mv)) {
        printf("{err}o1: sigma or Delta differs between a scalar map and a 1-vector\n")
        nf++
    }

    st = to_st()
    sv = dwlsmv_stage1(1)
    sv[1] = st
    f = dwls_fit(sv, pv)
    if (length(f.stage1) != 1) {
        printf("{err}o2: fit.stage1 has length %f, wanted 1\n", length(f.stage1))
        nf++
    }
    else if (f.stage1[1].s != st.s | f.stage1[1].Gamma != st.Gamma |
             f.stage1[1].N != st.N | f.stage1[1].map.K != st.map.K) {
        printf("{err}o2: fit.stage1[1] is not the Stage 1 passed in\n")
        nf++
    }
    if (!f.converged) {
        printf("{err}o2: the fit did not converge\n")
        nf++
    }
    return(nf)
}

real scalar to_id2()
{
    struct dwlsmv_statmap vector mv
    mv = dwlsmv_statmap(2)
    mv[1] = to_map()
    mv[2] = to_map()
    (void) dwls_identify(to_spec(), mv)
    return(0)
}

real scalar to_id0()
{
    struct dwlsmv_statmap vector mv
    mv = dwlsmv_statmap(0)
    (void) dwls_identify(to_spec(), mv)
    return(0)
}

real scalar to_model2(real scalar delta)
{
    struct dwlsmv_statmap vector mv
    struct dwlsmv_parmap scalar p
    real colvector th
    mv = dwlsmv_statmap(2)
    mv[1] = to_map()
    mv[2] = to_map()
    p  = dwls_identify(to_spec(), to_map())
    th = J(p.npar, 1, 0.3)
    if (delta) (void) dwls_delta(th, p, mv)
    else       (void) dwls_model(th, p, mv)
    return(0)
}

real scalar to_fit2()
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar p
    sv = dwlsmv_stage1(2)
    sv[1] = to_st()
    sv[2] = to_st()
    p = dwls_identify(to_spec(), to_map())
    (void) dwls_fit(sv, p)
    return(0)
}

real scalar to_inf2()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar p
    struct dwlsmv_fit scalar f
    st = to_st()
    p  = dwls_identify(to_spec(), st.map)
    f  = dwls_fit(st, p)
    sv = dwlsmv_stage1(2)
    sv[1] = st
    sv[2] = st
    f.stage1 = sv
    dwls_inference(f)
    return(0)
}

end

local nfail 0

mata: st_local("nok", strofreal(to_ok()))
if `nok' {
    display as error "test_o14: o1/o2 `nok' failure(s)"
    local nfail = `nfail' + `nok'
}
else display as text "  o1, o2 ok: 1-vector = scalar bitwise; fit.stage1 is the 1-vector"

local c1 to_id2()
local l1 o3 dwls_identify, two maps
local c2 to_model2(0)
local l2 o3 dwls_model, two maps
local c3 to_model2(1)
local l3 o3 dwls_delta, two maps
local c4 to_fit2()
local l4 o4 dwls_fit, two Stage 1 results and a one-group table
local c5 to_inf2()
local l5 o5 dwls_inference, a one-group fit carrying two Stage 1 results
local c6 to_id0()
local l6 o6 dwls_identify, zero maps
forvalues i = 1/6 {
    capture mata: (void) `c`i''
    local rc = _rc
    display as text "  `l`i'': rc `rc'"
    if `rc' != 3498 {
        display as error "test_o14: `l`i'' gave rc `rc', wanted 3498"
        local nfail = `nfail' + 1
    }
}

if `nfail' {
    display as error _n "test_o14: `nfail' failure(s)"
    exit 9
}

display as result "test_o14: one group is a 1-vector at every public entry point, bitwise the old scalar call; a group count that does not match the table or fit is refused with rc 3498"
exit 0
