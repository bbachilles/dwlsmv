*! test_mg_model.do -- Modules 5 and 6 at G > 1 (spec 4.8 O14 contract C2;
*! derivation @sec-mg-sigma-delta), from hand-built tables: no grammar.
*!
*! The reference values are embedded in tests/unit/mg_model_ref.mata, written
*! by tools/export_mg2_ref.R from the derivation's mg2_sigma and mg2_delta
*! (themselves asserted there against central differences and lavaan 0.7.2).
*! Cases: S (two factors, two groups, every row type 1-7, two lambda ties, a
*! beta, K differing between groups) and design E at its four levels.
*!
*! Checks (tolerances fixed BEFORE the first run):
*!   m1  sigma = mg2_sigma, all five cases, max abs <= 1e-13 (two codings of
*!       one closed form on O(1) values)
*!   m2  Delta = mg2_delta, all five cases, max abs <= 1e-9 (lambda, beta and
*!       psi columns are dwlsmv_delta_cols' central differences, h = 1e-5:
*!       rounding ~ eps/h; Lemma MS2.1, no truncation error)
*!   m3  E configural: Delta on the threshold rows (a selection, G6) and the
*!       tau columns is EXACTLY the identity, and exactly 0 on every other
*!       column of those rows
*!   m4  E configural: each group's block of sigma and Delta is bitwise the
*!       single-group dwls_model/dwls_delta on that group's rows (MS4)
*!   m5  S with its table rows reversed: sigma and Delta bitwise unchanged
*!       (Modules 5/6 do not depend on row order)
*!   m6  refusals, each _rc == 3498: group key shorter than the table; a key
*!       outside 1..G; three maps for a two-group table; a nu row at G = 1
*!       (dwls_model and dwls_delta); no map at all
*!
*! Mutants, run 2026-10-10: eight killed, one equivalent.
*!   X1 rho not scaled by omega_i omega_j (mg_sigma)        m1
*!   X2 nu enters thresholds with + (mg_sigma)              m1
*!   X3 absent omega read as 0 (mg_parts)                   m1
*!   X4 alpha column without L (mg_dcols: -omega_j)          m2
*!   X5 tied columns assigned, not summed (dwls_delta)      SURVIVES, EQUIVALENT:
*!      each group writes only its own block of rows, and a tie has at most
*!      one member per group (G7 refuses two), so assignment and summation
*!      write the same numbers. It would matter only if one group could hold
*!      two rows with one fidx, which Module 3 refuses.
*!   X6 group offset not advanced (dwls_delta)              m2
*!   X7 omega column's pair entry uses V_ii (mg_dcols)      m2
*!   X8 tau column omega_i replaced by 1 (mg_dcols)         m2
*!   X9 mg_check skipped (dwls_model)                       m6
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
local load `load' tests/unit/mg_model_ref.mata
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_mg_model: `f' failed to load"
        exit 9
    }
}

mata:

struct dwlsmv_statmap vector tm_maps(real matrix K)
{
    struct dwlsmv_statmap vector m
    real scalar g
    m = dwlsmv_statmap(cols(K))
    for (g = 1; g <= cols(K); g++) {
        m[g] = dwlsmv_statmap_build(K[., g])
    }
    return(m)
}

void tm_case(string scalar nm, struct dwlsmv_parmap scalar par, real colvector th,
             real matrix K, real colvector sig, real matrix Dl)
{
    if (nm == "S")                    mgref_S(par, th, K, sig, Dl)
    else if (nm == "E_configural")    mgref_E_configural(par, th, K, sig, Dl)
    else if (nm == "E_thresholds")    mgref_E_thresholds(par, th, K, sig, Dl)
    else if (nm == "E_loadings")      mgref_E_loadings(par, th, K, sig, Dl)
    else if (nm == "E_intercepts")    mgref_E_intercepts(par, th, K, sig, Dl)
    else _error("tm_case: unknown case")
}

// m1, m2. Returns the number of failures.
real scalar tm_m12()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    string rowvector cs
    real colvector th, sig, s
    real matrix K, Dl, D
    real scalar c, nf, e1, e2

    nf = 0
    cs = ("S", "E_configural", "E_thresholds", "E_loadings", "E_intercepts")
    for (c = 1; c <= cols(cs); c++) {
        tm_case(cs[c], par, th, K, sig, Dl)
        maps = tm_maps(K)
        s = dwls_model(th, par, maps)
        D = dwls_delta(th, par, maps)
        e1 = (rows(s) == rows(sig) ? max(abs(s - sig)) : .)
        e2 = (rows(D) == rows(Dl) & cols(D) == cols(Dl) ? max(abs(D - Dl)) : .)
        printf("{txt}  %-14s sigma %9.2e   Delta %9.2e   (%f x %f)\n", cs[c], e1, e2, rows(D), cols(D))
        if (missing(e1) | e1 > 1e-13) {
            printf("{err}m1: %s sigma differs from mg2_sigma by %g\n", cs[c], e1)
            nf++
        }
        if (missing(e2) | e2 > 1e-9) {
            printf("{err}m2: %s Delta differs from mg2_delta by %g\n", cs[c], e2)
            nf++
        }
    }
    return(nf)
}

// m3
real scalar tm_m3()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    real colvector th, sig, tr, tc
    real matrix K, Dl, D, X, O
    real scalar i, nf

    nf = 0
    tm_case("E_configural", par, th, K, sig, Dl)
    maps = tm_maps(K)
    D = dwls_delta(th, par, maps)
    tr = (1::maps[1].T) \ (maps[1].D :+ (1::maps[2].T))
    tc = J(0, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 4) tc = tc \ par.fidx[i]
    }
    X = D[tr, tc]
    if (rows(X) != cols(X) | X != I(rows(X))) {
        printf("{err}m3: Delta on the threshold rows and tau columns is not exactly the identity\n")
        nf++
    }
    O = D[tr, .]
    O[., tc] = J(rows(tr), rows(tc), 0)
    if (O != J(rows(O), cols(O), 0)) {
        printf("{err}m3: the threshold rows are not exactly 0 off the tau columns\n")
        nf++
    }
    return(nf)
}

// m4
real scalar tm_m4()
{
    struct dwlsmv_parmap scalar par, sp
    struct dwlsmv_statmap vector maps
    real colvector th, sig, idx, fr, s, s1
    real matrix K, Dl, D, D1
    real scalar g, i, n, off, nf

    nf = 0
    tm_case("E_configural", par, th, K, sig, Dl)
    maps = tm_maps(K)
    s = dwls_model(th, par, maps)
    D = dwls_delta(th, par, maps)
    off = 0
    for (g = 1; g <= 2; g++) {
        sp = dwlsmv_mg_sub(par, g, idx)
        // renumber the group's free rows 1..n in table order (no ties at configural)
        fr = J(0, 1, .)
        n = 0
        for (i = 1; i <= sp.nall; i++) {
            if (sp.free[i]) {
                n++
                fr = fr \ sp.fidx[i]
                sp.fidx[i] = n
            }
        }
        sp.npar = n
        s1 = dwls_model(th[fr], sp, maps[g])
        D1 = dwls_delta(th[fr], sp, maps[g])
        if (s[(off + 1)::(off + maps[g].D)] != s1 |
            D[(off + 1)::(off + maps[g].D), fr] != D1) {
            printf("{err}m4: group %f's block is not bitwise the single-group sigma/Delta\n", g)
            nf++
        }
        off = off + maps[g].D
    }
    return(nf)
}

// m5
real scalar tm_m5()
{
    struct dwlsmv_parmap scalar par, rv
    struct dwlsmv_statmap vector maps
    real colvector th, sig, o
    real matrix K, Dl

    tm_case("S", par, th, K, sig, Dl)
    maps = tm_maps(K)
    o = par.nall::1
    rv = par
    rv.name = par.name[o]
    rv.ptype = par.ptype[o]
    rv.row = par.row[o]
    rv.col = par.col[o]
    rv.free = par.free[o]
    rv.value = par.value[o]
    rv.fidx = par.fidx[o]
    rv.grp = par.grp[o]
    if (dwls_model(th, par, maps) != dwls_model(th, rv, maps) |
        dwls_delta(th, par, maps) != dwls_delta(th, rv, maps)) {
        printf("{err}m5: reversing the table rows changes sigma or Delta\n")
        return(1)
    }
    return(0)
}

// m6 cases: each is called under capture and must abort with rc 3498
real scalar tm_bad(real scalar which)
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps, m3, m0
    real colvector th, sig
    real matrix K, Dl

    tm_case("S", par, th, K, sig, Dl)
    maps = tm_maps(K)
    if (which == 1) {
        par.grp = par.grp[1::(par.nall - 1)]
        (void) dwls_model(th, par, maps)
    }
    else if (which == 2) {
        par.grp[par.nall] = 3
        (void) dwls_delta(th, par, maps)
    }
    else if (which == 3) {
        m3 = dwlsmv_statmap(3)
        m3[1] = maps[1]
        m3[2] = maps[2]
        m3[3] = maps[2]
        (void) dwls_model(th, par, m3)
    }
    else if (which == 4 | which == 5) {
        // group 1's rows alone, plus group 2's first nu row, as a single group
        tm_onegroup_nu(par, th, maps[1], which == 5)
    }
    else if (which == 6) {
        m0 = dwlsmv_statmap(0)
        (void) dwls_model(th, par, m0)
    }
    return(0)
}

void tm_onegroup_nu(struct dwlsmv_parmap scalar par, real colvector th,
                    struct dwlsmv_statmap scalar m1, real scalar delta)
{
    struct dwlsmv_parmap scalar sp
    real colvector idx, t1
    real scalar i, n

    sp = dwlsmv_mg_sub(par, 1, idx)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 5) break
    }
    sp.nall  = sp.nall + 1
    sp.name  = sp.name  \ par.name[i]
    sp.ptype = sp.ptype \ 5
    sp.row   = sp.row   \ par.row[i]
    sp.col   = sp.col   \ 1
    sp.free  = sp.free  \ 0
    sp.value = sp.value \ 0.1
    sp.fidx  = sp.fidx  \ 0
    sp.grp   = sp.grp   \ 1
    // renumber free rows 1..n; ties do not matter for a refusal
    n = 0
    for (i = 1; i <= sp.nall; i++) {
        if (sp.free[i]) {
            n++
            sp.fidx[i] = n
        }
    }
    sp.npar = n
    t1 = J(n, 1, 0.3)
    if (delta) (void) dwls_delta(t1, sp, m1)
    else       (void) dwls_model(t1, sp, m1)
}

end

local nfail 0
foreach t in m12 m3 m4 m5 {
    mata: st_local("k", strofreal(tm_`t'()))
    if `k' {
        display as error "test_mg_model: `t' `k' failure(s)"
        local nfail = `nfail' + `k'
    }
    else display as text "  `t' ok"
}

local l1 group key shorter than the table
local l2 a group key outside 1..G
local l3 three maps for a two-group table
local l4 a nu row at G = 1, dwls_model
local l5 a nu row at G = 1, dwls_delta
local l6 no map at all
forvalues i = 1/6 {
    capture mata: (void) tm_bad(`i')
    local rc = _rc
    display as text "  m6 `l`i'': rc `rc'"
    if `rc' != 3498 {
        display as error "test_mg_model: m6 `l`i'' gave rc `rc', wanted 3498"
        local nfail = `nfail' + 1
    }
}

if `nfail' {
    display as error _n "test_mg_model: `nfail' failure(s)"
    exit 9
}

display as result "test_mg_model: stacked sigma and Delta at G = 2 match the derivation on five tables, the configural threshold block is exactly the identity, each group's block is bitwise the single-group pair, row order does not matter, and malformed group keys are refused with rc 3498"
exit 0
