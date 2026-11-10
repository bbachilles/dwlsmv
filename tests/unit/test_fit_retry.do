*! Module 7: optimization recovery in dwls_fit (spec Module 7, 2026-11-01;
*! E-013).
*!
*! Executed under StataMP 18 on 2026-11-01 and it passes. Mutants: E-074, 2026-11-01.
*!
*!   1  the case that motivated it (mg_B, F -> y1-y3, H -> y4 y5, thresholds):
*!      the command converges and says it retried; test_standardized_g.do check 2
*!      compares that fit with lavaan (the exclusion made for this cell is removed),
*!      and tools outside this test measured F = .000576643 at the same optimum;
*!   2  no trigger, G = 2 and G = 1: dwls_fit equals ULS then dwlsmv_fit_from,
*!      theta exactly, and no note;
*!   3  a failed retry, G = 2: two different ridge starts; the first run's result
*!      is returned, converged 0, with the failure note;
*!   4  G = 1, forced: a first start on a ridge and Module 3's second; the
*!      retried fit equals the direct fit from Module 3's starts exactly.
*!   5  a marker fit takes no retry: the replicate where a retry converged and
*!      erased D11's trigger (E-012, b2f n = 200, replicate
*!      16, its draws replayed here). The first run fails, the retry alone would
*!      converge, and the command reports e(recovery) recovered with no retry note.
*!   6  G = 1, unit variance, a natural failure a retry rescues, through the command.
*!   7  the same data under marker: D11's unit-variance run takes no retry, so
*!      e(recovery) is failed where check 6 converges. This pins the behavior as built.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
quietly do tests/dwlsmv_assert.mata
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

local nfail 0
scalar TFRF = 9
local m2 (F -> y1 y2 y3) (H -> y4 y5)
local ok  "it converged from"
local bad "or, on a retry"

capture program drop tfr_log
program tfr_log
    args cmd
    tempfile lg
    capture log close tfr
    quietly log using "`lg'", text replace name(tfr)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tfr
    mata: TFRLOG = invtokens(cat(st_local("lg"))')
    c_local tfr_rc `rc'
end

capture program drop tfr_has
program tfr_has
    args s
    mata: st_local("h", strofreal(strpos(TFRLOG, st_local("s")) > 0))
    c_local hit `h'
end

mata:
struct dwlsmv_stage1 vector tfr_sts(real matrix D, real scalar G)
{
    struct dwlsmv_stage1 vector sts
    real scalar g
    real matrix Y

    sts = dwlsmv_stage1(G, 1)
    for (g = 1; g <= G; g++) {
        Y = select(D[., 2..6], D[., 1] :== g)
        sts[g] = dwls_sample_stats(Y, J(5, 1, 3))
    }
    return(sts)
}

struct dwlsmv_parmap scalar tfr_par(struct dwlsmv_stage1 vector sts)
{
    struct dwlsmv_statmap vector maps
    real scalar g

    maps = dwlsmv_statmap(length(sts), 1)
    for (g = 1; g <= length(sts); g++) {
        maps[g] = sts[g].map
    }
    return(dwls_identify(subinstr(st_global("e(canon)"), ";", char(10)) + char(10), maps, sts))
}

real colvector tfr_start(struct dwlsmv_parmap scalar par)
{
    real colvector start
    real scalar i

    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            start[par.fidx[i]] = par.value[i]
        }
    }
    return(start)
}

real colvector tfr_uls(struct dwlsmv_stage1 vector sts, struct dwlsmv_parmap scalar par)
{
    struct dwlsmv_statmap vector maps
    real colvector s, w, thu
    real scalar N, fu, iu

    dwlsmv_fit_stack(sts, s, w, N, maps)
    (void) dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, tfr_start(par), "d2", "nr", 300, thu, fu, iu)
    return(thu)
}

void tfr_set(real scalar v)
{
    st_numscalar("TFRF", v)
}

// 5: replicate 16 of the d11r.do stream, left in memory as the data
void tfr_replay(real scalar nt)
{
    real matrix Y
    real colvector ix
    real scalar r

    rseed(20261101)
    for (r = 1; r <= 60; r++) {
        (void) runiform(300, 1)
    }
    Y = st_data(., "y1 y2 y3 y4 y5")
    for (r = 1; r <= 16; r++) {
        ix = ceil(runiform(200, 1) :* rows(Y))
    }
    Y = Y[ix, .]
    st_dropvar(.)
    st_addobs(200 - st_nobs())
    st_store(., st_addvar("byte", ("y1", "y2", "y3", "y4", "y5")), Y)
}

// 5: the precondition: on these data the first marker run fails, a retry converges
real scalar tfr_pre(string scalar canon)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar a, b
    real colvector thu, start

    sts = dwls_sample_stats(st_data(., "y1 y2 y3 y4 y5"), J(5, 1, 3))
    par = dwls_identify(subinstr(canon, ";", char(10)), sts[1].map)
    dwlsmv_fit_uls(sts, par, start, thu)
    a = dwlsmv_fit_from(sts, par, thu)
    b = dwlsmv_fit_from(sts, par, start)
    return(a.converged != 0 | b.converged != 1)
}

// 6, 7: a natural G = 1 unit-variance failure a retry rescues (replicate 1 of
// E-012, scan_g1, n = 150 from mg_B), left in memory as the data
void tfr_replay2()
{
    real matrix Y
    real colvector ix

    Y = st_data(., "y1 y2 y3 y4 y5")
    rseed(20261102)
    ix = ceil(runiform(150, 1) :* rows(Y))
    Y = Y[ix, .]
    st_dropvar(.)
    st_addobs(150)
    st_store(., st_addvar("byte", ("y1", "y2", "y3", "y4", "y5")), Y)
}

// 2: dwls_fit against the direct path
real scalar tfr_direct(real matrix D, real scalar G)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar a, b

    sts = tfr_sts(D, G)
    par = tfr_par(sts)
    a   = dwls_fit(sts, par)
    b   = dwlsmv_fit_from(sts, par, tfr_uls(sts, par))
    return(!a.converged | a.theta != b.theta | a.iterations != b.iterations)
}

// 3: both starts on the ridge
real scalar tfr_fail(real matrix D)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar a, b, c
    real colvector r, r2

    sts = tfr_sts(D, 2)
    par = tfr_par(sts)
    r   = tfr_uls(sts, par)
    r2  = r :* 1.001
    c   = dwlsmv_fit_from(sts, par, r2)
    a   = dwlsmv_fit_retry(sts, par, r, r2)
    b   = dwlsmv_fit_from(sts, par, r)
    if (c.converged | c.theta == b.theta) {
        printf("{err}  3: the second start does not separate the two runs\n")
        return(1)
    }
    return(a.converged != 0 | b.converged != 0 | a.theta != b.theta)
}

// 4: G = 1, first start on a ridge by construction, the loadings of H set
// as the G = 2 ridge left them
real scalar tfr_g1(real matrix D)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar first, a, b
    real colvector r
    real scalar i

    sts = tfr_sts(D, 1)
    par = tfr_par(sts)
    r   = tfr_start(par)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.ptype[i] == 1 & par.col[i] == 2) {
            r[par.fidx[i]] = (par.row[i] == 4 ? -1e-3 : 12)
        }
    }
    first = dwlsmv_fit_from(sts, par, r)
    printf("{txt}  4 first run from the ridge: converged %f, F %g\n", first.converged, first.fmin)
    if (first.converged) {
        printf("{err}  4 FAILED: the forced start did not fail, so nothing was tested\n")
        return(1)
    }
    a = dwlsmv_fit_retry(sts, par, r, tfr_start(par))
    b = dwlsmv_fit_from(sts, par, tfr_start(par))
    return(!a.converged | a.theta != b.theta)
}
end

import delimited using "fixtures/data/mg_B.csv", clear

* ---- 1 ---------------------------------------------------------------------
tfr_log "dwlsmv `m2', group(g) invariance(thresholds) metric(unitvariance)"
tfr_has "`ok'"
local h1 `hit'
if `tfr_rc' | !`h1' | e(converged) != 1 {
    display as error "  1 FAILED (rc `tfr_rc', note `h1', converged " e(converged) ")"
    local ++nfail
}
else display as text "  1 ok, " e(iterations) " iterations, T3 " e(T3)
mata: TFRD = st_data(., "g y1 y2 y3 y4 y5")

* ---- 3 (same e(canon)) -----------------------------------------------------
tfr_log "mata: tfr_set(tfr_fail(TFRD))"
tfr_has "`bad'"
if `tfr_rc' | TFRF != 0 | !`hit' {
    display as error "  3 FAILED (rc `tfr_rc', result " TFRF ", note `hit')"
    local ++nfail
}
else display as text "  3 ok"

* ---- 2, G = 2 --------------------------------------------------------------
quietly dwlsmv `m2', group(g) invariance(metric) metric(unitvariance)
tfr_log "mata: tfr_set(tfr_direct(TFRD, 2))"
tfr_has "ULS starting values"
if `tfr_rc' | TFRF != 0 | `hit' {
    display as error "  2 G = 2 FAILED (rc `tfr_rc', result " TFRF ", note `hit')"
    local ++nfail
}
else display as text "  2 G = 2 ok"

* ---- 2 and 4, G = 1 --------------------------------------------------------
keep if g == 1
mata: TFRD = st_data(., "g y1 y2 y3 y4 y5")
quietly dwlsmv `m2', metric(unitvariance)
tfr_log "mata: tfr_set(tfr_direct(TFRD, 1))"
tfr_has "ULS starting values"
if `tfr_rc' | TFRF != 0 | `hit' {
    display as error "  2 G = 1 FAILED (rc `tfr_rc', result " TFRF ", note `hit')"
    local ++nfail
}
else display as text "  2 G = 1 ok"
tfr_log "mata: tfr_set(tfr_g1(TFRD))"
tfr_has "`ok'"
if `tfr_rc' | TFRF != 0 | !`hit' {
    display as error "  4 FAILED (rc `tfr_rc', result " TFRF ", note `hit')"
    local ++nfail
}
else display as text "  4 ok"

* ---- 5, marker: no retry ---------------------------------------------------
import delimited using "fixtures/data/marker_trap.csv", clear
local nt = _N
import delimited using "fixtures/data/mg_B.csv", clear
mata: tfr_replay(`nt')
quietly _dwlsmv_translate , paths(`m2') metric(marker)
local canon5 `"`r(canon)'"'
mata: tfr_set(tfr_pre(st_local("canon5")))
tfr_log "dwlsmv `m2', metric(marker)"
tfr_has "ULS starting values"
local h5 `hit'
tfr_has "restarted from the unit-variance"
if `tfr_rc' | TFRF != 0 | `h5' | !`hit' | "`e(recovery)'" != "recovered" {
    display as error "  5 FAILED (rc `tfr_rc', precondition " TFRF ", retry note `h5', D11 note `hit', `e(recovery)')"
    local ++nfail
}
else display as text "  5 ok, e(recovery) `e(recovery)'"

* ---- 6, unit variance, G = 1, through the command: the retry applies -------
import delimited using "fixtures/data/mg_B.csv", clear
mata: tfr_replay2()
quietly _dwlsmv_translate , paths(`m2') metric(unitvariance)
local canon6 `"`r(canon)'"'
mata: tfr_set(tfr_pre(st_local("canon6")))
tfr_log "dwlsmv `m2', metric(unitvariance)"
tfr_has "`ok'"
if `tfr_rc' | TFRF != 0 | !`hit' | e(converged) != 1 {
    display as error "  6 FAILED (rc `tfr_rc', precondition " TFRF ", note `hit', converged " e(converged) ")"
    local ++nfail
}
else display as text "  6 ok"

* ---- 7, same data, marker: D11's unit-variance run takes no retry -----------
tfr_log "dwlsmv `m2', metric(marker)"
tfr_has "ULS starting values"
if `tfr_rc' | `hit' | "`e(recovery)'" != "failed" | e(converged) != 0 {
    display as error "  7 FAILED (rc `tfr_rc', retry note `hit', `e(recovery)', converged " e(converged) ")"
    local ++nfail
}
else display as text "  7 ok, e(recovery) `e(recovery)'"

adopath - "`here'/src/ado"
adopath - "`here'/build"
if `nfail' {
    display as error "test_fit_retry: `nfail' failure(s)"
    exit 9
}
display as text "test_fit_retry: all checks passed"
