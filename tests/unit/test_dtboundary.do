*! Module 12, DIFFTEST validation check 11: dwlsmv_dtboundary() refuses H0
*! restrictions whose FIXED values force the boundary of the parameter space.
*!
*! Executed under StataMP 18 on 2026-10-08 and it passes: 20 cases, every
*! precondition holding, each refusal from the intended rule ((a) for r1-r5
*! and r11, (b) for r6-r10, read from the messages). Drafted by a read-only
*! agent; reviewed, run and landed by the main session. run_all.do discovers
*! it by name.
*!
*! MUTANTS, 2026-10-08, each turning this test red (killed by):
*!   Psi read as Phi, A = I (r8); absent lambda free (r6-r10); free beta,
*!   psi or lambda ignored (g5, g6, g7); free beta read at its start (g5,
*!   g9); (a) on 2x2 blocks only (r5); det == 0 on maximal blocks (r4, r11);
*!   a pivot or a Theta cutoff of 1e-2 (g3); th < 0 for th <= 0 (r6); (a)
*!   disabled (r1-r5, r11).
*!
*! Rule: spec Module 12, Validation contract, check 11 (decided 2026-10-08):
*!   (a) a principal submatrix of Psi whose entries are all fixed is singular
*!       or not positive semidefinite;
*!   (b) a Theta_jj = 1 - (Lambda Phi Lambda')_jj that is a function of fixed
*!       values only is <= 0.
*! Exact, no ratio cutoff. Absent rows count as fixed at 0.
*!
*! Every case is a hand-written canonical specification (spec 4.8), so this
*! file waits on no surface grammar. Each case runs in two steps:
*!   PRECONDITION  dwls_identify() alone must accept the specification (rc 0).
*!                 Otherwise a refusal could come from Module 3 (the rank
*!                 check, say) and an error case would pass for the wrong
*!                 reason. A precondition failure is a test-design failure.
*!   CHECK         dwls_identify() then dwlsmv_dtboundary(): rc 0 for an
*!                 accepted case, rc 3498 exactly for a refused one.
*!
*! Accepted (positive controls; each also names the mutant it kills):
*!   g1  Case E, marker metric, Psi free: nothing determined.
*!   g2  Case E, unit variance, everything free; psi21 absent = fixed 0, so
*!       {F1,F2} is an all-fixed block equal to I.
*!   g3  NEAR-BOUNDARY BUT LEGAL: psi21 fixed 0.999 and a loading fixed 0.999
*!       under unit variance (pivot 0.001999, Theta 0.001999). Kills a check
*!       with a cutoff.
*!   g4  Endogenous, legal: beta fixed 0.3, lambda[y4,F2] fixed 0.9;
*!       Theta = 1 - 0.81*1.09 = 0.1171.
*!   g5  Free beta into Anc: beta free, psi21 fixed -0.5, lambda fixed 1.05.
*!       Theta is -0.1025 at beta = 0 but 0.173 at beta = 0.5: not forced.
*!       Kills a check that zeroes (ignores) a free beta.
*!   g6  Free psi: y7 loads 0.8/0.8 (fixed) on F1/F2, psi21 free. Theta is
*!       -1.048 at psi21 = 0.6 but 0.36 at -0.5. Kills ignoring a free psi.
*!   g7  Free loading in the row: lambda[y7,F1] fixed 1.05, lambda[y7,F2]
*!       free, psi21 fixed -0.5. Theta -0.1025 at 0, 0.1725 at 0.5. Kills
*!       ignoring a free lambda.
*!   g8  Three factors, all Psi fixed, positive definite (0.5, 0.5, 0.4):
*!       control for r5.
*!   g9  As g4 but beta FREE with start 2: kills a check that reads a free
*!       parameter's start value (it would get Theta = 1 - 0.81*5 < 0).
*!
*! Refused, rc 3498:
*!   r1  DISCRIMINANT VALIDITY: unit variance, psi21 fixed at 1. (a)
*!   r2  psi21 fixed at -1. (a)
*!   r3  A variance fixed at 0 (marker-metric F2). (a)
*!   r4  A variance fixed at 0 with a covariance fixed at 0.3: the block
*!       [[1,.3],[.3,0]] is indefinite but its 1x1 block is singular. Kills
*!       a "det of the maximal block == 0" check.
*!   r5  Three unit-variance factors, psi fixed 0.5, 0.5, -0.5: singular PSD
*!       (eigenvalues 1.5, 1.5, 0) with NO singular 2x2 block. Kills a check
*!       of pairs only. LDL' pivots are exact: 1, 0.75, 0.
*!   r6  Unit variance, a loading fixed at 1: Theta = 0. (b), the spec's
*!       example.
*!   r7  Unit variance, a loading fixed at 1.2. (b)
*!   r8  ENDOGENOUS (case E shape): beta fixed 0.5, lambda[y4,F2] fixed 0.9;
*!       Phi22 = 1.25, Theta = -0.0125. A check that reads psi22 as F2's
*!       total variance gets 0.19 and passes: this kills the Psi/Phi
*!       confusion. Also kills treating the ABSENT lambda[y4,F1] as free.
*!   r9  y7 loads 0.8/0.8 (fixed), psi21 fixed 0.6: Theta = -1.048. (b)
*!       through a fixed covariance.
*!   r10 Marker metric, psi[F1,F1] fixed 1.2: Theta(y1) = -0.2. (b) only;
*!       the 1x1 block is positive definite.
*!   r11 Unit variances, psi21 fixed 1.2: indefinite with no singular
*!       principal submatrix, refused by (a) (user, 2026-10-08).
*!
*! Messages (spec check 11, decided 2026-10-13; added the same day, exemption
*! gate step 2a; the 20 cases above unchanged, each call now passes a
*! consequence sentence):
*!   m0  an empty consequence is refused, rc 3498, on an accepted table (g1)
*!   m1  r6's message names item y1, says "zero or negative" and ends with
*!       the caller's sentence; r1's states rule (a)'s condition and the
*!       caller's sentence
*! Mutants, predicted | observed 2026-10-13 (scratch worktree, diff printed):
*!   D1 the routine drops the caller's sentence   m1 15, 10 | m1 15, 10
*!   D2 the empty-consequence guard removed       m0        | m0
*!   D3 rule (b) names the row index, not item    m1 15     | m1 15
*!
*! Also: the call does not modify the table it is given (Mata passes by
*! reference).
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
        display as error "test_dtboundary: `f' failed to load"
        exit 9
    }
}

mata:

string scalar tdb_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

void tdb_add(string colvector WHY, string colvector SPEC, real colvector NIND,
             real colvector EXPECT, string scalar why, string rowvector lines,
             real scalar nind, real scalar expect)
{
    WHY    = WHY \ why
    SPEC   = SPEC \ tdb_spec(lines)
    NIND   = NIND \ nind
    EXPECT = EXPECT \ expect
}

void tdb_cases(string colvector WHY, string colvector SPEC, real colvector NIND,
               real colvector EXPECT)
{
    // ---------------- accepted ------------------------------------------
    tdb_add(WHY, SPEC, NIND, EXPECT, "g1 case E, marker metric, Psi free", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric marker F1 y1", "metric marker F2 y4",
        "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y5 free", "lambda F2 y6 free",
        "beta F2 F1 free", "psi F1 F1 free", "psi F2 F2 free"), 6, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g2 case E, unit variance, all free", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "beta F2 F1 free"), 6, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g3 near boundary: psi21 = 0.999, loading 0.999", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 fixed 0.999", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "psi F2 F1 fixed 0.999"), 6, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g4 endogenous, beta 0.3, loading 0.9: Theta 0.1171", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 fixed 0.9", "lambda F2 y5 free", "lambda F2 y6 free",
        "beta F2 F1 fixed 0.3"), 6, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g5 free beta into the ancestors: not determined", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 fixed 1.05", "lambda F2 y5 free", "lambda F2 y6 free",
        "beta F2 F1 free", "psi F2 F1 fixed -0.5"), 6, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g6 free psi21 under fixed cross-loadings", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6 y7",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "lambda F1 y7 fixed 0.8", "lambda F2 y7 fixed 0.8",
        "psi F2 F1 free"), 7, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g7 a free loading in the row: not determined", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6 y7",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "lambda F1 y7 fixed 1.05", "lambda F2 y7 free",
        "psi F2 F1 fixed -0.5"), 7, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g8 three factors, Psi fixed and positive definite", (
        "factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
        "metric unitvariance F1", "metric unitvariance F2", "metric unitvariance F3",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
        "psi F2 F1 fixed 0.5", "psi F3 F1 fixed 0.5", "psi F3 F2 fixed 0.4"), 9, 0)

    tdb_add(WHY, SPEC, NIND, EXPECT, "g9 a free beta's start value is never read", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 fixed 0.9", "lambda F2 y5 free", "lambda F2 y6 free",
        "beta F2 F1 free 2"), 6, 0)

    // ---------------- refused -------------------------------------------
    tdb_add(WHY, SPEC, NIND, EXPECT, "r1 discriminant validity: correlation fixed at 1", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "psi F2 F1 fixed 1"), 6, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r2 correlation fixed at -1", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "psi F2 F1 fixed -1"), 6, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r3 a variance fixed at 0", (
        "factor F1 F2", "ind y1 y2 y3 y4",
        "metric unitvariance F1", "metric marker F2 y4",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "psi F2 F2 fixed 0"), 4, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r4 variance 0 with a covariance 0.3 (indefinite block)", (
        "factor F1 F2", "ind y1 y2 y3 y4",
        "metric unitvariance F1", "metric marker F2 y4",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "psi F2 F2 fixed 0", "psi F2 F1 fixed 0.3"), 4, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r5 singular 3x3 Psi with no singular 2x2 block", (
        "factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
        "metric unitvariance F1", "metric unitvariance F2", "metric unitvariance F3",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
        "psi F2 F1 fixed 0.5", "psi F3 F1 fixed 0.5", "psi F3 F2 fixed -0.5"), 9, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r6 unit variance, loading fixed at 1: Theta = 0", (
        "factor F1", "ind y1 y2 y3",
        "metric unitvariance F1",
        "lambda F1 y1 fixed 1", "lambda F1 y2 free", "lambda F1 y3 free"), 3, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r7 unit variance, loading fixed at 1.2", (
        "factor F1", "ind y1 y2 y3",
        "metric unitvariance F1",
        "lambda F1 y1 fixed 1.2", "lambda F1 y2 free", "lambda F1 y3 free"), 3, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r8 endogenous: beta 0.5, loading 0.9, Theta -0.0125", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 fixed 0.9", "lambda F2 y5 free", "lambda F2 y6 free",
        "beta F2 F1 fixed 0.5"), 6, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r9 fixed cross-loadings 0.8/0.8, psi21 fixed 0.6", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6 y7",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "lambda F1 y7 fixed 0.8", "lambda F2 y7 fixed 0.8",
        "psi F2 F1 fixed 0.6"), 7, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r10 marker metric, factor variance fixed at 1.2", (
        "factor F1", "ind y1 y2 y3",
        "metric marker F1 y1",
        "lambda F1 y2 free", "lambda F1 y3 free",
        "psi F1 F1 fixed 1.2"), 3, 3498)

    tdb_add(WHY, SPEC, NIND, EXPECT, "r11 covariance 1.2 between unit variances (indefinite)", (
        "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
        "metric unitvariance F1", "metric unitvariance F2",
        "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
        "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
        "psi F2 F1 fixed 1.2"), 6, 3498)
}

// Structs stay inside functions (a project coding rule): the map is rebuilt per call.
void tdb_ident(string scalar spec, real scalar nind)
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar map

    map = dwlsmv_statmap_build(J(nind, 1, 3))
    par = dwls_identify(spec, map)
}

void tdb_run(string scalar spec, real scalar nind)
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar map

    map = dwlsmv_statmap_build(J(nind, 1, 3))
    par = dwls_identify(spec, map)
    dwlsmv_dtboundary(par, "H0", "TDB caller consequence.")
}

// As tdb_run with the caller's consequence sentence given (m0, m1).
void tdb_runc(string scalar spec, real scalar nind, string scalar cons)
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar map

    map = dwlsmv_statmap_build(J(nind, 1, 3))
    par = dwls_identify(spec, map)
    dwlsmv_dtboundary(par, "H0", cons)
}

// 1 if the table is unchanged by the call, on an accepted case (g4).
real scalar tdb_unchanged(string scalar spec, real scalar nind)
{
    struct dwlsmv_parmap  scalar par
    struct dwlsmv_statmap scalar map
    real colvector v0, f0

    map = dwlsmv_statmap_build(J(nind, 1, 3))
    par = dwls_identify(spec, map)
    v0  = par.value
    f0  = par.free
    dwlsmv_dtboundary(par, "H0", "TDB caller consequence.")
    return(v0 == par.value & f0 == par.free)
}

end

mata: WHY = J(0, 1, "")
mata: SPEC = J(0, 1, "")
mata: NIND = J(0, 1, .)
mata: EXPECT = J(0, 1, .)
mata: tdb_cases(WHY, SPEC, NIND, EXPECT)
mata: st_local("ncase", strofreal(rows(SPEC)))

local nfail 0
local npre 0

* Guard against an empty or truncated case list passing vacuously.
if `ncase' != 20 {
    display as error "test_dtboundary: expected 20 cases, built `ncase'"
    exit 9
}

forvalues i = 1/`ncase' {
    mata: st_local("why", WHY[`i'])
    mata: st_local("want", strofreal(EXPECT[`i']))

    * PRECONDITION: Module 3 alone accepts the specification.
    capture noisily mata: tdb_ident(SPEC[`i'], NIND[`i'])
    local rc = _rc
    if `rc' {
        display as error "    case `i' PRECONDITION FAILED (dwls_identify rc `rc'): `why'"
        local npre = `npre' + 1
        continue
    }

    capture noisily mata: tdb_run(SPEC[`i'], NIND[`i'])
    local rc = _rc
    if `rc' == `want' {
        display as text "    case `i' (rc `rc', as required): `why'"
    }
    else {
        display as error "    case `i' FAILED with rc `rc', wanted `want': `why'"
        local nfail = `nfail' + 1
    }
}

* The table is not modified through the reference (case g4 is the 4th).
mata: st_local("same", strofreal(tdb_unchanged(SPEC[4], NIND[4])))
if `same' != 1 {
    display as error "    dwlsmv_dtboundary() modified the parameter table it was given"
    local nfail = `nfail' + 1
}

* m0: an empty consequence is refused, even on an accepted table (g1, case 1).
capture noisily mata: tdb_runc(SPEC[1], NIND[1], " ")
if _rc != 3498 {
    display as error "    m0 FAILED: empty consequence gave rc `=_rc', wanted 3498"
    local nfail = `nfail' + 1
}
else {
    display as text "    m0 ok: an empty consequence is refused"
}

* m1: rule (b)'s message on r6 (case 15) names the item, says "zero or
* negative", and ends with the caller's sentence; rule (a)'s on r1 (case 10)
* states the condition and the caller's sentence (spec check 11, 2026-10-13).
mata: st_local("w15", WHY[15])
mata: st_local("w10", WHY[10])
if substr("`w15'", 1, 3) != "r6 " | substr("`w10'", 1, 3) != "r1 " {
    display as error "    m1 FAILED: case order changed (15: `w15'; 10: `w10')"
    local nfail = `nfail' + 1
}
tempfile lg
foreach c in 15 10 {
    quietly log using "`lg'", text replace name(tdbm)
    capture noisily mata: tdb_runc(SPEC[`c'], NIND[`c'], "TDB consequence `c' here.")
    local rc = _rc
    quietly log close tdbm
    mata: TDBLOG = invtokens(cat(st_local("lg"))')
    if `c' == 15 {
        local need `""residual variance of y1" "which is zero or negative" "TDB consequence 15 here.""'
    }
    else {
        local need `""the fixed values in Psi force a boundary" "two factors made one" "TDB consequence 10 here.""'
    }
    local bad 0
    foreach s of local need {
        mata: st_local("h", strofreal(strpos(TDBLOG, st_local("s")) > 0))
        if !`h' {
            display as error "    m1 case `c' FAILED: message lacks: `s'"
            local bad 1
        }
    }
    if `rc' != 3498 | `bad' {
        local nfail = `nfail' + 1
    }
    else {
        display as text "    m1 case `c' ok: condition and the caller's consequence"
    }
}

local ntot = `nfail' + `npre'
if `ntot' {
    display as error _n "test_dtboundary: `nfail' check failure(s), `npre' precondition failure(s)"
    exit 9
}

display as result "test_dtboundary: check 11 refuses every all-fixed singular Psi block and every fixed-determined Theta_jj <= 0 with rc 3498, and accepts near-boundary, free-dependent and endogenous-legal tables"
exit 0
