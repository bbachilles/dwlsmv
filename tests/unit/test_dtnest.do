*! Module 12, DIFFTEST validation checks 8 and 9: dwlsmv_dtnest() refuses a pair
*! of parameter tables in which H0 is not nested in H1, and dwlsmv_dtchain()
*! orders the models and tests every adjacent pair.
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9; r12 and N14 added by the
*! main session. Executed 2026-10-08: green; all 21 mutants below each turned it red.
*!
*! Rule: spec Module 12, Validation contract, checks 8 and 9 (decided
*! 2026-10-08). Rows are (ptype, row, col); a row absent from a table is fixed at
*! 0; every parameter free in H0 is free in H1; every parameter fixed in both is
*! fixed at the same value; d = q1 - q0 >= 1 and d = df0 - df1.
*!
*! Every case is a hand-written canonical specification, so this file waits on no
*! surface grammar. Each case runs in two steps (as test_dtboundary.do):
*!   PRECONDITION  dwls_identify() alone must accept BOTH specifications (rc 0).
*!                 A precondition failure is a test-design failure, not a
*!                 finding about dwlsmv_dtnest.
*!   CHECK         dwlsmv_dtnest(): rc 0 and the returned d for an accepted
*!                 pair, rc 3498 exactly for a refused one.
*!
*! Shorthand. U0 = two factors, unit variance, three free loadings each,
*! psi21 ABSENT (fixed 0). U1 = U0 + "psi F2 F1 free". M0 = two factors, marker
*! metric (y1, y4), the other loadings free, psi F1 F1 and psi F2 F2 free, psi21
*! absent. All indicators have 3 categories. d is the expected return.
*!
*! Accepted (positive controls; each also names the mutant it kills):
*!   a1 U0 in U1 (free a covariance that is absent in H0).               d 1
*!   a2 psi21 fixed 0.3 in H0, free in H1. Kills a check that compares the
*!      values of a fixed row with a FREE row's start value.             d 1
*!   a3 psi21 fixed 0 (explicit) in H0, free in H1. Kills ordering or d by
*!      nall (the row counts tie, the free counts do not).               d 1
*!   a4 U1 in U1 + lambda[y4,F1] free (an absent loading freed).         d 1
*!   a5 U1 in U1 + lambda[y4,F1] fixed 0 + lambda[y1,F2] free: an explicit
*!      fixed 0 in H1 equals the absent row in H0. Kills "absent is not 0". d 1
*!   a6 three factors, psi fixed 0.5, 0.5, 0.4; H1 frees psi21 only: the
*!      OTHER fixed entries must still be matched, not skipped.          d 1
*!   a7 H0 gives lambda[y2,F1] the start 0.9, H1 the default. Kills a check
*!      that compares a free row's start.                                d 1
*!   a8 M0 in M0 + psi21 free (marker metric).                           d 1
*!   a9 lambda[y1,F1] fixed 0.8 in both. Kills "fixed in both must be 0". d 1
*!
*! Refused, rc 3498:
*!   r1  Reversed a1: H0 = U1, H1 = U0.
*!   r2  Counts agree on d = 1 but the sets are not nested: H0 = U1, H1 = U0 +
*!       lambda[y4,F1] free + lambda[y1,F2] free (psi21 free in H0, absent in
*!       H1). Kills a check of the counts alone.
*!   r3  Fixed in both at different values: psi21 0.3 against 0.5 (H1 also
*!       frees lambda[y4,F1], so d = 1). Kills a check that ignores values.
*!   r4  psi21 ABSENT in H0 (0), fixed 0.5 in H1. The row exists only in H1:
*!       kills a loop over H0's rows alone, and "absent is not 0".
*!   r5  psi21 fixed 0.5 in H0, ABSENT in H1. The row exists only in H0: kills a
*!       loop over H1's rows alone.
*!   r6  REPARAMETERIZATION ONLY: H0 = M0, H1 = U1 (d = 1). The marker model
*!       is the unit-variance orthogonal model rescaled, hence "nested" in U1
*!       only through a reparameterization: psi F1 F1 is free in H0 and fixed
*!       at 1 in H1.
*!   r7  H0 = H1 = U1: d = 0.
*!   r8  H0 = U0, H1 = U1 with H1's df shifted by 2: the tables give d = 1 and
*!       the df difference disagrees. Kills dropping rule 4.
*!   r9  Index collision: H0 = U0 + lambda[y4,F1] free, H1 = U0 + lambda[y1,F2]
*!       free + lambda[y5,F1] free (d = 1). H0's extra row is absent from H1,
*!       but H1 has lambda[y4,F2] free. Kills a key that drops the column.
*!   r10 Different indicator sets (3 against 4 indicators, one factor): nesting
*!       looks satisfied row by row; d = 3 but the df difference is -2. Backstop
*!       only: checks 4-5 own this refusal.
*!   r11 DIFFERENT MARKER (the spec's own example): H0 = M0 (marker y1), H1 =
*!       M0 with marker y2 and psi21 free (d = 1): lambda[y2,F1] is free in H0,
*!       fixed at 1 in H1.
*!   r12 FACTOR ORDER (user, 2026-10-08): H0 = U0; H1 declares "factor F2 F1"
*!       with F2 on y1-y3, F1 on y4-y6, psi21 free. Index-identical to U1, so
*!       only rule 0 (matched rows carry the same name) refuses it.
*!
*! Check 9 (dwlsmv_dtchain), models given in the order listed; ord is the
*! returned permutation, "least restricted last":
*!   c1 (U1 + lambda[y4,F1] free), U0, U1               ord 2 3 1  (q+2, q, q+1)
*!   c2 U1, U0 (two models, reversed)                    ord 2 1
*!   c3 U1, (U0 + psi21 fixed 0)    ord 2 1  -- nall ties, npar does not
*!   c4 U1, (U0 + lambda[y4,F1] free): equal npar, not nested   -> 3498
*!   c5 U0, (U0 + lambda[y4,F1] free), (U0 + psi21 free + lambda[y1,F2] free):
*!      the middle model is nested in the third? NO (its extra loading is absent
*!      from the third), while the endpoints ARE nested (U0 in the third). An
*!      endpoint-only test accepts it                           -> 3498
*!   c6 one model; c7 df of the wrong length; c8 names of the wrong length
*!                                                              -> 3498
*!   (s1, the stub's abort, retired 2026-10-13: the exemption is built and
*!   tested in test_dtgenerated.do.)
*!
*! Also: neither call modifies the tables it is given (checked on a2).
*!
*! MUTANTS to apply to dwlsmv_dtnest / dwlsmv_dtchain, each of which must turn
*! this file red (killed by):
*!   N1  counts only (rules 1-2 removed)                       r2, r3, r4, r5, r9
*!   N2  fixed values never compared                           r3, r4, r5
*!   N3  absent rows skipped instead of read as 0              r4, r5
*!   N4  only H0's rows visited                                r4
*!   N5  only H1's rows visited                                r5, r1
*!   N6  d >= 0 instead of d >= 1                              r7
*!   N7  rule 4 (df) removed                                   r8, r10
*!   N8  a FREE row's start value compared in rule 2           a7
*!   N9  value equality demanded also when H1's row is free    a2
*!   N10 key (ptype, row) without col                          r9
*!   N11 d from par.nall instead of par.npar                   a3 (d would be 0)
*!   N12 free-in-H0 / fixed-in-H1 accepted                     r6, r11
*!   N13 the call assigns to an argument (table changed)       unchanged check
*!   N14 rule 0 removed (names of matched rows not compared)   r12
*!   C1  dtchain does not sort                                 c1, c2
*!   C2  dtchain sorts descending                              c1, c2
*!   C3  only the first and last models compared               c5
*!   C4  sorted by par.nall                                     c3
*!   C5  n < 2 not refused                                      c6
*!   C6  df or names length unchecked                           c7, c8
*!   (C7, the stub returning 0, retired with s1 on 2026-10-13.)
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
        display as error "test_dtnest: `f' failed to load"
        exit 9
    }
}

mata:

string scalar tdn_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// ---- building blocks, as canonical lines (spec 4.8) ------------------------
string rowvector tdn_U0()
{
    return(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "metric unitvariance F1", "metric unitvariance F2",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free"))
}

// a7's H0: lambda[y2,F1] carries a start value
string rowvector tdn_U0s()
{
    return(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "metric unitvariance F1", "metric unitvariance F2",
            "lambda F1 y1 free", "lambda F1 y2 free 0.9", "lambda F1 y3 free",
            "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free"))
}

// a9: lambda[y1,F1] fixed at 0.8
string rowvector tdn_U0fx()
{
    return(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "metric unitvariance F1", "metric unitvariance F2",
            "lambda F1 y1 fixed 0.8", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free"))
}

string rowvector tdn_M0()
{
    return(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "metric marker F1 y1", "metric marker F2 y4",
            "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y5 free", "lambda F2 y6 free",
            "psi F1 F1 free", "psi F2 F2 free"))
}

// r11's H1: marker y2 for F1, psi21 free
string rowvector tdn_M0b()
{
    return(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "metric marker F1 y2", "metric marker F2 y4",
            "lambda F1 y1 free", "lambda F1 y3 free",
            "lambda F2 y5 free", "lambda F2 y6 free",
            "psi F1 F1 free", "psi F2 F2 free", "psi F2 F1 free"))
}

// test_dtboundary.do's g8, and the same with psi21 freed
string rowvector tdn_G8()
{
    return(("factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
            "metric unitvariance F1", "metric unitvariance F2", "metric unitvariance F3",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
            "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
            "psi F2 F1 fixed 0.5", "psi F3 F1 fixed 0.5", "psi F3 F2 fixed 0.4"))
}

string rowvector tdn_G8b()
{
    return(("factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
            "metric unitvariance F1", "metric unitvariance F2", "metric unitvariance F3",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
            "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
            "psi F2 F1 free", "psi F3 F1 fixed 0.5", "psi F3 F2 fixed 0.4"))
}

// r12's H1: factors declared F2 F1, F2 on y1-y3, F1 on y4-y6, psi21 free
string rowvector tdn_U1r()
{
    return(("factor F2 F1", "ind y1 y2 y3 y4 y5 y6",
            "metric unitvariance F1", "metric unitvariance F2",
            "lambda F2 y1 free", "lambda F2 y2 free", "lambda F2 y3 free",
            "lambda F1 y4 free", "lambda F1 y5 free", "lambda F1 y6 free",
            "psi F1 F2 free"))
}

// r10: one factor, 3 and 4 indicators
string rowvector tdn_T3()
{
    return(("factor F1", "ind y1 y2 y3", "metric unitvariance F1",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free"))
}

string rowvector tdn_T4()
{
    return(("factor F1", "ind y1 y2 y3 y4", "metric unitvariance F1",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F1 y4 free"))
}

// NUM columns: nind0, nind1, df shift added to H1's df, expected rc, expected d
void tdn_add(string colvector WHY, string colvector SP0, string colvector SP1,
             real matrix NUM, string scalar why, string rowvector l0,
             string rowvector l1, real rowvector num)
{
    WHY = WHY \ why
    SP0 = SP0 \ tdn_spec(l0)
    SP1 = SP1 \ tdn_spec(l1)
    NUM = NUM \ num
}

void tdn_cases(string colvector WHY, string colvector SP0, string colvector SP1,
               real matrix NUM)
{
    // ---------------- accepted ------------------------------------------
    tdn_add(WHY, SP0, SP1, NUM, "a1 free a covariance absent in H0",
            tdn_U0(), (tdn_U0(), "psi F2 F1 free"), (6, 6, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a2 psi21 fixed 0.3 freed",
            (tdn_U0(), "psi F2 F1 fixed 0.3"), (tdn_U0(), "psi F2 F1 free"), (6, 6, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a3 psi21 explicit fixed 0 freed",
            (tdn_U0(), "psi F2 F1 fixed 0"), (tdn_U0(), "psi F2 F1 free"), (6, 6, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a4 an absent loading freed",
            (tdn_U0(), "psi F2 F1 free"),
            (tdn_U0(), "psi F2 F1 free", "lambda F1 y4 free"), (6, 6, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a5 explicit fixed 0 in H1 equals absent in H0",
            (tdn_U0(), "psi F2 F1 free"),
            (tdn_U0(), "psi F2 F1 free", "lambda F1 y4 fixed 0", "lambda F2 y1 free"),
            (6, 6, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a6 three factors, other fixed psi still matched",
            tdn_G8(), tdn_G8b(), (9, 9, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a7 a free row's start differs",
            tdn_U0s(), (tdn_U0(), "psi F2 F1 free"), (6, 6, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a8 marker metric, psi21 freed",
            tdn_M0(), (tdn_M0(), "psi F2 F1 free"), (6, 6, 0, 0, 1))
    tdn_add(WHY, SP0, SP1, NUM, "a9 lambda fixed 0.8 in both",
            tdn_U0fx(), (tdn_U0fx(), "psi F2 F1 free"), (6, 6, 0, 0, 1))

    // ---------------- refused -------------------------------------------
    tdn_add(WHY, SP0, SP1, NUM, "r1 reversed a1",
            (tdn_U0(), "psi F2 F1 free"), tdn_U0(), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r2 d = 1 but not nested (psi21 free only in H0)",
            (tdn_U0(), "psi F2 F1 free"),
            (tdn_U0(), "lambda F1 y4 free", "lambda F2 y1 free"), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r3 fixed in both at 0.3 against 0.5",
            (tdn_U0(), "psi F2 F1 fixed 0.3"),
            (tdn_U0(), "psi F2 F1 fixed 0.5", "lambda F1 y4 free"), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r4 absent in H0, fixed 0.5 in H1",
            tdn_U0(),
            (tdn_U0(), "psi F2 F1 fixed 0.5", "lambda F1 y4 free"), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r5 fixed 0.5 in H0, absent in H1",
            (tdn_U0(), "psi F2 F1 fixed 0.5"),
            (tdn_U0(), "lambda F1 y4 free"), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r6 reparameterization only: marker model in unit-variance U1",
            tdn_M0(), (tdn_U0(), "psi F2 F1 free"), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r7 identical models, d = 0",
            (tdn_U0(), "psi F2 F1 free"), (tdn_U0(), "psi F2 F1 free"), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r8 df difference disagrees with d",
            tdn_U0(), (tdn_U0(), "psi F2 F1 free"), (6, 6, 2, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r9 lambda[y4,F1] against lambda[y4,F2]: same row, other column",
            (tdn_U0(), "lambda F1 y4 free"),
            (tdn_U0(), "lambda F2 y1 free", "lambda F1 y5 free"), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r10 3 against 4 indicators (backstop)",
            tdn_T3(), tdn_T4(), (3, 4, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r12 factors declared in a different order",
            tdn_U0(), tdn_U1r(), (6, 6, 0, 3498, .))
    tdn_add(WHY, SP0, SP1, NUM, "r11 different marker (spec's example)",
            tdn_M0(), tdn_M0b(), (6, 6, 0, 3498, .))
}

// Structs stay inside functions (a project coding rule).
void tdn_ident(string scalar spec0, string scalar spec1, real scalar n0, real scalar n1)
{
    struct dwlsmv_parmap  scalar p0, p1
    struct dwlsmv_statmap scalar m0, m1

    m0 = dwlsmv_statmap_build(J(n0, 1, 3))
    m1 = dwlsmv_statmap_build(J(n1, 1, 3))
    p0 = dwls_identify(spec0, m0)
    p1 = dwls_identify(spec1, m1)
}

real scalar tdn_run(string scalar spec0, string scalar spec1, real scalar n0,
                    real scalar n1, real scalar dfshift)
{
    struct dwlsmv_parmap  scalar p0, p1
    struct dwlsmv_statmap scalar m0, m1

    m0 = dwlsmv_statmap_build(J(n0, 1, 3))
    m1 = dwlsmv_statmap_build(J(n1, 1, 3))
    p0 = dwls_identify(spec0, m0)
    p1 = dwls_identify(spec1, m1)
    return(dwlsmv_dtnest(p0, p1, m0.D - p0.npar, m1.D - p1.npar + dfshift, "H0", "H1"))
}

// 1 if neither table is changed by the call (a2).
real scalar tdn_unchanged(string scalar spec0, string scalar spec1, real scalar n)
{
    struct dwlsmv_parmap  scalar p0, p1
    struct dwlsmv_statmap scalar m
    real colvector v0, f0, v1, f1
    real scalar d

    m  = dwlsmv_statmap_build(J(n, 1, 3))
    p0 = dwls_identify(spec0, m)
    p1 = dwls_identify(spec1, m)
    v0 = p0.value
    f0 = p0.free
    v1 = p1.value
    f1 = p1.free
    d  = dwlsmv_dtnest(p0, p1, m.D - p0.npar, m.D - p1.npar, "H0", "H1")
    return(v0 == p0.value & f0 == p0.free & v1 == p1.value & f1 == p1.free)
}

// ---- check 9 ---------------------------------------------------------------
string colvector tdn_chain_specs(real scalar which)
{
    if (which == 1) {
        return(tdn_spec((tdn_U0(), "psi F2 F1 free", "lambda F1 y4 free")) \
               tdn_spec(tdn_U0()) \ tdn_spec((tdn_U0(), "psi F2 F1 free")))
    }
    if (which == 2) {
        return(tdn_spec((tdn_U0(), "psi F2 F1 free")) \ tdn_spec(tdn_U0()))
    }
    if (which == 3) {
        return(tdn_spec((tdn_U0(), "psi F2 F1 free")) \
               tdn_spec((tdn_U0(), "psi F2 F1 fixed 0")))
    }
    if (which == 4) {
        return(tdn_spec((tdn_U0(), "psi F2 F1 free")) \
               tdn_spec((tdn_U0(), "lambda F1 y4 free")))
    }
    if (which == 5) {
        return(tdn_spec(tdn_U0()) \ tdn_spec((tdn_U0(), "lambda F1 y4 free")) \
               tdn_spec((tdn_U0(), "psi F2 F1 free", "lambda F2 y1 free")))
    }
    return(tdn_spec(tdn_U0()))
}

// mode 0: as given; 1: df one element short; 2: names one element short
real colvector tdn_chain(string colvector SP, real scalar mode)
{
    struct dwlsmv_parmap  vector P
    struct dwlsmv_statmap scalar m
    real colvector df
    string colvector nm
    real scalar k, n

    n  = rows(SP)
    m  = dwlsmv_statmap_build(J(6, 1, 3))
    P  = dwlsmv_parmap(n, 1)
    df = J(n, 1, .)
    nm = J(n, 1, "")
    for (k = 1; k <= n; k++) {
        P[k]  = dwls_identify(SP[k], m)
        df[k] = m.D - P[k].npar
        nm[k] = "M" + strofreal(k)
    }
    if (mode == 1) {
        df = df[(1::(n - 1))]
    }
    if (mode == 2) {
        nm = nm[(1::(n - 1))]
    }
    return(dwlsmv_dtchain(P, df, nm))
}

void tdn_chain_ident(string colvector SP)
{
    struct dwlsmv_statmap scalar m
    struct dwlsmv_parmap  scalar p
    real scalar k

    m = dwlsmv_statmap_build(J(6, 1, 3))
    for (k = 1; k <= rows(SP); k++) {
        p = dwls_identify(SP[k], m)
    }
}

end

mata: WHY = J(0, 1, "")
mata: SP0 = J(0, 1, "")
mata: SP1 = J(0, 1, "")
mata: NUM = J(0, 5, .)
mata: tdn_cases(WHY, SP0, SP1, NUM)
mata: st_local("ncase", strofreal(rows(SP0)))

local nfail 0
local npre 0

* Guard against an empty or truncated case list passing vacuously.
if `ncase' != 21 {
    display as error "test_dtnest: expected 21 cases, built `ncase'"
    exit 9
}

forvalues i = 1/`ncase' {
    mata: st_local("why", WHY[`i'])
    mata: st_local("want", strofreal(NUM[`i', 4]))
    mata: st_local("wantd", strofreal(NUM[`i', 5]))

    * PRECONDITION: Module 3 alone accepts both specifications.
    capture noisily mata: tdn_ident(SP0[`i'], SP1[`i'], NUM[`i', 1], NUM[`i', 2])
    local rc = _rc
    if `rc' {
        display as error "    case `i' PRECONDITION FAILED (dwls_identify rc `rc'): `why'"
        local npre = `npre' + 1
        continue
    }

    local dgot .
    capture noisily mata: st_local("dgot", strofreal(tdn_run(SP0[`i'], SP1[`i'], NUM[`i', 1], NUM[`i', 2], NUM[`i', 3])))
    local rc = _rc
    if `rc' == `want' & (`want' != 0 | "`dgot'" == "`wantd'") {
        display as text "    case `i' (rc `rc', d `dgot', as required): `why'"
    }
    else {
        display as error "    case `i' FAILED with rc `rc' and d `dgot', wanted rc `want' and d `wantd': `why'"
        local nfail = `nfail' + 1
    }
}

* Neither table is modified through the reference (a2 is the 2nd case).
mata: st_local("same", strofreal(tdn_unchanged(SP0[2], SP1[2], NUM[2, 1])))
if `same' != 1 {
    display as error "    dwlsmv_dtnest() modified a parameter table it was given"
    local nfail = `nfail' + 1
}

* ---- check 9: dwlsmv_dtchain ------------------------------------------------
* accepted: the returned permutation, least restricted last
local want1 2 3 1
local want2 2 1
local want3 2 1
forvalues c = 1/3 {
    capture noisily mata: tdn_chain_ident(tdn_chain_specs(`c'))
    if _rc {
        display as error "    chain c`c' PRECONDITION FAILED (dwls_identify rc `=_rc')"
        local npre = `npre' + 1
        continue
    }
    local got ""
    capture noisily mata: st_local("got", invtokens(strofreal(tdn_chain(tdn_chain_specs(`c'), 0))'))
    local rc = _rc
    if `rc' == 0 & "`got'" == "`want`c''" {
        display as text "    chain c`c' (order `got', as required)"
    }
    else {
        display as error "    chain c`c' FAILED with rc `rc', order `got', wanted `want`c''"
        local nfail = `nfail' + 1
    }
}

* refused: c4 and c5 (nesting); c6-c8 (arguments); all rc 3498 exactly
foreach c in 4 5 {
    capture noisily mata: tdn_chain_ident(tdn_chain_specs(`c'))
    if _rc {
        display as error "    chain c`c' PRECONDITION FAILED (dwls_identify rc `=_rc')"
        local npre = `npre' + 1
        continue
    }
    capture mata: (void) tdn_chain(tdn_chain_specs(`c'), 0)
    local rc = _rc
    display as text "    chain c`c': rc `rc'"
    if `rc' != 3498 {
        capture noisily mata: (void) tdn_chain(tdn_chain_specs(`c'), 0)
        display as error "test_dtnest: chain c`c' gave rc `rc', not 3498"
        local nfail = `nfail' + 1
    }
}
capture mata: (void) tdn_chain(tdn_chain_specs(6), 0)
local rc = _rc
display as text "    chain c6 (one model): rc `rc'"
if `rc' != 3498 {
    display as error "test_dtnest: one model gave rc `rc', not 3498"
    local nfail = `nfail' + 1
}
capture mata: (void) tdn_chain(tdn_chain_specs(1), 1)
local rc = _rc
display as text "    chain c7 (df short): rc `rc'"
if `rc' != 3498 {
    display as error "test_dtnest: df of the wrong length gave rc `rc', not 3498"
    local nfail = `nfail' + 1
}
capture mata: (void) tdn_chain(tdn_chain_specs(1), 2)
local rc = _rc
display as text "    chain c8 (names short): rc `rc'"
if `rc' != 3498 {
    display as error "test_dtnest: names of the wrong length gave rc `rc', not 3498"
    local nfail = `nfail' + 1
}

local ntot = `nfail' + `npre'
if `ntot' {
    display as error _n "test_dtnest: `nfail' check failure(s), `npre' precondition failure(s)"
    exit 9
}

display as result "test_dtnest: check 8 accepts nested pairs and refuses every non-nested, reparameterized, d <= 0 or df-inconsistent pair with rc 3498; check 9 orders by npar and tests adjacent pairs"
exit 0
