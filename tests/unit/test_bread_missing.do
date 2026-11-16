*! Executed by the main session 2026-11-15 (StataMP 18): passes on its first run; mutants M1-M6 killed (E-097).
*! The caller-level "bread not evaluated" rule: dwlsmv_m11_bread(Delta, Gamma) (src/mata/dwlsmv_m11_bread.mata), used by
*! dwlsmv_m11_estimate() at both call sites. dwlsmv_breaddiag() keeps its strict contract (rc 3498 on a missing Delta).
*! Mata level, hand-built inputs, no fit: a real missing-Delta fit is not known to be reachable (see the probe draft), so
*! there is NO end-to-end check here that e(bread_ratio) is missing through the command; checks 3 and 5 stand in for the two
*! consumers (the D11 trigger; dwlsmv_m11_fit's bratio expression and dwlsmv_diagadvice), and check 6 guards the call sites.
*!
*! Checks
*!   1  A Delta with one missing value, ordinary-shaped (D = 6, npar = 3) and conditional-shaped (D = 11, npar = 5): the helper
*!      returns exactly (J(4, 2, .), J(4, 2, 0)) (4 x 4, == elementwise: every flag is ., every i, j is 0), no abort. The
*!      missing sits in a first, a middle and a last cell in turn, and once as .a.
*!   2  A finite Delta: the helper equals dwlsmv_breaddiag() with == and mreldif 0, for both shapes, and is NOT the
*!      not-evaluated shape, so a helper that always returns "not evaluated" fails.
*!   3  dwlsmv_m11_rectrigger() with a converged fit, an admissible A (J(6, 4, 0), as check 8 of test_condx_diag.do) and the
*!      not-evaluated Bd returns ""; with converged = 0 it returns only "did not converge", so Bd adds no cause.
*!   4  dwlsmv_breaddiag() itself still aborts with _rc == 3498 on the missing Delta.
*!   5  Bd[4, 1] is missing, so dwlsmv_m11_fit's expression (Bd[4, 1] == . ? . : Bd[4, 2]) is missing; and
*!      dwlsmv_diagadvice() on it prints no bread note.
*!   6  Source guard: dwlsmv_m11_estimate.mata has no dwlsmv_breaddiag( call outside comments and exactly two
*!      dwlsmv_m11_bread( calls. Text-level, since no fit reaches a missing Delta here.
*!
*! MUTANTS (draft's mutants.md)
*!  M1  helper returns J(0, 0, .) for a missing Delta.          Killed by 1 (shape) and 3/5 setup.
*!  M2  helper drops the missing(Delta) test.                   Killed by 1 (the call aborts, rc 3498).
*!  M3  helper's not-evaluated Bd has row 3 flag 1.             Killed by 1 and 3.
*!  M4  helper always returns the not-evaluated Bd.             Killed by 2.
*!  M5  estimate.mata: one call site back to breaddiag.         Killed by 6 only.
*!
*! Return codes: 0 pass, 9 fail.
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

global TB_NFAIL 0

capture program drop tb_cond
program tb_cond
    args ok tag
    if `ok' != 1 {
        display as error "test_bread_missing: FAIL `tag'"
        global TB_NFAIL = ${TB_NFAIL} + 1
    }
end

mata:

// a deterministic Delta (D x np, full column rank for these sizes) and a strictly positive diagonal Gamma
real matrix tb_delta(real scalar D, real scalar np)
{
    real matrix M
    real scalar i, j

    M = J(D, np, 0)
    for (i = 1; i <= D; i++) {
        for (j = 1; j <= np; j++) {
            M[i, j] = sin(0.7 * i * j + 0.3 * j) + (i == j)
        }
    }
    return(M)
}

real matrix tb_gamma(real scalar D)
{
    return(diag(1 :+ 0.1 * (1::D)))
}

real scalar tb_notev(real matrix B)
{
    if (rows(B) != 4 | cols(B) != 4) {
        return(0)
    }
    return(all(B :== (J(4, 2, .), J(4, 2, 0))))
}

// 1, 2: the helper on both shapes
real scalar tb_helper()
{
    real matrix Dl, Dm, G, B, B0
    real scalar nf, s, D, np, k, r, c

    nf = 0
    for (s = 1; s <= 2; s++) {
        if (s == 1) {
            D  = 6
            np = 3
        }
        else {
            D  = 11
            np = 5
        }
        G  = tb_gamma(D)
        Dl = tb_delta(D, np)

        // finite: identical to dwlsmv_breaddiag()
        B  = dwlsmv_m11_bread(Dl, G)
        B0 = dwlsmv_breaddiag(Dl, G)
        if (rows(B) != 4 | cols(B) != 4) {
            nf++
            printf("{err}  2: finite Delta (shape %f): the helper is %f x %f\n", s, rows(B), cols(B))
        }
        else if (!all(B :== B0) | mreldif(B, B0) != 0) {
            nf++
            printf("{err}  2: finite Delta (shape %f): the helper is not dwlsmv_breaddiag()\n", s)
        }
        if (tb_notev(B)) {
            nf++
            printf("{err}  2: finite Delta (shape %f) returned the not-evaluated shape\n", s)
        }

        // one missing value, in three places
        for (k = 1; k <= 3; k++) {
            if (k == 1) {
                r = 1
                c = 1
            }
            else if (k == 2) {
                r = trunc(D / 2)
                c = 2
            }
            else {
                r = D
                c = np
            }
            Dm = Dl
            Dm[r, c] = .
            if (!tb_notev(dwlsmv_m11_bread(Dm, G))) {
                nf++
                printf("{err}  1: shape %f, missing at (%f, %f): not exactly (J(4, 2, .), J(4, 2, 0))\n", s, r, c)
            }
        }
        // an extended missing value counts as missing too
        Dm = Dl
        Dm[1, 1] = .a
        if (!tb_notev(dwlsmv_m11_bread(Dm, G))) {
            nf++
            printf("{err}  1: shape %f, .a is not treated as missing\n", s)
        }
    }
    return(nf)
}

// 3, 5: the consumers of the not-evaluated Bd
real scalar tb_consumers()
{
    struct dwlsmv_fit scalar fit
    real matrix A, Bd, Dl, G, C
    string scalar why, msg
    string colvector ind, pname, slab
    real scalar nf, D

    nf = 0
    D  = 6
    Dl = tb_delta(D, 3)
    Dl[2, 2] = .
    Bd = dwlsmv_m11_bread(Dl, tb_gamma(D))
    if (!tb_notev(Bd)) {
        nf++
        printf("{err}  3: setup: Bd is not the not-evaluated shape\n")
        return(nf)
    }

    // 3: the trigger, built as check 8 of test_condx_diag.do builds it
    fit.converged = 1
    A   = J(6, 4, 0)
    why = dwlsmv_m11_rectrigger(fit, A, Bd)
    if (why != "") {
        nf++
        printf("{err}  3: not-evaluated Bd fired the trigger: '%s'\n", why)
    }
    fit.converged = 0
    why = dwlsmv_m11_rectrigger(fit, A, Bd)
    if (why != "did not converge") {
        nf++
        printf("{err}  3: converged = 0 with not-evaluated Bd: wanted only 'did not converge', got '%s'\n", why)
    }

    // 5: dwlsmv_m11_fit posts bratio from Bd[4, 1] == . ? . : Bd[4, 2]
    if (Bd[4, 1] != .) {
        nf++
        printf("{err}  5: Bd[4, 1] is not missing\n")
    }
    if ((Bd[4, 1] == . ? . : Bd[4, 2]) != .) {
        nf++
        printf("{err}  5: the bratio expression is not missing on the not-evaluated Bd\n")
    }
    // 5: dwlsmv_diagadvice() reads flags only for its notes; all-zero G flags and evaluated, unflagged counts
    // give no note at all, so any bread note would come from Bd.
    G     = J(7, 4, 0)
    C     = (0, 20, 1, 1 \ 0, 20, 1, 2)
    ind   = ("y1" \ "y2" \ "y3")
    slab  = ("a" \ "b" \ "c")
    pname = ("p1" \ "p2" \ "p3")
    msg = dwlsmv_diagadvice(G, C, Bd, ind, slab, pname)
    if (strpos(msg, "Identification of the parameters at the estimate") > 0) {
        nf++
        printf("{err}  5: dwlsmv_diagadvice() printed a bread note for the not-evaluated Bd\n")
    }
    return(nf)
}

// 4: the strict contract, so that the abort is a return code in Stata
void tb_strict()
{
    real matrix Dl

    Dl = tb_delta(6, 3)
    Dl[3, 1] = .
    (void) dwlsmv_breaddiag(Dl, tb_gamma(6))
}

// 6: the source guard
real scalar tb_source(string scalar fn)
{
    string colvector L
    string scalar t
    real scalar i, nb, nm, nf

    L  = cat(fn)
    nb = 0
    nm = 0
    nf = 0
    for (i = 1; i <= rows(L); i++) {
        t = strtrim(L[i])
        if (substr(t, 1, 2) == "//" | substr(t, 1, 2) == "*!") {
            continue
        }
        if (strpos(t, "dwlsmv_breaddiag(") > 0) {
            nb++
        }
        if (strpos(t, "dwlsmv_m11_bread(") > 0) {
            nm++
        }
    }
    if (nb != 0 | nm != 2) {
        nf++
        printf("{err}  6: dwlsmv_m11_estimate.mata has %f breaddiag call(s) and %f dwlsmv_m11_bread call(s); wanted 0 and 2\n", nb, nm)
    }
    return(nf)
}

end

mata: st_local("nh", strofreal(tb_helper()))
tb_cond `=(`nh' == 0)' "1-2 helper (`nh' failures)"

mata: st_local("nc", strofreal(tb_consumers()))
tb_cond `=(`nc' == 0)' "3, 5 consumers (`nc' failures)"

* 4: asserted == 3498, never != 0
capture noisily mata: tb_strict()
local rc4 = _rc
tb_cond `=(`rc4' == 3498)' "4 dwlsmv_breaddiag() on a missing Delta: wanted rc 3498, got `rc4'"

mata: st_local("ns", strofreal(tb_source("src/mata/dwlsmv_m11_estimate.mata")))
tb_cond `=(`ns' == 0)' "6 call sites (`ns' failures)"

if ${TB_NFAIL} {
    display as error "test_bread_missing: ${TB_NFAIL} check(s) failed"
    exit 9
}
display as text "test_bread_missing: all checks passed"
exit 0
