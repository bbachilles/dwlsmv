*! dwlsmv test harness -- Mata comparison primitives
*!
*! Run tests/unit/test_harness.do first: it self-tests the comparator below
*! against cases with known answers, so a syntax or logic error here surfaces
*! immediately rather than as a mysterious downstream failure. Do not trust any
*! other test in this suite until that one passes.
*!
*! Verified under StataMP 18 on Windows, 2026-09-17: every function here is
*! exercised by that self-test and it passes.
*!
*! Mata note: statements are terminated by a newline OR a semicolon, so the
*! last statement before a closing brace on the SAME line needs its own
*! semicolon. Every block here is written multi-line to sidestep that -- the
*! one-line form is what failed on the first run.
*!
*! This file does NOT issue `mata clear`, and must not. A definition file that
*! clears can only ever be loaded first: whatever was loaded before it is
*! silently wiped. Clearing belongs to the caller, and every test here opens
*! with `clear all`, which clears Mata too -- verified, not assumed. See
*! E-072, 2026-09-18.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Spec Section 7.2's acceptance rule, as a single number:
//      |a - b| <= eps_abs + eps_rel * |b|
// Returned as the ratio of gap to tolerance, so a value <= 1 is a pass
// regardless of which tier's tolerances were supplied. Reporting the ratio
// rather than a boolean means a near-miss is visible as 1.03 instead of
// "failed", which matters when calibrating a tolerance.
real scalar dwlsmv_gapratio(real matrix a, real matrix b,
                            real scalar eps_abs, real scalar eps_rel)
{
    real matrix tol

    if (rows(a) != rows(b) | cols(a) != cols(b)) return(.)
    if (eps_abs <= 0 & eps_rel <= 0)             return(.)

    tol = eps_abs :+ eps_rel :* abs(b)
    return(max(abs(a :- b) :/ tol))
}

// Location of the worst element, for diagnosis. A scalar max tells you a
// comparison failed; this tells you where, which is the difference between a
// five-minute fix and an afternoon.
real rowvector dwlsmv_worstcell(real matrix a, real matrix b,
                                real scalar eps_abs, real scalar eps_rel)
{
    real matrix tol, ratio
    real scalar i, j, best, bi, bj

    tol   = eps_abs :+ eps_rel :* abs(b)
    ratio = abs(a :- b) :/ tol

    best = -1
    bi   = 1
    bj   = 1
    for (i = 1; i <= rows(ratio); i++) {
        for (j = 1; j <= cols(ratio); j++) {
            if (ratio[i, j] > best) {
                best = ratio[i, j]
                bi   = i
                bj   = j
            }
        }
    }
    return((bi, bj, best, a[bi, bj], b[bi, bj]))
}

// Structural checks that need no reference values at all. These catch a large
// share of Gamma bugs before any comparison is made, and cost nothing.
real scalar dwlsmv_symmetry(real matrix g)
{
    return(max(abs(g :- g')))
}

real scalar dwlsmv_mineig(real matrix g)
{
    real matrix sym
    real colvector lambda

    sym    = (g :+ g') :/ 2
    lambda = symeigenvalues(sym)'
    return(min(lambda))
}

end
