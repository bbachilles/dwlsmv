*! Self-test of the comparison harness itself.
*!
*! Run this before any other test. It has no dependency on dwlsmv code, so it
*! is runnable the moment Stata is available -- before Module 0 exists.
*!
*! Verified passing under StataMP 18 on Windows, 2026-09-17.
*!
*! The point: a test suite whose comparator is untested reports whatever the
*! comparator happens to do. Each case below has an answer that can be
*! computed by hand, so a broken comparator cannot quietly pass.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

do "tests/dwlsmv_assert.mata"

local failures 0

mata:
    // ---- exact equality passes at any tolerance --------------------------
    a = (1, 2 \ 3, 4)
    r = dwlsmv_gapratio(a, a, 1e-12, 0)
    printf("{txt}identical matrices, ratio %g (expect 0)\n", r)
    ok1 = (r == 0)

    // ---- a gap exactly at the tolerance is a pass, just barely -----------
    // |1.0 - 1.0| = 0 elsewhere; the (1,1) cell is off by exactly eps_abs
    b = (1.001, 2 \ 3, 4)
    r = dwlsmv_gapratio(b, a, 0.001, 0)
    printf("{txt}gap exactly at tolerance, ratio %g (expect 1)\n", r)
    ok2 = (reldif(r, 1) < 1e-10)

    // ---- twice the tolerance is a clear fail ----------------------------
    r = dwlsmv_gapratio(b, a, 0.0005, 0)
    printf("{txt}gap at twice tolerance, ratio %g (expect 2)\n", r)
    ok3 = (reldif(r, 2) < 1e-10)

    // ---- the relative term must actually bite on large values -----------
    // Expected 1000, actual 1001, so the absolute gap is 1. The tolerance is
    // eps_rel * |expected| = eps_rel * 1000 -- note it keys off the EXPECTED
    // value, the second argument, not the actual.
    //   eps_rel = 2e-3 -> tolerance 2.0 -> ratio 0.5, comfortably a pass
    //   eps_rel = 1e-4 -> tolerance 0.1 -> ratio 10,  comfortably a fail
    // Both are deliberately off the boundary; ok2 above tests the boundary.
    c = (1000)
    d = (1001)
    r1 = dwlsmv_gapratio(d, c, 0, 2e-3)
    r2 = dwlsmv_gapratio(d, c, 0, 1e-4)
    printf("{txt}relative term, ratios %g and %g (expect 0.5 then 10)\n", r1, r2)
    ok4 = (reldif(r1, 0.5) < 1e-10 & reldif(r2, 10) < 1e-10)

    // ---- shape mismatch must not silently broadcast ---------------------
    // Mata would happily broadcast a conformable scalar; returning missing
    // instead is the whole reason for the dimension guard.
    r = dwlsmv_gapratio((1, 2), a, 1e-6, 1e-6)
    printf("{txt}shape mismatch, ratio %g (expect missing)\n", r)
    ok5 = (r == .)

    // ---- worst-cell locator points at the right element -----------------
    e = (1, 2 \ 3, 9)
    w = dwlsmv_worstcell(e, a, 1e-6, 0)
    printf("{txt}worst cell (%g,%g), ratio %g\n", w[1], w[2], w[3])
    ok6 = (w[1] == 2 & w[2] == 2)

    // ---- structural checks ----------------------------------------------
    sym = (2, 1 \ 1, 2)
    printf("{txt}symmetry of a symmetric matrix: %g (expect 0)\n", dwlsmv_symmetry(sym))
    printf("{txt}min eigenvalue of [[2,1],[1,2]]: %g (expect 1)\n", dwlsmv_mineig(sym))
    ok7 = (dwlsmv_symmetry(sym) == 0 & reldif(dwlsmv_mineig(sym), 1) < 1e-10)

    asym = (2, 1 \ 0, 2)
    printf("{txt}symmetry of an asymmetric matrix: %g (expect 1)\n", dwlsmv_symmetry(asym))
    ok8 = (reldif(dwlsmv_symmetry(asym), 1) < 1e-10)

    st_local("allok", strofreal(ok1 & ok2 & ok3 & ok4 & ok5 & ok6 & ok7 & ok8))
end

if `allok' {
    display as result "PASS: comparison harness self-test"
}
else {
    display as error "FAIL: comparison harness self-test -- do not trust other tests"
    exit 9
}
