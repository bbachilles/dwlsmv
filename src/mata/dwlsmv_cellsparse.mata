*! dwlsmv -- classify a pair's contingency table for sparsity (Module 2,
*! spec Section 4.5.6)
*!
*! Executed under StataMP 18 on 2026-09-19. Driven by tests/unit/
*! test_diagnostics.do, whose pinned table is the one measured that day that
*! produced a clean, unflagged rho = 0.436 with an empty cell in it.
*!
*! Returns (code, nzero, nsparse, mincount) where
*!
*!     code = (nsparse > 0) + 2*(nzero > 0)
*!          = 0 none, 1 sparse only, 2 an empty cell, 3 both
*!
*! mirroring dwlsmv_catdegen()'s encoding so the two read the same way.
*!
*! WHY THIS EXISTS. Spec 4.5.6 lists three failure modes and Module 1 detects
*! only the first two, because those are properties of a single variable's
*! MARGIN. The third -- "category present but very sparse" -- is a property of
*! the PAIR's cell counts, and its symptom is the one the spec singles out as
*! dangerous: psi and hence Gamma are inflated and "no error is raised at all".
*!
*! That was not hypothetical. Measured 2026-09-19 on a 2000-observation 4x3
*! design whose top category held 11 observations, dwls_polychoric() returned
*! rho = 0.436 with converged = 1 -- a clean, plausible, entirely unflagged
*! answer -- from a table containing an EMPTY cell and a one-observation cell:
*!
*!         171    61    12
*!         383   405   102
*!         204   412   239
*!           0     1    10
*!
*! Module 1 passes this correctly: the category is sparse, not empty, so the
*! threshold is finite and there is nothing for dwlsmv_catdegen() to report.
*! The degeneracy is at the cell level, which nothing inspected until now.
*!
*! THIS CLASSIFIES AND DOES NOT ACT, exactly as dwlsmv_catdegen() does. Whether
*! a sparse table collapses, warns, aborts, or takes a continuity correction is
*! the open sparse-cell policy question (E-001), and it needs
*! Module 10's rank check before it can be answered. A routine that repaired
*! anything here would settle that question by accident.
*!
*! `minok` IS A REQUIRED ARGUMENT AND HAS NO DEFAULT. Where to put the line
*! between "sparse" and "fine" is precisely the open policy, so this routine
*! refuses to own it. A caller that has not decided should pass 1 and read
*! nzero, which needs no threshold to interpret: an empty cell is empty.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_cellsparse(real matrix n, real scalar minok)
{
    real scalar r, s, nzero, nsparse, mincount

    if (rows(n) < 2 | cols(n) < 2) {
        _error("dwlsmv_cellsparse: a contingency table needs at least 2x2 cells")
    }
    if (minok < 0 | minok >= .) {
        _error("dwlsmv_cellsparse: minok must be a non-negative number")
    }
    if (min(n) < 0) {
        _error("dwlsmv_cellsparse: counts cannot be negative")
    }

    nzero    = 0
    nsparse  = 0
    mincount = n[1, 1]

    for (r = 1; r <= rows(n); r++) {
        for (s = 1; s <= cols(n); s++) {
            if (n[r, s] < mincount) {
                mincount = n[r, s]
            }
            if (n[r, s] == 0) {
                nzero = nzero + 1
            }
            else if (n[r, s] < minok) {
                nsparse = nsparse + 1
            }
        }
    }

    return(((nsparse > 0) + 2 * (nzero > 0), nzero, nsparse, mincount))
}

end
