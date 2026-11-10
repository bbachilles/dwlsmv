*! Module 1-2: the diagnostics a user actually sees when data is degenerate.
*!
*! This test exists because of a specific question that turned out to have an
*! uncomfortable answer: when a pair has a sparse-cell problem, would the user
*! know? Measured before this file was written, the answer was no. A
*! 2000-observation 4x3 design whose top category held 11 observations
*! produced rho = 0.436 with converged = 1 -- a clean, plausible, entirely
*! unflagged result -- from a table containing an empty cell and a
*! one-observation cell. Module 1 passed it correctly, because the category
*! was sparse rather than empty. Nothing inspected the cells.
*!
*! So the checks here are about DETECTION AND MESSAGES, not about numerical
*! agreement. Five things:
*!
*!   1. dwlsmv_cellsparse() classifies the three conditions -- clean, sparse,
*!      empty-cell -- and its codes mirror dwlsmv_catdegen()'s so the two read
*!      alike. The pinned table is the one measured above, so a future change
*!      that stops detecting it fails here.
*!
*!   2-3. REMOVED 2026-10-07 (the user's): dwlsmv_celladvice(), which these
*!      checked, was never called by the command, and its text ("dwlsmv adds
*!      nothing", stale since O6's repair) was deleted with it. The command's
*!      own notes are tested in test_diagadvice.do and test_dwlsmv*.do.
*!
*!   4. The divergence itself is pinned numerically, so it is documented as
*!      KNOWN rather than discovered later. Module 2's own suite does not
*!      cover it: case_A has no empty cell.
*!
*!   5. The error messages are actionable. Asserting on message text is
*!      usually brittle and usually not worth it; here it is, because the
*!      whole point of the change was that the messages tell a user what to
*!      do. A message that silently reverts to "a boundary category is empty"
*!      would pass every other test in this suite.
*!
*! What this test deliberately does NOT do is assert any policy. Sparse cells
*! do not abort, are not repaired, and no continuity correction is applied.
*! That decision is open (E-001) and needs Module 10's rank check.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in src/mata/dwlsmv_types.mata                ///
             src/mata/dwlsmv_catcounts.mata            ///
             src/mata/dwlsmv_catdegen.mata             ///
             src/mata/dwlsmv_cumprop.mata              ///
             src/mata/dwlsmv_thresholds.mata           ///
             src/mata/dwlsmv_cellprob.mata             ///
             src/mata/dwlsmv_dcellprob_rho.mata        ///
             src/mata/dwlsmv_wmat.mata                 ///
             src/mata/dwlsmv_crosstab.mata             ///
             src/mata/dwlsmv_cellsparse.mata           ///
             src/mata/dwls_polychoric.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_diagnostics: `f' failed to load"
        exit 9
    }
}

mata:

real scalar dwlsmv_t_diagnostics()
{
    string colvector FAIL
    real matrix      clean, sparse, empty2, big
    real rowvector   z
    real scalar      i

    FAIL = J(0, 1, "")

    // ---- 1. classification ------------------------------------------------
    clean = (50, 60 \ 70, 80)
    z = dwlsmv_cellsparse(clean, 5)
    if (z[1] != 0 | z[2] != 0 | z[3] != 0 | z[4] != 50) {
        FAIL = FAIL \ "cellsparse flagged a healthy table"
    }

    // The table measured on 2026-09-19 that produced rho = 0.436 with
    // converged = 1 and no flag of any kind. Pinned so a regression is loud.
    big = (171,  61,  12 \
           383, 405, 102 \
           204, 412, 239 \
             0,   1,  10)
    z = dwlsmv_cellsparse(big, 5)
    if (z[1] != 3) {
        FAIL = FAIL \ ("cellsparse gave code " + strofreal(z[1]) +
                       " for the measured sparse table, expected 3 (both modes)")
    }
    if (z[2] != 1) {
        FAIL = FAIL \ ("cellsparse counted " + strofreal(z[2]) + " empty cells, expected 1")
    }
    // At minok = 5 exactly one non-empty cell qualifies: the count of 1. The
    // 10 and the 12 do not. Stated explicitly because the first version of
    // this assertion said 2, having counted the 10 without checking it
    // against the threshold.
    if (z[3] != 1) {
        FAIL = FAIL \ ("cellsparse counted " + strofreal(z[3]) +
                       " sparse cells at minok = 5, expected 1 (the count of 1; the 10 and 12 are above it)")
    }
    if (z[4] != 0) {
        FAIL = FAIL \ "cellsparse did not report a minimum count of 0"
    }

    // minok must actually move the answer, or it is decoration. At 11 the
    // count of 10 joins the count of 1, and the 12 still does not.
    z = dwlsmv_cellsparse(big, 11)
    if (z[3] != 2) {
        FAIL = FAIL \ ("cellsparse counted " + strofreal(z[3]) +
                       " sparse cells at minok = 11, expected 2 -- minok is not being applied")
    }

    // Sparse but no empty cell: code 1, not 3.
    sparse = (50, 60 \ 70,  2)
    z = dwlsmv_cellsparse(sparse, 5)
    if (z[1] != 1 | z[2] != 0 | z[3] != 1) {
        FAIL = FAIL \ ("cellsparse misclassified a sparse-but-nonempty table, code " +
                       strofreal(z[1]))
    }

    // Empty but nothing else sparse: code 2.
    empty2 = (50, 60 \ 70,  0)
    z = dwlsmv_cellsparse(empty2, 5)
    if (z[1] != 2 | z[2] != 1 | z[3] != 0) {
        FAIL = FAIL \ ("cellsparse misclassified an empty-cell table, code " +
                       strofreal(z[1]))
    }

    printf("{txt}test_diagnostics: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_diagnostics()))

* ---------------------------------------------------------------------------
* 5. The messages are actionable.
*
* Asserting on message text is usually brittle and usually not worth doing.
* It is worth doing here: the entire point of the change was that the message
* tells the user what to do, and a silent revert to the old wording would pass
* every other test in this suite.
* ---------------------------------------------------------------------------
* EVERY PATH IS CHECKED FOR rc 3498 SPECIFICALLY, not merely for aborting.
* That is the structural point of this block, and it is not fussiness.
*
* Mata's _error() accepts a message of at most 100 characters. Past that it
* raises rc 3300, "argument out of range", and the message you wrote is
* discarded -- so the error a user sees is not merely less helpful than
* intended, it is worse than the one you replaced, and it points at _error()
* rather than at the function that failed. Measured 2026-09-19: 100 characters
* works, 101 does not.
*
* An assertion of "it aborted" passes either way. Asserting the code is what
* distinguishes "the user got the guidance" from "the guidance was silently
* thrown away", which is exactly the failure this whole test file exists to
* prevent from recurring in a different form.
local nerr 0

* Indexed locals rather than a single list: `local probes "a b" "c d"` strips
* the quotes at assignment, so the loop then splits these on their internal
* spaces and backslashes and tests nonsense. Found by the first run reporting
* `dwlsmv_thresholds((1 raised 3000`.
local p1 "dwlsmv_thresholds((1 \ 1 \ 1 \ 1 \ 2 \ 2), 3)"
local p2 "dwlsmv_thresholds((0 \ 0 \ 2 \ 2), 3)"
local p3 "dwlsmv_catcounts((1 \ 2 \ 2 \ 1), 2)"
local p4 "dwlsmv_catcounts((0 \ 1 \ -1), 2)"
local p5 "dwlsmv_catcounts((0 \ 1 \ 1.5), 2)"
local p6 "dwlsmv_catcounts((0 \ 1 \ .), 2)"
local p7 "dwlsmv_cellsparse((1, 2 \ 3, 4), -1)"
local p8 "dwlsmv_cellsparse((1, 2, 3), 5)"

forvalues i = 1/8 {
    local p : copy local p`i'
    capture mata: `p'
    if _rc == 0 {
        display as error "test_diagnostics: `p' did not abort"
        local nerr = `nerr' + 1
    }
    else if _rc == 3300 {
        display as error "test_diagnostics: `p' raised 3300 -- its _error() message"
        display as error "  exceeds 100 characters and was DISCARDED. Move the detail into"
        display as error "  errprintf() and leave a short summary in _error()."
        local nerr = `nerr' + 1
    }
    else if _rc != 3498 {
        display as error "test_diagnostics: `p' raised `=_rc', expected 3498"
        local nerr = `nerr' + 1
    }
}

if `nfail' | `nerr' {
    display as error _n "test_diagnostics: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_diagnostics: sparse cells are detected, named, and explained"
exit 0
