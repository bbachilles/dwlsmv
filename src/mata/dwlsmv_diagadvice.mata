*! dwlsmv -- Module 10's rows as notes a user can read (spec Section 5)
*!
*! Executed under StataMP 18 on 2026-10-07 (first 2026-09-26).
*! tests/unit/test_diagadvice.do: case D gives no note, no ratio and the
*! category counts; every flag its note in the practitioner wording; counts
*! never worded as Gamma; asymmetry is an internal-check report; rc 3498 on
*! bad input.
*!
*!     string scalar dwlsmv_diagadvice(real matrix G, real matrix C,
*!                                     real matrix B, string colvector ind,
*!                                     string colvector slab,
*!                                     string colvector pname
*!                                     | real matrix Cg, string scalar xname)
*!
*! CONDITIONAL FIT (decided 2026-11-15, the owner's; one binary covariate). Cg is
*! dwlsmv_countcells()'s table and xname the covariate's name; both or Cg alone (xname
*! then reads "x"). Given Cg, the notes change in four ways and nothing else:
*!   - the pattern note counts "(the set of responses together with the value of x)" and
*!     says "thresholds, slopes and correlations" (a pattern is a row of (y, x), and
*!     the slopes are statistics);
*!   - the centring note says "estimate" where it says "correlation" (an item row can
*!     name a threshold or a slope; slab carries the label, "y1~x" for a slope);
*!   - the pooled category counts are headed as pooled over both groups, and a block of
*!     counts WITHIN each value of x follows: the two group sizes, the smallest boundary
*!     and interior cell, and one sentence on what a slope's standard error rests on;
*!   - nothing is said about threshold aliasing, whose row is not evaluated (flag `.`).
*! G, C and B keep their shapes; a caller with no bread to report passes B as
*! (J(4, 2, .), J(4, 2, 0)), four rows not evaluated.
*!
*! G, C and B are dwlsmv_gammadiag()'s 7 x 4, dwlsmv_countdiag()'s 2 x 4 and
*! dwlsmv_breaddiag()'s 4 x 4. ind names the indicators in declaration order,
*! slab the sample statistics (map.label) and pname the free parameters in
*! theta order: the rows carry positions, and a note names what is there.
*! The same flags/prose split as dwlsmv_admissadvice(): a program branches on
*! the matrices, a person reads this. It neither prints nor aborts on a
*! condition; Module 11 decides what to show and when.
*!
*! TWO PARTS.
*!   Notes, one per flagged row, under its family's heading: the covariance of
*!   the sample statistics (Gamma), or the sandwich bread. Nothing when
*!   nothing flags.
*!   Category counts, ALWAYS: the smallest boundary and interior category
*!   counts, as descriptive diagnostics -- no threshold, no warning, not an
*!   admissibility criterion (decided 2026-10-07, the user's; minok is open).
*!   Until 2026-10-07 a values block printed Gamma's and the bread's
*!   singular-value ratios on every fit (decided 2026-09-26). They have no
*!   cutoff and no practitioner reading, so they are no longer printed:
*!   Module 11 stores them as e(gamma_ratio), e(gamma_ratio_raw),
*!   e(gamma_rank) and e(bread_ratio).
*!   Wording rewritten for practitioners 2026-10-07 (the user's).
*!
*! THE FAMILIES ARE NEVER MIXED. The category counts sit under their own
*! heading and no count line mentions Gamma, conditioning or singularity:
*! small-count failure is a separate family that routes to the sparse-cell
*! policy (spec Section 5). Admissibility is a separate channel altogether,
*! Module 7's dwlsmv_admissadvice(); Module 11 may show both.
*!
*! WHAT THE NOTES DO NOT SAY: what to do about an exactly singular Gamma or a
*! sparse table. Both policies are open (E-001); the
*! notes state the condition and what it touches, and stop.
*!
*! Row states: a flag of 1 gives a note; `.` (not evaluated) gives "not
*! evaluated" in the values block; .n values are printed as values.
*!
*! CONTRACT (rc 3498): G 7 x 4, C 2 x 4, B 4 x 4, and every non-zero location
*! inside its labels. Newlines are char(10) (a project coding rule).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Private to this file. Ratios in scientific notation, three decimals; the
// same format test_diagadvice.do reads back.
string scalar dwlsmv_diag_num(real scalar x)
{
    return(strtrim(strofreal(x, "%9.3e")))
}

// Private to this file: a label, or rc 3498 if the position is outside it.
string scalar dwlsmv_diag_lab(string colvector L, real scalar k, string scalar what)
{
    if (k < 1 | k > rows(L) | k != floor(k)) {
        errprintf("dwlsmv_diagadvice: %s position %f is outside its %f labels.\n", what, k, rows(L))
        _error("dwlsmv_diagadvice: a row names a position outside its labels")
    }
    return(L[k])
}

string scalar dwlsmv_diagadvice(real matrix G, real matrix C, real matrix B,
                                string colvector ind, string colvector slab,
                                string colvector pname,
                                | real matrix Cg, string scalar xname)
{
    string scalar nl, msg, gnote, bnote, where, xn, what
    real scalar r, cx, Kj, bmin, imin, bi, bc, bg, ii, ic, ig

    nl = char(10)
    cx = (args() >= 7)
    xn = "x"
    if (args() >= 8) {
        xn = xname
    }
    what = "correlation"
    if (cx) {
        what = "estimate"
    }
    if (rows(G) != 7 | cols(G) != 4 | rows(C) != 2 | cols(C) != 4 |
        rows(B) != 4 | cols(B) != 4) {
        errprintf("dwlsmv_diagadvice: expected G 7 x 4, C 2 x 4, B 4 x 4; got %f x %f, %f x %f, %f x %f.\n",
                  rows(G), cols(G), rows(C), cols(C), rows(B), cols(B))
        _error("dwlsmv_diagadvice: G, C or B is not a Module 10 row matrix")
    }
    // Every non-zero location must name something, flagged or not.
    for (r = 1; r <= 2; r++) {
        if (C[r, 3] != 0) (void) dwlsmv_diag_lab(ind, C[r, 3], "indicator")
    }
    if (G[1, 3] != 0) (void) dwlsmv_diag_lab(slab, G[1, 3], "statistic")
    if (G[1, 4] != 0) (void) dwlsmv_diag_lab(slab, G[1, 4], "statistic")
    if (G[4, 3] != 0) (void) dwlsmv_diag_lab(slab, G[4, 3], "statistic")
    if (G[4, 4] != 0) (void) dwlsmv_diag_lab(slab, G[4, 4], "statistic")
    if (G[5, 3] != 0) (void) dwlsmv_diag_lab(slab, G[5, 3], "statistic")
    if (B[2, 3] != 0) (void) dwlsmv_diag_lab(pname, B[2, 3], "parameter")
    // Category positions: countdiag's row 1 is a boundary category (position
    // >= 1), row 2 an interior one (>= 2); a location of 0 means none (D4).
    if (C[1, 3] != 0 & (C[1, 4] < 1 | C[1, 4] != floor(C[1, 4]))) {
        _error("dwlsmv_diagadvice: C's boundary row names an invalid category position")
    }
    if (C[2, 3] != 0 & (C[2, 4] < 2 | C[2, 4] != floor(C[2, 4]))) {
        _error("dwlsmv_diagadvice: C's interior row names an invalid category position")
    }
    // the per-group counts: dwlsmv_countcells()'s table, every item inside its labels
    if (cx) {
        if (cols(Cg) != 4 | rows(Cg) < 6) {
            errprintf("dwlsmv_diagadvice: Cg is %f x %f; expected the 4-column table of dwlsmv_countcells().\n",
                      rows(Cg), cols(Cg))
            _error("dwlsmv_diagadvice: Cg is not a dwlsmv_countcells() table")
        }
        if (missing(Cg) > 0 | Cg[1, 4] < 1 | Cg[2, 4] < 1) {
            _error("dwlsmv_diagadvice: Cg has missing values or an empty group")
        }
        for (r = 3; r <= rows(Cg); r++) {
            (void) dwlsmv_diag_lab(ind, Cg[r, 1], "indicator")
            if (Cg[r, 2] < 1 | Cg[r, 2] != floor(Cg[r, 2]) | (Cg[r, 3] != 0 & Cg[r, 3] != 1)) {
                _error("dwlsmv_diagadvice: Cg names an invalid category position or group")
            }
        }
    }

    // Wording 2026-10-07 (the user's): a technical term is kept where it
    // helps and explained at once. The notes on an exactly singular Gamma
    // state the condition and its consequence and stop: the remedy is an open
    // policy.

    // ---- notes: Gamma --------------------------------------------------------
    // Rows 1-2 cannot come from Stage 1's construction: they are internal
    // checks, and say so.
    gnote = ""
    if (G[1, 1] == 1) {
        gnote = gnote + "    An internal numerical check failed: Gamma, the sampling covariance, is not" + nl +
                "    symmetric (largest scale-free asymmetry " + dwlsmv_diag_num(G[1, 2]) + " at (" +
                slab[G[1, 3]] + ", " + slab[G[1, 4]] + ")). The remaining Gamma" + nl +
                "    checks were not run. Please report this result to the package developers." + nl
    }
    if (G[2, 1] == 1) {
        gnote = gnote + "    An internal numerical check failed: Gamma, the sampling covariance, is not" + nl +
                "    positive semidefinite (smallest/largest eigenvalue " + dwlsmv_diag_num(G[2, 2]) + ")." + nl +
                "    Please report this result to the package developers." + nl
    }
    if (G[3, 1] == 1) {
        if (cx) {
            gnote = gnote + "    Too few distinct response patterns: " + strofreal(G[3, 2]) +
                    " patterns (a pattern is the set of" + nl +
                    "    responses together with the value of " + xn + ") for " + strofreal(G[3, 3]) +
                    " thresholds, slopes and correlations." + nl
        }
        else {
            gnote = gnote + "    Too few distinct response patterns: " + strofreal(G[3, 2]) +
                    " patterns for " + strofreal(G[3, 3]) + " thresholds and correlations." + nl
        }
        gnote = gnote + "    Gamma is then exactly singular, so standard errors and the WLSMV test rest" + nl +
                "    on less information than the model uses." + nl
    }
    if (G[4, 1] == 1) {
        where = G[4, 4] == 0 ? slab[G[4, 3]] + " is constant (an empty category)" :
                slab[G[4, 3]] + " and " + slab[G[4, 4]] + " move together on every observation"
        gnote = gnote + "    " + where + nl +
                "    Gamma is then exactly singular (" + strofreal(G[4, 2]) + " dimension(s)), so standard" + nl +
                "    errors and the WLSMV test rest on less information than the model uses." + nl
    }
    if (G[5, 1] == 1) {
        gnote = gnote + "    The " + what + " " + slab[G[5, 3]] + " was not estimated precisely enough (its" + nl +
                "    estimating equation was not solved to tolerance: " + dwlsmv_diag_num(G[5, 2]) + ")." + nl +
                "    Gamma, the standard errors and the WLSMV test can be wrong even though" + nl +
                "    nothing else looks unusual." + nl
    }
    if (G[3, 1] == 1 | G[4, 1] == 1) {
        gnote = gnote + "    Point estimates are not affected; they use only Gamma's diagonal." + nl
    }

    // ---- notes: identification at the estimate (the bread) -------------------
    bnote = ""
    if (B[1, 1] == 1) {
        bnote = bnote + "    More free parameters (" + strofreal(B[1, 2]) + ") than sample statistics (" +
                strofreal(B[1, 3]) + "): the model is not identified." + nl
    }
    if (B[2, 1] == 1) {
        bnote = bnote + "    " + pname[B[2, 3]] + " has no effect on any fitted threshold or correlation at" + nl +
                "    the estimate (its column of the derivative matrix, which describes how the" + nl +
                "    model changes as parameters change, is zero), so its standard error" + nl +
                "    cannot be computed. " + strofreal(B[2, 2]) + " parameter(s) like this." + nl
    }
    if (B[3, 1] == 1) {
        bnote = bnote + "    The information in the fit is insufficient to determine all parameters" + nl +
                "    reliably (the derivative matrix has rank " + strofreal(B[3, 2]) + " of " +
                strofreal(B[3, 3]) + " columns at the" + nl +
                "    estimate). Some parameters are not identified there, and their standard" + nl +
                "    errors are not meaningful." + nl
    }

    msg = ""
    if (gnote != "") {
        msg = msg + "Sampling covariance of the thresholds and correlations (Gamma):" + nl + gnote
    }
    if (bnote != "") {
        msg = msg + "Identification of the parameters at the estimate:" + nl + bnote
    }

    // ---- category counts, on every fit (2026-10-07, the user's) ---------------
    // Descriptive: no threshold, no warning, not an admissibility criterion
    // (minok is open). The singular-value ratios are no longer printed; Module
    // 11 stores them in e(). "Not evaluated" is read from the FLAG, `.`.
    if (cx) {
        msg = msg + "Category counts (the fewest observations in any one category, both groups together):" + nl
    }
    else {
        msg = msg + "Category counts (the fewest observations in any one category):" + nl
    }
    if (C[1, 1] == .) {
        msg = msg + "    lowest or highest category: not evaluated" + nl
    }
    else {
        msg = msg + "    lowest or highest category: " + strofreal(C[1, 2]) + " (" + ind[C[1, 3]] +
                ", " + (C[1, 4] == 1 ? "lowest" : "highest") + ")" + nl
    }
    if (C[2, 1] == .) {
        msg = msg + "    middle category:            not evaluated" + nl
    }
    else if (C[2, 2] == .) {
        msg = msg + "    middle category:            none (every item is binary)" + nl
    }
    else {
        msg = msg + "    middle category:            " + strofreal(C[2, 2]) + " (" + ind[C[2, 3]] +
                ", category " + strofreal(C[2, 4]) + ")" + nl
    }

    // ---- conditional fit: the counts within each value of x (2026-11-15, the owner's) ----
    // Beside the pooled rows, never in place of them. Descriptive, like the pooled block: no
    // threshold, no warning. The last sentence says what a slope's standard error rests on.
    if (cx) {
        msg = msg + "Category counts within each value of " + xn +
                " (the fewest observations in any one category and group):" + nl
        msg = msg + "    group sizes:                " + xn + " = 0: " + strofreal(Cg[1, 4]) + ", " +
                xn + " = 1: " + strofreal(Cg[2, 4]) + nl
        bmin = .
        imin = .
        bi = 0
        bc = 0
        bg = 0
        ii = 0
        ic = 0
        ig = 0
        for (r = 3; r <= rows(Cg); r++) {
            Kj = max(select(Cg[., 2], Cg[., 1] :== Cg[r, 1]))
            if (Cg[r, 2] == 1 | Cg[r, 2] == Kj) {
                if (bmin == . | Cg[r, 4] < bmin) {
                    bmin = Cg[r, 4]
                    bi = Cg[r, 1]
                    bc = Cg[r, 2]
                    bg = Cg[r, 3]
                }
            }
            else if (imin == . | Cg[r, 4] < imin) {
                imin = Cg[r, 4]
                ii = Cg[r, 1]
                ic = Cg[r, 2]
                ig = Cg[r, 3]
            }
        }
        msg = msg + "    lowest or highest category: " + strofreal(bmin) + " (" + ind[bi] + ", " +
                (bc == 1 ? "lowest" : "highest") + ", " + xn + " = " + strofreal(bg) + ")" + nl
        if (imin == .) {
            msg = msg + "    middle category:            none (every item is binary)" + nl
        }
        else {
            msg = msg + "    middle category:            " + strofreal(imin) + " (" + ind[ii] +
                    ", category " + strofreal(ic) + ", " + xn + " = " + strofreal(ig) + ")" + nl
        }
        msg = msg + "    The standard error of a slope (the effect of " + xn + " on an item) depends on the" + nl +
                "    counts within each value of " + xn + ". A category with many observations overall" + nl +
                "    can have few in one group, which the pooled counts above do not show." + nl
    }
    return(msg)
}

end
