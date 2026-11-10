*! dwlsmv -- turn dwlsmv_admissible()'s flags into notes a user can act on
*! (Module 7, spec Section 5)
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_fit.do check 12 prints the rendered notes to the log and asserts their content; the first run printed -.21, which is why values go through dwlsmv_admiss_num().
*!
*!     string scalar dwlsmv_admissadvice(real matrix A,
*!                                       string colvector ind,
*!                                       string colvector fac)
*!
*! A is dwlsmv_admissible()'s 6 x 4 result; ind and fac are the indicator and
*! factor names in declaration order, which Module 11 holds and the parameter
*! table does not. Returns a multi-line note, or "" when every condition
*! passes. The caller decides whether and at what severity to show it -- this
*! routine neither prints nor aborts, the same split as dwlsmv_diagadvice(): a
*! program branches on the flags, a person reads this. Wording rewritten for
*! practitioners 2026-10-07 (the user's): no matrix notation, each technical
*! term explained where it appears.
*!
*! One line per failed condition, naming the indicator or factor pair and the
*! value, then what it means for the output. Point estimates are still the
*! DWLS minimum and are still reported -- the spec asks for a diagnostic, not a
*! refusal -- but Section 4.10's sandwich and Section 4.9's statistic are
*! asymptotics about an interior point, and the note says so.
*!
*! Values are written by dwlsmv_admiss_num(), a helper private to this file
*! and defined above the routine that uses it.
*!
*! Newlines are char(10). A backslash-n inside a Mata string is two characters,
*! not an escape; an escape only inside a printf FORMAT. A project coding rule.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Private to this file: how a value is written into a note. Four decimals
// with a leading zero -- Stata's %g writes -.21 -- and scientific notation
// below 1e-4, where four decimals would print an eigenvalue of -1e-6 as
// -0.0000 and hide the very thing the note is about.
string scalar dwlsmv_admiss_num(real scalar x)
{
    if (abs(x) >= 1e-4) {
        return(strtrim(strofreal(x, "%12.4f")))
    }
    return(strtrim(strofreal(x, "%10.2e")))
}

string scalar dwlsmv_admissadvice(real matrix A,
                                  string colvector ind,
                                  string colvector fac)
{
    string scalar msg, nl

    if ((rows(A) != 6 & rows(A) != 7) | cols(A) != 4) {
        _error("dwlsmv_admissadvice: A is not a 6 x 4 or 7 x 4 admissibility result")
    }
    if (!any(A[., 1] :== 1)) {
        return("")
    }

    // Wording 2026-10-07 (the user's): a technical term is kept where it
    // helps and explained at once; no matrix notation.
    nl  = char(10)
    msg = "The fitted solution is not admissible (some estimates are outside the range" + nl +
          "a valid model allows):" + nl

    if (A[5, 1] == 1) {
        msg = msg + "    the structural paths form a feedback loop with no solution, so the" + nl +
              "    factor variances and covariances, and everything checked from them," + nl +
              "    could not be computed" + nl
    }
    if (A[1, 1] == 1) {
        msg = msg + "    " + ind[A[1, 3]] + ": residual variance " +
              dwlsmv_admiss_num(A[1, 2]) + ", below 0 (a Heywood case: the model" + nl +
              "    implies the factors explain more than all of " + ind[A[1, 3]] +
              "'s variance)" + nl
    }
    if (A[2, 1] == 1) {
        msg = msg + "    " + fac[A[2, 3]] + " <-> " + fac[A[2, 4]] +
              ": implied factor correlation " +
              dwlsmv_admiss_num(A[2, 2]) + ", outside -1 to 1" + nl
    }
    if (A[3, 1] == 1) {
        msg = msg + "    the factors' implied variances and covariances are impossible" + nl +
              "    (their matrix is not positive definite; smallest eigenvalue " +
              dwlsmv_admiss_num(A[3, 2]) + ")" + nl
    }
    if (A[4, 1] == 1) {
        msg = msg + "    the implied correlations among the items' underlying continuous" + nl +
              "    responses are impossible (their matrix is not positive definite;" + nl +
              "    smallest eigenvalue " + dwlsmv_admiss_num(A[4, 2]) + ")" + nl
    }
    if (A[6, 1] == 1) {
        msg = msg + "    " + ind[A[6, 3]] +
              ": thresholds are not in increasing order (smallest gap " +
              dwlsmv_admiss_num(A[6, 2]) + ")" + nl
    }
    if (rows(A) == 7) {
        if (A[7, 1] == 1) {
            msg = msg + "    " + ind[A[7, 3]] + ": scale factor " +
                  dwlsmv_admiss_num(A[7, 2]) + ", not above 0 (the scale of the item's" + nl +
                  "    underlying continuous response must be positive)" + nl
        }
    }

    msg = msg + "  Point estimates are reported. Standard errors and the WLSMV test assume a" + nl +
          "  solution inside the allowed range and should not be relied on." + nl

    return(msg)
}

end
