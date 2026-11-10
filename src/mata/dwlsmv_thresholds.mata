*! dwlsmv -- Stage 1 thresholds for one variable (Module 1, spec Section 4.5.1)
*!
*! Executed under StataMP 18 on 2026-09-19. Its output on case_D.csv agrees
*! with lavaan's `th` to 3.3e-15 across all ten thresholds -- the same figure
*! the R-side comparison found, now measured from inside Stata. Both
*! degenerate-category refusals are exercised.
*!
*!     tau-hat[j,c] = invnormal( q-hat[j,c] )
*!
*! and NOTHING more. That is not an inference from the spec, it is a
*! measurement: lavaan's `th` was compared against qnorm(empirical cumulative
*! proportion) computed directly from case_D.csv and agrees to 3.3e-15 across
*! all ten thresholds, with Mplus agreeing to every digit it prints. See
*! E-072, 2026-09-18. A first stage that needs an optimizer, a
*! starting value or an iteration is wrong, not merely slower.
*!
*! The ordering constraint tau[1] < ... < tau[C] holds automatically, with no
*! constrained optimization, because q-hat is a cumulative proportion and
*! invnormal is monotone. It holds STRICTLY provided no interior category is
*! empty -- which is exactly the degeneracy classified below.
*!
*! WHY THIS ABORTS ON AN EMPTY CATEGORY. An empty boundary category makes the
*! threshold infinite and an empty interior category makes two thresholds
*! collide; in neither case is there a number to return that is not a lie. Spec
*! Section 4.5.6 and Module 2's stated v1 default agree on the behaviour: fail
*! or warn, do not silently repair. The error names the mode so the caller can
*! tell the loud failure from the silent one.
*!
*! This is NOT the sparse-cell policy decision, which is open (minimum cell
*! counts, collapse vs error vs warn, continuity corrections -- see
*! E-001). A caller that wants to decide for itself calls
*! dwlsmv_catcounts() and dwlsmv_catdegen() first; they are separate routines
*! precisely so that is possible without catching an error.
*!
*! Requires dwlsmv_catcounts.mata, dwlsmv_catdegen.mata and dwlsmv_cumprop.mata
*! to be loaded. No definition file here issues `mata clear`, so order does not
*! matter.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_thresholds(real colvector y, real scalar K)
{
    real colvector n
    real scalar d, c
    string scalar msg, empty

    n = dwlsmv_catcounts(y, K)
    d = dwlsmv_catdegen(n)

    // The message names the offending categories and what to do about them.
    // "A boundary category is empty" tells a user what happened; it does not
    // tell them which category, or that collapsing is the usual remedy, and
    // an error they cannot act on is only marginally better than none.
    if (d) {
        empty = ""
        for (c = 1; c <= K; c++) {
            if (n[c] == 0) {
                if (empty != "") {
                    empty = empty + ", "
                }
                empty = empty + strofreal(c - 1)
            }
        }

        // The guidance goes through errprintf, NOT through _error. Mata's
        // _error() accepts at most 100 characters and raises rc 3300
        // "argument out of range" past that -- so a long, helpful message is
        // silently replaced by a useless one that points at _error() instead
        // of at this function. Measured 2026-09-19: 100 characters works, 101
        // does not. A project coding rule.
        errprintf("dwlsmv_thresholds: no observations in categor%s%s of %s" +
                  " (N = %s).\n",
                  (strpos(empty, ",") ? "ies " : "y "), empty,
                  strofreal(K), strofreal(rows(y)))

        if (d == 1) {
            errprintf("  This is a BOUNDARY category, so its threshold is" +
                      " infinite and cannot be estimated.\n")
            errprintf("  Collapse it into the adjacent category, or drop" +
                      " this indicator from the model.\n")
            msg = "dwlsmv_thresholds: empty boundary category -- collapse it or drop the indicator"
        }
        else if (d == 2) {
            errprintf("  This is an INTERIOR category, so two thresholds" +
                      " collide at the same finite value\n")
            errprintf("  and Gamma becomes rank-deficient. That leaves the" +
                      " point estimates correct while the\n")
            errprintf("  standard errors and the test statistic are wrong," +
                      " which is why it is refused here\n")
            errprintf("  rather than warned about. Collapse it into an" +
                      " adjacent category.\n")
            msg = "dwlsmv_thresholds: empty interior category -- collapse it into a neighbour"
        }
        else {
            errprintf("  Both boundary and interior categories are empty." +
                      " Collapse them into their\n")
            errprintf("  neighbours, or drop this indicator from the model.\n")
            msg = "dwlsmv_thresholds: empty boundary and interior categories"
        }

        _error(msg)
    }

    return(invnormal(dwlsmv_cumprop(n)))
}

end
