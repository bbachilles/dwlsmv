*! dwlsmv -- classify empty categories by failure mode (Module 1, spec 4.5.6)
*!
*! Executed under StataMP 18 on 2026-09-19: tests/unit/test_thresholds.do
*! asserts all four return codes on hand-built count vectors, so the boundary
*! and interior modes are checked as DISTINGUISHED, not merely as caught.
*!
*! Spec Section 4.5.6 is explicit that "Module 1's category validation must
*! distinguish empty INTERIOR categories from empty BOUNDARY categories rather
*! than treating 'degenerate category' as one class". They are not one class
*! because they fail differently:
*!
*!   BOUNDARY (first or last category empty)   q-hat is 0 or 1, so tau-hat is
*!       +/- infinity and phi(tau-hat) is 0. Loud. It kills that threshold
*!       outright and nothing downstream can proceed past it.
*!
*!   INTERIOR (some middle category empty)     q-hat[c] == q-hat[c+1], so two
*!       thresholds collide at a FINITE value and nothing overflows. Two psi
*!       columns become identical and Gamma is rank-deficient -- silently. DWLS
*!       point estimation is unaffected, because W_D = diag(Gamma)^-1 never
*!       inverts the full matrix, so Module 7 converges and returns a plausible
*!       theta-hat. Modules 8 and 9 are the ones that are wrong. That is the
*!       "runs without error but is quietly wrong" mode, and it is why this
*!       distinction is carried rather than collapsed.
*!
*! Return codes: 0 none, 1 boundary only, 2 interior only, 3 both.
*!
*! This routine DECIDES NOTHING. It reports. Whether a rank-deficient Gamma
*! blocks estimation entirely or only suppresses the robust corrections is an
*! open item (E-054) and belongs to Module 10, not here. The
*! near-boundary third mode of Section 4.5.6 -- a category present but very
*! sparse -- is deliberately absent: it has no threshold at which it becomes a
*! different object, so it is a judgement about counts that the caller makes
*! from the count vector itself, not a classification this routine can make.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_catdegen(real colvector n)
{
    real scalar K, c, boundary, interior

    if (cols(n) != 1) {
        _error("dwlsmv_catdegen: n must be a column vector")
    }
    K = rows(n)
    if (K < 2) {
        _error("dwlsmv_catdegen: n must name at least two categories")
    }

    boundary = (n[1] == 0 | n[K] == 0)

    interior = 0
    for (c = 2; c <= K - 1; c++) {
        if (n[c] == 0) {
            interior = 1
        }
    }

    return(boundary + 2 * interior)
}

end
