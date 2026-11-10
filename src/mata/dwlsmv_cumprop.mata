*! dwlsmv -- cumulative sample proportions from a category count vector
*! (Module 1, spec Section 4.5.1)
*!
*! Executed under StataMP 18 on 2026-09-19: tests/unit/test_thresholds.do and
*! tests/unit/test_loading.do both drive it, the latter on a hand-built count
*! vector whose answer is obvious by inspection.
*!
*!     q[c] = P-hat(Y < c) = (n[1] + ... + n[c]) / N,   c = 1 .. K-1
*!
*! Takes COUNTS, not data. Two reasons, both about avoiding a second pass over
*! something that is already known. dwlsmv_thresholds() needs the counts anyway
*! for the degeneracy classification, so passing them in means the data is
*! walked once; and a function of counts alone is a pure function that the test
*! can drive with hand-built vectors where the right answer is obvious.
*!
*! The running total is built by hand rather than with runningsum(), matching
*! dwlsmv_statmap_build()'s prefix sum: the file then leans on as little as
*! possible, which is worth more at the version floor than the line it saves.
*!
*! Precision note. n and N are exact integers in double, so each n/N is the
*! correctly-rounded proportion and the sequence is bitwise reproducible across
*! processor counts. This is one of the accumulations the project's quadcross()
*! rule does NOT apply to -- there is no reduction whose order can vary.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_cumprop(real colvector n)
{
    real colvector q
    real scalar K, c, N, run

    if (cols(n) != 1) {
        _error("dwlsmv_cumprop: n must be a column vector")
    }
    K = rows(n)
    if (K < 2) {
        _error("dwlsmv_cumprop: n must name at least two categories")
    }
    if (min(n) < 0 | sum(n :!= floor(n)) > 0) {
        _error("dwlsmv_cumprop: counts must be non-negative integers")
    }
    N = sum(n)
    if (N < 1) {
        _error("dwlsmv_cumprop: counts are all zero")
    }

    q   = J(K - 1, 1, 0)
    run = 0
    for (c = 1; c <= K - 1; c++) {
        run  = run + n[c]
        q[c] = run / N
    }

    return(q)
}

end
