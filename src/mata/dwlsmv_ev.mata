*! dwlsmv -- Module 8: e(V) from ACOV
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_acov.do checks 2, 4 and 6 reject /N here: replacing N - 1 with N fails all three (mutation, 2026-09-24). Check 2b does not call this routine; it shows the fixture can tell the two divisors apart.
*!
*! Spec Section 4.10, decided 2026-09-23:
*!
*!     e(V) = ACOV / (N - 1)
*!
*! and this routine implements literally that and nothing else. ACOV is
*! dwlsmv_acov()'s output, the asymptotic covariance of sqrt(N)(theta-hat -
*! theta_0); this is the one place a sample-size divisor enters the SEs.
*!
*! N - 1, NOT N, IS A COMPATIBILITY CONVENTION AND NOT A RESULT -- the same one
*! as T = (N-1) F in Section 4.9. lavaan 0.7.2 divides by ntotal - ngroups in
*! lav_model_vcov() with an N-scaled Gamma. /N would miss lavaan's SEs by
*! exactly sqrt(N/(N-1)) - 1, 1.0e-4 at N = 5000, which was INSIDE Section
*! 7.2's SE tier until 2026-09-24 (the current 1e-5 tier sees it by 9.2-9.7x, a one-off R figure), so
*! test_acov.do requires N - 1 and rejects N explicitly regardless.
*!
*! It takes N rather than a divisor so that the convention lives in exactly
*! one line. Its caller is the code that posts e(V) -- Module 11 -- and not
*! dwls_inference(): the fit stores ACOV, never e(V), because e(V) is
*! derivable from acov and N and a stored copy could disagree.
*!
*! G > 1 (decided 2026-10-09, check 10; spec D8d): e(V) = ACOV / (N - G),
*! lavaan's ntotal - ngroups, through the optional third argument; omitted,
*! G = 1 and the divisor is N - 1 as before (2026-10-12, invariance() posting;
*! test_invariance_post.do check 5, mutant M-b red).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_ev(real matrix acov, real scalar N, | real scalar G)
{
    real scalar ng

    if (rows(acov) != cols(acov)) {
        _error("dwlsmv_ev: acov is not square")
    }
    if (missing(N) | N < 2 | N != floor(N)) {
        errprintf("dwlsmv_ev: N = %f; the e(V) divisor N - 1 needs an integer N >= 2\n", N)
        _error("dwlsmv_ev: N must be an integer >= 2")
    }
    // an omitted optional argument is never read
    ng = 1
    if (args() >= 3) {
        ng = G
    }
    if (missing(ng) | ng < 1 | ng != floor(ng) | ng >= N) {
        errprintf("dwlsmv_ev: G = %f; the number of groups must be an integer in 1..N-1\n", ng)
        _error("dwlsmv_ev: G must be an integer >= 1 and < N")
    }

    return(acov :/ (N - ng))
}

end
