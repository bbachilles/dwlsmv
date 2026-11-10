*! dwlsmv -- DIFFTEST validation check 6, the comparison: a refit reproduces the
*! stored fit (Module 12)
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9; the missing-T3 and df
*! rules rewritten by the main session (user, 2026-10-08). Executed 2026-10-08:
*! tests/unit/test_dtrefit.do 22/22, mutants R1-R12 (R5, R8 as sharpened) killed.
*!
*!     void dwlsmv_dtrefit(real matrix b0, real scalar T0, real matrix bound0,
*!                         real scalar df0,
*!                         real matrix b1, real scalar T1, real matrix bound1,
*!                         real scalar df1,
*!                         real scalar tol, string scalar name)
*!
*! Spec Module 12, Validation contract, check 6. (b0, T0, bound0, df0) are the
*! STORED e(b), e(T3), e(bound) and e(df_T3); (b1, T1, bound1, df1) the
*! refit's. name labels the messages. Returns silently, or aborts rc 3498.
*!
*! THE COMPARISON (decided, user 2026-10-08):
*!   - e(b) and e(T3): the measure m(x, y) = |x - y| / (1 + |y|), y the stored
*!     value, maximum over elements, must be <= tol;
*!   - e(bound), a 0/1 flag vector, must be IDENTICAL: a discrete status, never
*!     loosened by the tolerance;
*!   - e(df_T3) must be equal: the refit is the same model;
*!   - a missing e(T3) in either result is refused. dwlsmv never posts one: a
*!     saturated fit (D = 0) posts T3 = T_raw, finite (dwls_inference, spec
*!     4.9, decided 2026-09-24), and the measure compares it near 0 as an
*!     absolute difference. Whether a DIFFTEST against a saturated model is
*!     meaningful is the front end's, not this check's.
*! tol is an ARGUMENT; the caller passes dwlsmv_dt_reftol(), which returns
*! missing until the class (k) measurement sets it, and this routine refuses a
*! missing tolerance.
*!
*! WHAT IT DOES NOT DO. It does not refit anything, read e(), compare
*! e(converged) (check 7), the data (checks 2-3) or the version (check 1). It
*! does not choose a tolerance.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_dtrefit(real matrix b0, real scalar T0, real matrix bound0,
                    real scalar df0,
                    real matrix b1, real scalar T1, real matrix bound1,
                    real scalar df1,
                    real scalar tol, string scalar name)
{
    real scalar gb, gt, flip

    if (missing(tol)) {
        errprintf("dwlsmv check 6 (%s): the refit-reproduction tolerance has not been set.\n", name)
        errprintf("dwlsmv_dt_reftol() returns missing until the class (k) measurement is done.\n")
        _error("dwlsmv check 6: the reproduction tolerance is not set")
    }
    if (tol <= 0) {
        errprintf("dwlsmv check 6 (%s): tolerance %g; it must be positive.\n", name, tol)
        _error("dwlsmv check 6: the tolerance must be positive")
    }
    if (rows(b0) != rows(b1) | cols(b0) != cols(b1) |
        rows(bound0) != rows(bound1) | cols(bound0) != cols(bound1) |
        rows(b0) * cols(b0) == 0) {
        errprintf("dwlsmv check 6 (%s): e(b) is %g x %g stored and %g x %g refit; e(bound) %g x %g and %g x %g.\n",
                  name, rows(b0), cols(b0), rows(b1), cols(b1),
                  rows(bound0), cols(bound0), rows(bound1), cols(bound1))
        _error("dwlsmv check 6: stored and refit results are not conformable")
    }
    if (missing(b0) > 0 | missing(b1) > 0 | missing(bound0) > 0 | missing(bound1) > 0) {
        errprintf("dwlsmv check 6 (%s): a stored or refit e(b) or e(bound) value is missing.\n", name)
        _error("dwlsmv check 6: a stored or refit value is missing")
    }
    if (missing(T0) | missing(T1)) {
        errprintf("dwlsmv check 6 (%s): e(T3) is %g stored and %g refit. dwlsmv posts a T3 for\n",
                  name, T0, T1)
        errprintf("every fit, a saturated one included (T3 = T_raw), so a missing one is malformed.\n")
        _error("dwlsmv check 6: e(T3) is missing")
    }
    if (missing(df0) | missing(df1)) {
        errprintf("dwlsmv check 6 (%s): e(df_T3) is %g stored and %g refit.\n", name, df0, df1)
        _error("dwlsmv check 6: e(df_T3) is missing")
    }
    if (df0 != df1) {
        errprintf("dwlsmv check 6 (%s): e(df_T3) is %g stored and %g refit; the refit is not\n",
                  name, df0, df1)
        errprintf("the same model.\n")
        _error("dwlsmv check 6: the refit has different degrees of freedom")
    }

    gb   = max(abs(b1 :- b0) :/ (1 :+ abs(b0)))
    gt   = abs(T1 - T0) / (1 + abs(T0))
    flip = any(bound1 :!= bound0)

    if (gb > tol | gt > tol | flip) {
        errprintf("dwlsmv check 6 (%s): the refit does not reproduce the stored fit.\n", name)
        errprintf("    e(b): worst difference %g; e(T3): %g; tolerance %g;\n", gb, gt, tol)
        if (flip) {
            errprintf("    e(bound) differs.\n")
        }
        else {
            errprintf("    e(bound) identical.\n")
        }
        errprintf("Estimation is deterministic, so a larger difference means the data, the\n")
        errprintf("specification, the version or the machine's arithmetic is not the one that\n")
        errprintf("produced the stored result. Nothing is repaired or substituted.\n")
        _error("dwlsmv check 6: the refit does not reproduce the stored fit")
    }
}

end
