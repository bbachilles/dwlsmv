*! dwlsmv -- conditional covariate, one 0/1 x, one group: the fit and inference acceptance check
*!
*! Executed 2026-11-13 by the main session (Stata 18 MP), mutant-tested. Drafted 2026-11-13 by a read-only agent; the main session loads and
*! runs it (tests/unit/test_condx_fit.do) and removes this banner only then.
*!
*!     real scalar dwlsmv_cx1_fit_check(struct dwlsmv_stage1 vector sts, string scalar who)
*!
*! Replaces the dwlsmv_s2_guard_st() call in dwlsmv_fit_stack() (for dwls_fit(), through
*! dwlsmv_fit_uls() and dwlsmv_fit_from(), which pass the acceptance flag) and in
*! dwls_inference(). It returns q, 0 or 1, as dwlsmv_cx1_check() does for the maps:
*!   - several Stage 1 results: dwlsmv_s2_guard_st() decides, so a conditional map is refused at
*!     rc 3498 with its group number, as before (several groups with a covariate are not
*!     specified);
*!   - one result: q is read from its map (dwlsmv_statmap_q); q > 1 is refused (several
*!     covariates are not specified, spec 4.8); q = 0 and q = 1 pass. An element never filled
*!     (D, T, P or p absent or missing) has no q and is left to the caller's own checks, as in
*!     dwlsmv_s2_guard(). A model table with gamma or kappa rows is checked against the map by
*!     dwlsmv_cx1_check(), which dwlsmv_fit_uls() and dwlsmv_fit_from() call as well, since this
*!     function does not see the table.
*! The diagnostics that rebuild the stack (dwlsmv_sens*) call dwlsmv_fit_stack() without the flag
*! and keep the old refusal. An unconditional result (q = 0) is a no-op, as the guard was.
*!
*! Refuses with rc 3498: errprintf() carries the guidance, _error() a short literal.
*!
*! Requires dwlsmv_types, dwlsmv_statmap_q, dwlsmv_s2_guard.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_cx1_fit_check(struct dwlsmv_stage1 vector sts, string scalar who)
{
    struct dwlsmv_statmap scalar m
    real scalar q, built

    if (length(sts) < 1) {
        return(0)
    }

    if (length(sts) > 1) {
        // several groups: a conditional map stays refused (not specified yet)
        dwlsmv_s2_guard_st(sts, who)
        return(0)
    }

    m = sts[1].map
    q = 0
    built = (length(m.D) == 1) + (length(m.T) == 1) +
            (length(m.P) == 1) + (length(m.p) == 1)
    if (built == 4) {
        if (missing(m.D) + missing(m.T) + missing(m.P) + missing(m.p) == 0) {
            q = dwlsmv_statmap_q(m)
        }
    }

    if (q > 1) {
        errprintf("%s: the index map has q = %g covariates; Stage 2 accepts one 0/1 covariate (spec 4.8).\n", who, q)
        _error("several covariates are not specified yet")
    }

    return(q)
}

end
