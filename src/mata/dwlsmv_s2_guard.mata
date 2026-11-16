*! dwlsmv -- Stage 2's refusal of a conditional index map (spec 5.1)
*!
*!     void dwlsmv_s2_guard(struct dwlsmv_statmap vector maps, string scalar who)
*!     void dwlsmv_s2_guard_st(struct dwlsmv_stage1 vector sts, string scalar who)
*!
*! Stage 2 (Modules 3 to 9) indexes the correlations at map.T + k, which is
*! wrong once slopes exist (spec 4.2, 5.1). A conditional map (q > 0, from
*! dwls_sample_stats_cx()) must therefore never reach it. Every Stage 2 entry
*! point that receives a map, a vector of maps, or a Stage 1 result or vector
*! of them calls one of these two functions before it reads the map; the
*! unconditional path (q = 0) is a no-op and is unchanged.
*!
*! Callers (2026-11-12): dwls_identify, dwlsmv_cx1_check (for dwls_model and dwls_delta, which
*! accept one group with one covariate and pass several maps here),
*! dwlsmv_fit_stack (every caller that does not pass its cxok flag, i.e. the dwlsmv_sens*
*! diagnostics; dwls_fit's own paths pass it and use dwlsmv_cx1_fit_check(), which falls back
*! to dwlsmv_s2_guard_st() for several groups), dwlsmv_cx1_fit_check (several groups, for
*! dwls_fit and dwls_inference), dwlsmv_sens_hess, dwlsmv_mg_admissible,
*! dwlsmv_fitidx_base.
*!
*! q is read from the map by dwlsmv_statmap_q(), never recomputed here. An
*! element that was never filled (D, T, P or p absent or missing) has no q and
*! is left for the caller's own checks, which refuse it as before; a
*! conditional map always has all four.
*!
*! Refuses with rc 3498: errprintf() carries the guidance, _error() a short
*! literal.
*!
*! Drafted by a read-only agent 2026-11-12; loaded and run by the main session
*! the same day (tests/unit/test_s2_guard.do).
*!
*! Requires dwlsmv_types, dwlsmv_statmap_q.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_s2_guard(struct dwlsmv_statmap vector maps, string scalar who,
                     | real scalar qok)
{
    real scalar g, q, built, qx

    // qx: the number of covariates the caller's specification declares. Only
    // dwls_identify() passes it (spec 5.1); every other entry point passes
    // two arguments, so qx is 0 and a conditional map is refused as before.
    if (args() < 3) {
        qx = 0
    }
    else {
        qx = qok
    }

    for (g = 1; g <= length(maps); g++) {
        built = (length(maps[g].D) == 1) + (length(maps[g].T) == 1) +
                (length(maps[g].P) == 1) + (length(maps[g].p) == 1)
        if (built == 4) {
            if (missing(maps[g].D) + missing(maps[g].T) +
                missing(maps[g].P) + missing(maps[g].p) == 0) {
                q = dwlsmv_statmap_q(maps[g])
                if (q > 0 & qx == 0) {
                    errprintf("%s: the index map of group %g is conditional (q = %g).\n",
                              who, g, q)
                    errprintf("Stage 2 does not yet accept a conditional Stage 1 result: its slope rows\n")
                    errprintf("shift the correlations away from map.T + k (spec 5.1).\n")
                    _error("a conditional index map is not accepted by Stage 2 yet")
                }
                if (q != qx & qx > 0) {
                    errprintf("%s: the index map of group %g has q = %g covariate(s), but the model\n",
                              who, g, q)
                    errprintf("specification declares %g. Build the map with dwls_sample_stats_cx() for exactly\n",
                              qx)
                    errprintf("the declared covariates (spec 5.1, Conditional covariates).\n")
                    _error("the index map and the declared covariates disagree")
                }
            }
        }
    }
}

// The same refusal for Stage 1 results: only each result's map is read. The
// maps are copied into a vector so the group number in the message is right.
void dwlsmv_s2_guard_st(struct dwlsmv_stage1 vector sts, string scalar who,
                        | real scalar qok)
{
    struct dwlsmv_statmap vector maps
    real scalar g, qx

    if (args() < 3) {
        qx = 0
    }
    else {
        qx = qok
    }

    if (length(sts) > 0) {
        maps = dwlsmv_statmap(length(sts), 1)
        for (g = 1; g <= length(sts); g++) {
            maps[g] = sts[g].map
        }
        dwlsmv_s2_guard(maps, who, qx)
    }
}

end
