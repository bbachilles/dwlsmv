*! dwlsmv -- conditional-covariate mode: Stage 1 as a whole, one binary x
*! (Modules 1-4; spec 5.1, 4.2 conditional paragraph, 4.8 *Conditional
*! covariates*; derivation @sec-condx, @sec-condx-norm, @sec-condx-hess,
*! @sec-condx-policy)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*! tests/unit/test_condx_stage1.do is its development test. The
*! recorded link 2 run (prereg.md, docs/notes/condx-validation) needs its own
*! do-file, hashed in a dated addendum first ([V1]).
*!
*! TAKES (Y, K, X) AND NOTHING MODEL-IMPLIED (spec 4.8: "X is data, and
*! nothing model-implied enters"). No model specification, no starting
*! values, no theta. This is a NEW entry point beside dwls_sample_stats(),
*! which is unchanged: a model without covariates never reaches this file.
*!
*! SCOPE: one covariate, coded 0/1, one group (@sec-condx). Anything else is
*! refused at rc 3498 by dwlsmv_cx_check().
*!
*! The order of operations is the two-step estimator:
*!
*!     0. refusals: X's shape, missing, coding; Y's codes and pooled
*!        categories; separation by x (Result P1)
*!     1. each item's thresholds and slope, by ordered probit on x (Module 1)
*!     2. each pair's residual rho given x, with the item fits HELD FIXED
*!        (Module 2), by v1's root-finder on the pooled score
*!     3. Gamma, in the Hessian form, assembled from 1 and 2 (Module 4)
*!
*! s = (thresholds item-major at x = 0, slopes covariate-major, residual
*! correlations lower-triangle), D = T + p + P (spec 4.2).
*!
*! WHICH PAIRS ENTER (decided 2026-10-25, P2): status 1, and status 6 reached
*! by the endpoint check (five or more evaluations), as in v1. A status-6
*! pair that enters is named in a printed note (the user, 2026-11-13): it is
*! an endpoint-sign-certification limitation, not ordinary convergence. v1
*! prints no such note. Every other
*! pair aborts Stage 1 at rc 3498 after all pairs are classified and
*! reported. There is no one-empty-cell repair and no keep_bound in this mode:
*! lavaan applies no repair when x is present, and v1's per-table
*! classification does not transfer (a table can look like a boundary in one
*! group and not be one pooled). stage1.bound is therefore all zeros.
*!
*! WHAT IS NOT DONE HERE. Listwise deletion (the caller deletes, jointly over
*! Y and x, and hands over a matrix whose row count is the N it will report).
*! The (A1) per-item diagnostic (P4, a non-gating reported table) and Module
*! 10's rank flag (P3): output-layer work, not this build.
*!
*! Requires dwlsmv_types, dwlsmv_statmap_build, dwlsmv_statmap_q,
*! dwlsmv_pos_tau, dwlsmv_pos_slope, dwlsmv_pos_rho, dwlsmv_catcounts,
*! dwlsmv_catdegen, dwlsmv_cumprop, dwlsmv_thresholds, dwlsmv_cellprob,
*! dwlsmv_dcellprob_rho, dwlsmv_wmat, dwlsmv_crosstab, dwlsmv_pairinfo,
*! dwls_polychoric, dwls_sample_stats (dwlsmv_ss_statusword and
*! dwlsmv_ss_readings), dwlsmv_cx_check, dwlsmv_cx_item, dwlsmv_cx_pair,
*! dwlsmv_cx_psi_matrix, dwlsmv_cx_psibar and dwls_gamma_cx.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

struct dwlsmv_stage1 scalar dwls_sample_stats_cx(real matrix Y, real colvector K,
                                                 real matrix X)
{
    struct dwlsmv_stage1 scalar st
    real matrix    R
    real colvector th, thi, thj, cls
    real rowvector r
    real scalar    j, c, i, k, m, pos, it, nbad, nbad2
    string scalar  line

    st.map = dwlsmv_statmap_build(K, 1)

    if (cols(Y) != st.map.p) {
        errprintf("dwls_sample_stats_cx: Y has %s columns but K names %s variables.\n",
                  strofreal(cols(Y)), strofreal(st.map.p))
        _error("dwls_sample_stats_cx: Y and K disagree on the number of variables")
    }
    st.N = rows(Y)
    if (st.N < 2) {
        _error("dwls_sample_stats_cx: Y needs at least two observations")
    }

    // ---- 0. refusals ---------------------------------------------------------
    dwlsmv_cx_check(Y, K, X)

    st.s = J(st.map.D, 1, .)

    // ---- 1. items: thresholds at x = 0 and slopes ---------------------------
    for (j = 1; j <= st.map.p; j++) {
        th = dwlsmv_cx_item_fit(Y[., j], X, st.map.K[j], j, it)
        for (c = 1; c <= st.map.C[j]; c++) {
            st.s[dwlsmv_pos_tau(st.map, j, c)] = th[c]
        }
        st.s[dwlsmv_pos_slope(st.map, j, 1)] = th[st.map.K[j]]
    }

    // ---- 2. residual correlations, items held fixed --------------------------
    // R: one row per pair in the map's order: (i, j, rho-hat, status, evals).
    R   = J(st.map.P, 5, .)
    cls = J(st.map.P, 1, 0)
    m   = 0
    for (i = 1; i <= st.map.p - 1; i++) {
        for (k = i + 1; k <= st.map.p; k++) {
            m   = m + 1
            pos = dwlsmv_pos_rho(st.map, i, k)
            thi = J(st.map.K[i], 1, .)
            thj = J(st.map.K[k], 1, .)
            for (c = 1; c <= st.map.C[i]; c++) {
                thi[c] = st.s[dwlsmv_pos_tau(st.map, i, c)]
            }
            for (c = 1; c <= st.map.C[k]; c++) {
                thj[c] = st.s[dwlsmv_pos_tau(st.map, k, c)]
            }
            thi[st.map.K[i]] = st.s[dwlsmv_pos_slope(st.map, i, 1)]
            thj[st.map.K[k]] = st.s[dwlsmv_pos_slope(st.map, k, 1)]

            r = dwlsmv_cx_pair_fit(Y[., i], Y[., k], X, thi, thj, (i, k))
            R[m, .]    = (i, k, r)
            st.s[pos]  = r[1]

            // v1's entry rule: status 1, or status 6 from the endpoint check.
            if (r[2] == 1) {
                continue
            }
            if (r[2] == 6) {
                if (r[3] != 3) {
                    // Not silent (the user, 2026-11-13): the pair enters, with a note.
                    printf("{txt}note: dwls_sample_stats_cx: the pair y%g~~y%g enters at status 6, residual rho-hat %9.4f.\n",
                           i, k, r[1])
                    printf("{txt}  The sign of the score could not be certified at an end of the range of rho.\n")
                    printf("{txt}  This is a limit of the endpoint sign check, not ordinary convergence (status 1).\n")
                    continue
                }
            }
            cls[m] = 2
        }
    }

    nbad = sum(cls :== 2)
    if (nbad) {
        nbad2 = 0
        errprintf("dwls_sample_stats_cx: %s of %s pairs cannot enter the sample statistics:\n",
                  strofreal(nbad), strofreal(st.map.P))
        for (m = 1; m <= st.map.P; m++) {
            if (cls[m] != 2) {
                continue
            }
            line = "  y" + strofreal(R[m, 1]) + "~~y" + strofreal(R[m, 2]) +
                   ": residual rho-hat " + strofreal(R[m, 3], "%9.4f") +
                   ", status " + strofreal(R[m, 4]) +
                   " (" + dwlsmv_ss_statusword(R[m, 4], R[m, 5]) + ")"
            errprintf("%s\n", line)
            if (R[m, 4] == 2) {
                nbad2 = nbad2 + 1
            }
        }
        if (nbad2) {
            dwlsmv_ss_readings()
        }
        errprintf("  With a covariate no cell repair is applied and no pair is kept at the\n")
        errprintf("  bound. The statuses are defined in spec Module 2.\n")
        exit(3498)
    }

    // ---- 3. Gamma (Hessian form) and the centring diagnostic -----------------
    st.Gamma  = dwls_gamma_cx(Y, X, st.map, st.s)
    st.psibar = dwlsmv_cx_psibar(Y, X, st.map, st.s)
    st.bound  = J(st.map.P, 1, 0)

    return(st)
}

end
