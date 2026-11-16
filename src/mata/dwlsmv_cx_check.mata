*! dwlsmv -- conditional mode: the covariate's refusals and item separation
*! (spec 4.8 *Conditional covariates*, Coding and Missing x; @sec-condx (c),
*! (g); @sec-condx-policy Result P1; docs/decisions.md 2026-10-25 policy row)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session). Every refusal
*! below is meant to reach the caller as rc 3498 with its message; that is
*! tests/unit/test_condx_stage1.do check 4.
*!
*! What is refused, in this order, each at rc 3498:
*!
*!   1. X not N x 1. This build takes exactly one covariate. Several
*!      covariates are specified only once their derivations land (E-006).
*!   2. A missing x. Stage 1 refuses missing values rather than assuming
*!      deletion happened, as dwlsmv_catcounts() does for Y. The rule is
*!      listwise JOINTLY with the indicators (@sec-condx (f)), so that every
*!      statistic and Gamma refer to one N; it belongs to the caller.
*!   3. x taking one value: no slope is estimable (Result P1's companion).
*!   4. x taking two values other than exactly 0 and 1 (spec 4.8, Coding).
*!      +-1 coding halves gamma and centring moves every threshold, silently.
*!      The tau_U bound also assumes gradient entries 0 and +-1, which needs
*!      0/1 (prereg.md, the phi2 clause).
*!   5. x taking more than two values: continuous x is not derived
*!      (@sec-condx (b): the pattern lemma needs x discrete; E-006).
*!   6. Y's own checks (codes 0..K-1, no missing, every pooled category
*!      observed), by dwlsmv_thresholds(), which aborts with v1's messages.
*!   7. Separation (Result P1, accepted 2026-10-25): item j has a finite MLE
*!      of (tau_j, b_j) iff  max{c : x = 0} > min{c : x = 1}  and
*!      max{c : x = 1} > min{c : x = 0}. Data-only and exact. Every
*!      separated item is named with each group's observed categories, then
*!      Stage 1 aborts.
*!
*! Requires dwlsmv_catcounts, dwlsmv_catdegen, dwlsmv_cumprop and
*! dwlsmv_thresholds.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The covariate alone (refusals 1-5). Separate from the items so that
// dwls_gamma_cx(), which also takes X, can apply the same rule.
void dwlsmv_cx_check_x(real matrix X, real scalar N)
{
    real colvector v
    real scalar nm

    if (cols(X) != 1) {
        errprintf("dwls_sample_stats_cx: X has %s columns; this build takes exactly one covariate.\n",
                  strofreal(cols(X)))
        errprintf("  Several covariates are specified only once their derivations land.\n")
        _error("dwls_sample_stats_cx: X must have exactly one column")
    }
    if (rows(X) != N) {
        errprintf("dwls_sample_stats_cx: X has %s rows but Y has %s.\n",
                  strofreal(rows(X)), strofreal(N))
        _error("dwls_sample_stats_cx: X and Y have different numbers of rows")
    }
    nm = missing(X)
    if (nm) {
        errprintf("dwls_sample_stats_cx: %s of %s values of the covariate are missing.\n",
                  strofreal(nm), strofreal(N))
        errprintf("  Stage 1 refuses missing values rather than assuming deletion happened.\n")
        errprintf("  Delete listwise across the covariate AND all indicators together, so\n")
        errprintf("  that every statistic and Gamma refer to one complete-case sample.\n")
        _error("dwls_sample_stats_cx: missing covariate values -- delete listwise first")
    }

    v = uniqrows(X)
    if (rows(v) == 1) {
        errprintf("dwls_sample_stats_cx: the covariate takes the single value %s.\n",
                  strofreal(v[1]))
        errprintf("  A constant covariate has no estimable effect; drop it from the model.\n")
        _error("dwls_sample_stats_cx: the covariate is constant")
    }
    if (rows(v) > 2) {
        errprintf("dwls_sample_stats_cx: the covariate takes %s distinct values.\n",
                  strofreal(rows(v)))
        errprintf("  This build takes one binary covariate coded 0/1. A continuous covariate\n")
        errprintf("  is not yet derived; an ordinal predictor enters as 0/1 dummies once\n")
        errprintf("  several covariates are supported.\n")
        _error("dwls_sample_stats_cx: the covariate is not binary")
    }
    if (v[1] != 0 | v[2] != 1) {
        errprintf("dwls_sample_stats_cx: the covariate takes the values %s and %s, not 0 and 1.\n",
                  strofreal(v[1]), strofreal(v[2]))
        errprintf("  A binary covariate must be coded exactly 0/1. Other codings change what\n")
        errprintf("  the slopes mean: +-1 coding halves them, and centring moves every\n")
        errprintf("  reported threshold, since thresholds are reported at x = 0. Recode it.\n")
        _error("dwls_sample_stats_cx: a binary covariate must be coded 0/1")
    }
}

void dwlsmv_cx_check(real matrix Y, real colvector K, real matrix X)
{
    real colvector y0, y1, n0, n1, np
    real scalar j, c, nsep
    string scalar s0, s1

    dwlsmv_cx_check_x(X, rows(Y))

    // Pooled categories: v1's checks and messages (empty categories abort).
    // dwlsmv_thresholds() refuses exactly when dwlsmv_catdegen(dwlsmv_catcounts())
    // is not 0; the same two calls name the item first, on that condition only,
    // and the shared function then refuses with v1's own message (rc 3498).
    for (j = 1; j <= cols(Y); j++) {
        np = dwlsmv_catcounts(Y[., j], K[j])
        if (dwlsmv_catdegen(np)) {
            errprintf("dwls_sample_stats_cx: item y%s has an empty category in the pooled sample.\n",
                      strofreal(j))
        }
        (void) dwlsmv_thresholds(Y[., j], K[j])
    }

    // Result P1. Both groups are non-empty (x is 0/1 and not constant).
    nsep = 0
    for (j = 1; j <= cols(Y); j++) {
        y0 = select(Y[., j], X :== 0)
        y1 = select(Y[., j], X :== 1)
        if (max(y0) > min(y1)) {
            if (max(y1) > min(y0)) {
                continue
            }
        }
        if (nsep == 0) {
            errprintf("dwls_sample_stats_cx: separation by the covariate; no finite estimate exists for:\n")
        }
        nsep = nsep + 1
        n0 = dwlsmv_catcounts(y0, K[j])
        n1 = dwlsmv_catcounts(y1, K[j])
        s0 = ""
        s1 = ""
        for (c = 1; c <= K[j]; c++) {
            if (n0[c] > 0) {
                s0 = s0 + " " + strofreal(c - 1)
            }
            if (n1[c] > 0) {
                s1 = s1 + " " + strofreal(c - 1)
            }
        }
        errprintf("  y%s: categories observed at x = 0:%s; at x = 1:%s\n",
                  strofreal(j), s0, s1)
    }
    if (nsep) {
        errprintf("  The two groups' observed categories overlap in at most one category, so\n")
        errprintf("  the likelihood rises without bound along the slope and no estimate of\n")
        errprintf("  that item's thresholds and slope exists.\n")
        _error("dwls_sample_stats_cx: an item is separated by the covariate")
    }
}

end
