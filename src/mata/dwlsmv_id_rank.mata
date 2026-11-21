*! dwlsmv -- Module 3: the rank check at the starting values (spec 4.8)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_identify.do (unedited), tests/unit/test_identify_groups.do.
*!
*! Moved out of dwls_identify() on 2026-10-10 so that the single-group table
*! and the stacked multiple-group table go through one check (spec 4.8 G6:
*! "the Delta rank check runs once, on the stacked table"; contract C3: never
*! separately by group). Delta is dwls_delta(theta, par, maps), which stacks
*! the groups when more than one map is passed.
*!
*! The second point (spec 4.8, A second point before aborting; O7 under
*! groups): every free lambda, beta and psi entry (ptypes 1-3) is shifted by
*! 0.1*(frac(i*phi) - 1/2), i its theta position, a tie being one theta
*! element and so one shift; thresholds, nu, omega and alpha are not moved.
*! Until 2026-10-10 the condition read "not a threshold", which is the same
*! set for a single group (no ptypes 5-7 exist there) and is wrong under
*! groups, where it would also move nu, omega and alpha.
*!
*! Passing at theta2 writes theta2 into par.value (the optimizer's start).
*! A deficiency at both points aborts with theta2's null directions, labelled a
*! diagnostic at that point and not a proof of nonidentification (B3.13, owner,
*! 2026-11-20; before then theta0's were printed).
*!
*! Load order: dwlsmv_types, dwls_delta (and what it needs), dwlsmv_numrank.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_id_rank(struct dwlsmv_parmap scalar par,
                    struct dwlsmv_statmap vector maps)
{
    real colvector theta0, theta2, sv, sv2, nullv
    real matrix Dm, V, V2
    real scalar npar, nall, Dtot, rank, rank2, big, phi, i, j, k, g, hint

    npar = par.npar
    nall = par.nall
    Dtot = 0
    for (g = 1; g <= length(maps); g++) {
        Dtot = Dtot + maps[g].D
    }

    // A backstop, and not reachable from any canonical specification: every
    // threshold becomes a FREE tau row and dwlsmv_statmap_build() guarantees
    // T >= p >= 1, so npar >= 1 however much the user fixes. Kept because
    // npar == 0 would otherwise reach svd() with a D x 0 matrix.
    if (npar == 0) {
        errprintf("dwls_identify: the specification has no free parameters.\n")
        _error("dwls_identify: no free parameters")
    }
    if (npar > Dtot) {
        errprintf("dwls_identify: %f free parameter(s) against %f sample statistic(s), so df = %f.\n",
                  npar, Dtot, Dtot - npar)
        _error("dwls_identify: more free parameters than sample statistics")
    }

    theta0 = J(npar, 1, 0)
    for (i = 1; i <= nall; i++) {
        if (par.free[i]) {
            theta0[par.fidx[i]] = par.value[i]
        }
    }

    V    = J(0, 0, 0)
    sv   = J(0, 1, 0)
    Dm   = dwls_delta(theta0, par, maps)
    rank = dwlsmv_numrank(Dm, sv, V)

    // A backstop, not a live branch. Delta's threshold rows carry an identity
    // (a selection under groups) and dwlsmv_statmap_build() refuses K_j < 2, so
    // sv[1] >= 1 on every reachable input. It is kept because the cost of a
    // division by zero here is a matrix of missing values that looks like a
    // rank-0 model -- but do not read it as a path that has been exercised.
    big = sv[1]
    if (big <= 0) {
        errprintf("dwls_identify: Delta is identically zero at the starting values.\n")
        _error("dwls_identify: Delta is identically zero")
    }

    // The criterion, 1e-8 of the largest singular value, is dwlsmv_numrank()'s
    // and only its; the sizing argument lives in that file's header.
    //
    // Section 4.8's second point (derivation Finding 5): a deficiency at the
    // symmetric starts alone is not a redundancy. Re-check at theta2 and abort
    // only if that is deficient too; the message then reports theta2's sv2
    // and V2 (B3.13). No note when theta2 passes.
    if (rank < npar) {
        phi    = (sqrt(5) - 1) / 2
        theta2 = theta0
        for (i = 1; i <= nall; i++) {
            if (par.free[i] & par.ptype[i] <= 3) {
                k         = par.fidx[i]
                theta2[k] = theta0[k] + 0.1 * (mod(k * phi, 1) - 0.5)
            }
        }
        V2    = J(0, 0, 0)
        sv2   = J(0, 1, 0)
        rank2 = dwlsmv_numrank(dwls_delta(theta2, par, maps), sv2, V2)
        // Passing at theta2 makes theta2 the optimizer's start: Delta'W Delta is
        // singular at theta0, and a fit started there did not converge
        // (tools/probe_rank_starts.do, 2026-10-04). par.value carries starts.
        if (rank2 == npar) {
            for (i = 1; i <= nall; i++) {
                if (par.free[i]) {
                    par.value[i] = theta2[par.fidx[i]]
                }
            }
            return
        }
    }
    if (rank < npar) {
        // Wording for an applied reader (owner, 2026-11-20): plain language
        // first, the linear algebra in one closing "Technical detail" line.
        errprintf("dwls_identify: the model appears not to be identified. It has %f free\n", npar)
        errprintf("parameters, but the data seem to determine only %f of them.\n", rank2)
        // O7 leaves nu, omega and alpha out of theta2, so a free one at its
        // reference value was never moved, and the deficiency may be tied to it.
        hint = 0
        for (i = 1; i <= nall; i++) {
            if (par.free[i] & par.ptype[i] >= 5) {
                if (par.value[i] == (par.ptype[i] == 6 ? 1 : 0)) {
                    hint = 1
                }
            }
        }
        if (hint) {
            errprintf("The check was repeated at changed starting values, but only loadings, paths\n")
            errprintf("and (co)variances were changed. These free parameters stayed at their\n")
            errprintf("reference values, and the problem may hold only at those values:\n")
            for (i = 1; i <= nall; i++) {
                if (par.free[i] & par.ptype[i] >= 5) {
                    if (par.value[i] == (par.ptype[i] == 6 ? 1 : 0)) {
                        errprintf("    %s\n", par.name[i])
                    }
                }
            }
            errprintf("A different start is checked at that value instead. That can pass this check,\n")
            errprintf("but it does not identify the model if the true value is the reference value.\n")
        }
        else {
            errprintf("The problem remained when the starting values were changed, so it is unlikely\n")
            errprintf("to come from the starting values alone. This is a numerical check at two sets\n")
            errprintf("of trial values, not a proof that the model is unidentified.\n")
        }

        // Section 4.8: "a deficiency of k means k redundant parameters and the
        // NULL SPACE says which". Every right singular vector past the rank
        // spans it, and they all come from the decomposition already done, so
        // report all of them rather than only the last one. B3.13 (owner,
        // 2026-11-20): the directions are theta2's, and are labelled a
        // diagnostic at that point. theta0's mislead when B = 0 adds
        // deficiencies of its own (G over two, @sec-second-order correction 6).
        errprintf("Where to look: the parameters below appear to trade off against one another.\n")
        errprintf("At the trial values, changing them together in approximately the proportions\n")
        errprintf("shown has almost no effect on the model-implied quantities. This often\n")
        errprintf("indicates that an identifying constraint is missing or that these parameters\n")
        errprintf("should not all be free.\n")
        for (k = rank2 + 1; k <= npar; k++) {
            errprintf("Set %f of %f:\n", k - rank2, npar - rank2)
            nullv = V2[k, .]'
            for (i = 1; i <= npar; i++) {
                if (abs(nullv[i]) > 0.1) {
                    for (j = 1; j <= nall; j++) {
                        if (par.free[j] & par.fidx[j] == i) {
                            errprintf("    %-22s %9.4f\n", par.name[j], nullv[i])
                        }
                    }
                }
            }
        }

        // The redundancy this names exists only when the model frees both
        // elements, between one pair of factors, in one group.
        hint = 0
        for (i = 1; i <= nall; i++) {
            if (par.free[i] & par.ptype[i] == 2 & par.row[i] != par.col[i]) {
                for (j = 1; j <= nall; j++) {
                    if (par.free[j] & par.ptype[j] == 3 & par.grp[j] == par.grp[i] &
                        min((par.row[j], par.col[j])) == min((par.row[i], par.col[i])) &
                        max((par.row[j], par.col[j])) == max((par.row[i], par.col[i]))) {
                        hint = 1
                    }
                }
            }
        }
        if (hint) {
            errprintf("This model frees both a path and a disturbance covariance between the same\n")
            errprintf("two factors. The path already accounts for that covariance, so the two cannot\n")
            errprintf("both be estimated. Loadings listed in the same set belong to that problem.\n")
        }
        errprintf("Technical detail: rank(Delta) is %f at the starting values and %f at a second\n",
                  rank, rank2)
        errprintf("point, with %f free parameters. The smallest singular value there is %s\n",
                  npar, strtrim(sprintf("%9.3g", sv2[npar])))
        errprintf("(largest: %s). Singular values at or below 1e-8 times the largest are\n",
                  strtrim(sprintf("%9.3g", sv2[1])))
        errprintf("treated as zero.\n")
        _error("dwls_identify: the model appears not to be identified (Delta is rank deficient)")
    }
}

end
