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
*! A deficiency at both points aborts with theta0's null directions.
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
    // only if that is deficient too. theta0's sv and V are kept for the
    // message: its null directions are the clean ones. No note when theta2
    // passes.
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
        errprintf("dwls_identify: rank(Delta) is %f at the starting values but %f parameter(s) are free.\n",
                  rank, npar)
        if (rank2 > rank) {
            errprintf("%f parameter(s) are redundant here and %f at a second point. Singular values run %g down to %g.\n",
                      npar - rank, npar - rank2, big, sv[npar])
        }
        else {
            errprintf("%f parameter(s) are redundant. Singular values run %g down to %g.\n",
                      npar - rank, big, sv[npar])
        }
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
            errprintf("At a second, perturbed point the rank is %f, still deficient. That point moves%s", rank2, char(10))
            errprintf("lambda, beta and psi only. These free parameters stayed at\n")
            errprintf("their reference values, and the deficiency may hold only at those values:\n")
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
            errprintf("At a second, perturbed point the rank is %f, still deficient, so this is not\n",
                      rank2)
            errprintf("an artefact of the symmetric starting values.\n")
        }

        // Section 4.8: "a deficiency of k means k redundant parameters and the
        // NULL SPACE says which". Every right singular vector past the rank
        // spans it, and they all come from the decomposition already done, so
        // report all k rather than only the last one.
        for (k = rank + 1; k <= npar; k++) {
            errprintf("Null direction %f of %f loads on:\n", k - rank, npar - rank)
            nullv = V[k, .]'
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
            errprintf("Freeing a psi covariance alongside a beta path between the same two factors does\n")
            errprintf("this: Phi = (I-B)^-1 Psi (I-B)^-T already implies that covariance.\n")
        }
        _error("dwls_identify: Delta is rank deficient at the starting values")
    }
}

end
