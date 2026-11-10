*! dwlsmv -- Module 3: identification and the parameter map
*!
*! Executed under StataMP 18 on 2026-09-22. tests/unit/test_identify.do: nine table checks and fifteen rejections, every one of them at rc 3498.
*!
*! The locked Section 5.1 entry point:
*!
*!     struct dwlsmv_parmap scalar
*!         dwls_identify(string scalar modelspec,
*!                       struct dwlsmv_statmap vector maps,
*!                       | struct dwlsmv_stage1 vector st)
*!
*! maps holds one index map per group (spec 4.8 O14); a single group is a
*! 1-vector, the only way to pass it. A text with no 'groups' line is a single
*! group and parses exactly as before (G0); st is then not read. A 'groups'
*! text goes to dwlsmv_id_groups() (spec 4.8 G1-G7, built 2026-10-10), which
*! needs st, one Stage 1 result per group, for the threshold starts (O5) and
*! builds one stacked table with one rank check (contract C3).
*!
*! 'modelspec' IS THE CANONICAL MODEL SPECIFICATION OF SPEC SECTION 4.8, NEVER
*! THE USER'S TEXT. Module 11 owns the sem-style surface syntax and normalizes
*! into this form; this routine parses only the canonical form. The seam is as
*! much part of the lock as the type is. It buys three things: Stage 2 is not
*! built against a grammar that is still moving, malformed user input is
*! diagnosed at the ado level rather than inside Mata where _error() truncates
*! destructively at 100 characters, and Modules 3, 5 and 6 are testable before
*! Module 11 exists at all -- a canonical specification is writable by hand.
*!
*! THE GRAMMAR (Section 4.8). One directive per line, '#' comments, blank lines
*! ignored, and order within the file is irrelevant:
*!
*!     factor  <F> [<F> ...]         declaration order fixes B's and Psi's
*!                                   row/column order
*!     ind     <y> [<y> ...]         declaration order is the observed-variable
*!                                   order behind the Section 4.2 index map
*!     lambda  <F> <y>      free [<start>] | fixed <value>
*!     beta    <to> <from>  free [<start>] | fixed <value>
*!     psi     <F> <F>      free [<start>] | fixed <value>
*!     metric  marker <F> <y>  |  unitvariance <F>
*!
*! WHAT IS GENERATED AND WHAT IS NOT. Threshold rows are generated here, one
*! free tau per threshold, from the index map -- which knows K. Everything else
*! must be written: the form is a literal serialization of the parameter table,
*! so an element with no directive is not a parameter and is zero. In
*! particular a factor under a MARKER metric needs an explicit 'psi F F free',
*! and its absence is an error rather than a default, because a default that
*! lives in both Module 3 and Module 11 is a default that will eventually
*! differ between them.
*!
*! 'metric' is the one structural directive, and it is expanded HERE rather
*! than in Module 11. What 'unitvariance F' fixes is psi_ff, the DISTURBANCE
*! variance, whether or not F is endogenous -- for an exogenous factor that is
*! also its total variance, for a predicted one it is not, and the two differ
*! by the explained variance in a way that rescales every loading on that
*! factor. Expanding it upstream would put spec Section 4.8's hardest-won
*! result in the command layer, out of reach of this module's own tests. See
*! tests/unit/test_structural_metric.do for the lavaan measurement behind it
*! and tests/unit/test_identify.do for the assertion that we implement it.
*!
*! THETA'S ORDERING: thresholds first, in Section 4.2's index-map order, then
*! lambda, beta and psi, each ascending by (row, col). The table itself is
*! built in that order with fixed rows interleaved in place, and fidx numbers
*! only the free ones. Thresholds first is what makes Delta's leading block a
*! literal identity of order T, needing no computation and no permutation to
*! align with s -- the tell that this has been violated is a permutation
*! appearing anywhere between theta and sigma(theta).
*!
*! Fixed rows are CARRIED, not dropped. Sorbom modification indices (Module 12,
*! v2) differentiate with respect to currently-fixed parameters, and reshaping
*! a free-only table later is the expensive kind of change.
*!
*! ONE THING THIS DELIBERATELY DOES NOT REJECT, recorded so that nobody adds
*! it as "missing validation":
*!
*!   - An INDICATOR with no lambda directive. Its row of Lambda is zero, so
*!     every model-implied correlation involving it is fixed at 0. That is a
*!     strange model but a well-formed and identified one, and rank(Delta) can
*!     still equal npar, so the rank check does not see it either. The
*!     converse -- a FACTOR with no loadings -- is rejected below, because
*!     nothing observed measures it. lavaan warns here; v1 does not, and the
*!     restraint is deliberate rather than an oversight.
*!
*! A CYCLIC B IS REJECTED, since 2026-09-23 (spec Section 4.8): v1 fits
*! recursive models only. Until then this header listed it as deliberately
*! allowed, which was wrong in two ways. With (I - B) nonsingular a cycle is
*! computable, but at the default starts (B = 0, Psi = I) dPhi = dB + dB', so a
*! reciprocal pair's two Delta columns are identical and the rank check called
*! EVERY such model overparameterized, identified or not. And a feedback loop
*! needs a stability condition (spectral radius of B below 1) that nothing
*! checks. So the model is refused up front with the real reason, via a
*! transitive closure of B's directed graph below.
*!
*! TAU ROWS CARRY THE INDEX MAP'S NAMES, not the user's. par.name for a
*! threshold is map.label -- "y3|t2" -- whatever the 'ind' directive called
*! variable 3, because the map is keyed on position and K alone. The rank
*! -check diagnostic below can therefore print "y3|t2" beside
*! "lambda[gender,F1]" in the same list. Harmless in v1, where Module 11 owns
*! surface naming and nothing user-facing reads these names; noted so that it
*! is recognised as a naming seam rather than "fixed" by having Module 3
*! rewrite the map's labels, which would put two spellings of one statistic
*! into circulation.
*!
*! THE RANK CHECK, and the one subtlety in it. Section 4.8 requires
*! rank(Delta) to equal the number of free parameters, checked at the STARTING
*! values before any optimization, and names the case it exists for: freeing
*! psi_12 alongside a beta path between the same two factors, when
*! Phi = (I-B)^-1 Psi (I-B)^-T already implies that covariance.
*!
*! The derivation's Finding 4 (rewritten 2026-09-23 under the settled
*! disturbance-variance metric) derives the mechanism, and it decides what the
*! evaluation point changes:
*!
*!   - with psi_22 fixed at 1, Phi_12 = beta + psi_12 and
*!     Phi_22 = beta^2 + 2 beta psi_12 + 1. Freeing psi_12 makes Phi_22 vary,
*!     which re-opens the factor-scale indeterminacy on F2 that the metric
*!     constraint existed to remove, so the deficiency is still exactly one --
*!     but its direction now also rescales F2's LOADINGS, with beta and psi_12
*!     absorbing the rescaling.
*!   - AT psi_12 = 0 the two collapse. d Phi_22 vanishes along
*!     (d beta = -d psi_12) and the null direction becomes exactly
*!     (beta, psi_12). The two Delta columns then agree analytically, and
*!     numerically to a few ulps -- the measured ratio is ~2e-17 rather than
*!     identically 0, so the central differences are not *bitwise* identical,
*!     only equal to rounding.
*!
*! psi_12 = 0 is where this runs, because 0 is a disturbance covariance's
*! default start. Measured on case E, the singular-value ratio is 6.0e-17 at
*! the defaults and 2.1e-17 at the oracle's beta = 0.4, both with psi_12 = 0.
*! AWAY FROM psi_12 = 0 THE CHECK STILL FIRES: the deficiency is still exactly
*! one, at ratios of 4.8e-13 and 8.9e-13 at psi_12 = 0.1 and 0.3, but the null
*! direction rescales F2's loadings instead of naming (beta, psi_12) alone
*! (test_delta.do check 8; derivation, Gate 2 Finding 4). What checking at the
*! starting values buys is that clean null direction, which is the one the
*! error message reports. Until 2026-09-23 this header said a nonzero psi_12
*! gives full rank; that was never measured and is false.
*!
*! A SECOND POINT (decided 2026-10-04; spec 4.8, derivation Finding 5). The
*! starts are symmetric, and a model of full rank almost everywhere can lose
*! rank there (a random intercept with psi_21 free, at equal loadings; two
*! factors of two indicators, at psi_21 = 0). When the starts are deficient
*! the check is repeated at theta2, every free lambda/beta/psi entry shifted
*! by 0.1*(frac(i*phi) - 1/2), and only a deficiency at both aborts. The
*! message reports theta0's null directions. When theta2 is what passes, it
*! becomes the optimizer's start (written into par.value): at theta0
*! Delta'W Delta is singular, and a fit started there did not converge.
*!
*! The 1e-8 threshold is therefore NOT sized against that gap, which is
*! machine zero. It is sized against the numerical Delta's own noise floor,
*! ~1e-10 relative (see dwlsmv_delta_cols on the step size), while the
*! smallest genuinely nonzero singular value measured across the four
*! identified configurations of the tier-2 oracle is 0.226. Two orders of
*! margin below, seven above. test_identify.do and test_delta.do both print
*! the measured values. Since 2026-09-26 the criterion is dwlsmv_numrank(),
*! the one definition Module 10's bread diagnostic also uses; this file no
*! longer carries the constant.
*!
*! Load order: dwlsmv_types, dwlsmv_statmap_build, dwlsmv_pos_tau,
*! dwlsmv_pos_rho, dwlsmv_nfac, dwlsmv_parvals, dwlsmv_matrices,
*! dwlsmv_rystar, dwlsmv_corrvec, dwlsmv_delta_cols, dwlsmv_mg_* (five),
*! dwls_delta, dwlsmv_canon_lines, dwlsmv_canon_find, dwlsmv_canon_status,
*! dwlsmv_numrank, dwlsmv_id_parse, dwlsmv_id_rank, dwlsmv_id_tie,
*! dwlsmv_id_groups. (Mata links every referenced function on first call.)
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

struct dwlsmv_parmap scalar dwls_identify(string scalar modelspec,
                                          struct dwlsmv_statmap vector maps,
                                          | struct dwlsmv_stage1 vector st)
{
    struct dwlsmv_statmap scalar map
    struct dwlsmv_stage1 vector st0
    struct dwlsmv_parmap scalar par
    string colvector lines, fac, ind
    string rowvector w
    real colvector eptype, erow, ecol, efree, eval
    real scalar nl, i, r, c
    real scalar nfac, nel, nall, npar, pos, hit

    lines = dwlsmv_canon_lines(modelspec)
    if (rows(lines) == 0) {
        errprintf("dwls_identify: the canonical model specification has no directives.\n")
        _error("dwls_identify: the model specification is empty")
    }

    // A 'groups' text (spec 4.8 G1-G7) is the group path: one stacked table,
    // one rank check (contract C3). Its threshold starts come from Stage 1.
    for (nl = 1; nl <= rows(lines); nl++) {
        w = tokens(lines[nl])
        if (w[1] == "groups") {
            if (args() >= 3) {
                return(dwlsmv_id_groups(lines, maps, st, 1))
            }
            st0 = dwlsmv_stage1(0)
            return(dwlsmv_id_groups(lines, maps, st0, 0))
        }
    }

    // No 'groups' line: a single group, parsed exactly as before (G0). One map;
    // st, if passed, is not read.
    if (length(maps) != 1) {
        errprintf("dwls_identify: %f index maps passed for a text with no 'groups' line; a single\n",
                  length(maps))
        errprintf("group is a 1-vector (spec 4.8 O14, G7).\n")
        _error("dwls_identify: index maps do not match the number of groups")
    }
    map = maps[1]

    // Passes 1-3 and the section checks: dwlsmv_id_parse(), shared with each
    // group section, in single-group mode.
    dwlsmv_id_parse(lines, map, 1, 0, fac, ind, eptype, erow, ecol, efree, eval)
    nfac = rows(fac)

    // =====================================================================
    // Build the table, in theta order: tau, then lambda, beta, psi.
    // =====================================================================
    nel  = rows(eptype)
    nall = map.T + nel

    par.nall  = nall
    par.name  = J(nall, 1, "")
    par.ptype = J(nall, 1, 0)
    par.row   = J(nall, 1, 0)
    par.col   = J(nall, 1, 0)
    par.free  = J(nall, 1, 0)
    par.value = J(nall, 1, 0)
    par.fidx  = J(nall, 1, 0)
    par.grp   = J(nall, 1, 1)   // spec 4.8 G6 group key; single group

    pos = 0

    // Thresholds. Positions 1..T of the index map ARE the threshold block, in
    // Section 4.2's order, so this loop needs no sort and no lookup. All free:
    // v1 has no constrained thresholds.
    for (i = 1; i <= map.T; i++) {
        pos            = pos + 1
        par.name[pos]  = map.label[i]
        par.ptype[pos] = 4
        par.row[pos]   = map.var_i[i]
        par.col[pos]   = map.thr_k[i]
        par.free[pos]  = 1
        par.value[pos] = 0
    }

    // Lambda, ascending by (row, col) = (indicator, factor).
    for (r = 1; r <= map.p; r++) {
        for (c = 1; c <= nfac; c++) {
            hit = 0
            for (i = 1; i <= nel; i++) {
                if (eptype[i] == 1 & erow[i] == r & ecol[i] == c) {
                    hit = i
                }
            }
            if (hit > 0) {
                pos            = pos + 1
                par.name[pos]  = "lambda[" + ind[r] + "," + fac[c] + "]"
                par.ptype[pos] = 1
                par.row[pos]   = r
                par.col[pos]   = c
                par.free[pos]  = efree[hit]
                par.value[pos] = eval[hit]
            }
        }
    }

    // Beta, ascending by (row, col) = (to, from).
    for (r = 1; r <= nfac; r++) {
        for (c = 1; c <= nfac; c++) {
            hit = 0
            for (i = 1; i <= nel; i++) {
                if (eptype[i] == 3 & erow[i] == r & ecol[i] == c) {
                    hit = i
                }
            }
            if (hit > 0) {
                pos            = pos + 1
                par.name[pos]  = "beta[" + fac[r] + "," + fac[c] + "]"
                par.ptype[pos] = 3
                par.row[pos]   = r
                par.col[pos]   = c
                par.free[pos]  = efree[hit]
                par.value[pos] = eval[hit]
            }
        }
    }

    // Psi, lower triangle, ascending by (row, col).
    for (r = 1; r <= nfac; r++) {
        for (c = 1; c <= r; c++) {
            hit = 0
            for (i = 1; i <= nel; i++) {
                if (eptype[i] == 2 & erow[i] == r & ecol[i] == c) {
                    hit = i
                }
            }
            if (hit > 0) {
                pos            = pos + 1
                par.name[pos]  = "psi[" + fac[r] + "," + fac[c] + "]"
                par.ptype[pos] = 2
                par.row[pos]   = r
                par.col[pos]   = c
                par.free[pos]  = efree[hit]
                par.value[pos] = eval[hit]
            }
        }
    }

    if (pos != nall) {
        _error("dwls_identify: the parameter table lost a row while being ordered")
    }

    npar = 0
    for (i = 1; i <= nall; i++) {
        if (par.free[i]) {
            npar         = npar + 1
            par.fidx[i]  = npar
        }
    }
    par.npar = npar

    // Section 4.8's admissibility check: rank(Delta) at the starting values,
    // with the second point. Shared with the group path (dwlsmv_id_rank).
    dwlsmv_id_rank(par, maps)
    return(par)
}

end
