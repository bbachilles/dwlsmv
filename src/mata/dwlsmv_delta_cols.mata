*! dwlsmv -- Delta columns for a chosen SET of parameter rows (Module 6)
*!
*! Executed under StataMP 18 on 2026-09-22. tests/unit/test_delta.do check 1 asserts both structural blocks at exactly 0.0, and check 4 puts the computed block against lavaan's analytic Jacobian at ~3e-12.
*!
*! 'sel' indexes rows of the parameter table -- free or FIXED -- and the result
*! has one column per selected row, in the order given. dwls_delta() calls this
*! with the free rows in theta order, which is the v1 use. Spec Section 5.1's
*! design note requires the shape: Sorbom modification indices (Module 12, v2)
*! need Delta columns for parameters that are fixed in the current model, and
*! retrofitting that into a routine built around "the free ones" is the
*! expensive kind of change. It costs one argument now. Nothing in v1 passes a
*! fixed row, and nothing in v1 should.
*!
*! TWO OF THE THREE BLOCKS ARE ASSERTED, NOT COMPUTED, and that is Gate 2's
*! Finding 3 rather than an optimization:
*!
*!   - a THRESHOLD column is a unit vector. d(tau_jc)/d(tau_jc) = 1 and every
*!     other entry is exactly zero, including the entire correlation block --
*!     thresholds do not enter the population polychoric (Finding 2).
*!   - a STRUCTURAL column has exactly zero in every threshold row.
*!
*! Delta is allocated as zeros and those blocks are simply left alone. Spec
*! Section 4.8 orders theta with thresholds first precisely so the leading
*! block comes out as a literal identity of order T, needing no computation and
*! no permutation to line up with s. If you ever find yourself BUILDING that
*! block rather than asserting it, something upstream has reordered theta.
*!
*! THE STEP SIZE. Central differences at h = 1e-5 * max(1, |v|). The balance
*! that matters is between truncation, O(h^2) ~ 1e-10, and cancellation,
*! O(eps/h) ~ 2e-11; the minimizing h for a central difference is around
*! eps^(1/3) = 6e-6 and 1e-5 is the round number next to it. This is not a free
*! choice: Module 3's rank check reads the singular values of this matrix, and
*! a null direction shows up at the noise floor rather than at zero. At 1e-5
*! that floor is ~1e-10 against a threshold of 1e-8 x sigma_max, roughly two
*! orders of margin. The derivation used h = 1e-6, where the floor is ~2e-10
*! and the margin is ten times thinner. Do not "tighten" h.
*!
*! Spec Section 5 prescribes starting numerically and validating any analytic
*! version against this one. tests/unit/test_delta.do compares whole columns
*! against closed forms that hold for these designs -- d(R_jk)/d(lambda_j) =
*! lambda_k on a single-factor model, and the entire beta column on case E --
*! so "numerical" here does not mean "unchecked". Whether analytic derivatives
*! are worth the cost stays Module 6's open call; they are a performance
*! question, not a correctness one, and nothing downstream is waiting on it.
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_pos_rho, dwlsmv_nfac,
*! dwlsmv_matrices, dwlsmv_rystar, dwlsmv_corrvec.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_delta_cols(real colvector vals,
                              struct dwlsmv_parmap  scalar par,
                              struct dwlsmv_statmap scalar map,
                              real colvector sel)
{
    real matrix Dm
    real colvector v, rp, rm
    real scalar k, i, h

    Dm = J(map.D, rows(sel), 0)

    for (k = 1; k <= rows(sel); k++) {
        i = sel[k]

        if (i < 1 | i > par.nall) {
            _error("dwlsmv_delta_cols: sel names a row outside the table")
        }

        if (par.ptype[i] == 4) {
            Dm[dwlsmv_pos_tau(map, par.row[i], par.col[i]), k] = 1
            continue
        }

        if (map.P == 0) {
            continue
        }

        h    = 1e-5 * max((1, abs(vals[i])))
        v    = vals
        v[i] = vals[i] + h
        rp   = dwlsmv_corrvec(v, par, map)
        v[i] = vals[i] - h
        rm   = dwlsmv_corrvec(v, par, map)

        Dm[(map.T + 1)::map.D, k] = (rp - rm) / (2 * h)
    }

    return(Dm)
}

end
