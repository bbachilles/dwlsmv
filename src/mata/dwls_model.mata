*! dwlsmv -- Module 5: the model-implied moment vector sigma(theta)
*!
*! Executed under StataMP 18 on 2026-09-22. sigma(theta) agrees with lavaan's model-implied vector to 5.6e-17 absolute on four configurations: tests/unit/test_model.do check 1.
*!
*! The locked Section 5.1 entry point:
*!
*!     real colvector dwls_model(real colvector theta,
*!                               struct dwlsmv_parmap  scalar par,
*!                               struct dwlsmv_statmap vector maps)
*!
*! maps: one index map per group (spec 4.8 O14), a single group a 1-vector.
*! At G > 1 (built 2026-10-10, contract C2) the result is the stacked
*! (sigma_1', ..., sigma_G')', from dwlsmv_mg_sigma() on each group's rows; the
*! table needs its group key (dwlsmv_mg_check). G = 1 is the path below,
*! unchanged; it never reads the group key and refuses nu, omega, alpha rows.
*!
*! theta is taken EXPLICITLY rather than read from a fitted object, because
*! Module 7 calls this at trial values and Module 3 calls it (through
*! dwls_delta) at STARTING values for the rank check. A signature that
*! presumed a converged fit would support neither caller.
*!
*! Spec Section 4.6's pipeline, and the partitioned vector it ends with:
*!
*!     sigma(theta) = [ tau(theta) ; vech_off(R_y*(theta)) ]
*!
*! indexed exactly as s is. The threshold entries are theta's own values,
*! placed by dwlsmv_pos_tau(); the correlation entries are the corresponding
*! entries of R_y*, placed by dwlsmv_pos_rho() inside dwlsmv_corrvec().
*! Neither block searches the map -- Gate 1 made both positions closed form.
*!
*! WHY THE THRESHOLD ENTRIES ARE A COPY AND NOT A COMPUTATION. Thresholds are
*! direct parameters of sigma. That is the same fact whose derivative form
*! makes Delta's leading block a literal identity (Gate 2, Finding 3), and it
*! is why spec Section 4.8 orders theta with thresholds first: no permutation
*! is needed anywhere between theta and sigma(theta). A permutation appearing
*! between the two is the tell that something has been reordered that should
*! not have been.
*!
*! THE MISSING CHECK AT THE END IS NOT PARANOIA. sigma is allocated missing and
*! every entry must be written. A parameter table carrying fewer tau rows than
*! the index map has thresholds -- which is what a hand-written canonical
*! specification produces if the map was built for a different K -- otherwise
*! leaves a missing entry that travels silently into the DWLS objective and
*! emerges as a non-converging optimizer with no diagnosis attached.
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_pos_rho, dwlsmv_nfac,
*! dwlsmv_parvals, dwlsmv_matrices, dwlsmv_rystar, dwlsmv_corrvec.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwls_model(real colvector theta,
                          struct dwlsmv_parmap  scalar par,
                          struct dwlsmv_statmap vector maps)
{
    struct dwlsmv_statmap scalar map
    struct dwlsmv_parmap  scalar sub
    real colvector vals, sigma, idx
    real scalar i, g

    // O14: one map per group, a single group a 1-vector.
    if (length(maps) < 1) {
        errprintf("dwls_model: no index map passed; a single group is a 1-vector.\n")
        _error("dwls_model: no index map passed")
    }

    // G > 1 (contract C2): sigma_g from group g's rows, stacked group by group
    // in the order of s (spec 4.8 G6; derivation @sec-mg-sigma-delta MS1).
    if (length(maps) > 1) {
        dwlsmv_mg_check(par, length(maps), "dwls_model")
        vals  = dwlsmv_parvals(theta, par)
        sigma = J(0, 1, .)
        for (g = 1; g <= length(maps); g++) {
            sub   = dwlsmv_mg_sub(par, g, idx)
            sigma = sigma \ dwlsmv_mg_sigma(vals[idx], sub, maps[g])
        }
        return(sigma)
    }

    // G = 1: the single-group path, unchanged. It knows no nu, omega or alpha.
    map = maps[1]
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] >= 5) {
            errprintf("dwls_model: row %f (%s) is a nu, omega or alpha row, which exists only with more than one group.\n",
                      i, par.name[i])
            _error("dwls_model: nu, omega and alpha need more than one group")
        }
    }

    vals  = dwlsmv_parvals(theta, par)
    sigma = J(map.D, 1, .)

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 4) {
            sigma[dwlsmv_pos_tau(map, par.row[i], par.col[i])] = vals[i]
        }
    }

    if (map.P > 0) {
        sigma[(map.T + 1)::map.D] = dwlsmv_corrvec(vals, par, map)
    }

    if (missing(sigma)) {
        errprintf("dwls_model: %f of sigma's %f entries were never written.\n",
                  missing(sigma), map.D)
        errprintf("The parameter table carries %f threshold row(s); the index map has %f.\n",
                  sum(par.ptype :== 4), map.T)
        _error("dwls_model: sigma has entries no parameter row filled")
    }

    return(sigma)
}

end
