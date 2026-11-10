*! dwlsmv -- Module 6: Delta = d sigma(theta) / d theta'
*!
*! Executed under StataMP 18 on 2026-09-22. tests/unit/test_delta.do.
*!
*! The locked Section 5.1 entry point:
*!
*!     real matrix dwls_delta(real colvector theta,
*!                            struct dwlsmv_parmap  scalar par,
*!                            struct dwlsmv_statmap vector maps)
*!
*! maps: one index map per group (spec 4.8 O14), a single group a 1-vector.
*! At G > 1 (built 2026-10-10, contract C2) the result is the stacked
*! (sum_g D_g) x npar Delta from dwlsmv_mg_dcols(), tied rows summed into one
*! column. G = 1 is the path below, unchanged; it never reads the group key
*! and refuses nu, omega, alpha rows.
*!
*! D x npar, with column k belonging to theta[k]. Like dwls_model() it takes
*! theta explicitly, because its two callers evaluate it at different points:
*! Module 7 at trial values, and Module 3 at STARTING values for the Section
*! 4.8 rank check, before any optimization has happened.
*!
*! This routine is a thin wrapper. It resolves theta and the table's fixed
*! values into one vector, works out which table rows the free parameters are,
*! and hands both to dwlsmv_delta_cols(), which owns the differentiation and
*! the two structural blocks. The split exists so that the parameter-SET shape
*! Module 12 needs is available without this signature -- which is locked --
*! having to grow an argument.
*!
*! 'sel' is built by inverting fidx rather than by filtering in table order.
*! Those two orders are the same today, because spec Section 4.8 orders theta
*! exactly as the table is ordered, but column k must belong to theta[k] even
*! if that ever stops being true. Inverting fidx says so; filtering assumes it.
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_pos_rho, dwlsmv_nfac,
*! dwlsmv_parvals, dwlsmv_matrices, dwlsmv_rystar, dwlsmv_corrvec,
*! dwlsmv_delta_cols.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwls_delta(real colvector theta,
                       struct dwlsmv_parmap  scalar par,
                       struct dwlsmv_statmap vector maps)
{
    struct dwlsmv_statmap scalar map
    struct dwlsmv_parmap  scalar sub
    real matrix Dm, Dg
    real colvector vals, sel, idx, hit
    real scalar i, g, k, off, Dtot

    // O14: one map per group, a single group a 1-vector.
    if (length(maps) < 1) {
        errprintf("dwls_delta: no index map passed; a single group is a 1-vector.\n")
        _error("dwls_delta: no index map passed")
    }

    // G > 1 (contract C2): rows of Delta in the order of the stacked s, group
    // by group; column t sums the columns of every row with fidx t, so a tie
    // is one shared column (derivation @sec-mg-sigma-delta MS3).
    if (length(maps) > 1) {
        dwlsmv_mg_check(par, length(maps), "dwls_delta")
        vals = dwlsmv_parvals(theta, par)
        Dtot = 0
        for (g = 1; g <= length(maps); g++) {
            Dtot = Dtot + maps[g].D
        }
        Dm  = J(Dtot, par.npar, 0)
        hit = J(par.npar, 1, 0)
        off = 0
        for (g = 1; g <= length(maps); g++) {
            sub = dwlsmv_mg_sub(par, g, idx)
            Dg  = dwlsmv_mg_dcols(vals[idx], sub, maps[g])
            for (k = 1; k <= sub.nall; k++) {
                if (sub.free[k]) {
                    i = sub.fidx[k]
                    Dm[(off + 1)::(off + maps[g].D), i] =
                        Dm[(off + 1)::(off + maps[g].D), i] + Dg[., k]
                    hit[i] = 1
                }
            }
            off = off + maps[g].D
        }
        if (min(hit) < 1) {
            _error("dwls_delta: fidx does not number the free rows 1..npar")
        }
        return(Dm)
    }

    // G = 1: the single-group path, unchanged. It knows no nu, omega or alpha.
    map = maps[1]
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] >= 5) {
            errprintf("dwls_delta: row %f (%s) is a nu, omega or alpha row, which exists only with more than one group.\n",
                      i, par.name[i])
            _error("dwls_delta: nu, omega and alpha need more than one group")
        }
    }

    vals = dwlsmv_parvals(theta, par)

    sel = J(par.npar, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            sel[par.fidx[i]] = i
        }
    }

    if (min(sel) < 1) {
        _error("dwls_delta: fidx does not number the free rows 1..npar")
    }

    return(dwlsmv_delta_cols(vals, par, map, sel))
}

end
