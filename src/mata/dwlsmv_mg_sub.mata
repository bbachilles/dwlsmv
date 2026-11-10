*! dwlsmv -- one group's rows of a multiple-group parameter table (Modules 5, 6)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_mg_model.do.
*!
*! Returns the rows of par whose group key is g, in table order, as a table of
*! their own, and the positions of those rows in par (idx). Group g's sigma
*! reads only group g's rows (spec 4.8 G6; derivation @sec-mg-sigma-delta MS1),
*! so the single-group helpers (dwlsmv_matrices, dwlsmv_corrvec,
*! dwlsmv_delta_cols) can run on it unchanged: they ignore ptypes 5-7.
*!
*! npar and fidx are copied as they are, so a tied row keeps its pooled theta
*! position; the sub-table is never handed to dwlsmv_parvals.
*!
*! src/mata/dwlsmv_types.mata must be loaded first.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

struct dwlsmv_parmap scalar dwlsmv_mg_sub(struct dwlsmv_parmap scalar par,
                                          real scalar g,
                                          real colvector idx)
{
    struct dwlsmv_parmap scalar s
    real scalar i, n

    n = 0
    for (i = 1; i <= par.nall; i++) {
        n = n + (par.grp[i] == g)
    }
    idx = J(n, 1, 0)
    n = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.grp[i] == g) {
            n++
            idx[n] = i
        }
    }

    s.npar  = par.npar
    s.nall  = n
    s.name  = par.name[idx]
    s.ptype = par.ptype[idx]
    s.row   = par.row[idx]
    s.col   = par.col[idx]
    s.free  = par.free[idx]
    s.value = par.value[idx]
    s.fidx  = par.fidx[idx]
    s.grp   = par.grp[idx]
    return(s)
}

end
