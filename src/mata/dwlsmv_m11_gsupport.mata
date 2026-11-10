*! dwlsmv -- Module 11: D9's category recoding per group (invariance() flow)
*!
*!     void dwlsmv_m11_gsupport(touse, indlist, gvar, gvalslist)
*!
*! Build step 3 of the group surface (2026-10-12). D9 under groups, option (b)
*! (decided 2026-10-12, user; spec 4.8 after F4): each indicator is recoded on
*! each group's own observed support, so K is G x p, rows the groups in O13
*! order (gvalslist ascending, group 1 the lowest value), columns the
*! indicators in ind-line order. Nothing is checked here; dwlsmv_m11_gcheck()
*! compares the supports, and dwlsmv_m11_stage1_g() runs Stage 1.
*!
*! Externals, removed by dwlsmv_m11_clear(): dwlsmv_m11_Yg (pointer colvector,
*! one recoded Y per group), dwlsmv_m11_K (G x p) and dwlsmv_m11_sup (G x p
*! string, each support's original values, space-separated, ascending).
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_invariance_flow.do.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_gsupport(string scalar touse, string scalar indlist,
                         string scalar gvar, string scalar gvalslist)
{
    string rowvector ind
    string matrix sup
    real rowvector gv
    real matrix X, Y, K
    real colvector gcol, x, u, r
    real scalar G, g, j, k
    pointer(real matrix) colvector Yg
    string scalar s
    pointer() scalar p

    ind  = tokens(indlist)
    gv   = strtoreal(tokens(gvalslist))
    G    = cols(gv)
    X    = st_data(., ind, touse)
    gcol = st_data(., gvar, touse)
    K    = J(G, cols(ind), 0)
    sup  = J(G, cols(ind), "")
    Yg   = J(G, 1, NULL)
    for (g = 1; g <= G; g++) {
        Y = X[selectindex(gcol :== gv[g]), .]
        for (j = 1; j <= cols(Y); j++) {
            x = Y[., j]
            u = uniqrows(x)
            r = J(rows(x), 1, 0)
            s = strofreal(u[1], "%21.0g")
            for (k = 2; k <= rows(u); k++) {
                r = r + (x :>= u[k])
                s = s + " " + strofreal(u[k], "%21.0g")
            }
            Y[., j]   = r
            K[g, j]   = rows(u)
            sup[g, j] = s
        }
        // a fresh matrix per group: &Y would point every group at one Y
        Yg[g]  = &(J(0, 0, .))
        *Yg[g] = Y
    }
    p = crexternal("dwlsmv_m11_Yg")
    if (p == NULL) {
        p = findexternal("dwlsmv_m11_Yg")
    }
    *p = Yg
    p = crexternal("dwlsmv_m11_K")
    if (p == NULL) {
        p = findexternal("dwlsmv_m11_K")
    }
    *p = K
    p = crexternal("dwlsmv_m11_sup")
    if (p == NULL) {
        p = findexternal("dwlsmv_m11_sup")
    }
    *p = sup
}

end
