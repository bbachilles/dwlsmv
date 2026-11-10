*! dwlsmv -- Module 11: translate sem's multiple-group reading into the configural text
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_translate_groups.do.
*!
*!   dwlsmv_translate_g(final, gname, gvals)
*!
*! Spec principle 4 (group surface, decided 2026-10-10), the 2026-10-12 decisions
*! "Group-specific paths under invariance()" (a), "invariance() requires an
*! explicit metric(unitvariance)" (b) and spec 4.8 G5, G8. _dwlsmv_translate
*! runs sem <paths>, group(gname) ginvariant(none) noestimate on dummy data whose
*! group variable takes the user's values gvals (sorted, O13); this function
*! reads e() into the multiple-group configural canonical text that
*! dwlsmv_levelgen takes (L1 (b)): "groups G", then per group "group g",
*! "factor", "ind" and that group's lines, group 1 scaled by
*! "metric unitvariance F", later groups by "psi F F fixed 1".
*!
*! Each group's columns are found from the group tag in sem's stripe
*! (y1:2.grp#c.F, y1:2.grp, /var(F)#2.grp; tools/probe_sem_group_surface.do,
*! cases 5.1-5.7), the tag removed, and the single-group translation
*! dwlsmv_tr_core run on that group's columns and Cns rows, so every D1-D3
*! refusal holds per group. sem's own factor means (group 1 fixed at 0, later
*! groups free, ginvariant(none)) have no counterpart in the configural text
*! (means enter only as alpha, at the scalar level, from the generator) and are
*! skipped; a mean the user constrained otherwise is refused.
*!
*! final = 0 returns the text unchecked (D2's first pass, for the injection).
*! final = 1 refuses at rc 198, in this order (most specific first, as G7):
*!   a loading fixed at 1 (point 8's message, per factor);
*!   a structural path between factors (G5);
*!   any difference between a later group and group 1 (L1, Module 11 first;
*!   dwlsmv_levelgen's rc 3498 is the backstop), every difference listed.
*! A constraint tying parameters of different groups is refused on either pass.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string colvector dwlsmv_translate_g(real scalar final, string scalar gname,
    real rowvector gvals)
{
    string matrix st, sg
    real matrix C, Cg
    string colvector fac, ind, oxv, out, L
    real colvector gi, skip, keep, J, rr
    real scalar G, j, r, k, g, gv, nzg, sumg
    string scalar co, rx1, rx2, f

    st  = st_matrixcolstripe("e(b)")
    C   = st_matrix("e(Cns)")
    fac = tokens(st_global("e(lxvars)") + " " + st_global("e(lyvars)"))'
    ind = tokens(st_global("e(oyvars)"))'
    oxv = tokens(st_global("e(oxvars)"))'
    G   = cols(gvals)

    // group of every column, tag removed
    rx1 = "^([0-9]+)(bn)?\." + gname + "(#c\.(.+))?$"
    rx2 = "^(.+)#([0-9]+)(bn)?\." + gname + "$"
    sg   = st
    gi   = J(rows(st), 1, 0)
    skip = J(rows(st), 1, 0)
    for (j = 1; j <= rows(st); j++) {
        co = st[j, 2]
        if (st[j, 1] != "/" & ustrregexm(co, rx1)) {
            gv = strtoreal(ustrregexs(1))
            sg[j, 2] = (ustrregexs(4) == "" ? "_cons" : ustrregexs(4))
        }
        else if (st[j, 1] == "/" & ustrregexm(co, rx2)) {
            gv = strtoreal(ustrregexs(2))
            sg[j, 2] = ustrregexs(1)
        }
        else {
            dwlsmv_tr_refuse("dwlsmv: sem reported a parameter dwlsmv cannot read (" + co + ").")
        }
        for (g = 1; g <= G; g++) {
            if (gvals[g] == gv) {
                gi[j] = g
            }
        }
        if (gi[j] == 0) {
            dwlsmv_tr_refuse("dwlsmv: sem reported a group value dwlsmv did not expect (" + co + ").")
        }
        if (substr(sg[j, 2], 1, 5) == "mean(") {
            skip[j] = 1
        }
    }

    // Cns rows: all of a row's coefficients in one group; sem's default
    // means (group 1 fixed at 0, later groups free) and nothing else on them
    if (rows(C) == 0) {
        C = J(0, rows(st) + 1, 0)
    }
    keep = J(rows(C), 1, 1)
    for (r = 1; r <= rows(C); r++) {
        nzg = 0
        sumg = 0
        for (j = 1; j <= rows(st); j++) {
            if (C[r, j] != 0) {
                if (nzg == 0) {
                    nzg = gi[j]
                }
                else if (gi[j] != nzg) {
                    errprintf("dwlsmv: a constraint ties parameters of different groups (for example a\n")
                    errprintf("shared @name). Equality constraints across groups are not available in\n")
                    errprintf("this version; the level named in invariance() sets which parameters are\n")
                    errprintf("equal across groups. No model was fitted.\n")
                    exit(198)
                }
                if (skip[j]) {
                    sumg++
                    if (gi[j] != 1 | C[r, cols(C)] != 0) {
                        dwlsmv_tr_mean(sg[j, 2])
                    }
                    keep[r] = 0
                }
            }
        }
        if (sumg > 0 & sumg != 1) {
            dwlsmv_tr_mean("mean")
        }
    }
    for (j = 1; j <= rows(st); j++) {
        if (skip[j] & gi[j] == 1) {
            // group 1's mean must be fixed (at 0, checked above)
            k = 0
            for (r = 1; r <= rows(C); r++) {
                if (C[r, j] != 0) {
                    k = 1
                }
            }
            if (!k) {
                dwlsmv_tr_mean(sg[j, 2])
            }
        }
    }

    out = "groups " + strofreal(G)
    for (g = 1; g <= G; g++) {
        J = selectindex((gi :== g) :& !skip)
        rr = J(0, 1, .)
        for (r = 1; r <= rows(C); r++) {
            if (keep[r]) {
                if (any(C[r, J'] :!= 0)) {
                    rr = rr \ r
                }
            }
        }
        if (rows(rr)) {
            Cg = C[rr, (J', cols(C))]
        }
        else {
            Cg = J(0, rows(J) + 1, 0)
        }
        L = dwlsmv_tr_core(sg[J, .], Cg, fac, ind, oxv, final)
        if (g > 1) {
            for (k = 3; k <= rows(L); k++) {
                if (substr(L[k], 1, 20) == "metric unitvariance ") {
                    f = substr(L[k], 21, .)
                    L[k] = "psi " + f + " " + f + " fixed 1"
                }
            }
            L = L[1 :: 2] \ sort(L[3 :: rows(L)], 1)
        }
        out = out \ ("group " + strofreal(g)) \ L
    }
    if (final) {
        dwlsmv_tr_gcheck(out, gname, gvals)
    }
    return(out)
}

void dwlsmv_tr_mean(string scalar m)
{
    errprintf("dwlsmv: a factor mean is constrained (%s). With invariance(), factor means\n", m)
    errprintf("are set by the invariance level: 0 in every group, except at invariance(scalar),\n")
    errprintf("where they are estimated in groups after the first. No model was fitted.\n")
    exit(198)
}

end
