*! dwlsmv -- Module 11: the invariance() table checks before any fit
*!
*!     void dwlsmv_m11_gcheck(canon, level, indlist, gvar, gvalslist)
*!
*! Build step 3 of the group surface (2026-10-12). Past configural only, in
*! F4's order (spec 4.8, decided 2026-10-12, user):
*!
*! 1. Equal categories (G3, "No silent partial tie"; D9 under groups (b)).
*!    Every level past configural ties every indicator's thresholds, so every
*!    indicator's observed original values must be the same set in every
*!    group. One line per offending indicator, every group named by position
*!    and value: the count line when the counts differ, else the values line;
*!    then the shared closing sentence once. Printed here; the local gref is
*!    set to 1 and dwlsmv.ado exits 198 (the wording rows: display as error and
*!    exit 198, no _error()).
*! 2. The binary scalar refusal (T0, T2, F4), only when 1 passed: at
*!    invariance(scalar), a factor is all-binary iff every indicator of it has
*!    K = 2 on the now common K. A factor's indicators are its group 1 lambda
*!    lines free or fixed nonzero (dwlsmv_levelgen's rule). Message G-A
*!    (E-079) by errprintf, then
*!    E1, exit(198): a refusal, so no Mata traceback.
*!
*! Reads dwlsmv_m11_K and dwlsmv_m11_sup (dwlsmv_m11_gsupport()). Nothing is
*! fitted; at configural nothing is checked (no tie).
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_invariance_flow.do.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_gcheck(string scalar canon, string scalar level, string scalar indlist,
                       string scalar gvar, string scalar gvalslist)
{
    string rowvector ind, gv, w, fac
    string colvector L
    string matrix sup
    real matrix K
    real scalar G, g, j, i, f, bad, nb, allb, has, on
    string scalar line, bin

    st_local("gref", "0")
    if (level == "configural") {
        return
    }
    ind = tokens(indlist)
    gv  = tokens(gvalslist)
    K   = *findexternal("dwlsmv_m11_K")
    sup = *findexternal("dwlsmv_m11_sup")
    G   = rows(K)

    // ---- 1. the same categories in every group --------------------------
    bad = 0
    for (j = 1; j <= cols(ind); j++) {
        if (any(sup[., j] :!= sup[1, j])) {
            bad++
            if (any(K[., j] :!= K[1, j])) {
                line = ind[j] + " has " + strofreal(K[1, j]) + " categories in group 1 (" +
                       gvar + " = " + gv[1] + ")"
                for (g = 2; g <= G; g++) {
                    line = line + (g < G ? ", " : " and ") + strofreal(K[g, j]) +
                           " in group " + strofreal(g) + " (" + gvar + " = " + gv[g] + ")"
                }
            }
            else {
                line = ind[j] + " takes the values " + sup[1, j] + " in group 1 (" +
                       gvar + " = " + gv[1] + ")"
                for (g = 2; g <= G; g++) {
                    line = line + (g < G ? ", " : " and ") + sup[g, j] +
                           " in group " + strofreal(g) + " (" + gvar + " = " + gv[g] + ")"
                }
            }
            errprintf("%s.\n", line)
        }
    }
    if (bad) {
        errprintf("Every invariance level past configural ties each indicator%sthresholds\n", char(39) + "s ")
        errprintf("across groups, which needs the same categories in every group.\n")
        errprintf("No model was fitted. invariance(configural) is not affected; partial\n")
        errprintf("threshold invariance is not available in this version.\n")
        st_local("gref", "1")
        return
    }
    if (level != "scalar") {
        return
    }

    // ---- 2. the binary scalar refusal, on the common K -------------------
    L   = dwlsmv_canon_lines(subinstr(canon, ";", char(10)))
    fac = J(1, 0, "")
    for (i = 1; i <= rows(L); i++) {
        w = tokens(L[i])
        if (w[1] == "factor") {
            fac = w[2..cols(w)]
            break
        }
    }
    bin = ""
    nb  = 0
    for (f = 1; f <= cols(fac); f++) {
        allb = 1
        has  = 0
        for (i = 1; i <= rows(L); i++) {
            w = tokens(L[i])
            // nested, not &: Mata's & does not short-circuit, and a free
            // lambda line has no fifth token
            if (w[1] == "group") {
                if (w[2] == "2") {
                    break
                }
            }
            if (w[1] == "lambda") {
                if (w[2] != fac[f]) {
                    continue
                }
                on = (w[4] == "free")
                if (w[4] == "fixed") {
                    on = (strtoreal(w[5]) != 0)
                }
                if (on) {
                    has = 1
                    for (j = 1; j <= cols(ind); j++) {
                        if (ind[j] == w[3] & K[1, j] != 2) {
                            allb = 0
                        }
                    }
                }
            }
        }
        if (has & allb) {
            nb++
            bin = bin + (nb > 1 ? ", " : "") + fac[f]
        }
    }
    if (nb) {
        errprintf("Factor(s) %s: every indicator has two categories.\n", bin)
        errprintf("For such a factor the scalar level is close to unidentified when the items are\n")
        errprintf("about equally difficult, and no validated rule yet tells safe data from harmful,\n")
        errprintf("so this version does not fit it. No model was fitted.\n")
        errprintf("You can fit invariance(configural), invariance(thresholds) or invariance(metric);\n")
        errprintf("factors with at least one indicator of three or more categories are not affected.\n")
        errprintf("invariance(scalar) is not available for an all-binary factor in this version.\n")
        exit(198)
    }
}

end
