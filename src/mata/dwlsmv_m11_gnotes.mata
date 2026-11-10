*! dwlsmv -- Module 11: the notes block of an invariance() fit (step 4)
*!
*!     void dwlsmv_m11_gnotes(string scalar canon, string scalar indlist,
*!                            string scalar gvar, string scalar gvalslist)
*!
*! Build step 4 (2026-10-12). S4-3 (decided 2026-10-12, user, (c) narrowed):
*! only notes whose groupwise meaning is established, each printed per group
*! under "Group g (gvar = v):" and only for a group where it fires, with
*! dwlsmv_m11_fit()'s G = 1 wording: pairs included at the bound under
*! keepbound; the one-empty-cell repair; polychoric correlations at or beyond
*! 0.99; admissibility (O15, decided 2026-10-31, user: group g on its own
*! scale, omega entering, by dwlsmv_mg_admissible(), worded by
*! dwlsmv_admissadvice()). Module 10's notes are deferred at G > 1 (they read
*! the pooled Delta or SEs). Then the d = 0 note, if any
*! (dwlsmv_m11_fit_g_note()). Nothing is stored and nothing is refused.
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_invariance_post.do.
*! Admissibility added 2026-10-31: tests/unit/test_mg_admissible.do.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_gnotes(string scalar canon, string scalar indlist, string scalar gvar,
                       string scalar gvalslist)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    pointer(real matrix) colvector Yg
    string rowvector ind, gv, fac, w
    string colvector L
    string matrix bnm
    real colvector rep, rho
    real scalar G, g, k, T, P, rpt, anyr, anyl, i
    string scalar txt, nl, line

    nl  = char(10)
    sts = *findexternal("dwlsmv_m11_st")
    Yg  = *findexternal("dwlsmv_m11_Yg")
    par = *findexternal("dwlsmv_m11_gpar")
    fit = *findexternal("dwlsmv_m11_gfit")
    ind = tokens(indlist)
    gv  = tokens(gvalslist)
    G   = rows(sts)
    // factor names, as dwlsmv_m11_post_g() reads them
    L   = dwlsmv_canon_lines(subinstr(canon, ";", char(10)))
    fac = J(1, 0, "")
    for (i = 1; i <= rows(L); i++) {
        w = tokens(L[i])
        if (w[1] == "factor") {
            fac = w[2..cols(w)]
            break
        }
    }
    anyr = 0
    anyl = 0
    for (g = 1; g <= G; g++) {
        T   = sts[g].map.T
        P   = sts[g].map.P
        bnm = J(P, 2, "")
        for (k = 1; k <= P; k++) {
            bnm[k, 1] = ind[sts[g].map.var_i[T + k]]
            bnm[k, 2] = ind[sts[g].map.var_j[T + k]]
        }
        txt = ""
        if (any(sts[g].bound)) {
            txt = txt + "Pairs included at the bound under keepbound:" + nl
            for (k = 1; k <= P; k++) {
                if (sts[g].bound[k]) {
                    txt = txt + "    " + bnm[k, 1] + "~~" + bnm[k, 2] + nl
                }
            }
        }
        rep = dwlsmv_emptycell_pairs(*Yg[g], sts[g].map)
        rpt = 0
        if (any(rep)) {
            rpt = 1
            txt = txt + "One-empty-cell repair applied to:" + nl
            for (k = 1; k <= P; k++) {
                if (rep[k]) {
                    txt = txt + "    " + bnm[k, 1] + "~~" + bnm[k, 2] + "   rho = " +
                          strtrim(strofreal(sts[g].s[T + k], "%9.4f")) + nl
                }
            }
        }
        rho  = sts[g].s[(T + 1)::(T + P)]
        line = ""
        for (k = 1; k <= P; k++) {
            if (!sts[g].bound[k]) {
                if (abs(rho[k]) >= 0.99) {
                    line = line + "    " + bnm[k, 1] + "~~" + bnm[k, 2] + "   rho = " +
                           strtrim(strofreal(rho[k], "%9.4f")) + nl
                }
            }
        }
        if (line != "") {
            txt = txt + "Polychoric correlations at or beyond 0.99 in absolute value:" + nl + line
        }
        txt = txt + dwlsmv_admissadvice(dwlsmv_mg_admissible(fit.theta, par, g, sts[g].map),
                                        ind', fac')
        if (txt != "") {
            printf("{txt}Group %f (%s = %s):\n", g, gvar, gv[g])
            dwlsmv_m11_gnotes_print(txt)
        }
        anyr = anyr | rpt
        anyl = anyl | (line != "")
    }
    // the explanations once, after every group, in dwlsmv_m11_fit()'s G = 1 words
    if (anyr) {
        dwlsmv_m11_gnotes_print(
          "0.5 was added to the empty cell and its diagonal partner and taken from" + nl +
          "the other two; margins and thresholds are unchanged. lavaan and Mplus apply" + nl +
          "the same repair and the same Gamma. Its effect on standard errors and T3 at" + nl +
          "small N is not established." + nl)
    }
    if (anyl) {
        dwlsmv_m11_gnotes_print(
          "A correlation this close to the bound has several possible readings: near-" + nl +
          "duplicate (redundant) items, strong skew collapsing the table, or too few" + nl +
          "observations. The data cannot tell them apart; this is not a diagnosis." + nl)
    }
    dwlsmv_m11_fit_g_note()
}

// one printf per line (a %s holding several char(10)s lets SMCL wrap)
void dwlsmv_m11_gnotes_print(string scalar s0)
{
    string scalar s
    real scalar k

    s = s0
    while ((k = strpos(s, char(10))) > 0) {
        printf("{txt}%s\n", substr(s, 1, k - 1))
        s = substr(s, k + 1, .)
    }
    if (s != "") {
        printf("{txt}%s\n", s)
    }
}

end
