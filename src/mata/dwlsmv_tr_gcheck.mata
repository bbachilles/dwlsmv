*! dwlsmv -- Module 11: the invariance() checks on the configural text
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_translate_groups.do.
*!
*!   dwlsmv_tr_gcheck(text, gname, gvals)
*!
*! Called by dwlsmv_translate_g on its final pass. Refuses at rc 198, in this
*! order (most specific first, as G7), each with every instance listed:
*!   1  a loading fixed at 1 (decided 2026-10-12, point 8's second message);
*!   2  a structural path between factors (spec 4.8 G5: refused with
*!      invariance levels in 0.3.0; Module 11 first, Module 3 the backstop);
*!   3  a later group that does not mirror group 1 (decided 2026-10-12,
*!      "Group-specific paths under invariance()" (a); dwlsmv_levelgen's rc 3498
*!      is the backstop). A loading fixed at 0 counts as absent, as sem's
*!      e(Cns) records a path written for another group only (probe case 5.4).
*! Groups are named by position and value (O13), e.g. "group 2 (grp = 1)".
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_tr_gcheck(string colvector text, string scalar gname, real rowvector gvals)
{
    real colvector at
    string colvector L1, Lg, msg, seen
    string rowvector t
    real scalar G, g, i, n
    string scalar key, s1, sg

    G  = cols(gvals)
    at = J(G + 1, 1, rows(text) + 1)
    for (i = 1; i <= rows(text); i++) {
        t = tokens(text[i])
        if (t[1] == "group") {
            at[strtoreal(t[2])] = i
        }
    }

    // 1: a loading fixed at 1, one message per factor
    seen = J(0, 1, "")
    for (i = 1; i <= rows(text); i++) {
        t = tokens(text[i])
        if (t[1] == "metric" & cols(t) == 4) {
            t = ("lambda", t[3], t[4], "fixed", "1")
        }
        if (t[1] == "lambda") {
            if (t[4] == "fixed") {
                if (strtoreal(t[5]) == 1 & !anyof(seen, t[2])) {
                    seen = seen \ t[2]
                    errprintf("Factor %s has a loading fixed at 1 (%s -> %s@1). With invariance(), each factor's\n",
                              t[2], t[2], t[3])
                    errprintf("scale is set by metric(unitvariance) instead; remove the fixed loading and\n")
                    errprintf("specify metric(unitvariance). No model was fitted.\n")
                }
            }
        }
    }
    if (rows(seen)) {
        exit(198)
    }

    // 2: structural paths
    msg = J(0, 1, "")
    for (i = 1; i <= rows(text); i++) {
        t = tokens(text[i])
        if (t[1] == "beta") {
            if (!anyof(msg, t[2] + " <- " + t[3])) {
                msg = msg \ (t[2] + " <- " + t[3])
            }
        }
    }
    if (rows(msg)) {
        errprintf("Structural paths between factors (%s) are not available with\n", invtokens(msg', ", "))
        errprintf("invariance() in this version. No model was fitted.\n")
        exit(198)
    }

    // 3: every later group mirrors group 1
    L1  = dwlsmv_tr_glines(text, at[1] + 1, at[2] - 1)
    msg = J(0, 1, "")
    for (g = 2; g <= G; g++) {
        Lg = dwlsmv_tr_glines(text, at[g] + 1, at[g + 1] - 1)
        n  = rows(L1)
        for (i = 1; i <= rows(L1) + rows(Lg); i++) {
            key = (i <= n ? L1[i] : Lg[i - n])
            key = invtokens(tokens(key)[1..3])
            s1 = dwlsmv_tr_gstat(L1, key)
            sg = dwlsmv_tr_gstat(Lg, key)
            if (s1 == sg) {
                continue
            }
            if (i > n & s1 != "") {
                continue
            }
            msg = msg \ dwlsmv_tr_gdiff(key, s1, sg, g, gname, gvals)
        }
    }
    if (rows(msg)) {
        for (i = 1; i <= rows(msg); i++) {
            errprintf("%s\n", msg[i])
        }
        errprintf("With invariance(), every group must have the same model. No model was fitted.\n")
        exit(198)
    }
}

// a group's lines from the third on (after factor and ind), with group 1's
// metric line written as later groups write it, and loadings fixed at 0 dropped
string colvector dwlsmv_tr_glines(string colvector text, real scalar a, real scalar b)
{
    string colvector out
    string rowvector t
    real scalar i

    out = J(0, 1, "")
    for (i = a + 2; i <= b; i++) {
        t = tokens(text[i])
        if (t[1] == "metric") {
            out = out \ ("psi " + t[3] + " " + t[3] + " fixed 1")
        }
        else if (t[1] == "lambda" & t[4] == "fixed") {
            if (strtoreal(t[5]) != 0) {
                out = out \ text[i]
            }
        }
        else {
            out = out \ text[i]
        }
    }
    return(out)
}

// the status of a key ("lambda F y4") in a group: "" if absent, else the rest
string scalar dwlsmv_tr_gstat(string colvector L, string scalar key)
{
    real scalar i
    string rowvector t

    for (i = 1; i <= rows(L); i++) {
        t = tokens(L[i])
        if (invtokens(t[1..3]) == key) {
            return(invtokens(t[4..cols(t)]))
        }
    }
    return("")
}

string scalar dwlsmv_tr_gdiff(string scalar key, string scalar s1, string scalar sg,
    real scalar g, string scalar gname, real rowvector gvals)
{
    string rowvector t
    string scalar what, n1, ng

    t = tokens(key)
    if (t[1] == "lambda") {
        what = "Path " + t[2] + " -> " + t[3]
    }
    else if (t[2] == t[3]) {
        what = "The variance of " + t[2]
    }
    else {
        what = "The covariance of " + t[2] + " and " + t[3]
    }
    n1 = "group 1 (" + gname + " = " + strofreal(gvals[1]) + ")"
    ng = "group " + strofreal(g) + " (" + gname + " = " + strofreal(gvals[g]) + ")"
    if (s1 == "") {
        return(what + " is in " + ng + " but not in " + n1 + ".")
    }
    if (sg == "") {
        return(what + " is in " + n1 + " but not in " + ng + ".")
    }
    return(what + " is " + dwlsmv_tr_gword(sg) + " in " + ng + " but " +
           dwlsmv_tr_gword(s1) + " in " + n1 + ".")
}

string scalar dwlsmv_tr_gword(string scalar s)
{
    string rowvector t

    t = tokens(s)
    if (t[1] == "fixed") {
        return("fixed at " + t[2])
    }
    return(s)
}

end
