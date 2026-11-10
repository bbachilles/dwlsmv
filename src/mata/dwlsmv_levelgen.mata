*! dwlsmv -- the invariance level generator (Module 12)
*!
*!     string scalar dwlsmv_levelgen(string scalar conf, real matrix K,
*!                                   string scalar level)
*!
*! Spec 4.8 G8 (decided 2026-10-08/09, location and contract 2026-10-12, user).
*! A canonical-to-canonical function: from the configural canonical text, the
*! category counts K and a level name it returns that level's canonical text,
*! ties included, in G8's normative emission order. No data, no fit. Called by
*! Module 11's invariance() and, later, by check 8's generated-chain exemption,
*! which regenerates a level and compares it exactly with the stored text.
*!
*! INPUT CONTRACT (L1 (b)). conf is the full multi-group configural text. It is
*! validated: every later group must mirror group 1 (same factors, indicators
*! and order; the same lambda and off-diagonal psi lines, status, value and
*! start, Q-5); group 1 scales every factor by "metric unitvariance F"
*! and later groups by "psi F F fixed 1". Nothing else
*! may differ. No beta, nu, omega, alpha or tie line is accepted (beta is never
*! in a generated chain, scope (B); the others belong to levels, not to the
*! configural input). Anything else aborts, rc 3498.
*!
*! K is G x p, rows the groups, columns the indicators in ind-line order. At any
*! level past configural an indicator's K must be equal in every group (G3: a
*! tie of thresholds across different category counts is an incompatible
*! restriction; Module 11 refuses it first at rc 198, this is the backstop for
*! a direct caller). Polytomous means K > 2.
*!
*! LEVELS (L2 (a)): configural, thresholds, metric, scalar, spelled in full.
*! In each later group (G8's table; metric = loadings, scalar = intercepts):
*!   configural  psi F F fixed 1;  no nu, omega, alpha;           no ties
*!   thresholds  psi F F fixed 1;  nu all; omega polytomous;      tau ties
*!   metric      psi F F free iff F has a polytomous indicator, else fixed 1;
*!               nu all; omega all;                               tau, lambda
*!   scalar      psi F F free;     omega all; alpha every factor; tau, lambda
*! A factor's indicators are its group-1 lambda lines free or fixed nonzero.
*! Only free lambda lines are tied (a fixed row is never tied, G2). Scalar with
*! an all-binary factor IS emitted: the 0.3.0 refusal of it is Module 11's
*! policy (T0, T2), not this function's mathematics.
*!
*! EMISSION ORDER (normative; spec 4.8 G8). groups G; then each group in order:
*! group g, one factor line, one ind line, the lambda lines in group 1's
*! input order; then group 1's metric unitvariance F lines in factor order, or
*! a later group's psi F F lines in factor order; then the off-diagonal psi
*! lines in group 1's input order; nu (indicator order), omega (indicator
*! order), alpha (factor order). Then the tie lines: tau, by indicator then
*! threshold; then lambda, in group 1's lambda order. Tokens are single-spaced,
*! comments and blank lines dropped, every line ended by char(10). Copied lines
*! keep their tokens as written (values are never reformatted), so a configural
*! text that is the generator's own output regenerates byte for byte.
*!
*! Executed 2026-10-12: tests/unit/test_levelgen.do green; mutants M1-M5 red.
*!
*! Load order: dwlsmv_canon_lines.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Lines of a vector joined with char(10) after each.
string scalar dwlsmv_lg_join(string colvector L)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= rows(L); i++) {
        s = s + L[i] + char(10)
    }
    return(s)
}

// Position of s in v, or 0.
real scalar dwlsmv_lg_pos(string rowvector v, string scalar s)
{
    real scalar i

    for (i = 1; i <= cols(v); i++) {
        if (v[i] == s) {
            return(i)
        }
    }
    return(0)
}

void dwlsmv_lg_fail(string scalar msg)
{
    errprintf("dwlsmv_levelgen: %s\n", msg)
    _error("dwlsmv_levelgen: the configural text or K is not valid for a level")
}

string scalar dwlsmv_levelgen(string scalar conf, real matrix K, string scalar level)
{
    string colvector raw, lam, pcv, out, srt0, srt1
    string rowvector w, fac, ind
    string scalar s
    real scalar G, ng, i, j, k, g, f, cur, p, nf, hasf, poly
    real rowvector seen, fpoly
    pointer(string colvector) colvector Plam, Ppcv, Pmet, Ppdg
    pointer(string rowvector) colvector Pfac, Pind

    if (level != "configural" & level != "thresholds" & level != "metric" &
        level != "scalar") {
        errprintf("dwlsmv_levelgen: unknown level %s; use configural, thresholds, metric or scalar.\n",
                  level)
        _error("dwlsmv_levelgen: unknown level")
    }

    raw = dwlsmv_canon_lines(conf)
    for (i = 1; i <= rows(raw); i++) {
        raw[i] = invtokens(tokens(raw[i]))
    }

    // the groups line: exactly one, G an integer >= 2
    G  = .
    ng = 0
    for (i = 1; i <= rows(raw); i++) {
        w = tokens(raw[i])
        if (w[1] == "groups") {
            ng++
            if (cols(w) != 2) dwlsmv_lg_fail("a groups line must be groups G.")
            G = strtoreal(w[2])
        }
    }
    if (ng != 1) dwlsmv_lg_fail("the configural text needs exactly one groups line.")
    if (missing(G) | G < 2 | G != floor(G)) {
        dwlsmv_lg_fail("groups must give an integer of at least 2.")
    }

    Plam = J(G, 1, NULL)
    Ppcv = J(G, 1, NULL)
    Pmet = J(G, 1, NULL)
    Ppdg = J(G, 1, NULL)
    Pfac = J(G, 1, NULL)
    Pind = J(G, 1, NULL)
    for (g = 1; g <= G; g++) {
        Plam[g] = &J(0, 1, "")
        Ppcv[g] = &J(0, 1, "")
        Pmet[g] = &J(0, 1, "")
        Ppdg[g] = &J(0, 1, "")
        Pfac[g] = &J(1, 0, "")
        Pind[g] = &J(1, 0, "")
    }
    seen = J(1, G, 0)
    cur  = 0
    for (i = 1; i <= rows(raw); i++) {
        w = tokens(raw[i])
        if (w[1] == "groups") {
            continue
        }
        if (w[1] == "group") {
            g = (cols(w) == 2 ? strtoreal(w[2]) : .)
            if (missing(g) | g < 1 | g > G | g != floor(g)) {
                dwlsmv_lg_fail("a group line must name a group 1.." + strofreal(G) + ": " + raw[i])
            }
            if (seen[g]) dwlsmv_lg_fail("group " + strofreal(g) + " appears twice.")
            seen[g] = 1
            cur = g
            continue
        }
        if (w[1] == "tie" | w[1] == "nu" | w[1] == "omega" | w[1] == "alpha" |
            w[1] == "beta") {
            dwlsmv_lg_fail("a configural input has no " + w[1] +
                           " line (levels add nu, omega, alpha and ties; beta is never in a generated chain): " +
                           raw[i])
        }
        if (cur == 0) dwlsmv_lg_fail("a line lies outside any group section: " + raw[i])
        if ((w[1] == "factor" | w[1] == "ind") & cols(w) < 2) {
            dwlsmv_lg_fail("a factor or ind line names nothing: " + raw[i])
        }
        if (w[1] == "factor") {
            *Pfac[cur] = *Pfac[cur], w[|2 \ cols(w)|]
        }
        else if (w[1] == "ind") {
            *Pind[cur] = *Pind[cur], w[|2 \ cols(w)|]
        }
        else if (w[1] == "lambda") {
            *Plam[cur] = *Plam[cur] \ raw[i]
        }
        else if (w[1] == "metric") {
            *Pmet[cur] = *Pmet[cur] \ raw[i]
        }
        else if (w[1] == "psi") {
            if (cols(w) < 4) dwlsmv_lg_fail("a psi line is incomplete: " + raw[i])
            if (w[2] == w[3]) *Ppdg[cur] = *Ppdg[cur] \ raw[i]
            else              *Ppcv[cur] = *Ppcv[cur] \ raw[i]
        }
        else {
            dwlsmv_lg_fail("unrecognized directive for a configural input: " + raw[i])
        }
    }
    if (sum(seen) != G) dwlsmv_lg_fail("every group 1.." + strofreal(G) + " needs a section.")

    fac = *Pfac[1]
    ind = *Pind[1]
    nf  = cols(fac)
    p   = cols(ind)
    if (nf == 0 | p == 0) dwlsmv_lg_fail("group 1 declares no factor or no indicator.")
    lam = *Plam[1]
    pcv = *Ppcv[1]

    // group 1: one metric unitvariance line per factor, no psi variance line
    if (rows(*Ppdg[1]) > 0) {
        dwlsmv_lg_fail("group 1 scales its factors by metric unitvariance, not by a psi line: " +
                       (*Ppdg[1])[1])
    }
    if (rows(*Pmet[1]) != nf) dwlsmv_lg_fail("group 1 needs one metric unitvariance line per factor.")
    for (f = 1; f <= nf; f++) {
        if (!anyof(*Pmet[1], "metric unitvariance " + fac[f])) {
            dwlsmv_lg_fail("group 1 lacks metric unitvariance " + fac[f] + ".")
        }
    }

    // later groups mirror group 1
    srt0 = sort(lam, 1) \ "|" \ sort(pcv, 1)
    for (g = 2; g <= G; g++) {
        if (*Pfac[g] != fac | *Pind[g] != ind) {
            dwlsmv_lg_fail("group " + strofreal(g) + " must declare the factors and indicators of group 1, in the same order.")
        }
        if (rows(*Pmet[g]) > 0) dwlsmv_lg_fail("a metric line lies in group " + strofreal(g) + ".")
        if (rows(*Ppdg[g]) != nf) {
            dwlsmv_lg_fail("group " + strofreal(g) + " needs psi F F fixed 1 for every factor.")
        }
        for (f = 1; f <= nf; f++) {
            if (!anyof(*Ppdg[g], "psi " + fac[f] + " " + fac[f] + " fixed 1")) {
                dwlsmv_lg_fail("group " + strofreal(g) + " lacks psi " + fac[f] + " " + fac[f] + " fixed 1.")
            }
        }
        srt1 = sort(*Plam[g], 1) \ "|" \ sort(*Ppcv[g], 1)
        if (rows(srt1) != rows(srt0)) {
            dwlsmv_lg_fail("group " + strofreal(g) + " does not mirror group 1 (loadings or factor covariances differ).")
        }
        if (srt1 != srt0) {
            dwlsmv_lg_fail("group " + strofreal(g) + " does not mirror group 1 (loadings or factor covariances differ).")
        }
    }

    // K: G x p, integers >= 2; equal across groups past configural
    if (rows(K) != G | cols(K) != p) {
        dwlsmv_lg_fail("K must be " + strofreal(G) + " x " + strofreal(p) + " (groups by indicators).")
    }
    if (missing(K)) dwlsmv_lg_fail("K has a missing entry.")
    if (any(K :< 2) | any(K :!= floor(K))) dwlsmv_lg_fail("every K must be an integer of at least 2.")
    if (level != "configural") {
        for (j = 1; j <= p; j++) {
            for (g = 2; g <= G; g++) {
                if (K[g, j] != K[1, j]) {
                    dwlsmv_lg_fail(ind[j] + " has " + strofreal(K[1, j]) + " categories in group 1 and " +
                                   strofreal(K[g, j]) + " in group " + strofreal(g) +
                                   ": its thresholds cannot be tied, so no level past configural exists.")
                }
            }
        }
    }

    // a factor is polytomous if one of its indicators (lambda free or fixed nonzero) is
    fpoly = J(1, nf, 0)
    for (i = 1; i <= rows(lam); i++) {
        w = tokens(lam[i])
        if (cols(w) < 4) dwlsmv_lg_fail("a lambda line is incomplete: " + lam[i])
        f = dwlsmv_lg_pos(fac, w[2])
        j = dwlsmv_lg_pos(ind, w[3])
        if (f == 0 | j == 0) dwlsmv_lg_fail("a lambda line names an undeclared factor or indicator: " + lam[i])
        if (w[4] == "free") {
            hasf = 1
        }
        else if (w[4] == "fixed" & cols(w) == 5) {
            hasf = (strtoreal(w[5]) != 0)
        }
        else {
            dwlsmv_lg_fail("a lambda line must be free [start] or fixed value: " + lam[i])
        }
        if (hasf & K[1, j] > 2) fpoly[f] = 1
    }

    // emission
    out = "groups " + strofreal(G)
    for (g = 1; g <= G; g++) {
        out = out \ ("group " + strofreal(g)) \ ("factor " + invtokens(fac)) \ ("ind " + invtokens(ind))
        out = out \ lam
        if (g == 1) {
            for (f = 1; f <= nf; f++) {
                out = out \ ("metric unitvariance " + fac[f])
            }
        }
        else {
            for (f = 1; f <= nf; f++) {
                s = "psi " + fac[f] + " " + fac[f]
                if (level == "scalar" | (level == "metric" & fpoly[f])) {
                    out = out \ (s + " free")
                }
                else {
                    out = out \ (s + " fixed 1")
                }
            }
        }
        out = out \ pcv
        if (g > 1) {
            if (level == "thresholds" | level == "metric") {
                for (j = 1; j <= p; j++) {
                    out = out \ ("nu " + ind[j] + " free")
                }
            }
            if (level != "configural") {
                for (j = 1; j <= p; j++) {
                    poly = (K[1, j] > 2)
                    if (level != "thresholds" | poly) {
                        out = out \ ("omega " + ind[j] + " free")
                    }
                }
            }
            if (level == "scalar") {
                for (f = 1; f <= nf; f++) {
                    out = out \ ("alpha " + fac[f] + " free")
                }
            }
        }
    }
    if (level != "configural") {
        for (j = 1; j <= p; j++) {
            for (k = 1; k <= K[1, j] - 1; k++) {
                s = "tie"
                for (g = 1; g <= G; g++) {
                    s = s + " tau " + strofreal(g) + " " + ind[j] + " " + strofreal(k)
                }
                out = out \ s
            }
        }
    }
    if (level == "metric" | level == "scalar") {
        for (i = 1; i <= rows(lam); i++) {
            w = tokens(lam[i])
            if (w[4] != "free") continue
            s = "tie"
            for (g = 1; g <= G; g++) {
                s = s + " lambda " + strofreal(g) + " " + w[2] + " " + w[3]
            }
            out = out \ s
        }
    }
    return(dwlsmv_lg_join(out))
}

end
