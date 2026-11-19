*! dwlsmv -- Module 11: the compact grouped coefficient table (release A, owner 2026-11-17 and 2026-11-18)
*!
*! Drafted by a read-only agent 2026-11-18; run by the main session the same day (tests/unit/test_ctab.do;
*! fixes on the first run: e(canon) split on ";", C2 at every G, the driver named after the file).
*!
*!     void dwlsmv_m11_ctab(string scalar bn, string scalar vn, string scalar cn, real scalar showthr,
*!                         real scalar std, real scalar level, real scalar ls, string scalar tn)
*!         The driver dwlsmv.ado calls. bn, vn, cn name Stata matrices: b (1 x k, with its column
*!         stripe), V ("" when there is no inference to show) and the fixed-value rows ("" for none).
*!         Reads e(metric), e(canon), e(indicators), e(groupvar), e(invariance) and e(K) itself, so
*!         that no model text passes through a Stata macro. Prints the table, writes the full
*!         9 x k matrix r(table) into the Stata matrix tn, and sets the locals ct_nshow, ct_nomit,
*!         ct_nthr (rows hidden as thresholds) and ct_nlev (rows omitted as metric fixings or level
*!         conventions at G > 1).
*!
*!     void dwlsmv_ct_plan(stp, b, V, hasv, C, hasc, metric, canon, ind, gv, gvar, showthr, std, level,
*!                         ptxt, pnum, otxt, oidx, tab)
*!         The pure planner: no printing, no Stata state. The row plan:
*!           ptxt (n x 3)  heading of the block, group title ("" at G = 1 and for thresholds), label
*!           pnum (n x 11) b index, kind (1 free, 2 marker, 3 constrained), est, se, z, p, ll, ul,
*!                         block number, group key (G + 1 for thresholds), section rank
*!         the omitted rows, otxt (m x 2) the parameter name eq:name and the reason (thresholds,
*!         metric, level, standardized), oidx (m x 1) their b index, and tab (9 x k), r(table)'s
*!         layout (b se z pvalue ll ul df crit eform) for EVERY column of b.
*!         A row's numbers are read from tab at the row's own b index, so reordering the rows for
*!         display cannot move a standard error away from its parameter.
*!
*!     void dwlsmv_ct_print(ptxt, pnum, hasv, level, ls)
*!         The printer. 78 columns at linesize 78: label column (16, or the longest label plus 2
*!         up to linesize - 62), then Estimate, Std. err., z, P>|z| and the interval.
*!
*! Row labels. Within a block the stripe name without the equation the block heading already
*! states: a loading y1 under "Loadings on F"; F2 <- F1 under "Structural"; F <- x, y4 <- x under
*! "Covariates"; var(F), cov(F1,F2), scale(y), mean(F) as posted; a threshold cutk, followed by
*! " (Group g)" at G > 1. The leading "/" of the stripe is not shown.
*!
*! Which fixed rows are omitted (the planner's rule, one place):
*!   G = 1, e(metric) unitvariance: the factor variance fixed at 1 of each factor named by e(canon)'s
*!       "metric unitvariance F" line (var(F), or var(e.F) for an endogenous factor): reason metric.
*!   G > 1: a factor variance fixed at 1 (group 1: reason metric; later groups: level) and the nu
*!       at 0, omega at 1 and alpha at 0 that the level fixes: reason level. Nothing else is
*!       omitted: a loading, covariance or other row fixed by the user is shown as constrained.
*!   standardized (C2 (a), owner 2026-11-18): a constant row that is an exogenous var(F), at every G,
*!       and at G > 1 a scale, an intercept or a mean: reason standardized. Any other constant row
*!       (a user-fixed covariance, say) is shown as constrained.
*!   thresholds: every threshold row when showthr is 0: reason thresholds.
*! The marker (G = 1, raw): the loading e(canon)'s "metric marker F y" line names (without e(canon),
*! the first loading fixed at 1 of each factor under e(metric) marker). Any other fixed row is
*! constrained, including a user loading fixed at 1 elsewhere.
*!
*! Malformed input (a stripe that is not e(b)'s, a fixed-value row that does not select one
*! parameter, a missing value, a V or C of the wrong shape) is rc 3498: errprintf() says what,
*! _error() is short.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_ct_bad(string scalar why)
{
    errprintf("dwlsmv: the coefficient table could not be built: %s\n", why)
    _error("dwlsmv: internal error, coefficient table input is malformed")
}

// the first position of s in v, or 0
real scalar dwlsmv_ct_pos(string rowvector v, string scalar s)
{
    real scalar i

    for (i = 1; i <= cols(v); i++) {
        if (v[i] == s) {
            return(i)
        }
    }
    return(0)
}

// "head#tag" (the tag after the #): the group, or 0; head returns the part before the #.
// At one group (tags has one column) there must be no #.
real scalar dwlsmv_ct_tagafter(string scalar nm, string rowvector tags, string scalar head)
{
    real scalar p, g

    head = nm
    p = strpos(nm, "#")
    if (cols(tags) == 1) {
        if (p == 0) {
            return(1)
        }
        return(0)
    }
    if (p == 0) {
        return(0)
    }
    head = substr(nm, 1, p - 1)
    for (g = 1; g <= cols(tags); g++) {
        if (substr(nm, p + 1, .) == tags[g]) {
            return(g)
        }
    }
    return(0)
}

// "tag#c.rest": the group, or 0; rest returns what follows the c.
real scalar dwlsmv_ct_tagbefore(string scalar nm, string rowvector tags, string scalar rest)
{
    real scalar p, g

    rest = nm
    p = strpos(nm, "#c.")
    if (p == 0) {
        return(0)
    }
    rest = substr(nm, p + 3, .)
    for (g = 1; g <= cols(tags); g++) {
        if (substr(nm, 1, p - 1) == tags[g]) {
            return(g)
        }
    }
    return(0)
}

// r(table): b se z pvalue ll ul df crit eform for every column of b. Normal theory, as ereturn
// display computes it (dwlsmv posts no e(df_r)): z = b/se, p = 2 normal(-|z|), b -/+ crit se.
// se is the square root of V's diagonal at the column's own index. A zero or missing se gives
// missing z, p, ll and ul. Without V (hasv 0) only b, df, crit and eform are filled.
real matrix dwlsmv_ct_tab(real rowvector b, real matrix V, real scalar hasv, real scalar level)
{
    real matrix T
    real scalar nb, i, crit, se, z

    nb   = cols(b)
    crit = invnormal(.5 + level / 200)
    T    = J(9, nb, .)
    for (i = 1; i <= nb; i++) {
        T[1, i] = b[i]
        T[8, i] = crit
        T[9, i] = 0
        if (hasv) {
            se      = sqrt(V[i, i])
            T[2, i] = se
            if (se > 0 & se < .) {
                z       = b[i] / se
                T[3, i] = z
                T[4, i] = 2 * normal(-abs(z))
                T[5, i] = b[i] - crit * se
                T[6, i] = b[i] + crit * se
            }
        }
    }
    return(T)
}

void dwlsmv_ct_plan(string matrix stp, real rowvector b, real matrix V, real scalar hasv,
                    real matrix C, real scalar hasc, string scalar metric, string scalar canon,
                    string rowvector ind, string rowvector gv, string scalar gvar,
                    real scalar showthr, real scalar std, real scalar level,
                    string matrix ptxt, real matrix pnum, string matrix otxt, real colvector oidx,
                    real matrix tab)
{
    real scalar nb, G, i, j, k, g, f, it, nz, cj, p, nf, nsh, nom, blk, bi
    real scalar fv, isfac
    real rowvector fixed, fixval, kind, grp, sub, rk, omit, rkd, seenmk
    real colvector o
    real matrix K
    string scalar eq, nm, head, inner, rest, reason
    string rowvector fac, tags, mkf, mky, uvf, w
    string colvector lines, lab, why

    // ---- shapes ------------------------------------------------------------
    nb = cols(b)
    if (nb == 0 | rows(b) != 1) {
        dwlsmv_ct_bad("e(b) is not a nonempty row vector")
    }
    if (rows(stp) != nb | cols(stp) != 2) {
        dwlsmv_ct_bad("the column stripe does not fit e(b)")
    }
    if (missing(b)) {
        dwlsmv_ct_bad("e(b) has a missing value")
    }
    if (hasv) {
        if (rows(V) != nb | cols(V) != nb) {
            dwlsmv_ct_bad("e(V) does not fit e(b)")
        }
        if (missing(diagonal(V))) {
            dwlsmv_ct_bad("e(V) has a missing diagonal entry")
        }
    }
    if (level <= 0 | level >= 100 | level == .) {
        dwlsmv_ct_bad("the confidence level is not between 0 and 100")
    }
    G = cols(gv)
    if (G < 2) {
        G = 1
    }
    if (G > 1 & gvar == "") {
        dwlsmv_ct_bad("group values without the group variable")
    }
    tags = J(1, G, "")
    if (G > 1) {
        for (g = 1; g <= G; g++) {
            if (gv[g] == "") {
                dwlsmv_ct_bad("an empty group value")
            }
            tags[g] = gv[g] + (g == 1 ? "bn." : ".") + gvar
        }
    }
    for (i = 1; i <= nb; i++) {
        for (j = 1; j < i; j++) {
            if (stp[i, 1] == stp[j, 1] & stp[i, 2] == stp[j, 2]) {
                dwlsmv_ct_bad("a parameter appears twice in e(b): " + stp[i, 1] + ":" + stp[i, 2])
            }
        }
    }

    // ---- fixed rows --------------------------------------------------------
    fixed  = J(1, nb, 0)
    fixval = J(1, nb, .)
    if (hasc) {
        if (cols(C) != nb + 1) {
            dwlsmv_ct_bad("the fixed-value rows do not fit e(b)")
        }
        if (missing(C)) {
            dwlsmv_ct_bad("the fixed-value rows have a missing value")
        }
        for (k = 1; k <= rows(C); k++) {
            nz = 0
            cj = 0
            for (j = 1; j <= nb; j++) {
                if (C[k, j] != 0) {
                    nz++
                    cj = j
                }
            }
            if (nz != 1) {
                dwlsmv_ct_bad("a fixed-value row does not select exactly one parameter")
            }
            if (fixed[cj]) {
                dwlsmv_ct_bad("a parameter is fixed twice")
            }
            fixed[cj]  = 1
            fixval[cj] = C[k, nb + 1] / C[k, cj]
        }
    }

    // ---- the numbers, by b index -------------------------------------------
    tab = dwlsmv_ct_tab(b, V, hasv, level)

    // ---- factors, from the variance rows (group 1's at G > 1); an indicator's residual variance
    // var(e.y) of the standardized table is not a factor
    fac = J(1, 0, "")
    for (i = 1; i <= nb; i++) {
        if (stp[i, 1] == "/") {
            g = dwlsmv_ct_tagafter(stp[i, 2], tags, head)
            if (g == 1) {
                if (substr(head, 1, 4) == "var(" & substr(head, strlen(head), 1) == ")") {
                    inner = substr(head, 5, strlen(head) - 5)
                    if (substr(inner, 1, 2) == "e.") {
                        inner = substr(inner, 3, .)
                    }
                    if (!anyof(ind, inner) & !anyof(fac, inner)) {
                        fac = fac, inner
                    }
                }
            }
        }
    }
    nf = cols(fac)
    if (nf == 0) {
        dwlsmv_ct_bad("no factor variance row in e(b)")
    }

    // ---- the metric, from e(canon) at one group ------------------------------
    mkf = J(1, 0, "")
    mky = J(1, 0, "")
    uvf = J(1, 0, "")
    if (G == 1 & canon != "") {
        // e(canon) is posted with ";" between lines (one Stata macro); dwlsmv_canon_lines splits on
        // newlines
        lines = dwlsmv_canon_lines(subinstr(canon, ";", char(10), .))
        for (k = 1; k <= rows(lines); k++) {
            w = tokens(lines[k])
            if (cols(w) >= 3) {
                if (w[1] == "metric" & w[2] == "marker" & cols(w) == 4) {
                    mkf = mkf, w[3]
                    mky = mky, w[4]
                }
                if (w[1] == "metric" & w[2] == "unitvariance" & cols(w) == 3) {
                    uvf = uvf, w[3]
                }
            }
        }
    }
    else if (G == 1 & metric == "unitvariance") {
        uvf = fac
    }
    seenmk = J(1, nf, 0)

    // ---- classify every column ---------------------------------------------
    // kind: 1 loading, 2 path, 3 gamma, 4 kappa, 5 variance, 6 covariance, 7 threshold,
    //       8 intercept (nu), 9 scale (omega), 10 mean (alpha)
    kind = J(1, nb, 0)
    grp  = J(1, nb, 1)
    sub  = J(1, nb, 0)
    rk   = J(1, nb, 0)
    lab  = J(nb, 1, "")
    for (i = 1; i <= nb; i++) {
        eq = stp[i, 1]
        nm = stp[i, 2]
        if (eq == "/") {
            g = dwlsmv_ct_tagafter(nm, tags, head)
            if (g == 0) {
                dwlsmv_ct_bad("a row without its group tag: " + eq + ":" + nm)
            }
            grp[i] = g
            lab[i] = head
            if (substr(head, strlen(head), 1) != ")") {
                dwlsmv_ct_bad("an unreadable row: " + eq + ":" + nm)
            }
            if (substr(head, 1, 4) == "var(") {
                kind[i] = 5
                rk[i]   = 4
            }
            else if (substr(head, 1, 4) == "cov(") {
                kind[i] = 6
                rk[i]   = 5
            }
            else if (substr(head, 1, 6) == "scale(") {
                kind[i] = 9
                rk[i]   = 7
            }
            else if (substr(head, 1, 5) == "mean(") {
                kind[i] = 10
                rk[i]   = 8
            }
            else {
                dwlsmv_ct_bad("an unreadable row: " + eq + ":" + nm)
            }
            if (G == 1 & (kind[i] == 9 | kind[i] == 10)) {
                dwlsmv_ct_bad("a scale or mean row at one group: " + eq + ":" + nm)
            }
        }
        else if (substr(eq, 1, 1) == "/") {
            it = dwlsmv_ct_pos(ind, substr(eq, 2, .))
            if (it == 0) {
                dwlsmv_ct_bad("a threshold of a variable that is not an indicator: " + eq)
            }
            kind[i] = 7
            rk[i]   = 9
            sub[i]  = it
            if (G == 1) {
                if (strpos(nm, "#") > 0) {
                    dwlsmv_ct_bad("a group tag at one group: " + eq + ":" + nm)
                }
                grp[i] = 1
                lab[i] = nm
            }
            else {
                g = dwlsmv_ct_tagbefore(nm, tags, rest)
                if (g == 0) {
                    dwlsmv_ct_bad("a threshold without its group tag: " + eq + ":" + nm)
                }
                grp[i] = g
                lab[i] = rest + " (Group " + strofreal(g) + ")"
            }
        }
        else if (anyof(ind, eq)) {
            if (G == 1) {
                f = dwlsmv_ct_pos(fac, nm)
                if (f > 0) {
                    kind[i] = 1
                    rk[i]   = 1
                    sub[i]  = f
                    lab[i]  = eq
                }
                else {
                    if (strpos(nm, "#") > 0) {
                        dwlsmv_ct_bad("a group tag at one group: " + eq + ":" + nm)
                    }
                    kind[i] = 4
                    rk[i]   = 3
                    lab[i]  = eq + " <- " + nm
                }
            }
            else {
                g = dwlsmv_ct_pos(tags, nm)
                if (g > 0) {
                    kind[i] = 8
                    rk[i]   = 6
                    grp[i]  = g
                    lab[i]  = eq
                }
                else {
                    g = dwlsmv_ct_tagbefore(nm, tags, rest)
                    f = dwlsmv_ct_pos(fac, rest)
                    if (g == 0 | f == 0) {
                        dwlsmv_ct_bad("a loading without its group tag or factor: " + eq + ":" + nm)
                    }
                    kind[i] = 1
                    rk[i]   = 1
                    sub[i]  = f
                    grp[i]  = g
                    lab[i]  = eq
                }
            }
        }
        else if (anyof(fac, eq) & G == 1) {
            kind[i] = (anyof(fac, nm) ? 2 : 3)
            rk[i]   = (anyof(fac, nm) ? 2 : 3)
            lab[i]  = eq + " <- " + nm
        }
        else {
            dwlsmv_ct_bad("a row that is not in e(b)'s grammar: " + eq + ":" + nm)
        }
    }

    // ---- which rows are shown, which are omitted, and how a fixed row is labelled ----------
    // rkd: 0 omitted, 1 free, 2 marker, 3 constrained
    rkd  = J(1, nb, 1)
    why  = J(nb, 1, "")
    nsh  = 0
    nom  = 0
    for (i = 1; i <= nb; i++) {
        nm = stp[i, 2]
        if (kind[i] == 7 & !showthr) {
            rkd[i] = 0
            why[i] = "thresholds"
        }
        else if (fixed[i]) {
            rkd[i] = 3
            fv = fixval[i]
            if (kind[i] == 5) {
                // a factor variance: var(F) or var(e.F); a residual variance var(e.y) is not one
                inner = substr(lab[i], 5, strlen(lab[i]) - 5)
                isfac = 1
                if (substr(inner, 1, 2) == "e.") {
                    inner = substr(inner, 3, .)
                    isfac = 0
                }
                f = dwlsmv_ct_pos(fac, inner)
            }
            if (std) {
                // C2 (a), owner 2026-11-18, at every G: an exogenous factor variance (and, at G > 1,
                // a scale, intercept or mean) that the standardization makes constant
                if (kind[i] == 5) {
                    if (isfac & f > 0) {
                        rkd[i] = 0
                        why[i] = "standardized"
                    }
                }
                else if (G > 1) {
                    if (kind[i] == 8 | kind[i] == 9 | kind[i] == 10) {
                        rkd[i] = 0
                        why[i] = "standardized"
                    }
                }
            }
            else if (G == 1) {
                if (kind[i] == 5) {
                    if (f > 0 & fv == 1) {
                        if (anyof(uvf, fac[f])) {
                            rkd[i] = 0
                            why[i] = "metric"
                        }
                    }
                }
                else if (kind[i] == 1) {
                    f = sub[i]
                    if (fv == 1) {
                        if (canon != "") {
                            // the marker is the pair e(canon)'s metric marker line names
                            j = 0
                            for (k = 1; k <= cols(mkf); k++) {
                                if (mkf[k] == fac[f] & mky[k] == stp[i, 1]) {
                                    j = 1
                                }
                            }
                            if (j) {
                                rkd[i] = 2
                            }
                        }
                        else if (metric == "marker") {
                            if (!seenmk[f]) {
                                seenmk[f] = 1
                                rkd[i]    = 2
                            }
                        }
                    }
                }
            }
            else {
                if (kind[i] == 5) {
                    if (f > 0 & fv == 1) {
                        rkd[i] = 0
                        why[i] = (grp[i] == 1 ? "metric" : "level")
                    }
                }
                else if (kind[i] == 8 | kind[i] == 10) {
                    if (fv == 0) {
                        rkd[i] = 0
                        why[i] = "level"
                    }
                }
                else if (kind[i] == 9) {
                    if (fv == 1) {
                        rkd[i] = 0
                        why[i] = "level"
                    }
                }
            }
        }
        if (rkd[i] == 0) {
            nom++
        }
        else {
            nsh++
        }
    }

    // ---- the omitted rows ----------------------------------------------------
    otxt = J(nom, 2, "")
    oidx = J(nom, 1, .)
    k = 0
    for (i = 1; i <= nb; i++) {
        if (rkd[i] == 0) {
            k++
            otxt[k, 1] = stp[i, 1] + ":" + stp[i, 2]
            otxt[k, 2] = why[i]
            oidx[k]    = i
        }
    }

    // ---- the shown rows, in display order -------------------------------------
    // keys: group (thresholds after every group), section, factor or item, group within an
    // item's thresholds, the e(b) position (which makes the order total)
    ptxt = J(nsh, 3, "")
    pnum = J(nsh, 11, .)
    if (nsh == 0) {
        return
    }
    K = J(nsh, 5, 0)
    j = 0
    for (i = 1; i <= nb; i++) {
        if (rkd[i] > 0) {
            j++
            K[j, 1] = (kind[i] == 7 ? G + 1 : grp[i])
            K[j, 2] = rk[i]
            K[j, 3] = sub[i]
            K[j, 4] = (kind[i] == 7 ? grp[i] : 0)
            K[j, 5] = i
        }
    }
    o   = order(K, (1, 2, 3, 4, 5))
    blk = 0
    for (j = 1; j <= nsh; j++) {
        bi = K[o[j], 5]
        if (j == 1) {
            blk = 1
        }
        else if (K[o[j], 1] != K[o[j - 1], 1] | K[o[j], 2] != K[o[j - 1], 2] | K[o[j], 3] != K[o[j - 1], 3]) {
            blk++
        }
        k = K[o[j], 2]
        if (k == 1) {
            reason = "Loadings on " + fac[K[o[j], 3]]
        }
        else if (k == 2) {
            reason = "Structural"
        }
        else if (k == 3) {
            reason = "Covariates"
        }
        else if (k == 4) {
            reason = "Variances"
        }
        else if (k == 5) {
            reason = "Covariances"
        }
        else if (k == 6) {
            reason = "Intercepts"
        }
        else if (k == 7) {
            reason = "Scales"
        }
        else if (k == 8) {
            reason = "Factor means"
        }
        else {
            reason = "Thresholds: " + ind[K[o[j], 3]]
        }
        ptxt[j, 1] = reason
        if (G > 1 & K[o[j], 1] <= G) {
            ptxt[j, 2] = "Group " + strofreal(K[o[j], 1]) + " (" + gvar + " = " + gv[K[o[j], 1]] + ")"
        }
        ptxt[j, 3] = lab[bi]
        pnum[j, 1]  = bi
        pnum[j, 2]  = rkd[bi]
        pnum[j, 3]  = tab[1, bi]
        pnum[j, 4]  = tab[2, bi]
        pnum[j, 5]  = tab[3, bi]
        pnum[j, 6]  = tab[4, bi]
        pnum[j, 7]  = tab[5, bi]
        pnum[j, 8]  = tab[6, bi]
        pnum[j, 9]  = blk
        pnum[j, 10] = K[o[j], 1]
        pnum[j, 11] = K[o[j], 2]
    }
}

// a fixed value as the table shows it: 1, 0.8, -0.3, 0.123457
string scalar dwlsmv_ct_val(real scalar x)
{
    string scalar s

    s = strtrim(sprintf("%12.6g", x))
    if (substr(s, 1, 1) == ".") {
        s = "0" + s
    }
    else if (substr(s, 1, 2) == "-.") {
        s = "-0" + substr(s, 2, .)
    }
    return(s)
}

// a number in a column of width wd with d decimals; a large one in g form so the column holds
string scalar dwlsmv_ct_f(real scalar x, real scalar wd, real scalar d)
{
    if (x < . & abs(x) >= 10^(wd - d - 2)) {
        return(sprintf("%" + strofreal(wd) + "." + strofreal(d) + "g", x))
    }
    return(sprintf("%" + strofreal(wd) + "." + strofreal(d) + "f", x))
}

// a label cut to at most nl bytes, with a ~ in the last place when cut
string scalar dwlsmv_ct_cut(string scalar s, real scalar nl)
{
    if (strlen(s) <= nl) {
        return(s)
    }
    return(substr(s, 1, nl - 1) + "~")
}

void dwlsmv_ct_print(string matrix ptxt, real matrix pnum, real scalar hasv, real scalar level,
                     real scalar ls)
{
    real scalar n, r, longest, w, nl, ng, needrule, prevgrp, prevblk, prevrk, body
    string scalar pad, rule, line

    n = rows(ptxt)
    printf("\n")
    if (n == 0) {
        return
    }
    longest = 0
    for (r = 1; r <= n; r++) {
        longest = max((longest, strlen(ptxt[r, 1]), strlen(ptxt[r, 2]), strlen(ptxt[r, 3]) + 2))
    }
    w   = max((16, longest + 2))
    w   = min((w, max((ls - 62, 16))))
    nl  = w - 1
    pad = "%-" + strofreal(nl) + "s"
    body = (hasv ? 61 : 12)
    rule = "{txt}{hline " + strofreal(w) + "}{c +}{hline " + strofreal(body) + "}\n"

    // the column titles
    if (hasv) {
        line = sprintf("%11s%10s%8s%8s", "Estimate", "Std. err.", "z", "P>|z|") +
               sprintf("%24s", "[" + strofreal(level) + "% conf. interval]")
    }
    else {
        line = sprintf("%11s", "Estimate")
    }
    printf("{txt}" + pad + " {c |}%s\n", "", line)
    printf(rule)

    needrule = 0
    prevgrp  = 0
    prevblk  = 0
    prevrk   = 0
    for (r = 1; r <= n; r++) {
        // a new group: a rule, then its title (the first block after it takes no rule of its own)
        if (pnum[r, 10] != prevgrp) {
            if (ptxt[r, 2] != "") {
                if (needrule) {
                    printf(rule)
                }
                printf("{txt}" + pad + " {c |}\n", dwlsmv_ct_cut(ptxt[r, 2], nl))
                needrule = 0
            }
            prevgrp = pnum[r, 10]
        }
        // a new block: a rule, except between two loading blocks or two threshold blocks
        if (pnum[r, 9] != prevblk) {
            if (needrule) {
                if (!(pnum[r, 11] == prevrk & (prevrk == 1 | prevrk == 9) & pnum[r, 10] == prevgrp)) {
                    printf(rule)
                }
            }
            printf("{txt}" + pad + " {c |}\n", dwlsmv_ct_cut(ptxt[r, 1], nl))
            prevblk  = pnum[r, 9]
            prevrk   = pnum[r, 11]
            needrule = 1
        }
        if (pnum[r, 2] == 1) {
            if (hasv) {
                line = dwlsmv_ct_f(pnum[r, 3], 11, 4) + dwlsmv_ct_f(pnum[r, 4], 10, 4) +
                       dwlsmv_ct_f(pnum[r, 5], 8, 2) + dwlsmv_ct_f(pnum[r, 6], 8, 3) +
                       dwlsmv_ct_f(pnum[r, 7], 12, 4) + dwlsmv_ct_f(pnum[r, 8], 12, 4)
            }
            else {
                line = dwlsmv_ct_f(pnum[r, 3], 11, 4)
            }
        }
        else if (pnum[r, 2] == 2) {
            line = sprintf("%11s", dwlsmv_ct_val(pnum[r, 3])) + "  (marker)"
        }
        else {
            line = sprintf("%11s", dwlsmv_ct_val(pnum[r, 3])) + "  (constrained)"
        }
        printf("{txt}" + pad + " {c |}{res}%s\n", dwlsmv_ct_cut("  " + ptxt[r, 3], nl), line)
    }
    // the closing rule, as Stata's own table has (the veto package's mock)
    printf("{txt}{hline " + strofreal(w + 1 + body) + "}\n")
}

void dwlsmv_m11_ctab(string scalar bn, string scalar vn, string scalar cn, real scalar showthr,
                    real scalar std, real scalar level, real scalar ls, string scalar tn)
{
    real rowvector b
    real matrix V, C, tab, pnum, K
    string matrix stp, ptxt, otxt, K_s
    string rowvector ind, gv
    string scalar metric, canon, gvar
    real colvector oidx
    real scalar hasv, hasc, i, nthr, nlev

    b      = st_matrix(bn)
    stp    = st_matrixcolstripe(bn)
    hasv   = (vn != "")
    hasc   = (cn != "")
    V      = J(0, 0, .)
    C      = J(0, 0, .)
    if (hasv) {
        V = st_matrix(vn)
    }
    if (hasc) {
        C = st_matrix(cn)
    }
    metric = st_global("e(metric)")
    canon  = st_global("e(canon)")
    gvar   = st_global("e(groupvar)")
    ind    = tokens(st_global("e(indicators)"))
    gv     = J(1, 0, "")
    if (st_global("e(invariance)") != "") {
        K = st_matrix("e(K)")
        if (rows(K) < 2) {
            dwlsmv_ct_bad("e(K) is missing under invariance()")
        }
        K_s = st_matrixrowstripe("e(K)")
        gv  = K_s[., 2]'
    }
    dwlsmv_ct_plan(stp, b, V, hasv, C, hasc, metric, canon, ind, gv, gvar, showthr, std, level,
                   ptxt, pnum, otxt, oidx, tab)
    dwlsmv_ct_print(ptxt, pnum, hasv, level, ls)

    nthr = 0
    nlev = 0
    for (i = 1; i <= rows(otxt); i++) {
        if (otxt[i, 2] == "thresholds") {
            nthr++
        }
        if (otxt[i, 2] == "metric" | otxt[i, 2] == "level") {
            if (cols(gv) > 1) {
                nlev++
            }
        }
    }
    st_matrix(tn, tab)
    st_matrixcolstripe(tn, stp)
    st_matrixrowstripe(tn, (J(9, 1, ""), ("b" \ "se" \ "z" \ "pvalue" \ "ll" \ "ul" \ "df" \ "crit" \ "eform")))
    st_local("ct_nshow", strofreal(rows(ptxt)))
    st_local("ct_nomit", strofreal(rows(otxt)))
    st_local("ct_nthr", strofreal(nthr))
    st_local("ct_nlev", strofreal(nlev))
}

end
