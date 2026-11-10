*! dwlsmv -- Module 11: translate sem's reading of the model into the canonical form
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_translate.do.
*!
*! Spec Module 11 D1-D3, D9. sem parses (D1): dwlsmv.ado runs
*! sem <paths>, <options> noestimate on a dummy dataset and this file reads
*! what sem left in e() -- e(b)'s stripe, e(Cns), e(lxvars), e(lyvars),
*! e(oyvars) -- into Section 4.8's canonical form. Nothing here parses the
*! user's text.
*!
*!   dwlsmv_translate(final)  canonical lines, one directive per row. 'ind' is
*!                            e(oyvars) in sem's order (D9). final = 0 is D2's
*!                            forcenoanchor pass, which reports and skips the
*!                            D3 scale checks (an unanchored factor is what that
*!                            pass is for); final = 1 applies them.
*! Its companions, one file each: dwlsmv_inject (D2's injection),
*! dwlsmv_semnames (sem's stripe for D6's check; uses dwlsmv_tr_pos from
*! this file) and dwlsmv_benign (D1's dummy-data assertion).
*!
*! What sem accepts and v1a does not is refused here at rc 198 (D8a), with the
*! reason through errprintf(): a constraint tying several coefficients;
*! covariances among indicator residuals (v1b); a fixed indicator residual
*! variance (Theta is derived); a constrained intercept; observed exogenous
*! variables or a path between observed variables; more than one group; and,
*! on the final pass, D3's non-unit scale fixes and factors with no observed
*! indicator. A factor with no constraint at all (sem leaves it unanchored,
*! D3's random-intercept case) is NOT refused here: Module 3's rc 3498 says so.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_tr_pos(string colvector v, string scalar s)
{
    real scalar i
    for (i = 1; i <= rows(v); i++) {
        if (v[i] == s) {
            return(i)
        }
    }
    return(0)
}

// a fixed value as canonical text: the first of the short form, %21.0g and
// %21x that reads back as the identical double. %21.0g alone does not: it
// failed on 22% of random doubles near 1e-3 and 92% at 1e-4 and below
// (measured 2026-10-07, StataMP 18); %21x read back on every value (spec
// D11, precision; fixed here by the user's decision 2026-10-07)
string scalar dwlsmv_tr_num(real scalar x)
{
    string scalar s
    s = strofreal(x)
    if (strtoreal(s) != x) {
        s = strtrim(strofreal(x, "%21.0g"))
    }
    if (strtoreal(s) != x) {
        s = strtrim(strofreal(x, "%21x"))
    }
    return(s)
}

void dwlsmv_tr_refuse(string scalar msg)
{
    errprintf("%s\n", msg)
    exit(198)
}

string colvector dwlsmv_translate(real scalar final)
{
    string matrix st
    string colvector fac, ind, oxv

    st  = st_matrixcolstripe("e(b)")
    fac = tokens(st_global("e(lxvars)") + " " + st_global("e(lyvars)"))'
    ind = tokens(st_global("e(oyvars)"))'
    oxv = tokens(st_global("e(oxvars)"))'

    if (rows(st_numscalar("e(N_groups)")) == 1) {
        if (st_numscalar("e(N_groups)") > 1) {
            dwlsmv_tr_refuse("dwlsmv: more than one group; dwlsmv fits a single group.")
        }
    }
    return(dwlsmv_tr_core(st, st_matrix("e(Cns)"), fac, ind, oxv, final))
}

// the translation proper, from one group's stripe (eq, name) and its Cns rows
// (columns as st, then the right-hand side); dwlsmv_translate_g calls it once
// per group with the group tag removed from the names
string colvector dwlsmv_tr_core(string matrix st, real matrix C, string colvector fac,
    string colvector ind, string colvector oxv, real scalar final)
{
    string colvector out
    real colvector fixed, fixval, isvarfix, haslam, hasmet, fixnon1
    real scalar j, r, nz, cj, f, fa, fz
    string scalar eq, co, a, z, t

    if (rows(oxv) > 0) {
        dwlsmv_tr_refuse("dwlsmv: observed exogenous variables (" + invtokens(oxv') +
            ") are not supported: dwlsmv fits latent factors measured by ordinal indicators.")
    }
    if (rows(fac) == 0) {
        dwlsmv_tr_refuse("dwlsmv: the model has no latent factor.")
    }

    // fixed[j] = 1 if a Cns row constrains column j alone; fixval its value
    fixed  = J(rows(st), 1, 0)
    fixval = J(rows(st), 1, .)
    for (r = 1; r <= rows(C); r++) {
        nz = 0
        for (j = 1; j < cols(C); j++) {
            if (C[r, j] != 0) {
                nz++
                cj = j
            }
        }
        if (nz != 1) {
            dwlsmv_tr_refuse("dwlsmv: a constraint ties several parameters together " +
                "(for example a shared @name). dwlsmv accepts only constraints that fix " +
                "one parameter at a value.")
        }
        fixed[cj]  = 1
        fixval[cj] = C[r, cols(C)] / C[r, cj]
    }

    out      = J(0, 1, "")
    isvarfix = J(rows(fac), 1, 0)
    fixnon1  = J(rows(fac), 1, 0)
    for (j = 1; j <= rows(st); j++) {
        co = st[j, 2]
        for (r = 1; r <= rows(fac); r++) {
            if (co == "var(" + fac[r] + ")" | co == "var(e." + fac[r] + ")") {
                if (fixed[j] & fixval[j] == 1) {
                    isvarfix[r] = 1
                    out = out \ ("metric unitvariance " + fac[r])
                }
                else if (fixed[j]) {
                    fixnon1[r] = 1
                }
            }
        }
    }

    hasmet = isvarfix
    haslam = J(rows(fac), 1, 0)
    for (j = 1; j <= rows(st); j++) {
        eq = st[j, 1]
        co = st[j, 2]
        if (co == "_cons") {
            if (fixed[j]) {
                dwlsmv_tr_refuse("dwlsmv: an intercept (_cons) is constrained; thresholds " +
                    "replace intercepts for ordinal indicators.")
            }
            continue
        }
        if (eq != "/") {
            if (dwlsmv_tr_pos(ind, eq)) {
                r = dwlsmv_tr_pos(fac, co)
                if (r == 0) {
                    dwlsmv_tr_refuse("dwlsmv: a path between observed variables (" + eq +
                        " <- " + co + ") is not supported by dwlsmv.")
                }
                haslam[r] = 1
                if (fixed[j] & fixval[j] == 1 & !hasmet[r]) {
                    out = out \ ("metric marker " + co + " " + eq)
                    hasmet[r] = 1
                }
                else if (fixed[j]) {
                    out = out \ ("lambda " + co + " " + eq + " fixed " + dwlsmv_tr_num(fixval[j]))
                    if (fixval[j] != 1) {
                        fixnon1[r] = 1
                    }
                }
                else {
                    out = out \ ("lambda " + co + " " + eq + " free")
                }
            }
            else {
                if (!dwlsmv_tr_pos(fac, eq) | !dwlsmv_tr_pos(fac, co)) {
                    dwlsmv_tr_refuse("dwlsmv: a path involving an observed variable (" + eq +
                        " <- " + co + ") is not supported by dwlsmv.")
                }
                if (fixed[j]) {
                    out = out \ ("beta " + eq + " " + co + " fixed " + dwlsmv_tr_num(fixval[j]))
                }
                else {
                    out = out \ ("beta " + eq + " " + co + " free")
                }
            }
            continue
        }
        // equation "/": variances and covariances
        if (substr(co, 1, 6) == "var(e.") {
            a = substr(co, 7, strlen(co) - 7)
            if (dwlsmv_tr_pos(ind, a)) {
                if (fixed[j]) {
                    dwlsmv_tr_refuse("dwlsmv: var(e." + a + ") is constrained, but an " +
                        "indicator's residual variance is not a parameter here: it is " +
                        "1 - (Lambda Phi Lambda')_jj, derived from the model.")
                }
                continue
            }
            z = a
        }
        else if (substr(co, 1, 4) == "var(") {
            a = substr(co, 5, strlen(co) - 5)
            z = a
        }
        else if (substr(co, 1, 4) == "cov(") {
            t = substr(co, 5, strlen(co) - 5)
            a = substr(t, 1, strpos(t, ",") - 1)
            z = substr(t, strpos(t, ",") + 1, .)
            if (dwlsmv_tr_pos(ind, subinstr(a, "e.", "")) | dwlsmv_tr_pos(ind, subinstr(z, "e.", ""))) {
                dwlsmv_tr_refuse("dwlsmv: covariances among indicator residuals are not " +
                    "supported (planned for a later version).")
            }
            a = subinstr(a, "e.", "")
            z = subinstr(z, "e.", "")
            // the lower triangle, in factor declaration order
            if (dwlsmv_tr_pos(fac, a) < dwlsmv_tr_pos(fac, z)) {
                t = a
                a = z
                z = t
            }
        }
        else {
            dwlsmv_tr_refuse("dwlsmv: sem reported a parameter dwlsmv cannot read (" +
                co + ").")
        }
        fa = dwlsmv_tr_pos(fac, a)
        fz = dwlsmv_tr_pos(fac, z)
        if (fa == 0 | fz == 0) {
            dwlsmv_tr_refuse("dwlsmv: a variance or covariance of an observed variable (" +
                co + ") is not supported by dwlsmv.")
        }
        if (a == z & fixed[j] & fixval[j] == 1 & isvarfix[fa]) {
            continue
        }
        if (fixed[j]) {
            out = out \ ("psi " + a + " " + z + " fixed " + dwlsmv_tr_num(fixval[j]))
        }
        else {
            out = out \ ("psi " + a + " " + z + " free")
        }
    }

    if (final) {
        for (f = 1; f <= rows(fac); f++) {
            if (!haslam[f]) {
                errprintf("dwlsmv: factor %s has no observed indicator. Second-order factors\n", fac[f])
                errprintf("are not supported: every factor must load on at least one indicator.\n")
                exit(198)
            }
            if (!hasmet[f] & fixnon1[f]) {
                errprintf("dwlsmv: factor %s is scaled by a loading or variance fixed at a value\n", fac[f])
                errprintf("other than 1. A factor's scale is set by a loading fixed at 1 (the marker)\n")
                errprintf("or by its variance fixed at 1 (the disturbance variance, var(e.%s@1), if\n", fac[f])
                errprintf("%s is endogenous); see metric(unitvariance). Other scale values are not\n", fac[f])
                errprintf("supported.\n")
                exit(198)
            }
        }
    }

    return(("factor " + invtokens(fac')) \ ("ind " + invtokens(ind')) \ sort(out, 1))
}

end
