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
*! variance (Theta is derived); a constrained intercept; more than one observed
*! exogenous variable (one binary covariate is read as covariate, gamma and kappa
*! lines, dwlsmv_tr_cxterm); a path between observed variables; more than one
*! group; and,
*! on the final pass, D3's non-unit scale fixes, a factor with no observed
*! indicator and no path to another factor, a nonzero fixed path out of a factor
*! with no indicator, and such a factor with a covariate (release B, B3.1/B3.5;
*! a second-order factor is otherwise accepted, scaled by 'metric
*! unitvariance'). A factor with no constraint at all (sem leaves it unanchored,
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
    string colvector out, head
    real colvector fixed, fixval, isvarfix, haslam, hasmet, fixnon1, sopath, sofix
    real scalar j, r, nz, cj, f, fa, fz, kind
    string scalar eq, co, a, z, t, cvn

    // DECIDED P1 (owner, 2026-11-14, docs/decisions.md "Stage 2 decisions package"): exactly one observed exogenous variable is
    // accepted, as the covariate of a conditional fit (one binary 0/1 x, one group; spec 4.8,
    // Conditional covariates). Two or more are refused here, explicitly and with their own message (not as a side
    // effect of another check). dwlsmv_translate_g refuses any (a covariate with group() is not in the owner's package).
    cvn = ""
    if (rows(oxv) > 1) {
        dwlsmv_tr_refuse("dwlsmv: more than one observed exogenous variable (" + invtokens(oxv') +
            "). A model can have one binary 0/1 covariate in this version; the others are not available.")
    }
    if (rows(oxv) == 1) {
        cvn = oxv[1]
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
    // release B (B3.1): sopath[f] = a beta row out of f that is free or fixed nonzero;
    // sofix[f] = such a row fixed at a nonzero value (a beta fixed at 0 is no path)
    sopath   = J(rows(fac), 1, 0)
    sofix    = J(rows(fac), 1, 0)
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
                    // (y <- x): the covariate's direct effect on the indicator's latent response.
                    // P3 (@c semantics), approved by the owner on 2026-11-14: a value fixed with @ is written literally (fixed 0 too).
                    if (co == cvn) {
                        out = out \ ("kappa " + eq + " " + co +
                            (fixed[j] ? " fixed " + dwlsmv_tr_num(fixval[j]) : " free"))
                        continue
                    }
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
                // (F <- x): the covariate's effect on the factor. P3 as for kappa: semantics approved by the owner on 2026-11-14.
                if (co == cvn) {
                    if (dwlsmv_tr_pos(fac, eq)) {
                        out = out \ ("gamma " + eq + " " + co +
                            (fixed[j] ? " fixed " + dwlsmv_tr_num(fixval[j]) : " free"))
                        continue
                    }
                }
                if (!dwlsmv_tr_pos(fac, eq) | !dwlsmv_tr_pos(fac, co)) {
                    dwlsmv_tr_refuse("dwlsmv: a path involving an observed variable (" + eq +
                        " <- " + co + ") is not supported by dwlsmv.")
                }
                r = dwlsmv_tr_pos(fac, co)
                if (fixed[j]) {
                    out = out \ ("beta " + eq + " " + co + " fixed " + dwlsmv_tr_num(fixval[j]))
                    if (fixval[j] != 0) {
                        sopath[r] = 1
                        sofix[r]  = 1
                    }
                }
                else {
                    out = out \ ("beta " + eq + " " + co + " free")
                    sopath[r] = 1
                }
            }
            continue
        }
        // equation "/": variances and covariances
        kind = dwlsmv_tr_cxterm(co, cvn)
        if (kind == 1 | kind == 2) {
            // var(x) and mean(x) listed unconstrained: x is data in a conditional fit, so they
            // are skipped. A constrained one is the user asking for a parameter that is not one.
            if (fixed[j]) {
                dwlsmv_tr_refuse("dwlsmv: the variance or mean of the covariate " + cvn + " is constrained (" +
                    co + "). In a conditional fit x is data; its variance and mean are not parameters.")
            }
            continue
        }
        if (kind == 3) {
            // DECIDED P2 (owner, 2026-11-14: the refusal is retained): sem lets an observed exogenous x covary
            // with every exogenous factor by default. In dwlsmv's conditional model that
            // covariance is not a parameter: the factor either depends on x (gamma) or does not.
            // Only a covariance the user fixed at zero is accepted (and then gives no gamma row);
            // a free one is refused with the two ways out. A covariance directive is not reinterpreted as a regression path.
            if (fixed[j]) {
                if (fixval[j] == 0) {
                    continue
                }
            }
            t = substr(co, 5, strlen(co) - 5)
            if (substr(t, 1, strlen(cvn) + 1) == cvn + ",") {
                a = substr(t, strlen(cvn) + 2, .)
            }
            else {
                a = substr(t, 1, strlen(t) - strlen(cvn) - 1)
            }
            if (fixed[j]) {
                errprintf("dwlsmv: the covariance of the covariate %s with %s is fixed at a value other than 0.\n", cvn, a)
                errprintf("Only a covariance fixed at zero is accepted: cov(%s*%s@0).\n", cvn, a)
                exit(198)
            }
            errprintf("dwlsmv: sem lets the covariate %s correlate with %s (%s), and in a conditional\n", cvn, a, co)
            errprintf("fit that covariance is not a parameter. To let %s depend on %s, write the path (%s <- %s).\n", a, cvn, a, cvn)
            errprintf("To say it does not, fix the covariance at zero with the sem option cov(%s*%s@0).\n", cvn, a)
            exit(198)
        }
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
                // release B (B3.1, B3.5): a factor with no indicator and a path out of it is a
                // second-order factor, scaled by its variance ('metric unitvariance')
                if (!sopath[f]) {
                    errprintf("dwlsmv: factor %s has no observed indicator and no path to another factor.\n", fac[f])
                    errprintf("A factor needs an observed indicator, or paths to the factors it is\n")
                    errprintf("measured by (a second-order factor).\n")
                    exit(198)
                }
                if (sofix[f]) {
                    errprintf("dwlsmv: factor %s has no observed indicator, and a path out of it is fixed at a nonzero value.\n", fac[f])
                    errprintf("A factor with no indicator is scaled by its variance; fixing a path as well adds a second\n")
                    errprintf("restriction. Free the path and write var(%s@1), or leave the scale to dwlsmv.\n", fac[f])
                    exit(198)
                }
                if (cvn != "") {
                    errprintf("dwlsmv: a second-order factor (%s) is not available with a covariate (%s) in this version.\n", fac[f], cvn)
                    errprintf("No model was fitted.\n")
                    exit(198)
                }
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

    head = ("factor " + invtokens(fac')) \ ("ind " + invtokens(ind'))
    if (cvn != "") {
        head = head \ ("covariate " + cvn)
    }
    return(head \ sort(out, 1))
}

// release B: the factors of canonical text P (lines, as dwlsmv_tr_core returns them) that have
// no observed indicator (no lambda row, no 'metric marker' row) and a path out of them: a row
// 'beta child F' that is free or fixed nonzero. A beta fixed at 0 is no path (spec 4.8).
string rowvector dwlsmv_tr_sofac(string colvector P)
{
    string rowvector fac, out, t
    real scalar k, i, hasind, outp

    fac = tokens(substr(P[1], 8, .))
    out = J(1, 0, "")
    for (k = 1; k <= cols(fac); k++) {
        hasind = 0
        outp   = 0
        for (i = 3; i <= rows(P); i++) {
            t = tokens(P[i])
            if (t[1] == "lambda") {
                if (t[2] == fac[k]) {
                    hasind = 1
                }
            }
            else if (t[1] == "metric") {
                if (t[2] == "marker") {
                    if (t[3] == fac[k]) {
                        hasind = 1
                    }
                }
            }
            else if (t[1] == "beta") {
                if (t[3] == fac[k]) {
                    if (t[4] == "free") {
                        outp = 1
                    }
                    else if (strtoreal(t[5]) != 0) {
                        outp = 1
                    }
                }
            }
        }
        if (!hasind) {
            if (outp) {
                out = out, fac[k]
            }
        }
    }
    return(out)
}

// the same, for the ";"-joined text _dwlsmv_translate holds in a local
string rowvector dwlsmv_tr_sofac_t(string scalar canon)
{
    string colvector P
    string scalar s
    real scalar p

    P = J(0, 1, "")
    s = canon
    while (strlen(s) > 0) {
        p = strpos(s, ";")
        if (p == 0) {
            P = P \ s
            s = ""
        }
        else {
            P = P \ substr(s, 1, p - 1)
            s = substr(s, p + 1, .)
        }
    }
    return(dwlsmv_tr_sofac(P))
}

end
