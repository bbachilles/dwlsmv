*! dwlsmv -- Module 11 D2 probe: metric(unitvariance) by injection
*!
*! A probe, not a test (tools/). Module 11's D2 (2026-10-04, the user's): the
*! default metric is sem's marker (metric(marker)); metric(unitvariance)
*! switches every factor the USER did not normalise to unit variance. Which
*! factors those are comes from sem itself: under forcenoanchor, e(Cns) holds
*! only the user's own constraints. Module 11 then adds var(F@1), or
*! var(e.F@1) for an endogenous F (one with a beta line), and runs sem again.
*! Reuses tools/probe_d1_sem.do's prototype (d1_T, d1_semdummy), which is
*! scratch, not the Module 11 implementation.
*!
*! Checks: forcenoanchor's e(Cns) lists exactly the user's constraints; the
*! injected models' canonical text. Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_d2_metric.do    (log: probe_d2_metric.log)
*!     grep -E "^D2(ENV|NA|TXT|INJ|RC)|^r\(" probe_d2_metric.log | diff - tools/probe_d2_metric_stata18.txt
*! Only the D2ENV line is expected to differ.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 200
display "D2ENV stata_version " c(stata_version) " processors " c(processors)
mata:
// ---------------- route T -------------------------------------------------
string colvector d1_T()
{
    real matrix b, C
    string matrix st, cs
    string colvector out, fac, ind, endo
    real scalar j, r, nz, cj, v
    string scalar eq, co, f, y, a, z
    real colvector fixed
    real rowvector fixval, isvarfix

    b  = st_matrix("e(b)")
    st = st_matrixcolstripe("e(b)")
    C  = st_matrix("e(Cns)")
    fac  = tokens(st_global("e(lxvars)") + " " + st_global("e(lyvars)"))'
    endo = tokens(st_global("e(lyvars)"))'
    ind  = tokens(st_global("e(oyvars)"))'
    // fixed[j] = 1 if a Cns row constrains column j alone; fixval its value
    fixed  = J(cols(b), 1, 0)
    fixval = J(1, cols(b), .)
    for (r = 1; r <= rows(C); r++) {
        nz = 0
        for (j = 1; j <= cols(b); j++) {
            if (C[r, j] != 0) {
                nz++
                cj = j
            }
        }
        if (nz != 1) _error("d1_T: a multi-coefficient constraint")
        fixed[cj] = 1
        fixval[cj] = C[r, cols(C)] / C[r, cj]
    }
    out = J(0, 1, "")
    // which factors have their variance fixed to 1 -> unitvariance
    isvarfix = J(1, rows(fac), 0)
    for (j = 1; j <= cols(b); j++) {
        co = st[j, 2]
        for (r = 1; r <= rows(fac); r++) {
            if ((co == "var(" + fac[r] + ")" | co == "var(e." + fac[r] + ")") & fixed[j] & fixval[j] == 1) {
                isvarfix[r] = 1
                out = out \ ("metric unitvariance " + fac[r])
            }
        }
    }
    for (j = 1; j <= cols(b); j++) {
        eq = st[j, 1]
        co = st[j, 2]
        if (co == "_cons") continue
        if (eq != "/") {
            if (anyof(ind, eq)) {
                // loading lambda[co, eq]
                r = d1_pos(fac, co)
                if (fixed[j] & fixval[j] == 1 & !isvarfix[r] & !d1_hasmarker(out, co)) {
                    out = out \ ("metric marker " + co + " " + eq)
                }
                else if (fixed[j]) out = out \ ("lambda " + co + " " + eq + " fixed " + strofreal(fixval[j]))
                else out = out \ ("lambda " + co + " " + eq + " free")
            }
            else {
                if (fixed[j]) out = out \ ("beta " + eq + " " + co + " fixed " + strofreal(fixval[j]))
                else out = out \ ("beta " + eq + " " + co + " free")
            }
            continue
        }
        // eq "/": variances and covariances
        if (substr(co, 1, 6) == "var(e.") {
            f = substr(co, 7, strlen(co) - 7)
            if (anyof(ind, f)) continue            // Theta: derived, never free
            a = f
            z = f
        }
        else if (substr(co, 1, 4) == "var(") {
            a = substr(co, 5, strlen(co) - 5)
            z = a
        }
        else if (substr(co, 1, 4) == "cov(") {
            f = substr(co, 5, strlen(co) - 5)
            a = substr(f, 1, strpos(f, ",") - 1)
            z = substr(f, strpos(f, ",") + 1, .)
            if (substr(a, 1, 2) == "e." | substr(z, 1, 2) == "e.") {
                if (anyof(ind, subinstr(a, "e.", "")) | anyof(ind, subinstr(z, "e.", ""))) {
                    _error("d1_T: residual covariance among indicators (v1b)")
                }
                a = subinstr(a, "e.", "")
                z = subinstr(z, "e.", "")
            }
            // lower triangle in declaration order
            if (d1_pos(fac, a) < d1_pos(fac, z)) {
                f = a
                a = z
                z = f
            }
        }
        if (a == z & fixed[j] & fixval[j] == 1 & isvarfix[d1_pos(fac, a)]) continue
        if (fixed[j]) out = out \ ("psi " + a + " " + z + " fixed " + strofreal(fixval[j]))
        else out = out \ ("psi " + a + " " + z + " free")
    }
    out = ("factor " + invtokens(fac')) \ ("ind " + invtokens(ind')) \ sort(out, 1)
    return(out)
}

real scalar d1_pos(string colvector v, string scalar s)
{
    real scalar i
    for (i = 1; i <= rows(v); i++) if (v[i] == s) return(i)
    _error("d1_pos: " + s + " not found")
}

real scalar d1_hasmarker(string colvector out, string scalar f)
{
    real scalar i
    for (i = 1; i <= rows(out); i++) {
        if (substr(out[i], 1, 14 + strlen(f)) == "metric marker " + f) return(1)
    }
    return(0)
}

// ---------------- route H -------------------------------------------------
// paths: "(F -> y1 y2@.5 y3)" measurement, "(F2 <- F1 [F3])" structural.
// opts: var(F@1) var(e.F@1) cov(A*B) cov(A*B@v). No other syntax.
string colvector d1_H(string scalar paths, string scalar opts)
{
    string rowvector grp, tk, ot
    string colvector out, fac, ind, endo, exo
    string matrix lam, cv
    real scalar g, i, j, at, nf
    string scalar s, lhs, rhs, it, v, a, z, t
    real rowvector uv, hasfixl

    fac = J(0, 1, "")
    ind = J(0, 1, "")
    endo = J(0, 1, "")
    lam = J(0, 3, "")             // F, y, value or ""
    out = J(0, 1, "")
    s = paths
    while ((i = strpos(s, "(")) > 0) {
        j = strpos(s, ")")
        it = strtrim(substr(s, i + 1, j - i - 1))
        s = substr(s, j + 1, .)
        if (strpos(it, "->")) {
            lhs = strtrim(substr(it, 1, strpos(it, "->") - 1))
            rhs = substr(it, strpos(it, "->") + 2, .)
            if (!anyof(fac, lhs)) fac = fac \ lhs
            tk = tokens(rhs)
            for (g = 1; g <= cols(tk); g++) {
                at = strpos(tk[g], "@")
                if (at) lam = lam \ (lhs, substr(tk[g], 1, at - 1), substr(tk[g], at + 1, .))
                else    lam = lam \ (lhs, tk[g], "")
                if (!anyof(ind, lam[rows(lam), 2])) ind = ind \ lam[rows(lam), 2]
            }
        }
        else if (strpos(it, "<-")) {
            lhs = strtrim(substr(it, 1, strpos(it, "<-") - 1))
            tk = tokens(substr(it, strpos(it, "<-") + 2, .))
            if (!anyof(endo, lhs)) endo = endo \ lhs
            for (g = 1; g <= cols(tk); g++) out = out \ ("beta " + lhs + " " + tk[g] + " free")
        }
        else _error("d1_H: unsupported path")
    }
    nf = rows(fac)
    // options
    uv = J(1, nf, 0)
    cv = J(0, 3, "")
    ot = tokens(opts)
    for (g = 1; g <= cols(ot); g++) {
        t = ot[g]
        if (substr(t, 1, 4) == "var(") {
            it = substr(t, 5, strlen(t) - 5)
            a = substr(it, 1, strpos(it, "@") - 1)
            v = substr(it, strpos(it, "@") + 1, .)
            if (substr(a, 1, 2) == "e.") a = substr(a, 3, .)
            if (anyof(endo, a) & substr(substr(t, 5, .), 1, 2) != "e.") _error("d1_H: var(F) of an endogenous factor is its total variance")
            if (v == "1") uv[d1_pos(fac, a)] = 1
            else out = out \ ("psi " + a + " " + a + " fixed " + v)
        }
        else if (substr(t, 1, 4) == "cov(") {
            it = substr(t, 5, strlen(t) - 5)
            v = ""
            if (strpos(it, "@")) {
                v = substr(it, strpos(it, "@") + 1, .)
                it = substr(it, 1, strpos(it, "@") - 1)
            }
            a = substr(it, 1, strpos(it, "*") - 1)
            z = substr(it, strpos(it, "*") + 1, .)
            if (substr(a, 1, 2) == "e." | substr(z, 1, 2) == "e.") _error("d1_H: covariance with an error term")
            cv = cv \ (a, z, v)
        }
        else _error("d1_H: unsupported option " + t)
    }
    // metric: unitvariance where var fixed at 1; else marker = first-listed
    // indicator, unless one of F's loadings is already fixed (sem's rule)
    hasfixl = J(1, nf, 0)
    for (i = 1; i <= rows(lam); i++) if (lam[i, 3] != "") hasfixl[d1_pos(fac, lam[i, 1])] = 1
    for (g = 1; g <= nf; g++) {
        if (uv[g]) out = out \ ("metric unitvariance " + fac[g])
    }
    for (i = 1; i <= rows(lam); i++) {
        g = d1_pos(fac, lam[i, 1])
        if (lam[i, 3] == "1" & !uv[g] & !d1_hasmarker(out, fac[g])) {
            out = out \ ("metric marker " + fac[g] + " " + lam[i, 2])
        }
        else if (lam[i, 3] != "") out = out \ ("lambda " + fac[g] + " " + lam[i, 2] + " fixed " + lam[i, 3])
        else if (!uv[g] & !hasfixl[g] & !d1_hasmarker(out, fac[g])) {
            out = out \ ("metric marker " + fac[g] + " " + lam[i, 2])
        }
        else out = out \ ("lambda " + fac[g] + " " + lam[i, 2] + " free")
    }
    // variances: every factor's psi free unless unitvariance or fixed above
    for (g = 1; g <= nf; g++) {
        if (!uv[g] & !d1_haspsi(out, fac[g], fac[g])) out = out \ ("psi " + fac[g] + " " + fac[g] + " free")
    }
    // covariances: exogenous pairs free by default (sem's default), cov() overrides
    exo = J(0, 1, "")
    for (g = 1; g <= nf; g++) if (!anyof(endo, fac[g])) exo = exo \ fac[g]
    for (i = 2; i <= rows(exo); i++) {
        for (j = 1; j < i; j++) {
            a = exo[i]
            z = exo[j]
            v = "free"
            for (g = 1; g <= rows(cv); g++) {
                if ((cv[g, 1] == a & cv[g, 2] == z) | (cv[g, 1] == z & cv[g, 2] == a)) {
                    v = (cv[g, 3] == "" ? "free" : "fixed " + cv[g, 3])
                }
            }
            if (d1_pos(fac, a) < d1_pos(fac, z)) {
                t = a
                a = z
                z = t
            }
            if (v != "fixed 0") out = out \ ("psi " + a + " " + z + " " + v)
        }
    }
    return(("factor " + invtokens(fac')) \ ("ind " + invtokens(ind')) \ sort(out, 1))
}

real scalar d1_haspsi(string colvector out, string scalar a, string scalar z)
{
    real scalar i
    for (i = 1; i <= rows(out); i++) if (substr(out[i], 1, 5 + strlen(a) + strlen(z)) == "psi " + a + " " + z) return(1)
    return(0)
}

void d1_cmp(string colvector T, string colvector H)
{
    real scalar i
    if (rows(T) == rows(H)) {
        if (all(T :== H)) {
            printf("D1CMP identical, %f lines\n", rows(T))
            for (i = 1; i <= rows(T); i++) printf("   %s\n", T[i])
            return
        }
    }
    printf("D1CMP DIFFER\n")
    printf("  T:\n")
    for (i = 1; i <= rows(T); i++) printf("   %s\n", T[i])
    printf("  H:\n")
    for (i = 1; i <= rows(H); i++) printf("   %s\n", H[i])
}
end
program d1_semdummy, rclass
    args paths opts
    * the model's observed variables: tokens of the text that unab in the data
    local txt `"`paths' `opts'"'
    foreach ch in "(" ")" "<" ">" "-" "@" "*" "," "=" {
        if "`ch'" != "-" local txt : subinstr local txt "`ch'" " ", all
    }
    * "-" separates arrows from varlist ranges; handle arrows only
    local txt : subinstr local txt "->" " ", all
    local txt : subinstr local txt "<-" " ", all
    local txt : subinstr local txt "(" " ", all
    local txt : subinstr local txt ")" " ", all
    local obs
    foreach t of local txt {
        capture unab u : `t'
        if !_rc local obs `obs' `u'
    }
    local obs : list uniq obs
    local p : word count `obs'
    if `p' == 0 | `p' > 100 {
        display as error "d1_semdummy: `p' observed variables (1 to 100 supported)"
        exit 198
    }
    local rng `c(rngstate)'
    tempname est
    capture _estimates hold `est', restore nullok
    preserve
    capture noisily {
        drop _all
        quietly set obs 200
        set seed 20261003
        tempvar f
        quietly generate double `f' = rnormal()
        foreach v of local obs {
            quietly generate double `v' = 0.6 * `f' + 0.8 * rnormal()
        }
        drop `f'
        mata: d1_benign(tokens(st_local("obs")))
        quietly sem `paths', `opts' noestimate
        mata: D1OUT = d1_T()
    }
    local rc = _rc
    restore
    set rngstate `rng'
    if `rc' exit `rc'
end

mata:
// asserted, not assumed: no missing, no zero variance, |r| < 0.9, and the
// correlation matrix well conditioned
void d1_benign(string rowvector v)
{
    real matrix X, R
    real colvector ev
    X = st_data(., v)
    if (missing(X)) _error("d1_benign: missing values")
    R = correlation(X)
    if (any(diagonal(variance(X)) :<= 0)) _error("d1_benign: zero variance")
    if (cols(X) > 1) {
        if (max(abs(R - I(cols(X)))) >= 0.9) _error("d1_benign: a correlation at 0.9 or more")
    }
    ev = symeigenvalues(R)'
    if (min(ev) / max(ev) < 1e-3) _error("d1_benign: ill-conditioned correlation matrix")
}
end
import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
destring y*, replace

mata:
// factors the user normalised, from the forcenoanchor pass's canonical text:
// a metric line, a fixed loading, or a fixed variance on F
real scalar d2_usernorm(string colvector P, string scalar f)
{
    real scalar i
    string rowvector t
    for (i = 3; i <= rows(P); i++) {
        // nested, not &: Mata's & does not short-circuit, and a metric
        // line has three tokens where lambda and psi lines have four
        t = tokens(P[i])
        if (t[1] == "metric") {
            if (t[3] == f) return(1)
        }
        else if (t[1] == "lambda" | t[1] == "psi") {
            if (t[2] == f & t[4] == "fixed") {
                if (t[1] == "lambda" | t[3] == f) return(1)
            }
        }
    }
    return(0)
}

string scalar d2_inject(string colvector P)
{
    string rowvector fac, endo
    string scalar inj
    real scalar i
    fac  = tokens(substr(P[1], 8, .))
    endo = J(1, 0, "")
    for (i = 3; i <= rows(P); i++) {
        if (substr(P[i], 1, 5) == "beta ") endo = endo, tokens(P[i])[2]
    }
    inj = ""
    for (i = 1; i <= cols(fac); i++) {
        if (d2_usernorm(P, fac[i])) continue
        inj = inj + (anyof(endo, fac[i]) ? " var(e." : " var(") + fac[i] + "@1)"
    }
    return(strtrim(inj))
}
end

program d2uv
    args paths opts
    display as text "==== `paths' , `opts'"
    * the user's own constraints only
    quietly sem `paths', `opts' forcenoanchor noestimate
    local ncns 0
    capture confirm matrix e(Cns)
    if !_rc {
        tempname C
        matrix `C' = e(Cns)
        local ncns = rowsof(`C')
    }
    display as result "D2NA `ncns' user constraint(s)"
    capture noisily d1_semdummy "`paths'" "`opts' forcenoanchor"
    if _rc {
        display as result "D2RC forcenoanchor pass " _rc
        exit
    }
    mata: st_local("inj", d2_inject(D1OUT))
    display as result "D2INJ `inj'"
    capture noisily d1_semdummy "`paths'" "`opts' `inj'"
    if _rc {
        display as result "D2RC injected pass " _rc
        exit
    }
    mata: for (i = 1; i <= rows(D1OUT); i++) printf("D2TXT %s\n", D1OUT[i])
end

d2uv "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" ""
d2uv "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)" ""
d2uv "(F1 -> y1@1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" ""
d2uv "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)" "var(F1@1)"
d2uv "(F1 -> y1 y2@.5 y3) (F2 -> y4 y5 y6)" ""
d2uv "(y1 y2 y3 <- F1) (y4 y5 y6 <- F2)" ""
d2uv "(F1 -> y1-y3) (F2 -> y4-y6)" "cov(F1*F2@.3)"
d2uv "(F1 -> y1@1 y2 y3) (F2 -> y4 y5 y6)" "var(F1@1)"
