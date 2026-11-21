*! Executed under StataMP 18 on 2026-11-20: all checks pass on the first run. Every fit maps
*! all its rows (69, 52, 74, 52, 105); the largest gap is 0.29 of the SE tier (J_V), V against U
*! agree to 2e-12 (estimates) and 6e-11 (covariance), no lavaan fit needed reflecting.
*!
*! Release B, unit B5: the standardized solution of a model with a second-order factor
*! (derivation @sec-so-standardize SO5; SO3.2: the standardized vector is the same in every
*! metric), dwlsmv_m11_std against fixtures/expected/lavaan_so_standardized.csv (tier 2).
*!
*! Builds build/ldwlsmv.mlib first and runs against it and src/ado, as test_standardized.do.
*!
*! Tolerances, as test_standardized.do's (spec 7.2) and set before any run:
*!   values 1e-5 abs + 1e-4 rel; SEs 1e-8 abs + 1e-5 rel (lavaan's forward-difference SE error
*!   is ~1.5e-8 relative); a row lavaan gives SE 0 must have SE exactly 0 here, and only those.
*!   The lavaan references stop at their own optimizer floor (about 2e-6 on Case J, E-114);
*!   that is inside the value tier and, on SEs (1.3e-7 absolute against 4e-2), inside the SE tier.
*!   Nothing is loosened if a comparison fails.
*!   Check 3: 1e-5 (the B4 cross-metric tolerance, chosen after measuring, E-115).
*!   Check 4: 1e-12 for the algebraic identity psi*_k = 1 - b*_k^2.
*!
*! Where the standardized results live: nothing is stored (D8d). Every check calls
*! dwlsmv_m11_std("<b>","<V>","<C>") on the current e(), as test_standardized.do does; the
*! posted matrices are tg_<fit> (1 x m), tV_<fit> (m x m), tC_<fit>. Each fit is also run with
*! the standardized option (rc 0 asserted) so the display path is exercised on these models.
*!
*! Fit labels (fixtures/generate/99s_lavaan_second_order.R):
*!   H4_V H4_U  case_H  y1-y12  G over F1-F4            metric(marker) / metric(unitvariance)
*!   H3_V H3_U  case_H  y1-y9   G over F1-F3
*!   HD_M HD_U  case_H  y1-y13  G -> y13 F1-F4; the item is listed first so it anchors G
*!                              (lavaan HD_M: y13 is G's marker). HD_M is dwlsmv metric(marker),
*!                              HD_U metric(unitvariance).
*!   I_V        case_I  y1-y9   G over F1-F3, improper (psi_F1 < 0)  metric(marker)
*!   J_V  J_U   case_J  y1-y18  G1 over F1-F3, G2 over F4-F6; G1, G2 covary by sem's default
*! Not compared: the *_cor fits (correlated factors, not second-order), I_U (did not converge).
*!
*! Row map, lavaan row -> dwlsmv e(b) name (every row must map to exactly one and every dwlsmv
*! row must be used; lavaan has no '~' rows in these fits, a second-order path is G =~ F):
*!   P =~ C       -> C:P                (item loading and second-order path alike, a stripe y:F or F:G)
*!   y | tK       -> /y:cutK
*!   A ~~ B       -> /cov(A,B) or /cov(B,A)
*!   y ~~ y       -> /var(e.y)
*!   F ~~ F       -> /var(e.F) (endogenous; the disturbance variance psi*_k) else /var(G) (exogenous, = 1)
*!
*! Checks:
*!   1  EVERY STANDARDIZED ESTIMATE AND SE vs lavaan, nine fits, at the tiers above. Kills S1
*!      (path standardization written s[to]/s[from]: every F:G value moves), S4 (ds_t term of the
*!      path Jacobian dropped: path SEs move), S2 (psi_ff of an exogenous factor not reset to 1:
*!      var(G) rows lose the exact SE 0 that lavaan has; whether the unreset row is exactly 0 by
*!      rounding is not known, so S2 may be killed by check 1 only if the residue is nonzero).
*!   2  SIGN (B3.10). Orientation of lavaan's rows: a first-order factor is reflected only if its
*!      first item loading is negative, a factor without an item (G) only if the sum over its children
*!      of (child's reflection * std path) is negative, a factor with an item (MD's G) by its first item;
*!      reflected rows are loading/path rows * r_parent * r_child and covariances * r_a * r_b. The fits
*!      reflected are printed (by the CSV none is, so this is a guard on the fixture). On dwlsmv's side,
*!      in every fit each factor's first item loading is positive and each G's sum of standardized
*!      paths is positive, without reflecting anything. Kills S3 (pass 2 of dwlsmv_signcanon removed:
*!      H4_U's paths came out negative, E-116), through checks 1 and 2.
*!   3  METRIC INVARIANCE (SO3.2, SO5.4): dwlsmv's standardized vector of the V fit equals the U fit's
*!      (mreldif, 1e-5) and its covariance the U fit's (max abs gap over max |V_U|, 1e-5), by name,
*!      H4, H3, J. HD (information only, no assertion). Kills S1 only if it differs by metric; it is
*!      the one check that sees a standardization bug that depends on psi_G, and a dwlsmv_signcanon bug
*!      that leaves G reflected in one metric.
*!   4  IMPROPER I_V: standardized path of F1 on G above 1, standardized disturbance variance of F1
*!      negative, psi*_1 = 1 - b*_1^2 (1e-12), and both equal lavaan's (read from the fixture) at the
*!      tiers. SO5's psi_k < 0, phi_kk > 0 case.
*!   Printed, no assertion: the smallest SE among each fit's standardized second-order paths.
*!
*! Mutants: S1 "b * s[from] / s[to]" inverted, S2 the last loop of dwlsmv_m11_std, S3 pass 2 of
*! dwlsmv_signcanon removed, S4 the ds_t term dropped from the path row of Jg.
*! Run 2026-11-20, each in a scratch worktree, all four killed (exit 9):
*!   S1 check 1, value gap 4.7e+03 of tier (H4_V first).
*!   S2 check 1 only through a rounding residue: HD_M's /var(G) SE is 1.6e-17, not 0. Every other
*!      exogenous psi_ff row came out exactly 0 unreset, so this kill depends on rounding: killed,
*!      but not a guarantee.
*!   S3 checks 1 and 2 (H4_U, H3_V: orientation violation, value gap 1.8e+04 of tier).
*!   S4 check 1, SE gap 3.1e+04 of tier (H4_V first).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 200

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

* ---- the fixture, into Mata ---------------------------------------------------
import delimited using "fixtures/expected/lavaan_so_standardized.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: tst_fx = st_sdata(., ("fit", "lhs", "op", "rhs", "est_std", "se"))
mata: tst_nfail = 0

mata:
// the e(b)-style full name of a stripe row
string colvector tst_names(string matrix stp)
{
    string colvector nm
    real scalar i
    nm = J(rows(stp), 1, "")
    for (i = 1; i <= rows(stp); i++) {
        nm[i] = (stp[i, 1] == "/" ? "/" + stp[i, 2] : stp[i, 1] + ":" + stp[i, 2])
    }
    return(nm)
}

// candidate dwlsmv names of one lavaan row (fit, lhs, op, rhs, est, se)
string rowvector tst_cand(string rowvector r)
{
    string scalar l, o, h
    l = r[2]
    o = r[3]
    h = r[4]
    if (o == "=~") {
        return((h + ":" + l))
    }
    if (o == "|") {
        return(("/" + l + ":cut" + substr(h, 2, .)))
    }
    if (l != h) {
        return(("/cov(" + l + "," + h + ")", "/cov(" + h + "," + l + ")"))
    }
    if (substr(l, 1, 1) == "y") {
        return(("/var(e." + l + ")"))
    }
    return(("/var(e." + l + ")", "/var(" + l + ")"))
}

// the reflection of a factor (1 if it is not in the list)
real scalar tst_r(string scalar f, string colvector facs, real colvector r)
{
    real colvector k
    k = selectindex(facs :== f)
    if (length(k) == 0) {
        return(1)
    }
    return(r[k[1]])
}

// check 1, 2 for one fit: dwlsmv_m11_std on the current e() against the fixture
void tst_vs_lavaan(string scalar fitn)
{
    external string matrix tst_fx
    external real scalar tst_nfail
    string matrix stp, fx
    string colvector nm, lhs, op, rhs, facs
    string rowvector cand
    real rowvector g
    real matrix V
    real colvector est, used, r, kk, idx, hit, gv, sg, cnt
    real scalar i, j, k, n, a, b, se, lse, re, rs, bad, nref, s, f, minse, bado, nm1, sse
    string scalar refl

    dwlsmv_m11_std("tg_" + fitn, "tV_" + fitn, "tC_" + fitn)
    g   = st_matrix("tg_" + fitn)
    V   = st_matrix("tV_" + fitn)
    stp = st_matrixcolstripe("tg_" + fitn)
    nm  = tst_names(stp)
    fx  = select(tst_fx, tst_fx[., 1] :== fitn)
    n   = rows(fx)
    lhs = fx[., 2]
    op  = fx[., 3]
    rhs = fx[., 4]
    est = strtoreal(fx[., 5])

    // map every lavaan row to exactly one dwlsmv row
    kk   = J(n, 1, .)
    gv   = J(n, 1, .)
    used = J(rows(nm), 1, 0)
    bad  = 0
    for (i = 1; i <= n; i++) {
        cand = tst_cand(fx[i, .])
        hit  = J(0, 1, .)
        for (j = 1; j <= cols(cand); j++) {
            idx = selectindex(nm :== cand[j])
            hit = hit \ idx
        }
        if (rows(hit) != 1) {
            printf("{err}    %s: lavaan row %s %s %s maps to %f dwlsmv rows\n", fitn, lhs[i], op[i], rhs[i], rows(hit))
            bad++
            continue
        }
        kk[i]      = hit[1]
        gv[i]      = g[hit[1]]
        used[hit[1]] = used[hit[1]] + 1
    }

    // lavaan's orientation (B3.10): pass 1 a factor with an item by its first item, pass 2 a factor
    // without one by the sum of (child's reflection) * (its std path)
    facs = uniqrows(select(lhs, op :== "=~"))
    r    = J(rows(facs), 1, 1)
    for (i = 1; i <= rows(facs); i++) {
        idx = selectindex((lhs :== facs[i]) :& (op :== "=~") :& (substr(rhs, 1, 1) :== "y"))
        if (length(idx) > 0) {
            if (est[idx[1]] < 0) {
                r[i] = -1
            }
        }
    }
    for (i = 1; i <= rows(facs); i++) {
        idx = selectindex((lhs :== facs[i]) :& (op :== "=~") :& (substr(rhs, 1, 1) :== "y"))
        if (length(idx) == 0) {
            idx = selectindex((lhs :== facs[i]) :& (op :== "=~"))
            s = 0
            for (j = 1; j <= length(idx); j++) {
                s = s + tst_r(rhs[idx[j]], facs, r) * est[idx[j]]
            }
            if (s < 0) {
                r[i] = -1
            }
        }
    }
    refl = ""
    nref = 0
    for (i = 1; i <= rows(facs); i++) {
        if (r[i] < 0) {
            refl = refl + " " + facs[i]
            nref++
        }
    }

    // compare
    sg  = J(n, 1, 1)
    for (i = 1; i <= n; i++) {
        if (op[i] == "=~") {
            sg[i] = tst_r(lhs[i], facs, r) * tst_r(rhs[i], facs, r)
        }
        else if (op[i] == "~~" & lhs[i] != rhs[i]) {
            sg[i] = tst_r(lhs[i], facs, r) * tst_r(rhs[i], facs, r)
        }
    }
    re = 0
    rs = 0
    minse = .
    for (i = 1; i <= n; i++) {
        if (kk[i] >= .) {
            continue
        }
        k   = kk[i]
        a   = g[k]
        b   = sg[i] * strtoreal(fx[i, 5])
        se  = sqrt(V[k, k])
        lse = strtoreal(fx[i, 6])
        // max() skips missing: count a missing value as a failure explicitly
        if (missing(a) | missing(b) | missing(se) | missing(lse)) {
            printf("{err}    %s: %s has a missing value or SE\n", fitn, nm[k])
            bad++
            continue
        }
        re  = max((re, abs(a - b) / (1e-5 + 1e-4 * abs(b))))
        if (lse == 0) {
            if (V[k, k] != 0) {
                printf("{err}    %s: %s has SE 0 in lavaan, %g here\n", fitn, nm[k], se)
                bad++
            }
        }
        else {
            if (V[k, k] == 0) {
                printf("{err}    %s: %s has SE 0 here, %g in lavaan\n", fitn, nm[k], lse)
                bad++
            }
            rs = max((rs, abs(se - lse) / (1e-8 + 1e-5 * abs(lse))))
        }
        if (op[i] == "=~" & substr(rhs[i], 1, 1) == "F" & substr(lhs[i], 1, 1) == "G") {
            minse = min((minse, se))
        }
    }
    printf("{txt} 1. %-7s %f rows matched of %f (lavaan %f); value gap %8.2e of tier, SE gap %8.2e of tier\n",
           fitn, sum(used :> 0), rows(nm), n, re, rs)
    if (bad | sum(used :> 0) != rows(nm) | max(used) > 1 | n != rows(nm) | re > 1 | rs > 1 | re >= . | rs >= .) {
        tst_nfail++
    }
    printf("{txt}    %-7s reflected in lavaan's rows by the sign rule: %s; smallest SE of a second-order path %9.5f (information only)\n",
           fitn, (nref ? refl : "none"), minse)

    // 2: dwlsmv's own orientation, unreflected: a first item loading positive, a G's path sum positive
    bado = 0
    for (i = 1; i <= rows(facs); i++) {
        idx = selectindex((lhs :== facs[i]) :& (op :== "=~") :& (substr(rhs, 1, 1) :== "y"))
        if (length(idx) > 0) {
            if (!(gv[idx[1]] > 0)) {
                printf("{err}    %s: %s first item loading standardizes to %g, not positive\n", fitn, facs[i], gv[idx[1]])
                bado++
            }
        }
        else {
            idx = selectindex((lhs :== facs[i]) :& (op :== "=~"))
            sse = sum(gv[idx])
            if (!(sse > 0)) {
                printf("{err}    %s: the standardized paths of %s sum to %g, not positive\n", fitn, facs[i], sse)
                bado++
            }
        }
    }
    printf("{txt} 2. %-7s dwlsmv's orientation: %f violation(s)\n", fitn, bado)
    if (bado) {
        tst_nfail++
    }
}

// check 3: the standardized vector and covariance of fit a equal fit b's, matched by name
void tst_cross(string scalar a, string scalar b, real scalar assert1)
{
    external real scalar tst_nfail
    real rowvector ga, gb
    real matrix Va, Vb
    string colvector na, nb
    real colvector p, k
    real scalar i, gap1, gap2

    ga = st_matrix("tg_" + a)
    gb = st_matrix("tg_" + b)
    Va = st_matrix("tV_" + a)
    Vb = st_matrix("tV_" + b)
    na = tst_names(st_matrixcolstripe("tg_" + a))
    nb = tst_names(st_matrixcolstripe("tg_" + b))
    if (rows(na) != rows(nb)) {
        printf("{err} 3. %s vs %s: %f against %f standardized rows\n", a, b, rows(na), rows(nb))
        if (assert1) {
            tst_nfail++
        }
        return
    }
    p = J(rows(na), 1, .)
    for (i = 1; i <= rows(na); i++) {
        k = selectindex(nb :== na[i])
        if (length(k) != 1) {
            printf("{err} 3. %s vs %s: row %s is not in both\n", a, b, na[i])
            if (assert1) {
                tst_nfail++
            }
            return
        }
        p[i] = k[1]
    }
    gap1 = mreldif(ga, gb[., p'])
    gap2 = max(abs(Va - Vb[p, p])) / max(abs(Vb))
    printf("{txt} 3. %-4s vs %-4s standardized estimates mreldif %9.2e, covariance gap over max|V| %9.2e (1e-5)%s\n",
           a, b, gap1, gap2, (assert1 ? "" : "; information only"))
    if (assert1) {
        if (gap1 >= 1e-5 | gap2 >= 1e-5 | gap1 >= . | gap2 >= .) {
            tst_nfail++
        }
    }
}

// check 4: Case I, V
void tst_improper(string scalar fitn)
{
    external string matrix tst_fx
    external real scalar tst_nfail
    string matrix stp, fx
    string colvector nm
    real rowvector g
    real matrix V
    real scalar kb, kp, bstar, pstar, lb, lp, lsb, lsp, ok, d

    g   = st_matrix("tg_" + fitn)
    V   = st_matrix("tV_" + fitn)
    nm  = tst_names(st_matrixcolstripe("tg_" + fitn))
    fx  = select(tst_fx, tst_fx[., 1] :== fitn)
    kb  = selectindex(nm :== "F1:G")
    kp  = selectindex(nm :== "/var(e.F1)")
    lb  = strtoreal(select(fx[., 5], fx[., 2] :== "G" :& fx[., 3] :== "=~" :& fx[., 4] :== "F1"))
    lp  = strtoreal(select(fx[., 5], fx[., 2] :== "F1" :& fx[., 3] :== "~~" :& fx[., 4] :== "F1"))
    lsb = strtoreal(select(fx[., 6], fx[., 2] :== "G" :& fx[., 3] :== "=~" :& fx[., 4] :== "F1"))
    lsp = strtoreal(select(fx[., 6], fx[., 2] :== "F1" :& fx[., 3] :== "~~" :& fx[., 4] :== "F1"))
    bstar = g[kb]
    pstar = g[kp]
    ok = 1
    if (!(bstar > 1)) {
        printf("{err}    F1:G standardizes to %g, not above 1\n", bstar)
        ok = 0
    }
    if (!(pstar < 0)) {
        printf("{err}    var(e.F1) standardizes to %g, not negative\n", pstar)
        ok = 0
    }
    d = abs(pstar - (1 - bstar^2))
    if (!(d < 1e-12)) {
        printf("{err}    psi* - (1 - b*^2) = %g\n", d)
        ok = 0
    }
    if (!(abs(bstar - lb) <= 1e-5 + 1e-4 * abs(lb)) | !(abs(pstar - lp) <= 1e-5 + 1e-4 * abs(lp))) {
        printf("{err}    against lavaan: b* %g (%g), psi* %g (%g)\n", bstar, lb, pstar, lp)
        ok = 0
    }
    if (!(abs(sqrt(V[kb, kb]) - lsb) <= 1e-8 + 1e-5 * lsb) | !(abs(sqrt(V[kp, kp]) - lsp) <= 1e-8 + 1e-5 * lsp)) {
        printf("{err}    SEs %g (%g), %g (%g)\n", sqrt(V[kb, kb]), lsb, sqrt(V[kp, kp]), lsp)
        ok = 0
    }
    printf("{txt} 4. %s improper: b*(F1<-G) %9.5f (lavaan %9.5f), psi*(F1) %9.5f (lavaan %9.5f), psi* - (1 - b*^2) %8.1e\n",
           fitn, bstar, lb, pstar, lp, d)
    if (!ok) {
        tst_nfail++
    }
}
end

* ---- one fit: option runs, then check 1, 2 ---------------------------------------
capture program drop tt_fit
program tt_fit
    args tag model met
    capture noisily quietly dwlsmv `model', metric(`met') standardized
    if _rc | e(converged) != 1 {
        display as error "  `tag': dwlsmv rc " _rc ", converged " e(converged)
        mata: tst_nfail++
        exit
    }
    mata: tst_vs_lavaan("`tag'")
end

local m4  "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (F4 -> y10 y11 y12) (G -> F1 F2 F3 F4)"
local m3  "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (G -> F1 F2 F3)"
local mD  "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (F4 -> y10 y11 y12) (G -> y13 F1 F2 F3 F4)"
local mJ  "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (F4 -> y10 y11 y12) (F5 -> y13 y14 y15) (F6 -> y16 y17 y18) (G1 -> F1 F2 F3) (G2 -> F4 F5 F6)"

* ---- Case H ------------------------------------------------------------------------
import delimited using "fixtures/data/case_H.csv", clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
tt_fit H4_V `"`m4'"' marker
tt_fit H4_U `"`m4'"' unitvariance
tt_fit H3_V `"`m3'"' marker
tt_fit H3_U `"`m3'"' unitvariance
tt_fit HD_M `"`mD'"' marker
tt_fit HD_U `"`mD'"' unitvariance

* ---- Case I (improper) ---------------------------------------------------------------
import delimited using "fixtures/data/case_I.csv", clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
tt_fit I_V `"`m3'"' marker

* ---- Case J ----------------------------------------------------------------------------
import delimited using "fixtures/data/case_J.csv", clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y*, replace
tt_fit J_V `"`mJ'"' marker
tt_fit J_U `"`mJ'"' unitvariance

* ---- 3. metric invariance ----------------------------------------------------------------
mata: tst_cross("H4_V", "H4_U", 1)
mata: tst_cross("H3_V", "H3_U", 1)
mata: tst_cross("J_V",  "J_U",  1)
mata: tst_cross("HD_M", "HD_U", 0)

* ---- 4. improper I_V -----------------------------------------------------------------------
mata: tst_improper("I_V")

adopath - "`here'/src/ado"
adopath - "`here'/build"

mata: st_local("nfail", strofreal(tst_nfail))
if `nfail' {
    display as error _n "test_standardized_so: `nfail' failure(s)"
    exit 9
}
display as result "test_standardized_so: the standardized solution of nine fits, Case I's improper one among them, matches lavaan, dwlsmv's own sign rule holds, and the V and U standardized vectors and covariances agree"
exit 0
