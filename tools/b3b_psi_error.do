*! B5 / gate B3b: numerical error in a fitted psi-hat (release B). A diagnostic,
*! not a pass/fail test. Plan: E-114 (the psi measurement),
*! written before the first run.
*!
*! Executed under StataMP 18 on 2026-11-20, at six processors and at one; results
*! in E-114.
*!
*! Usage, from the repository root:  do tools/b3b_psi_error.do [tag]
*! Writes docs/notes/second-order-id/b5/psi_error[tag].csv. Run once at the
*! default processor count (no tag) and once at one processor (tag p1, from a
*! %TEMP% wrapper that does 'set processors 1') and diff the two files.
*!
*! Every fit is from hand-written canonical text: Module 11 still refuses a
*! factor with no observed indicator, and this goes through Module 3 directly.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255
args tag

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    quietly do "`f'"
}

foreach c in E H I J {
    import delimited using "fixtures/data/case_`c'.csv", clear varnames(1) ///
        delimiters(",") stringcols(_all)
    destring _all, replace
    mata: Y`c' = st_data(., .)
}
local out "docs/notes/second-order-id/b5/psi_error`tag'.csv"

mata:

string scalar bp_join(string colvector L)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= rows(L); i++) {
        s = s + L[i] + char(10)
    }
    return(s)
}

string scalar bp_g(real scalar g, real scalar nG)
{
    if (nG == 1) {
        return("G")
    }
    return("G" + strofreal(g))
}

// Metric V: first item of each F_k its marker, psi_kk free, beta F_k G free;
// each G by 'metric unitvariance' unless gmark (then the caller anchors it).
string colvector bp_v(real rowvector nper, real rowvector parent, real scalar nG,
                      real scalar xitems, real scalar gmark)
{
    string colvector L
    string scalar fl, il, f
    real scalar k, j, y, g

    fl = "factor"
    for (k = 1; k <= cols(nper); k++) {
        fl = fl + " F" + strofreal(k)
    }
    for (g = 1; g <= nG; g++) {
        fl = fl + " " + bp_g(g, nG)
    }
    il = "ind"
    for (y = 1; y <= sum(nper) + xitems; y++) {
        il = il + " y" + strofreal(y)
    }
    L = (fl \ il)
    y = 0
    for (k = 1; k <= cols(nper); k++) {
        f = "F" + strofreal(k)
        for (j = 1; j <= nper[k]; j++) {
            y++
            if (j == 1) {
                L = L \ ("metric marker " + f + " y" + strofreal(y))
            }
            else {
                L = L \ ("lambda " + f + " y" + strofreal(y) + " free")
            }
        }
        L = L \ ("psi " + f + " " + f + " free")
        if (parent[k] > 0) {
            L = L \ ("beta " + f + " " + bp_g(parent[k], nG) + " free")
        }
    }
    if (!gmark) {
        for (g = 1; g <= nG; g++) {
            L = L \ ("metric unitvariance " + bp_g(g, nG))
        }
    }
    return(L)
}

// Central difference of sigma at step hb * max(1, |v|): exact up to rounding,
// sigma being at most quadratic in each coordinate (SO1.4).
real matrix bp_dref(real colvector th, struct dwlsmv_parmap scalar par,
                    struct dwlsmv_statmap vector maps, real scalar hb)
{
    real matrix R
    real colvector tp, tm
    real scalar j, h

    R = J(rows(dwls_model(th, par, maps)), rows(th), .)
    for (j = 1; j <= rows(th); j++) {
        h = hb * max((1, abs(th[j])))
        tp = th
        tm = th
        tp[j] = th[j] + h
        tm[j] = th[j] - h
        R[., j] = (dwls_model(tp, par, maps) - dwls_model(tm, par, maps)) / (2 * h)
    }
    return(R)
}

// Fit, measure, write one row per free diagonal psi. Returns theta-hat.
real colvector bp_measure(string scalar tag, struct dwlsmv_stage1 scalar st,
                          string scalar canon, real scalar fh,
                          real scalar nrest, real colvector truth)
{
    struct dwlsmv_fit scalar f, fk
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    real colvector th, s, w, sig, r, g, e1, e1g, e2, se, tp, tm, st0, dev
    string scalar sdev, strue
    real matrix D, HGN, HGNi, Hex, Hexi, C, Dr1, Dr2, E
    real scalar n, j, h, q, qg, rho, k, i, nconv

    sts     = dwlsmv_stage1(1, 1)
    sts[1]  = st
    maps    = dwlsmv_statmap(1, 1)
    maps[1] = st.map
    par     = dwls_identify(canon, st.map)

    f = dwls_fit(sts, par)
    dwls_inference(f)
    th = f.theta
    n  = rows(th)
    se = sqrt(diagonal(f.acov) :/ (st.N - 1))

    s   = st.s
    w   = 1 :/ diagonal(st.Gamma)
    sig = dwls_model(th, par, maps)
    r   = s - sig
    D   = dwls_delta(th, par, maps)
    g   = -2 * D' * (w :* r)
    HGN = 2 * cross(D, w, D)

    // The exact Hessian: 2 Delta'W Delta - 2 sum_d w_d r_d d2sigma_d.
    C = J(n, n, 0)
    for (j = 1; j <= n; j++) {
        h  = 1e-4 * max((1, abs(th[j])))
        tp = th
        tm = th
        tp[j] = th[j] + h
        tm[j] = th[j] - h
        C[., j] = (dwls_delta(tp, par, maps) - dwls_delta(tm, par, maps))' * (w :* r) / (2 * h)
    }
    C   = (C + C') / 2
    Hex = HGN - 2 * C
    Hexi = luinv(Hex)
    HGNi = luinv(HGN)

    e1  = Hexi * g
    q   = g' * Hexi * g
    e1g = HGNi * g
    qg  = g' * HGNi * g

    Dr1 = bp_dref(th, par, maps, 1e-3)
    Dr2 = bp_dref(th, par, maps, 2e-3)
    E   = D - Dr1
    e2  = 2 * Hexi * (E' * (w :* r))
    rho = max(abs(eigenvalues(I(n) - HGNi * Hex)))

    dev   = J(n, 1, 0)
    nconv = 0
    for (k = 1; k <= nrest; k++) {
        st0 = th + 0.04 * (runiform(n, 1) :- 0.5)
        fk  = dwlsmv_fit_from(sts, par, st0)
        nconv = nconv + fk.converged
        dev = rowmax((dev, abs(fk.theta - th)))
    }

    printf("{txt}%s: converged %f, fmin %g, ||r||_W %g, rho %g, max|E| %g, Dref check %g, restarts converged %f of %f\n",
           tag, f.converged, f.fmin, sqrt(r' * (w :* r)), rho,
           max(abs(E)), max(abs(Dr1 - Dr2)), nconv, nrest)

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.row[i] == par.col[i] & par.free[i]) {
            j = par.fidx[i]
            // if/else, not ?:, so dev[j] and truth[j] are never indexed when
            // they do not apply (whether Mata's ?: evaluates only one branch
            // is not verified here).
            sdev = "NA"
            if (nrest > 0) {
                sdev = strofreal(dev[j], "%12.4e")
            }
            strue = "NA"
            if (rows(truth) > 0) {
                strue = strofreal(abs(th[j] - truth[j]), "%12.4e")
            }
            // The name is quoted: psi[F1,F1] carries a comma.
            fput(fh, tag + "," + char(34) + par.name[i] + char(34) + "," +
                 strofreal(th[j], "%21.15g") + "," + strofreal(se[j], "%12.6g") + "," +
                 strofreal(abs(e1[j]), "%12.4e") + "," +
                 strofreal(sqrt(abs(Hexi[j, j] * q)), "%12.4e") + "," +
                 strofreal(abs(e1g[j]), "%12.4e") + "," +
                 strofreal(sqrt(abs(HGNi[j, j] * qg)), "%12.4e") + "," +
                 strofreal(abs(e2[j]), "%12.4e") + "," +
                 strtrim(sdev) + "," + strtrim(strue) + "," +
                 strofreal(f.converged) + "," + strofreal(sqrt(r' * (w :* r)), "%12.4e") + "," +
                 strofreal(rho, "%12.4e") + "," + strofreal(max(abs(E)), "%12.4e"))
            printf("{txt}  %-10s psi-hat %12.6g (SE %8.4g)  stop %9.2e [bound %9.2e; GN %9.2e]  fd %9.2e  restarts %s  truth %s\n",
                   par.name[i], th[j], se[j], abs(e1[j]), sqrt(abs(Hexi[j, j] * q)),
                   abs(e1g[j]), abs(e2[j]), strtrim(sdev), strtrim(strue))
        }
    }
    return(th)
}

void bp_main(string scalar out)
{
    external real matrix YE, YH, YI, YJ
    struct dwlsmv_stage1 scalar sE, sH4, sH3, sMD, sI, sJ, sZ
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    real colvector thH, thZ, none
    real scalar fh, i, z, kpsi
    string scalar cH4
    real rowvector zval
    string rowvector ztag

    rseed(20261120)
    if (fileexists(out)) {
        unlink(out)
    }
    fh = fopen(out, "w")
    fput(fh, "fit,param,psihat,se,stop_exact,stop_bound_exact,stop_gn,stop_bound_gn,fd_term,restart_maxdev,truth_dev,converged,rnorm_w,rho,max_abs_E")
    none = J(0, 1, .)

    // Case E, marker metric: an ordinary model's endogenous psi, for comparison.
    sE = dwls_sample_stats(YE, J(6, 1, 3))
    (void) bp_measure("E_marker", sE, bp_join(("factor F1 F2" \ "ind y1 y2 y3 y4 y5 y6" \
        "metric marker F1 y1" \ "lambda F1 y2 free" \ "lambda F1 y3 free" \
        "metric marker F2 y4" \ "lambda F2 y5 free" \ "lambda F2 y6 free" \
        "psi F1 F1 free" \ "psi F2 F2 free" \ "beta F2 F1 free")), fh, 10, none)

    sH4 = dwls_sample_stats(YH[., 1..12], J(12, 1, 4))
    cH4 = bp_join(bp_v((3, 3, 3, 3), (1, 1, 1, 1), 1, 0, 0))
    thH = bp_measure("H4_V", sH4, cH4, fh, 10, none)

    sH3 = dwls_sample_stats(YH[., 1..9], J(9, 1, 4))
    (void) bp_measure("H3_V", sH3, bp_join(bp_v((3, 3, 3), (1, 1, 1), 1, 0, 0)), fh, 10, none)

    // MD: G anchored on y13, as sem does, so psi_G is free too.
    sMD = dwls_sample_stats(YH, J(13, 1, 4))
    (void) bp_measure("MD", sMD, bp_join(bp_v((3, 3, 3, 3), (1, 1, 1, 1), 1, 1, 1) \
        "metric marker G y13" \ "psi G G free"), fh, 10, none)

    sJ = dwls_sample_stats(YJ, J(18, 1, 4))
    (void) bp_measure("J_V", sJ, bp_join(bp_v((3, 3, 3, 3, 3, 3), (1, 1, 1, 2, 2, 2), 2, 0, 0) \
        "psi G2 G1 free"), fh, 10, none)

    sI = dwls_sample_stats(YI, J(9, 1, 4))
    (void) bp_measure("I_V", sI, bp_join(bp_v((3, 3, 3), (1, 1, 1), 1, 0, 0)), fh, 10, none)

    // Known truth: s = sigma(theta*) exactly, theta* = H4's estimate with
    // psi_F1 set to 0, -0.01, +0.01; Gamma (so W) from Case H.
    maps    = dwlsmv_statmap(1, 1)
    maps[1] = sH4.map
    par     = dwls_identify(cH4, sH4.map)
    kpsi    = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.row[i] == 1 & par.col[i] == 1) {
            kpsi = par.fidx[i]
        }
    }
    zval = (0, -0.01, 0.01)
    ztag = ("Z0_psi0", "Zm_psi-.01", "Zp_psi+.01")
    for (z = 1; z <= 3; z++) {
        thZ       = thH
        thZ[kpsi] = zval[z]
        sZ        = sH4
        sZ.s      = dwls_model(thZ, par, maps)
        (void) bp_measure(ztag[z], sZ, cH4, fh, 0, thZ)
    }
    fclose(fh)
    printf("{txt}written: %s\n", out)
}

end

mata: bp_main(st_local("out"))
