*! dwlsmv -- check 6 sensitivity diagnostic, step 1: the gradient noise at a saved theta
*!
*! Executed under StataMP 18 on 2026-10-14 at six processors (default count only);
*! output and reading in E-082.
*!
*! E-027, 2026-10-14 (the controlled diagnostic, step 1). Fixture
*! E (consistency data) at the metric level, one processor count (the
*! default). No processor-count comparison is made here. This step needs no
*! frozen protocol (the same day's decision row).
*!
*! 1. theta-hat saved bitwise (fputmatrix), reloaded, compared.
*! 2. The reference derivative: dgs_delta_acc (H = 1e-2) against the closed
*!    form dgs_delta_an on the central-difference entries; acc at H = 1e-2
*!    against H = 5e-3; production against the closed form.
*! 3. The gradient noise e = (Delta_prod - Delta_ref)'W r, with Delta_ref the
*!    closed form, at theta-hat: its size, where it is largest, its Gauss-
*!    Newton step P^-1 e, and the same for the accurate derivative.
*! 4. Around element 17 and nearby: theta-hat moved by k ulps in element 17
*!    (k = -3..3), by delta in 17 (1e-12, 1e-10, 1e-8, 1e-6), and by ten
*!    random 1-ulp patterns on every element.
*! 5. S5's bound with epsV = 2u at every point of 4 (|E| <= bound), and the
*!    measured epsV / u.
*!
*! Run from the repository root:
*!     StataMP-64.exe /e do tools/diag_sens_step1.do
*! Writes SENS1| lines to the log and c(tmpdir)/sens_step1.csv.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index
quietly do tools/diag_sens_lib.do
display "SENS1|env|stata_version|`c(stata_version)'|processors|`c(processors)'"

mata:

real colvector ds1_ebound(real colvector th, struct dwlsmv_parmap scalar par,
                          struct dwlsmv_statmap vector maps, real matrix D, real scalar epsV,
                          real matrix Eb)
{
    real matrix Mk
    real colvector sg, h
    real scalar j, u

    u  = epsilon(1) / 2
    sg = dwls_model(th, par, maps)
    Mk = dgs_cdmask(par, maps)
    Eb = J(rows(D), cols(D), 0)
    h  = 1e-5 :* rowmax((J(rows(th), 1, 1), abs(th)))
    for (j = 1; j <= cols(D); j++) {
        Eb[., j] = Mk[., j] :* (epsV :* abs(sg) + u * (abs(th[j]) + h[j]) :* abs(D[., j])) / h[j]
    }
    return(sg)
}

void ds1_point(string scalar tag, real colvector th, real colvector s, real colvector w,
               struct dwlsmv_parmap scalar par, struct dwlsmv_statmap vector maps,
               real colvector th0, real scalar fh)
{
    real matrix Dp, Dr, Da, Eb, P, Mk, Sc
    real colvector r, e, ea, st, sg, h
    real scalar u, ratio, over, epsm, j

    u  = epsilon(1) / 2
    Dp = dwls_delta(th, par, maps)
    Dr = dgs_delta_an(th, par, maps)
    Da = dgs_delta_acc(th, par, maps, 1e-2)
    r  = s - dwls_model(th, par, maps)
    P  = quadcross(Dr, w, Dr)
    e  = quadcross(Dp - Dr, w :* r)
    ea = quadcross(Da - Dr, w :* r)
    st = cholsolve(P, e)
    sg = ds1_ebound(th, par, maps, Dp, 2 * u, Eb)
    Mk = dgs_cdmask(par, maps)
    // |E| against S5's bound at epsV = 2u, and E in units of the bound's scale, as mg-sens-d
    ratio = max(abs(Dp - Dr) :/ (Eb + (Eb :== 0)) :* Mk)
    over  = max(abs(Dp - Dr) - Eb)
    h  = 1e-5 :* rowmax((J(rows(th), 1, 1), abs(th)))
    Sc = J(rows(Dp), cols(Dp), 0)
    for (j = 1; j <= cols(Dp); j++) {
        Sc[., j] = (abs(sg) + (abs(th[j]) + h[j]) :* abs(Dp[., j])) / h[j]
    }
    epsm = max(abs(Dp - Dr) :/ (Sc + (Sc :== 0)) :* Mk) / u
    printf("SENS1|%s|maxabs_e|%9.3e|argmax_e|%f|e17|%10.3e|maxabs_step|%9.3e|argmax_step|%f|maxabs_e_acc|%9.3e|E_over_bound2u|%6.3f|eps_over_u|%6.3f|bound_excess|%9.2e|move|%9.2e\n",
           tag, max(abs(e)), order(-abs(e), 1)[1], e[17], max(abs(st)), order(-abs(st), 1)[1],
           max(abs(ea)), ratio, epsm, over, max(abs(th - th0)))
    fput(fh, tag + "," + strofreal(max(abs(e)), "%21x") + "," + strofreal(e[17], "%21x") + "," +
         strofreal(max(abs(st)), "%21x") + "," + strofreal(max(abs(ea)), "%21x") + "," + strofreal(ratio, "%21x") + "," + strofreal(epsm, "%21x"))
}

void ds1_run()
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    string scalar cc, fn
    real matrix K, Y, Da, Da2, Dr, Dp, Mk
    real colvector gv, s, w, th, thr, t2, d
    real scalar G, g, N, fh, k, j, i

    K   = st_matrix("e(K)")
    cc  = st_global("e(canon_configural)")
    G   = rows(K)
    gv  = uniqrows(st_data(., "g"))
    sts = dwlsmv_stage1(G, 1)
    maps = dwlsmv_statmap(G, 1)
    for (g = 1; g <= G; g++) {
        Y       = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
        sts[g]  = dwls_sample_stats(Y, K[g, .]')
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(cc, K, "metric"), maps, sts)
    fit = dwls_fit(sts, par)
    dwlsmv_fit_stack(sts, s, w, N, maps)

    // 1. the saved parameter vector
    fn = pathjoin(c("tmpdir"), "sens_step1_theta.mmat")
    if (fileexists(fn)) unlink(fn)
    fh = fopen(fn, "w")
    fputmatrix(fh, fit.theta)
    fclose(fh)
    fh  = fopen(fn, "r")
    thr = fgetmatrix(fh)
    fclose(fh)
    th = thr
    printf("SENS1|saved|converged|%f|iterations|%f|npar|%f|reload_identical|%f\n",
           fit.converged, fit.iterations, rows(th), mreldif(thr, fit.theta) == 0 & thr == fit.theta)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.fidx[i] == 17) printf("SENS1|element17|row|%f|name|%s|grp|%f|value|%10.6f\n", i, par.name[i], par.grp[i], th[17])
    }

    // 2. the reference derivative
    Mk  = dgs_cdmask(par, maps)
    Dr  = dgs_delta_an(th, par, maps)
    Da  = dgs_delta_acc(th, par, maps, 1e-2)
    Da2 = dgs_delta_acc(th, par, maps, 5e-3)
    Dp  = dwls_delta(th, par, maps)
    printf("SENS1|reference|cd_entries|%f|acc_vs_closed|%9.3e|acc_H_vs_H2|%9.3e|prod_vs_closed|%9.3e|offmask_prod_vs_closed|%9.3e\n",
           sum(Mk), max(abs(Da - Dr)), max(abs(Da - Da2)), max(abs(Dp - Dr)), max(abs((Dp - Dr) :* (1 :- Mk))))

    // 3.-5. the gradient noise
    fn = pathjoin(c("tmpdir"), "sens_step1.csv")
    if (fileexists(fn)) unlink(fn)
    fh = fopen(fn, "w")
    fput(fh, "tag,maxabs_e,e17,maxabs_step,maxabs_e_acc,E_over_bound2u,eps_over_u")
    ds1_point("theta_hat", th, s, w, par, maps, th, fh)
    for (k = -3; k <= 3; k++) {
        if (k == 0) continue
        t2 = th
        for (j = 1; j <= abs(k); j++) t2[17] = t2[17] + sign(k) * epsilon(t2[17])
        ds1_point("ulp17_" + strofreal(k), t2, s, w, par, maps, th, fh)
    }
    d = (1e-12 \ 1e-10 \ 1e-8 \ 1e-6)
    for (k = 1; k <= rows(d); k++) {
        t2 = th
        t2[17] = th[17] + d[k]
        ds1_point("delta17_" + strofreal(d[k], "%5.0e"), t2, s, w, par, maps, th, fh)
    }
    for (k = 1; k <= 10; k++) {
        t2 = th + sign(sin((1::rows(th)) :* (k + 0.5))) :* epsilon(1) :* abs(th)
        ds1_point("ulpall_" + strofreal(k), t2, s, w, par, maps, th, fh)
    }
    fclose(fh)
}

end

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
mata: ds1_run()
display "SENS1|done"
