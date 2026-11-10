* export_sens_t3_in.do -- write dwlsmv_sens_t3()'s inputs and output, in %21x,
* for tools/export_sens_t3_ref.R (the R reference of tests/unit/test_sens_t3.do)
*
* Executed under StataMP 18 on 2026-10-16 (six processors), writing the inputs of tests/unit/sens_t3_ref.mata.
*
* Run from the repository root, with the output directory in $ST3DIR:
*     global ST3DIR "C:/Users/<you>/AppData/Local/Temp/st3ref"
*     do tools/export_sens_t3_in.do
* Optional: global ST3PATH "exact" exports the exact path's inputs (dwlsmv_sens_x,
* dwlsmv_sens_jmat_x, dwlsmv_sens_ddelta_x, dwlsmv_sens_t3_x; Result MS6). Unset or
* "fd", the finite-difference path, as before (2026-10-19, plan step 5).
*
* Writes $ST3DIR/in_<case>.csv (name,i,j,hex) for the four cases of
* tests/unit/sens_t3_setup.mata, with epsV = 2u and rhoT = 0.01. Design E is
* consistency-only data; this measures no processor-count residue.
*
* SPDX-License-Identifier: GPL-3.0-or-later
* SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

if `"$ST3DIR"' == "" {
    display as error "export_sens_t3_in: set the global ST3DIR to the output directory"
    exit 198
}
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index
do tests/unit/sens_unpolished.mata
do tests/unit/sens_t3_setup.mata

mata:

void st3x_put(real scalar fh, string scalar nm, real matrix x)
{
    real scalar i, j

    for (i = 1; i <= rows(x); i++) {
        for (j = 1; j <= cols(x); j++) {
            fput(fh, nm + "," + strofreal(i) + "," + strofreal(j) + "," + strofreal(x[i, j], "%21x"))
        }
    }
}

void st3x_export(string scalar cs, string scalar fn)
{
    struct dwlsmv_fit scalar f
    struct dwlsmv_statmap vector maps
    real colvector step, ds, eps, s, wf, sg
    real matrix dOm, sens, A, Om
    real rowvector out
    real scalar fh, N, c, nu, epsV, rhoT, full, exact

    exact = (st_global("ST3PATH") == "exact")
    epsV = epsilon(1)
    rhoT = 0.01
    full = st3_inputs(cs, f, step, ds, eps, dOm)
    dwlsmv_fit_stack(f.stage1, s, wf, N, maps)
    sg = dwls_model(f.theta, f.par, maps)
    Om = dwlsmv_omega(f.stage1)
    if (full) {
        if (exact) sens = dwlsmv_sens_x(f, step, epsV, A, c, nu, ds, eps)
        else       sens = dwlsmv_sens(f, step, epsV, A, c, nu, ds, eps)
    }
    else {
        ds   = J(rows(s), 1, 0)
        eps  = J(rows(s), 1, 0)
        dOm  = J(rows(s), rows(s), 0)
        if (exact) sens = dwlsmv_sens_x(f, step, epsV, A, c, nu)
        else       sens = dwlsmv_sens(f, step, epsV, A, c, nu)
        sens[., 5] = J(rows(sens), 1, 0)          // missing without (ds, eps); the routine reads it as 0
    }
    out = (exact ? st3_case_x(cs, epsV, rhoT) : st3_case(cs, epsV, rhoT))

    unlink(fn)
    fh = fopen(fn, "w")
    st3x_put(fh, "th", f.theta)
    st3x_put(fh, "Dl", f.Delta)
    st3x_put(fh, "Om", Om)
    st3x_put(fh, "s", s)
    st3x_put(fh, "sg", sg)
    st3x_put(fh, "wf", wf)
    st3x_put(fh, "N", N)
    st3x_put(fh, "G", length(maps))
    st3x_put(fh, "T3", f.T3)
    st3x_put(fh, "T_raw", f.T_raw)
    st3x_put(fh, "shift_a", f.shift_a)
    st3x_put(fh, "step", step)
    st3x_put(fh, "K", f.iterations)
    st3x_put(fh, "epsV", epsV)
    st3x_put(fh, "rhoT", rhoT)
    st3x_put(fh, "full", full)
    st3x_put(fh, "ds", ds)
    st3x_put(fh, "eps", eps)
    st3x_put(fh, "dOm", dOm)
    st3x_put(fh, "A", A)
    st3x_put(fh, "c", c)
    st3x_put(fh, "sens", sens)
    st3x_put(fh, "Eb", dwlsmv_sens_ebar(f, epsV))
    st3x_put(fh, "Jm", (exact ? dwlsmv_sens_jmat_x(f) : dwlsmv_sens_jmat(f)))
    st3x_put(fh, "dD", (exact ? dwlsmv_sens_ddelta_x(f) : dwlsmv_sens_ddelta(f)))
    st3x_put(fh, "out", out)
    fclose(fh)
    printf("{txt}export_sens_t3_in: %s -> %s\n", cs, fn)
}

end

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
mata: st3_setup()
foreach cs in metric_full scalar_full metric_none metric_bigeps {
    mata: st3x_export("`cs'", "$ST3DIR/in_`cs'.csv")
}
display as result "export_sens_t3_in: done"
