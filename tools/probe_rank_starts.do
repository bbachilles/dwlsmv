*! dwlsmv -- the Delta rank check at default starting values: a random-intercept case
*!
*! A probe, not a test (tools/). Records an open item (E-061,
*! 2026-10-04, the user's to decide): spec 4.8 checks rank(Delta) at the STARTING
*! values, and at dwlsmv's default starts (every free loading 0.7, psi
*! covariances 0) this model, of full rank away from that point, reads rank 19
*! of 20. Model: F1, F2 unit variance, three items each; G a random intercept
*! (all six loadings fixed at 1, variance free), orthogonal to F1 and F2.
*! Four runs with psi[F2,F1] free (default starts; psi21 start .3; distinct
*! loading starts; both), and one with psi[F2,F1] fixed at 0 (default starts).
*! Case E data (generated without G): the runs show behaviour, not fit.
*! RESOLVED 2026-10-04 (spec 4.8, derivation Finding 5): dwls_identify() now
*! re-checks at a second point and starts the fit there when that is what
*! passes. The reference file is re-recorded under that rule: all five runs
*! converge, admissible, to one solution. The original record (default and
*! psi21 .3 starts at rank 19 of 20, rc 3498) is in commit d087927. Keeping
*! theta0 as the start was measured too: identification passed and the fit
*! did not converge (psi21 -> -25.5).
*! Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_rank_starts.do
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 200
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
quietly do src/mata/dwlsmv_types.mata
foreach f of local srcs {
    quietly do "src/mata/`f'"
}
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata:
string scalar ri_spec(string scalar lamF, string scalar psi21)
{
    string scalar s, nl
    nl = char(10)
    s = "factor F1 F2 G" + nl + "ind y1 y2 y3 y4 y5 y6" + nl
    s = s + "lambda F1 y1 free" + lamF + nl + "lambda F1 y2 free" + nl + "lambda F1 y3 free" + nl
    s = s + "lambda F2 y4 free" + lamF + nl + "lambda F2 y5 free" + nl + "lambda F2 y6 free" + nl
    s = s + "metric marker G y1" + nl
    s = s + "lambda G y2 fixed 1" + nl + "lambda G y3 fixed 1" + nl + "lambda G y4 fixed 1" + nl
    s = s + "lambda G y5 fixed 1" + nl + "lambda G y6 fixed 1" + nl
    s = s + "metric unitvariance F1" + nl + "metric unitvariance F2" + nl
    if (psi21 != "omit") s = s + "psi F2 F1 free" + psi21 + nl
    s = s + "psi G G free" + nl
    return(s)
}
void ri_try(string scalar label, string scalar lamF, string scalar psi21)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    real matrix A
    real scalar rc, i
    st = dwls_sample_stats(st_data(., .), J(6, 1, 3))
    printf("RI %s\n", label)
    par = dwls_identify(ri_spec(lamF, psi21), st.map)
    printf("RI   identify ok, %f free\n", par.npar)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.ptype[i] != 4) printf("RI   start %-14s %8.4f\n", par.name[i], par.value[i])
    }
    fit = dwls_fit(st, par)
    printf("RI   converged %f, T3 %9.4f, df %f\n", fit.converged, fit.T3, fit.T3_df)
    A = dwlsmv_admissible(fit.theta, fit.par, fit.stage1[1].map)
    printf("RI   inadmissible flags %f %f %f %f %f %f; min Theta_jj %8.4f\n",
           A[1,1], A[2,1], A[3,1], A[4,1], A[5,1], A[6,1], A[1,2])
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i] & par.ptype[i] != 4) printf("RI   est   %-14s %8.4f\n", par.name[i], fit.theta[par.fidx[i]])
    }
}
end
capture noisily mata: ri_try("default starts", "", "")
display "RIRC " _rc
capture noisily mata: ri_try("psi21 start .3 only", "", " .3")
display "RIRC " _rc
capture noisily mata: ri_try("distinct loading starts only", " .9", "")
display "RIRC " _rc
capture noisily mata: ri_try("both", " .9", " .3")
display "RIRC " _rc
capture noisily mata: ri_try("psi21 fixed at 0, default starts", "", "omit")
display "RIRC " _rc
