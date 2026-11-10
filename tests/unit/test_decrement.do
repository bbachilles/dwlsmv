*! The Newton decrement at the returned theta-hat: g' H^-1 g <= nrtol (1e-12)
*!
*! DIFFTEST check 10's band eps = (N-1) nrtol rests on dwlsmv_minimize()
*! stopping on optimize()'s Newton-decrement criterion (derivation
*! @sec-difftest, *The precision band eps*). E-003 (2026-10-08) recorded
*! as NOT established that optimize() evaluates that decrement at the theta-hat
*! it returns. This test asserts it: for every fit, dwlsmv_fdwls() at
*! fit.theta (todo = 2: the analytic gradient and the Gauss-Newton H that
*! optimize() itself used) gives g H^-1 g' <= 1e-12, the nrtol set in
*! dwlsmv_minimize(). No tolerance is chosen here; 1e-12 is the optimizer's own.
*!
*! Fits: test_fit.do's five (case D marker; case E correlated and structural,
*! unit variance and marker) and the four DIFFTEST models of
*! fixtures/generate/97_lavaan_difftest.R (beta0 and xload2, H0 and H1, case
*! E, marker metric; cor_marker and str_marker are two of them): seven fits.
*! Sign canonicalization after convergence is a reflection,
*! under which g H^-1 g' is invariant.
*!
*! Executed under StataMP 18 on 2026-10-08 and it passes: every decrement at the
*! returned theta-hat is <= 2.0e-18 (beta0_h0 8.5e-32), F there equals fit.fmin.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_decrement: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

string scalar tdc_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

string scalar tdc_canon(string scalar fit)
{
    string rowvector head, lam6, mk

    head = ("factor f1 f2", "ind y1 y2 y3 y4 y5 y6")
    lam6 = ("lambda f1 y1 free", "lambda f1 y2 free", "lambda f1 y3 free",
            "lambda f2 y4 free", "lambda f2 y5 free", "lambda f2 y6 free")
    mk   = ("lambda f1 y2 free", "lambda f1 y3 free",
            "lambda f2 y5 free", "lambda f2 y6 free",
            "psi f1 f1 free", "psi f2 f2 free",
            "metric marker f1 y1", "metric marker f2 y4")

    if (fit == "D_marker") {
        return(tdc_spec(("factor f", "ind y1 y2 y3 y4",
                         "lambda f y2 free", "lambda f y3 free",
                         "lambda f y4 free", "psi f f free",
                         "metric marker f y1")))
    }
    if (fit == "cor_unitvar") {
        return(tdc_spec((head, lam6, "psi f1 f2 free",
                         "metric unitvariance f1", "metric unitvariance f2")))
    }
    if (fit == "str_unitvar") {
        return(tdc_spec((head, lam6, "beta f2 f1 free",
                         "metric unitvariance f1", "metric unitvariance f2")))
    }
    // cor_marker is also xload2's H0; str_marker is beta0's H1
    if (fit == "cor_marker") {
        return(tdc_spec((head, mk, "psi f1 f2 free")))
    }
    if (fit == "str_marker") {
        return(tdc_spec((head, mk, "beta f2 f1 free")))
    }
    // beta0's H0: f2 ~ 0*f1, no covariance (absent = fixed 0)
    if (fit == "beta0_h0") {
        return(tdc_spec((head, mk)))
    }
    // xload2's H1: + y4 on f1 and y3 on f2
    if (fit == "xload2_h1") {
        return(tdc_spec((head, mk, "psi f1 f2 free",
                         "lambda f1 y4 free", "lambda f2 y3 free")))
    }
    _error("tdc_canon: unknown fit")
}

real scalar tdc_main(real matrix YD, real matrix YE)
{
    struct dwlsmv_stage1 scalar stD, stE, st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    string colvector FAIL
    string rowvector fits
    real colvector w
    real rowvector g
    real matrix H
    real scalar k, v, dec

    FAIL = J(0, 1, "")
    stD  = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    stE  = dwls_sample_stats(YE, J(6, 1, 3))
    fits = ("D_marker", "cor_unitvar", "str_unitvar", "cor_marker", "str_marker",
            "beta0_h0", "xload2_h1")

    for (k = 1; k <= cols(fits); k++) {
        if (fits[k] == "D_marker") {
            st = stD
        }
        else {
            st = stE
        }
        par = dwls_identify(tdc_canon(fits[k]), st.map)
        fit = dwls_fit(st, par)
        w   = 1 :/ diagonal(st.Gamma)

        dwlsmv_fdwls(2, fit.theta', st.s, w, par, st.map, v, g, H)
        dec = g * invsym(H) * g'
        printf("{txt}  %-11s converged %g in %2.0f it; npar %2.0f; g H^-1 g' = %9.2e; max|g| %9.2e; F %21.15g vs %21.15g\n",
               fits[k], fit.converged, fit.iterations, par.npar, dec, max(abs(g)), v, fit.fmin)
        if (fit.converged != 1) FAIL = FAIL \ (fits[k] + ": not converged")
        if (!(abs(dec) <= 1e-12)) FAIL = FAIL \ (fits[k] + ": decrement " + strofreal(dec) + " > nrtol")
    }

    for (k = 1; k <= rows(FAIL); k++) {
        printf("{err}%s\n", FAIL[k])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(tdc_main(YD, YE)))

if `nfail' {
    display as error _n "test_decrement: `nfail' failure(s)"
    exit 9
}
display as result "test_decrement: g H^-1 g' <= nrtol = 1e-12 at every returned theta-hat"
exit 0
