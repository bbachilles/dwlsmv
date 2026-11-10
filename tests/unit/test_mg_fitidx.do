*! test_mg_fitidx.do -- G > 1 fit indices and difftest's dCFI, dRMSEA
*! against lavaan 0.7.2 (tier 2; study E-018)
*!
*! Written 2026-10-31. Reference fixtures/expected/lavaan_mg_fitidx_fits.csv
*! and _steps.csv (99j_lavaan_mg_fitidx.R): lavaan's T3 inputs converted to
*! dwlsmv's T = (N - G) F, FI1 on that convention. Designs E, D (committed
*! mg_E.csv, mg_D.csv), P, N (mg_P.csv, mg_N.csv), one factor on y1-y5, four
*! levels each, metric(unitvariance).
*!
*! Decisive, in this order (components before differences; prereg 3a.2):
*!   a  each fit: e(rmsea_s), e(cfi_s), e(tli_s) against rmsea_NG, cfi_NG,
*!      tli_NG at the test-statistic tier (1e-8 abs + 1e-5 rel); missing on
*!      both sides agrees, on one side fails; e(df_T3) = df, e(df_T3_b) = DB.
*!   b  difftest on the four: r(models) rows named by model, rmsea and cfi
*!      as in a; r(steps) rows H0:H1 the expected pair, d and status exact.
*!   c  r(steps) dCFI, dRMSEA (H0 minus H1) within the sum of the two
*!      components' tolerances; missing exactly where the fixture is (F6).
*! Diagnostics, printed only: e(T3) against T3_NG, e(T3_b) against T3B_NG,
*! the largest gap ratio per quantity.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
quietly do tests/dwlsmv_assert.mata
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

import delimited using "fixtures/expected/lavaan_mg_fitidx_fits.csv", ///
    clear varnames(1) case(lower) delimiters(",") stringcols(_all)
mata: FK = st_sdata(., ("design", "level"))
mata: FV = strtoreal(st_sdata(., ("df", "db", "rmsea_ng", "cfi_ng", "tli_ng", "t3_ng", "t3b_ng")))
import delimited using "fixtures/expected/lavaan_mg_fitidx_steps.csv", ///
    clear varnames(1) case(lower) delimiters(",") stringcols(_all)
mata: SK = st_sdata(., ("design", "h0", "h1"))
mata: SV = strtoreal(st_sdata(., ("d", "dcfi_ng", "drmsea_ng")))

mata:
// gap ratio with explicit missing handling (max() skips missing): both
// missing 0, one missing ., else |a - b| / (ea + er |b|)
real scalar tmf_gap(real scalar a, real scalar b, real scalar ea, real scalar er)
{
    if (missing(a) & missing(b)) {
        return(0)
    }
    if (missing(a) | missing(b)) {
        return(.)
    }
    return(abs(a - b) / (ea + er * abs(b)))
}

real scalar tmf_tol(real scalar x)
{
    return(1e-8 + 1e-5 * abs(x))
}

// a: the active fit against its fixture row; returns failures
real scalar tmf_fit(string scalar ds, string scalar lv, real rowvector worst)
{
    real scalar i, bad, k, r
    real rowvector ours

    external real matrix FV
    external string matrix FK
    i = selectindex(FK[., 1] :== ds :& FK[., 2] :== lv)
    if (rows(i) != 1) {
        printf("{err}check a FAILED: %s %s has %f fixture rows\n", ds, lv, rows(i))
        return(1)
    }
    bad  = 0
    ours = (st_numscalar("e(rmsea_s)"), st_numscalar("e(cfi_s)"), st_numscalar("e(tli_s)"))
    for (k = 1; k <= 3; k++) {
        r = tmf_gap(ours[k], FV[i, 2 + k], 1e-8, 1e-5)
        if (!(r <= 1)) {
            printf("{err}check a FAILED: %s %s index %f: %21.15g vs %21.15g (ratio %g)\n",
                   ds, lv, k, ours[k], FV[i, 2 + k], r)
            bad++
        }
        else if (r > worst[k]) {
            worst[k] = r
        }
    }
    if (st_numscalar("e(df_T3)") != FV[i, 1] | st_numscalar("e(df_T3_b)") != FV[i, 2]) {
        printf("{err}check a FAILED: %s %s df %g / %g, df_b %g / %g\n", ds, lv,
               st_numscalar("e(df_T3)"), FV[i, 1], st_numscalar("e(df_T3_b)"), FV[i, 2])
        bad++
    }
    printf("{txt}diag %s %-10s T3 rel %9.2e  T3_b rel %9.2e  rmsea %8.5f cfi %8.5f tli %8.5f\n",
           ds, lv, reldif(st_numscalar("e(T3)"), FV[i, 6]), reldif(st_numscalar("e(T3_b)"), FV[i, 7]),
           ours[1], ours[2], ours[3])
    return(bad)
}

// b and c: difftest's r(models), r(steps) for design ds; returns failures
real scalar tmf_dt(string scalar ds, real rowvector worst)
{
    real scalar bad, k, i, j, r, t, lo, hi
    real matrix M, S
    string matrix mr, sr

    external real matrix FV, SV
    external string matrix FK, SK
    M  = st_matrix("TM")
    S  = st_matrix("TS")
    mr = st_matrixrowstripe("TM")
    sr = st_matrixrowstripe("TS")
    bad = 0
    // b: models, most restricted first, each row named by its model
    for (k = 1; k <= rows(M); k++) {
        i = selectindex(FK[., 1] :== ds :& ("m_" :+ FK[., 2]) :== mr[k, 2])
        if (rows(i) != 1) {
            printf("{err}check b FAILED: %s r(models) row %f named %s\n", ds, k, mr[k, 2])
            bad++
            continue
        }
        for (j = 1; j <= 2; j++) {
            // r(models) cols 4 rmsea, 5 cfi; fixture cols 3 rmsea, 4 cfi
            if (!(tmf_gap(M[k, 3 + j], FV[i, 2 + j], 1e-8, 1e-5) <= 1)) {
                printf("{err}check b FAILED: %s r(models) %s col %f %g vs %g\n", ds, mr[k, 2], 3 + j,
                       M[k, 3 + j], FV[i, 2 + j])
                bad++
            }
        }
    }
    // b and c: each step against the fixture row of its own pair
    if (rows(S) != 3) {
        printf("{err}check b FAILED: %s r(steps) has %f rows\n", ds, rows(S))
        return(bad + 1)
    }
    for (k = 1; k <= 3; k++) {
        i = selectindex(SK[., 1] :== ds :& ("m_" :+ SK[., 2]) :== sr[k, 1] :& ("m_" :+ SK[., 3]) :== sr[k, 2])
        if (rows(i) != 1) {
            printf("{err}check b FAILED: %s r(steps) row %f is %s:%s, no such fixture step\n",
                   ds, k, sr[k, 1], sr[k, 2])
            bad++
            continue
        }
        // the row must be the k-th step in the fixture's order as well
        if (i != selectindex(SK[., 1] :== ds)[k]) {
            printf("{err}check b FAILED: %s r(steps) row %f %s:%s out of order\n", ds, k, sr[k, 1], sr[k, 2])
            bad++
        }
        t = (SV[i, 1] == 0 ? 0 : 1)
        if (S[k, 1] != SV[i, 1] | S[k, 6] != t) {
            printf("{err}check b FAILED: %s %s:%s d %g (fixture %g), status %g\n", ds, sr[k, 1], sr[k, 2],
                   S[k, 1], SV[i, 1], S[k, 6])
            bad++
        }
        // c: the tolerance is the two components' (H0 and H1 rows of r(models))
        lo = selectindex(mr[., 2] :== sr[k, 1])
        hi = selectindex(mr[., 2] :== sr[k, 2])
        for (j = 1; j <= 2; j++) {
            // dCFI: r(steps) col 4, r(models) col 5; dRMSEA: col 5, col 4
            r = (j == 1 ? 5 : 4)
            if (missing(SV[i, 1 + j]) & missing(S[k, 3 + j])) {
                continue
            }
            if (missing(SV[i, 1 + j]) | missing(S[k, 3 + j]) | rows(lo) != 1 | rows(hi) != 1) {
                printf("{err}check c FAILED: %s %s:%s diff %f %g vs fixture %g\n", ds, sr[k, 1], sr[k, 2], j,
                       S[k, 3 + j], SV[i, 1 + j])
                bad++
                continue
            }
            t = abs(S[k, 3 + j] - SV[i, 1 + j]) / (tmf_tol(M[lo, r]) + tmf_tol(M[hi, r]))
            if (!(t <= 1)) {
                printf("{err}check c FAILED: %s %s:%s diff %f %21.15g vs %21.15g (ratio %g)\n", ds,
                       sr[k, 1], sr[k, 2], j, S[k, 3 + j], SV[i, 1 + j], t)
                bad++
            }
            else if (t > worst[3 + j]) {
                worst[3 + j] = t
            }
            printf("{txt}diag %s %-11s:%-11s d%s %12.8f (fixture %12.8f)\n", ds, sr[k, 1], sr[k, 2],
                   (j == 1 ? "CFI  " : "RMSEA"), S[k, 3 + j], SV[i, 1 + j])
        }
    }
    return(bad)
}
end

local nfail 0
mata: W = J(1, 5, 0)
foreach ds in E D P N {
    import delimited using "fixtures/data/mg_`ds'.csv", clear varnames(1) delimiters(",")
    foreach L in configural thresholds metric scalar {
        capture noisily dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(`L') metric(unitvariance)
        if _rc | e(converged) != 1 {
            display as error "check a FAILED: `ds' `L' did not fit (rc `=_rc')"
            local ++nfail
            continue
        }
        mata: st_local("b", strofreal(tmf_fit("`ds'", "`L'", W)))
        local nfail = `nfail' + `b'
        estimates store m_`L'
    }
    capture noisily difftest m_configural m_thresholds m_metric m_scalar
    if _rc {
        display as error "check b FAILED: `ds' difftest refused, rc `=_rc'"
        local ++nfail
        continue
    }
    matrix TM = r(models)
    matrix TS = r(steps)
    mata: st_local("b", strofreal(tmf_dt("`ds'", W)))
    local nfail = `nfail' + `b'
    estimates clear
}
mata: printf("{txt}largest gap ratios: rmsea %g cfi %g tli %g dCFI %g dRMSEA %g\n", W[1], W[2], W[3], W[4], W[5])

if `nfail' {
    display as error _n "test_mg_fitidx: `nfail' failure(s)"
    exit 9
}
display _n "test_mg_fitidx: all checks passed"
