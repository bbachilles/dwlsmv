*! Drafted 2026-09-28 by a read-only agent. Executed unchanged by the main
*! session on StataMP 18 with polychoric v1.5, 2026-09-28 (nrep 5, treps 5).
*! Results: E-089.
*!
*! dwlsmv -- spec 7.5's Stata-side comparator: Kolenikov's -polychoric- (v1.5)
*! against dwlsmv's Stage 1 thresholds and per-pair rho-hat
*!
*! A diagnostic, not a test: no pass/fail, and the suite runner ignores tools/.
*! SCOPE: thresholds and each pair's rho-hat ONLY. Never Gamma, never SEs.
*!
*!     cd <repo root>
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/compare_polychoric.do [nrep [treps]]
*!     grep "^KOL|" compare_polychoric.log
*!
*! nrep  = replicates per high-rho condition compared (default 5)
*! treps = timed repetitions per timed design (default 5), after one warm-up
*!
*! CLEAN ROOM. His command is called as a user calls it and only its returned
*! r(R) is read. No code from his files is in this driver. His behaviour is
*! taken from E-081.
*!
*! HIS BY-DESIGN DIFFERENCES (E-081), stated
*! before any comparison and printed as KOL|design| lines before any result:
*!   D1 thresholds are invnormal((cumulative count - 1/2)/N), not the ML
*!      margins invnormal(cumulative count / N) that dwlsmv, lavaan and Mplus
*!      use. O(1/N) by design.
*!   D2 he returns no thresholds. The "his" threshold below is his documented
*!      rule computed by this driver from the counts, not his output; it is
*!      checked only indirectly, through rho-hat (D3).
*!   D3 thresholds enter each pair, so his rho-hat moves by O(1/N) too.
*!   D4 rho optimized directly by -ml- (lf, numerical derivatives), no tanh,
*!      no 0.999 cap; start rho = 0; 50 iterations then a random search.
*!      -ml-'s convergence tolerance adds its own small difference.
*!   D5 boundary pairs (an empty cell, MLE |rho| = 1) usually come back
*!      MISSING in r(R), without a warning. dwlsmv reports them at the 0.999
*!      bound, flagged converged = 0. Neither is compared, neither resolved:
*!      the 0.999 bound and boundary pairs are open decisions.
*!   D6 listwise deletion by default; variables with 10 or more distinct
*!      values are treated as continuous. The designs here have no missing
*!      values and at most 4 categories (checked and printed).
*!   D7 only the LAST pair's scalars are stored; r(R) is the only matrix.
*!      His displayed X2/G2 labels are swapped (not used here).
*!   D8 his random-search fallback draws random numbers, so the seed is set
*!      before every call.
*!
*! WHAT "AGREEMENT" MEANS HERE. Not spec 7.2's tiers. Every difference is
*! printed with the size first-order theory predicts and the ratio
*! observed / predicted. No pass/fail.
*!   Threshold: d = invnormal(q - 1/(2N)) - invnormal(q)
*!              predicted  -1 / (2 N phi(tau))
*!   rho-hat:   rho(tau) solves s(rho, tau) = 0, so
*!              d rho = -(I_ra' d_a + I_rb' d_b) / I_rr
*!              with dwlsmv's own pair information (I_ra = E[s_rho s_a] =
*!              -E[d s_rho / d a], the information identity; pihat weights).
*!              The residual (observed - predicted) is also printed, against 1/N.
*!   This first-order expansion is written here, not in the derivation file:
*!   a diagnostic's expected size, not production mathematics.
*!
*! DESIGNS. case_D (fixtures/data/case_D.csv) and replicates 1..nrep of the
*! high-rho Monte Carlo conditions mc_c_98, mc_c_99, mc_c_995, mc_c_998,
*! mc_c_999, mc_B_99, mc_B_995, mc_B_9975 (fixtures/replication/, gitignored;
*! regenerate with fixtures/generate/50_replicates.R). An absent file is
*! skipped with a printed note. K = (4, 4, 3, 3) throughout (case D's items).
*!
*! TIMING. Median wall clock over treps repetitions after one untimed
*! warm-up, data loading excluded, on case_D and replicate 1 of each
*! condition. His: -quietly polychoric varlist, nolog-. Ours: the SAME scope,
*! every variable's thresholds plus every pair's dwls_polychoric() -- not
*! dwls_sample_stats(), which also builds Gamma. RUN ON AN IDLE MACHINE.
*!
*! Every result line starts KOL|.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
args nrep treps
if "`nrep'"  == "" local nrep  5
if "`treps'" == "" local treps 5
clear all
set more off
set linesize 255

capture which polychoric
if _rc {
    display "KOL|absent|polychoric is not installed: net from https://staskolenikov.net/stata/ (not on SSC)"
    exit
}

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
quietly do src/mata/dwlsmv_types.mata
foreach f of local srcs {
    quietly do "src/mata/`f'"
}

display "KOL|machine|stata `c(stata_version)' `c(flavor)' MP=`c(MP)'|processors `c(processors)' of licensed `c(processors_lic)', machine `c(processors_mach)'|`c(os)' `c(machine_type)'|host `c(hostname)'|`c(current_date)' `c(current_time)'"
display "KOL|design|D1 his thresholds: invnormal((cum count - .5)/N); dwlsmv: invnormal(cum count / N). O(1/N) by design"
display "KOL|design|D2 he returns no thresholds; his column is his documented rule computed here, checked only through rho"
display "KOL|design|D3 threshold difference propagates to rho at O(1/N); predicted by -(I_ra'da + I_rb'db)/I_rr"
display "KOL|design|D4 rho by ml directly, no tanh, no 0.999 cap, start 0, numerical derivatives; ml tolerance adds noise"
display "KOL|design|D5 boundary pairs: his rho usually missing; dwlsmv at 0.999 flagged converged 0. Reported, not compared, not resolved"
display "KOL|design|D6 listwise deletion; >= 10 distinct values means continuous. Missing values and category counts printed per design"
display "KOL|design|D7 only r(R) is read; his last-pair scalars and swapped X2/G2 labels are not used"
display "KOL|design|D8 seed set before every call of his command (random-search fallback)"
display "KOL|design|agreement is O(1/N), not spec 7.2 tiers: observed, predicted, ratio; never pass/fail"

mata:

real scalar kd_median(real colvector t)
{
    real colvector s
    real scalar n

    s = sort(t, 1)
    n = rows(s)
    return(mod(n, 2) ? s[(n + 1) / 2] : (s[n / 2] + s[n / 2 + 1]) / 2)
}

// His documented threshold rule (D1), from the counts.
real colvector kd_tauhis(real colvector y, real scalar K)
{
    real colvector n, cn

    n  = dwlsmv_catcounts(y, K)
    cn = runningsum(n)
    return(invnormal((cn[|1 \ K - 1|] :- 0.5) :/ rows(y)))
}

// Load Y, set up both sets of thresholds, print the threshold lines.
void kd_prep(string scalar tag)
{
    external real matrix KD_Y
    external real colvector KD_K
    external string rowvector KD_V
    external pointer(real colvector) colvector KD_TAU, KD_TAUH
    real colvector t, h
    real scalar j, c, N, d, pr

    KD_V = st_varname(1..st_nvar())
    KD_Y = st_data(., .)
    KD_K = (4 \ 4 \ 3 \ 3)
    N    = rows(KD_Y)
    printf("KOL|%s|N %g|missing values %g|distinct values per variable", tag, N, missing(KD_Y))
    for (j = 1; j <= cols(KD_Y); j++) printf(" %g", rows(uniqrows(KD_Y[., j])))
    printf("|1/N %9.3e\n", 1 / N)
    KD_TAU  = J(cols(KD_Y), 1, NULL)
    KD_TAUH = J(cols(KD_Y), 1, NULL)
    for (j = 1; j <= cols(KD_Y); j++) {
        t = dwlsmv_thresholds(KD_Y[., j], KD_K[j])
        h = kd_tauhis(KD_Y[., j], KD_K[j])
        KD_TAU[j]  = &(t[.])
        KD_TAUH[j] = &(h[.])
        for (c = 1; c <= rows(t); c++) {
            d  = h[c] - t[c]
            pr = -1 / (2 * N * normalden(t[c]))
            printf("KOL|%s|threshold %s t%g|dwlsmv %13.10f|his rule %13.10f|diff %10.3e|predicted %10.3e|ratio %8.4f\n",
                   tag, KD_V[j], c, t[c], h[c], d, pr, d / pr)
        }
    }
    displayflush()
}

// One pair: compare, or report why not.
void kd_pair(string scalar tag, real scalar a, real scalar b)
{
    external real matrix KD_Y
    external string rowvector KD_V
    external pointer(real colvector) colvector KD_TAU, KD_TAUH
    struct dwlsmv_pairfit scalar f
    real matrix R
    real scalar rk, N, pred, d, atb
    string scalar lab

    N   = rows(KD_Y)
    lab = tag + "|pair " + KD_V[a] + " " + KD_V[b]
    R   = st_matrix("KOL_R")
    rk  = R[a, b]
    f   = dwls_polychoric(KD_Y[., a], KD_Y[., b], *KD_TAU[a], *KD_TAU[b])
    atb = abs(f.rho) >= 0.9989
    if (f.converged == 0 | atb | missing(rk)) {
        printf("KOL|%s|not compared|dwlsmv rho %12.9f converged %g%s|his rho %s\n",
               lab, f.rho, f.converged, atb ? " at 0.999 bound" : "",
               missing(rk) ? "missing" : strofreal(rk, "%12.9f"))
        displayflush()
        return
    }
    pred = -(f.I_ra' * (*KD_TAUH[a] - *KD_TAU[a]) + f.I_rb' * (*KD_TAUH[b] - *KD_TAU[b])) / f.I_rr
    d    = rk - f.rho
    printf("KOL|%s|dwlsmv %12.9f|his %12.9f|diff %10.3e|predicted %10.3e|ratio %8.4f|residual %10.3e|residual*N %8.4f\n",
           lab, f.rho, rk, d, pred, d / pred, d - pred, (d - pred) * N)
    displayflush()
}

// Our side of the timing: thresholds plus every pair, the scope his command has.
void kd_ours()
{
    external real matrix KD_Y
    external real colvector KD_K
    pointer(real colvector) colvector th
    real scalar a, b

    th = J(cols(KD_Y), 1, NULL)
    for (a = 1; a <= cols(KD_Y); a++) {
        th[a] = &(dwlsmv_thresholds(KD_Y[., a], KD_K[a]))
    }
    for (a = 1; a < cols(KD_Y); a++) {
        for (b = a + 1; b <= cols(KD_Y); b++) {
            (void) dwls_polychoric(KD_Y[., a], KD_Y[., b], *th[a], *th[b])
        }
    }
}

void kd_time_ours(real scalar reps)
{
    external real colvector KD_TO
    real scalar r

    kd_ours()
    KD_TO = J(reps, 1, .)
    for (r = 1; r <= reps; r++) {
        timer_clear(92); timer_on(92)
        kd_ours()
        timer_off(92); KD_TO[r] = timer_value(92)[1]
    }
}

void kd_time_report(string scalar tag, real scalar reps)
{
    external real colvector KD_TO, KD_TH
    external real matrix KD_Y

    printf("KOL|%s|timing|p %g N %g|reps %g after 1 warm-up|dwlsmv %9.4f s|polychoric %9.4f s|ratio his/ours %8.2f\n",
           tag, cols(KD_Y), rows(KD_Y), reps, kd_median(KD_TO), kd_median(KD_TH),
           kd_median(KD_TH) / kd_median(KD_TO))
    displayflush()
}

end

* ---- one design: the data in memory, compare and optionally time -------------
program define kol_design
    args tag dotime treps
    unab vars : _all
    capture noisily mata: kd_prep("`tag'")
    if _rc {
        display "KOL|`tag'|not compared|dwlsmv thresholds aborted rc " _rc
        exit
    }
    set seed 20260928
    capture quietly polychoric `vars', nolog
    if _rc {
        display "KOL|`tag'|not compared|polychoric failed rc " _rc
        exit
    }
    matrix KOL_R = r(R)
    local rn : rownames KOL_R
    if "`rn'" != "`vars'" {
        display "KOL|`tag'|not compared|r(R) rows (`rn') are not the varlist (`vars')"
        exit
    }
    local p : word count `vars'
    forvalues a = 1/`p' {
        forvalues b = `=`a'+1'/`p' {
            capture noisily mata: kd_pair("`tag'", `a', `b')
            if _rc {
                local va : word `a' of `vars'
                local vb : word `b' of `vars'
                display "KOL|`tag'|pair `va' `vb'|not compared|dwlsmv aborted rc " _rc
            }
        }
    }
    if `dotime' {
        capture noisily mata: kd_time_ours(`treps')
        if _rc {
            display "KOL|`tag'|timing|dwlsmv aborted rc " _rc
            exit
        }
        set seed 20260928
        quietly polychoric `vars', nolog
        mata: KD_TH = J(`treps', 1, .)
        forvalues r = 1/`treps' {
            set seed 20260928
            timer clear 91
            timer on 91
            quietly polychoric `vars', nolog
            timer off 91
            quietly timer list 91
            mata: KD_TH[`r'] = st_numscalar("r(t91)")
        }
        mata: kd_time_report("`tag'", `treps')
    }
end

* ---- case D -------------------------------------------------------------------
import delimited using "fixtures/data/case_D.csv", clear varnames(1) ///
    delimiters(",") stringcols(_all)
unab cvars : _all
mata: KD_S = strtoreal(st_sdata(., .))
clear
mata: st_addobs(rows(KD_S))
mata: (void) st_addvar("double", tokens(st_local("cvars")))
mata: st_store(., ., KD_S)
kol_design case_D 1 `treps'

* ---- the high-rho Monte Carlo conditions --------------------------------------
foreach c in mc_c_98 mc_c_99 mc_c_995 mc_c_998 mc_c_999 mc_B_99 mc_B_995 mc_B_9975 {
    forvalues r = 1/`nrep' {
        local f : display "fixtures/replication/`c'/`c'_" %04.0f `r' ".dta"
        capture confirm file "`f'"
        if _rc {
            display "KOL|`c' rep `r'|absent: `f' (gitignored; regenerate with 50_replicates.R)"
            continue
        }
        quietly use "`f'", clear
        kol_design "`c' rep `r'" `=(`r' == 1)' `treps'
    }
}
display "KOL|done"
