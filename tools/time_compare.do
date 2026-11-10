*! ======================================================================
*! Drafted 2026-09-28 by a read-only agent. Executed unchanged by the main
*! session on StataMP 18, 2026-09-28 (writedata, smoke, then the full run at
*! 6 and at 4 processors). Results: E-089.
*! Re-run 2026-10-05 in nogsem mode at 6 and 4 processors:
*! E-089.
*! ======================================================================
*!
*! dwlsmv -- spec 7.5's timing, Stata side: dwlsmv's TOTAL beside its Stage 1 /
*! Module 2 / Stage 2 breakdown, and gsem (ordinal probit) as a total.
*!
*! Binding protocol: E-008, 2026-09-28, "Spec 7.5's timing
*! protocol, agreed before any driver." A diagnostic, not a test: no
*! pass/fail, and the suite runner ignores tools/.
*!
*! TWO MODES
*!   writedata  simulate the 12 grid cells exactly as tools/time_baseline.do
*!              does (same functions, seed 20260926, same draw order) and write
*!              each to fixtures/replication/timing/grid_p<p>_K<K>_N<N>.csv
*!              (gitignored: fixtures/replication/ is). Each file is read back
*!              and compared with the matrix that was written. Replicate 1 of
*!              the five Monte Carlo conditions is copied from its .dta to
*!              mc_<cond>.csv beside them. Run ONCE; every tool (this file,
*!              tools/time_lavaan.R, tools/time_mplus.R) then reads these files.
*!   (none)     time. Reads the CSVs; never simulates.
*!   nogsem     the same, dwlsmv rows only; gsem is not run (added 2026-10-05).
*!
*!     cd <repo root>
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/time_compare.do writedata
*!     then, for 6 and for 4 processors, the %TEMP% wrappers:
*!     a one-line-per-command do-file that cds to the
*!     repo, sets processors, and runs  do tools/time_compare.do  ; its log is
*!     <wrapper>.log in %TEMP%, so the two runs do not overwrite each other.
*!     grep "^TIME|" <wrapper>.log
*!
*! Same data as the 2026-09-26 baseline: nothing in src/mata draws random
*! numbers (checked by grep for rseed/runiform/rnormal, 2026-09-28), so drawing
*! all 12 cells up front gives the same matrices as the baseline's interleaved
*! draw-then-time loop. The writedata log's cell list is the evidence.
*!
*! RUN IT ON AN IDLE MACHINE.
*!
*! WHAT IS TIMED
*!   dwlsmv  stage1 = dwls_sample_stats(); stage2 = dwls_identify() through
*!           dwls_inference(); module2 = every pair's dwls_polychoric() at the
*!           sample thresholds (as the baseline); TOTAL = one timer spanning
*!           stage1 and stage2 in the same repetition, reported as the median
*!           of per-repetition totals (not the sum of the two medians).
*!           Excluded: Mata compilation, reading the CSV into Y, one untimed
*!           warm-up per design, Module 10 diagnostics and display.
*!   gsem    one gsem call per fit, as a total; N = 1000 cells only (p 10, 20,
*!           40, i.e. 1, 2, 4 factors; both K). Excluded: import of the CSV,
*!           ONE untimed warm-up fit on the first cell only (see OPEN
*!           QUESTIONS), output display (quietly).
*!   Repetitions: median of 5 (3 at p 40, N 10000), as agreed.
*!
*! MODEL MATCH
*!   Grid: the true structure, unit-variance metric, factor correlations free,
*!   thresholds free. gsem: (Ff -> ten items), oprobit, var(Ff@1). gsem's
*!   default covstructure for latent exogenous variables leaves their
*!   covariances free. With a factor variance fixed, gsem should NOT also
*!   anchor the first loading; this is asserted per fit (y1's loading must
*!   carry a positive standard error), and e(rank) is printed beside dwlsmv's
*!   npar so a mismatch in the parameter count is visible.
*!   MC rows: the one-factor marker model of tools/paired_mc.do. No gsem row.
*!
*! WHERE THE TOTALS CANNOT BE MATCHED EXACTLY (state in the performance note)
*!   - gsem is full-information ML with numerical integration and a theta-type
*!     parameterization; dwlsmv is DWLS, delta. Different estimators; the
*!     total answers "how long does each take on the same data and model".
*!   - gsem's total includes its own data setup from the Stata dataset in
*!     memory (marking the sample, building its equations); dwlsmv's starts
*!     from Y already in Mata. Neither can be separated through gsem's public
*!     interface.
*!   - gsem computes starting values by fitting simpler models first; that is
*!     inside its total, as dwlsmv's own start values are inside Stage 2.
*!
*! gsem TIME CAP. Stata cannot interrupt a running estimation from a do-file.
*! The cap is therefore enforced two ways, both honest about what they bound:
*!   iterate(GITER) bounds iterations, not seconds; a fit that hits it is
*!     reported converged 0, never dropped.
*!   GCAP seconds: after any timed repetition longer than GCAP, the remaining
*!     repetitions of that design are skipped and the row says how many ran.
*!     A single fit can still exceed GCAP once.
*! GITER 300 (Stata's default maxiter) and GCAP 1800 are chosen BEFORE any
*! measurement; revise only with that label.
*!
*! OPEN QUESTIONS (for the user; not settled here)
*!   - gsem warm-up once (first cell) rather than per design, because a
*!     per-design warm-up doubles the cost of the 4-factor fits. The warm-up's
*!     purpose, loading gsem's ado and Mata code, is met by one fit.
*!   - intmethod(mvaghermite) intpoints(7) are gsem's defaults, pinned
*!     explicitly so they are recorded; whether to add a second, cheaper
*!     setting (e.g. intpoints(3)) is not decided.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
args mode
clear all
set more off
set linesize 255

local tdir "fixtures/replication/timing"
local mcs  "mc_d mc_a1p_16 mc_b1_m1 mc_c_998 mc_B_9975"
local giter 300
local gcap  1800
local gpts  7

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
quietly do src/mata/dwlsmv_types.mata
foreach f of local srcs {
    quietly do "src/mata/`f'"
}

display "TIME|machine|stata `c(stata_version)' `c(flavor)' MP=`c(MP)'|processors `c(processors)' of licensed `c(processors_lic)', machine `c(processors_mach)'|`c(os)' `c(machine_type)'|host `c(hostname)'|`c(current_date)' `c(current_time)'"

mata:

// ---- copied unchanged from tools/time_baseline.do ----------------------
string scalar tb_spec(real scalar p, real scalar nf, real scalar marker)
{
    string scalar s
    real scalar j, f, g

    s = "factor"
    for (f = 1; f <= nf; f++) s = s + " f" + strofreal(f)
    s = s + char(10) + "ind"
    for (j = 1; j <= p; j++) s = s + " y" + strofreal(j)
    s = s + char(10)
    for (j = 1; j <= p; j++) {
        f = ceil(j / (p / nf))
        if (!(marker & j == 1)) {
            s = s + "lambda f" + strofreal(f) + " y" + strofreal(j) + " free" + char(10)
        }
    }
    if (marker) {
        return(s + "psi f1 f1 free" + char(10) + "metric marker f1 y1" + char(10))
    }
    for (f = 1; f <= nf; f++) {
        for (g = f + 1; g <= nf; g++) {
            s = s + "psi f" + strofreal(f) + " f" + strofreal(g) + " free" + char(10)
        }
        s = s + "metric unitvariance f" + strofreal(f) + char(10)
    }
    return(s)
}

real matrix tb_sim(real scalar p, real scalar K, real scalar N)
{
    real matrix Phi, Lam, Z, Y, eta
    real colvector lam, tau
    real scalar nf, j, f, c

    nf  = p / 10
    Phi = J(nf, nf, 0.3) + 0.7 * I(nf)
    lam = J(p, 1, .)
    Lam = J(p, nf, 0)
    for (j = 1; j <= p; j++) {
        f         = ceil(j / 10)
        lam[j]    = 0.5 + 0.3 * mod(j - 1, 10) / 9
        Lam[j, f] = lam[j]
    }
    tau = invnormal((1::(K - 1)) / K)
    eta = rnormal(N, nf, 0, 1) * cholesky(Phi)'
    Z   = eta * Lam' + rnormal(N, p, 0, 1) :* sqrt(1 :- lam' :^ 2)
    Y   = J(N, p, 0)
    for (c = 1; c <= rows(tau); c++) {
        Y = Y + (Z :> tau[c])
    }
    return(Y)
}

real scalar tb_median(real colvector t)
{
    real colvector s
    real scalar n

    s = sort(t, 1)
    n = rows(s)
    return(mod(n, 2) ? s[(n + 1) / 2] : (s[n / 2] + s[n / 2 + 1]) / 2)
}

void tb_module2(real matrix Y, real colvector K)
{
    pointer(real colvector) colvector th
    real scalar a, b

    th = J(cols(Y), 1, NULL)
    for (a = 1; a <= cols(Y); a++) {
        th[a] = &(dwlsmv_thresholds(Y[., a], K[a]))
    }
    for (a = 1; a < cols(Y); a++) {
        for (b = a + 1; b <= cols(Y); b++) {
            (void) dwls_polychoric(Y[., a], Y[., b], *th[a], *th[b])
        }
    }
}
// ---- end of the copy ------------------------------------------------------

string scalar tc_gridfile(real scalar p, real scalar K, real scalar N)
{
    return(st_local("tdir") + "/grid_p" + strofreal(p) + "_K" + strofreal(K)
           + "_N" + strofreal(N) + ".csv")
}

// the data in memory, y1..yp, as a matrix
void tc_put(real matrix Y)
{
    real scalar j

    stata("clear")
    st_addobs(rows(Y))
    for (j = 1; j <= cols(Y); j++) {
        (void) st_addvar("byte", "y" + strofreal(j))
    }
    st_store(., ., Y)
}

real matrix tc_read(string scalar f)
{
    stata("quietly import delimited using " + char(34) + f + char(34)
          + ", delimiters(" + char(34) + "," + char(34) + ") varnames(1) clear")
    return(st_data(., .))
}

void tc_write(real matrix Y, string scalar f)
{
    tc_put(Y)
    stata("quietly export delimited using " + char(34) + f + char(34) + ", replace nolabel")
    if (tc_read(f) != Y) {
        errprintf("%s does not read back as the matrix written\n", f)
        _error(3498, "timing data did not round-trip")
    }
    printf("TIME|data|%s|%g x %g|written and read back\n", f, rows(Y), cols(Y))
}

real colvector tc_mck(string scalar c)
{
    return(c == "mc_a1p_16" ? J(16, 1, 3) : (4 \ 4 \ 3 \ 3))
}

void tc_writedata()
{
    real rowvector ps, Ks, Ns
    string rowvector mc
    string scalar f
    real scalar p, K, N, i

    rseed(20260926)
    ps = (10, 20, 40); Ks = (3, 5); Ns = (1000, 10000)
    for (p = 1; p <= cols(ps); p++) {
        for (K = 1; K <= cols(Ks); K++) {
            for (N = 1; N <= cols(Ns); N++) {
                tc_write(tb_sim(ps[p], Ks[K], Ns[N]), tc_gridfile(ps[p], Ks[K], Ns[N]))
            }
        }
    }
    mc = tokens(st_local("mcs"))
    for (i = 1; i <= cols(mc); i++) {
        f = "fixtures/replication/" + mc[i] + "/" + mc[i] + "_0001.dta"
        if (!fileexists(f)) {
            errprintf("absent: %s (gitignored; regenerate with 50_replicates.R)\n", f)
            _error(3498, "a Monte Carlo replicate file is absent")
        }
        stata("quietly use " + char(34) + f + char(34) + ", clear")
        tc_write(st_data(., .), st_local("tdir") + "/" + mc[i] + ".csv")
    }
}

void tc_one(string scalar tag, real matrix Y, real colvector K, string scalar spec,
            real scalar reps)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    real colvector t1, t2, tm, tt
    real scalar r

    // warm-up, untimed
    st  = dwls_sample_stats(Y, K)
    par = dwls_identify(spec, st.map)
    fit = dwls_fit(st, par)
    dwls_inference(fit)

    t1 = J(reps, 1, .); t2 = J(reps, 1, .); tm = J(reps, 1, .); tt = J(reps, 1, .)
    for (r = 1; r <= reps; r++) {
        timer_clear(84); timer_on(84)
        timer_clear(81); timer_on(81)
        st = dwls_sample_stats(Y, K)
        timer_off(81); t1[r] = timer_value(81)[1]

        timer_clear(82); timer_on(82)
        par = dwls_identify(spec, st.map)
        fit = dwls_fit(st, par)
        dwls_inference(fit)
        timer_off(82); t2[r] = timer_value(82)[1]
        timer_off(84); tt[r] = timer_value(84)[1]

        timer_clear(83); timer_on(83)
        tb_module2(Y, K)
        timer_off(83); tm[r] = timer_value(83)[1]
    }
    printf("TIME|dwlsmv|%s|p %g K %g N %g|D %g npar %g|reps %g|total %9.4f s|stage1 %9.4f s|module2 %9.4f s (%5.1f%% of stage1)|stage2 %9.4f s|converged %g|iterations %g\n",
           tag, cols(Y), max(K), rows(Y), st.map.D, par.npar, reps,
           tb_median(tt), tb_median(t1), tb_median(tm),
           100 * tb_median(tm) / tb_median(t1), tb_median(t2), fit.converged, fit.iterations)
    displayflush()
}

void tc_main()
{
    real rowvector ps, Ks, Ns
    string rowvector mc
    string scalar f
    real matrix Y
    real scalar p, K, N, i, reps, smoke

    smoke = st_local("mode") == "smoke"
    ps = (10, 20, 40); Ks = (3, 5); Ns = (1000, 10000)
    if (smoke) {
        ps = 10; Ks = 3; Ns = 1000
    }
    for (p = 1; p <= cols(ps); p++) {
        for (K = 1; K <= cols(Ks); K++) {
            for (N = 1; N <= cols(Ns); N++) {
                f = tc_gridfile(ps[p], Ks[K], Ns[N])
                if (!fileexists(f)) {
                    errprintf("absent: %s; run this file with writedata first\n", f)
                    _error(3498, "timing data absent")
                }
                Y    = tc_read(f)
                reps = smoke ? 1 : ((ps[p] == 40 & Ns[N] == 10000) ? 3 : 5)
                tc_one("grid", Y, J(ps[p], 1, Ks[K]), tb_spec(ps[p], ps[p] / 10, 0), reps)
            }
        }
    }
    mc = tokens(st_local("mcs"))
    for (i = 1; i <= cols(mc); i++) {
        f = st_local("tdir") + "/" + mc[i] + ".csv"
        if (!fileexists(f)) {
            errprintf("absent: %s; run this file with writedata first\n", f)
            _error(3498, "timing data absent")
        }
        Y = tc_read(f)
        tc_one(mc[i], Y, tc_mck(mc[i]), tb_spec(cols(Y), 1, 1), smoke ? 1 : 5)
    }
}

end

if "`mode'" == "writedata" {
    capture mkdir "fixtures/replication"
    capture mkdir "`tdir'"
    mata: tc_writedata()
    exit
}

mata: tc_main()

* nogsem (added 2026-10-05): dwlsmv rows only, for a re-time after src/mata
* changes, compared with an earlier run's gsem rows (code unchanged there)
if "`mode'" == "nogsem" {
    exit
}

* ---- gsem, ordinal probit, total only: N = 1000, 1 / 2 / 4 factors --------
* Stata-level loop: gsem is a Stata command. Local macros here are ordinary
* Stata macros; no Mata string below contains a quoted word.
local plist 10 20 40
local klist 3 5
if "`mode'" == "smoke" {
    local plist 10
    local klist 3
}
local warm 0
foreach p of local plist {
    foreach K of local klist {
        quietly import delimited using "`tdir'/grid_p`p'_K`K'_N1000.csv", delimiters(",") varnames(1) clear
        local nf = `p' / 10
        local paths
        local vars
        forvalues f = 1/`nf' {
            local a = (`f' - 1) * 10 + 1
            local b = `f' * 10
            local paths `paths' (F`f' -> y`a'-y`b')
            local vars `vars' F`f'@1
        }
        local gcmd gsem `paths', oprobit var(`vars') intmethod(mvaghermite) intpoints(`gpts') iterate(`giter')
        if !`warm' {
            quietly `gcmd'
            local warm 1
        }
        local reps = cond("`mode'" == "smoke", 1, 5)
        local ts
        local ran 0
        forvalues r = 1/`reps' {
            timer clear 91
            timer on 91
            capture quietly `gcmd'
            local rc = _rc
            timer off 91
            quietly timer list 91
            local ts `ts' `r(t91)'
            local ++ran
            if `rc' {
                display "TIME|gsem|grid|p `p' K `K' N 1000|rep `r' failed rc `rc'"
                continue, break
            }
            if `r(t91)' > `gcap' {
                continue, break
            }
        }
        if `rc' {
            continue
        }
        * y1's loading must be free: a positive, nonmissing standard error
        capture assert _se[y1:F1] > 0 & !missing(_se[y1:F1])
        local free1 = (_rc == 0)
        mata: st_local("med", strofreal(tb_median(strtoreal(tokens(st_local("ts")))'), "%9.4f"))
        display "TIME|gsem|grid|p `p' K `K' N 1000|factors `nf'|rank `e(rank)'|reps `ran' of `reps'|total `med' s|intmethod `e(intmethod)'|intpoints `e(n_quad)'|converged `e(converged)'|iterations `e(ic)'|iterate `giter'|cap `gcap' s|y1 loading free `free1'"
    }
}
