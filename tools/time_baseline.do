*! dwlsmv -- spec 7.5's timing, Stata side only: the pre-Module-2-fix baseline
*!
*! Executed under StataMP 18 on 2026-09-26 at e28b5eb, on the idle 6-core
*! development machine. Results: E-089.
*!
*! A diagnostic, not a test: no pass/fail, and the suite runner ignores tools/.
*! Decided with the user 2026-09-26: a Stata-only baseline now, the full
*! comparison against lavaan and Mplus after the Module 2 fix, so the pair
*! shows what the fix cost and feeds spec 6's Module 2 row (profile before
*! switching optimize() to an analytic Newton-Raphson).
*!
*!     cd <repo root>
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/time_baseline.do
*!     grep "^TIME|" time_baseline.log
*!
*! RUN IT ON AN IDLE MACHINE. Wall-clock timing moves with load; a result taken
*! beside another Stata or R job is not a baseline.
*!
*! WHAT IS TIMED (spec 7.5). Stage 1, dwls_sample_stats(); Stage 2,
*! dwls_identify() through dwls_inference(); and, inside Stage 1's time,
*! Module 2 alone: every pair's dwls_polychoric() at the sample thresholds,
*! which is where Stage 1's cost is known to sit. Each is the MEDIAN over
*! `reps' timed repetitions, after one untimed warm-up call, so one-time Mata
*! loading and first-call overhead are excluded.
*!
*! DESIGNS. A grid over p in {10, 20, 40}, K in {3, 5}, N in {1000, 10000},
*! simulated as tools/time_fit.do does (one factor per ten items, factor
*! correlations .3, loadings .5 to .8, fixed seed), fitted with the true
*! structure under unit-variance metric. It includes the quadcross() benchmark
*! size, p 40, K 5, N 10000. Plus replicate 1 of three Monte Carlo
*! conditions (mc_d, mc_a1p_16, mc_b1_m1) under the one-factor marker model
*! the paired run fits, when their gitignored .dta files are present.
*!
*! FAIRNESS RECORD (spec 7.5): Stata version and flavour, processors in use,
*! processors licensed and on the machine, host. Every line of output starts
*! TIME| so the performance note can be built from a grep.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
args mode
clear all
set more off
set linesize 255

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

// Module 2 alone: every pair's polychoric at the sample thresholds.
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

void tb_one(string scalar tag, real matrix Y, real colvector K, string scalar spec,
            real scalar reps)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    real colvector t1, t2, tm
    real scalar r

    // warm-up, untimed
    st  = dwls_sample_stats(Y, K)
    par = dwls_identify(spec, st.map)
    fit = dwls_fit(st, par)
    dwls_inference(fit)

    t1 = J(reps, 1, .); t2 = J(reps, 1, .); tm = J(reps, 1, .)
    for (r = 1; r <= reps; r++) {
        timer_clear(81); timer_on(81)
        st = dwls_sample_stats(Y, K)
        timer_off(81); t1[r] = timer_value(81)[1]

        timer_clear(82); timer_on(82)
        par = dwls_identify(spec, st.map)
        fit = dwls_fit(st, par)
        dwls_inference(fit)
        timer_off(82); t2[r] = timer_value(82)[1]

        timer_clear(83); timer_on(83)
        tb_module2(Y, K)
        timer_off(83); tm[r] = timer_value(83)[1]
    }
    printf("TIME|%s|p %g K %g N %g|D %g npar %g|reps %g|stage1 %9.4f s|module2 %9.4f s (%5.1f%% of stage1)|stage2 %9.4f s|converged %g\n",
           tag, cols(Y), max(K), rows(Y), st.map.D, par.npar, reps,
           tb_median(t1), tb_median(tm), 100 * tb_median(tm) / tb_median(t1),
           tb_median(t2), fit.converged)
    displayflush()
}

void tb_main()
{
    real matrix Y
    real scalar p, K, N, reps
    real rowvector ps, Ks, Ns
    string rowvector mc
    string scalar mcf
    real scalar smoke
    real colvector mcK
    real scalar i

    smoke = st_local("mode") == "smoke"
    rseed(20260926)
    ps = (10, 20, 40); Ks = (3, 5); Ns = (1000, 10000)
    if (smoke) {
        ps = 10; Ks = 3; Ns = 1000
    }
    for (p = 1; p <= cols(ps); p++) {
        for (K = 1; K <= cols(Ks); K++) {
            for (N = 1; N <= cols(Ns); N++) {
                Y    = tb_sim(ps[p], Ks[K], Ns[N])
                reps = smoke ? 1 : ((ps[p] == 40 & Ns[N] == 10000) ? 3 : 5)
                tb_one("grid", Y, J(ps[p], 1, Ks[K]), tb_spec(ps[p], ps[p] / 10, 0), reps)
            }
        }
    }

    mc  = ("mc_d", "mc_a1p_16", "mc_b1_m1")
    for (i = 1; i <= cols(mc); i++) {
        mcf = "fixtures/replication/" + mc[i] + "/" + mc[i] + "_0001.dta"
        if (!fileexists(mcf)) {
            printf("TIME|%s|absent: %s (gitignored; regenerate with 50_replicates.R)\n", mc[i], mcf)
            continue
        }
        stata("quietly use " + char(34) + mcf + char(34) + ", clear")
        Y   = st_data(., .)
        mcK = (mc[i] == "mc_a1p_16") ? J(16, 1, 3) : (4 \ 4 \ 3 \ 3)
        tb_one(mc[i], Y, mcK, tb_spec(cols(Y), 1, 1), smoke ? 1 : 5)
    }
}

tb_main()

end
