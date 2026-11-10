*! dwlsmv -- spec 7.5's timing at G > 1, Stata side (E-089)
*!
*! Executed 2026-11-02 (E-089). A diagnostic, not a test: no pass/fail, and the runner ignores tools/.
*!
*! MODES (run from a repo copy, on an idle machine for time):
*!   writedata  simulate the 24 datasets (seed 20261103), write them to
*!              fixtures/replication/timing/mg_G<G>_p<p>_K<K>_N<Ng>.csv, read each back.
*!   match      one untimed fit per dataset and level: writes match_dwlsmv.csv and
*!              match_theta.csv there, for tools/time_mg_lavaan.R and time_mg_mplus.R.
*!   time       the two dwlsmv rows, median of 5 after one warm-up:
*!              estimator total  Stage 1 per group, then dwlsmv_levelgen through
*!                               dwls_inference (the API, no Module 10, no baseline, no
*!                               fit indices, no display), with Module 2 beside it;
*!              command total    quietly dwlsmv ..., group() invariance() metric(unitvariance).
*!   smoke      as time, on the smallest dataset only.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
args mode
clear all
set more off
set linesize 255

local tdir "fixtures/replication/timing"
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

display "TIME|machine|stata `c(stata_version)' MP=`c(MP)'|processors `c(processors)'|`c(current_date)' `c(current_time)'"

mata:

string scalar tm_file(real scalar G, real scalar p, real scalar K, real scalar N)
{
    return(st_local("tdir") + "/mg_G" + strofreal(G) + "_p" + strofreal(p) + "_K" +
           strofreal(K) + "_N" + strofreal(N) + ".csv")
}

// One factor per ten items, loadings .5 to .8, factor correlations .3. Later groups:
// factor means .3, variances 1.2. Loadings, thresholds and residual variances are as in
// group 1, so the scalar model holds in the population.
real matrix tm_sim(real scalar G, real scalar p, real scalar K, real scalar Ng)
{
    real matrix Phi, Lam, Z, Y, eta, out
    real colvector lam, tau
    real scalar nf, j, f, c, g

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
    out = J(0, p + 1, .)
    for (g = 1; g <= G; g++) {
        eta = rnormal(Ng, nf, 0, 1) * cholesky(Phi)'
        if (g > 1) {
            eta = sqrt(1.2) * eta :+ 0.3
        }
        Z = eta * Lam' + rnormal(Ng, p, 0, 1) :* sqrt(1 :- lam' :^ 2)
        Y = J(Ng, p, 0)
        for (c = 1; c <= rows(tau); c++) {
            Y = Y + (Z :> tau[c])
        }
        for (j = 1; j <= p; j++) {
            if (rows(uniqrows(Y[., j])) != K) {
                _error(3498, "tm_sim: a category is empty in a group; the design needs review")
            }
        }
        out = out \ (J(Ng, 1, g), Y)
    }
    return(out)
}

void tm_put(real matrix D)
{
    real scalar j

    stata("clear")
    st_addobs(rows(D))
    (void) st_addvar("byte", "g")
    for (j = 2; j <= cols(D); j++) {
        (void) st_addvar("byte", "y" + strofreal(j - 1))
    }
    st_store(., ., D)
}

real matrix tm_read(string scalar f)
{
    stata("quietly import delimited using " + char(34) + f + char(34) +
          ", delimiters(" + char(34) + "," + char(34) + ") varnames(1) clear")
    return(st_data(., .))
}

string scalar tm_conf(real scalar G, real scalar p)
{
    string scalar s, nl, fl, il
    real scalar g, j, f, h, nf

    nl = char(10)
    nf = p / 10
    fl = "factor"
    for (f = 1; f <= nf; f++) {
        fl = fl + " F" + strofreal(f)
    }
    il = "ind"
    for (j = 1; j <= p; j++) {
        il = il + " y" + strofreal(j)
    }
    s = "groups " + strofreal(G) + nl
    for (g = 1; g <= G; g++) {
        s = s + "group " + strofreal(g) + nl + fl + nl + il + nl
        for (j = 1; j <= p; j++) {
            s = s + "lambda F" + strofreal(ceil(j / 10)) + " y" + strofreal(j) + " free" + nl
        }
        for (f = 1; f <= nf; f++) {
            for (h = f + 1; h <= nf; h++) {
                s = s + "psi F" + strofreal(f) + " F" + strofreal(h) + " free" + nl
            }
        }
        for (f = 1; f <= nf; f++) {
            s = s + (g == 1 ? "metric unitvariance F" + strofreal(f) :
                              "psi F" + strofreal(f) + " F" + strofreal(f) + " fixed 1") + nl
        }
    }
    return(s)
}

string scalar tm_cmd(real scalar p, string scalar level)
{
    string scalar s
    real scalar f

    s = "quietly dwlsmv"
    for (f = 1; f <= p / 10; f++) {
        s = s + " (F" + strofreal(f) + " -> y" + strofreal((f - 1) * 10 + 1) + "-y" + strofreal(f * 10) + ")"
    }
    return(s + ", group(g) invariance(" + level + ") metric(unitvariance)")
}

real scalar tm_median(real colvector t)
{
    real colvector s
    real scalar n

    s = sort(t, 1)
    n = rows(s)
    return(mod(n, 2) ? s[(n + 1) / 2] : (s[n / 2] + s[n / 2 + 1]) / 2)
}

void tm_module2(real matrix Y, real colvector K)
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

struct dwlsmv_stage1 vector tm_stage1(real matrix D, real scalar G, real scalar K)
{
    struct dwlsmv_stage1 vector sts
    real scalar g

    sts = dwlsmv_stage1(G, 1)
    for (g = 1; g <= G; g++) {
        sts[g] = dwls_sample_stats(select(D[., 2..cols(D)], D[., 1] :== g), J(cols(D) - 1, 1, K))
    }
    return(sts)
}

struct dwlsmv_fit scalar tm_stage2(struct dwlsmv_stage1 vector sts, string scalar conf,
                                   real matrix KM, string scalar level)
{
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    real scalar g

    maps = dwlsmv_statmap(rows(sts), 1)
    for (g = 1; g <= rows(sts); g++) {
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(conf, KM, level), maps, sts)
    fit = dwls_fit(sts, par)
    dwls_inference(fit)
    return(fit)
}

// retry flag: the same fit with no retry, untimed; 1 if its first run failed
real scalar tm_retried(struct dwlsmv_stage1 vector sts, string scalar conf, real matrix KM,
                       string scalar level)
{
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    real scalar g

    maps = dwlsmv_statmap(rows(sts), 1)
    for (g = 1; g <= rows(sts); g++) {
        maps[g] = sts[g].map
    }
    par = dwls_identify(dwlsmv_levelgen(conf, KM, level), maps, sts)
    fit = dwlsmv_fit_noretry(sts, par)
    return(!fit.converged)
}

void tm_writedata()
{
    real rowvector Gs, ps, Ks, Ns
    real scalar a, b, c, d
    real matrix D
    string scalar f

    Gs = (2, 3); ps = (10, 20); Ks = (3, 5); Ns = (500, 5000)
    for (a = 1; a <= 2; a++) for (b = 1; b <= 2; b++) for (c = 1; c <= 2; c++) for (d = 1; d <= 2; d++) {
        D = tm_sim(Gs[a], ps[b], Ks[c], Ns[d])
        f = tm_file(Gs[a], ps[b], Ks[c], Ns[d])
        tm_put(D)
        stata("quietly export delimited using " + char(34) + f + char(34) + ", replace nolabel")
        if (tm_read(f) != D) {
            _error(3498, "timing data did not round-trip")
        }
        printf("TIME|data|%s|%g x %g|written and read back\n", f, rows(D), cols(D))
    }
}

void tm_match()
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_fit scalar fit
    real rowvector Gs, ps, Ks, Ns
    real scalar a, b, c, d, l, fh, ft, j, re, N
    real matrix D, KM
    string rowvector levs
    string scalar conf, cell

    Gs = (2, 3); ps = (10, 20); Ks = (3, 5); Ns = (500, 5000)
    levs = ("configural", "scalar")
    unlink(st_local("tdir") + "/match_dwlsmv.csv")
    unlink(st_local("tdir") + "/match_theta.csv")
    fh = fopen(st_local("tdir") + "/match_dwlsmv.csv", "w")
    ft = fopen(st_local("tdir") + "/match_theta.csv", "w")
    fput(fh, "cell,level,npar,df,T3,Tmplus,converged,retried")
    fput(ft, "cell,level,abs_theta")
    for (a = 1; a <= 2; a++) for (b = 1; b <= 2; b++) for (c = 1; c <= 2; c++) for (d = 1; d <= 2; d++) {
        D    = tm_read(tm_file(Gs[a], ps[b], Ks[c], Ns[d]))
        sts  = tm_stage1(D, Gs[a], Ks[c])
        conf = tm_conf(Gs[a], ps[b])
        KM   = J(Gs[a], ps[b], Ks[c])
        cell = "G" + strofreal(Gs[a]) + "_p" + strofreal(ps[b]) + "_K" + strofreal(Ks[c]) + "_N" + strofreal(Ns[d])
        N    = rows(D)
        for (l = 1; l <= 2; l++) {
            fit = tm_stage2(sts, conf, KM, levs[l])
            re  = tm_retried(sts, conf, KM, levs[l])
            fput(fh, cell + "," + levs[l] + "," + strofreal(rows(fit.theta)) + "," + strofreal(fit.T3_df) + "," +
                 strofreal(fit.T3, "%21.0g") + "," +
                 strofreal(fit.shift_a * fit.T_raw * N / (N - Gs[a]) + fit.shift_b, "%21.0g") + "," +
                 strofreal(fit.converged) + "," + strofreal(re))
            for (j = 1; j <= rows(fit.theta); j++) {
                fput(ft, cell + "," + levs[l] + "," + strofreal(abs(fit.theta[j]), "%21.0g"))
            }
            printf("TIME|match|dwlsmv|%s|%s|npar %g|df %g|T3 %12.6f|converged %g|retried %g\n",
                   cell, levs[l], rows(fit.theta), fit.T3_df, fit.T3, fit.converged, re)
            displayflush()
        }
    }
    fclose(fh)
    fclose(ft)
}

void tm_time(real scalar smoke)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_fit scalar fit
    real rowvector Gs, ps, Ks, Ns
    real scalar a, b, c, d, l, r, g, R
    real matrix D, KM
    real colvector t1, t2, tt, tm, tc
    string rowvector levs
    string scalar conf, cell, cmd, dv

    Gs = (2, 3); ps = (10, 20); Ks = (3, 5); Ns = (500, 5000)
    levs = ("configural", "scalar")
    R = 5
    for (a = 1; a <= 2; a++) for (b = 1; b <= 2; b++) for (c = 1; c <= 2; c++) for (d = 1; d <= 2; d++) {
        if (smoke & (a + b + c + d > 4)) {
            continue
        }
        D    = tm_read(tm_file(Gs[a], ps[b], Ks[c], Ns[d]))
        conf = tm_conf(Gs[a], ps[b])
        KM   = J(Gs[a], ps[b], Ks[c])
        cell = "G" + strofreal(Gs[a]) + "_p" + strofreal(ps[b]) + "_K" + strofreal(Ks[c]) + "_N" + strofreal(Ns[d])
        tm = J(R, 1, .)
        for (r = 1; r <= R; r++) {
            timer_clear(4)
            timer_on(4)
            for (g = 1; g <= Gs[a]; g++) {
                tm_module2(select(D[., 2..cols(D)], D[., 1] :== g), J(cols(D) - 1, 1, Ks[c]))
            }
            timer_off(4)
            tm[r] = timer_value(4)[1]
        }
        for (l = 1; l <= 2; l++) {
            sts = tm_stage1(D, Gs[a], Ks[c])
            fit = tm_stage2(sts, conf, KM, levs[l])
            t1 = J(R, 1, .); t2 = J(R, 1, .); tt = J(R, 1, .)
            for (r = 1; r <= R; r++) {
                timer_clear(1); timer_clear(2); timer_clear(3)
                timer_on(3)
                timer_on(1)
                sts = tm_stage1(D, Gs[a], Ks[c])
                timer_off(1)
                timer_on(2)
                fit = tm_stage2(sts, conf, KM, levs[l])
                timer_off(2)
                timer_off(3)
                t1[r] = timer_value(1)[1]; t2[r] = timer_value(2)[1]; tt[r] = timer_value(3)[1]
            }
            printf("TIME|dwlsmv|estimator|%s|%s|npar %g|reps %g|stage1 %9.4f|module2 %9.4f|stage2 %9.4f|total %9.4f s|converged %g|iterations %g\n",
                   cell, levs[l], rows(fit.theta), R, tm_median(t1), tm_median(tm), tm_median(t2), tm_median(tt),
                   fit.converged, fit.iterations)
            displayflush()

            tm_put(D)
            cmd = tm_cmd(ps[b], levs[l])
            stata(cmd)
            tc = J(R, 1, .)
            for (r = 1; r <= R; r++) {
                timer_clear(5)
                timer_on(5)
                stata(cmd)
                timer_off(5)
                tc[r] = timer_value(5)[1]
            }
            dv = strofreal(st_numscalar("e(T3)"), "%12.6f")
            printf("TIME|dwlsmv|command|%s|%s|reps %g|total %9.4f s|T3 %s\n", cell, levs[l], R, tm_median(tc), dv)
            displayflush()
        }
    }
}
end

if "`mode'" == "writedata" {
    set seed 20261103
    mata: tm_writedata()
}
else if "`mode'" == "match" {
    mata: tm_match()
}
else if "`mode'" == "time" {
    mata: tm_time(0)
}
else if "`mode'" == "smoke" {
    mata: tm_time(1)
}
else {
    display as error "mode must be writedata, match, time or smoke"
    exit 198
}
display "TIME|done"
