*! dwlsmv -- the production paired run over the Monte Carlo replicates
*!
*! Executed under StataMP 18 on 2026-09-26 at f1c95c6, all 42 conditions in four
*! shards (about 15 minutes each), after a 5-replicate smoke test on mc_d,
*! mc_b1_m1 and mc_c_999. Results: E-071, 2026-09-26.
*!
*! A diagnostic tool, not a pass/fail test (tools/). Spec 7.1's paired design
*! on the Module 10 Monte Carlo reference: production is fitted to exactly
*! the replicate .dta files lavaan was fitted to by
*! fixtures/generate/70_lavaan_mc_reference.R, and one slim row per replicate
*! is written to fixtures/replication/<condition>/production_fits.csv
*! (gitignored, like lavaan_fits.csv beside it). tools/paired_mc_summary.R
*! compares the two. Decided with the user 2026-09-26: all 42 conditions in
*! four shards; the primary comparison is production's own z_var; Module 10's
*! rows are recorded per replicate.
*!
*! FAILURES ARE RECORDED, NEVER EXCLUDED. Every step runs under -capture-, and
*! a failure writes its row with status "error", the stage it died in
*! (stage1, identify, fit, inference, diag) and its rc. A replicate whose
*! Stage 1 aborts while countdiag() reads an empty declared category is
*! labelled "empty_category" -- Module 1's deliberate abort, the same class
*! lavaan did not fit -- but it is still ATTEMPTED, so the label is observed
*! and not assumed. A non-converged fit is "nonconverged" and keeps its
*! estimates and SEs. Nothing here depends on the open Module 2 near +/-1
*! defect: a stalled pair builds its Gamma anyway, and gammadiag's centring
*! row is recorded next to the SEs to show whether it catches it.
*!
*! THE MODEL is each condition's, as 70_ fits it: one factor, every indicator
*! loading, marker identification on y1 (lambda_y1 = 1), psi free, the
*! thresholds free. K comes from the committed lavaan_mc_conditions.csv.
*!
*! USAGE, from the repository root. Stata batch names its log after the
*! do-file invoked, so run each shard through its own one-line wrapper:
*!     wrapper paired_mc_s1.do:  cd <repo>  /  do tools/paired_mc.do 1 4
*!     StataMP-64.exe /e do paired_mc_s1.do        (log: paired_mc_s1.log)
*! Arguments: shard, nshards [, one condition [, max replicates [, keepbound]]].
*! keepbound (added 2026-10-05) runs Stage 1 with keep_bound = 1 and writes
*! production_fits_keepbound.csv instead. Conditions
*! go to shards round-robin in the committed file's order.
*!
*! Columns: condition, replicate, seed, status, stage, rc, N, converged,
*! iterations, T3, T3_df, then Module 10's rows as g<r>f / g<r>v (gammadiag
*! flag and value, r = 1..7), g_rstar, c1v, c2v (countdiag), ps1..ps8 and
*! ps2_empty (Module 2 pair statuses, from 2026-10-01), b<r>f / b<r>v
*! (breaddiag, r = 1..4), then est:<name> and se:<name> in lavaan's names
*! (f=~y2, y1|t1, f~~f). SE is sqrt(diag(e(V))), e(V) = ACOV / (N - 1)
*! (dwlsmv_ev). Numbers are written %21.17g; a missing value as Stata writes
*! it (".", ".n").
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
args shard nshards onecond nmax kb
if "`shard'" == "" | "`nshards'" == "" {
    display as error "usage: do tools/paired_mc.do shard nshards [condition [max [keepbound]]]"
    exit 198
}
if "`nmax'" == "" local nmax = .
local kbsfx = cond("`kb'" == "keepbound", "_keepbound", "")

clear all
set more off
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
quietly do src/mata/dwlsmv_types.mata
foreach f of local srcs {
    quietly do "src/mata/`f'"
}

mata:

// the canonical specification of the one-factor marker model on p items
string scalar pm_spec(real scalar p)
{
    string scalar s, ind
    real scalar j

    ind = "ind"
    for (j = 1; j <= p; j++) {
        ind = ind + " y" + strofreal(j)
    }
    s = "factor f" + char(10) + ind + char(10)
    for (j = 2; j <= p; j++) {
        s = s + "lambda f y" + strofreal(j) + " free" + char(10)
    }
    return(s + "psi f f free" + char(10) + "metric marker f y1" + char(10))
}

// par.name in lavaan's spelling: thresholds already are; lambda[yj,f] is
// f=~yj; psi[f,f] is f~~f.
string scalar pm_lavname(string scalar nm)
{
    if (substr(nm, 1, 7) == "lambda[") {
        return("f=~" + substr(nm, 8, strpos(nm, ",") - 8))
    }
    if (substr(nm, 1, 4) == "psi[") {
        return("f~~f")
    }
    return(nm)
}

// K from lavaan_mc_conditions.csv's thresholds field: variables separated
// by "/", each variable's thresholds by ";". K_j = thresholds + 1.
real colvector pm_kof(string scalar thr)
{
    string rowvector t, parts
    real colvector K
    real scalar j

    t = tokens(thr, "/")
    K = J(0, 1, .)
    for (j = 1; j <= cols(t); j++) {
        if (t[j] != "/") {
            parts = tokens(t[j], ";")
            K = K \ (sum(parts :!= ";") + 1)
        }
    }
    return(K)
}

string scalar pm_num(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

// The row buffer: fields set as they become known, so a failure part-way
// still writes everything that was computed before it.
void pm_begin(string scalar cond, real scalar rep, real scalar seed, real scalar p)
{
    external string rowvector PM_ROW, PM_HEAD
    external string scalar PM_STAGE
    external real scalar PM_P
    external real colvector PM_K
    external real scalar PM_KB
    string rowvector nm
    real scalar j, c

    PM_P = p
    // parameter names in theta order: thresholds, lambdas y2..yp, psi
    nm = J(1, 0, "")
    for (j = 1; j <= p; j++) {
        for (c = 1; c <= PM_K[j] - 1; c++) {
            nm = nm, "y" + strofreal(j) + "|t" + strofreal(c)
        }
    }
    for (j = 2; j <= p; j++) {
        nm = nm, "f=~y" + strofreal(j)
    }
    nm = nm, "f~~f"
    PM_HEAD = ("condition", "replicate", "seed", "status", "stage", "rc", "N", "converged",
               "iterations", "T3", "T3_df")
    for (j = 1; j <= 7; j++) PM_HEAD = PM_HEAD, "g" + strofreal(j) + "f", "g" + strofreal(j) + "v"
    PM_HEAD = PM_HEAD, "g_rstar", "c1v", "c2v"
    for (j = 1; j <= 8; j++) PM_HEAD = PM_HEAD, "ps" + strofreal(j)
    PM_HEAD = PM_HEAD, "ps2_empty"
    for (j = 1; j <= 4; j++) PM_HEAD = PM_HEAD, "b" + strofreal(j) + "f", "b" + strofreal(j) + "v"
    PM_HEAD = PM_HEAD, ("est:" :+ nm), ("se:" :+ nm)
    PM_ROW = J(1, cols(PM_HEAD), ".")
    PM_ROW[1..3] = (cond, strofreal(rep), strofreal(seed, "%12.0f"))
    PM_ROW[4] = "ok"
    PM_ROW[5] = ""
    PM_ROW[6] = "0"
    PM_STAGE = "load"
}

real scalar pm_col(string scalar h)
{
    external string rowvector PM_HEAD
    real scalar j

    for (j = 1; j <= cols(PM_HEAD); j++) {
        if (PM_HEAD[j] == h) return(j)
    }
    _error("pm_col: no column " + h)
}

void pm_one()
{
    external string rowvector PM_ROW
    external string scalar PM_STAGE
    external real colvector PM_K
    external real scalar PM_KB
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    struct dwlsmv_pairfit scalar f
    real matrix Y, G, C, B, eV
    real colvector se
    real rowvector ps
    real scalar i, j, k, pos, vi, vj, ps2e

    Y = st_data(., .)
    PM_ROW[pm_col("N")] = strofreal(rows(Y))

    PM_STAGE = "countdiag"
    C = dwlsmv_countdiag(Y, PM_K)
    PM_ROW[pm_col("c1v")] = pm_num(C[1, 2])
    PM_ROW[pm_col("c2v")] = pm_num(C[2, 2])

    PM_STAGE = "stage1"
    st = dwls_sample_stats(Y, PM_K, PM_KB)

    PM_STAGE = "diag"
    G = dwlsmv_gammadiag(st, Y)
    for (j = 1; j <= 7; j++) {
        PM_ROW[pm_col("g" + strofreal(j) + "f")] = pm_num(G[j, 1])
        PM_ROW[pm_col("g" + strofreal(j) + "v")] = pm_num(G[j, 2])
    }
    PM_ROW[pm_col("g_rstar")] = pm_num(G[7, 3])

    // Module 2's per-pair statuses (added 2026-10-01, option B). The stage1
    // struct does not carry them, so each pair is refitted exactly as
    // dwls_sample_stats() fits it (same thresholds, sliced from s; the fit is
    // deterministic). ps2_empty: status-2 pairs whose table has an empty
    // cell -- for a 2x2 table that means no interior root; for K > 2 it is
    // only a count, not an explanation (E-056, 2026-09-30).
    PM_STAGE = "pairs"
    ps = J(1, 8, 0)
    ps2e = 0
    for (k = 1; k <= st.map.P; k++) {
        pos = st.map.T + k
        vi  = st.map.var_i[pos]
        vj  = st.map.var_j[pos]
        f = dwls_polychoric(Y[., vi], Y[., vj],
                st.s[|st.map.Coff[vi] + 1 \ st.map.Coff[vi] + st.map.C[vi]|],
                st.s[|st.map.Coff[vj] + 1 \ st.map.Coff[vj] + st.map.C[vj]|])
        if (f.rho != st.s[pos]) _error("paired_mc: a refitted pair differs from Stage 1's")
        ps[f.status] = ps[f.status] + 1
        if (f.status == 2) {
            if (any(dwlsmv_crosstab(Y[., vi], Y[., vj], PM_K[vi], PM_K[vj]) :== 0)) ps2e++
        }
    }
    for (j = 1; j <= 8; j++) PM_ROW[pm_col("ps" + strofreal(j))] = pm_num(ps[j])
    PM_ROW[pm_col("ps2_empty")] = pm_num(ps2e)

    PM_STAGE = "identify"
    par = dwls_identify(pm_spec(rows(PM_K)), st.map)

    PM_STAGE = "fit"
    fit = dwls_fit(st, par)
    PM_ROW[pm_col("converged")]  = pm_num(fit.converged)
    PM_ROW[pm_col("iterations")] = pm_num(fit.iterations)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            PM_ROW[pm_col("est:" + pm_lavname(par.name[i]))] = pm_num(fit.theta[par.fidx[i]])
        }
    }

    PM_STAGE = "diag"
    B = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    for (j = 1; j <= 4; j++) {
        PM_ROW[pm_col("b" + strofreal(j) + "f")] = pm_num(B[j, 1])
        PM_ROW[pm_col("b" + strofreal(j) + "v")] = pm_num(B[j, 2])
    }

    PM_STAGE = "inference"
    dwls_inference(fit)
    PM_ROW[pm_col("T3")]    = pm_num(fit.T3)
    PM_ROW[pm_col("T3_df")] = pm_num(fit.T3_df)
    eV = dwlsmv_ev(fit.acov, st.N)
    se = sqrt(diagonal(eV))
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            PM_ROW[pm_col("se:" + pm_lavname(par.name[i]))] = pm_num(se[par.fidx[i]])
        }
    }
    PM_ROW[4] = fit.converged == 1 ? "ok" : "nonconverged"
    PM_STAGE = "done"
}

// After a failure: which status, from the stage and countdiag's reading.
void pm_fail(real scalar rc)
{
    external string rowvector PM_ROW
    external string scalar PM_STAGE

    PM_ROW[5] = PM_STAGE
    PM_ROW[6] = strofreal(rc)
    if (PM_STAGE == "stage1" & PM_ROW[pm_col("c1v")] == "0") PM_ROW[4] = "empty_category"
    else if (PM_STAGE == "stage1" & PM_ROW[pm_col("c2v")] == "0") PM_ROW[4] = "empty_category"
    else PM_ROW[4] = "error"
}

void pm_open(string scalar f)
{
    external real scalar PM_FH
    external string rowvector PM_HEAD

    if (fileexists(f)) unlink(f)
    PM_FH = fopen(f, "w")
}

void pm_write()
{
    external real scalar PM_FH, PM_NROW
    external string rowvector PM_ROW, PM_HEAD
    string scalar s
    real scalar j

    if (PM_NROW == 0) {
        s = char(34) + PM_HEAD[1] + char(34)
        for (j = 2; j <= cols(PM_HEAD); j++) s = s + "," + char(34) + PM_HEAD[j] + char(34)
        fput(PM_FH, s)
    }
    s = PM_ROW[1]
    for (j = 2; j <= cols(PM_ROW); j++) s = s + "," + PM_ROW[j]
    fput(PM_FH, s)
    PM_NROW = PM_NROW + 1
}

void pm_close()
{
    external real scalar PM_FH
    fclose(PM_FH)
}

end

* ---- the condition list, from the committed tier-2 file --------------------------
import delimited using "fixtures/expected/lavaan_mc_conditions.csv", clear ///
    varnames(1) delimiters(",") stringcols(_all)
quietly keep if arm == "ref"
local nc = _N
local todo
forvalues i = 1/`nc' {
    local c = condition[`i']
    if "`onecond'" != "" {
        if "`c'" == "`onecond'" local todo `todo' `i'
    }
    else if mod(`i' - 1, `nshards') + 1 == `shard' {
        local todo `todo' `i'
    }
}
tempfile conds
quietly save `conds'

foreach i of local todo {
    use `conds', clear
    local c   = condition[`i']
    local thr = thresholds[`i']
    * K from the thresholds string: variables split by "/", thresholds by ";"
    mata: PM_K = pm_kof(st_local("thr"))
    mata: PM_KB = ("`kb'" == "keepbound")
    import delimited using "fixtures/replication/`c'/index.csv", clear ///
        varnames(1) delimiters(",") stringcols(_all)
    local nr = min(_N, `nmax')
    tempfile idx
    quietly save `idx'
    mata: PM_NROW = 0
    mata: pm_open("fixtures/replication/`c'/production_fits`kbsfx'.csv")
    timer clear 1
    timer on 1
    local nok 0
    local nnc 0
    local nec 0
    local ner 0
    forvalues r = 1/`nr' {
        use `idx', clear
        local rep  = replicate[`r']
        local seed = seed[`r']
        local file = file[`r']
        mata: pm_begin("`c'", `rep', `seed', rows(PM_K))
        capture {
            quietly use "`file'", clear
            mata: pm_one()
        }
        local rc = _rc
        if `rc' {
            mata: pm_fail(`rc')
        }
        mata: pm_write()
        mata: st_local("st", PM_ROW[4])
        if "`st'" == "ok" local ++nok
        else if "`st'" == "nonconverged" local ++nnc
        else if "`st'" == "empty_category" local ++nec
        else local ++ner
    }
    mata: pm_close()
    timer off 1
    quietly timer list 1
    display as result "paired_mc: `c'  ok `nok'  nonconverged `nnc'  empty_category `nec'  error `ner'  (" %6.1f r(t1) " s)"
}
display as result "paired_mc: shard `shard' of `nshards' done"
