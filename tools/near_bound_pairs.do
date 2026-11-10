*! dwlsmv -- per-pair record of the paired Monte Carlo, for the near-bound note
*!
*! Executed under StataMP 18 on 2026-10-03, all 42 conditions in four shards
*! (about 75 minutes; the a1p designs dominate), after smoke runs on mc_d,
*! mc_c_999, mc_a2n_40 and mc_b1_m1 (an empty-category replicate). Summary:
*! tools/near_bound_summary.R. Result and the decision it informed:
*! E-009, 2026-10-03 (the near-bound note).
*!
*! A diagnostic tool, not a pass/fail test (tools/). Evidence for the deferred
*! note on status-1 pairs near the bound (spec Module 2, Stage 1 reporting;
*! E-009, 2026-10-03). Its threshold is open and is the user's; this
*! file records, it does not flag. Designed with the user 2026-10-03: all 42
*! conditions, every pair of every replicate, the score at the nearer bound,
*! distributions only (no thresholds set before the run).
*!
*! It reads exactly the replicate .dta files tools/paired_mc.do reads and
*! generates no data. Every pair is refitted with dwls_polychoric() from
*! Module 1's thresholds, so rho-hat and status are recorded whether or not
*! Stage 1 accepts the replicate. Stage 1 is then called twice, default and
*! keep_bound = 1, and each call's rc is recorded; Gamma's pair entry comes
*! from whichever call succeeded (they agree when both do: keep_bound only
*! changes which status-2 pairs are allowed). Where Stage 1 succeeded, the
*! refitted rho-hat is asserted equal to its s entry.
*!
*! FAILURES ARE RECORDED, NEVER EXCLUDED. A replicate whose thresholds abort
*! (an empty declared category) writes one row with k missing and rc_thr set.
*!
*! B is written into every row (0.999, dwls_polychoric()'s literal; spec
*! Module 2, The bound). The tool repeats the literal; if B changes, change
*! it here too. score_* is dwlsmv_rf_eval() at sign(rho-hat) * B, sign(0) = +1:
*! the evaluation the endpoint check makes. usable 0 means +-B is not an
*! evaluable point for that table (an observed cell's pi is at or below
*! delta_P there, production's own test) and S/N is missing: structurally
*! undefined, not a missing observation, and distinct from a row with k
*! missing (Module 1 aborted). B - epsilon is deliberately not tried; it would
*! be a different diagnostic from the endpoint score. Summaries report the
*! usable share by design before any quantile of score_SN.
*!
*! Columns: condition, replicate, k, i, j, B, rho, status, iterations, N,
*! gamma_rr (Gamma[T+k, T+k], N times the asymptotic variance), gamma_from
*! (default, keep_bound, none), empty_cell (any empty cell in the pair's
*! table), score_x, score_usable, score_SN, score_tauS, rc_thr, rc_default,
*! rc_keep. Reals are written %21.17g, which round-trips a double exactly.
*!
*! Output: <c(tmpdir)>/gammaC/nb/pairs/<condition>.csv, outside the repo.
*!
*! USAGE, from the repository root, one wrapper per shard (Stata batch names
*! its log after the invoked do-file):
*!     wrapper nb_s1.do:  cd <repo>  /  do tools/near_bound_pairs.do 1 4
*! Arguments: shard, nshards [, one condition [, max replicates]].
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
args shard nshards onecond nmax
if "`shard'" == "" | "`nshards'" == "" {
    display as error "usage: do tools/near_bound_pairs.do shard nshards [condition [max]]"
    exit 198
}
if "`nmax'" == "" local nmax = .

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

real colvector nb_kof(string scalar thr)
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

string scalar nb_num(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

void nb_line(string scalar s)
{
    external real scalar NB_FH
    fput(NB_FH, s)
}

void nb_open(string scalar f)
{
    external real scalar NB_FH

    if (fileexists(f)) unlink(f)
    NB_FH = fopen(f, "w")
    nb_line("condition,replicate,k,i,j,B,rho,status,iterations,N,gamma_rr,gamma_from," +
            "empty_cell,score_x,score_usable,score_SN,score_tauS,rc_thr,rc_default,rc_keep")
}

void nb_close()
{
    external real scalar NB_FH
    fclose(NB_FH)
}

// One replicate: the data in memory, K in NB_K. Writes one row per pair, or
// one row with k missing if the thresholds abort.
void nb_one(string scalar cond, real scalar rep)
{
    external real colvector NB_K
    struct dwlsmv_statmap scalar map
    struct dwlsmv_stage1 scalar sd, sk
    struct dwlsmv_pairfit scalar f
    real matrix Y, n, tau
    real rowvector e
    real colvector ti, tj
    real scalar B, j, k, pos, vi, vj, rct, rcd, rck, g, x, ok
    string scalar from, pre

    B = 0.999
    Y = st_data(., .)
    map = dwlsmv_statmap_build(NB_K)
    pre = cond + "," + strofreal(rep) + ","

    // Module 1's thresholds, one column per variable (padded with missing)
    tau = J(max(NB_K) - 1, map.p, .)
    rct = 0
    for (j = 1; j <= map.p; j++) {
        rct = _nb_thr(Y[., j], NB_K[j], ti)
        if (rct) break
        tau[|1, j \ rows(ti), j|] = ti
    }

    rcd = _nb_ss(Y, NB_K, 0, sd)
    rck = _nb_ss(Y, NB_K, 1, sk)

    if (rct) {
        nb_line(pre + ".,.,.," + nb_num(B) + ",.,.,.," + strofreal(rows(Y)) +
                ",.,none,.,.,.,.,.," + strofreal(rct) + "," + strofreal(rcd) + "," +
                strofreal(rck))
        return
    }

    for (k = 1; k <= map.P; k++) {
        pos = map.T + k
        vi  = map.var_i[pos]
        vj  = map.var_j[pos]
        ti  = tau[|1, vi \ NB_K[vi] - 1, vi|]
        tj  = tau[|1, vj \ NB_K[vj] - 1, vj|]
        f   = dwls_polychoric(Y[., vi], Y[., vj], ti, tj)
        n   = dwlsmv_crosstab(Y[., vi], Y[., vj], NB_K[vi], NB_K[vj])

        g = .
        from = "none"
        if (rcd == 0) {
            if (f.rho != sd.s[pos]) _error("near_bound_pairs: refit differs from Stage 1 (default)")
            g = sd.Gamma[pos, pos]
            from = "default"
        }
        else if (rck == 0) {
            if (f.rho != sk.s[pos]) _error("near_bound_pairs: refit differs from Stage 1 (keep_bound)")
            g = sk.Gamma[pos, pos]
            from = "keep_bound"
        }

        x = (f.rho < 0 ? -B : B)
        e = dwlsmv_rf_eval(x, n, ti, tj)
        ok = e[2]

        nb_line(pre + strofreal(k) + "," + strofreal(vi) + "," + strofreal(vj) + "," +
                nb_num(B) + "," + nb_num(f.rho) + "," + strofreal(f.status) + "," +
                strofreal(f.iterations) + "," + strofreal(rows(Y)) + "," + nb_num(g) + "," +
                from + "," + strofreal(any(n :== 0)) + "," + nb_num(x) + "," +
                strofreal(ok) + "," + nb_num(ok ? e[3] : .) + "," + nb_num(ok ? e[5] : .) +
                ",0," + strofreal(rcd) + "," + strofreal(rck))
    }
}

// Module 1 under a captured call: returns the rc, ti set on success.
real scalar _nb_thr(real colvector y, real scalar K, real colvector ti)
{
    external real colvector NB_Y1
    external real scalar NB_K1
    real scalar rc

    NB_Y1 = y
    NB_K1 = K
    rc = _stata("mata: NB_T1 = dwlsmv_thresholds(NB_Y1, NB_K1)", 1)
    if (rc == 0) {
        ti = _nb_get_t1()
    }
    return(rc)
}

real colvector _nb_get_t1()
{
    external real colvector NB_T1
    return(NB_T1)
}

// Stage 1 under a captured call: returns the rc, st set on success.
real scalar _nb_ss(real matrix Y, real colvector K, real scalar kb,
                   struct dwlsmv_stage1 scalar st)
{
    external real matrix NB_YY
    external real colvector NB_KK
    external real scalar NB_KB
    external struct dwlsmv_stage1 scalar NB_ST
    real scalar rc

    NB_YY = Y
    NB_KK = K
    NB_KB = kb
    rc = _stata("mata: NB_ST = dwls_sample_stats(NB_YY, NB_KK, NB_KB)", 1)
    if (rc == 0) {
        st = NB_ST
    }
    return(rc)
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

local outdir "`c(tmpdir)'/gammaC/nb/pairs"
capture mkdir "`c(tmpdir)'/gammaC"
capture mkdir "`c(tmpdir)'/gammaC/nb"
capture mkdir "`outdir'"

foreach i of local todo {
    use `conds', clear
    local c   = condition[`i']
    local thr = thresholds[`i']
    mata: NB_K = nb_kof(st_local("thr"))
    import delimited using "fixtures/replication/`c'/index.csv", clear ///
        varnames(1) delimiters(",") stringcols(_all)
    local nr = min(_N, `nmax')
    tempfile idx
    quietly save `idx'
    mata: nb_open("`outdir'/`c'.csv")
    timer clear 1
    timer on 1
    forvalues r = 1/`nr' {
        use `idx', clear
        local rep  = replicate[`r']
        local file = file[`r']
        quietly use "`file'", clear
        mata: nb_one("`c'", `rep')
    }
    mata: nb_close()
    timer off 1
    quietly timer list 1
    display as result "near_bound_pairs: `c'  `nr' replicates  (" %6.1f r(t1) " s)"
}
display as result "near_bound_pairs: shard `shard' of `nshards' done"
