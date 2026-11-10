*! Module 12, fit indices, option A (family S, the c1 baseline, dCFI/dRMSEA).
*! Executed 2026-10-21 under StataMP 16.1, 4 processors: passes. Mutants
*! (2026-10-21, applied to src/ in place, each checked applied, then
*! restored): M1 CFI/TLI conditions merged (lavaan's form) red at 1; M2 TLI
*! with max(E_L, 0) red at 1, 5; M3 no sqrt(G) red at 1; M4 T_B with N not
*! N - G red at 4, 5, 6; M5 pi_g dropped red at 4, 6; M6 delta against the
*! previous row, not the last fitted one, red at 3 (after check 3's data were
*! changed so rows 2 and 3 differ: with them equal M6 survived); M7 an
*! invalid baseline blanking RMSEA red at 1; M8 n = N red at 5, 6; M9 (check 7,
*! added with the display) the RMSEA line removed red at 7, fit and replay.
*!
*! Decisions under test: option A (user, 2026-10-10); F1-F6 frozen by the
*! user 2026-10-21, F4 as amended (model side and baseline side independent).
*! Derivation @sec-fit-indices (Results FI1-FI3, chunks fi-c1, fi-frozen).
*!
*! Unit checks (Mata):
*!   1  dwlsmv_fitidx, the four sign combinations of (E_L, E_B) at D = 20,
*!      DB = 30, n = 998, G = 2: RMSEA = sqrt(G) sqrt(max(E_L,0)/(n D)) in
*!      every case and identical to the no-baseline call; CFI/TLI defined
*!      only when E_B > 0; E_L <= 0 with E_B > 0 gives CFI = 1 and TLI > 1
*!      (unclipped); E_L <= 0 with E_B <= 0 gives CFI MISSING (the merged
*!      mutant gives 1). F3: D = 0 or T3 missing, all missing. sqrt(G):
*!      RMSEA at G = 2 is sqrt(2) times G = 1's. Values to 1e-14 relative.
*!   2  dwlsmv_fitidx refusals, _rc == 3498: n = 0, n missing, G = 1.5.
*!   3  dwlsmv_fitidx_delta: a four-level chain with level 3 skipped: row 1
*!      and row 3 missing, row 2 = X2 - X1, row 4 = X4 - X2 (not X4 - X3).
*!      Refusals, _rc == 3498: first level skipped, skip = 2, three columns.
*!   4  dwlsmv_fitidx_base on design E's two groups (mg_E.csv, unequal N,
*!      Stage 1 by dwls_sample_stats): T_B against the generic stacked
*!      quadratic form (N - G) r'Wr with W = blockdiag(pi_g / diag Gamma_g)
*!      and Delta_B the threshold selection; Tr(M) and Tr(M^2) against
*!      U Omega built generically; a_B, T3_B from those traces. 1e-10
*!      relative (two orderings of the same sums). Refusal: a non-positive
*!      polychoric variance, _rc == 3498.
*! End to end (dwlsmv, e()):
*!   5  case E, single group, against fixtures/expected/lavaan_fitindices.csv
*!      (99f_, lavaan 0.7.2): two-factor (E_2f, at the clamps) and one-factor
*!      (E_1f, away from them): e(T3_b), e(rmsea_s), e(cfi_s), e(tli_s) at
*!      1e-8 abs + 1e-5 rel (spec 7.2's test-statistic tier, chosen before
*!      the first run); e(df_T3_b) exact.
*!   6  design E, two groups, configural and scalar: e(T3_b) against 99f_'s
*!      T3B_NG (dwlsmv's T) at the same tier, identical at both levels;
*!      e(df_T3_b) = 20; the three indices equal dwlsmv_fitidx(e(T3),
*!      e(df_T3), N - 2, 2, e(T3_b), e(df_T3_b)) exactly.
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

local nfail 0

import delimited using "fixtures/expected/lavaan_fitindices.csv", ///
    clear varnames(1) stringcols(_all) delimiters(",")
mata: FX = st_sdata(., .)
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
mata: YALL = st_data(., .)

mata:

real scalar tfi_fx(string scalar cs, string scalar q)
{
    external string matrix FX
    real colvector k

    k = selectindex((FX[., 1] :== cs) :& (FX[., 3] :== q))
    if (rows(k) != 1) {
        _error("test_fitidx: fixture row not found")
    }
    return(strtoreal(FX[k, 4]))
}

real scalar tfi_rel(real scalar a, real scalar b)
{
    return(abs(a - b) / max((abs(b), 1e-300)))
}

// check 7: exactly one log line holds the label at column 51 and, after
// "= ", the value as %9.4f prints it
real scalar tfi_line(string scalar f, string scalar lab, string scalar val)
{
    string colvector L
    real scalar k, n, h

    L = cat(f)
    n = 0
    for (k = 1; k <= rows(L); k++) {
        if (substr(L[k], 51, strlen(lab)) != lab) continue
        h = strpos(L[k], "= ")
        if (h == 0) continue
        if (strtrim(substr(L[k], h + 2, .)) == val) n++
    }
    return(n == 1)
}

// check 1
real scalar tfi_c1()
{
    real scalar D, DB, n, G, bad, k, EL, EB, r0
    real rowvector f, f0, f1
    real matrix C

    D = 20; DB = 30; n = 998; G = 2; bad = 0
    C = (15, 300 \ -4, 300 \ 15, -2 \ -4, -2)
    for (k = 1; k <= 4; k++) {
        EL = C[k, 1]; EB = C[k, 2]
        f  = dwlsmv_fitidx(D + EL, D, n, G, DB + EB, DB)
        f0 = dwlsmv_fitidx(D + EL, D, n, G, ., .)
        r0 = sqrt(G) * sqrt(max((EL, 0)) / (n * D))
        if (f[1] != f0[1]) bad++
        if (tfi_rel(f[1], r0) > 1e-14 & r0 > 0) bad++
        if (r0 == 0 & f[1] != 0) bad++
        if (EB > 0) {
            if (missing(f[2]) | missing(f[3])) bad++
            if (tfi_rel(f[2], 1 - max((EL, 0)) / max((EL, EB))) > 1e-14) bad++
            if (tfi_rel(f[3], 1 - EL * DB / (EB * D)) > 1e-14) bad++
        }
        else {
            if (!missing(f[2]) | !missing(f[3])) bad++
        }
    }
    f = dwlsmv_fitidx(D - 4, D, n, G, DB + 300, DB)
    if (f[2] != 1 | !(f[3] > 1)) bad++
    f = dwlsmv_fitidx(D - 4, D, n, G, DB - 2, DB)
    if (!missing(f[2])) bad++
    if (missing(dwlsmv_fitidx(5, 0, n, G, 330, 30)) != 3) bad++
    if (missing(dwlsmv_fitidx(., D, n, G, 330, 30)) != 3) bad++
    f1 = dwlsmv_fitidx(D + 15, D, n, 1, DB + 300, DB)
    f  = dwlsmv_fitidx(D + 15, D, n, 2, DB + 300, DB)
    if (tfi_rel(f[1], sqrt(2) * f1[1]) > 1e-14) bad++
    if (f[2] != f1[2] | f[3] != f1[3]) bad++
    // survivor ledger (2026-10-31): DB = 0 and E_B = 0 leave CFI and TLI
    // missing, RMSEA not (FI-F2, FI-F3); E_L > E_B > 0 takes CFI's max arm,
    // CFI 0 (FI-F5)
    f = dwlsmv_fitidx(D + 15, D, n, G, 5, 0)
    if (missing(f[1]) | !missing(f[2]) | !missing(f[3])) bad++
    f = dwlsmv_fitidx(D + 15, D, n, G, DB, DB)
    if (missing(f[1]) | !missing(f[2]) | !missing(f[3])) bad++
    f = dwlsmv_fitidx(D + 50, D, n, G, DB + 10, DB)
    if (f[2] != 0 | tfi_rel(f[3], 1 - 50 * DB / (10 * D)) > 1e-14) bad++
    return(bad)
}

// check 3
real scalar tfi_c3()
{
    real matrix X, d
    real scalar bad

    bad = 0
    X = (.90, .05 \ .92, .04 \ .95, .02 \ .89, .06)
    d = dwlsmv_fitidx_delta(X, (0 \ 0 \ 1 \ 0))
    if (missing(d[1, .]) != 2 | missing(d[3, .]) != 2) bad++
    if (d[2, .] != X[2, .] - X[1, .]) bad++
    if (d[4, .] != X[4, .] - X[2, .]) bad++
    return(bad)
}

// check 4: closed form against the generic stacked construction
real scalar tfi_c4()
{
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    real colvector KE, sel, s, wd, isr
    real matrix Om, W, Dl, U, M
    real rowvector b
    real scalar g, G, N, o, k, bad, TB, tr1, tr2, a, ng

    KE  = (4 \ 4 \ 4 \ 2 \ 2)
    G   = 2
    sts = dwlsmv_stage1(G, 1)
    for (g = 1; g <= G; g++) {
        sel    = selectindex(YALL[., 1] :== g)
        sts[g] = dwls_sample_stats(YALL[sel, 2..6], KE)
    }
    N = sts[1].N + sts[2].N
    if (sts[1].N == sts[2].N) return(99)
    b = dwlsmv_fitidx_base(sts)

    // stacked s (thresholds kept, residual 0 at the optimum), W, Omega, Delta_B
    k   = rows(sts[1].s) + rows(sts[2].s)
    s   = J(0, 1, .); wd = J(0, 1, .); isr = J(0, 1, .)
    Om  = J(k, k, 0)
    ng  = sum(sts[1].map.stat_type :== 1) + sum(sts[2].map.stat_type :== 1)
    Dl  = J(k, ng, 0)
    o   = 0
    for (g = 1; g <= G; g++) {
        s   = s \ (sts[g].s :* (sts[g].map.stat_type :== 2))
        wd  = wd \ ((sts[g].N / N) :/ diagonal(sts[g].Gamma))
        isr = isr \ (sts[g].map.stat_type :== 1)
        Om[|o + 1, o + 1 \ o + rows(sts[g].s), o + rows(sts[g].s)|] = sts[g].Gamma :/ (sts[g].N / N)
        o = o + rows(sts[g].s)
    }
    o = 0
    for (k = 1; k <= rows(isr); k++) {
        if (isr[k]) {
            o++
            Dl[k, o] = 1
        }
    }
    W   = diag(wd)
    U   = W - W * Dl * invsym(Dl' * W * Dl) * Dl' * W
    M   = U * Om
    TB  = (N - G) * (s' * W * s)
    tr1 = trace(M)
    tr2 = trace(M * M)
    a   = sqrt(tr1 / tr2)
    bad = 0
    if (b[2] != 20 | abs(tr1 - 20) > 1e-10) bad++
    if (tfi_rel(b[3], TB) > 1e-10) bad++
    if (tfi_rel(b[4], a) > 1e-10) bad++
    if (tfi_rel(b[1], a * TB + 20 - sqrt(20 * tr1^2 / tr2)) > 1e-10) bad++
    printf("check 4: T_B %21.15g generic %21.15g; a_B %g / %g; T3_B %g\n", b[3], TB, b[4], a, b[1])
    return(bad)
}

void tfi_c4_bad()
{
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    real colvector KE, j

    KE  = (4 \ 4 \ 4 \ 2 \ 2)
    sts = dwlsmv_stage1(1, 1)
    sts[1] = dwls_sample_stats(YALL[selectindex(YALL[., 1] :== 1), 2..6], KE)
    j = selectindex(sts[1].map.stat_type :== 2)
    sts[1].Gamma[j[1], j[1]] = 0
    (void) dwlsmv_fitidx_base(sts)
}

// survivor ledger (2026-10-31), FI-F13: a missing polychoric aborts
void tfi_c4_miss()
{
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    real colvector KE, j

    KE  = (4 \ 4 \ 4 \ 2 \ 2)
    sts = dwlsmv_stage1(1, 1)
    sts[1] = dwls_sample_stats(YALL[selectindex(YALL[., 1] :== 1), 2..6], KE)
    j = selectindex(sts[1].map.stat_type :== 2)
    sts[1].s[j[1]] = .
    (void) dwlsmv_fitidx_base(sts)
}

// survivor ledger (2026-10-31), FI-F14: no polychoric row gives (., 0, ., .);
// the map is relabelled so that no row is a correlation
real scalar tfi_c4_nopair()
{
    external real matrix YALL
    struct dwlsmv_stage1 vector sts
    real colvector KE
    real rowvector b

    KE  = (4 \ 4 \ 4 \ 2 \ 2)
    sts = dwlsmv_stage1(1, 1)
    sts[1] = dwls_sample_stats(YALL[selectindex(YALL[., 1] :== 1), 2..6], KE)
    sts[1].map.stat_type = J(rows(sts[1].map.stat_type), 1, 1)
    b = dwlsmv_fitidx_base(sts)
    return(!(missing(b[1]) & b[2] == 0 & missing(b[3]) & missing(b[4])))
}

end

* ---- check 1 ---------------------------------------------------------------
mata: st_local("b1", strofreal(tfi_c1()))
if `b1' {
    display as error "check 1 FAILED: `b1' violation(s)"
    local ++nfail
}
else display "check 1 passed"

* ---- check 2 ---------------------------------------------------------------
foreach call in "25, 20, 0, 2, 330, 30" "25, 20, ., 2, 330, 30" "25, 20, 998, 1.5, 330, 30" {
    capture mata: (void) dwlsmv_fitidx(`call')
    if _rc != 3498 {
        display as error "check 2 FAILED: dwlsmv_fitidx(`call') rc " _rc
        local ++nfail
    }
}
display "check 2 done"

* ---- check 3 ---------------------------------------------------------------
mata: st_local("b3", strofreal(tfi_c3()))
if `b3' {
    display as error "check 3 FAILED: `b3' violation(s)"
    local ++nfail
}
capture mata: (void) dwlsmv_fitidx_delta((.9, .05 \ .91, .04), (1 \ 0))
if _rc != 3498 {
    display as error "check 3 FAILED: first level skipped, rc " _rc
    local ++nfail
}
capture mata: (void) dwlsmv_fitidx_delta((.9, .05 \ .91, .04), (0 \ 2))
if _rc != 3498 {
    display as error "check 3 FAILED: skip = 2, rc " _rc
    local ++nfail
}
capture mata: (void) dwlsmv_fitidx_delta((.9, .05, 1 \ .91, .04, 1), (0 \ 0))
if _rc != 3498 {
    display as error "check 3 FAILED: three columns, rc " _rc
    local ++nfail
}
* survivor ledger (2026-10-31), FI-F16: skip of another length
capture mata: (void) dwlsmv_fitidx_delta((.9, .05 \ .91, .04), (0))
if _rc != 3498 {
    display as error "check 3 FAILED: short skip, rc " _rc
    local ++nfail
}
display "check 3 done"

* ---- check 4 ---------------------------------------------------------------
mata: st_local("b4", strofreal(tfi_c4()))
if `b4' {
    display as error "check 4 FAILED: `b4' violation(s)"
    local ++nfail
}
capture mata: tfi_c4_bad()
if _rc != 3498 {
    display as error "check 4 FAILED: non-positive variance, rc " _rc
    local ++nfail
}
capture mata: tfi_c4_miss()
if _rc != 3498 {
    display as error "check 4 FAILED: missing polychoric, rc " _rc
    local ++nfail
}
mata: st_local("b4n", strofreal(tfi_c4_nopair()))
if `b4n' {
    display as error "check 4 FAILED: no polychoric row does not give (., 0, ., .)"
    local ++nfail
}
display "check 4 done"

* ---- check 5 ---------------------------------------------------------------
foreach cs in E_2f E_1f {
    import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
    if "`cs'" == "E_2f" quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6)
    else                quietly dwlsmv (F -> y1 y2 y3 y4 y5 y6)
    foreach pr in "T3_b baseline.chisq.scaled" "rmsea_s rmsea.scaled" "cfi_s cfi.scaled" "tli_s tli.scaled" "df_T3 df" {
        gettoken e q : pr
        local q = strtrim("`q'")
        mata: st_local("r", strofreal(dwlsmv_gapratio(st_numscalar("e(`e')"), tfi_fx("`cs'", "`q'"), 1e-8, 1e-5)))
        if !(`r' <= 1) {
            display as error "check 5 FAILED: `cs' e(`e') " e(`e') " vs lavaan `q' (ratio `r')"
            local ++nfail
        }
    }
    mata: st_local("r", strofreal(st_numscalar("e(df_T3_b)") == tfi_fx("`cs'", "baseline.df")))
    if !`r' {
        display as error "check 5 FAILED: `cs' e(df_T3_b) " e(df_T3_b)
        local ++nfail
    }
    display "check 5 `cs': rmsea " e(rmsea_s) " cfi " e(cfi_s) " tli " e(tli_s) " T3_b " e(T3_b)
}

* ---- check 6 ---------------------------------------------------------------
scalar tfi_T3b0 = .
foreach lev in configural scalar {
    import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
    quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(`lev') metric(unitvariance)
    mata: st_local("r", strofreal(dwlsmv_gapratio(st_numscalar("e(T3_b)"), tfi_fx("mgE", "T3B_NG"), 1e-8, 1e-5)))
    if !(`r' <= 1) | e(df_T3_b) != 20 {
        display as error "check 6 FAILED: `lev' e(T3_b) " e(T3_b) " (ratio `r'), df " e(df_T3_b)
        local ++nfail
    }
    if "`lev'" == "scalar" & e(T3_b) != tfi_T3b0 {
        display as error "check 6 FAILED: e(T3_b) differs between levels"
        local ++nfail
    }
    scalar tfi_T3b0 = e(T3_b)
    mata: f = dwlsmv_fitidx(st_numscalar("e(T3)"), st_numscalar("e(df_T3)"), st_numscalar("e(N)") - 2, 2, ///
        st_numscalar("e(T3_b)"), st_numscalar("e(df_T3_b)"))
    mata: st_local("r", strofreal(f == (st_numscalar("e(rmsea_s)"), st_numscalar("e(cfi_s)"), st_numscalar("e(tli_s)"))))
    if !`r' | ("`lev'" == "configural" & !(e(rmsea_s) > 0)) {
        display as error "check 6 FAILED: `lev' indices are not dwlsmv_fitidx of e()"
        local ++nfail
    }
    display "check 6 `lev': T3 " e(T3) " df " e(df_T3) " rmsea " e(rmsea_s) " cfi " e(cfi_s) " tli " e(tli_s)
}

* ---- check 7: the printed lines, at estimation and on replay -----------------
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
foreach mode in fit replay {
    tempfile lg
    capture log close tfi
    quietly log using "`lg'", text replace name(tfi)
    if "`mode'" == "fit" dwlsmv (F -> y1 y2 y3 y4 y5 y6)
    else                 dwlsmv
    quietly log close tfi
    foreach pr in "RMSEA (scaled) rmsea_s" "CFI (scaled) cfi_s" "TLI (scaled) tli_s" {
        local e : word 3 of `pr'
        local lab = substr("`pr'", 1, strpos("`pr'", ")"))
        mata: st_local("ok", strofreal(tfi_line(st_local("lg"), st_local("lab"), ///
            strtrim(sprintf("%9.4f", st_numscalar("e(`e')"))))))
        if !`ok' {
            display as error "check 7 FAILED (`mode'): no line '`lab' ... = " %9.4f e(`e') "'"
            local ++nfail
        }
    }
}
display "check 7 done"

if `nfail' {
    display as error _n "test_fitidx: `nfail' failure(s)"
    exit 9
}
display _n "test_fitidx: all checks passed"
