*! Executed by the main session 2026-11-15 (StataMP 18): passes on its first run; mutants in docs/notes/condx-stage2/m10-build/mutants/.
*! Module 10 diagnostics for a conditional fit (one binary x, one group), decisions of 2026-11-15
*! (docs/decisions.md, "Module 10 diagnostics in conditional mode: five points accepted"; "Module 10 aliasing row in
*! conditional mode: (a) now, (b) later"). Needs the drafted src/mata/dwlsmv_countcells.mata, the edits of
*! dwlsmv_npatterns, dwlsmv_gammadiag, dwlsmv_diagadvice, dwlsmv_m11_{fit,slab,stage1,clear}.mata and dwlsmv.ado.
*! Data: fixtures/data/condx_A.csv and condx_D.csv (y1 y2 y3 x; K = (3, 3, 2); N = 1500; D = 5 + 3 + 3 = 11), read with
*! explicit delimiters. Every expected number comes from Stata (egen tag, count, -duplicates-style tags), never from the
*! Mata under test, except where a check says it is a Mata-to-Mata consistency check.
*!
*! Checks, per dataset
*!   1  dwlsmv_npatterns(Y, K, x) equals the number of distinct (y1 y2 y3 x) rows counted by egen tag(); the 2-argument
*!      count equals the y-only tag count, and the two differ (so a count that ignores x is caught); a hand-built 8-row design
*!      (4 y-patterns, 8 (y,x)-patterns) fixes both numbers by inspection.
*!   2  dwlsmv_countcells(): the whole table equals Stata -count- cell by cell; per (item, category) the two group counts
*!      sum to the pooled catcounts; per item each group's counts sum to its group size; rows 1-2 are the group sizes.
*!   3  dwlsmv_gammadiag(st, Y, x): 7 x 4; row 3 = (n_pat - 1 < D, n_pat, D, 0) with n_pat from check 1's Stata count; row 4
*!      is (., ., 0, 0) EXACTLY (not 0: missing, and not .n); rows 6-7 have r* = min(D, n_pat - 1) and the scale-free ratio
*!      equals sigma_r*/sigma_1 of svd() of G :/ (s s') to 1e-10 relative; rows 1-2 unflagged; row 5 unflagged (healthy fixtures:
*!      z is expected far below 1e-6; a failure here is a FINDING about the pooled-score pair bound, not a typo).
*!   4  rc 3498 (asserted == 3498, never != 0): 2-argument call on a conditional result; 3-argument call on an unconditional
*!      result; x of the wrong length; npatterns with x coded 2, with a missing x, with the wrong number of rows;
*!      countcells with a constant x and with a missing x.
*!   5  Regression, ordinary fit: on fixtures/data/case_D.csv (K = 4,4,3,3) dwlsmv_npatterns(Y, K) = 139 (the recorded count),
*!      dwlsmv_gammadiag(st, Y) row 4 is evaluated (flag 0, not missing), rows 6-7 .n; dwlsmv_countdiag unchanged shape 2 x 4; and
*!      dwlsmv_m11_slab() on an ordinary map (no third argument) returns the labels it returned (thresholds and pairs only).
*!   6  Through the command (model (F -> y1 y2 y3) (F <- x)): e(gamma_ratio), e(gamma_ratio_raw), e(gamma_rank) are posted and
*!      equal the direct Mata values (1e-15 relative: the %21.0g local); e(bread_ratio) is posted and equals dwlsmv_breaddiag()'s row 4
*!      of the same fit recomputed here (1e-12 relative; the fit is rerun with dwlsmv_fit_noretry, so the check applies only when
*!      e(recovery) is "none"), with row 3 unflagged; the printed output has
*!      the heading "Category counts within each value of x", the group sizes line with Stata's counts, the lowest-or-highest
*!      group-cell line the Stata counts predict, the slope sentence, the "not carried out" sentence with its limitation
*!      wording; the pooled heading names "both groups together"; and it has no "Too few distinct" note. A negative control: an
*!      ordinary fit (no covariate) posts e(bread_ratio), prints none of the new lines, and its pooled heading is the old one.
*!   7  The slope labels: dwlsmv_m11_slab(map, ind, "x") names slope rows "y1~x", "y2~x", "y3~x" and no row reads "~~" with a
*!      repeated indicator.
*!   8  (Mata level, no fit) dwlsmv_m11_rectrigger() with a hand-built converged fit, admissible A and Bd: an unflagged Bd gives "";
*!      Bd row 3 = 1 gives the bread cause ("could not determine all parameters reliably") and no other; Bd row 2 = 1 alone gives "";
*!      an empty Bd (0 x 0) gives "" and does not abort. This is the conditional trigger's bread term; the same function serves
*!      ordinary fits, so the check also guards them.
*!   9  (executed 2026-11-15) A SPARSE design built in this file from condx_A (nothing committed): the (y3 = 1, x = 1) cell cut to its first
*!      two rows in file order (2 observations; pooled y3 = 1 stays 343). Checks 1-3 (td_check, so countcells against -count- cell by
*!      cell, the sums, gammadiag rows) run on it; countcells reads 2 in that cell; the log's within-group "lowest or highest category"
*!      line is "2 (y3, highest, x = 1)" and the pooled line "343 (y3, highest)", both from Stata's counts, and the two differ.
*!  10  (executed 2026-11-15) dwlsmv_gammadiag() row 5 in conditional mode excludes nothing. Mata level, on a hand-built Stage 1 result
*!      (K = 2, 2, 3; Gamma = I or diagonal; planted psibar) over a Y whose pooled (y1,y2) table has exactly one empty cell: a value
*!      planted at each of the D = 10 entries is row 5's (1, value, entry, 0); a mixed psibar and diagonal Gamma give the direct
*!      maximum over all D entries with its first location; an ORDINARY map on the same Y does exclude the pair's rho entry (the
*!      control that the Y holds an excludable pair). SCOPE: the excludable pair here is two BINARY items, and for such a pair
*!      on the tested conditional path (one binary x, one group) the Y cannot come out of dwls_sample_stats_cx(): an empty
*!      pooled cell is empty in both groups, the pair has no interior root, status 2 (agent's argument; one design in R). Not
*!      established for items with three or more categories, where Stage 1 may admit a pair with one empty pooled cell. Check
*!      10 verifies the downstream handling (row 5 excludes nothing) whatever the source; it does not show the state unreachable.
*!   9R  (added 2026-11-15) the sparse fit of check 9 is a recorded real fit that reaches the bread cause: asserted
*!      e(converged) 0, e(recovery) "failed", e(bread_ratio) < 1e-8, the "rank 8 of 9" bread note, the trigger line naming all
*!      three causes, and no exposed inference (e(V) not posted, e(T3)/e(p_T3) not posted, no "Std. err." or "Prob > chi2"
*!      printed, the header's "No standard errors, T3 or p-value" line). Recorded on StataMP 18 at default processors only.
*!   The printed bread notes for a conditional fit are covered by 9R (the rank-8-of-9 note); other bread notes are not.
*!
*! WHAT THIS CERTIFIES AND DOES NOT. The output and posting of Module 10 for a conditional fit, and the exactness of the
*! counts. It certifies no accuracy of Gamma-hat, says nothing about the mixed-orientation aliasing limitation beyond the
*! printed wording (that is probe_cx_delta.do part B), and nothing about the bread ratio beyond its posting and the trigger term (check 8);
*! the 1e-8 criterion itself is measured by the main session's probe (docs/notes/condx-stage2/m10-draft/probe/), not here.
*!
*! MUTANTS the file must kill (edit, run, revert; each must turn this file red). Files under src/.
*!  N1  dwlsmv_npatterns.mata: `return(rows(uniqrows((Y, x))))` -> `return(rows(uniqrows(Y)))`.          Killed by 1 (and 3, 6).
*!  N2  dwlsmv_npatterns.mata: delete the `x :!= 0 & x :!= 1` block.                                       Killed by 4 (x coded 2 gives rc 0).
*!  N3  dwlsmv_npatterns.mata: delete the `rows(x) != rows(Y)` block.                                      Killed by 4.
*!  G1  dwlsmv_gammadiag.mata: cx branch `npat = dwlsmv_npatterns(Y, st.map.K, x)` -> `(Y, st.map.K)`.     Killed by 3 (row 3 value).
*!  G2  dwlsmv_gammadiag.mata: cx branch `da = 0` -> `da = 1`.                                             Killed by 3 (r* = D - 1 instead of D).
*!  G3  dwlsmv_gammadiag.mata: cx branch set A[4, .] = (0, 0, 0, 0)  (an evaluated "no aliasing").         Killed by 3 (row 4 must be `.`).
*!  G4  dwlsmv_gammadiag.mata: delete the `cx & args() < 3` refusal.  SURVIVED (run 2026-11-15): the `rows(x) != st.N` guard
*!      below still refuses the 2-argument call at rc 3498; only the message differs. No contract violation, no code change.
*!  G5  dwlsmv_gammadiag.mata: delete the `!cx & args() >= 3` refusal.                                     Killed by 4.
*!  G6  dwlsmv_gammadiag.mata: cx branch calls dwlsmv_emptycell_pairs / the old thralias on conditional Y. SURVIVED the first run (2026-11-15, as
*!      `if (!cx)` -> `if (1)`): condx_A and condx_D have no one-empty-cell pair, so the exclusion zeroes nothing. A coverage
*!      gap, not a demonstrated contract violation. KILLED by 10 on the rerun (2026-11-15): on the hand-built Y the mutant's
*!      exclusion zeroes entry T + 2 = 6 (the slope of y2; emptycell_pairs reads (var_i, var_j) = (2, 1) there as the pair y2~y1),
*!      so a value planted at 6 is lost.
*!  C1  dwlsmv_countcells.mata: `(x :== g)` -> `(x :== 0)` in the cell count.                              Killed by 2.
*!  C2  dwlsmv_countcells.mata: group sizes swapped in rows 1-2.                                           Killed by 2 and 6.
*!  C3  dwlsmv_countcells.mata: loop `c <= K[j]` -> `c < K[j]` (top category dropped).                     Killed by 2 (row count / sums).
*!  A1  dwlsmv_diagadvice.mata: drop the within-group block.                                               Killed by 6.
*!  A2  dwlsmv_diagadvice.mata: pooled heading identical for cx and ordinary.                              Killed by 6 (both groups together; negative control).
*!  A3  dwlsmv_diagadvice.mata: the smallest-cell search uses `>` for `<`.                                 Killed by 6 (lowest-or-highest line).
*!  S1  dwlsmv_m11_slab.mata: delete the stat_type == 3 branch.                                            Killed by 7.
*!  F1  dwlsmv_m11_fit.mata: cx branch posts gratio from Gd[6, 2] (raw) instead of Gd[7, 2].              Killed by 6 (e(gamma_ratio) vs direct).
*!  F2  dwlsmv_m11_fit.mata: cx branch does not set `grank`.                                               Killed by 6 (e(gamma_rank) missing/empty).
*!  F3  dwlsmv_m11_fit.mata: pass Y instead of the stored x to countcells/gammadiag (xv = Y[., 1]).        Killed by 2 (via cell counts in the printed lines) and 6.
*!  E1  (retired: the guard around e(bread_ratio) no longer exists; see B4.)
*!  E2  dwlsmv.ado: gamma_* posting left inside an `if "`cov'" == ""` block (the pre-change state).        Killed by 6.
*!  X1  dwlsmv_m11_stage1.mata: the stored x not written (px block removed).                              Killed by 6 (findexternal on a missing external: rc 3499).
*!  X2  dwlsmv_m11_clear.mata: rmexternal("dwlsmv_m11_x") removed.                                         Killed by 6 (the external survives the fit).
*!  B1  dwlsmv_m11_estimate.mata: `bread_cx = 1` -> `bread_cx = 0`.        Killed by 6 (Bd is 0 x 0: the fit aborts reading Bd[4, 1]) and by 8.
*!  B2  dwlsmv_m11_rectrigger.mata: (a) disable the Bd[3, 1] term; (b) drop its `rows(Bd) >= 3` guard. Killed by 8 ((b) by the abort
*!      on the empty Bd, rc 3301 from the file).
*!  B3  dwlsmv_m11_fit.mata: `bratio` from Bd[4, 3], or not set.            Killed by 6 (e(bread_ratio) against the direct row 4).
*!  B4  dwlsmv.ado: e(bread_ratio) left inside an `if "`cov'" == ""` block (the earlier state).     Killed by 6.
*!  B5  dwlsmv_m11_fit.mata: restore the old "Not computed ... the diagnostics of the weight and derivative matrices" text.   Killed by 6.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

global TD_NFAIL 0

capture program drop td_cond
program td_cond
    args ok tag
    if `ok' != 1 {
        display as error "test_condx_diag: FAIL `tag'"
        global TD_NFAIL = ${TD_NFAIL} + 1
    }
end

* run a command with its output written to log file fn; r(rc) is its return code
capture program drop td_run
program td_run, rclass
    gettoken fn 0 : 0
    quietly log using "`fn'", text replace name(tdlog)
    capture noisily `0'
    local rc = _rc
    quietly log close tdlog
    return scalar rc = `rc'
end

* does the log have phrase?  r(hit)
capture program drop td_has
program td_has, rclass
    args fn phrase
    mata: st_local("hit", strofreal(td_has("`fn'", st_local("phrase"))))
    return scalar hit = `hit'
end

mata:

real scalar td_has(string scalar fn, string scalar ph)
{
    string colvector L
    real scalar i

    L = cat(fn)
    for (i = 1; i <= rows(L); i++) {
        if (strpos(L[i], ph) > 0) {
            return(1)
        }
    }
    return(0)
}

// a Y, K, x from the data in memory (values are already 0 .. K-1)
void td_data(real matrix Y, real colvector K, real colvector x)
{
    Y = st_data(., ("y1", "y2", "y3"))
    K = (3 \ 3 \ 2)
    x = st_data(., "x")
}

// expected line for the smallest boundary group cell, from the Stata-built table HC (same row layout as countcells)
string scalar td_expect_bline(real matrix HC, real colvector K)
{
    real scalar r, j, c, g, n, bn, bj, bc, bg

    bn = .
    bj = 0
    bc = 0
    bg = 0
    for (r = 3; r <= rows(HC); r++) {
        j = HC[r, 1]
        c = HC[r, 2]
        g = HC[r, 3]
        n = HC[r, 4]
        if (c == 1 | c == K[j]) {
            if (bn == . | n < bn) {
                bn = n
                bj = j
                bc = c
                bg = g
            }
        }
    }
    return("lowest or highest category: " + strofreal(bn) + " (y" + strofreal(bj) + ", " +
           (bc == 1 ? "lowest" : "highest") + ", x = " + strofreal(bg) + ")")
}

// checks 1-3 on the data in memory; HC is the Stata-built count table, npx and npy the Stata pattern counts
real scalar td_check(real scalar npx, real scalar npy, real matrix HC)
{
    struct dwlsmv_stage1 scalar st
    real matrix Y, A, G, Gs, C
    real colvector K, x, s, sv, pooled
    real scalar nf, D, rs, ratio, j, c

    nf = 0
    td_data(Y, K, x)
    D = 11

    // ---- 1 patterns
    if (dwlsmv_npatterns(Y, K, x) != npx) {
        nf++
        printf("{err}  1: npatterns with x is %f, Stata counts %f\n", dwlsmv_npatterns(Y, K, x), npx)
    }
    if (dwlsmv_npatterns(Y, K) != npy) {
        nf++
        printf("{err}  1: npatterns without x is %f, Stata counts %f\n", dwlsmv_npatterns(Y, K), npy)
    }
    if (npx <= npy) {
        nf++
        printf("{err}  1: premise failed: (y,x) patterns %f do not exceed y patterns %f on this fixture\n", npx, npy)
    }

    // ---- 2 countcells
    C = dwlsmv_countcells(Y, K, x)
    if (rows(C) != rows(HC) | cols(C) != 4) {
        nf++
        printf("{err}  2: countcells is %f x %f, wanted %f x 4\n", rows(C), cols(C), rows(HC))
    }
    else if (any(C :!= HC)) {
        nf++
        printf("{err}  2: countcells differs from Stata -count- in %f entries\n", sum(C :!= HC))
    }
    for (j = 1; j <= cols(Y); j++) {
        pooled = dwlsmv_catcounts(Y[., j], K[j])
        for (c = 1; c <= K[j]; c++) {
            if (sum(select(C[., 4], (C[., 1] :== j) :& (C[., 2] :== c))) != pooled[c]) {
                nf++
                printf("{err}  2: item %f category %f: group counts do not sum to the pooled count\n", j, c)
            }
        }
        if (sum(select(C[., 4], (C[., 1] :== j) :& (C[., 3] :== 0))) != C[1, 4] |
            sum(select(C[., 4], (C[., 1] :== j) :& (C[., 3] :== 1))) != C[2, 4]) {
            nf++
            printf("{err}  2: item %f: a group's counts do not sum to its size\n", j)
        }
    }

    // ---- 3 gammadiag with x
    st = dwls_sample_stats_cx(Y, K, x)
    A  = dwlsmv_gammadiag(st, Y, x)
    if (rows(A) != 7 | cols(A) != 4) {
        return(nf + 1)
    }
    if (A[1, 1] != 0 | A[2, 1] != 0) {
        nf++
        printf("{err}  3: asymmetry or PSD flagged on a healthy fit: %f %f\n", A[1, 1], A[2, 1])
    }
    if (A[3, 1] != (npx - 1 < D) | A[3, 2] != npx | A[3, 3] != D | A[3, 4] != 0) {
        nf++
        printf("{err}  3: row 3 is (%f, %f, %f, %f), wanted (%f, %f, %f, 0)\n", A[3, 1], A[3, 2], A[3, 3], A[3, 4],
               (npx - 1 < D), npx, D)
    }
    // row 4: not evaluated; `.` exactly, neither 0 nor .n
    if (!(A[4, 1] == . & !(A[4, 1] > .) ) | !(A[4, 2] == . & !(A[4, 2] > .)) | A[4, 3] != 0 | A[4, 4] != 0) {
        nf++
        printf("{err}  3: row 4 must be (., ., 0, 0); it is (%f, %f, %f, %f)\n", A[4, 1], A[4, 2], A[4, 3], A[4, 4])
    }
    if (A[5, 1] != 0 | A[5, 2] >= 1e-6) {
        nf++
        printf("{err}  3: centring row (%f, z = %g, k = %f): FINDING, not a typo\n", A[5, 1], A[5, 2], A[5, 3])
    }
    rs = min((D, npx - 1))
    if (A[7, 3] != rs | A[6, 3] != rs | A[6, 1] != .n | A[7, 1] != .n) {
        nf++
        printf("{err}  3: rows 6-7 r* is (%f, %f), wanted %f; flags must be .n\n", A[6, 3], A[7, 3], rs)
    }
    G  = st.Gamma
    s  = sqrt(diagonal(G))
    Gs = G :/ (s * s')
    sv = svdsv(Gs)
    sv = sort(sv, -1)
    ratio = sv[rs] / sv[1]
    if (abs(A[7, 2] - ratio) > 1e-10 * ratio) {
        nf++
        printf("{err}  3: scale-free ratio %g against svd() %g\n", A[7, 2], ratio)
    }
    return(nf)
}

// the 8-row hand design of check 1: y-patterns {00, 01, 10, 11}, each in both groups
real scalar td_hand()
{
    real matrix Y
    real colvector x
    real scalar nf

    nf = 0
    Y  = (0, 0 \ 0, 0 \ 0, 1 \ 0, 1 \ 1, 0 \ 1, 0 \ 1, 1 \ 1, 1)
    x  = (0 \ 1 \ 0 \ 1 \ 0 \ 1 \ 0 \ 1)
    if (dwlsmv_npatterns(Y, (2 \ 2)) != 4) {
        nf++
        printf("{err}  1: hand design, y-only count is %f, by inspection 4\n", dwlsmv_npatterns(Y, (2 \ 2)))
    }
    if (dwlsmv_npatterns(Y, (2 \ 2), x) != 8) {
        nf++
        printf("{err}  1: hand design, (y,x) count is %f, by inspection 8\n", dwlsmv_npatterns(Y, (2 \ 2), x))
    }
    // x constant: the count is the y-only count
    if (dwlsmv_npatterns(Y, (2 \ 2), J(8, 1, 1)) != 4) {
        nf++
        printf("{err}  1: hand design with constant x, count is %f, wanted 4\n", dwlsmv_npatterns(Y, (2 \ 2), J(8, 1, 1)))
    }
    return(nf)
}

// ---- 4: refusals, each run under capture noisily by the caller, so the rc is read there
void td_r1()
{
    struct dwlsmv_stage1 scalar st
    real matrix Y, A
    real colvector K, x

    td_data(Y, K, x)
    st = dwls_sample_stats_cx(Y, K, x)
    A  = dwlsmv_gammadiag(st, Y)
}

void td_r2()
{
    struct dwlsmv_stage1 scalar st
    real matrix Y, A
    real colvector K, x

    // an unconditional Stage 1 result: runs on case_D (y1-y4, in memory), where Stage 1 is known to succeed,
    // so the rc 3498 read by the caller can only come from the x-given refusal
    Y = st_data(., ("y1", "y2", "y3", "y4"))
    K = (4 \ 4 \ 3 \ 3)
    x = J(rows(Y), 1, 0)
    st = dwls_sample_stats(Y, K)
    A  = dwlsmv_gammadiag(st, Y, x)
}

void td_r3()
{
    struct dwlsmv_stage1 scalar st
    real matrix Y, A
    real colvector K, x

    td_data(Y, K, x)
    st = dwls_sample_stats_cx(Y, K, x)
    A  = dwlsmv_gammadiag(st, Y, x[1::(rows(x) - 1)])
}

void td_r4()
{
    real matrix Y
    real colvector K, x

    td_data(Y, K, x)
    x[1] = 2
    (void) dwlsmv_npatterns(Y, K, x)
}

void td_r5()
{
    real matrix Y
    real colvector K, x

    td_data(Y, K, x)
    x[1] = .
    (void) dwlsmv_npatterns(Y, K, x)
}

void td_r6()
{
    real matrix Y
    real colvector K, x

    td_data(Y, K, x)
    (void) dwlsmv_npatterns(Y, K, x[1::10])
}

void td_r7()
{
    real matrix Y
    real colvector K, x

    td_data(Y, K, x)
    (void) dwlsmv_countcells(Y, K, J(rows(x), 1, 0))
}

void td_r8()
{
    real matrix Y
    real colvector K, x

    td_data(Y, K, x)
    x[1] = .
    (void) dwlsmv_countcells(Y, K, x)
}

// ---- 5: ordinary regression on case_D (data in memory: y1 y2 y3 y4)
real scalar td_ordinary()
{
    struct dwlsmv_stage1 scalar st
    real matrix Y, A, C
    real colvector K
    real scalar nf

    nf = 0
    Y  = st_data(., ("y1", "y2", "y3", "y4"))
    K  = (4 \ 4 \ 3 \ 3)
    if (dwlsmv_npatterns(Y, K) != 139) {
        nf++
        printf("{err}  5: case_D npatterns %f, recorded 139\n", dwlsmv_npatterns(Y, K))
    }
    st = dwls_sample_stats(Y, K)
    A  = dwlsmv_gammadiag(st, Y)
    if (rows(A) != 7 | cols(A) != 4) {
        return(nf + 1)
    }
    if (A[4, 1] != 0 | A[4, 2] != 0) {
        nf++
        printf("{err}  5: ordinary row 4 must stay evaluated (0, 0), is (%f, %f)\n", A[4, 1], A[4, 2])
    }
    if (A[6, 1] != .n | A[7, 1] != .n) {
        nf++
        printf("{err}  5: ordinary rows 6-7 flags must be .n\n")
    }
    C = dwlsmv_countdiag(Y, K)
    if (rows(C) != 2 | cols(C) != 4) {
        nf++
        printf("{err}  5: countdiag shape %f x %f\n", rows(C), cols(C))
    }
    return(nf)
}

// ---- 7: slab labels, conditional map and ordinary map
real scalar td_slab()
{
    struct dwlsmv_statmap scalar m, m0
    string colvector L, L0
    real scalar nf, k

    nf = 0
    m  = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
    L  = dwlsmv_m11_slab(m, ("y1", "y2", "y3"), "x")
    if (L[6] != "y1~x" | L[7] != "y2~x" | L[8] != "y3~x") {
        nf++
        printf("{err}  7: slope labels are %s %s %s\n", L[6], L[7], L[8])
    }
    for (k = 1; k <= rows(L); k++) {
        if (strpos(L[k], "y1~~y1") | strpos(L[k], "y2~~y1")) {
            nf++
            printf("{err}  7: a slope row reads as a pair: %s\n", L[k])
        }
    }
    m0 = dwlsmv_statmap_build((3 \ 3 \ 2), 0)
    L0 = dwlsmv_m11_slab(m0, ("y1", "y2", "y3"))
    if (rows(L0) != 8 | L0[1] != "y1|t1" | L0[6] != "y1~~y2") {
        nf++
        printf("{err}  7: ordinary labels changed: %s ... %s\n", L0[1], L0[6])
    }
    return(nf)
}

// ---- 6: e() against the direct values, on the data in memory
// canon and metric are the fit's e(canon) and e(metric): the same first fit as dwlsmv_m11_estimate's, rerun here
// (marker: no retry; else dwls_fit). Valid for comparison only when no recovery replaced the fit (e(recovery) "none").
real rowvector td_direct(string scalar canon, string scalar metric)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    real matrix Y, A, Bd
    real colvector K, x

    td_data(Y, K, x)
    st = dwls_sample_stats_cx(Y, K, x)
    A  = dwlsmv_gammadiag(st, Y, x)
    par = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    if (metric == "marker") {
        fit = dwlsmv_fit_noretry(st, par)
    }
    else {
        fit = dwls_fit(st, par)
    }
    dwls_inference(fit)
    Bd = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    return((A[7, 2], A[6, 2], A[7, 3], Bd[4, 2], Bd[3, 1]))
}

// ---- 8: the conditional D11 trigger reads the bread row (Mata level, no fit): hand-built fit, A and Bd
real scalar td_trigger()
{
    struct dwlsmv_fit scalar fit
    real matrix A, Bd
    real scalar nf
    string scalar why

    nf = 0
    fit.converged = 1
    A  = J(6, 4, 0)
    // a converged, admissible fit with an unflagged bread: no trigger
    Bd = (J(4, 2, 0), J(4, 2, 0))
    why = dwlsmv_m11_rectrigger(fit, A, Bd)
    if (why != "") {
        nf++
        printf("{err}  8: unflagged bread fired the trigger: %s\n", why)
    }
    // bread row 3 = 1: fires, with the bread's cause and no other
    Bd[3, 1] = 1
    why = dwlsmv_m11_rectrigger(fit, A, Bd)
    if (strpos(why, "could not determine all parameters reliably") == 0) {
        nf++
        printf("{err}  8: Bd row 3 = 1 did not fire the trigger: '%s'\n", why)
    }
    if (strpos(why, "did not converge") > 0 | strpos(why, "was not admissible") > 0) {
        nf++
        printf("{err}  8: Bd row 3 = 1 alone also reports another cause: '%s'\n", why)
    }
    // bread row 2 (a zero column) alone does not fire (decided 2026-10-07)
    Bd[3, 1] = 0
    Bd[2, 1] = 1
    if (dwlsmv_m11_rectrigger(fit, A, Bd) != "") {
        nf++
        printf("{err}  8: Bd row 2 = 1 fired the trigger\n")
    }
    // an empty Bd (0 x 0) evaluates no bread term and must not abort
    if (dwlsmv_m11_rectrigger(fit, A, J(0, 0, .)) != "") {
        nf++
        printf("{err}  8: an empty Bd fired the trigger\n")
    }
    return(nf)
}

// ======================================================================================================================
// checks 9 and 10: drafted by a read-only agent 2026-11-15; executed by the main session the same day (StataMP 18), first run passed
// ======================================================================================================================

// ---- 9: the sparse design. Pooled and per-group expected lines, from the Stata-built table HC (rows as countcells')
// pooled count of (item j, category position c) from HC
real scalar td_pooled(real matrix HC, real scalar j, real scalar c)
{
    return(sum(select(HC[., 4], (HC[., 1] :== j) :& (HC[., 2] :== c))))
}

// the pooled "lowest or highest category" line dwlsmv_countdiag() prints: first strict minimum in (item, category) order
string scalar td_expect_pline(real matrix HC, real colvector K)
{
    real scalar j, c, n, bn, bj, bc

    bn = .
    bj = 0
    bc = 0
    for (j = 1; j <= rows(K); j++) {
        for (c = 1; c <= K[j]; c++) {
            if (c == 1 | c == K[j]) {
                n = td_pooled(HC, j, c)
                if (bn == . | n < bn) {
                    bn = n
                    bj = j
                    bc = c
                }
            }
        }
    }
    return("lowest or highest category: " + strofreal(bn) + " (y" + strofreal(bj) + ", " +
           (bc == 1 ? "lowest" : "highest") + ")")
}

// the sparse design's own premises and the exact cell, on the data in memory (K = 3, 3, 2); HC as above.
// The one thin cell is (y3 = 1, x = 1) with 2 observations while the pooled y3 = 1 count is large.
real scalar td_sparse(real matrix HC)
{
    real matrix Y, C
    real colvector K, x
    real scalar nf, r

    nf = 0
    td_data(Y, K, x)
    C = dwlsmv_countcells(Y, K, x)
    for (r = 1; r <= rows(C); r++) {
        if (C[r, 1] == 3 & C[r, 2] == 2 & C[r, 3] == 1) {
            if (C[r, 4] != 2) {
                nf++
                printf("{err}  9: countcells says %f in (y3 = 1, x = 1); the design has 2\n", C[r, 4])
            }
        }
    }
    if (td_pooled(HC, 3, 2) < 100) {
        nf++
        printf("{err}  9: premise failed: the pooled y3 = 1 count is %f, not large\n", td_pooled(HC, 3, 2))
    }
    return(nf)
}

// ---- 10: row 5 of dwlsmv_gammadiag() in conditional mode excludes nothing (G6)
// A hand-built Stage 1 result (identity or diagonal Gamma, a planted psibar) on a Y whose (y1, y2) pair is binary with
// EXACTLY ONE empty cell in the pooled table. Such a Y cannot come out of dwls_sample_stats_cx(): with K = 2 the item
// fits are saturated in each x-group, an empty pooled cell is empty in both groups, and the common-rho likelihood then
// increases to the bound (status 2, rc 3498; precheck_out.txt, part B). So the test is at the gammadiag level,
// on a struct built here, which gammadiag accepts (its checks are shape and codes). K = (2, 2, 3); T = 4, p = 3, P = 3.
// Conditional map: D = 10; thresholds 1-4, slopes 5-7, pairs 8-10 (y1~~y2 at 8). Ordinary map: D = 7, pairs 5-7.
real scalar td_g6()
{
    struct dwlsmv_stage1 scalar st1, st0
    struct dwlsmv_statmap scalar m1, m0
    real matrix Y, A, A0, G, tab
    real colvector K, x, cnt, pb, rep0, zk
    real scalar nf, c, k, D1, D0, xx, y1, y2, y3, yy, z, kz

    nf = 0
    // the design: 18 cells (x = c mod 2; y3 = floor(c/2) mod 3; (y1, y2) = (0,0), (0,1), (1,1) by floor(c/6)); N = 358
    cnt = (30 \ 22 \ 26 \ 28 \ 18 \ 24 \ 14 \ 20 \ 16 \ 12 \ 18 \ 16 \ 22 \ 16 \ 20 \ 24 \ 14 \ 18)
    Y = J(0, 3, .)
    x = J(0, 1, .)
    for (c = 0; c < 18; c++) {
        xx = mod(c, 2)
        y3 = mod(floor(c / 2), 3)
        yy = floor(c / 6)
        y1 = (yy == 2)
        y2 = (yy >= 1)
        Y = Y \ (J(cnt[c + 1], 1, 1) * (y1, y2, y3))
        x = x \ J(cnt[c + 1], 1, xx)
    }
    K = (2 \ 2 \ 3)

    // premises: the pooled (y1, y2) table has exactly one empty cell; the ordinary pair index says so
    tab = dwlsmv_crosstab(Y[., 1], Y[., 2], 2, 2)
    if (sum(tab :== 0) != 1) {
        nf++
        printf("{err}  10: premise failed: the pooled (y1,y2) table has %f empty cells, wanted 1\n", sum(tab :== 0))
    }
    m1 = dwlsmv_statmap_build(K, 1)
    m0 = dwlsmv_statmap_build(K, 0)
    D1 = m1.D
    D0 = m0.D
    if (D1 != 10 | D0 != 7 | m1.T != 4) {
        nf++
        printf("{err}  10: premise failed: D = %f (conditional), %f (ordinary), T = %f; wanted 10, 7, 4\n", D1, D0, m1.T)
        return(nf)
    }
    rep0 = dwlsmv_emptycell_pairs(Y, m0)
    if (any(rep0 :!= (1 \ 0 \ 0))) {
        nf++
        printf("{err}  10: premise failed: dwlsmv_emptycell_pairs() on the ordinary map is (%f, %f, %f), wanted (1, 0, 0)\n",
               rep0[1], rep0[2], rep0[3])
    }

    st1.map    = m1
    st1.N      = rows(Y)
    st1.s      = J(D1, 1, 0)
    st1.Gamma  = I(D1)
    st1.bound  = J(m1.P, 1, 0)
    st0.map    = m0
    st0.N      = rows(Y)
    st0.s      = J(D0, 1, 0)
    st0.Gamma  = I(D0)
    st0.bound  = J(m0.P, 1, 0)

    // (a) CONDITIONAL: a value planted at ANY position k is reported at k, unexcluded; Gamma = I, so z = 1e-3 exactly.
    // k = 6 is the slope of y2, the entry a mis-indexed exclusion (pos = T + k for k the pair index) zeroes on this Y;
    // k = 8 is the rho entry of the one-empty pair itself.
    for (k = 1; k <= D1; k++) {
        pb = J(D1, 1, 0)
        pb[k] = 1e-3
        st1.psibar = pb
        A = dwlsmv_gammadiag(st1, Y, x)
        if (A[5, 1] != 1 | A[5, 2] != 1e-3 | A[5, 3] != k | A[5, 4] != 0) {
            nf++
            printf("{err}  10: conditional, value planted at entry %f: row 5 is (%f, %g, %f, %f), wanted (1, 0.001, %f, 0)\n",
                   k, A[5, 1], A[5, 2], A[5, 3], A[5, 4], k)
        }
    }

    // (b) CONDITIONAL, unplanted structure: row 5 equals the direct maximum over ALL D entries, with its first location
    G = diag((1 \ 2 \ 3 \ 4 \ 5 \ 6 \ 7 \ 8 \ 9 \ 10))
    pb = 1e-4 * (1 \ 2 \ 3 \ 4 \ 5 \ 6 \ 7 \ 8 \ 9 \ 10)
    pb[6] = 5e-3
    st1.Gamma  = G
    st1.psibar = pb
    A  = dwlsmv_gammadiag(st1, Y, x)
    zk = abs(pb) :/ sqrt(diagonal(G))
    z  = max(zk)
    for (kz = 1; zk[kz] != z; kz++) {
    }
    if (A[5, 1] != (z > 1e-6) | abs(A[5, 2] - z) > 1e-14 * z | A[5, 3] != kz | A[5, 4] != 0) {
        nf++
        printf("{err}  10: conditional, direct maximum: row 5 is (%f, %g, %f, %f), wanted (%f, %g, %f, 0)\n",
               A[5, 1], A[5, 2], A[5, 3], A[5, 4], (z > 1e-6), z, kz)
    }

    // (c) ORDINARY control, same Y: the one-empty pair's rho entry (5) IS excluded, the next pair's (6) is not. This is the
    // live premise: the Y does contain an excludable pair, so (a) is a test of the conditional branch, not of nothing.
    st0.psibar = (0 \ 0 \ 0 \ 0 \ 1e-3 \ 0 \ 0)
    A0 = dwlsmv_gammadiag(st0, Y)
    if (A0[5, 1] != 0 | A0[5, 2] != 0) {
        nf++
        printf("{err}  10: ordinary control: the one-empty pair's rho entry was not excluded: row 5 is (%f, %g, %f, %f)\n",
               A0[5, 1], A0[5, 2], A0[5, 3], A0[5, 4])
    }
    st0.psibar = (0 \ 0 \ 0 \ 0 \ 0 \ 1e-3 \ 0)
    A0 = dwlsmv_gammadiag(st0, Y)
    if (A0[5, 1] != 1 | A0[5, 2] != 1e-3 | A0[5, 3] != 6) {
        nf++
        printf("{err}  10: ordinary control: a planted value at entry 6 is row 5 (%f, %g, %f), wanted (1, 0.001, 6)\n",
               A0[5, 1], A0[5, 2], A0[5, 3])
    }
    return(nf)
}

end

* ======================================================================================================
* Check 1 hand design, check 7, check 5 (case_D), then the two conditional fixtures
* ======================================================================================================
mata: st_local("nfh", strofreal(td_hand()))
td_cond `=(`nfh' == 0)' "1 hand design (`nfh' failures)"
mata: st_local("nft", strofreal(td_trigger()))
td_cond `=(`nft' == 0)' "8 conditional D11 trigger reads the bread row (`nft' failures)"
mata: st_local("nfs", strofreal(td_slab()))
td_cond `=(`nfs' == 0)' "7 slab labels (`nfs' failures)"

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: st_local("nfo", strofreal(td_ordinary()))
* check 4, r2: x handed to an unconditional Stage 1 result (case_D is in memory)
capture noisily mata: td_r2()
td_cond `=(_rc == 3498)' "case_D: refusal r2 (wanted rc 3498)"
td_cond `=(`nfo' == 0)' "5 ordinary regression (`nfo' failures)"

foreach t in A D {
    import delimited using "fixtures/data/condx_`t'.csv", clear varnames(1) delimiters(",")

    * ---- Stata's own counts: patterns, then the cell table in countcells' row order
    tempvar tx ty
    quietly egen `tx' = tag(y1 y2 y3 x)
    quietly count if `tx'
    local npx = r(N)
    quietly egen `ty' = tag(y1 y2 y3)
    quietly count if `ty'
    local npy = r(N)
    tempname HC
    matrix `HC' = J(2 + 2 * 8, 4, .)
    quietly count if x == 0
    matrix `HC'[1, 1] = (0, 0, 0, r(N))
    quietly count if x == 1
    matrix `HC'[2, 1] = (0, 0, 1, r(N))
    local n0 = `HC'[1, 4]
    local n1 = `HC'[2, 4]
    local r 2
    local Kv 3 3 2
    forvalues j = 1/3 {
        local Kj : word `j' of `Kv'
        forvalues c = 1/`Kj' {
            forvalues g = 0/1 {
                local ++r
                quietly count if y`j' == `c' - 1 & x == `g'
                matrix `HC'[`r', 1] = (`j', `c', `g', r(N))
            }
        }
    }
    mata: st_matrix("td_HC", st_matrix("`HC'"))
    mata: st_local("nfc", strofreal(td_check(`npx', `npy', st_matrix("td_HC"))))
    td_cond `=(`nfc' == 0)' "`t': checks 1-3 (`nfc' failures; patterns (y,x) `npx', y `npy')"

    * ---- 4 refusals, each asserted == 3498
    foreach k in 1 3 4 5 6 7 8 {
        capture noisily mata: td_r`k'()
        td_cond `=(_rc == 3498)' "`t': refusal r`k' (wanted rc 3498)"
    }

    * ---- 6 through the command
    tempfile lg
    td_run "`lg'" dwlsmv (F -> y1 y2 y3) (F <- x)
    td_cond `=(r(rc) == 0)' "`t': dwlsmv with a covariate returned nonzero rc"
    if r(rc) == 0 {
        mata: st_matrix("td_dir", td_direct(st_global("e(canon)"), st_global("e(metric)")))
        local ok = !missing(e(gamma_ratio)) & !missing(e(gamma_ratio_raw)) & !missing(e(gamma_rank))
        td_cond `ok' "`t': e(gamma_ratio), e(gamma_ratio_raw), e(gamma_rank) posted"
        if `ok' {
            td_cond `=abs(e(gamma_ratio) - td_dir[1, 1]) <= 1e-15 * abs(td_dir[1, 1])' "`t': e(gamma_ratio) `=e(gamma_ratio)' vs direct `=td_dir[1, 1]'"
            td_cond `=abs(e(gamma_ratio_raw) - td_dir[1, 2]) <= 1e-15 * abs(td_dir[1, 2])' "`t': e(gamma_ratio_raw) vs direct"
            td_cond `=(e(gamma_rank) == td_dir[1, 3])' "`t': e(gamma_rank) `=e(gamma_rank)' vs direct `=td_dir[1, 3]'"
        }
        * e(bread_ratio) is posted (as for an ordinary fit) and equals dwlsmv_breaddiag()'s row 4 of the same first fit;
        * a healthy fit: row 3 unflagged and the ratio far above 1e-8. The rerun is the same fit only if no recovery ran.
        td_cond `=!missing(e(bread_ratio))' "`t': e(bread_ratio) is not posted for a conditional fit"
        if "`e(recovery)'" == "none" & !missing(e(bread_ratio)) {
            td_cond `=abs(e(bread_ratio) - td_dir[1, 4]) <= 1e-12 * abs(td_dir[1, 4])' "`t': e(bread_ratio) `=e(bread_ratio)' vs direct breaddiag `=td_dir[1, 4]'"
            td_cond `=(td_dir[1, 5] == 0 & td_dir[1, 4] > 1e-4)' "`t': healthy bread: flag `=td_dir[1, 5]', ratio `=td_dir[1, 4]'"
        }
        else if "`e(recovery)'" != "none" {
            display as text "  note: `t' e(recovery) = `e(recovery)': the bread comparison is skipped (the fit was replaced)"
        }
        mata: st_numscalar("tdx", findexternal("dwlsmv_m11_x") == NULL)
        td_cond `=(scalar(tdx) == 1)' "`t': the external dwlsmv_m11_x is removed after the fit"

        td_has "`lg'" "Category counts within each value of x"
        td_cond `=(r(hit) == 1)' "`t': within-group heading missing"
        td_has "`lg'" "both groups together"
        td_cond `=(r(hit) == 1)' "`t': pooled heading does not say both groups together"
        td_has "`lg'" "group sizes:                x = 0: `n0', x = 1: `n1'"
        td_cond `=(r(hit) == 1)' "`t': group sizes line (Stata counts `n0', `n1')"
        mata: st_local("bl", td_expect_bline(st_matrix("td_HC"), (3 \ 3 \ 2)))
        td_has "`lg'" "`bl'"
        td_cond `=(r(hit) == 1)' "`t': lowest-or-highest group-cell line: `bl'"
        td_has "`lg'" "The standard error of a slope (the effect of x on an item) depends on the"
        td_cond `=(r(hit) == 1)' "`t': slope sentence missing"
        td_has "`lg'" "that move together on every observation. An item that is an exact function of another"
        td_cond `=(r(hit) == 1)' "`t': aliasing 'not computed' sentence missing"
        td_has "`lg'" "value near 1e-16 of e(gamma_ratio)"
        td_cond `=(r(hit) == 1)' "`t': aliasing limitation wording missing"
        td_has "`lg'" "when Gamma has lost rank in this way has not been established."
        td_cond `=(r(hit) == 1)' "`t': the sentence on the unestablished accuracy under a rank-lost Gamma is missing"
        td_has "`lg'" "Too few distinct response patterns"
        td_cond `=(r(hit) == 0)' "`t': a 'too few patterns' note on a healthy fit"
        td_has "`lg'" "Not computed for a model with a covariate in this version: the check for thresholds"
        td_cond `=(r(hit) == 1)' "`t': the 'not computed' sentence names the aliasing check"
        td_has "`lg'" "weight and derivative matrices"
        td_cond `=(r(hit) == 0)' "`t': the bread is computed now; the old 'not computed' wording must be gone"
    }
}

* ======================================================================================================================
* checks 9 and 10: drafted by a read-only agent 2026-11-15; executed by the main session the same day (StataMP 18), first run passed
* ======================================================================================================================
* ---- 9: a sparse conditional design, built here from condx_A (not a committed fixture): the (y3 = 1, x = 1) cell is cut
* to its first two rows in file order (the other 312 become y3 = 0); everything else is condx_A. The pooled y3 = 1 count
* stays large, so the within-group "lowest or highest category" line names a cell the pooled line does not.
* (i) the printed lines; (ii) check 2's comparison (checks 1-3 via td_check) on it.
import delimited using "fixtures/data/condx_A.csv", clear varnames(1) delimiters(",")
quietly generate byte tdhit = (x == 1 & y3 == 1)
quietly generate long tdcum = sum(tdhit)
quietly replace y3 = 0 if tdhit & tdcum > 2
quietly drop tdhit tdcum
quietly count if x == 1 & y3 == 1
td_cond `=(r(N) == 2)' "S: premise, the thin cell has `=r(N)' observations, wanted 2"
quietly count if x == 0
local n0 = r(N)
quietly count if x == 1
local n1 = r(N)

tempvar tx ty
quietly egen `tx' = tag(y1 y2 y3 x)
quietly count if `tx'
local npx = r(N)
quietly egen `ty' = tag(y1 y2 y3)
quietly count if `ty'
local npy = r(N)
tempname HS
matrix `HS' = J(2 + 2 * 8, 4, .)
matrix `HS'[1, 1] = (0, 0, 0, `n0')
matrix `HS'[2, 1] = (0, 0, 1, `n1')
local r 2
local Kv 3 3 2
forvalues j = 1/3 {
    local Kj : word `j' of `Kv'
    forvalues c = 1/`Kj' {
        forvalues g = 0/1 {
            local ++r
            quietly count if y`j' == `c' - 1 & x == `g'
            matrix `HS'[`r', 1] = (`j', `c', `g', r(N))
        }
    }
}
mata: st_matrix("td_HS", st_matrix("`HS'"))
* checks 1-3 (countcells against -count- cell by cell, the sums, gammadiag rows) on the sparse data
mata: st_local("nfc", strofreal(td_check(`npx', `npy', st_matrix("td_HS"))))
td_cond `=(`nfc' == 0)' "S: checks 1-3 on the sparse design (`nfc' failures; patterns (y,x) `npx', y `npy')"
mata: st_local("nfs2", strofreal(td_sparse(st_matrix("td_HS"))))
td_cond `=(`nfs2' == 0)' "S: the thin cell (`nfs2' failures)"
mata: st_local("bl", td_expect_bline(st_matrix("td_HS"), (3 \ 3 \ 2)))
mata: st_local("pl", td_expect_pline(st_matrix("td_HS"), (3 \ 3 \ 2)))
td_cond `=("`bl'" != "`pl'")' "S: the within-group line equals the pooled line: `bl'"
td_cond `=(strpos("`bl'", ": 2 (y3, highest, x = 1)") > 0)' "S: expected within-group line is not the thin cell: `bl'"

tempfile lgs
td_run "`lgs'" dwlsmv (F -> y1 y2 y3) (F <- x)
td_cond `=(r(rc) == 0)' "S: dwlsmv on the sparse design returned nonzero rc"
if r(rc) == 0 {
    td_has "`lgs'" "Category counts within each value of x"
    td_cond `=(r(hit) == 1)' "S: within-group heading missing"
    td_has "`lgs'" "group sizes:                x = 0: `n0', x = 1: `n1'"
    td_cond `=(r(hit) == 1)' "S: group sizes line (Stata counts `n0', `n1')"
    td_has "`lgs'" "`bl'"
    td_cond `=(r(hit) == 1)' "S: within-group line missing: `bl'"
    td_has "`lgs'" "`pl'"
    td_cond `=(r(hit) == 1)' "S: pooled line missing or changed: `pl'"

    * ---- 9R (added 2026-11-15, the owner's review): this real fit REACHES the conditional bread cause, with the other two.
    * Recorded on its first run (StataMP 18, default processors): the marker fit diverges (lambda_3 about 3028, var(e.F) about
    * 1e-4), the bread has rank 8 of 9 at the last iterate, all three trigger causes fire, recovery is attempted and fails.
    * Asserted, not printed: a change in any of these is a FINDING to read, not a typo (the iterate is a nonconverged one).
    local rk = `"(the derivative matrix has rank 8 of 9 columns at the"'
    td_cond `=(e(converged) == 0)' "9R: the sparse marker fit converged (e(converged) = `e(converged)'); the recorded fit did not"
    td_cond `=("`e(recovery)'" == "failed")' "9R: e(recovery) is '`e(recovery)'', recorded 'failed'"
    td_cond `=(!missing(e(bread_ratio)) & e(bread_ratio) < 1e-8)' "9R: e(bread_ratio) = `=e(bread_ratio)' is not below the 1e-8 criterion"
    td_has "`lgs'" "`rk'"
    td_cond `=(r(hit) == 1)' "9R: the bread note does not say rank 8 of 9"
    td_has "`lgs'" "did not converge, was not admissible and could not determine all parameters reliably (the information in the fit was insufficient);"
    td_cond `=(r(hit) == 1)' "9R: the recovery note does not name all three trigger causes"
    td_has "`lgs'" "Recovery attempt failed; first fit retained:"
    td_cond `=(r(hit) == 1)' "9R: the recovery note does not say the attempt failed"
    * the failed recovery exposes no standard errors and no test (spec D8, Exception): nothing posted, nothing printed
    mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
    td_cond `=(`nv' == 0)' "9R: e(V) is posted (`nv' rows) for a nonconverged conditional fit"
    td_cond `=("`e(T3)'" == "" & "`e(p_T3)'" == "")' "9R: e(T3) or e(p_T3) is posted ('`e(T3)'', '`e(p_T3)'')"
    td_has "`lgs'" "No standard errors, T3 or p-value are reported for them."
    td_cond `=(r(hit) == 1)' "9R: the header does not say no standard errors, T3 or p-value are reported"
    td_has "`lgs'" "Scaled-and-shifted test (T3): not reported (not a solution)"
    td_cond `=(r(hit) == 1)' "9R: the T3 line is not 'not reported'"
    td_has "`lgs'" "Std. err."
    td_cond `=(r(hit) == 0)' "9R: a standard-error column is printed"
    td_has "`lgs'" "Prob > chi2"
    td_cond `=(r(hit) == 0)' "9R: a p-value is printed"
}

* ---- 10: gammadiag row 5 in conditional mode excludes no entry (mutant G6); Mata-level, hand-built Stage 1 result
mata: st_local("nfg", strofreal(td_g6()))
td_cond `=(`nfg' == 0)' "10 gammadiag row 5, conditional mode excludes nothing (`nfg' failures)"

* ---- 6 negative control: an ordinary fit
import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
tempfile lg0
td_run "`lg0'" dwlsmv (F -> y1 y2 y3 y4)
td_cond `=(r(rc) == 0)' "control: ordinary fit returned nonzero rc"
if r(rc) == 0 {
    td_cond `=!missing(e(bread_ratio)) & !missing(e(gamma_ratio))' "control: ordinary fit posts e(bread_ratio) and e(gamma_ratio)"
    td_has "`lg0'" "Category counts (the fewest observations in any one category):"
    td_cond `=(r(hit) == 1)' "control: the ordinary pooled heading changed"
    td_has "`lg0'" "within each value of"
    td_cond `=(r(hit) == 0)' "control: an ordinary fit prints a within-group block"
    td_has "`lg0'" "both groups together"
    td_cond `=(r(hit) == 0)' "control: an ordinary fit's heading says both groups together"
}

if ${TD_NFAIL} {
    display as error "test_condx_diag: ${TD_NFAIL} check(s) failed"
    exit 9
}
display as text "test_condx_diag: all checks passed"
exit 0
