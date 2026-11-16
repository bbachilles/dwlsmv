*! Executed by the main session (Stata 18 MP) and mutant-tested, 2026-11-14 at the latest.
*! Conditional mode, output path through the command: dwlsmv with a binary 0/1 covariate
*! ((F <- x) and (y <- x)) runs end to end, displays and posts what decision (a) of 2026-11-13
*! allows, and nothing it does not. Drafted by a read-only agent 2026-11-14; remove this banner
*! only after the file has been run and what broke has been fixed.
*!
*! Needs the input side too (translator accepts the covariate paths; dwlsmv.ado routes Stage 1 to
*! dwls_sample_stats_cx and passes the covariate name on as the local cov) and the edits of
*! edits.md (header, posting, refusals). REVISED 2026-11-14 with the input-side draft: the skip
*! (exit 8 while the translator refused a covariate) is removed, so that a regression of the
*! input side fails here instead of being skipped; check 9 is now the command's own refusal of a
*! covariate not coded 0/1 (rc 198, P5), since the command checks before Stage 1.
*!
*! Data (generated here, fixed seed, no fixture files): N = 1500, x = Bernoulli(.5), one factor
*! with loadings .7 and conditional residual variance .51 (Var(y*|x) = 1), the factor shifted .4
*! by x, four items with thresholds -.5 and .5 (K = 3 each), and a direct effect of +.3 of x on
*! y4*. D = T + p + P = 8 + 4 + 6 = 18 statistics. Expected degrees of freedom, derived by hand
*! and not read from the code: model A, (F -> y1 y2 y3 y4) (F <- x): free parameters 8
*! thresholds, 3 loadings (marker), gamma, the disturbance variance = 13, so df = 5. Model B adds
*! (y4 <- x): 14 parameters, df = 4. metric(unitvariance) trades the marker for a fixed variance
*! and a fourth loading, the same count.
*!
*! Checks
*!   1  model A runs (rc 0), converged, e(N) = 1500, e(covariates) = "x".
*!   2  e(T3), e(df_T3), e(p_T3) posted and finite, df = 5; p = chi2tail(df, T3).
*!   3  e(b): thresholds first, then loadings, then F:x, under sem's names.
*!   4  no fit index or baseline is posted, omitted and not missing: e(rmsea_s), e(cfi_s),
*!      e(tli_s), e(T3_b), e(df_T3_b) read as empty strings (a posted missing reads "."), no
*!      e(bound); Module 10's e(gamma_ratio), e(gamma_ratio_raw), e(gamma_rank) and e(bread_ratio) ARE posted from
*!      2026-11-15 (test_condx_diag.do).
*!   5  the printed output: the header names the covariate, the note that no conditional CFI/TLI
*!      is defined is present, no "RMSEA (scaled)" / "CFI (scaled)" line is printed; the same on
*!      replay.
*!   2b model B with metric(unitvariance): F:x and y4:x rows, df = 4, no index.
*!   6  negative control: an unconditional fit still posts and prints the indices and e(bound).
*!   7  difftest on a conditional fit refuses at rc 3498, for the covariate and not another reason.
*!   8  standardized with a covariate is refused at rc 198, at estimation and on replay, and the
*!      fit stays posted after the refused replay.
*!   9  x coded 0/2 is refused at rc 198 by the command, before Stage 1 (P5, provisional); the
*!      message must say 0 and 1. Stage 1's own rc 3498 for the same coding is covered by
*!      test_condx_stage1.do and, for the command, test_condx_input.do check 7 (separation).
*!
*! WHAT THIS CERTIFIES AND DOES NOT. It certifies the output path (display, e(), refusals) of a
*! conditional fit. It says nothing about the estimates: the tier-1 reference check of fit and
*! inference on data that do not fit exactly is a separate required item (decisions 2026-11-13).
*! Model A is misspecified on purpose-free data (y4 carries a direct effect), so T3 is not
*! asserted; only the count of degrees of freedom is.
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

local nfail 0

* run a command capturing what it prints: r(rc), r(found) for the first phrase, r(found2) for
* the second
capture program drop tcx_run
program tcx_run, rclass
    gettoken phrase 0 : 0
    gettoken phrase2 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tcxlog)
    capture noisily `0'
    local rc = _rc
    quietly log close tcxlog
    mata: st_local("txt", invtokens(cat(st_local("lg"))', " "))
    mata: st_local("found",  strofreal(strpos(st_local("txt"), st_local("phrase"))  > 0))
    mata: st_local("found2", strofreal(strpos(st_local("txt"), st_local("phrase2")) > 0))
    return scalar rc     = `rc'
    return scalar found  = `found'
    return scalar found2 = `found2'
end

* ---- the data ---------------------------------------------------------------------------
set seed 20261114
quietly set obs 1500
quietly generate byte x = runiform() < .5
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
drop eta z1 z2 z3 z4

local mA "(F -> y1 y2 y3 y4) (F <- x)"
local mB "(F -> y1 y2 y3 y4) (F <- x) (y4 <- x)"

* ---- 1. model A runs, and the first print is captured for check 5 -----------------------
tcx_run "Conditional on covariate: x" "no CFI or TLI" dwlsmv `mA'
local rcA   = r(rc)
local hdrA  = r(found)
local noteA = r(found2)
if `rcA' {
    display as error " 1. model A: rc `rcA'"
    local ++nfail
}
else {
    if e(converged) != 1 | e(N) != 1500 | "`e(covariates)'" != "x" {
        display as error " 1. converged " e(converged) ", N " e(N) ", e(covariates) '`e(covariates)''"
        local ++nfail
    }
    display as text " 1. model A runs: converged " e(converged) ", N " e(N) ", covariates `e(covariates)'"
}

* ---- 2. statistic, df, p ----------------------------------------------------------------
if !`rcA' {
    local pwant = chi2tail(e(df_T3), e(T3))
    if missing(e(T3)) | missing(e(df_T3)) | missing(e(p_T3)) | e(df_T3) != 5 | abs(e(p_T3) - `pwant') > 1e-12 {
        display as error " 2. T3 " e(T3) " df " e(df_T3) " p " e(p_T3) " (want df 5, p " `pwant' ")"
        local ++nfail
    }
    else display as text " 2. T3 " %9.3f e(T3) ", df 5, p " %6.4f e(p_T3)
}

* ---- 3. e(b) names: thresholds, loadings, gamma ------------------------------------------
if !`rcA' {
    local names : colfullnames e(b)
    local want "/y1:cut1 /y1:cut2 /y2:cut1 /y2:cut2 /y3:cut1 /y3:cut2 /y4:cut1 /y4:cut2 y1:F y2:F y3:F y4:F F:x /var(e.F)"
    if "`names'" != "`want'" {
        display as error " 3. model A e(b) names: `names'"
        display as error "    wanted:              `want'"
        local ++nfail
    }
    else display as text " 3. model A e(b) names are as wanted"
}

* ---- 4. omitted, not missing ---------------------------------------------------------------
if !`rcA' {
    local bad
    foreach s in rmsea_s cfi_s tli_s T3_b df_T3_b {
        if "`e(`s')'" != "" local bad `bad' e(`s')
    }
    mata: st_local("nb", strofreal(rows(st_matrix("e(bound)"))))
    if "`bad'" != "" | `nb' != 0 {
        display as error " 4. posted but must be omitted: `bad'; rows of e(bound) = `nb'"
        local ++nfail
    }
    else display as text " 4. no fit index, baseline or e(bound) is posted"
}

* ---- 5. the printed output ------------------------------------------------------------------
if !`rcA' {
    if !`hdrA' | !`noteA' {
        display as error " 5. header names the covariate: `hdrA'; the no-CFI/TLI note: `noteA'"
        local ++nfail
    }
    tcx_run "RMSEA (scaled)" "CFI (scaled)" dwlsmv
    local rrc = r(rc)
    local rm  = r(found)
    local cf  = r(found2)
    tcx_run "Conditional on covariate: x" "no CFI or TLI" dwlsmv
    if `rrc' | `rm' | `cf' | !r(found) | !r(found2) {
        display as error " 5. replay: rc `rrc', RMSEA line `rm', CFI line `cf'; header " r(found) ", note " r(found2)
        local ++nfail
    }
    tcx_run "RMSEA (scaled)" "CFI (scaled)" dwlsmv `mA'
    if r(found) | r(found2) {
        display as error " 5. an RMSEA or CFI line is printed for a conditional fit"
        local ++nfail
    }
    else display as text " 5. header, note and replay as wanted; no index line"
}

* ---- 2b. model B, unitvariance, kappa row -------------------------------------------------
capture quietly dwlsmv `mB', metric(unitvariance)
local rcB = _rc
if `rcB' {
    display as error " 2b. model B, metric(unitvariance): rc `rcB'"
    local ++nfail
}
else {
    local namesB : colfullnames e(b)
    local okB = strpos(" `namesB'", " F:x") > 0 & strpos(" `namesB'", " y4:x") > 0
    if !`okB' | e(df_T3) != 4 | "`e(rmsea_s)'" != "" {
        display as error " 2b. model B names `namesB'; df " e(df_T3) "; rmsea_s '`e(rmsea_s)''"
        local ++nfail
    }
    else display as text " 2b. model B: F:x and y4:x rows, df 4, no index"
    * diagnostic only, nothing asserted: the generating gamma is .4 and kappa .3
    display as text "     (diagnostic) gamma " %6.3f _b[F:x] ", kappa(y4) " %6.3f _b[y4:x] ", T3 " %7.3f e(T3)
}

* ---- 6. negative control: an unconditional fit keeps its indices ---------------------------
quietly dwlsmv (F -> y1 y2 y3 y4)
mata: st_local("nb0", strofreal(rows(st_matrix("e(bound)"))))
if "`e(rmsea_s)'" == "" | "`e(cfi_s)'" == "" | "`e(covariates)'" != "" | `nb0' == 0 {
    display as error " 6. unconditional fit: rmsea_s '`e(rmsea_s)'', cfi_s '`e(cfi_s)'', covariates '`e(covariates)'', bound rows `nb0'"
    local ++nfail
}
else display as text " 6. an unconditional fit still posts the indices and e(bound), no e(covariates)"
tcx_run "RMSEA (scaled)" "Conditional on covariate" dwlsmv
if !r(found) | r(found2) {
    display as error " 6. unconditional replay: index line " r(found) ", covariate header " r(found2)
    local ++nfail
}

* ---- 7. difftest refuses a conditional fit --------------------------------------------------
quietly dwlsmv `mA'
quietly estimates store cx1
tcx_run "covariate" "difftest" difftest cx1
if r(rc) != 3498 | !r(found) | !r(found2) {
    display as error " 7. difftest on a conditional fit: rc " r(rc) ", names the covariate " r(found) ", says difftest " r(found2)
    local ++nfail
}
else display as text " 7. difftest refuses a conditional fit at rc 3498"
capture estimates drop cx1

* ---- 8. standardized is refused ---------------------------------------------------------------
tcx_run "standardized" "covariate" dwlsmv `mA', standardized
if r(rc) != 198 | !r(found) | !r(found2) {
    display as error " 8. standardized at estimation: rc " r(rc) ", message " r(found) r(found2)
    local ++nfail
}
quietly dwlsmv `mA'
tcx_run "standardized" "covariate" dwlsmv, standardized
local rr = r(rc)
if `rr' != 198 | "`e(covariates)'" != "x" | e(df_T3) != 5 {
    display as error " 8. standardized replay: rc `rr'; fit still posted: covariates '`e(covariates)'', df " e(df_T3)
    local ++nfail
}
else display as text " 8. standardized is refused at rc 198, at estimation and on replay"

* ---- 9. x coded 0/2: the command's refusal, rc 198, before Stage 1 (P5) -------------------------
preserve
quietly replace x = 2 * x
tcx_run "0 and 1" "covariate" dwlsmv `mA'
local rc9 = r(rc)
local f9  = r(found)
restore
if `rc9' != 198 | !`f9' {
    display as error " 9. x coded 0/2: rc `rc9' (want 198); message says 0 and 1: `f9'"
    local ++nfail
}
else display as text " 9. x coded 0/2 is refused at rc 198"

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_condx_output: `nfail' failure(s)"
    exit 9
}
display as result "test_condx_output: a conditional fit runs from the command, posts T3, df and p and no fit index or baseline, names the covariate and the missing CFI/TLI in the header and on replay, and difftest and standardized refuse it"
exit 0
