*! REVISED 2026-11-14 on the owner's decision (do not change the optimizer; keep check 5 as a characterization
*! fixture and document the limitation): check 5 is rewritten and check 5b added; the revision is NOT YET EXECUTED.
*! Executed 2026-11-14 by the main session (Stata 18 MP). Checks 1 and 4 expected no e(Cns) row for the marker
*! loading and were corrected (2 rows and 1 row); check 5 found a local optimum and no longer asserts the sign
*! of the loadings at gamma = -.4 (open finding, now check 5's characterization). Formerly: NOT YET EXECUTED. Drafted by a read-only agent 2026-11-14; target path tests/unit/test_condx_fixed.do.
*! Remove this banner only after the file has been run and what broke has been fixed.
*!
*! What a fixed value on a covariate path does, through the command: (F <- x@c) is gamma fixed at c,
*! (y <- x@c) is kappa fixed at c (spec 4.8, Conditional covariates), and (F <- x@0) is dropped by sem
*! and so writes no row at all. Pattern, loading and helper as in test_condx_input.do. No fixture files:
*! the data are simulated here with a fixed seed (the same design as test_condx_input.do).
*!
*! Hand counts (T 8 thresholds, p 4 slopes, P 6 correlations, D = 18). Model A = (F -> y1 y2 y3 y4)
*! (F <- x): free 8 thresholds + 3 loadings (marker) + gamma + disturbance variance = 13, df 5. Gamma
*! fixed: 12, df 6. Model B adds (y4 <- x): 14, df 4; with kappa fixed: 13, df 5.
*!
*! Checks (surface refusals rc 198; Module 3's refusal rc 3498)
*!   1  (F <- x@.4): _b[F:x] == .4 exactly, SE 0 or missing, e(Cns) has the row, canon has `gamma F x fixed .4`,
*!      df is the free model's + 1 (and 6), e(rank) the free model's - 1.
*!   2  (y4 <- x@-.3): the same with kappa and the value -.3 (the sign survives); df 5 against 4.
*!   3  gamma fixed at its own free-fit estimate: every other estimate equals the free fit's. The raw
*!      statistic T_raw is NOT posted by the command (only e(T3), whose scale and shift depend on df), so
*!      the equality of T_raw is implied by equal estimates (F is a function of theta) and not tested.
*!   4  (F <- x@0): runs, no F:x column and no e(Cns), df = free model + 1, canon and estimates identical to
*!      the same model written (F -> ...) (y4 <- x), cov(x*F@0), where sem leaves F exogenous (4a). 4b: with
*!      EVERY covariate path @0, the covariate is still in the model (a probe of the defect D-1 of the
*!      draft's atc_facts.md: if sem then omits x, the fit would be unconditional).
*!   5  CHARACTERIZATION of a known limitation, not desired behaviour: gamma fixed at -.4 under
*!      metric(unitvariance) converges to a local optimum (all loadings positive, T3 far above the +.4 fit),
*!      although the mirrored solution would fit as well. Pins what happens today; fails, on purpose, the day the
*!      limitation is lifted. 5b: under the marker metric the marker loading is fixed at 1 (the way out).
*!   6  free gamma with every kappa free: rc 3498 (screen N). The same model with gamma fixed runs; with
*!      every kappa fixed at a nonzero value (not free) and gamma free it runs too.
*!
*! TOLERANCES. Check 3: |difference| <= 1e-6 absolute on each estimate. Source: the optimizer stops at
*! ptol 1e-10 relative, vtol 1e-12 (dwlsmv_minimize.mata:120-122), and test_condx_fit.do check 6 accepts
*! 1e-9 relative between fits (marked "(agent)" there); 1e-6 is the T3 tier's absolute value, looser than
*! both. NO MEASURED SOURCE for two different start points reaching the same optimum: FLAGGED FOR THE
*! OWNER; the differences are printed, so the bound can be set from the first run, not before.
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

* run a command capturing what it prints: r(rc), r(found) for the first phrase, r(found2) for the second
capture program drop tf_run
program tf_run, rclass
    gettoken phrase 0 : 0
    gettoken phrase2 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(tflog)
    capture noisily `0'
    local rc = _rc
    quietly log close tflog
    mata: st_local("txt", invtokens(cat(st_local("lg"))', " "))
    mata: st_local("found",  strofreal(strpos(st_local("txt"), st_local("phrase"))  > 0))
    mata: st_local("found2", strofreal(strpos(st_local("txt"), st_local("phrase2")) > 0))
    return scalar rc     = `rc'
    return scalar found  = `found'
    return scalar found2 = `found2'
end

* ---- the data (the design of test_condx_input.do) -----------------------------------------------
set seed 20261114
quietly set obs 1000
quietly generate byte x = runiform() < .5
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
drop eta z1 z2 z3 z4

* ---- the free comparators, kept as matrices ------------------------------------------------------
* A: (F <- x) free.  B: (F <- x) (y4 <- x) free.
quietly dwlsmv (F -> y1 y2 y3 y4) (F <- x)
local dfA   = e(df_T3)
local rankA = e(rank)
matrix bA = e(b)
quietly dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- x)
local dfB   = e(df_T3)
local rankB = e(rank)
matrix bB = e(b)
local gB    = _b[F:x]

* ---- 1. gamma fixed at .4 -----------------------------------------------------------------------------
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@.4)
local rc1 = _rc
if `rc1' {
    display as error " 1. (F <- x@.4): rc `rc1'"
    local ++nfail
}
else {
    local canon `"`e(canon)'"'
    mata: st_local("hasf", strofreal(strpos(st_local("canon"), "gamma F x fixed .4") > 0))
    capture confirm matrix e(Cns)
    local hascns = (_rc == 0)
    local ncns   = cond(`hascns', rowsof(e(Cns)), 0)
    local se1    = _se[F:x]
    if _b[F:x] != .4 | !`hasf' | !`hascns' | `ncns' != 2 | !(`se1' == 0 | missing(`se1')) {
        display as error " 1. _b " _b[F:x] ", canon has fixed .4: `hasf', e(Cns) rows `ncns', SE `se1'"
        display as error "    canon: `canon'"
        local ++nfail
    }
    else if e(df_T3) != `dfA' + 1 | e(rank) != `rankA' - 1 | e(df_T3) != 6 | e(rank) != 12 {
        display as error " 1. df " e(df_T3) " (free model " `dfA' ", want +1 and 6), rank " e(rank) " (free " `rankA' ", want -1 and 12)"
        local ++nfail
    }
    else display as text " 1. gamma fixed at .4: exact in e(b), SE 0, in e(Cns) and e(canon), df +1, rank -1"
}

* ---- 2. kappa fixed at -.3 ----------------------------------------------------------------------------
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- x@-.3)
local rc2 = _rc
if `rc2' {
    display as error " 2. (y4 <- x@-.3): rc `rc2'"
    local ++nfail
}
else {
    local canon `"`e(canon)'"'
    mata: st_local("hask", strofreal(strpos(st_local("canon"), "kappa y4 x fixed -.3") > 0))
    capture confirm matrix e(Cns)
    local hascns = (_rc == 0)
    local se2    = _se[y4:x]
    if _b[y4:x] != -.3 | !`hask' | !`hascns' | !(`se2' == 0 | missing(`se2')) {
        display as error " 2. _b " _b[y4:x] ", canon has fixed -.3: `hask', e(Cns) `hascns', SE `se2'"
        display as error "    canon: `canon'"
        local ++nfail
    }
    else if e(df_T3) != `dfB' + 1 | e(rank) != `rankB' - 1 | e(df_T3) != 5 {
        display as error " 2. df " e(df_T3) " (free model " `dfB' ", want +1 and 5), rank " e(rank) " (free " `rankB' ")"
        local ++nfail
    }
    else display as text " 2. kappa fixed at -.3: exact (sign kept) in e(b), SE 0, in e(Cns) and e(canon), df +1, rank -1"
}

* ---- 3. gamma fixed at its own free-fit estimate --------------------------------------------------------
* model B's other estimates, with gamma fixed at B's own estimate (17 digits round-trip a double)
local gtxt = string(`gB', "%21.17g")
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@`gtxt') (y4 <- x)
local rc3 = _rc
if `rc3' {
    display as error " 3. gamma fixed at the free estimate: rc `rc3'"
    local ++nfail
}
else {
    matrix b3 = e(b)
    local maxd 0
    local worst ""
    local nc = colsof(b3)
    forvalues i = 1/`nc' {
        local nm : word `i' of `: colfullnames b3'
        local d = abs(b3[1, `i'] - bB[1, `i'])
        if `d' > `maxd' {
            local maxd = `d'
            local worst "`nm'"
        }
    }
    display as text "    3. largest |difference| to the free fit: " %9.2e `maxd' " at `worst'"
    if colsof(b3) != colsof(bB) | `maxd' > 1e-6 | reldif(_b[F:x], `gB') > 1e-14 | e(df_T3) != `dfB' + 1 {
        display as error " 3. other estimates differ from the free fit by " `maxd' " (bound 1e-6, a documented but unmeasured bound), gamma " _b[F:x] " vs " `gB' ", df " e(df_T3)
        local ++nfail
    }
    else display as text " 3. gamma fixed at its free estimate reproduces the other estimates (max diff above), df + 1"
}

* ---- 4. (F <- x@0) ------------------------------------------------------------------------------------------
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@0) (y4 <- x)
local rc4 = _rc
if `rc4' {
    display as error " 4. (F <- x@0): rc `rc4'"
    local ++nfail
}
else {
    local canon4 `"`e(canon)'"'
    matrix b4 = e(b)
    local df4 = e(df_T3)
    capture confirm matrix e(Cns)
    local hascns = (_rc == 0)
    capture local junk = _b[F:x]
    local hasFx = (_rc == 0)
    mata: st_local("nogam", strofreal(strpos(st_local("canon4"), "gamma") == 0))
    if `hasFx' | rowsof(e(Cns)) != 1 | !`nogam' | `df4' != `dfB' + 1 {
        display as error " 4. F:x in e(b) `hasFx', e(Cns) `hascns', no gamma line `nogam', df `df4' (want free model B " `dfB' " + 1)"
        display as error "    canon: `canon4'"
        local ++nfail
    }
    else {
        * 4a. the same model by the other route: F exogenous, its covariance with x fixed at 0
        capture noisily dwlsmv (F -> y1 y2 y3 y4) (y4 <- x), cov(x*F@0)
        local rc4a = _rc
        if `rc4a' {
            * if the translator refuses here instead, assert the refusal and document the asymmetry
            display as error " 4a. cov(x*F@0) route: rc `rc4a' (want 0; if 198, assert the refusal and record the asymmetry)"
            local ++nfail
        }
        else {
            local same = (`"`e(canon)'"' == `"`canon4'"')
            matrix b4a = e(b)
            local maxd 0
            forvalues i = 1/`=colsof(b4)' {
                local maxd = max(`maxd', abs(b4[1, `i'] - b4a[1, `i']))
            }
            if !`same' | colsof(b4a) != colsof(b4) | `maxd' > 1e-12 {
                display as error " 4a. canon identical `same', max |diff| of estimates " `maxd'
                local ++nfail
            }
            else display as text " 4. (F <- x@0): no F:x, no e(Cns), df free + 1; canon and estimates identical to cov(x*F@0) (max diff " %9.2e `maxd' ")"
        }
    }
}

* 4b. a probe: every covariate path @0. sem may then not list x at all (e(oxvars) empty).
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@0)
local rc4b = _rc
if `rc4b' {
    display as error " 4b. every path @0: rc `rc4b'"
    local ++nfail
}
else {
    local canon4b `"`e(canon)'"'
    mata: st_local("hascv", strofreal(strpos(st_local("canon4b"), "covariate x") > 0))
    * with x in the model: D = 18 statistics, 12 free parameters (8 thresholds, 3 loadings, psi), df 6;
    * without it the fit is unconditional of x: D = 14, df 14 - 12 = 2
    if !`hascv' | e(df_T3) != 6 {
        display as error " 4b. every covariate path @0: covariate line `hascv', df " e(df_T3) " (want x in the model, df 6). If df is 2 the covariate was dropped silently."
        local ++nfail
    }
    else display as text " 4b. every covariate path @0 still conditions on x (df 6)"
}

* ---- 5. CHARACTERIZATION of a known limitation (NOT a statement of desired behaviour) -------------------------
* Mirror algebra (spec 4.8, "Factor signs"): reversing the factor sends (lambda, gamma) to (-lambda, -gamma) and
* leaves thresholds, psi = 1 and b = lambda * gamma, hence sigma(theta), unchanged. So the GLOBAL optimum at
* gamma = -.4 has all loadings negative and the same fit as the +.4 model. What the command does TODAY (main
* session, 2026-11-14, Stata 18 MP, this file's data; true effect positive): at -.4 it reports convergence with all
* loadings POSITIVE (largest .68) and a T3 323 above the +.4 fit, i.e. it stops at a local optimum. The owner
* decided on 2026-11-14 not to change the optimizer on the strength of this one test; the limitation is documented
* (dwlsmv.sthlp, spec 4.8) and this check pins the behaviour. The bound of 10 on the T3 gap is COARSE ON PURPOSE and
* was CHOSEN AFTER MEASURING the gap of 323: it separates "a local optimum" from "the mirrored optimum" (gap near 0),
* it is not a tolerance on any comparison and is not derived.
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@.4), metric(unitvariance)
local rc5p = _rc
if !`rc5p' {
    local T3p   = e(T3)
    local convp = e(converged)
    local lpmin = min(_b[y1:F], _b[y2:F], _b[y3:F], _b[y4:F])
}
capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@-.4), metric(unitvariance)
local rc5 = _rc
if `rc5' | `rc5p' {
    display as error " 5. unitvariance with gamma fixed at -.4 / .4: rc `rc5' / `rc5p'"
    local ++nfail
}
else {
    local conv  = e(converged)
    local T3m   = e(T3)
    local lmin  = min(_b[y1:F], _b[y2:F], _b[y3:F], _b[y4:F])
    local lmax  = max(_b[y1:F], _b[y2:F], _b[y3:F], _b[y4:F])
    local gap   = `T3m' - `T3p'
    local canon `"`e(canon)'"'
    mata: st_local("hasn", strofreal(strpos(st_local("canon"), "gamma F x fixed -.4") > 0))
    display as text "    5. CHARACTERIZATION: T3(-.4) - T3(+.4) = " %9.3f `gap' "; loadings at -.4 from " %6.3f `lmin' " to " %6.3f `lmax' "; at +.4 smallest " %6.3f `lpmin'
    display as text "    5. If this check ever fails because the loadings came out NEGATIVE with an equal T3, the limitation has been lifted:"
    display as text "       update this check, dwlsmv.sthlp and spec 4.8 together."
    if !`hasn' | !`conv' | !`convp' | `lpmin' <= 0 {
        display as error " 5. settled facts broken: canon has fixed -.4: `hasn'; converged at -.4 / +.4: `conv' / `convp'; smallest loading at +.4: `lpmin' (want > 0)"
        local ++nfail
    }
    else if `lmax' < 0 & abs(`gap') <= 10 {
        display as error " 5. the limitation looks LIFTED: all loadings at -.4 are negative (largest `lmax') with T3 within 10 of the +.4 fit (gap `gap'). Update this check and the documentation together."
        local ++nfail
    }
    else if `lmin' <= 0 | `gap' <= 10 {
        display as error " 5. behaviour at -.4 changed: smallest loading `lmin' (characterized: all positive), T3 gap `gap' (characterized: above 10)"
        local ++nfail
    }
    else display as text " 5. CHARACTERIZED: gamma fixed at -.4 under unitvariance converges (local optimum): loadings all positive, T3 above the +.4 fit by more than 10"
}

* ---- 5b. the way out under the marker metric, derived and not measured ----------------------------------------
* Under the marker metric the first indicator's loading is FIXED at 1 (a fixed nonzero loading), so the factor is
* non-reflectable (dwlsmv_signcanon(): a fixed nonzero loading on the factor sets can = 0) and its orientation is
* the marker's. (F <- x@c) then fixes the effect of x on the marker's latent response, in the marker's units; the
* model at -.4 is a different model from the one at +.4, not a mirror image, so no T3 relation is asserted. Only
* what follows from the table: the fit converges and the marker loading is exactly 1 at both signs. (Fixing +.4
* under unitvariance and reading the sign, the other way out, is the +.4 fit of check 5: its loadings are
* positive; reversing the factor is the interpretation, with no further assertion derivable.)
local ok5b 1
foreach c in -.4 .4 {
    capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@`c')
    if _rc {
        display as error " 5b. marker metric, gamma fixed at `c': rc " _rc
        local ok5b 0
        continue
    }
    local m5 = .
    capture local m5 = _b[y1:F]
    if _rc | e(converged) != 1 | `m5' != 1 {
        display as error " 5b. marker metric, gamma fixed at `c': converged " e(converged) ", marker loading _b[y1:F] " cond(missing(`m5'), ".", string(`m5'))
        local ok5b 0
    }
}
if `ok5b' display as text " 5b. marker metric: gamma fixed at -.4 and at .4 both converge with the marker loading exactly 1 (orientation set by the marker)"
else local ++nfail

* ---- 6. identification ----------------------------------------------------------------------------------------
local mAll "(F -> y1 y2 y3 y4) (F <- x) (y1 <- x) (y2 <- x) (y3 <- x) (y4 <- x)"
tf_run "free gamma" "kappa" dwlsmv `mAll'
if r(rc) != 3498 {
    display as error " 6a. free gamma, every kappa free: rc " r(rc) " (want 3498)"
    local ++nfail
}
else display as text " 6a. free gamma with every kappa free: rc 3498"

capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x@.4) (y1 <- x) (y2 <- x) (y3 <- x) (y4 <- x)
local rc6b = _rc
if `rc6b' {
    display as error " 6b. gamma fixed at .4, every kappa free: rc `rc6b' (want 0: no free gamma, no anchor needed; the rank check passes since fixing removes a column)"
    local ++nfail
}
else display as text " 6b. the same model with gamma fixed runs (df " e(df_T3) ")"

capture noisily dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y1 <- x@.1) (y2 <- x@.1) (y3 <- x@.1) (y4 <- x@.1)
local rc6c = _rc
if `rc6c' {
    display as error " 6c. gamma free, every kappa fixed at .1: rc `rc6c' (want 0: a fixed kappa is not free)"
    local ++nfail
}
else {
    local canon `"`e(canon)'"'
    mata: st_local("hasq", strofreal(strpos(st_local("canon"), "kappa y1 x fixed .1") > 0))
    if !`hasq' {
        display as error " 6c. canon lacks kappa y1 x fixed .1: `canon'"
        local ++nfail
    }
    else display as text " 6c. gamma free, every kappa fixed at .1: runs"
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if `nfail' {
    display as error _n "test_condx_fixed: `nfail' failure(s)"
    exit 9
}
display as result "test_condx_fixed: a value fixed with @ on a covariate path is exact in e(b), e(Cns) and e(canon), adds one degree of freedom, and (@0) is the same model as no path; a negative fixed gamma under unitvariance is characterized as a known limitation; the anchor screen counts only free gamma and free kappa"
exit 0

* ======================================================================================================
* MUTANTS (apply one at a time, run this file, expect exit 9; restore after)
*  M1 src/mata/dwlsmv_translate.mata line 134
*       old: fixval[cj] = C[r, cols(C)] / C[r, cj]
*       new: fixval[cj] = abs(C[r, cols(C)] / C[r, cj])
*       killed by check 2 (_b[y4:x] == -.3, canon `kappa y4 x fixed -.3').
*  M2 src/mata/dwlsmv_translate.mata line 201
*       old (lines 200-201, the gamma line, not the kappa line at 174-175): the term (fixed[j] ? " fixed " + dwlsmv_tr_num(fixval[j]) : " free")
*       new: " free")
*       killed by check 1 (canon has `gamma F x fixed .4', _b == .4, df + 1).
*  M3 src/mata/dwlsmv_m11_fit.mata line 148
*       old: bvec[i] = (par.free[i] ? fit.theta[par.fidx[i]] : par.value[i])
*       new: bvec[i] = (par.free[i] ? fit.theta[par.fidx[i]] : 0)
*       killed by checks 1 and 2 (_b of a fixed row is c, not 0).
*  M4 src/mata/dwlsmv_m11_fit.mata line 190
*       old: st_local("rank",       strofreal(par.npar))
*       new: st_local("rank",       strofreal(par.nall))
*       killed by checks 1 and 2 (e(rank) would read 13 = nall, not 12).
*  M5 src/mata/dwlsmv_signcanon.mata line 100
*       old: if (!par.free[i] & par.value[i] != 0) {
*       new: if (!par.free[i] & par.value[i] != 0 & t != 8) {
*       NO LONGER killed by check 5 after the 2026-11-14 rewrite: the loadings at -.4 are already positive, so
*       there is nothing for a reflection to flip. Expected to survive this file; the non-reflectability is
*       covered by the signcanon unit test, if it has a fixed-gamma case (not checked by the drafting agent).
*  M6 src/mata/dwlsmv_id_cx.mata line 301
*       old: if (gmat[k, 1] == a & gmat[k, 3] == 1) {
*       new: if (gmat[k, 1] == a) {
*       killed by check 6b (a fixed gamma now triggers the anchor screen: rc 3498 instead of a fit).
*  Bonus: dwlsmv_id_cx.mata line 317 `kmat[k, 3] == 1' -> `kmat[k, 3] == 1 | kmat[k, 3] == 0' is NOT killed
*  by anything here (kappa fixed counts as non-free either way); `kmat[k, 3] >= 0' is killed by 6c.
