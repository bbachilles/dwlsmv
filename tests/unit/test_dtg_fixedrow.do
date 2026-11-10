*! Check 8's generated-chain exemption and the I2 note on a fixed nonzero row
*! (2026-11-02, user). With a fixed nonzero factor covariance or loading the
*! levels are not nested by construction: metric against thresholds was
*! demonstrated on simulated data, T3 24.83 on 28 df against 1524.80 on 24. Such
*! pairs now fail the exemption and reach check 8's tie refusal, and the I2 note
*! does not claim equivalence. A valid pair (scalar in metric) is refused too,
*! until a derivation covers these rows. Controls: free and zero covariance.
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

program tfr_log
    args cmd
    tempfile lg
    capture log close tfr
    quietly log using "`lg'", text replace name(tfr)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tfr
    mata: TFRLOG = invtokens(cat(st_local("lg"))')
    c_local tfr_rc `rc'
end

set seed 20261102
set obs 4000
gen g = 1 + (_n > 2000)
gen e1 = rnormal()
gen e2 = rnormal()
gen F1 = e1
gen F2 = .4*e1 + sqrt(1-.16)*e2
forvalues j = 1/6 {
    local F = cond(`j'<=3, "F1", "F2")
    gen s`j' = .65*`F' + sqrt(1-.65^2)*rnormal()
    egen y`j' = cut(s`j'), at(-100 -.6 .2 .9 100) icodes
}
local nfail 0
local m (F1 -> y1 y2 y3) (F2 -> y4 y5 y6)

* the fixed nonzero covariance: every pair refused at check 8, with the note
foreach L in configural thresholds metric scalar {
    quietly dwlsmv `m', cov(F1*F2@.4) group(g) invariance(`L') metric(unitvariance)
    estimates store a_`L'
}
foreach pair in "a_metric a_thresholds" "a_scalar a_metric" "a_thresholds a_configural" {
    tfr_log "difftest `pair'"
    mata: st_local("hit", strofreal(strpos(TFRLOG, "fixed nonzero") > 0))
    if `tfr_rc' != 3498 | !`hit' {
        display as error "  `pair' FAILED (rc `tfr_rc', note `hit')"
        local ++nfail
    }
    else display as text "  `pair': refused, with the note"
}

* the I2 note on a design it would otherwise match
preserve
forvalues j = 2/6 {
    if `j' != 4 replace y`j' = (y`j' > 1)
}
tfr_log "dwlsmv `m', cov(F1*F2@.4) group(g) invariance(metric) metric(unitvariance)"
mata: st_local("hit", strofreal(strpos(TFRLOG, "equivalence is established") > 0))
if `tfr_rc' | `hit' {
    display as error "  I2 note FAILED (rc `tfr_rc', claims equivalence `hit')"
    local ++nfail
}
else display as text "  I2 note: no equivalence claim with a fixed covariance"
restore

* a fixed nonzero loading: the surface refuses it before fitting, so the
* exemption's guard is reached only from canonical text
capture quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y3@.3 y5 y6), group(g) invariance(metric) metric(unitvariance)
if _rc != 198 {
    display as error "  surface cross-loading FAILED: rc `=_rc'"
    local ++nfail
}
else display as text "  surface cross-loading: refused before fitting"
mata: st_local("g1", strofreal(dwlsmv_dtg_fixedrow("lambda F2 y3 fixed .3" + char(10))))
mata: st_local("g2", strofreal(dwlsmv_dtg_fixedrow("lambda F2 y3 fixed 0" + char(10) + "lambda F2 y4 free" + char(10) + "psi F2 F2 fixed 1" + char(10))))
mata: st_local("g3", strofreal(dwlsmv_dtg_fixedrow("psi F1 F2 fixed .4" + char(10))))
if `g1' != 1 | `g2' != 0 | `g3' != 1 {
    display as error "  guard FAILED: `g1' `g2' `g3', wanted 1 0 1"
    local ++nfail
}
else display as text "  guard: nonzero rows caught; zero rows and fixed variances are not"

* controls: free and zero covariance stay exempt
foreach c in "cov(F1*F2)" "cov(F1*F2@0)" {
    foreach L in metric scalar {
        quietly dwlsmv `m', `c' group(g) invariance(`L') metric(unitvariance)
        estimates store c_`L'
    }
    tfr_log "difftest c_scalar c_metric"
    if `tfr_rc' {
        display as error "  control `c' FAILED: rc `tfr_rc'"
        local ++nfail
    }
    else display as text "  control `c': difftest runs"
}

if `nfail' {
    display as error "test_dtg_fixedrow: `nfail' failure(s)"
    exit 9
}
display as result "test_dtg_fixedrow: all checks passed"
exit 0
