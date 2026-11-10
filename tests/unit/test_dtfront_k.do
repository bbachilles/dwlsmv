*! difftest at G > 1 asserts e(K) against the refit's own Stage 1 (spec 4.8 G8,
*! implemented 2026-11-02, user). A review found the sentence unimplemented. Both
*! models are given the same forged e(K), so the exemption's own comparison of
*! the two cannot catch it. Control: the same fits with e(K) untouched.
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

program tfk_forge, eclass
    tempname K
    matrix `K' = e(K)
    matrix `K'[1, 1] = `K'[1, 1] + 1
    ereturn matrix K = `K'
end

program tfk_log
    args cmd
    tempfile lg
    capture log close tfk
    quietly log using "`lg'", text replace name(tfk)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tfk
    mata: TFKLOG = invtokens(cat(st_local("lg"))')
    c_local tfk_rc `rc'
end

set seed 20261102
set obs 2000
gen g = 1 + (_n > 1000)
gen F = rnormal()
forvalues j = 1/6 {
    gen s`j' = .65*F + sqrt(1-.65^2)*rnormal()
    egen y`j' = cut(s`j'), at(-100 -.6 .2 .9 100) icodes
}
local m (F -> y1 y2 y3 y4 y5 y6)
local nfail 0
foreach L in configural thresholds metric {
    quietly dwlsmv `m', group(g) invariance(`L') metric(unitvariance)
    estimates store a_`L'
    tfk_forge
    estimates store f_`L'
}

foreach pair in "a_thresholds a_configural" "a_metric a_thresholds" {
    tfk_log "difftest `pair'"
    if `tfk_rc' {
        display as error "  control `pair' FAILED: rc `tfk_rc'"
        local ++nfail
    }
    else display as text "  control `pair': difftest runs"
}

foreach pair in "f_thresholds f_configural" "f_metric f_thresholds" {
    tfk_log "difftest `pair'"
    mata: st_local("hit", strofreal(strpos(TFKLOG, "own Stage 1") > 0))
    if `tfk_rc' != 3498 | !`hit' {
        display as error "  forged `pair' FAILED (rc `tfk_rc', message `hit')"
        local ++nfail
    }
    else display as text "  forged `pair': refused, naming e(K)"
}

if `nfail' {
    display as error "test_dtfront_k: `nfail' failure(s)"
    exit 9
}
display as result "test_dtfront_k: all checks passed"
exit 0
