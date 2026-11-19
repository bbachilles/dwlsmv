*! dwlsmv -- A2b follow-up: does the null-direction behaviour of test/lincom
*! depend on the Stata version or on the sign of the rounding residue?
*!
*! Executed under StataMP 18 and StataMP 19 (6 processors) on 2026-11-17, main
*! session; record in tools/probe_a2b_sign_stata18.txt (Stata 19 identical apart from
*! the SGENV line). Written after tools/probe_a2b_singular_v.do gave different
*! null-direction results on Stata 16.1 (milestone machine, ARM) and Stata 18 (x86):
*! that probe rebuilds V on each machine, so the two runs did not test the same bits.
*! Here probe 1's V is rebuilt and given a rank-one correction along the null
*! direction k1, so that k1'V k1 is a small positive (Vp) or negative (Vm) residue,
*! then test and lincom of k1'b are run on each. Not a test; chooses no convention.
*! Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_a2b_sign.do      (log: probe_a2b_sign.log)
*!     grep -E "^SG|dropped" probe_a2b_sign.log | diff - tools/probe_a2b_sign_stata18.txt
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 200
display "SGENV stata_version " c(stata_version)
set obs 100
mata:
void sg_build()
{
    real matrix K, X, Q, R, Nn, V, Gm, Vp, Vm
    real colvector k
    real scalar q, r
    rseed(20261117)
    K  = rnormal(25, 2, 0, 1)
    X  = (K, rnormal(25, 23, 0, 1))
    qrd(X, Q, R)
    Gm = rnormal(23, 23, 0, 1)
    Nn = Q[., 3..25] * Gm * 0.02
    V  = Nn * Nn'
    V  = (V + V') / 2
    k  = K[., 1]
    q  = k' * V * k
    printf("SG k1'Vk1 as built %12.5e\n", q)
    r  = 2.5e-17
    // set the null-direction variance to +r and to -r exactly (k k' rank-one correction)
    Vp = V + ((r - q) / (k' * k)^2) * (k * k')
    Vm = V + ((-r - q) / (k' * k)^2) * (k * k')
    printf("SG k1'Vp k1 %12.5e   k1'Vm k1 %12.5e\n", k' * Vp * k, k' * Vm * k)
    st_matrix("SG_Vp", Vp)
    st_matrix("SG_Vm", Vm)
    st_matrix("SG_b", J(1, 25, 0.3))
    st_matrix("SG_k", k')
}
end
mata: sg_build()
local nm
forvalues i = 1/25 {
    local nm `nm' c`i'
}
local e
forvalues i = 1/25 {
    local kk = strofreal(SG_k[1, `i'], "%21.15g")
    local e `e' + (`kk')*_b[c`i']
}
local e = substr("`e'", 3, .)
program sgpost, eclass
    args vn nm
    tempname b V
    matrix `b' = SG_b
    matrix `V' = `vn'
    matrix colnames `b' = `nm'
    matrix rownames `V' = `nm'
    matrix colnames `V' = `nm'
    ereturn post `b' `V'
end
foreach s in p m {
    sgpost SG_V`s' "`nm'"
    capture noisily test `e' = 0
    display "SG `s' test rc " _rc " chi2 " r(chi2) " df " r(df) " drop " r(drop)
    quietly count
    capture noisily test `e' = .1
    display "SG `s' test_shift rc " _rc " chi2 " r(chi2) " df " r(df) " drop " r(drop)
    quietly count
    capture noisily lincom `e'
    display "SG `s' lincom rc " _rc " est " r(estimate) " se " r(se)
}
