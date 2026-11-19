*! dwlsmv -- A2b probe 1: a posted e(V) that is singular with no zero rows
*!
*! Executed under StataMP 18 (6 processors) on 2026-11-17; record in
*! tools/probe_a2b_singular_v_stata18.txt. Drafted by a read-only agent; on the first
*! runs the main session fixed two probe bugs: the helper read the touse variable that
*! ereturn post, esample() consumes (now quietly count), and test mixed [eq]coef with a
*! parenthesised constraint, read as a weight (now all parenthesised).
*!
*! A probe, not a test (tools/); the suite runner ignores it. Roadmap unit A2b
*! (E-108): what ereturn post, ereturn display, test, lincom and
*! e(rank) do with an e(V) that is singular by construction, with NO zero row, at the
*! size and rank of Case E's rotated EFA covariance. It reports what happens and
*! chooses no convention. Style and background: tools/probe_d5_post.do (ereturn post
*! does not set e(rank); a non-PD V posts with a warning and SEs blanked).
*!
*! Dimensions. Case E: p = 6 items of 3 categories (2 thresholds each, 12), m = 2
*! factors (derivation @sec-rot-se, rs-setup, and tests/unit/test_efa_unrotated.do).
*! n* = 12 loadings + m(m-1)/2 = 1 correlation + 12 thresholds = 25; n_u = 25 - m(m-1)
*! = 23, so V has rank 23 and a null space of dimension 2. The stripe follows
*! dwlsmv_m11_efa_rotate()'s e(b): loadings y1:F1 .. y6:F1, y1:F2 .. y6:F2, then
*! /:cov(F1,F2), then /y1:cut1 /y1:cut2 .. /y6:cut2.
*!
*! V = N N' (scaled), N = Q[,3..25] * G, Q the full orthogonal factor of a QR of
*! [K, random], K (25 x 2) the null basis, G a random 23 x 23: so K'V = 0 to rounding,
*! rank 23 (asserted), every row nonzero (asserted). seed 20261117. Control (d) is
*! V + .001 I (a full-rank V of the same scale; the brief's V + I would make every SE
*! about 1).
*!
*! Variants, each posted by its own eclass program:
*!   a  V as is, e(rank) as ereturn post leaves it
*!   b  V, then ereturn scalar rank = 23
*!   c  V, then ereturn scalar rank = 25
*!   d  V + .001 I (full rank), no e(rank) set
*! For each: post rc and e(rank); ereturn display (rc; the count of rows with missing
*! SE, z, p, CI from r(table) if the version stores it); _b[] and _se[] of a few
*! coefficients; test (single coefficient; all F1 loadings; one null-direction
*! constraint alone, exact and shifted; both null directions; a joint test mixing an
*! ordinary constraint with a null-direction one; the same with both null directions);
*! lincom (ordinary; along a null direction; shifted); estimates store / restore.
*! A warning from ereturn post, "constraint dropped" notes of test, and the table
*! itself are not on the A2BV lines: grep the FULL log, not only the A2BV lines.
*!
*! Run from the repository root (PowerShell):
*!     Start-Process -Wait -FilePath "C:\Program Files\Stata18\StataMP-64.exe" -ArgumentList "/e","do","tools/probe_a2b_singular_v.do"
*! (Git Bash: StataMP-64.exe //e do tools/probe_a2b_singular_v.do), then
*!     grep -E "^(====|A2BV| +Constraint)" probe_a2b_singular_v.log | diff - tools/probe_a2b_singular_v_stata18.txt
*! Only the A2BVENV line is expected to differ.
*! The roadmap wants StataMP 18; the same file is to be run on the Stata 16 milestone
*! machine if practical.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 200
display "A2BVENV stata_version " c(stata_version) " processors " c(processors)
set obs 100
generate byte touse = 1

// ---------------------------------------------------------------- the data
mata:
void a2bv_build()
{
    real matrix K, X, Q, R, Nn, V, Gm, ev, Xe
    real rowvector b
    real colvector sv, ssv
    real scalar n, rk, mm, c1, c2, j, i, nsv

    n  = 25
    mm = 2
    rk = 23
    rseed(20261117)

    K  = rnormal(n, mm, 0, 1)
    X  = (K, rnormal(n, n - mm, 0, 1))
    qrd(X, Q, R)
    Gm = rnormal(rk, rk, 0, 1)
    Nn = Q[., (mm + 1)..n] * Gm * 0.02
    V  = Nn * Nn'
    V  = (V + V') / 2

    // b: 12 loadings (F1 then F2), the correlation, 12 thresholds
    b = J(1, n, .)
    for (i = 1; i <= 12; i++) {
        b[i] = 0.2 + 0.6 * runiform(1, 1)
    }
    b[13] = 0.3
    for (j = 1; j <= 6; j++) {
        b[13 + 2 * j - 1] = -0.5 + 0.05 * j
        b[13 + 2 * j]     = 0.5 + 0.05 * j
    }

    // the facts about V, asserted
    if (min(rowsum(abs(V))) <= 0) {
        _error("A2BV: V has a zero row")
    }
    symeigensystem(V, Xe, ev)
    sv  = svdsv(V)
    nsv = sum(sv :> 1e-8 * max(sv))
    ssv = sort(sv, 1)
    printf("A2BV build n %g seeded 20261117; every row of V nonzero (min row abs sum %9.3e)\n",
           n, min(rowsum(abs(V))))
    printf("A2BV build eig min %12.5e max %12.5e; rank() %g; singular values above 1e-8 relative %g (expect 23)\n",
           min(ev), max(ev), rank(V), nsv)
    printf("A2BV build smallest 4 singular values %12.5e %12.5e %12.5e %12.5e\n",
           ssv[1], ssv[2], ssv[3], ssv[4])
    printf("A2BV build max|K'V| %9.3e (null basis K, 25 x 2); max|V-V'| %g\n",
           max(abs(K' * V)), max(abs(V - V')))
    if (nsv != rk) {
        _error("A2BV: the rank of V is not 23")
    }

    c1 = b * K[., 1]
    c2 = b * K[., 2]
    st_matrix("A2BV_b", b)
    st_matrix("A2BV_V", V)
    st_matrix("A2BV_K", K)
    st_numscalar("A2BV_c1", c1)
    st_numscalar("A2BV_c2", c2)
}
end

mata: a2bv_build()

// stripe names (ereturn post) and the _b[] spelling of each coefficient
global A2BV_NM
global A2BV_REF
foreach f in 1 2 {
    forvalues j = 1/6 {
        global A2BV_NM "$A2BV_NM y`j':F`f'"
        global A2BV_REF "$A2BV_REF y`j':F`f'"
    }
}
global A2BV_NM "$A2BV_NM /:cov(F1,F2)"
global A2BV_REF "$A2BV_REF /cov(F1,F2)"
forvalues j = 1/6 {
    forvalues c = 1/2 {
        global A2BV_NM "$A2BV_NM /y`j':cut`c'"
        global A2BV_REF "$A2BV_REF /y`j':cut`c'"
    }
}
matrix colnames A2BV_b = $A2BV_NM
matrix rownames A2BV_V = $A2BV_NM
matrix colnames A2BV_V = $A2BV_NM
display "A2BVNAMES " wordcount("$A2BV_NM") " coefficients: $A2BV_NM"

// the two null-direction linear forms K_k'b as _b[] expressions, and their constants
forvalues k = 1/2 {
    global A2BV_N`k'
    forvalues i = 1/25 {
        local ref : word `i' of $A2BV_REF
        local kk = strofreal(A2BV_K[`i', `k'], "%21.15g")
        if `i' == 1 {
            global A2BV_N`k' "(`kk')*_b[`ref']"
        }
        else {
            global A2BV_N`k' "${A2BV_N`k'} + (`kk')*_b[`ref']"
        }
    }
    global A2BV_C`k' = strofreal(A2BV_c`k', "%21.15g")
}
global A2BV_F1
global A2BV_F1P
forvalues j = 1/6 {
    global A2BV_F1 "$A2BV_F1 [y`j']F1"
    global A2BV_F1P "$A2BV_F1P (_b[y`j':F1] = 0)"
}

// ---------------------------------------------------------------- post
program a2bvpost, eclass
    args v
    tempname b V
    capture drop touse
    generate byte touse = 1
    matrix `b' = A2BV_b
    matrix `V' = A2BV_V
    if "`v'" == "d" {
        matrix `V' = A2BV_V + .001 * I(25)
        matrix rownames `V' = $A2BV_NM
        matrix colnames `V' = $A2BV_NM
    }
    ereturn post `b' `V', obs(100) esample(touse) properties(b V)
    if "`v'" == "b" {
        ereturn scalar rank = 23
    }
    if "`v'" == "c" {
        ereturn scalar rank = 25
    }
    ereturn local cmd "a2bvprobe"
end

// run one command, print rc and every r() value that matters, one A2BV line
program a2bvrun
    args v tag cmd
    quietly count
    capture noisily `cmd'
    local rc = _rc
    display as result "A2BV `v' `tag' rc `rc' chi2 " r(chi2) " df " r(df) " p " r(p) ///
        " drop " r(drop) " F " r(F) " est " r(estimate) " se " r(se) " z " r(z)
end

// ---------------------------------------------------------------- one variant
foreach v in a b c d {
    display as text "==== variant `v'"
    capture noisily a2bvpost `v'
    display as result "A2BV `v' post rc " _rc
    if _rc {
        continue
    }
    display as result "A2BV `v' e_rank " e(rank) " e_N " e(N) " e_df_r " e(df_r)

    // coefficient table
    capture noisily ereturn display
    display as result "A2BV `v' ereturn_display rc " _rc
    tempname T
    capture matrix `T' = r(table)
    if _rc {
        display as result "A2BV `v' r_table absent (rc " _rc ")"
    }
    else {
        local nse 0
        local nz 0
        local np 0
        local nll 0
        local nul 0
        forvalues i = 1/25 {
            if missing(`T'[2, `i']) local ++nse
            if missing(`T'[3, `i']) local ++nz
            if missing(`T'[4, `i']) local ++np
            if missing(`T'[5, `i']) local ++nll
            if missing(`T'[6, `i']) local ++nul
        }
        display as result "A2BV `v' table cols " colsof(`T') " missing_se `nse' missing_z `nz' missing_p `np' missing_ll `nll' missing_ul `nul'"
        display as result "A2BV `v' table se_first3 " `T'[2,1] " " `T'[2,2] " " `T'[2,3] " se_cov " `T'[2,13] " se_cut1 " `T'[2,14]
    }

    // _b[] and _se[]
    foreach r in y1:F1 y4:F2 /cov(F1,F2) /y1:cut1 {
        capture noisily display as result "A2BV `v' _b[`r'] " _b[`r'] " _se[`r'] " _se[`r']
        if _rc {
            display as result "A2BV `v' _b[`r'] rc " _rc
        }
    }

    // test
    a2bvrun `v' test_single "test _b[y1:F1] = 0"
    a2bvrun `v' test_F1_loadings "test $A2BV_F1"
    a2bvrun `v' test_null1_exact "test $A2BV_N1 = $A2BV_C1"
    a2bvrun `v' test_null1_shift "test $A2BV_N1 = $A2BV_C1 + .1"
    a2bvrun `v' test_null1_null2 "test ($A2BV_N1 = $A2BV_C1) ($A2BV_N2 = $A2BV_C2)"
    a2bvrun `v' test_ord_plus_null1 "test (_b[y1:F1] = 0) ($A2BV_N1 = $A2BV_C1)"
    a2bvrun `v' test_ord_plus_null1_shift "test (_b[y1:F1] = 0) ($A2BV_N1 = $A2BV_C1 + .1)"
    a2bvrun `v' test_ord_plus_both "test (_b[y1:F1] = 0) ($A2BV_N1 = $A2BV_C1) ($A2BV_N2 = $A2BV_C2)"
    a2bvrun `v' test_F1_plus_null1 "test $A2BV_F1P ($A2BV_N1 = $A2BV_C1)"

    // lincom
    a2bvrun `v' lincom_ord "lincom _b[y1:F1] - _b[y2:F1]"
    a2bvrun `v' lincom_null1 "lincom $A2BV_N1 - $A2BV_C1"
    a2bvrun `v' lincom_null1_shift "lincom $A2BV_N1 - $A2BV_C1 + .1"
    a2bvrun `v' lincom_null1_alone "lincom $A2BV_N1"

    // estimates store / restore
    capture noisily estimates store a2bvm_`v'
    display as result "A2BV `v' estimates_store rc " _rc
    capture noisily estimates restore a2bvm_`v'
    display as result "A2BV `v' estimates_restore rc " _rc " e_rank " e(rank) " e_cmd " e(cmd)
    capture noisily ereturn display
    display as result "A2BV `v' ereturn_display_after_restore rc " _rc
}
display as text "A2BV done"
