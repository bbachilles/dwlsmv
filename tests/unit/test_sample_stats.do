*! Stage 1's v1 pair policy (spec Module 2, *Which pairs enter s and Gamma*,
*! and *Reporting*; spec 5.1's keep_bound and stage1.bound; decided
*! 2026-10-03).
*!
*! Hand-built tables, each probed once with dwls_polychoric() before it was
*! chosen (2026-10-03):
*!   REM   (60,1,0 \ 1,60,1 \ 0,1,60): status 2 at 0.999, K = 3, not an empty-cell pair.
*!   KE1   L9d, (30,2 \ 0,8): one-empty-cell pair.
*!   KE2   (20,0 \ 0,20): perfect-association pair (an empty diagonal pair).
*!   OK    (30,5 \ 6,25): status 1.
*! A third variable, y3 = mod(i, K3), makes each a three-variable Stage 1.
*!
*!   1. Default (keep_bound absent or 0): REM and KE2 abort at rc 3498; OK
*!      passes with bound all zero and of length P. KE1 passes since the
*!      one-empty-cell repair (O6, 2026-10-06), with bound all zero; its values
*!      are tested in test_emptycell.do.
*!   2. keep_bound = 1: REM passes, bound marks exactly pair 1 (y1~~y2), and
*!      s carries 0.999 there; KE2 still aborts (the perfect-association pair
*!      is never kept); KE1 passes, repaired, with bound all zero (keep_bound
*!      does not touch the repair).
*!   3. The messages, read back from a log: every offending pair is named
*!      with its status in words; the perfect-association pair is named from
*!      its table (KE2), and KE1's repair is listed with its table; the
*!      readings appear once and end in "not a diagnosis"; the opt-in is
*!      named for REM and NOT for empty-cell pairs; under keep_bound the included pair
*!      is listed and the inclusion caveat is printed.
*!   4. Contract: keep_bound = 2 and keep_bound = . abort at rc 3498.
*!   5. Zero-probe status 6 is status 6 at exactly three evaluations, and
*!      status 6 by the endpoint check has five or more (the rule Stage 1
*!      classifies by). Forced with fake evaluators, as test_rootfind does.
*!   6. the empty-cell pair's assertable check (@sec-kindE, Corollary K4) in production
*!      arithmetic: (+-0.999, status 2, 11 evaluations) on L9d, its row
*!      mirror and a diagonal pair.
*!
*! Re-executed under StataMP 18 on 2026-10-06 for the one-empty-cell repair
*! (O6) and it passes; mutants M5 and M9 of test_emptycell.do's header fail
*! it. The failure printout was fixed the same day: it passed the whole FAIL
*! vector to printf(), which aborts with r(3204) instead of listing failures.
*! Executed under StataMP 18 on 2026-10-03 and it passes. Two mutants were
*! run against it and each fails it: keep_bound admitting empty-cell pairs (checks 1-2's
*! aborts) and the readings paragraph dropped (check 3).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

foreach f in src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_catdegen.mata          ///
             src/mata/dwlsmv_cumprop.mata           ///
             src/mata/dwlsmv_thresholds.mata        ///
             src/mata/dwlsmv_psi_tau_table.mata     ///
             src/mata/dwlsmv_cellprob.mata          ///
             src/mata/dwlsmv_dcellprob_rho.mata     ///
             src/mata/dwlsmv_wmat.mata              ///
             src/mata/dwlsmv_crosstab.mata          ///
             src/mata/dwlsmv_pairinfo.mata          ///
             src/mata/dwls_polychoric.mata          ///
             src/mata/dwlsmv_koff.mata              ///
             src/mata/dwlsmv_psi_tau_stack.mata     ///
             src/mata/dwlsmv_psi_rho_table.mata     ///
             src/mata/dwlsmv_psi_matrix.mata        ///
             src/mata/dwlsmv_psibar.mata            ///
             src/mata/dwls_gamma.mata               ///
             src/mata/dwlsmv_emptycell1.mata      ///
             src/mata/dwlsmv_emptycell_repair.mata ///
             src/mata/dwlsmv_emptycell_pairs.mata ///
             src/mata/dwls_sample_stats.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_sample_stats: `f' failed to load"
        exit 9
    }
}

mata:

// A two-variable table expanded to rows, plus y3 = mod(i, k3).
real matrix dwlsmv_t_ss_data(real matrix n, real scalar k3)
{
    real matrix Y
    real scalar r, c

    Y = J(0, 2, .)
    for (r = 1; r <= rows(n); r++) {
        for (c = 1; c <= cols(n); c++) {
            if (n[r, c]) {
                Y = Y \ J(n[r, c], 1, (r - 1, c - 1))
            }
        }
    }
    return((Y, mod((1::rows(Y)), k3)))
}

void dwlsmv_t_ss_put(string scalar name, real matrix n, real scalar k3)
{
    real matrix Y
    Y = dwlsmv_t_ss_data(n, k3)
    st_matrix(name, Y)
    st_matrix(name + "_K", (rows(n) \ cols(n) \ k3))
}

real rowvector dwlsmv_t_ss_z1(real scalar x, real matrix n, real colvector a, real colvector b)
{
    real rowvector e
    e = dwlsmv_rf_eval(x, n, a, b)
    if (e[2] & abs(x) < 0.01) {
        e[3] = 0
    }
    return(e)
}

real rowvector dwlsmv_t_ss_z2(real scalar x, real matrix n, real colvector a, real colvector b)
{
    real rowvector e
    e = dwlsmv_rf_eval(x, n, a, b)
    if (e[2]) {
        if (x == 0 | abs(x) >= 0.999) {
            e[3] = 0
        }
        else {
            e[3] = -sign(x) * max((abs(x) * e[4], 2 * e[5]))
        }
    }
    return(e)
}

string colvector dwlsmv_t_ss_ke(real matrix n, real scalar want, string scalar tag)
{
    struct dwlsmv_pairfit scalar f
    real matrix Y

    Y = dwlsmv_t_ss_data(n, 2)
    f = dwls_polychoric(Y[., 1], Y[., 2], dwlsmv_thresholds(Y[., 1], 2),
                        dwlsmv_thresholds(Y[., 2], 2))
    if (f.status != 2 | f.rho != want | f.iterations != 11) {
        return("6: " + tag + " gave (" + strofreal(f.rho) + ", " + strofreal(f.status) +
               ", " + strofreal(f.iterations) + "), not (" + strofreal(want) + ", 2, 11)")
    }
    return(J(0, 1, ""))
}

real scalar dwlsmv_t_ss_main()
{
    struct dwlsmv_stage1 scalar st
    string colvector FAIL
    real matrix Y, n
    real colvector a, b
    real rowvector r
    real scalar i

    FAIL = J(0, 1, "")

    // ---- 1 and 2: what enters, what is marked --------------------------
    Y  = dwlsmv_t_ss_data((30, 5 \ 6, 25), 2)
    st = dwls_sample_stats(Y, (2 \ 2 \ 2))
    if (rows(st.bound) != st.map.P | cols(st.bound) != 1) {
        FAIL = FAIL \ "1: bound is not P x 1"
    }
    if (any(st.bound :!= 0)) {
        FAIL = FAIL \ "1: bound is not all zero on a healthy Stage 1"
    }
    st = dwls_sample_stats(Y, (2 \ 2 \ 2), 0)
    if (any(st.bound :!= 0)) {
        FAIL = FAIL \ "1: bound is not all zero with keep_bound = 0"
    }

    Y  = dwlsmv_t_ss_data((60, 1, 0 \ 1, 60, 1 \ 0, 1, 60), 3)
    st = dwls_sample_stats(Y, (3 \ 3 \ 3), 1)
    if (st.bound != (1 \ 0 \ 0)) {
        FAIL = FAIL \ "2: keep_bound does not mark exactly y1~~y2"
    }
    if (st.s[st.map.T + 1] != 0.999) {
        FAIL = FAIL \ "2: the kept pair does not enter s at 0.999"
    }
    if (rows(st.Gamma) != st.map.D) {
        FAIL = FAIL \ "2: no Gamma under keep_bound"
    }

    // KE1, repaired: enters by default and under keep_bound, never in bound.
    Y  = dwlsmv_t_ss_data((30, 2 \ 0, 8), 2)
    st = dwls_sample_stats(Y, (2 \ 2 \ 2))
    if (st.bound != J(3, 1, 0)) {
        FAIL = FAIL \ "1: the repaired KE1 is marked in bound"
    }
    st = dwls_sample_stats(Y, (2 \ 2 \ 2), 1)
    if (st.bound != J(3, 1, 0)) {
        FAIL = FAIL \ "2: the repaired KE1 is marked in bound under keep_bound"
    }

    // ---- 5: zero-probe status 6 is (0, 6, 3) ---------------------------
    n = (30, 5 \ 6, 25)
    a = invnormal(35 / 66)
    b = invnormal(36 / 66)
    r = dwlsmv_rootfind(&dwlsmv_t_ss_z1(), n, a, b, (., .), 0.999, 40, 54)
    if (r != (0, 6, 3)) {
        FAIL = FAIL \ ("5: zero-probe ambiguity gave (" + invtokens(strofreal(r)) + "), not (0 6 3)")
    }
    r = dwlsmv_rootfind(&dwlsmv_t_ss_z2(), n, a, b, (., .), 0.999, 40, 54)
    if (r[1] != 0 | r[2] != 6 | r[3] < 5) {
        FAIL = FAIL \ ("5: endpoint-check status 6 gave (" + invtokens(strofreal(r)) + "), not (0 6 >=5)")
    }

    // ---- 6: the empty-cell pair's assertable check (@sec-kindE, Corollary K4) --------
    // Production's binormal() on L9d, its row mirror and a diagonal pair:
    // status 2, rho-hat = -sigma_e B exactly, 11 evaluations.
    FAIL = FAIL \ dwlsmv_t_ss_ke((30, 2 \ 0, 8), 0.999, "L9d")
    FAIL = FAIL \ dwlsmv_t_ss_ke((0, 8 \ 30, 2), -0.999, "L9d row mirror")
    FAIL = FAIL \ dwlsmv_t_ss_ke((20, 0 \ 0, 20), 0.999, "diagonal pair")

    // one line per failure: printf() takes a scalar, and a vector here
    // aborts with r(3204) instead of listing the failures
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}%s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// Called for its message only. Assigning the result keeps Mata from
// printing the returned struct's address, which differs run to run.
void dwlsmv_t_ss_keep()
{
    struct dwlsmv_stage1 scalar st
    st = dwls_sample_stats(st_matrix("REM"), st_matrix("REM_K"), 1)
}

dwlsmv_t_ss_put("REM", (60, 1, 0 \ 1, 60, 1 \ 0, 1, 60), 3)
dwlsmv_t_ss_put("KE1", (30, 2 \ 0, 8), 2)
dwlsmv_t_ss_put("KE2", (20, 0 \ 0, 20), 2)
dwlsmv_t_ss_put("OK",  (30, 5 \ 6, 25), 2)

end

mata: st_local("nfail", strofreal(dwlsmv_t_ss_main()))

local nerr 0

* ---- 1, 2, 4: the aborts, each at rc 3498 --------------------------------
foreach c in REM KE2 {
    capture mata: dwls_sample_stats(st_matrix("`c'"), st_matrix("`c'_K"))
    if _rc != 3498 {
        display as error "test_sample_stats: `c' by default gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}
foreach c in KE2 {
    capture mata: dwls_sample_stats(st_matrix("`c'"), st_matrix("`c'_K"), 1)
    if _rc != 3498 {
        display as error "test_sample_stats: `c' under keep_bound gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}
foreach v in 2 . {
    capture mata: dwls_sample_stats(st_matrix("OK"), st_matrix("OK_K"), `v')
    if _rc != 3498 {
        display as error "test_sample_stats: keep_bound = `v' gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}

* ---- 3: the messages, read back from a log -------------------------------
tempfile L
foreach c in REM KE1 KE2 KEEP {
    quietly log using "`L'_`c'.log", text replace name(ss)
    if "`c'" == "KEEP" {
        capture noisily mata: dwlsmv_t_ss_keep()
    }
    else {
        capture noisily mata: dwls_sample_stats(st_matrix("`c'"), st_matrix("`c'_K"))
    }
    quietly log close ss
}

mata:
real scalar dwlsmv_t_ss_msg(string scalar f, string rowvector has, string rowvector hasnt,
                            string scalar once, string scalar tag)
{
    string scalar txt
    real scalar i, bad, p

    txt = invtokens(cat(f)', char(10))
    bad = 0
    for (i = 1; i <= cols(has); i++) {
        if (!strpos(txt, has[i])) {
            printf("{err}test_sample_stats: %s message lacks: %s\n", tag, has[i])
            bad = bad + 1
        }
    }
    for (i = 1; i <= cols(hasnt); i++) {
        if (strpos(txt, hasnt[i])) {
            printf("{err}test_sample_stats: %s message should not contain: %s\n", tag, hasnt[i])
            bad = bad + 1
        }
    }
    if (once != "") {
        p = strpos(txt, once)
        if (!p | strpos(substr(txt, p + 1, .), once)) {
            printf("{err}test_sample_stats: %s message does not state once: %s\n", tag, once)
            bad = bad + 1
        }
    }
    if (strpos(txt, "\" + "n")) {
        printf("{err}test_sample_stats: %s message has a literal backslash-n\n", tag)
        bad = bad + 1
    }
    return(bad)
}
end

local tmpL "`L'"
mata: st_local("nmsg", strofreal(                                                     ///
    dwlsmv_t_ss_msg(st_local("tmpL") + "_REM.log", ("y1~~y2", "status 2 (no sign change up to the bound)", ///
        "cannot enter the sample statistics", "near-duplicate (redundant) items", "not a diagnosis",       ///
        "keep_bound = 1"), ("empty cell", "remove", "detected", "invalid"), "not a diagnosis", "REM") +    ///
    dwlsmv_t_ss_msg(st_local("tmpL") + "_KE2.log", ("y1~~y2", "2x2 table with two empty cells: no interior maximum, not repaired", ///
        "not a diagnosis"), ("keep_bound = 1", "repair applied"), "not a diagnosis", "KE2") +          ///
    dwlsmv_t_ss_msg(st_local("tmpL") + "_KE1.log", ("one-empty-cell repair applied to 1 of 3 pairs",        ///
        "y1~~y2: rho-hat 0.9656, table (30, 2 \ 0, 8)", "margins and thresholds are unchanged",           ///
        "is not established"), ("cannot enter", "keep_bound = 1"), "is not established", "KE1") +      ///
    dwlsmv_t_ss_msg(st_local("tmpL") + "_KEEP.log", ("included at the bound (keep_bound)", "y1~~y2",       ///
        "as lavaan does", "is not established", "not a diagnosis"), ("cannot enter"), "not a diagnosis", "KEEP")))

if `nfail' | `nerr' | `nmsg' {
    display as error _n "test_sample_stats: `nfail' value failure(s), `nerr' error-path failure(s), `nmsg' message failure(s)"
    exit 9
}

display as result "test_sample_stats: v1 pair policy verified -- aborts, repair, keep_bound, bound, messages, zero-probe status 6"
exit 0
