*! dwlsmv -- observed category counts for one variable (Module 1)
*!
*! Executed under StataMP 18 on 2026-09-19: tests/unit/test_thresholds.do
*! drives this on all four case_D columns and asserts every rejection path
*! below -- missing, non-integer, out of range above, and out of range below.
*!
*! The counting is trivial; the validation is the point. Every way the caller
*! can hand this routine something that is not an ordinal variable coded
*! 0..K-1 aborts here, because each of them produces a plausible number
*! further downstream rather than an error.
*!
*! Missing values in particular: Mata's `.` compares GREATER than every real,
*! so a missing code would sail past a lower-bound check and be counted into
*! whichever comparison it happened to satisfy. Listwise deletion belongs to
*! the caller (spec Section 3); arriving here with missings present means it
*! did not happen, which is worth an error rather than a silent recount.
*!
*! K is the DECLARED number of categories, not the number observed. That
*! distinction is the whole reason this routine returns counts rather than a
*! condensed table: an unobserved category still occupies a row, so
*! dwlsmv_catdegen() can see it and spec Section 4.5.6's two failure modes
*! stay distinguishable. A count vector that silently dropped empty categories
*! would erase exactly the information the classification needs.
*!
*! src/mata/dwlsmv_types.mata is NOT required -- this file returns a plain
*! vector and declares no struct.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_catcounts(real colvector y, real scalar K)
{
    real colvector n
    real scalar c

    if (cols(y) != 1) {
        _error("dwlsmv_catcounts: y must be a column vector")
    }
    if (rows(y) < 1) {
        _error("dwlsmv_catcounts: y has no observations")
    }
    if (K < 2 | K != floor(K)) {
        _error("dwlsmv_catcounts: K must be an integer >= 2")
    }
    // Each condition reports what was FOUND and what to do about it, not just
    // what was required. The detail goes through errprintf because _error()
    // truncates at 100 characters -- past that it raises rc 3300 and the
    // message is lost entirely. A project coding rule.
    if (sum(y :>= .) > 0) {
        errprintf("dwlsmv_catcounts: %s of %s observations are missing.\n",
                  strofreal(sum(y :>= .)), strofreal(rows(y)))
        errprintf("  Stage 1 refuses missing values rather than assuming" +
                  " deletion already happened.\n")
        errprintf("  Apply listwise deletion across ALL indicators first, so" +
                  " every pair is estimated\n")
        errprintf("  on the same sample -- deleting per variable would not be.\n")
        _error("dwlsmv_catcounts: missing values present -- apply listwise deletion first")
    }
    if (sum(y :!= floor(y)) > 0) {
        errprintf("dwlsmv_catcounts: %s observations are not whole numbers.\n",
                  strofreal(sum(y :!= floor(y))))
        errprintf("  v1 fits ordinal and binary indicators only. A continuous" +
                  " indicator must be\n")
        errprintf("  categorised before use, or left out of the model.\n")
        _error("dwlsmv_catcounts: non-integer codes -- v1 takes ordinal and binary only")
    }
    if (min(y) < 0) {
        errprintf("dwlsmv_catcounts: the lowest observed value is %s, but" +
                  " categories must be 0 .. K-1.\n", strofreal(min(y)))
        errprintf("  Recode so the lowest category is 0.\n")
        _error("dwlsmv_catcounts: negative category code -- recode from 0")
    }
    if (max(y) > K - 1) {
        errprintf("dwlsmv_catcounts: the highest observed value is %s, but" +
                  " K = %s declares 0 .. %s.\n",
                  strofreal(max(y)), strofreal(K), strofreal(K - 1))
        errprintf("  Either K is too small for this variable, or the data is" +
                  " coded 1 .. K and needs\n")
        errprintf("  shifting down by one. A 1..K coding is much the more" +
                  " common cause.\n")
        _error("dwlsmv_catcounts: category code exceeds K-1 -- check for 1..K coding")
    }

    n = J(K, 1, 0)
    for (c = 1; c <= K; c++) {
        n[c] = sum(y :== (c - 1))
    }

    return(n)
}

end
