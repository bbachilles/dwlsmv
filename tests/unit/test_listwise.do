*! Module 1: listwise deletion (spec Section 3).
*!
*! Executed under StataMP 18 on Windows, 2026-09-19, and it passes on the first
*! run. All three error paths fire.
*!
*! Every case here is hand-built, because the whole point of this routine is
*! that the right answer is obvious by inspection. A fixture would make the
*! expectation harder to check than the code.
*!
*! Five things are asserted, and the last two are the ones that matter:
*!
*!   1. A complete matrix comes back unchanged, and comes back IDENTICAL --
*!      the no-missing fast path returns Y itself, so this also checks that
*!      the shortcut did not quietly transpose, copy or reorder anything.
*!   2. Rows with a missing anywhere are dropped, and the survivors keep their
*!      original ORDER and CONTENT. A routine that dropped the right number of
*!      rows but permuted the rest would pass a count check and corrupt every
*!      downstream statistic.
*!   3. EXTENDED MISSINGS (.a through .z) count as missing. Stata has 27 of
*!      them and they are the ones a real dataset carries -- ".a" for "refused"
*!      is far more common in survey data than plain ".". A check written as
*!      `== .` would let all 26 extended codes through, and they would then be
*!      counted into whichever category comparison they happened to satisfy,
*!      because Mata orders every missing ABOVE every real number.
*!   4. Deletion composes with Module 1's own refusal to see a missing: the
*!      same matrix that dwlsmv_catcounts() rejects is accepted after passing
*!      through here. That is the contract between the two -- Stage 1 refuses
*!      missings rather than assuming deletion happened, so something must
*!      actually do the deleting.
*!   5. An all-missing matrix aborts rather than returning an empty one. An
*!      empty Y would otherwise surface as a division by N = 0 somewhere in
*!      Module 4, which is a long way from the cause.
*!
*! Module 1 exists, so this test never legitimately skips.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in src/mata/dwlsmv_listwise.mata      ///
             src/mata/dwlsmv_catcounts.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_listwise: `f' failed to load"
        exit 9
    }
}

mata:

real scalar dwlsmv_t_listwise()
{
    string colvector FAIL
    real matrix C, M, E, got
    real scalar i

    FAIL = J(0, 1, "")

    // ---- 1. a complete matrix survives intact ----------------------------
    C   = (0, 1 \ 1, 2 \ 2, 0)
    got = dwlsmv_listwise(C)
    if (rows(got) != 3 | cols(got) != 2) {
        FAIL = FAIL \ ("complete matrix came back " + strofreal(rows(got)) +
                       " x " + strofreal(cols(got)) + ", expected 3 x 2")
    }
    else if (got != C) {
        FAIL = FAIL \ "complete matrix came back altered"
    }

    // ---- 2. missings drop, survivors keep order and content --------------
    // Rows 2 and 4 are incomplete; rows 1, 3, 5 must survive in that order.
    M = (0, 1 \
         ., 2 \
         1, 1 \
         2, . \
         2, 0)
    E = (0, 1 \
         1, 1 \
         2, 0)
    got = dwlsmv_listwise(M)
    if (rows(got) != 3) {
        FAIL = FAIL \ ("expected 3 complete rows, got " + strofreal(rows(got)))
    }
    else if (got != E) {
        FAIL = FAIL \ "survivors are not the right rows, or not in the original order"
    }

    // ---- 3. extended missings count as missing ---------------------------
    // .a through .z are what a real dataset carries. Mata sorts every missing
    // above every real, so one that survived would be counted into the top
    // category by any comparison downstream.
    M   = (0, 1 \ .a, 1 \ 1, 0)
    got = dwlsmv_listwise(M)
    if (rows(got) != 2) {
        FAIL = FAIL \ ("extended missing .a was not treated as missing: " +
                       strofreal(rows(got)) + " rows survived, expected 2")
    }
    M   = (0, 1 \ 1, .z \ 1, 0)
    got = dwlsmv_listwise(M)
    if (rows(got) != 2) {
        FAIL = FAIL \ ("extended missing .z was not treated as missing: " +
                       strofreal(rows(got)) + " rows survived, expected 2")
    }

    // ---- 4. the contract with the rest of Module 1 ------------------------
    // dwlsmv_catcounts() refuses a column with a missing in it. After
    // deletion the same data is acceptable, and the counts are the ones the
    // surviving rows imply -- two 0s and one 1 in column 1 of E above.
    got = dwlsmv_listwise((0 \ . \ 0 \ 1))
    if (dwlsmv_catcounts(got, 2) != (2 \ 1)) {
        FAIL = FAIL \ "catcounts did not accept a listwise-deleted column, or counted it wrong"
    }

    printf("test_listwise: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_listwise()))

* ---- 5. error paths --------------------------------------------------------
local nerr 0

* An all-missing matrix must abort. Returning an empty Y would surface as a
* division by N = 0 inside Module 4, a long way from the cause.
capture mata: dwlsmv_listwise((., 1 \ 2, .))
if _rc != 3498 {
    display as error "test_listwise: listwise on a matrix with no complete cases gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_listwise(J(0, 0, .))
if _rc != 3498 {
    display as error "test_listwise: listwise on an empty matrix gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

* The contract, from the other side: Stage 1 must still refuse a missing that
* never went through deletion. If this ever stops erroring, the guarantee that
* a skipped deletion cannot pass silently is gone.
capture mata: dwlsmv_catcounts((0 \ . \ 1), 2)
if _rc != 3498 {
    display as error "test_listwise: catcounts on an undeleted missing gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_listwise: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_listwise: deletion verified -- order, content, extended missings, and the contract with Stage 1"
exit 0
