*! dwlsmv -- listwise deletion (Module 1, spec Section 3)
*!
*! Executed under StataMP 18 on 2026-09-19 and it passes: tests/unit/
*! test_listwise.do drives it on hand-built matrices, including the extended
*! missings .a and .z, and asserts that survivors keep their original order and
*! content rather than merely their count.
*!
*! Returns the rows of Y that are complete on every column. Spec Section 3
*! fixes listwise as v1's missing-data handling; this is that, and nothing
*! more. No pairwise deletion, no imputation, no per-pair sample size -- each
*! of those changes what Gamma is the asymptotic covariance OF, and none is in
*! v1 scope.
*!
*! WHY THIS IS A SEPARATE ROUTINE AND NOT A STEP INSIDE dwls_sample_stats().
*! Deleting rows changes N, and N is what every influence function is divided
*! by. A deletion buried inside Stage 1 would make the N that Gamma was built
*! with invisible to the caller reporting it. Keeping it separate means the
*! caller performs the deletion, sees the row count change, and hands Stage 1
*! a matrix whose rows() is the N it will report.
*!
*! THE WRAPPER STILL OWNS THE AUTHORITATIVE DELETION. A Stata estimation
*! command is expected to mark its estimation sample, and `e(sample)` cannot be
*! set from a Mata matrix that has already dropped the rows -- by then the
*! correspondence to observation numbers is gone. So Module 11 does the
*! deletion at the dataset level with a touse marker, and this routine is the
*! backstop for every other caller: tests, tools, and anything driving Stage 1
*! directly from a matrix. Both paths must agree, which they do trivially,
*! because "complete on every column" is not a judgement call.
*!
*! The complementary count is deliberately not returned. It is rows(Y) minus
*! rows of the result, which the caller has both of -- a second copy that can
*! disagree with the first is what the project's struct-field rule is about,
*! and it applies to return values for the same reason.
*!
*! select() and rowmissing() were checked against the version floor before
*! being used, by the whatsnew procedure now recorded in the allowlist header:
*! neither appears in Stata's whatsnew16to17 or whatsnew17to18, so neither was
*! introduced after the floor.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_listwise(real matrix Y)
{
    real matrix keep

    if (rows(Y) < 1 | cols(Y) < 1) {
        _error("dwlsmv_listwise: Y is empty")
    }

    // This is StataCorp's own documented idiom, not a construction of ours.
    // [M-5] select() example 5 reads, verbatim:
    //
    //     submat = select(X, rowmissing(X):==0)
    //
    // under the heading "To select rows of X for which there are no missing
    // values" -- which is listwise deletion. Using the documented form rather
    // than a hand-rolled loop means rowmissing() owns the definition of
    // "missing", so the extended codes .a through .z are handled by Stata
    // rather than by a comparison we wrote. tests/unit/test_listwise.do
    // asserts .a and .z anyway, because that is the assumption most worth
    // holding to a test rather than to a reading of the manual.
    keep = select(Y, rowmissing(Y) :== 0)

    if (rows(keep) == 0) {
        _error("dwlsmv_listwise: no complete cases")
    }

    return(keep)
}

end
