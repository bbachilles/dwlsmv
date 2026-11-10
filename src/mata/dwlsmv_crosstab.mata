*! dwlsmv -- observed contingency table for a pair (Module 2)
*!
*! Executed under StataMP 18 on 2026-09-19, on all six case_D pairs and on a
*! hand-built table in tests/unit/test_loading.do.
*!
*! Returns K1 x K2 counts, row r holding category r-1 of the first variable.
*! The declared K1, K2 govern the shape, exactly as the index map is keyed on
*! declared C and does not shrink for an unobserved category (spec 4.2). A
*! table that resized itself would make this pair's cell grid disagree with
*! every other object built from the same K.
*!
*! An empty CELL is fine and is not an error: it contributes nothing to the
*! log-likelihood, and the model probability of that cell is still positive.
*! An empty CATEGORY is a different thing entirely -- that is the Section 4.5.6
*! degeneracy, and dwlsmv_thresholds() has already refused it upstream, since
*! a variable with an empty category has no finite threshold to pass in here.
*!
*! Validation is delegated to dwlsmv_catcounts() rather than repeated. Running
*! it on each column costs one extra pass and means the coding rules -- integer,
*! in 0..K-1, no missings -- are defined in exactly one place.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_crosstab(real colvector yi, real colvector yj,
                            real scalar K1, real scalar K2)
{
    real matrix n
    real scalar N, i

    if (rows(yi) != rows(yj)) {
        _error("dwlsmv_crosstab: the two variables have different lengths")
    }

    // Validates coding and aborts on anything out of range or missing.
    (void) dwlsmv_catcounts(yi, K1)
    (void) dwlsmv_catcounts(yj, K2)

    N = rows(yi)
    n = J(K1, K2, 0)
    for (i = 1; i <= N; i++) {
        n[yi[i] + 1, yj[i] + 1] = n[yi[i] + 1, yj[i] + 1] + 1
    }

    return(n)
}

end
