*! dwlsmv -- fit indices, option A: dCFI and dRMSEA across generated levels
*! (Module 12)
*!
*! Written 2026-10-21 from derivation @sec-fit-indices, Result FI3, with F6
*! frozen by the user 2026-10-21. No user-facing chain exists yet; this is
*! the routine a chain front end calls.
*!
*!     real matrix dwlsmv_fitidx_delta(real matrix X, real colvector skip)
*!
*! X: one row per level in chain order, columns (CFI, RMSEA), each level's
*! dwlsmv_fitidx() values against the one c1 baseline. skip: 1 for a level
*! skipped because d = 0 (spec check 8), 0 otherwise. Returns rows(X) x 2,
*! (dCFI, dRMSEA): row 1 missing; a skipped row missing (F6: shown, the difference blank);
*! any other row is its value minus the LAST NON-SKIPPED earlier row's.
*! A missing index gives a missing difference. Refused, rc 3498: X not two
*! columns, skip not 0/1 or of another length, the first level skipped.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_fitidx_delta(real matrix X, real colvector skip)
{
    real scalar L, k, last
    real matrix out

    L = rows(X)
    if (cols(X) != 2) {
        _error("dwlsmv_fitidx_delta: X must have two columns, CFI and RMSEA")
    }
    if (rows(skip) != L) {
        _error("dwlsmv_fitidx_delta: skip must have one entry per level")
    }
    if (L == 0) {
        return(J(0, 2, .))
    }
    if (missing(skip)) {
        _error("dwlsmv_fitidx_delta: skip must be 0 or 1")
    }
    if (any((skip :!= 0) :& (skip :!= 1))) {
        _error("dwlsmv_fitidx_delta: skip must be 0 or 1")
    }
    if (skip[1]) {
        _error("dwlsmv_fitidx_delta: the first level cannot be skipped")
    }
    out  = J(L, 2, .)
    last = 1
    for (k = 2; k <= L; k++) {
        if (skip[k]) {
            continue
        }
        out[k, .] = X[k, .] - X[last, .]
        last = k
    }
    return(out)
}

end
