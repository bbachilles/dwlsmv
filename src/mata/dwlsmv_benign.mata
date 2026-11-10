*! dwlsmv -- Module 11: D1's dummy-data assertion
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_translate.do.
*! Called by _dwlsmv_translate.ado.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// D1's dummy data, asserted rather than assumed: no missing value, no zero
// variance, every |r| < 0.9, and a well-conditioned correlation matrix
void dwlsmv_benign(string rowvector v)
{
    real matrix X, R
    real colvector ev

    X = st_data(., v)
    if (missing(X)) {
        _error("dwlsmv_benign: missing values in the dummy data")
    }
    if (any(diagonal(variance(X)) :<= 0)) {
        _error("dwlsmv_benign: zero variance in the dummy data")
    }
    R = correlation(X)
    if (cols(X) > 1) {
        if (max(abs(R - I(cols(X)))) >= 0.9) {
            _error("dwlsmv_benign: a dummy correlation at 0.9 or more")
        }
    }
    ev = symeigenvalues(R)'
    if (min(ev) / max(ev) < 1e-3) {
        _error("dwlsmv_benign: ill-conditioned dummy correlation matrix")
    }
}


end
