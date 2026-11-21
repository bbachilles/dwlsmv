*! dwlsmv -- Module 11: D11's unit-variance to marker transform
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*! B3.3 change executed under StataMP 18 on 2026-11-20: tests/unit/test_recovery_so.do (all checks pass; mutants in its header).
*!
*!     real scalar dwlsmv_m11_rectransform(struct dwlsmv_parmap scalar parm,
*!                                         struct dwlsmv_parmap scalar paru,
*!                                         real colvector thetau,
*!                                         real colvector marker,
*!                                         real colvector thetam,
*!                                         real colvector c)
*!
*! Spec D11 step 3; derivation @sec-marker-recovery Result 1; Mata port of
*! mr_recover()'s transform (chunk mr-impl), table in, table out, rows matched
*! by key (ptype, row, col) -- Module 3 owns each theta order, so nothing here
*! assumes the two tables share one.
*!
*!   parm    the requested MARKER table (the first fit's par)
*!   paru    the UNIT-VARIANCE table of the same model (dwls_identify() on
*!           dwlsmv_m11_recuvtext()'s text)
*!   thetau  the unit-variance fit's theta, in paru's order
*!   marker  dwlsmv_m11_receligible()'s output: marker[f] = indicator row;
*!           marker[f] == 0 is a factor scaled by its variance (a second-order G, release B, B3.3, derivation
*!           @sec-so-recovery SO4.1): c_f = 1, no lookup. The transform is then the U -> V map: beta_fg *
*!           c_f / 1, psi_fg * c_f * 1, and no lambda row of f exists to divide. psi_ff of such a factor is
*!           fixed at 1 in BOTH tables, so it is never mapped.
*!
*! With c_f = the marker loading's value in the unit-variance point:
*!     lambda_jf -> lambda_jf / c_f
*!     beta_fg   -> beta_fg * c_f / c_g      (f = row = to, g = col = from)
*!     psi_fg    -> psi_fg * c_f * c_g       (psi_ff = 1 fixed in paru -> c_f^2)
*!     gamma_fk  -> gamma_fk * c_f           (ptype 8, row = factor f, col = covariate k; conditional fits,
*!                                            derivation @sec-condx-recovery Result C1; owner 2026-11-15)
*!     kappa_jk  -> unchanged                (ptype 9: a direct effect on y*, which is not rescaled)
*!     tau       -> unchanged
*! Outputs (overwritten): thetam, npar x 1 in parm's theta order; c, nfac x 1.
*!
*! Returns 0 ok; 1 some c_f is 0 or missing (Mata's non-finite: the map does
*! not exist, spec D11 Failure); 2 a transformed value is missing (overflow).
*! NO near-zero cutoff: that gate was closed 2026-10-07 (spec D11, derivation
*! c_f -> 0); only an exact 0 or a missing c_f stops it.
*!
*! A key present in parm's free rows and absent from paru is an internal
*! error, rc 3498: both tables come from the same model text.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_m11_rectransform(struct dwlsmv_parmap scalar parm,
                                    struct dwlsmv_parmap scalar paru,
                                    real colvector thetau,
                                    real colvector marker,
                                    real colvector thetam,
                                    real colvector c)
{
    real colvector vu
    real scalar i, k, j, f, u, r, cc

    if (rows(thetau) != paru.npar | parm.npar != paru.npar | parm.nall != paru.nall) {
        errprintf("dwlsmv_m11_rectransform: marker table %f/%f rows, unit-variance %f/%f, theta %f.\n",
                  parm.npar, parm.nall, paru.npar, paru.nall, rows(thetau))
        _error("dwlsmv: internal error, the two metric tables differ in size")
    }

    // every row's value at the unit-variance point; a copy, paru is not touched
    vu = paru.value
    for (i = 1; i <= paru.nall; i++) {
        if (paru.free[i]) {
            vu[i] = thetau[paru.fidx[i]]
        }
    }

    // c_f: the marker loading (ptype 1, row marker[f], col f) in the UV point
    c = J(rows(marker), 1, .)
    for (f = 1; f <= rows(marker); f++) {
        if (marker[f] == 0) {
            // scaled by its variance (B3.3): c_f = 1, exactly
            c[f] = 1
            continue
        }
        k = 0
        for (i = 1; i <= paru.nall; i++) {
            if (paru.ptype[i] == 1) {
                if (paru.row[i] == marker[f]) {
                    if (paru.col[i] == f) {
                        k = i
                    }
                }
            }
        }
        if (k == 0) {
            _error("dwlsmv: internal error, a marker loading is not in the unit-variance table")
        }
        if (!paru.free[k]) {
            _error("dwlsmv: internal error, a marker loading is fixed in the unit-variance table")
        }
        c[f] = vu[k]
    }
    thetam = J(parm.npar, 1, .)
    for (f = 1; f <= rows(c); f++) {
        if (missing(c[f])) {
            return(1)
        }
        if (c[f] == 0) {
            return(1)
        }
    }

    for (i = 1; i <= parm.nall; i++) {
        if (parm.free[i]) {
            k = 0
            for (j = 1; j <= paru.nall; j++) {
                if (paru.ptype[j] == parm.ptype[i]) {
                    if (paru.row[j] == parm.row[i]) {
                        if (paru.col[j] == parm.col[i]) {
                            k = j
                        }
                    }
                }
            }
            if (k == 0) {
                errprintf("dwlsmv_m11_rectransform: %s has no row in the unit-variance table.\n",
                          parm.name[i])
                _error("dwlsmv: internal error, a marker row is missing from the unit-variance table")
            }
            u  = vu[k]
            r  = parm.row[i]
            cc = parm.col[i]
            if (parm.ptype[i] == 4) {
                thetam[parm.fidx[i]] = u
            }
            else if (parm.ptype[i] == 1) {
                thetam[parm.fidx[i]] = u / c[cc]
            }
            else if (parm.ptype[i] == 3) {
                thetam[parm.fidx[i]] = u * c[r] / c[cc]
            }
            else if (parm.ptype[i] == 8) {
                thetam[parm.fidx[i]] = u * c[r]
            }
            else if (parm.ptype[i] == 9) {
                thetam[parm.fidx[i]] = u
            }
            else if (parm.ptype[i] == 2) {
                thetam[parm.fidx[i]] = u * c[r] * c[cc]
            }
            else {
                // ptypes 5-7 (nu, omega, alpha) cannot occur in a single-group table; none is mapped
                errprintf("dwlsmv_m11_rectransform: %s has ptype %f, which has no transform.\n",
                          parm.name[i], parm.ptype[i])
                _error("dwlsmv: internal error, a parameter type has no marker transform")
            }
        }
    }
    if (missing(thetam)) {
        return(2)
    }
    return(0)
}

end
