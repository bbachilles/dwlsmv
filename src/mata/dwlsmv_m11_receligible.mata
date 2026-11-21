*! dwlsmv -- Module 11: is the marker model eligible for D11's recovery?
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*! B3.3 change executed under StataMP 18 on 2026-11-20: tests/unit/test_recovery_so.do (all checks pass; mutants in its header).
*!
*!     real scalar dwlsmv_m11_receligible(struct dwlsmv_parmap scalar par,
*!                                        string colvector fac,
*!                                        real colvector marker,
*!                                        string colvector reasons)
*!
*! Spec D11, Eligibility; derivation @sec-marker-recovery Result 2 (E1-E3);
*! Mata port of mr_eligibility() (chunk mr-impl). All-or-nothing, the user's
*! rule: returns 1 only if
*!   E1  every factor has exactly one fixed nonzero loading, equal to 1;
*!   E2  every factor has one psi_ff row and it is free;
*!       EXCEPT (release B, B3.3, spec Module 12) a factor with NO lambda row at all (free or fixed), a
*!       second-order factor G scaled by 'metric unitvariance G': E1 and E2 do not apply to it; instead
*!       it must have exactly one psi_GG row, FIXED at exactly 1 (it is scaled by its variance, c_G = 1,
*!       derivation @sec-so-recovery SO4.1). Any other psi_GG on a factor with no lambda row (free, or
*!       fixed at another value) is its own reason. A factor with a lambda row, an item on G included,
*!       is an ordinary factor under E1/E2;
*!   E3 every other fixed lambda, beta or psi (off-diagonal) is exactly 0;
*!   E4  (conditional fits, 2026-11-15, owner; derivation @sec-condx-recovery Result C3) every fixed gamma
*!       (ptype 8, x -> factor) is exactly 0. A gamma fixed at v in the marker table is c_f*v at the
*!       unit-variance optimum, and c_f is unknown before the fit. A fixed kappa (ptype 9) of any value
*!       is the same parameter in both metrics and stays eligible.
*! Read from Module 3's MARKER parameter table (the first fit's par), never
*! from the canonical text. Every reason is collected, not the first, so the
*! note can name them all.
*!
*! Outputs (by reference, overwritten): marker[f] = the indicator row of
*! factor f's marker loading, 0 where E1 fails; reasons = one line each.
*! When the function returns 1, marker[f] == 0 means "scaled by its variance, c_f = 1" (B3.3): dwlsmv_m11_
*! rectransform() takes c_f = 1 there. When it returns 0, marker is all zeros and means nothing.
*! par is read only.
*!
*! Fixed values reach the table from canonical text written by
*! dwlsmv_tr_num() (exact round trip), so == 1 and != 0 are exact tests.
*! A missing fixed value is its own reason; it is excluded from E1/E3's
*! nonzero tests explicitly, since missing != 0 is true in Mata.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_m11_receligible(struct dwlsmv_parmap scalar par,
                                   string colvector fac,
                                   real colvector marker,
                                   string colvector reasons)
{
    real scalar f, i, nlam, vlam, rlam, npsi, psifree, psival, nlrow, offd

    marker  = J(rows(fac), 1, 0)
    reasons = J(0, 1, "")

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] != 4) {
            if (!par.free[i]) {
                if (missing(par.value[i])) {
                    reasons = reasons \ (par.name[i] + " is fixed at a missing value")
                }
            }
        }
    }

    // ---- E1, E2: per factor ------------------------------------------------
    for (f = 1; f <= rows(fac); f++) {
        nlam    = 0
        vlam    = .
        rlam    = 0
        npsi    = 0
        psifree = 0
        psival  = .
        nlrow   = 0
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] == 1) {
                if (par.col[i] == f) {
                    nlrow = nlrow + 1
                    if (!par.free[i]) {
                        if (!missing(par.value[i])) {
                            if (par.value[i] != 0) {
                                nlam = nlam + 1
                                vlam = par.value[i]
                                rlam = par.row[i]
                            }
                        }
                    }
                }
            }
            else if (par.ptype[i] == 2) {
                if (par.row[i] == f) {
                    if (par.col[i] == f) {
                        npsi    = npsi + 1
                        psifree = par.free[i]
                        psival  = par.value[i]
                    }
                }
            }
        }
        if (nlrow == 0) {
            // B3.3: no lambda row at all: scaled by its variance, c_f = 1, marker[f] stays 0. The variance
            // row must be fixed at exactly 1. A missing fixed value is already a reason (first loop), so
            // it is not tested here (missing != 1 is true in Mata).
            if (npsi != 1) {
                reasons = reasons \ ("factor " + fac[f] + " has no single variance row")
            }
            else if (psifree) {
                reasons = reasons \ ("factor " + fac[f] + " has no indicator and its variance is free")
            }
            else if (!missing(psival)) {
                if (psival != 1) {
                    reasons = reasons \ ("factor " + fac[f] + " has no indicator and its variance is fixed at " +
                                         strtrim(strofreal(psival, "%21.0g")) + ", not 1")
                }
            }
        }
        else {
            if (nlam == 0) {
                reasons = reasons \ ("factor " + fac[f] + " has no loading fixed at 1")
            }
            else if (nlam > 1) {
                reasons = reasons \ ("factor " + fac[f] + " has " + strofreal(nlam) +
                                     " loadings fixed at nonzero values")
            }
            else if (vlam != 1) {
                reasons = reasons \ ("factor " + fac[f] + " has its loading fixed at " +
                                     strtrim(strofreal(vlam, "%21.0g")) + ", not 1")
            }
            else {
                marker[f] = rlam
            }
            if (npsi != 1) {
                reasons = reasons \ ("factor " + fac[f] + " has no single variance row")
            }
            else if (!psifree) {
                reasons = reasons \ ("factor " + fac[f] + " has its variance fixed")
            }
        }
    }

    // ---- E4: a fixed gamma must be exactly 0 (ptype 8; a fixed kappa, ptype 9, is not tested) -------
    // A missing fixed value is its own reason in the first loop above, so it is not repeated here.
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 8) {
            if (!par.free[i]) {
                if (!missing(par.value[i])) {
                    if (par.value[i] != 0) {
                        reasons = reasons \ (par.name[i] + " is fixed at " +
                                             strtrim(strofreal(par.value[i], "%21.0g")) +
                                             ", not 0 (the effect of the covariate on a factor)")
                    }
                }
            }
        }
    }

    // ---- E3: fixed nonzero beta or off-diagonal psi ---------------------------
    // (a fixed nonzero lambda beyond the marker is already E1's count above)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) {
            offd = 0
            if (par.ptype[i] == 3) {
                offd = 1
            }
            else if (par.ptype[i] == 2) {
                if (par.row[i] != par.col[i]) {
                    offd = 1
                }
            }
            if (offd) {
                if (!missing(par.value[i])) {
                    if (par.value[i] != 0) {
                        reasons = reasons \ (par.name[i] + " is fixed at " +
                                             strtrim(strofreal(par.value[i], "%21.0g")) + ", not 0")
                    }
                }
            }
        }
    }

    if (rows(reasons) > 0) {
        marker = J(rows(fac), 1, 0)
    }
    return(rows(reasons) == 0)
}

end
