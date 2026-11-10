*! dwlsmv -- Module 11: is the marker model eligible for D11's recovery?
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
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
*!   E3  every other fixed lambda, beta or psi (off-diagonal) is exactly 0.
*! Read from Module 3's MARKER parameter table (the first fit's par), never
*! from the canonical text. Every reason is collected, not the first, so the
*! note can name them all.
*!
*! Outputs (by reference, overwritten): marker[f] = the indicator row of
*! factor f's marker loading, 0 where E1 fails; reasons = one line each.
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
    real scalar f, i, nlam, vlam, rlam, npsi, psifree, offd

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
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] == 1) {
                if (par.col[i] == f) {
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
                    }
                }
            }
        }
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
