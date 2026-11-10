*! dwlsmv -- Module 11: the fits of D11's recovery (dwlsmv_m11_recrun)
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*!
*!     void dwlsmv_m11_recrun()
*!
*! Spec D11, Procedure steps 1-4. Run ONLY by dwlsmv_m11_recover(), through
*! _stata("mata: dwlsmv_m11_recrun()"), so that an _error() anywhere inside
*! -- dwls_identify() refusing a text, dwls_fit(), dwls_inference(), an
*! internal error -- returns a nonzero rc to the caller instead of ending the
*! command (the user's decision 2026-10-07: a recovery-time abort is a failed
*! recovery, not a failed command). Mata has no catch; a nested Stata call is
*! the only boundary that stops an abort. Hence no arguments: everything
*! passes through the asarray external dwlsmv_m11_rec, which the caller
*! creates and dwlsmv_m11_clear() removes.
*!
*!   in   "canon"    the marker canonical text, ";"-separated
*!        "par"      the first fit's marker table
*!        "marker"   dwlsmv_m11_receligible()'s marker rows
*!   out  "stage"    the step being run, written BEFORE it runs, so that
*!                   after an abort the caller knows where it stopped
*!        "status"   "recovered" or "failed", written only on a return
*!        "why"      the failure's reason, for the note
*!        "parr", "fitr"  the final marker fit's table and fit (inference
*!                   done), written only on "recovered"
*!        "theta2"   1 if dwls_identify() replaced the transformed starts by
*!                   its second point (spec 4.8), 0 otherwise
*!
*! theta2 IS DETECTED, NOT ASSUMED (the user's decision 2026-10-07): the
*! requested starts are the transformed values written into the text; the
*! actual starts are par.value after dwls_identify(). They differ only when
*! Delta was rank deficient at the requested starts and theta2 passed.
*! dwlsmv_m11_recstarts() guarantees each written start reads back as the
*! identical double, so the comparison is exact. Thresholds are not written
*! (the grammar has no tau directive) and are not compared.
*!
*! The Stage 1 output is the session's (dwlsmv_m11_st): Stage 1 is not re-run.
*!
*! THE FINAL FIT SKIPS ULS (the user's decision 2026-10-07; spec D11 and
*! Module 7). It is dwlsmv_fit_from() from dwlsmv_m11_recstart(): park's own
*! starts (the transformed values, or theta2's) for lambda/beta/psi, and
*! Stage 1's sample thresholds for tau. Delta's threshold block is the
*! identity and sigma_tau = tau, so the DWLS optimum's thresholds ARE the
*! sample thresholds; measured 2026-10-07 on the trap and case E
*! (max|tau-hat - s_tau| = 0 exactly): 2 DWLS iterations against 7 after
*! ULS on the trap, F 1.3e-14 lower (E-016).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The final fit's start: par.value over the free rows, as dwls_fit() builds
// it, with each threshold at Stage 1's sample threshold. The table's first T
// rows are the tau rows in the index map's order (Section 4.2), checked.
real colvector dwlsmv_m11_recstart(struct dwlsmv_stage1 scalar st,
                                   struct dwlsmv_parmap  scalar par)
{
    real colvector start
    real scalar i

    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            start[par.fidx[i]] = par.value[i]
        }
    }
    if (par.nall < st.map.T) {
        _error("dwlsmv: internal error, the table has fewer rows than thresholds")
    }
    for (i = 1; i <= st.map.T; i++) {
        if (par.ptype[i] != 4 | par.free[i] != 1) {
            _error("dwlsmv: internal error, a threshold row is not where expected")
        }
        start[par.fidx[i]] = st.s[i]
    }
    if (missing(start)) {
        _error("dwlsmv: internal error, the final fit's start is incomplete")
    }
    return(start)
}

void dwlsmv_m11_recrun()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par, paru, park
    struct dwlsmv_fit    scalar fitu, fitk
    pointer() scalar R
    real colvector marker, thetam, c
    string scalar canon, nl, why
    string colvector fac
    real scalar rc, f, i, t2

    nl     = char(10)
    R      = findexternal("dwlsmv_m11_rec")
    st     = *findexternal("dwlsmv_m11_st")
    canon  = asarray(*R, "canon")
    par    = asarray(*R, "par")
    marker = asarray(*R, "marker")
    fac    = tokens(substr(tokens(canon, ";")[1], 8, .))'

    // 1-2. the unit-variance fit on the same Stage 1 output
    asarray(*R, "stage", "uvid")
    paru = dwls_identify(subinstr(dwlsmv_m11_recuvtext(canon), ";", nl), st.map)
    asarray(*R, "stage", "uvfit")
    fitu = dwlsmv_fit_noretry(st, paru)
    if (fitu.converged != 1) {
        asarray(*R, "why", "the unit-variance fit did not converge")
        asarray(*R, "status", "failed")
        return
    }

    // 3. the transform
    asarray(*R, "stage", "transform")
    rc = dwlsmv_m11_rectransform(par, paru, fitu.theta, marker, thetam, c)
    if (rc == 1) {
        why = ""
        for (f = 1; f <= rows(c); f++) {
            if (missing(c[f]) | c[f] == 0) {
                why = why + (why == "" ? "" : "; ") + "the marker loading of " + fac[f] +
                      " is " + (missing(c[f]) ? "not finite" : "0") + " under unit variance"
            }
        }
        asarray(*R, "why", why)
        asarray(*R, "status", "failed")
        return
    }
    if (rc == 2) {
        asarray(*R, "why", "a transformed starting value is not finite")
        asarray(*R, "status", "failed")
        return
    }

    // 4. the marker text with the transformed starts, through the normal path
    asarray(*R, "stage", "mid")
    park = dwls_identify(subinstr(dwlsmv_m11_recstarts(canon, par, thetam), ";", nl), st.map)
    t2 = 0
    for (i = 1; i <= park.nall; i++) {
        if (park.free[i]) {
            if (park.ptype[i] != 4) {
                if (park.value[i] != thetam[park.fidx[i]]) {
                    t2 = 1
                }
            }
        }
    }
    asarray(*R, "theta2", t2)
    asarray(*R, "stage", "mfit")
    fitk = dwlsmv_fit_from(st, park, dwlsmv_m11_recstart(st, park))
    if (fitk.converged != 1) {
        asarray(*R, "why", "the restarted marker-metric fit did not converge")
        asarray(*R, "status", "failed")
        return
    }
    asarray(*R, "stage", "minf")
    dwls_inference(fitk)
    asarray(*R, "parr", park)
    asarray(*R, "fitr", fitk)
    asarray(*R, "status", "recovered")
}

end
