*! dwlsmv -- Module 11: D11's marker-metric recovery (dwlsmv_m11_recover)
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*!
*!     string scalar dwlsmv_m11_recover(string scalar canon,
*!                                      struct dwlsmv_parmap scalar par,
*!                                      struct dwlsmv_fit    scalar fit,
*!                                      real matrix A, real matrix Bd,
*!                                      struct dwlsmv_parmap scalar parr,
*!                                      struct dwlsmv_fit    scalar fitr,
*!                                      string scalar note)
*!
*! Spec D11. Called by dwlsmv_m11_fit() under metric(marker) only, after the
*! first fit and its inference. par, fit: the first fit; A, Bd: its
*! dwlsmv_admissible() and dwlsmv_breaddiag() values, the caller's, reused.
*! Returns e(recovery):
*!   "none"        the trigger did not fire (dwlsmv_m11_rectrigger())
*!   "ineligible"  it fired; dwlsmv_m11_receligible() refused; no fit run
*!   "failed"      the recovery fits ran (dwlsmv_m11_recrun()) and failed,
*!                 or one of them stopped with an error
*!   "recovered"   parr, fitr = the final marker fit, inference done
*! parr and fitr are written only on "recovered"; otherwise the caller
*! reports its first fit unchanged (D11, Failure). note is the notes-block
*! text, "" on "none"; char(10)-terminated lines otherwise.
*!
*! THE ABORT BOUNDARY. The recovery fits run as _stata("mata:
*! dwlsmv_m11_recrun()", 1): an _error() inside returns its rc here and the
*! command continues with the first fit (decided by the user 2026-10-07).
*! Measured 2026-10-07 on StataMP 18 before this was built: a nested mata
*! call aborting with _error() returns 3498 to _stata() and the calling
*! function carries on. The abort's own message is not shown (nooutput);
*! the note names the step that stopped instead.
*!
*! A theta2 substitution in dwls_identify() (spec 4.8) is not a failure: the
*! note says the second identification start was used (decided by the user
*! 2026-10-07).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_m11_recover(string scalar canon,
                                 struct dwlsmv_parmap scalar par,
                                 struct dwlsmv_fit    scalar fit,
                                 real matrix A, real matrix Bd,
                                 struct dwlsmv_parmap scalar parr,
                                 struct dwlsmv_fit    scalar fitr,
                                 string scalar note)
{
    pointer() scalar R
    real colvector marker
    string colvector reasons, fac
    string scalar trig, nl, head, status, why, stage
    real scalar rc, i

    nl   = char(10)
    note = ""
    trig = dwlsmv_m11_rectrigger(fit, A, Bd)
    if (trig == "") {
        return("none")
    }
    head = "The marker-metric fit " + trig + ";"
    if (strlen(head) > 78) {
        head = "The marker-metric fit" + nl + trig + ";"
    }

    fac = tokens(substr(tokens(canon, ";")[1], 8, .))'
    if (!dwlsmv_m11_receligible(par, fac, marker, reasons)) {
        note = head + nl +
               "recovery through the unit-variance solution was not attempted, because" + nl +
               "the model is not eligible for it:" + nl
        for (i = 1; i <= rows(reasons); i++) {
            note = note + "    " + reasons[i] + nl
        }
        return("ineligible")
    }

    R = crexternal("dwlsmv_m11_rec")
    if (R == NULL) {
        R = findexternal("dwlsmv_m11_rec")
    }
    *R = asarray_create()
    asarray(*R, "canon",  canon)
    asarray(*R, "par",    par)
    asarray(*R, "marker", marker)
    asarray(*R, "stage",  "")
    asarray(*R, "status", "")
    asarray(*R, "why",    "")
    asarray(*R, "theta2", 0)

    rc     = _stata("mata: dwlsmv_m11_recrun()", 1)
    status = asarray(*R, "status")
    why    = asarray(*R, "why")
    if (rc) {
        stage  = asarray(*R, "stage")
        status = "failed"
        if (stage == "uvid") {
            why = "dwls_identify() refused the unit-variance specification"
        }
        else if (stage == "uvfit") {
            why = "the unit-variance fit stopped with an error"
        }
        else if (stage == "transform") {
            why = "the transform to the marker metric stopped with an error"
        }
        else if (stage == "mid") {
            why = "dwls_identify() refused the restarted marker-metric specification"
        }
        else if (stage == "mfit") {
            why = "the restarted marker-metric fit stopped with an error"
        }
        else if (stage == "minf") {
            why = "inference on the restarted fit stopped with an error"
        }
        else {
            why = "an internal error stopped it"
        }
        why = why + " (rc " + strofreal(rc) + ")"
    }
    if (status == "recovered") {
        parr = asarray(*R, "parr")
        fitr = asarray(*R, "fitr")
        note = head + nl +
               "it was restarted from the unit-variance solution, transformed to the" + nl +
               "marker metric." + nl
        if (asarray(*R, "theta2")) {
            note = note +
                   "Recovery used the second identification start because the transformed" + nl +
                   "start was rank deficient." + nl
        }
        rmexternal("dwlsmv_m11_rec")
        return("recovered")
    }
    if (status != "failed") {
        rmexternal("dwlsmv_m11_rec")
        _error("dwlsmv: internal error, the recovery returned no status")
    }
    note = head + nl +
           "a restart from the unit-variance solution, transformed to the marker" + nl +
           "metric, was attempted. Recovery attempt failed; first fit retained:" + nl +
           "    " + why + "." + nl
    rmexternal("dwlsmv_m11_rec")
    return("failed")
}

end
