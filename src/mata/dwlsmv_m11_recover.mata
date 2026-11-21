*! dwlsmv -- Module 11: D11's marker-metric recovery (dwlsmv_m11_recover)
*!
*! Executed under StataMP 18 on 2026-10-07: tests/unit/test_marker_recovery.do.
*! B3.3 change executed under StataMP 18 on 2026-11-20: tests/unit/test_recovery_so.do (all checks pass; mutants in its header).
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
*! A SECOND-ORDER MODEL (release B, B3.3). A factor scaled by its variance (marker[f] == 0 after an eligible
*! return) takes c_f = 1 in the transform. On a "failed" return, when such a factor exists, the first fit
*! converged, and the disturbance-variance check (dwlsmv_psi_check(), classes 2 and 3) finds a negative free
*! diagonal psi, the note is replaced by the owner's wording (2026-11-20): it names each negative variance,
*! says the restart failed and the first fit is retained, that the unit-variance parameterization cannot
*! represent it (SO4.4) although the first fit converged, and names the step that failed. Wrapped at 78
*! columns by dwlsmv_m11_recwrap(). The check is recomputed here from fit, with the same
*! band and classes as dwlsmv_psi_note()'s lines, so the signature is unchanged. "disturbance variance" when
*! the negative psi is an endogenous factor's (a beta or gamma row into it), "variance" when exogenous.
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

// 0: no free diagonal psi of fit is classed negative (class 2 or 3 of dwlsmv_psi_check()); 1: one is, and it
// is exogenous (a variance); 2: one is an endogenous factor's (a disturbance variance). Not a raw sign test:
// a psi inside its band (class 1) is zero within numerical error and does not count.
// Optional outputs (owner's wording, 2026-11-20: name the reason): given fac, the factor names in
// declaration order, what = one clause per negative psi joined by "; " ("F1 has a negative disturbance
// variance"), and n = how many.
real scalar dwlsmv_m11_recnegpsi(struct dwlsmv_fit scalar fit, | string colvector fac,
                                 string scalar what, real scalar n)
{
    struct dwlsmv_stage1 vector sv
    real colvector band
    real matrix C
    real scalar r, i, k, kind, endo

    sv   = fit.stage1
    band = dwlsmv_psi_band(fit.theta, fit.par, sv, fit.Delta)
    C    = dwlsmv_psi_check(fit.theta, fit.par, band)
    kind = 0
    what = ""
    n    = 0
    for (r = 1; r <= rows(C); r++) {
        if (C[r, 4] == 2 | C[r, 4] == 3) {
            k = C[r, 1]
            for (i = 1; i <= fit.par.nall; i++) {
                if (fit.par.ptype[i] == 2 & fit.par.row[i] == fit.par.col[i] & fit.par.free[i]) {
                    if (fit.par.fidx[i] == k) {
                        endo = any(((fit.par.ptype :== 3) :| (fit.par.ptype :== 8)) :& (fit.par.row :== fit.par.row[i]))
                        if (endo) {
                            kind = 2
                        }
                        else {
                            kind = max((kind, 1))
                        }
                        n = n + 1
                        if (args() >= 2) {
                            if (what != "") {
                                what = what + "; "
                            }
                            what = what + fac[fit.par.row[i]] + " has a negative " +
                                   (endo ? "disturbance variance" : "variance")
                        }
                        break
                    }
                }
            }
        }
    }
    return(kind)
}

// Word-wrap s at width w (no word is split); each line ends in char(10).
string scalar dwlsmv_m11_recwrap(string scalar s, real scalar w)
{
    string rowvector t
    string scalar out, line
    real scalar j

    t    = tokens(s)
    out  = ""
    line = ""
    for (j = 1; j <= cols(t); j++) {
        if (line == "") {
            line = t[j]
        }
        else if (strlen(line) + 1 + strlen(t[j]) <= w) {
            line = line + " " + t[j]
        }
        else {
            out  = out + line + char(10)
            line = t[j]
        }
    }
    if (line != "") {
        out = out + line + char(10)
    }
    return(out)
}

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
    real scalar rc, i, negpsi, nneg
    string scalar what, para, noun

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
    // B3.3: a second-order model (some marker[f] == 0), a converged first fit, a negative free diagonal psi.
    // The three tests are nested: Mata's & does not short-circuit and the check needs a converged fit.
    // The note is then the owner's wording (2026-11-20), naming the reason, in place of the generic one;
    // the step that stopped is still named (spec D11, Failure).
    negpsi = 0
    nneg   = 0
    what   = ""
    if (any(marker :== 0)) {
        if (fit.converged == 1) {
            negpsi = dwlsmv_m11_recnegpsi(fit, fac, what, nneg)
        }
    }
    if (negpsi > 0) {
        para = "The marker-metric fit produced an invalid solution: " + what + "."
        if (strpos(trig, "could not determine") > 0) {
            para = para + " It also could not determine all parameters reliably."
        }
        if (nneg == 1) {
            noun = "Its negative " + (negpsi == 2 ? "disturbance variance" : "variance")
        }
        else {
            noun = "Its negative variances"
        }
        para = para + " A restart from the unit-variance solution, transformed to the marker" +
               " metric, was attempted but failed. The first fit is retained. " + noun +
               " cannot be represented by the unit-variance parameterization, although the" +
               " first fit itself converged."
        note = dwlsmv_m11_recwrap(para, 78) +
               dwlsmv_m11_recwrap("The step that failed: " + why + ".", 78)
    }
    rmexternal("dwlsmv_m11_rec")
    return("failed")
}

end
