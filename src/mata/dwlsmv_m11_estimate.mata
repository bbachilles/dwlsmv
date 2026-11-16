*! dwlsmv -- Module 11: one single-group fit, recovery included (dwlsmv_m11_estimate)
*!
*! Factored out of dwlsmv_m11_fit() on 2026-10-23 so that difftest's refit
*! (spec Module 12 check 6) runs the same code as the fit it reproduces.
*!
*!     void dwlsmv_m11_estimate(string scalar canon, string scalar metric,
*!                              struct dwlsmv_stage1 scalar st,
*!                              struct dwlsmv_parmap scalar par,
*!                              struct dwlsmv_fit    scalar fit,
*!                              real matrix A, real matrix Bd,
*!                              string scalar rec, string scalar rnote)
*!
*! canon is the single-line, ";"-separated canonical text dwlsmv.ado passes.
*! Outputs: par, fit (inference done), A (dwlsmv_admissible), Bd
*! (dwlsmv_m11_bread: dwlsmv_breaddiag, or "not evaluated" for a Delta with a missing value), rec (e(recovery)) and rnote (the recovery note), all
*! of the reported fit: under metric "marker", D11's recovery replaces the
*! first fit when it returns "recovered". A conditional fit (one covariate) takes the same path (2026-11-15);
*! its Bd is dwlsmv_breaddiag() too, and the bread-rank term is part of its trigger (owner, 2026-11-15; see the body).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_estimate(string scalar canon, string scalar metric,
                         struct dwlsmv_stage1 scalar st,
                         struct dwlsmv_parmap scalar par,
                         struct dwlsmv_fit    scalar fit,
                         real matrix A, real matrix Bd,
                         string scalar rec, string scalar rnote)
{
    struct dwlsmv_parmap scalar parr
    struct dwlsmv_fit    scalar fitr
    real scalar cx, bread_cx

    par = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    // D11's trigger reads the first run, so a marker fit takes no retry (spec Module 7). That holds with a
    // covariate as well (owner, 2026-11-15, decision 3); a conditional fit under unitvariance takes
    // dwls_fit, as every unit-variance fit does.
    cx  = (dwlsmv_statmap_q(st.map) > 0)
    if (metric == "marker") {
        fit = dwlsmv_fit_noretry(st, par)
    }
    else {
        fit = dwls_fit(st, par)
    }
    dwls_inference(fit)

    // ===== DECIDED (owner, 2026-11-15, "Yes, both"): the bread-rank row in a CONDITIONAL fit ==============
    // dwlsmv_breaddiag() runs on a conditional Delta and Gamma, and its row 3 is part of the conditional D11
    // trigger as for an ordinary fit, and e(bread_ratio) is posted. Evidence: production's conditional bread
    // ratio 0.15 / 0.12 healthy against 4e-13 / 1.2e-13 for an exact scale twin (condx_A / condx_D, StataMP 18,
    // E-097); several factors, B and Stata 16.1 are not measured. The
    // switch stays a variable so that a measurement that fails can turn it off in one place.
    bread_cx = 1
    // =====================================================================================================

    // ---- D11: marker-metric recovery ------------------------------------
    // A and Bd at the first fit, computed once: the trigger reads them, and
    // they are the notes' values unless the recovered fit replaces it. A is
    // dwlsmv_admissible() at the REPORTED theta, converged or not (it accepts a conditional map of one group).
    A     = dwlsmv_admissible(fit.theta, par, st.map)
    Bd    = J(0, 0, .)
    if (!cx | bread_cx) {
        Bd = dwlsmv_m11_bread(fit.Delta, st.Gamma)
    }
    rec   = "none"
    rnote = ""
    if (metric == "marker") {
        rec = dwlsmv_m11_recover(canon, par, fit, A, Bd, parr, fitr, rnote)
        if (rec == "recovered") {
            par = parr
            fit = fitr
            A   = dwlsmv_admissible(fit.theta, par, st.map)
            Bd  = J(0, 0, .)
            if (!cx | bread_cx) {
                Bd = dwlsmv_m11_bread(fit.Delta, st.Gamma)
            }
        }
    }
}

end
