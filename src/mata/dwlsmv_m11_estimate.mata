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
*! (dwlsmv_breaddiag), rec (e(recovery)) and rnote (the recovery note), all
*! of the reported fit: under metric "marker", D11's recovery replaces the
*! first fit when it returns "recovered".
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

    par = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    // D11's trigger reads the first run, so a marker fit takes no retry (spec Module 7)
    fit = (metric == "marker" ? dwlsmv_fit_noretry(st, par) : dwls_fit(st, par))
    dwls_inference(fit)

    // ---- D11: marker-metric recovery ------------------------------------
    // A and Bd at the first fit, computed once: the trigger reads them, and
    // they are the notes' values unless the recovered fit replaces it.
    A     = dwlsmv_admissible(fit.theta, par, st.map)
    Bd    = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    rec   = "none"
    rnote = ""
    if (metric == "marker") {
        rec = dwlsmv_m11_recover(canon, par, fit, A, Bd, parr, fitr, rnote)
        if (rec == "recovered") {
            par = parr
            fit = fitr
            A   = dwlsmv_admissible(fit.theta, par, st.map)
            Bd  = dwlsmv_breaddiag(fit.Delta, st.Gamma)
        }
    }
}

end
