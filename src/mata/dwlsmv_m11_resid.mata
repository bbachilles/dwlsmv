*! dwlsmv -- Module 11: the residual diagnostics of one fit, for posting
*!
*! Written by the main session 2026-11-22 (derivation @sec-resid-diag SR5; owner decisions
*! E-117, package 2, 3b-3c: O3, O4, O6).
*!
*!     void dwlsmv_m11_resid(struct dwlsmv_stage1 scalar st, struct dwlsmv_fit scalar fit,
*!                           struct dwlsmv_parmap scalar par, string rowvector ind)
*!
*! Called by dwlsmv_m11_fit() on the single-group, unconditional path, after the fit, while the
*! Stage 1 struct and the recoded Y still exist (dwlsmv_m11_clear() removes them before the
*! command returns, so a replay cannot recompute these; SR5). For a converged fit it sets the
*! caller's locals
*!   srmr      SRMR (dwlsmv_resid())
*!   resid_m   the name of a Stata temporary matrix: E, p x p, signed r - rho, diagonal 0
*!   nresid_m  the name of a Stata temporary matrix: Z, p x p, normalized, diagonal missing,
*!             missing at repaired pairs and pairs held at the bound (O1, O2)
*! with rows and columns named by the indicators in e(indicators) order. dwlsmv.ado posts them
*! as e(srmr), e(Resid), e(Nresid) under the fit indices' rule: converged, no covariate, one
*! group (SR5). An inadmissible converged fit posts them (O3); df = 0 posts the computed values
*! (O4). For a fit that did not converge the three locals are set empty and nothing is computed:
*! the last iterate is not a solution (SR5).
*!
*! rho is sigma(theta-hat)'s polychoric rows from dwls_model(); for efa() this is the unrotated
*! fit's, which equals Lambda Phi Lambda' of any rotation (SR6). The repaired pairs come from
*! dwlsmv_emptycell_pairs() on the recoded Y, as Stage 1 and the notes derive them.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_resid(struct dwlsmv_stage1 scalar st, struct dwlsmv_fit scalar fit,
                      struct dwlsmv_parmap scalar par, string rowvector ind)
{
    real colvector sig, rho, rep
    real matrix E, Z
    real scalar srmr, T, P, p
    string scalar en, zn
    string matrix stp

    st_local("srmr", "")
    st_local("resid_m", "")
    st_local("nresid_m", "")
    if (fit.converged == 0) {
        return
    }
    p = st.map.p
    T = st.map.T
    P = st.map.P
    if (cols(ind) != p) {
        errprintf("dwlsmv_m11_resid: %g indicator names for p = %g.\n", cols(ind), p)
        _error("dwlsmv_m11_resid: indicator names do not match p")
    }
    sig  = dwls_model(fit.theta, par, st.map)
    rho  = sig[|T + 1 \ T + P|]
    rep  = dwlsmv_emptycell_pairs(*findexternal("dwlsmv_m11_Y"), st.map)
    srmr = dwlsmv_resid(st, rho, rep, E, Z)

    stp = (J(p, 1, ""), ind')
    en  = st_tempname()
    zn  = st_tempname()
    st_matrix(en, E)
    st_matrixrowstripe(en, stp)
    st_matrixcolstripe(en, stp)
    st_matrix(zn, Z)
    st_matrixrowstripe(zn, stp)
    st_matrixcolstripe(zn, stp)
    st_local("srmr",     strofreal(srmr, "%21.0g"))
    st_local("resid_m",  en)
    st_local("nresid_m", zn)
}

end
