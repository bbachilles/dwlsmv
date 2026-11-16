*! dwlsmv -- conditional mode: Gamma in the Hessian form (Module 4, spec 5.1;
*! spec 4.8 "Stage 1 covariance"; @sec-condx Item 3)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*! tests/unit/test_condx_stage1.do compares every entry with
*! fixtures/expected/condx_ref_{A,D}_gamma.csv at 1e-8 abs + 1e-8 rel (the
*! prereg's [V6] tier; spec 7.2).
*!
*!     Gamma = (1/N) sum_i psi_i psi_i'
*!
*! N-divided, as v1's dwls_gamma() and as the analytic R reference
*! (99n_condx_ref.R: cx_gamma, sum_k (n_k / N) psi_k psi_k', the manifest's
*! "Gamma-hat: cx_gamma, N-divided"). EMPIRICAL: an average of estimated
*! influence functions, never a plug-in (spec 4.5.4).
*!
*! quadcross(), NEVER cross(): cross() moves with the processor count
*! (dwls_gamma.mata's header has the measurement; prereg.md, the phi2 clause).
*!
*! WHAT THIS DOES NOT DO. It does not estimate. It takes s, which holds every
*! Stage 1 estimate (item thresholds and slopes, residual correlations), and
*! the blocks are EVALUATED there by Modules 1 and 2 inside
*! dwlsmv_cx_psi_matrix(), which is the one place the propagation terms are
*! assembled. Whether s came from a converged Stage 1 is the caller's
*! business: dwls_sample_stats_cx() calls this only after every item and pair
*! has been accepted.
*!
*! Requires dwlsmv_cx_psi_matrix (and what it requires).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwls_gamma_cx(real matrix Y, real matrix X,
                          struct dwlsmv_statmap scalar map,
                          real colvector s)
{
    real matrix Psi

    Psi = dwlsmv_cx_psi_matrix(Y, X, map, s)

    return(quadcross(Psi, Psi) :/ rows(Y))
}

end
