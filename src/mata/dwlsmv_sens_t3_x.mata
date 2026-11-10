*! dwlsmv -- B_T3 with the closed-form dDelta/dtheta and J of Result MS6: the opt-in exact path
*! f must be an UNPOLISHED fit; see dwlsmv_sens.mata.
*!
*! Executed under StataMP 16.1 on 2026-10-19: tests/unit/test_sens_hess.do green; suite 62/0/0 at four processors.
*!
*!     real rowvector dwlsmv_sens_t3_x(struct dwlsmv_fit scalar f, real colvector step,
*!                                     real scalar epsV, real scalar rhoT,
*!                                     | real colvector ds, real colvector eps,
*!                                     real matrix dOm)
*!
*! dwlsmv_sens_t3's arguments, columns and refusals (read its header). The one
*! difference: beta_2 takes Theta from dwlsmv_sens_x, v from
*! dwlsmv_sens_ddelta_x and J from dwlsmv_sens_jmat_x, all exact on Result MS6's
*! domain (@sec-mg-hessian; one factor, no beta, one row per theta element per
*! group), instead of their finite-difference versions. Outside the domain it
*! refuses at rc 3498, with no fallback. dwlsmv_sens_t3's default path is
*! unchanged. Diagnostic only; no production path calls it, and
*! test_sens_t3.do does not use it (plan steps 4-5 decide that).
*!
*! Load order: dwlsmv_sens_t3 (dwlsmv_sens_t3_core), dwlsmv_sens_x,
*! dwlsmv_sens_jmat_x, dwlsmv_sens_ddelta_x.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_sens_t3_x(struct dwlsmv_fit scalar f, real colvector step,
                                real scalar epsV, real scalar rhoT,
                                | real colvector ds, real colvector eps,
                                real matrix dOm)
{
    if (args() == 4) {
        return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 1, 4, J(0, 1, .), J(0, 1, .), J(0, 0, .)))
    }
    if (args() == 5) {
        return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 1, 5, ds, J(0, 1, .), J(0, 0, .)))
    }
    if (args() == 6) {
        return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 1, 6, ds, eps, J(0, 0, .)))
    }
    return(dwlsmv_sens_t3_core(f, step, epsV, rhoT, 1, 7, ds, eps, dOm))
}

end
