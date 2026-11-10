*! dwlsmv -- DIFFTEST: one adjacent pair, checks and statistic (Module 12)
*!
*!     real rowvector dwlsmv_dtpair(struct dwlsmv_parmap scalar par0,
*!                                  real colvector theta0,
*!                                  real scalar d, real scalar exempt,
*!                                  real scalar T0, real scalar T1,
*!                                  real scalar c, real scalar nrtol,
*!                                  real scalar df0,
*!                                  string scalar cc, real matrix K,
*!                                  string scalar inv0, string scalar inv1,
*!                                  string scalar name0, string scalar name1,
*!                                  real matrix Delta1, real matrix Delta0,
*!                                  real matrix Omega
*!                                  | string scalar note)
*!
*! Built 2026-10-14 (exemption gate step 3c; user, option (a), E-027,
*! 2026-10-14). dwlsmv_dtstep's inputs, then H1's and H0's stacked Jacobians
*! at their own estimates and Omega = dwlsmv_omega(st), the shared Stage 1's
*! blockdiag(Gamma_g / pi_g). Calls dwlsmv_dtstep; calls dwlsmv_difftest ONLY
*! when dwlsmv_dtstep returns status 1, with its Td (T = (N - G)F, the
*! convention of Omega's weighting, MG7(1)). Returns
*!
*!     (status, d, Td, eps, T_d3, p, a_d, b_d)
*!
*! the first four dwlsmv_dtstep's, the last four missing at status 0 (d = 0,
*! no test) and 2 (within the band; what is reported there is the front
*! end's). note is dwlsmv_dtstep's, passed through. Aborts are dwlsmv_dtstep's
*! and dwlsmv_difftest's, and one of its own: the statistic's d (q1 - q0 from
*! the Jacobians' columns) must be dwlsmv_dtstep's d, rc 3498 otherwise (the
*! caller passed Jacobians of another pair).
*!
*! Load order: dwlsmv_types, dwlsmv_dtstep and its requirements,
*! dwlsmv_difftest and its requirements.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_dtpair(struct dwlsmv_parmap scalar par0,
                             real colvector theta0,
                             real scalar d, real scalar exempt,
                             real scalar T0, real scalar T1,
                             real scalar c, real scalar nrtol,
                             real scalar df0,
                             string scalar cc, real matrix K,
                             string scalar inv0, string scalar inv1,
                             string scalar name0, string scalar name1,
                             real matrix Delta1, real matrix Delta0,
                             real matrix Omega,
                             | string scalar note)
{
    real rowvector st, r
    string scalar nt

    nt = ""
    st = dwlsmv_dtstep(par0, theta0, d, exempt, T0, T1, c, nrtol, df0, cc, K,
                       inv0, inv1, name0, name1, nt)
    if (args() >= 19) {
        note = nt
    }
    if (st[1] != 1) {
        return((st, J(1, 4, .)))
    }
    if (cols(Delta1) - cols(Delta0) != st[2]) {
        errprintf("dwlsmv_dtpair: the checks give d = %f but the Jacobians give q1 - q0 = %f - %f.\n",
                  st[2], cols(Delta1), cols(Delta0))
        _error("dwlsmv_dtpair: the Jacobians are not this pair's")
    }
    r = dwlsmv_difftest(Delta1, Delta0, Omega, st[3])
    return((st, r[1], r[3], r[4], r[5]))
}

end
