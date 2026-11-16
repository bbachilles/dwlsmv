*! dwlsmv -- conditional mode: the centring diagnostic psi-bar (spec 4.5.5)
*!
*! Executed 2026-11-10, suite green (drafted by a read-only agent, run by the
*! main session).
*!
*! psibar = (1/N) sum_i psi_i, D x 1, carried in struct dwlsmv_stage1 and
*! never asserted on. Every Stage 1 score equation is solved, so it is near 0
*! (@sec-condx: 3.7e-16 in the population). As in v1 it CANNOT detect a
*! missing propagation term: the naive influence function is also centred.
*! It catches derivative errors only. There is no O6 repair in this mode, so
*! no entry is nonzero by construction.
*!
*! Requires dwlsmv_cx_psi_matrix (and what it requires).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_cx_psibar(real matrix Y, real matrix X,
                                struct dwlsmv_statmap scalar map,
                                real colvector s)
{
    real matrix Psi

    Psi = dwlsmv_cx_psi_matrix(Y, X, map, s)

    return(quadcross(Psi, J(rows(Psi), 1, 1)) :/ rows(Y))
}

end
