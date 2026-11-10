*! dwlsmv -- the joint asymptotic covariance matrix Gamma (Module 4, spec 4.5)
*!
*! Executed under StataMP 18 on 2026-09-20 and it passes. On case_D it reproduces
*! lavInspect(fit, "gamma") to 1.4e-13 at lavaan's own rho-hat, and to 9.1e-10 end
*! to end -- the residual being lavaan's polychoric convergence, not arithmetic.
*!
*! The locked spec 5.1 entry point, and the project's highest-risk routine:
*! spec Section 6 rates a Gamma missing its threshold-propagation term as a
*! CRITICAL risk, because the result is not an error but a different estimator
*! that looks entirely plausible.
*!
*!     Gamma = (1/N) sum_i psi_i psi_i'
*!
*! under spec 4.4's normalisation, sqrt(N)(s-hat - s) -> N(0, Gamma), so
*! Var(s-hat) ~ Gamma/N. That matches lavInspect(fit, "gamma") -- N TIMES the
*! asymptotic covariance -- which is the cross-language target. Whichever
*! module assembles e(V) applies the 1/N explicitly, and must unit-test it
*! there; it is not applied here.
*!
*! EMPIRICAL, NEVER PLUG-IN (spec 4.5.4). It would be cheaper to evaluate the
*! closed forms of 4.5.1 and 4.5.3 at the estimates for the parts of Gamma that
*! have them. Do not. lavaan builds Gamma-hat by averaging estimated influence
*! functions, the two constructions differ at O(N^-1/2), and at realistic N a
*! plug-in disagrees with the primary compatibility target by an amount that
*! gets reported as a standard-error bug. The closed forms are unit-test
*! oracles and nothing else.
*!
*! quadcross(), NEVER cross(). This is the O(N*D^2) accumulation the project's
*! accumulator rule was written for. Measured on a bitwise-identical input,
*! cross() moves by up to ~126 ulps with the processor count and
*! non-monotonically, because MP varies reduction order with it; quadcross()
*! moves by at most about 1 ulp. The development machine has 6 cores and the
*! Stata 16 milestone machine has 4, so cross() here would make every milestone
*! run disagree with the reference for reasons unrelated to the version under
*! test -- and it would do so on spec 4.5.1's within-variable identity, which
*! is asserted at machine precision and documented as a coding error and never
*! noise. The cost is measured and is not worth re-arguing: ~4.1x cross(),
*! about one second on a large problem.
*!
*! quadcross(Psi, Psi) also computes only what symmetry requires, which is spec
*! 4.5.4's "compute only the lower triangle" -- taken by using the routine
*! built for it rather than by hand-rolling a triangular loop.
*!
*! WHAT THIS DOES NOT DO. It does not compute the pair fits -- spec 5.1 hands
*! them in. Module 2 estimates and Module 4 assembles, so the propagation term
*! appears exactly once, in dwlsmv_psi_rho_table(), and can be audited there.
*! It does not centre: see spec 4.5.5 and struct dwlsmv_stage1's psibar, which
*! carries the diagnostic rather than applying it. A non-zero psibar is a bug
*! report, not a reason to switch formulas -- except the rho entry of a
*! one-empty-cell pair (O6, spec Module 2, @sec-o6), nonzero at O(1) by
*! construction (a); Module 10 excludes it (dwlsmv_emptycell_pairs()).
*!
*! Requires dwlsmv_types, dwlsmv_koff, dwlsmv_pos_tau, dwlsmv_psi_rho_table and
*! dwlsmv_psi_matrix.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwls_gamma(real matrix Y,
                       struct dwlsmv_statmap scalar map,
                       struct dwlsmv_pairfit vector pairs,
                       real matrix psi_tau)
{
    real matrix Psi

    Psi = dwlsmv_psi_matrix(Y, map, pairs, psi_tau)

    return(quadcross(Psi, Psi) :/ rows(Y))
}

end
