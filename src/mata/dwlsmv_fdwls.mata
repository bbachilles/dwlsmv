*! dwlsmv -- the DWLS objective, its gradient and its Gauss-Newton Hessian (Mod. 7)
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_fit.do check 1 puts the gradient against a central difference of F at 2.1e-4 of a 1e-7 budget.
*!
*! An optimize() evaluator. Spec Section 5, Module 7:
*!
*!     F(theta) = (s - sigma(theta))' W (s - sigma(theta))
*!
*! with W = diag(w). dwls_fit() calls it twice: once with w = 1 (ULS, for
*! starting values) and once with w = 1 :/ diagonal(Gamma) (DWLS). Nothing here
*! knows which, and nothing here should.
*!
*! W IS A VECTOR AND IS NEVER FORMED AS A MATRIX. Every product below is a
*! weighted cross product, quadcross(X, w, Z) = X' diag(w) Z, so the D x D
*! diagonal matrix never exists. That is a memory argument at D ~ 1000, but it
*! is also the convention argument: 'w' is W_D = diag(Gamma)^-1 itself, the
*! spec's weight, and not its inverse. Muthen, du Toit & Spisic (1997) write
*! their fit function with W^-1, so their W is diag(Gamma); substituting one
*! for the other runs without complaint and fits a different estimator. The
*! caller forms w exactly once, in dwls_fit(), and this routine only multiplies
*! by it.
*!
*! THE DERIVATIVES. With r = s - sigma(theta) and Delta = d sigma / d theta':
*!
*!     g = dF/dtheta' = -2 r' W Delta              (exact, given Delta)
*!     H = 2 Delta' W Delta                        (Gauss-Newton, NOT exact)
*!
*! The gradient is the chain rule on Gate 2's Delta and nothing more. The
*! Hessian is labelled honestly: the exact Hessian carries a second term,
*! -2 sum_d w_d r_d d2 sigma_d / d theta d theta', which is dropped. That term
*! vanishes with the residual, and at a fitted model the residuals are
*! O(N^-1/2), so near the optimum H is the Hessian up to that order. It is
*! also the "bread" Delta' W_D Delta of Section 4.10's sandwich, so the matrix
*! Module 8 inverts is the one the optimizer steps with. Supplying it makes
*! optimize()'s "nr" a damped Gauss-Newton, which is the standard method for a
*! minimum-distance objective; the alternative, a "d1" evaluator, makes
*! optimize() build H by differentiating this gradient numerically, which costs
*! npar extra Delta evaluations per iteration. tools/time_fit.do measures both.
*!
*! What H being approximate does NOT affect: the point converged to. optimize()
*! stops on the gradient (nrtol measures g H^-1 g'), and the gradient is exact.
*! A wrong H costs iterations, never the answer.
*!
*! Every product goes through quadcross(), not cross() or '*', for the reason
*! The project's numerical rules record: cross() moves with the
*! processor count because MP varies reduction order, and an optimizer path is
*! the one place a last-ulp difference in F can compound into a different
*! stopping iterate. quadcross() is not immune in general: under StataMP 18 it
*! moved with the count on dwlsmv_traces's products (E-034,
*! 2026-10-28), though not on the polish's own operands in the 4 units probed.
*! Delta is D x npar, so the cost is small beside Delta
*! itself, which is 2 corrvec() evaluations per structural parameter.
*!
*! OUTPUTS: v, g and H are filled by reference -- v always, g when todo >= 1,
*! H when todo == 2 -- which is optimize()'s evaluator contract. A todo = 0
*! caller's g and H are left untouched.
*!
*! 'p' IS NOT WRITTEN. It is optimize()'s own parameter vector, passed by
*! address; theta is a copy. Neither is s, w, par nor map, and neither is
*! anything dwls_model() or dwls_delta() is handed -- checked routine by
*! routine on 2026-09-23 before the first loop was run over them.
*!
*! Signature (optimize()'s "d0/d1/d2" form, four extra arguments):
*!
*!     void dwlsmv_fdwls(real scalar todo, real rowvector p,
*!                       real colvector s, real colvector w,
*!                       struct dwlsmv_parmap  scalar par,
*!                       struct dwlsmv_statmap vector maps,
*!                       real scalar v, real rowvector g, real matrix H)
*!
*! MAPS (spec 4.8 O14, amended 2026-10-10): one index map per group. Under
*! groups s is stacked group by group and w is weighting (a), pi_g / diag(Gamma_g)
*! stacked, so v is the pooled F = sum_g pi_g F_g (@sec-multigroup-scaling MG2).
*! A single map passes where the vector is declared; G = 1 is unchanged.
*!
*! Load order: everything dwls_model() and dwls_delta() need, then those two.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_fdwls(real scalar todo, real rowvector p,
                  real colvector s, real colvector w,
                  struct dwlsmv_parmap  scalar par,
                  struct dwlsmv_statmap vector maps,
                  real scalar v, real rowvector g, real matrix H)
{
    real colvector theta, r
    real matrix Dm

    theta = p'
    r     = s - dwls_model(theta, par, maps)
    v     = quadcross(r, w, r)

    if (todo >= 1) {
        Dm = dwls_delta(theta, par, maps)
        g  = -2 * quadcross(w :* r, Dm)
        if (todo == 2) {
            H = 2 * quadcross(Dm, w, Dm)
        }
    }
}

end
