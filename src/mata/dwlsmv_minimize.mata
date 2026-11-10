*! dwlsmv -- minimize (s - sigma)' diag(w) (s - sigma) from a start (Module 7)
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_fit.do checks 2 and 8: convergence on five fits, and a one-iteration cap returning 0 with the last iterate.
*!
*! One weighted least-squares minimization, run by optimize() over
*! dwlsmv_fdwls(). dwls_fit() calls it twice -- w = 1 for ULS, then
*! w = diag(Gamma)^-1 for DWLS starting from the ULS solution -- and that
*! sequencing, which is spec Section 5's, lives in dwls_fit() and not here.
*!
*!     real scalar dwlsmv_minimize(real colvector s, real colvector w,
*!                                 struct dwlsmv_parmap  scalar par,
*!                                 struct dwlsmv_statmap vector maps,
*!                                 real colvector start,
*!                                 string scalar evaltype,
*!                                 string scalar technique,
*!                                 real scalar maxiter,
*!                                 real colvector theta,
*!                                 real scalar fmin,
*!                                 real scalar iter)
*!
*! Returns 1 if optimize() declared convergence, 0 otherwise. theta, fmin and
*! iter are filled by reference, and are filled IN BOTH CASES: on failure they
*! hold the last iterate, its objective and the iteration count, so a caller
*! can report where the optimizer stopped rather than only that it did.
*!
*! NON-CONVERGENCE IS A RETURN VALUE, NOT AN ABORT. _optimize() is used rather
*! than optimize() for the reason Module 2 uses it: the caller decides what an
*! unconverged fit means. An iteration limit, a flat region or a step that
*! cannot improve F all come back here as 0. What does still abort is an
*! evaluator that cannot evaluate at all -- dwls_model()'s unfilled-sigma check,
*! and dwlsmv_rystar()'s singular (I - B). The second can only be reached by a
*! non-recursive B whose trial step lands on a singular point; a recursive B has
*! a unit-triangular (I - B) at every theta. That abort names the cause, which
*! is better than a silent missing, but it is not a converged = 0 and nobody
*! should read it as one.
*!
*! 'evaltype' and 'technique' are arguments, not constants, so that
*! tools/time_fit.do can reproduce the measurement that chose dwls_fit()'s
*! pair. They are passed straight to optimize(): "d1" or "d2", and "nr",
*! "bfgs", ... . With "d2" the Hessian is dwlsmv_fdwls()'s Gauss-Newton one.
*!
*! CONVERGENCE TOLERANCES, and why they are tighter than optimize()'s defaults
*! (ptol 1e-6, vtol 1e-7, nrtol 1e-5). The default ptol is a RELATIVE change of
*! 1e-6 in theta between iterations, which is the same order as spec Section
*! 7.2's theta-hat tier against lavaan (1e-5 abs, 1e-4 rel). A comparison
*! against an external reference is only informative if our own stopping error
*! is far inside the tolerance, or a disagreement cannot be attributed. The
*! equivalent-models identity (spec 4.8) is sharper still: two parameterizations
*! of the same model, separately converged, compared with each other. So:
*!
*!     ptol  1e-10   relative change in theta
*!     vtol  1e-12   relative change in F
*!     nrtol 1e-12   g H^-1 g', the Newton decrement
*!
*! optimize() requires (ptol or vtol) AND nrtol. nrtol is the one that
*! certifies a stationary point, and with the Gauss-Newton H it is twice the
*! decrease in F a full step would predict. Their cost over the defaults has
*! NOT been measured -- they are hard-coded here, so tools/time_fit.do cannot
*! vary them. What it does show is that the fits converge in 4-7 iterations
*! under them.
*!
*! 'start', 'theta', 's', 'w', 'par' and 'map' are distinct variables: start is
*! handed to optimize() as a transposed copy, and theta is written only once,
*! at the end, from optimize()'s result. A caller may pass the same variable as
*! start and theta.
*!
*! MAPS: one per group (spec 4.8 O14, amended 2026-10-10); s and w are stacked
*! and D is the sum of the maps' D. A single map passes where the vector is
*! declared.
*!
*! Load order: as dwlsmv_fdwls, then dwlsmv_fdwls.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_minimize(real colvector s, real colvector w,
                            struct dwlsmv_parmap  scalar par,
                            struct dwlsmv_statmap vector maps,
                            real colvector start,
                            string scalar evaltype,
                            string scalar technique,
                            real scalar maxiter,
                            real colvector theta,
                            real scalar fmin,
                            real scalar iter)
{
    transmorphic S
    real scalar rc, conv, D, g

    if (rows(start) != par.npar) {
        errprintf("dwlsmv_minimize: start has %f row(s); the parameter table has %f free parameter(s)\n",
                  rows(start), par.npar)
        _error("dwlsmv_minimize: start is not npar rows long")
    }
    D = 0
    for (g = 1; g <= length(maps); g++) {
        D = D + maps[g].D
    }
    if (rows(s) != D | rows(w) != D) {
        errprintf("dwlsmv_minimize: s has %f row(s) and w %f; the index map(s) have D = %f\n",
                  rows(s), rows(w), D)
        _error("dwlsmv_minimize: s or w is not D rows long")
    }

    S = optimize_init()
    optimize_init_which(S, "min")
    optimize_init_evaluator(S, &dwlsmv_fdwls())
    optimize_init_evaluatortype(S, evaltype)
    optimize_init_technique(S, technique)
    optimize_init_argument(S, 1, s)
    optimize_init_argument(S, 2, w)
    optimize_init_argument(S, 3, par)
    optimize_init_argument(S, 4, maps)
    optimize_init_params(S, start')
    optimize_init_conv_maxiter(S, maxiter)
    optimize_init_conv_ptol(S, 1e-10)
    optimize_init_conv_vtol(S, 1e-12)
    optimize_init_conv_nrtol(S, 1e-12)
    optimize_init_conv_warning(S, "off")
    optimize_init_tracelevel(S, "none")

    rc = _optimize(S)

    theta = optimize_result_params(S)'
    fmin  = optimize_result_value(S)
    iter  = optimize_result_iterations(S)
    conv  = optimize_result_converged(S)

    if (rc != 0) {
        conv = 0
    }

    return(conv)
}

end
