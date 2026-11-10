*! dwlsmv -- Module 7: the DWLS fit
*!
*! Executed under StataMP 18 on 2026-09-23. theta-hat agrees with lavaan within 1.2e-3 of spec 7.2's tier on five fits: tests/unit/test_fit.do.
*!
*! The locked Section 5.1 entry point:
*!
*!     struct dwlsmv_fit scalar dwls_fit(struct dwlsmv_stage1 vector sts,
*!                                       struct dwlsmv_parmap scalar par)
*!
*! sts: one Stage 1 per group (spec 4.8 O14), a single group a 1-vector; the
*! fit's stage1 field is the same vector. Under groups (built 2026-10-10) the
*! objective is the pooled F = sum_g pi_g F_g, pi_g = N_g / N: s stacked group
*! by group and w = pi_g / diag(Gamma_g), weighting (a) (@sec-multigroup-scaling
*! MG2; spec Module 12 check 10). T_raw = (N - G) F-hat with N the total (MG3,
*! decided 2026-10-09), D_stat = sum_g D_g - npar. dwlsmv_fit_stack() builds the
*! stacked inputs; at G = 1 it returns s, W_D and N exactly as before.
*!
*! Spec Section 5, Module 7. Minimizes
*!
*!     F_DWLS(theta) = (s - sigma(theta))' W_D (s - sigma(theta)),
*!     W_D = diag(Gamma)^-1,
*!
*! starting from the ULS solution (W = I), which is itself a real nonlinear
*! optimization and not a closed form. Both runs are dwlsmv_minimize(), from
*! the table's own start values for ULS and from the ULS point for DWLS. The
*! DWLS run and everything after it are dwlsmv_fit_from(st, par, start),
*! defined here; D11's recovery calls it directly for its final fit, from the
*! transformed point with thresholds at Stage 1's sample thresholds, with no
*! ULS run (decided 2026-10-07, spec D11 and Module 7).
*!
*! W_D IS diag(Gamma)^-1 AND NOT diag(Gamma). Spec Section 4.9 records the trap:
*! Muthen, du Toit & Spisic (1997) write the fit function with W^-1, so their W
*! is this spec's W_D inverted, and the substitution fits a different estimator
*! without raising anything. It is formed exactly once, below, as a vector.
*! tests/unit/test_fit.do checks the convention against lavaan's unscaled
*! chi-square, which is proportional to F-hat: under the inverted weight both
*! theta-hat and F-hat move, and no constant reconciles them.
*!
*! WHAT IS FILLED, and what is deliberately left for Modules 8-9:
*!
*!     stage1, par    copies of the inputs, so the fit carries its own context
*!     theta          theta-hat, npar x 1, in Section 4.8's order, with factor
*!                    signs canonicalized by dwlsmv_signcanon() after
*!                    convergence: every reflectable factor's first free
*!                    loading is non-negative. The optimizer itself runs
*!                    unconstrained -- no sign bound, no artificial boundary
*!     Delta          dwls_delta() AT theta-hat -- Module 8's bread and Module
*!                    9's U both need it there, and nowhere else
*!     fmin           F_DWLS(theta-hat), exactly as defined above: no 1/2, no N
*!     converged      the DWLS run's; see below
*!     iterations     the DWLS run's
*!     T_raw          (N - G) * fmin; at G = 1 (N - 1) * fmin -- LAVAAN'S
*!                    single-group convention, decided
*!                    2026-09-23 (spec Section 4.9, derivation Gate 5). lavaan
*!                    0.7.2 sets fx = 1/2 (N-1)/N F in lav_model_objective()
*!                    and reports 2N fx, i.e. (N-1) F; Mplus uses N F. Asymptotically identical, and
*!                    1/N apart in finite samples -- 2e-4 at N = 5000, which
*!                    was inside the test-statistic tier until 2026-09-24
*!                    and invisible to it. test_fit.do check 5 tells them
*!                    apart on T, test_inference.do check 4 on T3.
*!     D_stat         sum_g D_g - npar (D - npar at G = 1)
*!
*!     acov, T3, T3_df, shift_a, shift_b     Modules 8-9: dwls_inference().
*!                                           Left at the struct's initial
*!                                           values here, never guessed.
*!
*! CONVERGENCE IS REPORTED, NEVER ENFORCED. 'converged' is 1 only if the DWLS
*! minimization declared convergence; theta then holds its optimum. If it did
*! not, theta holds the last iterate and converged is 0 -- non-convergence
*! itself never aborts; what an unconverged fit means is the caller's decision,
*! exactly as Module 2 leaves a non-converging pair to its caller. What DOES
*! abort is an evaluator that cannot evaluate at all: a singular (I - B) at a
*! trial point of a non-recursive model, or a sigma entry no parameter row
*! fills. dwlsmv_minimize()'s header lists both; neither is a converged = 0. The ULS run's outcome is not
*! carried: it only supplies a starting point, and a DWLS optimum reached from
*! an unconverged ULS iterate is still a DWLS optimum.
*!
*! ADMISSIBILITY IS NOT CHECKED HERE, deliberately. Spec Section 5 requires it
*! before a fitted theta-hat is accepted; every condition is a function of
*! (theta, par, map), so it is dwlsmv_admissible(fit.theta, fit.par,
*! fit.stage1.map) and not a field on this struct -- decided 2026-09-23, the
*! same rule as dwlsmv_nfac(). The caller asks, and dwlsmv_admissadvice()
*! turns the answer into notes. This routine prints nothing.
*!
*! THE OPTIMIZER is optimize()'s "nr" over a "d2" evaluator whose Hessian is
*! the Gauss-Newton 2 Delta' W Delta -- a damped Gauss-Newton. Chosen by
*! measurement over "d1" with "nr" and with "bfgs": tools/time_fit.do.
*!
*! Load order: everything dwls_model() and dwls_delta() need, then
*! dwlsmv_mg_check, dwlsmv_fdwls, dwlsmv_minimize, dwlsmv_signcanon. dwlsmv_fit_w,
*! dwlsmv_fit_stack and
*! dwlsmv_fit_from are defined above dwls_fit(), in this file.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// W_D = diag(Gamma)^-1 as a vector, after checking Stage 1's output. Formed
// here and nowhere else.
real colvector dwlsmv_fit_w(struct dwlsmv_stage1 scalar st)
{
    real colvector g

    g = diagonal(st.Gamma)

    if (rows(st.s) != st.map.D | rows(g) != st.map.D) {
        errprintf("dwls_fit: s has %f row(s) and Gamma %f; the index map has D = %f\n",
                  rows(st.s), rows(g), st.map.D)
        _error("dwls_fit: Stage 1 output does not match its own map")
    }
    if (missing(g) | min(g) <= 0) {
        errprintf("dwls_fit: diag(Gamma) has a zero, negative or missing entry, so W_D = diag(Gamma)^-1 does not exist.\n")
        errprintf("Every sample statistic needs a positive asymptotic variance; check Stage 1's output.\n")
        _error("dwls_fit: diag(Gamma) is not strictly positive")
    }
    return(1 :/ g)
}

// The stacked inputs of the objective (spec 4.8 O14; @sec-multigroup-scaling
// MG2). G = 1: s, W_D and N exactly as before. G > 1: s stacked group by group,
// w = pi_g W_D,g with pi_g = N_g / N (weighting (a), decided 2026-10-10), N the
// total. s, w, N and maps are filled by reference.
void dwlsmv_fit_stack(struct dwlsmv_stage1 vector sts,
                      real colvector s, real colvector w, real scalar N,
                      struct dwlsmv_statmap vector maps)
{
    real colvector wg
    real scalar G, g

    G = length(sts)
    if (G == 1) {
        s    = sts[1].s
        w    = dwlsmv_fit_w(sts[1])
        N    = sts[1].N
        maps = sts[1].map
        return
    }
    N = 0
    for (g = 1; g <= G; g++) {
        if (missing(sts[g].N) | sts[g].N < 1) {
            errprintf("dwls_fit: group %f's Stage 1 result has N = %f\n", g, sts[g].N)
            _error("dwls_fit: a group's N is missing or below 1")
        }
        N = N + sts[g].N
    }
    s    = J(0, 1, .)
    w    = J(0, 1, .)
    maps = dwlsmv_statmap(G, 1)
    for (g = 1; g <= G; g++) {
        wg      = dwlsmv_fit_w(sts[g])
        s       = s \ sts[g].s
        w       = w \ ((sts[g].N / N) :* wg)
        maps[g] = sts[g].map
    }
}

// dwlsmv_polish(s, w, par, maps, theta): the post-convergence polish (spec
// Module 7; @sec-polish; adopted 2026-10-23, E-014).
// optimize()'s "nr" extends a good step while F does not increase, and near
// theta-hat that comparison can be an ulp-level tie, so where it stops depends
// on the path by up to ~1e-7 (E-015). From a converged
// DWLS point this takes pure Gauss-Newton steps theta + P^-1 Delta'W r, with no
// extension, no halving and NO comparison of F, which near theta-hat contract
// at rate c = rho(P^-1 R) to a point that, in exact arithmetic, does not depend
// on the path; at machine level a floor remains (@sec-polish).
// Relative = check 6's unit, max_j |d_j| / (1 + |theta_j|). It stops after a
// step <= 1e-11, or after 20 steps (constants chosen after measuring: polish
// P0). The adopted rule also stops once a step <= tau_f fails to halve; with
// tau_f = tau_p = 1e-11 that test can never fire first, so it is not coded. It is
// ABANDONED, theta left exactly as passed, when a step cannot be computed
// (missing), the first step exceeds 1e-6, or a step grows while above 1e-11.
// theta is filled by reference. Returns (status, steps, last relative step):
// status 0 stopped by the threshold, 2 at the cap (theta kept; 1 unused),
// -1/-2/-3 abandoned as above. An evaluator that aborts still aborts.
real rowvector dwlsmv_polish(real colvector s, real colvector w,
                             struct dwlsmv_parmap  scalar par,
                             struct dwlsmv_statmap vector maps,
                             real colvector theta)
{
    real colvector th, d
    real rowvector g
    real matrix H
    real scalar v, k, rel, prev

    th   = theta
    prev = .
    rel  = .
    for (k = 1; k <= 20; k++) {
        dwlsmv_fdwls(2, th', s, w, par, maps, v, g, H)
        d = cholsolve(H, -g')
        if (missing(d) | missing(v)) {
            return((-1, k - 1, prev))
        }
        rel = max(abs(d) :/ (1 :+ abs(th)))
        if (k == 1 & rel > 1e-6) {
            return((-2, 0, rel))
        }
        if (k >= 2 & rel > prev & rel > 1e-11) {
            return((-3, k - 1, rel))
        }
        th = th + d
        if (rel <= 1e-11) {
            theta = th
            return((0, k, rel))
        }
        prev = rel
    }
    theta = th
    return((2, 20, rel))
}

// dwlsmv_fit_from(sts, par, start): the DWLS run and everything after it,
// from a given start. dwls_fit() is ULS from the table's starts, then this
// from the ULS point. D11's recovery calls it directly for its final fit,
// with no ULS run (spec D11, Module 7; decided 2026-10-07), passing one
// Stage 1 result, which is a 1-vector here. Not a Section 5.1 entry point.
// start is npar x 1 in theta order and is not modified.
struct dwlsmv_fit scalar dwlsmv_fit_from(struct dwlsmv_stage1 vector sts,
                                         struct dwlsmv_parmap  scalar par,
                                         real colvector start)
{
    struct dwlsmv_fit scalar f
    struct dwlsmv_statmap vector maps
    real colvector s, w, th
    real rowvector gr
    real matrix H
    real scalar fm, it, conv, N, G

    G = length(sts)
    dwlsmv_fit_stack(sts, s, w, N, maps)
    if (rows(start) != par.npar | missing(start)) {
        _error("dwlsmv_fit_from: start is not a complete npar x 1 vector")
    }

    // DWLS, W = W_D (pooled: weighting (a))
    conv = dwlsmv_minimize(s, w, par, maps, start,
                           "d2", "nr", 300, th, fm, it)

    // The post-convergence polish, only after a converged run (spec Module 7).
    // converged and iterations stay optimize()'s; an abandoned polish leaves
    // th exactly as optimize() returned it.
    if (conv) {
        (void) dwlsmv_polish(s, w, par, maps, th)
    }

    // Canonical factor signs, AFTER convergence and BEFORE anything is derived
    // from theta-hat, so Delta, F and everything Modules 8-9 build agree with
    // the reported theta. F is invariant under the reflection; it is
    // re-evaluated at the returned theta anyway, so fmin is F at theta by
    // construction rather than by an argument about rounding. Under groups a
    // factor is reflected in every group at once (spec 4.8 O16 (a)).
    th = dwlsmv_signcanon(th, par)
    dwlsmv_fdwls(0, th', s, w, par, maps, fm, gr, H)

    f.stage1     = sts
    f.par        = par
    f.theta      = th
    f.Delta      = dwls_delta(th, par, maps)
    f.fmin       = fm
    f.converged  = conv
    f.iterations = it
    f.T_raw      = (N - G) * fm         // spec 4.9 / MG3: T = (N - G) F; G = 1 is (N - 1) F
    f.D_stat     = rows(s) - par.npar

    return(f)
}

// DWLS from first; only if that does not converge, once more from second
// (spec Module 7, 2026-11-01). A converged first run is returned untouched. If
// the second also fails the first run's result is returned, as before the retry.
struct dwlsmv_fit scalar dwlsmv_fit_retry(struct dwlsmv_stage1 vector sts,
                                          struct dwlsmv_parmap  scalar par,
                                          real colvector first,
                                          real colvector second)
{
    struct dwlsmv_fit scalar f, f2

    f = dwlsmv_fit_from(sts, par, first)
    if (f.converged) {
        return(f)
    }
    f2 = dwlsmv_fit_from(sts, par, second)
    if (f2.converged) {
        printf("{txt}note: the fit did not converge from the ULS starting values; it converged from\n")
        printf("{txt}      the model's initial starting values instead.\n")
        return(f2)
    }
    printf("{txt}note: the fit did not converge from the ULS starting values or, on a retry, from\n")
    printf("{txt}      the model's initial starting values. The first attempt is reported.\n")
    return(f)
}

// The checks, the table's starts and the ULS point (spec Section 5, Module 7).
void dwlsmv_fit_uls(struct dwlsmv_stage1 vector sts,
                    struct dwlsmv_parmap  scalar par,
                    real colvector start, real colvector thu)
{
    struct dwlsmv_statmap vector maps
    real colvector s, w
    real scalar i, fu, iu, N

    if (length(sts) < 1) {
        _error("dwls_fit: no Stage 1 result was passed")
    }
    if (length(sts) > 1) {
        dwlsmv_mg_check(par, length(sts), "dwls_fit")
    }
    dwlsmv_fit_stack(sts, s, w, N, maps)

    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            start[par.fidx[i]] = par.value[i]
        }
    }
    if (missing(start)) {
        _error("dwls_fit: fidx does not number the free rows 1..npar")
    }

    // ULS, W = I: starting values only (spec Section 5, Module 7)
    (void) dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, start,
                           "d2", "nr", 300, thu, fu, iu)
}

// ULS, then DWLS, with no retry: the marker metric's first fit and D11's own
// unit-variance fit, so that D11 reads and runs as specified (user, 2026-11-01).
struct dwlsmv_fit scalar dwlsmv_fit_noretry(struct dwlsmv_stage1 vector sts,
                                            struct dwlsmv_parmap  scalar par)
{
    real colvector start, thu

    dwlsmv_fit_uls(sts, par, start, thu)
    return(dwlsmv_fit_from(sts, par, thu))
}

struct dwlsmv_fit scalar dwls_fit(struct dwlsmv_stage1 vector sts,
                                  struct dwlsmv_parmap  scalar par)
{
    real colvector start, thu

    dwlsmv_fit_uls(sts, par, start, thu)

    // DWLS, W = W_D (pooled: weighting (a)), from the ULS point; Module 3's
    // starts only if that fails
    return(dwlsmv_fit_retry(sts, par, thu, start))
}

end
