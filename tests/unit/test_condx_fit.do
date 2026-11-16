*! Conditional mode (one 0/1 covariate, one group): dwls_fit() and dwls_inference().
*!
*! Executed 2026-11-13 by the main session (Stata 18 MP). The condx_A and condx_D comparison (checks 1 to 8
*! below) was run once and then removed from this file: it is link 2 of the frozen preregistration.
*! fixes what breaks, and only then removes this banner.
*!
*! Needs fit_changes.md applied (dwlsmv_cx1_fit_check.mata, dwls_fit.mata, dwls_inference.mata,
*! dwlsmv_signcanon.mata). The tolerances below were fixed BEFORE any run and are not to be widened
*! to make a failure pass; a tolerance chosen by the agent without a measurement is marked (agent).
*!
*! REFERENCE. The committed tier-1 (base-R) reference of the Hessian-form fit, fixtures/generate/
*! 99n_condx_ref.R -> fixtures/MANIFEST_condx_ref.json, read with import delimited, stringcols(_all),
*! delimiters(",") and strtoreal():
*!     fixtures/data/condx_A.csv, condx_D.csv                          (y1 y2 y3 x; N = 1500; K = 3,3,2)
*!     fixtures/expected/condx_ref_A_theta.csv, condx_ref_D_theta.csv       (label, est, se)
*!     fixtures/expected/condx_ref_A_scalars.csv, condx_ref_D_scalars.csv   (quantity, value: T_raw,
*!                                                                           a, b, T3, df, ...)
*! The reference is the DWLS fit with W_D = diag(Gamma-hat)^-1, Gamma-hat in the decided Hessian
*! form, T = (N - 1) F, SE = sqrt(ACOV / (N - 1)); theta order tau, lambda, gamma is the table's.
*! The tiers are the preregistered link-2 table (docs/notes/condx-validation/prereg.md section 3,
*! set before any result): |a - b| <= e_abs + e_rel |b| with
*!     theta-hat  1e-5, 1e-4     SE  1e-8, 1e-5     T3  1e-6, 1e-5     a  none, 1e-5
*!     b  1e-6, none             df  exact           T_raw: the T3 tier (agent, not in the prereg)
*! THIS IS A DEVELOPMENT TEST, not link 2: link 2 is tests/unit/test_condx_link2.do with its own
*! hash (the prereg addendum). It does NOT use the lavaan conditional fit (lavaan_condx_fit_*.csv):
*! lavaan's Gamma is BHHH, and the decided Hessian form differs from it by up to 3.4% in an SE
*! (@sec-condx-full C2), so those files cannot be an equality target here.
*!
*! Checks, per dataset (condx_A, condx_D):
*!   1  converged; npar 9; D_stat == T3_df == df of the scalars file (2), exactly.
*!   2  theta-hat against the reference, every entry, with the first five labels matching the table's.
*!   3  SE = sqrt(diag(acov) / (N - 1)) against the reference.
*!   4  T3, a, b and T_raw against the reference.
*!   5  tau-hat equals the sample thresholds s_tau: |diff| <= 1e-12 (@sec-marker-recovery Result 5).
*!   6  stationarity: the Gauss-Newton step (Delta' W Delta)^-1 Delta' W r at theta-hat has
*!      max |d_j| / (1 + |theta_j|) <= 1e-9 (the polish stops at 1e-11; agent).
*!   7  independent dense algebra on the fit's own Delta: ACOV and (a, b) recomputed with explicit
*!      matrices (invsym, Wm = diag(w), U Gamma) against f.acov (relative to max |acov|, 1e-9) and
*!      f.shift_a, f.shift_b (1e-8 relative / absolute) (agent); T3 = a T_raw + b to 1e-14 relative.
*!   8  Delta against central differences of dwls_model() at theta-hat (step 1e-6 max(1,|theta|)):
*!      max |diff| <= 1e-8, and the ACOV built from the finite-difference Delta equals f.acov to
*!      1e-6 relative to max |acov| (agent).
*! Once:
*!   9  exact-fit populations (no data): s = sigma(theta*) and a symmetric positive definite Gamma
*!      constructed here. (i) 3 items, gamma; (ii) 4 items, gamma and one kappa. The fit must
*!      converge, return theta* to 1e-8, F <= 1e-16, sigma(theta-hat) - s <= 1e-9, D_stat = D - npar,
*!      and an ACOV that is symmetric with a positive minimum eigenvalue (agent tolerances).
*!  10  sign canonicalization with gamma: on a hand-built two-factor table with beta, gamma on F1 and
*!      kappa, reversing F1, F2 or both gives the expected theta EXACTLY and sigma(theta) unchanged
*!      to 1e-15; with gamma fixed nonzero F1 is not reflectable (theta returned unchanged).
*!  11  refusals, each rc == 3498 exactly, with the message read back from a log: G = 2 with a
*!      conditional second group (the stack call dwls_fit() makes); q = 2 in dwls_fit and in
*!      dwls_inference; a table with gamma and kappa rows against an unconditional Stage 1; the
*!      dwlsmv_sens* style call (no flag) on one conditional group.
*!
*! WHAT THIS CERTIFIES AND DOES NOT. It certifies that dwls_fit() and dwls_inference() reproduce a
*! tier-1 base-R Hessian-form fit on two datasets (K = 3,3,2, one factor, gamma only), to the link-2
*! tiers, and that the machinery is internally consistent on two more models. It does not certify
*! kappa, beta or a second factor against an external reference (no committed fixture has them; see
*! the report: generate one from fixtures/generate/99n_condx_ref.R), nor several covariates,
*! groups, admissibility, fit indices, DIFFTEST, or the model at |rho| near 1. A pass at the fit
*! layer is no evidence for Stage 1 or for the guarded layers.
*!
*! MUTANTS for the main session (edit, run, revert; each must turn this file red):
*!   M1  dwlsmv_fit_stack: ignore cxok (always dwlsmv_s2_guard_st): checks 1 to 9 abort at rc 3498.
*!   M2  dwlsmv_cx1_fit_check: `q > 1` -> `q > 0`: same as M1. `q > 1` -> `q > 2`: refusal 11b fails.
*!   M3  dwlsmv_cx1_fit_check: return(0) before the dwlsmv_s2_guard_st call (several groups): 11a.
*!   M4  dwls_inference: put dwlsmv_s2_guard_st back: inference checks abort at rc 3498.
*!   M5  dwlsmv_signcanon: delete the t == 8 branch of the apply loop: check 10 (sigma changes).
*!       Delete (t == 8 & r == f) from the condition: the fixed-gamma case of check 10 fails.
*!   M6  dwlsmv_fit_w: weight 1 on the slope rows (or only thresholds and correlations): 2, 3, 4.
*!   M7  dwls_fit: T_raw = N * fm: check 4 (T_raw, T3). D_stat = rows(s) - npar - map.p: 1.
*!   M8  dwlsmv_fit_uls: delete the dwlsmv_cx1_check call: 11c loses its "dwls_fit:" message.
*!   M9  dwls_inference: pass Gamma only on the threshold and correlation rows (zero the slope
*!       block of Omega): checks 3, 4, 7.
*!   M10 dwlsmv_fit_from: delete the dwlsmv_polish call: 5 and 6 may still pass at optimize()'s
*!       tolerances; a survivor here is a coverage gap, not a contract violation (mutant-survivor rule).
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

* Load every source, types first (as test_s2_guard.do does).
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_condx_fit: `f' failed to load"
        exit 9
    }
}

capture mata: (void) dwlsmv_cx1_fit_check(dwlsmv_stage1(0), "load check")
if _rc {
    display as error "test_condx_fit: dwlsmv_cx1_fit_check() is not defined (src/mata/dwlsmv_cx1_fit_check.mata missing?)"
    exit 9
}

* ---- no fixtures are read (2026-11-13) --------------------------------------------
* The comparison of production with the tier-1 reference on condx_A and condx_D is LINK 2 of the frozen
* preregistration docs/notes/condx-validation/prereg.md. It is not a unit test and must not run in the
* suite. The first draft of this file ran it once, before link 2's do-file was hashed
* (docs/notes/condx-validation/addendum-2026-11-13-early-run.md); the two blocks that did so are removed.
* cf_one() stays defined below, uncalled, as the record of what was compared.

mata:

// ---- helpers -------------------------------------------------------------------

real scalar cf_fail(string scalar what)
{
    printf("{err}test_condx_fit: %s\n", what)
    return(1)
}

real scalar cf_chk(real scalar ok, string scalar what)
{
    if (!ok) {
        return(cf_fail(what))
    }
    return(0)
}

// |a - b| <= eabs + erel |b|
real scalar cf_tier(real scalar a, real scalar b, real scalar eabs, real scalar erel)
{
    if (missing(a) | missing(b)) {
        return(0)
    }
    return(abs(a - b) <= eabs + erel * abs(b))
}

real scalar cf_look(string colvector QN, real colvector QV, string scalar nm)
{
    real scalar i

    for (i = 1; i <= rows(QN); i++) {
        if (QN[i] == nm) {
            return(QV[i])
        }
    }
    return(.)
}

string scalar cf_join(string rowvector a)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(a); i++) {
        s = s + a[i] + char(10)
    }
    return(s)
}

// condx_A / condx_D's model: y1..y3 on F, unit disturbance variance, f ~ x.
string scalar cf_spec3()
{
    return(cf_join(("factor F", "ind y1 y2 y3",
                    "lambda F y1 free", "lambda F y2 free", "lambda F y3 free",
                    "metric unitvariance F", "covariate x", "gamma F x free")))
}

// four items, K = (3,3,2,2): f ~ x and a direct effect of x on y3 (y4 is the anchor).
string scalar cf_spec4()
{
    return(cf_join(("factor F", "ind y1 y2 y3 y4",
                    "lambda F y1 free", "lambda F y2 free", "lambda F y3 free", "lambda F y4 free",
                    "metric unitvariance F", "covariate x", "gamma F x free", "kappa y3 x free")))
}

// A symmetric positive definite Gamma, exactly symmetric by construction (min eigenvalue >= 0.5).
real matrix cf_gamma(real scalar D)
{
    real colvector v
    real scalar i

    v = J(D, 1, 1)
    for (i = 1; i <= D; i++) {
        v[i] = (-1)^i
    }
    return(diag(0.5 :+ (1::D) :/ D) + 0.05 * (v * v'))
}

// An unconditional and a conditional Stage 1 result with placeholder statistics (refusal cases only).
struct dwlsmv_stage1 scalar cf_st(real colvector K, real scalar q)
{
    struct dwlsmv_stage1 scalar st

    st.map    = dwlsmv_statmap_build(K, q)
    st.s      = J(st.map.D, 1, 0.1)
    st.Gamma  = I(st.map.D)
    st.N      = 500
    st.psibar = J(st.map.D, 1, 0)
    st.bound  = J(st.map.P, 1, 0)
    return(st)
}

struct dwlsmv_statmap vector cf_maps(struct dwlsmv_statmap scalar m)
{
    struct dwlsmv_statmap vector mv

    mv = dwlsmv_statmap(1, 1)
    mv[1] = m
    return(mv)
}

// ---- checks 1 to 8: one dataset against the tier-1 reference ------------------------

real scalar cf_one(real matrix Y, real colvector X, string colvector TL, real colvector TE,
                   real colvector TS, string colvector QN, real colvector QV, string scalar tag)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f
    struct dwlsmv_statmap vector mv
    real colvector w, r, d, se, tp, tm
    real matrix Dfd, E, Ei, U, M, V, Vfd, Wm
    real scalar nf, i, k, h, N, df, trM, trM2, a, b, T, D, scl

    nf = 0
    st = dwls_sample_stats_cx(Y, (3 \ 3 \ 2), X)
    mv = cf_maps(st.map)
    par = dwls_identify(cf_spec3(), st.map)
    sv = dwlsmv_stage1(1, 1)
    sv[1] = st
    f = dwls_fit(sv, par)
    dwls_inference(f)

    N  = st.N
    D  = st.map.D
    df = cf_look(QN, QV, "df")

    // 1
    nf = nf + cf_chk(f.converged == 1, tag + ": the fit did not converge")
    nf = nf + cf_chk(par.npar == 9 & rows(f.theta) == 9, tag + ": npar is not 9")
    nf = nf + cf_chk(f.D_stat == df & f.T3_df == df & D - par.npar == df,
                     tag + ": D_stat is " + strofreal(f.D_stat) + ", wanted " + strofreal(df) + " (D = " +
                     strofreal(D) + ", npar " + strofreal(par.npar) + ")")
    if (rows(f.theta) != 9 | rows(f.acov) != 9) {
        printf("{err}test_condx_fit: %s: theta or acov has the wrong size; later checks skipped\n", tag)
        return(nf + 1)
    }

    // 2
    for (i = 1; i <= 5; i++) {
        nf = nf + cf_chk(par.name[i] == TL[i], tag + ": table label " + par.name[i] + " vs reference " + TL[i])
    }
    for (i = 1; i <= 9; i++) {
        if (!cf_tier(f.theta[i], TE[i], 1e-5, 1e-4)) {
            nf = nf + cf_fail(tag + ": theta[" + strofreal(i) + "] = " + strofreal(f.theta[i], "%21.15g") +
                              ", reference " + strofreal(TE[i], "%21.15g"))
        }
    }

    // 3
    se = sqrt(diagonal(f.acov) :/ (N - 1))
    for (i = 1; i <= 9; i++) {
        if (!cf_tier(se[i], TS[i], 1e-8, 1e-5)) {
            nf = nf + cf_fail(tag + ": SE[" + strofreal(i) + "] = " + strofreal(se[i], "%21.15g") +
                              ", reference " + strofreal(TS[i], "%21.15g"))
        }
    }

    // 4
    nf = nf + cf_chk(cf_tier(f.T3, cf_look(QN, QV, "T3"), 1e-6, 1e-5),
                     tag + ": T3 " + strofreal(f.T3, "%21.15g") + " vs " + strofreal(cf_look(QN, QV, "T3"), "%21.15g"))
    nf = nf + cf_chk(cf_tier(f.shift_a, cf_look(QN, QV, "a"), 0, 1e-5),
                     tag + ": a " + strofreal(f.shift_a, "%21.15g") + " vs " + strofreal(cf_look(QN, QV, "a"), "%21.15g"))
    nf = nf + cf_chk(cf_tier(f.shift_b, cf_look(QN, QV, "b"), 1e-6, 0),
                     tag + ": b " + strofreal(f.shift_b, "%21.15g") + " vs " + strofreal(cf_look(QN, QV, "b"), "%21.15g"))
    nf = nf + cf_chk(cf_tier(f.T_raw, cf_look(QN, QV, "T_raw"), 1e-6, 1e-5),
                     tag + ": T_raw " + strofreal(f.T_raw, "%21.15g") + " vs " + strofreal(cf_look(QN, QV, "T_raw"), "%21.15g"))

    // 5: tau-hat = s_tau
    if (max(abs(f.theta[1::5] - st.s[1::5])) > 1e-12) {
        nf = nf + cf_fail(tag + ": tau-hat differs from the sample thresholds by " +
                          strofreal(max(abs(f.theta[1::5] - st.s[1::5]))))
    }

    // 6: stationarity, in the unit of the polish's stopping rule
    w  = 1 :/ diagonal(st.Gamma)
    r  = st.s - dwls_model(f.theta, par, mv)
    d  = cholsolve(quadcross(f.Delta, w, f.Delta), quadcross(f.Delta, w :* r))
    scl = max(abs(d) :/ (1 :+ abs(f.theta)))
    nf = nf + cf_chk(scl <= 1e-9, tag + ": the Gauss-Newton step at theta-hat is " + strofreal(scl) + " (relative)")
    nf = nf + cf_chk(reldif(f.fmin, quadcross(r, w, r)) <= 1e-12, tag + ": fmin is not F at theta-hat")
    nf = nf + cf_chk(reldif(f.T_raw, (N - 1) * f.fmin) <= 1e-14, tag + ": T_raw is not (N - 1) fmin")

    // 7: dense algebra on the fit's own Delta
    Wm  = diag(w)
    E   = f.Delta' * Wm * f.Delta
    Ei  = invsym(E)
    V   = Ei * f.Delta' * Wm * st.Gamma * Wm * f.Delta * Ei
    U   = Wm - Wm * f.Delta * Ei * f.Delta' * Wm
    M   = U * st.Gamma
    trM  = trace(M)
    trM2 = trace(M * M)
    a   = sqrt(df / trM2)
    b   = df - a * trM
    nf = nf + cf_chk(max(abs(V - f.acov)) <= 1e-9 * max(abs(V)),
                     tag + ": acov differs from the dense sandwich by " + strofreal(max(abs(V - f.acov))))
    nf = nf + cf_chk(reldif(a, f.shift_a) <= 1e-8 & abs(b - f.shift_b) <= 1e-8,
                     tag + ": (a, b) differ from the dense traces: " + strofreal(a, "%12.8g") + ", " +
                     strofreal(b, "%12.8g") + " vs " + strofreal(f.shift_a, "%12.8g") + ", " +
                     strofreal(f.shift_b, "%12.8g"))
    nf = nf + cf_chk(reldif(f.T3, f.shift_a * f.T_raw + f.shift_b) <= 1e-14, tag + ": T3 is not a T_raw + b")

    // 8: Delta against central differences of dwls_model()
    Dfd = J(D, 9, .)
    for (k = 1; k <= 9; k++) {
        h     = 1e-6 * max((1, abs(f.theta[k])))
        tp    = f.theta
        tm    = f.theta
        tp[k] = f.theta[k] + h
        tm[k] = f.theta[k] - h
        Dfd[., k] = (dwls_model(tp, par, mv) - dwls_model(tm, par, mv)) / (2 * h)
    }
    nf = nf + cf_chk(max(abs(Dfd - f.Delta)) <= 1e-8,
                     tag + ": Delta differs from finite differences by " + strofreal(max(abs(Dfd - f.Delta))))
    E   = Dfd' * Wm * Dfd
    Ei  = invsym(E)
    Vfd = Ei * Dfd' * Wm * st.Gamma * Wm * Dfd * Ei
    nf = nf + cf_chk(max(abs(Vfd - f.acov)) <= 1e-6 * max(abs(Vfd)),
                     tag + ": acov differs from the finite-difference sandwich by " + strofreal(max(abs(Vfd - f.acov))))

    printf("%s: converged %g, iterations %g, T3 %12.8f, a %10.6f, b %10.6f\n", tag, f.converged, f.iterations,
           f.T3, f.shift_a, f.shift_b)
    return(nf)
}

// ---- check 9: exact-fit populations -------------------------------------------------

real scalar cf_pop(real scalar which)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f
    real colvector ths, ev
    real colvector K
    real scalar nf, D, tag

    nf = 0
    if (which == 1) {
        K     = (3 \ 3 \ 2)
        par   = dwls_identify(cf_spec3(), dwlsmv_statmap_build(K, 1))
        ths   = (-0.5 \ 0.5 \ -0.3 \ 0.6 \ 0.2 \ 0.8 \ 0.7 \ 0.6 \ 0.4)
    }
    else {
        K     = (3 \ 3 \ 2 \ 2)
        par   = dwls_identify(cf_spec4(), dwlsmv_statmap_build(K, 1))
        ths   = (-0.5 \ 0.5 \ -0.3 \ 0.6 \ 0.2 \ -0.2 \ 0.8 \ 0.7 \ 0.6 \ 0.75 \ 0.4 \ -0.3)
    }
    if (par.npar != rows(ths)) {
        return(cf_fail("population " + strofreal(which) + ": npar is " + strofreal(par.npar) +
                       ", the constructed theta has " + strofreal(rows(ths))))
    }

    st.map    = dwlsmv_statmap_build(K, 1)
    D         = st.map.D
    st.s      = dwls_model(ths, par, cf_maps(st.map))
    st.Gamma  = cf_gamma(D)
    st.N      = 1000
    st.psibar = J(D, 1, 0)
    st.bound  = J(st.map.P, 1, 0)
    sv = dwlsmv_stage1(1, 1)
    sv[1] = st

    f = dwls_fit(sv, par)
    dwls_inference(f)

    tag = which
    nf = nf + cf_chk(f.converged == 1, "population " + strofreal(tag) + ": the fit did not converge")
    nf = nf + cf_chk(max(abs(f.theta - ths)) <= 1e-8,
                     "population " + strofreal(tag) + ": theta-hat misses theta* by " + strofreal(max(abs(f.theta - ths))))
    nf = nf + cf_chk(f.fmin <= 1e-16, "population " + strofreal(tag) + ": F is " + strofreal(f.fmin) + ", not 0")
    nf = nf + cf_chk(max(abs(dwls_model(f.theta, par, cf_maps(st.map)) - st.s)) <= 1e-9,
                     "population " + strofreal(tag) + ": sigma(theta-hat) does not reproduce s")
    nf = nf + cf_chk(f.D_stat == D - par.npar & f.T3_df == D - par.npar,
                     "population " + strofreal(tag) + ": D_stat is not D - npar")
    nf = nf + cf_chk(max(abs(f.acov - f.acov')) <= 1e-12 * max(abs(f.acov)),
                     "population " + strofreal(tag) + ": acov is not symmetric")
    ev = symeigenvalues(f.acov)'
    nf = nf + cf_chk(min(ev) > 0, "population " + strofreal(tag) + ": acov is not positive definite")
    return(nf)
}

// ---- check 10: sign canonicalization with gamma ---------------------------------------
// A hand-built table (after test_condx_model_delta.do's cx_par2): K = (3,2,2,3,2,2), F1 -> y1..y3,
// F2 -> y4..y6, beta F2 <- F1, gamma F1 <- x, kappa on y2 and y5, psi_11 fixed 1, psi_22 free.
// theta = (tau: 8 | lambda: 6 | beta | gamma | kappa_2, kappa_5 | psi_22). gfree 0 fixes gamma at 0.45.

struct dwlsmv_parmap scalar cf_par2(real scalar gfree)
{
    struct dwlsmv_parmap scalar par
    real colvector pt, rw, cl, fr, vl
    real scalar i, n, k

    pt = (4\4\4\4\4\4\4\4 \ 1\1\1\1\1\1 \ 3 \ 8 \ 9\9 \ 2 \ 2)
    rw = (1\1\2\3\4\4\5\6 \ 1\2\3\4\5\6 \ 2 \ 1 \ 2\5 \ 1 \ 2)
    cl = (1\2\1\1\1\2\1\1 \ 1\1\1\2\2\2 \ 1 \ 1 \ 1\1 \ 1 \ 2)
    fr = (J(15, 1, 1) \ gfree \ 1 \ 1 \ 0 \ 1)
    vl = (-0.5\0.6\-0.4\0.3\-0.2\0.7\0.1\-0.3 \ 0.8\0.7\0.6\0.75\0.65\0.55 \ 0.4 \ 0.45 \ -0.3\0.25 \ 1 \ 0.8)

    n = rows(pt)
    par.nall  = n
    par.npar  = sum(fr)
    par.name  = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        par.name[i] = "row" + strofreal(i)
    }
    par.ptype = pt
    par.row   = rw
    par.col   = cl
    par.free  = fr
    par.value = vl
    par.fidx  = J(n, 1, 0)
    k = 0
    for (i = 1; i <= n; i++) {
        if (fr[i]) {
            k = k + 1
            par.fidx[i] = k
        }
    }
    par.grp = J(n, 1, 1)
    return(par)
}

real scalar cf_sign()
{
    struct dwlsmv_parmap scalar par, parf
    struct dwlsmv_statmap vector mv
    real colvector th, can, want, ths, thf
    real scalar nf, c

    nf = 0
    mv = cf_maps(dwlsmv_statmap_build((3 \ 2 \ 2 \ 3 \ 2 \ 2), 1))
    par = cf_par2(1)
    // positive reference theta: lambda 6, beta, gamma, kappa x 2, psi_22 after the 8 thresholds
    ths = (-0.45\0.55\-0.35\0.25\-0.15\0.65\0.05\-0.25 \ 0.78\0.66\0.58\0.72\0.61\0.52 \ 0.35 \ 0.45 \ -0.25\0.2 \ 0.85)

    for (c = 1; c <= 3; c++) {
        // th is the input (some factor's loadings negative); want is the canonical theta
        th   = ths
        want = ths
        if (c == 1) {                      // F1's loadings negative: F1 is reversed, so its
            th[9::11]   = -ths[9::11]      // loadings turn positive and beta (one end at F1) and
            want[15]    = -ths[15]         // gamma (F1's) change sign
            want[16]    = -ths[16]
        }
        else if (c == 2) {                 // F2 reversed: loadings 4-6 negative
            th          = ths
            th[12::14]  = -ths[12::14]
            want        = ths
            want[15]    = -ths[15]         // beta F2 <- F1 has one end at F2
        }
        else {                             // both reversed: beta joins two reversed factors
            th          = ths
            th[9::14]   = -ths[9::14]
            want        = ths
            want[16]    = -ths[16]         // gamma of F1 only
        }
        can = dwlsmv_signcanon(th, par)
        nf = nf + cf_chk(max(abs(can - want)) == 0,
                         "check 10, case " + strofreal(c) + ": canonical theta is not the expected one (max diff " +
                         strofreal(max(abs(can - want))) + ")")
        nf = nf + cf_chk(max(abs(dwls_model(can, par, mv) - dwls_model(th, par, mv))) <= 1e-15,
                         "check 10, case " + strofreal(c) + ": sigma(theta) changed under canonicalization by " +
                         strofreal(max(abs(dwls_model(can, par, mv) - dwls_model(th, par, mv)))))
    }

    // a point that is already canonical is returned exactly
    nf = nf + cf_chk(max(abs(dwlsmv_signcanon(ths, par) - ths)) == 0, "check 10: a canonical theta moved")

    // gamma FIXED nonzero: F1 is not reflectable, so a negative first loading stays (a fixed-gamma
    // theta has 18 entries: gamma is the 16th table row and not in theta)
    parf = cf_par2(0)
    thf  = (ths[1::15] \ ths[17::19])
    thf[9::11] = -thf[9::11]
    can  = dwlsmv_signcanon(thf, parf)
    nf = nf + cf_chk(max(abs(can - thf)) == 0, "check 10: a fixed nonzero gamma did not block the reflection")
    return(nf)
}

// ---- check 11: refusals, one scenario per function; the caller reads rc and log -------

real scalar cf_r1()               // the stack call dwls_fit() makes, a conditional second group
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_statmap vector maps
    real colvector s, w
    real scalar N

    sv = dwlsmv_stage1(2, 1)
    sv[1] = cf_st((3 \ 3 \ 2), 0)
    sv[2] = cf_st((3 \ 3 \ 2), 1)
    dwlsmv_fit_stack(sv, s, w, N, maps, 1)
    return(0)
}

real scalar cf_r2()               // q = 2 through dwls_fit()
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f

    par = dwls_identify(cf_spec3(), dwlsmv_statmap_build((3 \ 3 \ 2), 1))
    sv = dwlsmv_stage1(1, 1)
    sv[1] = cf_st((3 \ 3 \ 2), 2)
    f = dwls_fit(sv, par)
    return(0)
}

real scalar cf_r3()               // q = 2 through a hand-assembled fit into dwls_inference()
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_fit scalar f

    sv = dwlsmv_stage1(1, 1)
    sv[1] = cf_st((3 \ 3 \ 2), 2)
    f.stage1 = sv
    dwls_inference(f)
    return(0)
}

real scalar cf_r4()               // G = 2 hand-assembled fit, second group conditional
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_fit scalar f

    sv = dwlsmv_stage1(2, 1)
    sv[1] = cf_st((3 \ 3 \ 2), 0)
    sv[2] = cf_st((3 \ 3 \ 2), 1)
    f.stage1 = sv
    dwls_inference(f)
    return(0)
}

real scalar cf_r5()               // gamma and kappa rows in the table, an unconditional Stage 1
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f

    par = dwls_identify(cf_spec3(), dwlsmv_statmap_build((3 \ 3 \ 2), 1))
    sv = dwlsmv_stage1(1, 1)
    sv[1] = cf_st((3 \ 3 \ 2), 0)
    f = dwls_fit(sv, par)
    return(0)
}

real scalar cf_r6()               // the diagnostics' call: no flag, one conditional group
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_statmap vector maps
    real colvector s, w
    real scalar N

    sv = dwlsmv_stage1(1, 1)
    sv[1] = cf_st((3 \ 3 \ 2), 1)
    dwlsmv_fit_stack(sv, s, w, N, maps)
    return(0)
}

// Does the log at f contain the text has? Returns the number of failures.
real scalar cf_msg(string scalar f, string scalar has, string scalar tag)
{
    string scalar txt

    txt = invtokens(cat(f)', char(10))
    if (!strpos(txt, has)) {
        printf("{err}test_condx_fit: case %s: the log lacks: %s\n", tag, has)
        return(1)
    }
    return(0)
}

end

local nfail 0

* ---- 9, 10: populations and sign canonicalization -------------------------------------
foreach c in "cf_pop(1)" "cf_pop(2)" "cf_sign()" {
    capture noisily mata: st_local("nf", strofreal(`c'))
    if _rc {
        display as error "test_condx_fit: `c' aborted, rc " _rc
        local nfail = `nfail' + 1
    }
    else {
        display as text "  `c': `nf' failure(s)"
        local nfail = `nfail' + `nf'
    }
}

* ---- 11: refusals, rc exactly 3498 and the message on the log ----------------------------
tempfile L
local n 0
local ++n
local c`n' cf_r1()
local w`n' dwls_fit: the index map of group 2 is conditional
local ++n
local c`n' cf_r2()
local w`n' dwls_fit: the index map has q = 2
local ++n
local c`n' cf_r3()
local w`n' dwls_inference: the index map has q = 2
local ++n
local c`n' cf_r4()
local w`n' dwls_inference: the index map of group 2 is conditional
local ++n
local c`n' cf_r5()
local w`n' dwls_fit: the table has gamma or kappa rows but the index map has no covariate
local ++n
local c`n' cf_r6()
local w`n' dwls_fit: the index map of group 1 is conditional

forvalues i = 1/`n' {
    quietly log using "`L'_`i'.log", text replace name(cf)
    capture noisily mata: (void) `c`i''
    local rc = _rc
    quietly log close cf
    display as text "  `c`i'': rc `rc'"
    if `rc' != 3498 {
        display as error "test_condx_fit: `c`i'' gave rc `rc', wanted 3498"
        local nfail = `nfail' + 1
    }
    else {
        local tmpL "`L'_`i'.log"
        local want "`w`i''"
        capture noisily mata: st_local("nm", strofreal(cf_msg(st_local("tmpL"), st_local("want"), "`i'")))
        if _rc {
            display as error "test_condx_fit: reading the log of `c`i'' failed, rc " _rc
            local nfail = `nfail' + 1
        }
        else local nfail = `nfail' + `nm'
    }
}

if `nfail' {
    display as error _n "test_condx_fit: `nfail' failure(s)"
    exit 9
}

display as result "test_condx_fit: dwls_fit and dwls_inference on a conditional Stage 1 result are self-consistent on two exact-fit populations, canonicalize gamma with the factor, and refuse what is not specified at rc 3498"
exit 0
