*! test_sens_t3.do -- dwlsmv_sens_t3(): B_T3 of @sec-mg-sensitivity-t3 (TS1-TS7)
*!
*! Executed under StataMP 18 on 2026-10-16; green at six processors.
*! Executed under StataMP 16.1 on 2026-10-19 at 1, 2 and 4 processors against
*! the per-count references (results in E-047, 2026-10-19).
*!
*! EXACT PATH ADOPTED (user, plan step 5, 2026-10-19; E-047): r1x and
*! r2x run every case through dwlsmv_sens_t3_x / dwlsmv_sens_x (Result MS6's
*! closed-form dDelta/dtheta and J) against tests/unit/sens_t3_ref_s<v>_p<p>.mata,
*! now exported on the exact path (tools/export_sens_t3_in.do with ST3PATH exact).
*! r1 and r2 stay on the finite-difference path, against
*! tests/unit/sens_t3_fdref_s<v>_p<p>.mata (the earlier per-count references,
*! re-exported; their values are identical to the files they replace), so
*! dwlsmv_sens_t3's default path, dwlsmv_sens and the helpers keep their pins,
*! as do h1-h3 and a1-a5, unchanged. Both files are required for this version
*! and count, else the test skips (exit 8). Per-count references and 1e-12 are
*! kept. tests/unit/sens_t3_ref.mata (StataMP 18, six processors) is
*! finite-difference only and is no longer read: under Stata 18 the test skips
*! until an exact-path and FD pair is exported on the six-processor machine. It
*! stays in the repository as the provenance of the 2026-10-16 build.
*! Per-count green is not cross-count equivalence: green at 1, 2 and 4
*! processors means each count agrees with its OWN reference to 1e-12, not that
*! the counts agree with one another. The measured cross-count gaps on the exact
*! path stand (step 4, X7: case 2 at 2.7e-11 in beta_2 and 8.9e-12 in B_T3;
*! metric_bigeps 1.03e-12), as does the open T_g sensitivity at gbar ~ 0.
*!
*! Mutants (scratch worktree, each diff printed before its run,
*! %TEMP%\st3ref\mutate.sh, mutate2.sh and their logs): 31 of 32 red.
*!   M1  G_Delta without kappa G2 (r1)     M17 kappa without the 2 (r1)
*!   M2  G2 with -2, not -4 (r1)           M18 beta_elem 2u on Delta (r1)
*!   M3  G_Omega without the T part (r1)   M19 |T3| < 1e-2 guard removed (e4)
*!   M4  G1_Omega without dW (r1)          M20 partial-inputs guard removed (e7)
*!   M5  beta_1 not doubled (r1)           M21 dOm missing guard removed (e9)
*!   M6  Theta without S_C (r1)            M22 D_stat guard removed (e3)
*!   M7  beta_2 without Theta'|J|Theta     M23 inference guard on T3 only (e2)
*!   M8  beta_2 without |v| (r1)           M25 rhoT guard removed (e6)
*!   M9  ||eps||_inf forced to 0 (r1.4)    M26 beta_eval without X3's gq term
*!   M10 beta_3 without |g_s|'|ds| (r1)    M27 beta_eval without G_P (r1)
*!   M11 c_ch = 0 (r1)                     M28 Z2 without U Omega (r1)
*!   M12 beta_eval without 4u(...) (r1)    M29 Ebar without h_j in u(.) (r1, h1)
*!   M13 G_Gm without its C part (r1)      M30 J not symmetrized (h2)
*!   M14 beta_elem without sigma (r1)      M31 dDelta at h = 1e-5 (r1, h3.a)
*!   M15 rhoT on every term (r1, a1, a4)   M32 beta_1 from G1, not G_Delta (r1)
*!   M16 B_eval without r_T (r1, a4)
*! M9 survived the first three cases (eps 3e-16 is too small for its term to
*! reach 1e-12); case metric_bigeps was added for it, and it is then red.
*! M14 did not apply at first (pattern error) and was re-run. Survivor, judged
*! equivalent: M24, the epsV guard removed. dwlsmv_sens() is called first and
*! refuses a negative epsV at rc 3498 itself (e5 stays green).
*! Mutants after the exact path was adopted (2026-10-19, StataMP 16.1 at two
*! processors, scratch copy %TEMP%\t3mut, each diff printed before its run,
*! t3mut.log). None of the earlier red ones turned green: M5, M6, M7, M8 (r1 and
*! r1x), M29 (h1, r1, r1x), M30 (h2), M31 (h3.a, r1). New, for the exact
*! wiring: N1 default routed exact (r1), N2 t3_x's dD by FD (r1x), N4 t3_x's
*! Theta by FD (r1x), N6 jmat_x's R from r (r1x), N7 sens_x's R scaled (r1x,
*! r2x): red. N3 (t3_x's J by FD) and N5 (a wrong Hessian entry, {nu, omega})
*! survive here; test_sens_hess.do catches both (tj; x1). The references hold
*! Stata's own dD and J, so a Hessian entry is pinned in test_sens_hess.
*! Q4 ADOPTED in beta_eval (user, 2026-10-20; E-087): every
*! quadcross() output's eta is SE6's 2 gamma_(n-1) X, not 2u X. Both
*! StataMP 16.1 sets (exact and FD, 1, 2 and 4 processors) were re-exported
*! with tools/export_sens_t3_ref.R calling t3_eta_q4; only columns 5, 7 and 8
*! moved. Green at 1, 2 and 4 processors. Mutants on that run (scratch copy
*! %TEMP%\q4exp\mut, four processors, each diff printed before its run):
*! M33 hD = hq = 2u (Q1 restored): red (r1.1-4, r1x.1-4), which pins the Q4
*! form. M11, M12, M26 and M27 re-run because beta_eval grew: still red (same
*! checks). M5-M8 and M29-M31, N1-N7 do not touch beta_eval; not re-run.
*! Logs: E-087 (t_p1, t_p2, t_p4, m_M*).
*! Division of responsibility (clarified 2026-10-19): the exact references hold
*! Stata's own dDelta/dtheta and J, so r1x and r2x establish that the assembled
*! exact-path results reproduce those references; they do not independently
*! establish every Hessian entry. That is test_sens_hess.do's job (references
*! from Result MS6's closed form; its x1 and tj catch N5 and N3). N3 and N5
*! surviving here is that split, not a coverage gap.
*!
*! Fits: design E (fixtures/data/mg_E.csv, G = 2) at the metric and scalar
*! levels, from tests/unit/sens_t3_setup.mata. E is consistency-only data for
*! the check 6 investigation; nothing here measures a processor-count residue,
*! and no magnitude of B_T3 is asserted. epsV = 2u, rhoT = 0.01.
*!
*!   r1  every column, four cases (metric and scalar with the setup's
*!       synthetic ds, eps, dOm; metric with none; metric with eps at 1e-6,
*!       added after mutant M9 survived), against the subsection's own
*!       R functions on the same inputs (written by tools/export_sens_t3_in.do
*!       then tools/export_sens_t3_ref.R). Column 1 also shows the fit is the
*!       one the reference was computed on.
*!       PER VERSION AND PROCESSOR COUNT (policy, user, 2026-10-19): the inputs
*!       are one run's fit, and the scalar case's theta-hat moves with the count
*!       (E-047, 2026-10-19; the dDelta/dtheta diagnostic). The
*!       reference is tests/unit/sens_t3_ref_s<version>_p<processors>.mata,
*!       version = floor(c(stata_version)); StataMP 18 at six processors keeps
*!       tests/unit/sens_t3_ref.mata. With no reference for this version and
*!       count the test skips (exit 8). The tolerance is unchanged
*!   r2  Theta (Stata, from dwlsmv_sens()'s columns as the routine forms it)
*!       against R's t3_Theta() on Stata's P, A, c, r, step and K
*!   h1  dwlsmv_sens_ebar() is dwlsmv_sens()'s Eb: column 6 recomputed from it
*!       equals dwlsmv_sens()'s, bitwise
*!   h2  dwlsmv_sens_jmat() is dwlsmv_sens()'s J: cholsolve(J, (Delta :* w)')
*!       equals dwlsmv_sens()'s A, bitwise
*!   h3  dwlsmv_sens_ddelta(): block 1 is the central difference of
*!       dwls_delta() at h = 1e-4 max(1, |theta_1|), bitwise; and every block
*!       agrees with the same difference at 10 h (Lemma MS2.1: no truncation
*!       error, only Delta_c's noise)
*!   a1  rhoT enters column 7 only, linearly: columns 2-6 and 8 bitwise equal
*!       at rhoT = 0 and 0.01; B(0.01) - B(0) = 0.01 (c2 + c3 + c4)
*!   a2  no measured inputs: beta_3 (column 4) is exactly 0
*!   a3  beta_3 is linear in (ds, dOm) with eps = 0: doubling both doubles it
*!   a4  column 8 = column 5 + r_T and column 7 = TS7's formula, bitwise
*!   a5  epsV: columns 5 (beta_eval) bitwise independent of it; columns 2 and 6
*!       affine in it; column 3 at most quadratic (third difference over 0, u,
*!       2u, 3u)
*!   e1-e9 refusals, each rc 3498: not converged; inference not run;
*!       D_stat = 0; |T3| < 1e-2; epsV negative; rhoT negative; only ds and
*!       eps given; ds short; dOm with a missing entry
*!
*! Tolerances, chosen before the first run: 1e-12 relative on r1 and r2 (the
*! same inputs, R's crossprod against Mata's quadcross); bitwise on h1, h2,
*! h3's first part, a1's equalities, a2 and a4; 1e-6 of max |dD| on h3's
*! second part (O(|E| / h) noise, ~1e-7 estimated in TS4); 1e-12 relative on
*! a1's difference and a3; 1e-10 of the largest value on a5's differences.
*!
*! CHANGED AFTER THE FIRST RUN (green): r1's measure was normwise over the
*! row, so column 1 (T3 ~ 16) hid columns 2-8 (~1e-11), and the first run's
*! r1 figures (~3e-26) said nothing about them. r1 is now per column,
*! relative, and exact where R gives 0; the tolerance is unchanged. r2 stays
*! normwise over Theta, as written: elementwise, one threshold element with
*! Theta ~ 7e-22 (1e11 below the largest) differs by 6e-6 relative in the
*! no-input case, an absolute 4e-27 (tools/export_sens_t3_ref.R's log line).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index
do tests/unit/sens_unpolished.mata
do tests/unit/sens_t3_setup.mata
local st3v = floor(c(stata_version))
local st3p = c(processors)
local st3ref "tests/unit/sens_t3_ref_s`st3v'_p`st3p'.mata"
local st3fd  "tests/unit/sens_t3_fdref_s`st3v'_p`st3p'.mata"
capture confirm file "`st3ref'"
local st3m = _rc
capture confirm file "`st3fd'"
if _rc | `st3m' {
    display as text "test_sens_t3: no exact-path and FD reference pair for Stata `st3v' at `st3p' processor(s) (`st3ref', `st3fd'); skipped"
    exit 8
}
display as text "test_sens_t3: references `st3ref' (exact path) and `st3fd' (finite differences)"
* The FD reference's values are kept in externals; both files define st3_ref().
do "`st3fd'"
mata:
void st3t_keepfd()
{
    string rowvector cs
    real scalar k

    cs = ("metric_full", "scalar_full", "metric_none", "metric_bigeps")
    for (k = 1; k <= cols(cs); k++) {
        st3_put("ST3_FDREF_" + cs[k], st3_ref(cs[k]))
        st3_put("ST3_FDTH_" + cs[k], st3_ref_theta(cs[k]))
    }
}
st3t_keepfd()
mata drop st3t_keepfd() st3_ref() st3_ref_theta()
end
do "`st3ref'"

mata:

real scalar st3t_rel(real matrix x, real matrix y)
{
    return(max(abs(x - y)) / max(abs(y)))
}

// Per element: |x - y| / |y|, and exact equality where y = 0 (r1)
real scalar st3t_relel(real rowvector x, real rowvector y)
{
    real scalar j, m

    m = 0
    for (j = 1; j <= cols(y); j++) {
        if (y[j] == 0) {
            if (x[j] != 0) return(.)
        }
        else if (abs(x[j] - y[j]) / abs(y[j]) > m) {
            m = abs(x[j] - y[j]) / abs(y[j])
        }
    }
    return(m)
}

void st3t_chk(string colvector FAIL, string scalar nm, real scalar ok, real scalar val)
{
    if (ok) {
        printf("{txt}  %-6s ok   (%9.2e)\n", nm, val)
    }
    else {
        printf("{err}  %-6s FAIL (%9.2e)\n", nm, val)
        FAIL = FAIL \ nm
    }
}

real scalar st3t_run()
{
    struct dwlsmv_fit scalar f
    struct dwlsmv_statmap vector maps
    string colvector FAIL
    string rowvector cs
    real colvector step, ds, eps, s, w, th, r, sg, Th, h, z0
    real matrix dOm, sens, A, D, P, Pi, Eb, Jm, dD, c6, Dp, Dm, d10
    real rowvector o, o0, o1, oa, ob, f0, f1, f2, f3
    real scalar k, c, nu, N, q, u, epsV, rhoT, l, mx, rT, ok, val, x

    FAIL = J(0, 1, "")
    u    = epsilon(1) / 2
    epsV = 2 * u
    rhoT = 0.01
    cs   = ("metric_full", "scalar_full", "metric_none", "metric_bigeps")

    // ---- r1, r2 (finite-difference path) and r1x, r2x (exact path) ------------------
    for (k = 1; k <= cols(cs); k++) {
        for (x = 0; x <= 1; x++) {
            if (x) {
                o   = st3_case_x(cs[k], epsV, rhoT)
                val = st3t_relel(o, st3_ref(cs[k]))
            }
            else {
                o   = st3_case(cs[k], epsV, rhoT)
                val = st3t_relel(o, *findexternal("ST3_FDREF_" + cs[k]))
            }
            st3t_chk(FAIL, (x ? "r1x." : "r1.") + strofreal(k), val <= 1e-12, val)

            if (st3_inputs(cs[k], f, step, ds, eps, dOm)) {
                if (x) sens = dwlsmv_sens_x(f, step, epsV, A, c, nu, ds, eps)
                else   sens = dwlsmv_sens(f, step, epsV, A, c, nu, ds, eps)
            }
            else {
                if (x) sens = dwlsmv_sens_x(f, step, epsV, A, c, nu)
                else   sens = dwlsmv_sens(f, step, epsV, A, c, nu)
                sens[., 5] = J(rows(sens), 1, 0)
            }
            Th  = (1 :+ abs(f.theta)) :* (sens[., 5] + sens[., 3] + sens[., 6])
            val = st3t_rel(Th, (x ? st3_ref_theta(cs[k]) : *findexternal("ST3_FDTH_" + cs[k])))
            st3t_chk(FAIL, (x ? "r2x." : "r2.") + strofreal(k), val <= 1e-12, val)
        }
    }

    // ---- h1-h3, on the metric fit ---------------------------------------------
    f    = *findexternal("ST3_F_metric")
    step = *findexternal("ST3_STEP_metric")
    dwlsmv_fit_stack(f.stage1, s, w, N, maps)
    th   = f.theta
    q    = rows(th)
    sg   = dwls_model(th, f.par, maps)
    r    = s - sg
    D    = dwls_delta(th, f.par, maps)
    P    = quadcross(D, w, D)
    sens = dwlsmv_sens(f, step, epsV, A, c, nu)

    Eb  = dwlsmv_sens_ebar(f, epsV)
    Pi  = abs(cholinv(P))
    c6  = (2 * (Pi * quadcross(Eb, abs(w :* r))) / (1 - c)) :/ (1 :+ abs(th))
    val = max(abs(c6 - sens[., 6]))
    st3t_chk(FAIL, "h1", c6 == sens[., 6], val)

    Jm  = dwlsmv_sens_jmat(f)
    val = max(abs(cholsolve(Jm, (D :* w)') - A))
    st3t_chk(FAIL, "h2", cholsolve(Jm, (D :* w)') == A, val)

    dD   = dwlsmv_sens_ddelta(f)
    h    = J(q, 1, 0)
    h[1] = 1e-4 * max((1, abs(th[1])))
    Dp   = dwls_delta(th + h, f.par, maps)
    Dm   = dwls_delta(th - h, f.par, maps)
    val  = max(abs((Dp - Dm) / (2 * h[1]) - dD[., 1..q]))
    st3t_chk(FAIL, "h3.a", (Dp - Dm) / (2 * h[1]) == dD[., 1..q], val)
    val = 0
    mx  = max(abs(dD))
    for (l = 1; l <= q; l++) {
        h    = J(q, 1, 0)
        h[l] = 1e-3 * max((1, abs(th[l])))
        d10  = (dwls_delta(th + h, f.par, maps) - dwls_delta(th - h, f.par, maps)) / (2 * h[l])
        x    = max(abs(d10 - dD[., ((l - 1) * q + 1)..(l * q)]))
        if (x > val) val = x
    }
    st3t_chk(FAIL, "h3.b", val <= 1e-6 * mx, val / mx)

    // ---- a1-a5 ----------------------------------------------------------------
    o1 = st3_case("metric_full", epsV, 0.01)
    o0 = st3_case("metric_full", epsV, 0)
    ok = (o1[(1..6, 8)] == o0[(1..6, 8)])
    st3t_chk(FAIL, "a1.eq", ok, max(abs(o1[(1..6, 8)] - o0[(1..6, 8)])))
    val = abs((o1[7] - o0[7]) / (0.01 * (o1[2] + o1[3] + o1[4])) - 1)
    st3t_chk(FAIL, "a1.lin", val <= 1e-12, val)

    o = st3_case("metric_none", epsV, rhoT)
    st3t_chk(FAIL, "a2", o[4] == 0, o[4])

    st3_pert(f, ds, eps, dOm)
    z0 = J(rows(ds), 1, 0)
    oa = dwlsmv_sens_t3(f, step, epsV, rhoT, ds, z0, dOm)
    ob = dwlsmv_sens_t3(f, step, epsV, rhoT, 2 * ds, z0, 2 * dOm)
    val = abs(ob[4] / (2 * oa[4]) - 1)
    st3t_chk(FAIL, "a3", val <= 1e-12 & oa[4] > 0, val)

    o  = st3_case("metric_full", epsV, rhoT)
    rT = 2 * u * abs(o[1]) / (1 + abs(o[1]))
    ok = (o[8] == o[5] + rT) & (o[7] == (1 + rhoT) * (o[2] + o[3] + o[4]) + o[5] + o[6] + rT)
    st3t_chk(FAIL, "a4", ok, abs(o[8] - o[5] - rT))

    f0 = st3_case("metric_full", 0, rhoT)
    f1 = st3_case("metric_full", u, rhoT)
    f2 = st3_case("metric_full", 2 * u, rhoT)
    f3 = st3_case("metric_full", 3 * u, rhoT)
    st3t_chk(FAIL, "a5.ev", (f0[5] == f1[5]) & (f1[5] == f3[5]), abs(f0[5] - f3[5]))
    o   = f0 - 2 * f1 + f2
    val = max(abs(o[(2, 6)])) / max(abs(f2[(2, 6)]))
    st3t_chk(FAIL, "a5.aff", val <= 1e-10, val)
    val = abs(f3[3] - 3 * f2[3] + 3 * f1[3] - f0[3]) / abs(f3[3])
    st3t_chk(FAIL, "a5.quad", val <= 1e-10, val)

    if (rows(FAIL)) printf("{err}  failed: %s\n", invtokens(FAIL'))
    return(rows(FAIL))
}

// The refusals, each expected at rc 3498.
void st3t_err(real scalar which)
{
    struct dwlsmv_fit scalar f
    real colvector step, ds, eps
    real matrix dOm
    real rowvector o

    f    = *findexternal("ST3_F_metric")
    step = *findexternal("ST3_STEP_metric")
    st3_pert(f, ds, eps, dOm)
    if (which == 1) {
        f.converged = 0
        o = dwlsmv_sens_t3(f, step, epsilon(1), 0.01)
    }
    else if (which == 2) {
        f.acov = J(0, 0, .)
        o = dwlsmv_sens_t3(f, step, epsilon(1), 0.01)
    }
    else if (which == 3) {
        f.D_stat = 0
        o = dwlsmv_sens_t3(f, step, epsilon(1), 0.01)
    }
    else if (which == 4) {
        f.T3 = 0.001
        o = dwlsmv_sens_t3(f, step, epsilon(1), 0.01)
    }
    else if (which == 5) o = dwlsmv_sens_t3(f, step, -1, 0.01)
    else if (which == 6) o = dwlsmv_sens_t3(f, step, epsilon(1), -1)
    else if (which == 7) o = dwlsmv_sens_t3(f, step, epsilon(1), 0.01, ds, eps)
    else if (which == 8) o = dwlsmv_sens_t3(f, step, epsilon(1), 0.01, ds[2::rows(ds)], eps, dOm)
    else if (which == 9) {
        dOm[2, 3] = .
        o = dwlsmv_sens_t3(f, step, epsilon(1), 0.01, ds, eps, dOm)
    }
}

end

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F -> y1 y2 y3 y4 y5), group(g) invariance(configural) metric(unitvariance)
mata: st3_setup()

local nf 0
mata: st_local("nf", strofreal(st3t_run()))

forvalues k = 1/9 {
    capture noisily mata: st3t_err(`k')
    if _rc != 3498 {
        display as error "    e`k' FAILED: rc `=_rc', expected 3498"
        local ++nf
    }
    else {
        display as text "    e`k' refused, rc 3498"
    }
}

if `nf' {
    display as error "test_sens_t3: `nf' failure(s)"
    exit 9
}
display as result "test_sens_t3: dwlsmv_sens_t3() reproduces the subsection's R on its own inputs, its helpers are dwlsmv_sens()'s Ebar and J and an exact dDelta/dtheta, rhoT, the Stage 1 inputs and epsV enter as TS7 says, and every refusal is rc 3498"
