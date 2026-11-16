*! test_condx_rootcap.do -- the root-search cap in conditional mode (nine recorded pairs).
*!
*! Executed 2026-11-12 by the main session; drafted by a read-only agent. It passed on its first run with
*! remedy_applied 0. Re-run 2026-11-13 with remedy_applied 1, after dwlsmv_cx_pair_fit() was given the
*! budget 200 (the user, 2026-11-13; a margin chosen after measuring, @sec-condx-rootcap). P_ALT_K (78)
*! was filled from the first run; its structure, 8 + 2n, is derived in @sec-condx-rootcap and n = 35 is
*! computed, not closed-form. The rest of this banner is the agent's, written before the first run:
*! expected number marked "derived" below was worked out by hand from
*! dwlsmv_rootfind() in src/mata/dwls_polychoric.mata; the first run decides
*! whether the arithmetic was right. Remove this banner only after running the
*! file and fixing what broke.
*!
*! What this pins. dwlsmv_cx_pair_fit() passes the root search the cap 40; nine
*! constructed pairs end at status 7 (cap reached) although a root exists
*! (E-096). The decided remedy: conditional
*! mode passes a larger derived budget N_BUDGET (placeholder below, one line).
*! The shared root-finder, and v1's cap of 40, do not change.
*!
*! LEVEL 1 (no data, exact, fast): dwlsmv_rootfind() on SYNTHETIC linear scores
*!   S(x) = -h (x - root), reported information c, error bound tau_S, so the
*!   Fisher step is S/c = -(h/c)(x - root) and H/c = h/c is under control. The
*!   evaluator is usable everywhere, so the endpoint check evaluates f(+-B) once
*!   each (usable), reads S(-B) > tau_S and S(+B) < -tau_S (root inside (-B, B))
*!   and returns status 1: every converged synthetic run ends (.., 1, k) with
*!   k = loop evaluations + 2. The evaluator appends each x to rows 2.. of its
*!   argument n (by reference), so the test also reads the points visited; the
*!   number of evaluator calls must equal the returned evaluation count.
*!   P-alt   h/c = 0.1:  Fisher step undershoots, alternates with a bisection.
*!   P-one   h/c = 0.51: every Fisher step accepted, all on the same side; the
*!           bracket never narrows, the end comes from the far-side probe.
*!   P-near  h/c = 1.45: bracket inside tau_rho while |S| = 1.9 tau_S at call 40.
*!   P-tiny  h/c = 1e13, S and tau_S at 1e-30: every Fisher step rejected,
*!           pure bisection; stop decided by |S| alone (tau_rho = 1.6 > bracket).
*!   control h/c = 1:    converges in 9 evaluations, identical at both caps.
*! LEVEL 2 (production, nine recorded pairs; tables embedded, no reads from
*!   E-096): (i) cap 40 through dwlsmv_rootfind() with dwlsmv_cx_rf_eval():
*!   status 7, 40 evaluations (today's behaviour at the old cap; the cap is an
*!   argument, so this stays true after the remedy); (ii) cap N_BUDGET: the status,
*!   rho and evaluation count recorded at cap 1000 (mech/trace_cap1000.log), which
*!   requires N_BUDGET above the largest loop count (145 evaluations in all, case 9);
*!   (iii) the entry point dwls_sample_stats_cx(): BEFORE the remedy rc 3498 with
*!   "status 7" in the message; AFTER it rc 0 and rho-hat as recorded. Switch below.
*!
*! Mutants for the main session (edit, run, revert); see notes.md for which of the
*! checks catches each and which mutants are NOT caught:
*!   M1 budget argument in dwlsmv_cx_pair_fit reverted to 40 (or below 145)
*!   M2 half-step safeguard (abs(stp) > abs(dxold) / 2) removed from the loop
*!   M3a stop test: drop "R[1] - L[1] <= tr"      M3b drop "abs(x[3]) <= x[5]"
*!   M4 "k >= cap" changed to "k > cap"            M5 far-side probe sign flipped
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

* ---- the two knobs ------------------------------------------------------------
* N_BUDGET: 200, the conditional-mode budget (the user, 2026-11-13); v1 stays at 40.
*   Needs N_BUDGET > 145 for Level 2 (ii) (the recorded counts), > 85 for P-alt.
local N_BUDGET 200
* remedy_applied: 0 before dwlsmv_cx_pair_fit passes the new budget, 1 after.
local remedy_applied 1
* P-alt's exact evaluation count: "to be filled from the first run" (bounds hold
* meanwhile). A missing value means not yet pinned.
local P_ALT_K 78

* Load every source, types first, as test_s2_guard.do does.
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
        display as error "test_condx_rootcap: `f' failed to load"
        exit 9
    }
}

capture mata: st_local("chk", strofreal(dwlsmv_rf_ulp(1)))
if _rc {
    display as error "test_condx_rootcap: dwlsmv_rf_ulp() is not defined (src/mata/dwls_polychoric.mata missing?)"
    exit 9
}

mata:

real scalar rcap_nb()
{
    return(`N_BUDGET')
}

real scalar rcap_pk_alt()
{
    return(`P_ALT_K')
}

// 0 if ok, else prints what and returns 1.
real scalar rcap_chk(real scalar ok, string scalar what)
{
    if (ok) {
        return(0)
    }
    printf("{err}test_condx_rootcap: %s\n", what)
    return(1)
}

// ==== LEVEL 1: synthetic evaluator =============================================

// n = (root, h, c, tau_S) in row 1; rows 2.. receive the points visited (n is passed
// by reference through dwlsmv_rootfind, so the append is seen by the caller).
// Returns (x, usable = 1, S, c, tau_S) with the linear score S = -h (x - root).
real rowvector rcap_lin(real scalar x, real matrix n,
                        real colvector a, real colvector b)
{
    real scalar root, h, c, tauS

    root = n[1, 1]
    h    = n[1, 2]
    c    = n[1, 3]
    tauS = n[1, 4]
    n    = n \ (x, 0, 0, 0)
    return((x, 1, -h * (x - root), c, tauS))
}

// dwlsmv_rootfind() with the synthetic evaluator, B = 0.999, capc = 54, the
// production call's shape. tr receives the points visited, in order.
real rowvector rcap_run(real rowvector prm, real scalar cap, real colvector tr)
{
    real matrix    n
    real colvector z
    real rowvector r

    n = prm
    z = 0
    r = dwlsmv_rootfind(&rcap_lin(), n, z, z, (., .), 0.999, cap, 54)
    if (rows(n) > 1) {
        tr = n[|2, 1 \ rows(n), 1|]
    }
    else {
        tr = J(0, 1, .)
    }
    return(r)
}

// ---- P-alt: alternating step (cases 2, 3, 4, 5, 7, 8) ---------------------------
// prm = (root 0.9, h 0.1, c 1, tau_S 1e-14). Walk (S(0) > 0, d = +1): e0, then
// 0.5, 0.75, 0.875 (S > 0), 0.9375 (S < 0, sign change): k = 5, L = 0.875,
// R = 0.9375, W0 = 0.0625 exactly, x = L (distance 0.025). Loop, r = h/c = 0.1:
// first Fisher step 0.1 * 0.025 = 0.0025 <= W0/2: accepted, lands on L's side at
// 0.9 * 0.025. The second would be 0.1 * 0.0225 = 0.00225 > half of the previous
// step 0.00125: rejected, bisection (dxold becomes W/2). After a bisection the
// near end is within W/2, the Fisher step is at most 0.05 W <= W/4: accepted.
// So Fisher and bisection strictly alternate. Each bisection halves W (the new
// bracket is a half); a Fisher step shrinks W by at most 0.1 W.
//   tau_rho = 16 tau_S / c + 8 ulp(0.9) = 1.6e-13 + 8 * 2^-53.
//   Stop = |S| <= tau_S and W <= tau_rho; once W <= tau_rho the near end is within
//   W/2 of the root, |S| <= 0.1 * W/2 <= 8.5e-15 <= tau_S, so the stop is W <= tau_rho.
//   m = ceil(log2(W0 / tau_rho)) bisections are enough (W <= W0 / 2^m), and no
//   run can finish with fewer loop evaluations than m (an evaluation at best halves W).
//   Loop evaluations in [m, 2 m]; total k = 5 + loop + 2 (endpoint check, status 1).
//   Cap 40: after 35 loop evaluations W >= W0 / 2^35 = 1.8e-12 > tau_rho: status 7, k = 40.
real scalar rcap_l1_alt()
{
    real rowvector prm, r40, r40b, rb
    real colvector t40, t40b, tb
    real scalar nf, W0, trr, m, klo, khi, pk

    nf  = 0
    prm = (0.9, 0.1, 1, 1e-14)
    r40  = rcap_run(prm, 40, t40)
    r40b = rcap_run(prm, 40, t40b)
    nf = nf + rcap_chk(r40[2] == 7, "P-alt cap 40: status is " + strofreal(r40[2]) + ", wanted 7")
    nf = nf + rcap_chk(r40[3] == 40, "P-alt cap 40: evaluations " + strofreal(r40[3]) + ", wanted 40")
    nf = nf + rcap_chk(rows(t40) == r40[3], "P-alt cap 40: evaluator calls differ from the returned count")
    nf = nf + rcap_chk(r40 == r40b & t40 == t40b, "P-alt cap 40: two calls are not bitwise identical")

    W0  = 0.0625
    trr = 16 * prm[4] / prm[3] + 8 * epsilon(1) / 2
    m   = ceil(ln(W0 / trr) / ln(2))
    klo = 5 + m + 2
    khi = 5 + 2 * m + 2
    rb  = rcap_run(prm, rcap_nb(), tb)
    nf = nf + rcap_chk(rb[2] == 1, "P-alt budget: status is " + strofreal(rb[2]) + ", wanted 1 (endpoint check, usable everywhere)")
    nf = nf + rcap_chk(rb[3] >= klo & rb[3] <= khi,
                       "P-alt budget: evaluations " + strofreal(rb[3]) + " outside the derived [" +
                       strofreal(klo) + ", " + strofreal(khi) + "]")
    nf = nf + rcap_chk(rows(tb) == rb[3], "P-alt budget: evaluator calls differ from the returned count")
    nf = nf + rcap_chk(abs(rb[1] - prm[1]) <= trr,
                       "P-alt budget: |rho - root| = " + strofreal(abs(rb[1] - prm[1])) + " exceeds tau_rho " + strofreal(trr))
    pk = rcap_pk_alt()
    if (pk < .) {
        nf = nf + rcap_chk(rb[3] == pk, "P-alt budget: evaluations " + strofreal(rb[3]) + ", pinned " + strofreal(pk))
    }
    else {
        printf("{txt}  P-alt: exact evaluation count not yet pinned (observed %g; derived bounds %g .. %g)\n", rb[3], klo, khi)
    }
    return(nf)
}

// ---- P-one: one-sided convergence (case 6) ---------------------------------------
// prm = (root 0.9, h 0.51, c 1, tau_S 7.4e-15); same walk as P-alt: k = 5,
// L = 0.875 at distance d0 = 0.025, W0 = 0.0625. r = h/c = 0.51; a Fisher step moves
// the near end to (1 - r) d = 0.49 d, on the same side. Accepted iff
// |stp_j| <= |stp_(j-1)| / 2, stp_j = r d_(j-1), d_j = 0.49 d_(j-1): stp_j / stp_(j-1)
// = 0.49 <= 0.5 (and the first, 0.01275, is <= W0/2 = 0.03125): all accepted. R never
// moves (all points on L's side), so W stays 0.0625 > tau_rho = 1.18e-13 and the stop
// test cannot pass through W; the loop continues until |S| <= tau_S, i.e.
// d_j <= tau_S / h = 1.45e-14: d_j = 0.025 * 0.49^j: d_39 = 2.07e-14 (1.43 x), d_40 =
// 1.01e-14 (0.70 x): j = 40. Evaluation k = 5 + 40 = 45. At the top |S| <= tau_S but
// W > tau_rho: far-side probe (x == L): xn = L + tau_rho/2 (above the root), S < 0,
// R = xn (k = 46, W = tau_rho/2). Top: stop passes. Endpoint check +2 (status 1):
// k = 48. Cap 40: top of loop with k = 40 means 35 Fisher steps, d_35 = 3.5e-13 >> tol,
// status 7, k = 40.
real scalar rcap_l1_one()
{
    real rowvector prm, r40, rb
    real colvector t40, tb
    real scalar nf, root, tol, trr

    nf   = 0
    prm  = (0.9, 0.51, 1, 7.4e-15)
    root = prm[1]
    tol  = prm[4] / prm[2]
    trr  = 16 * prm[4] / prm[3]
    r40  = rcap_run(prm, 40, t40)
    nf = nf + rcap_chk(r40[2] == 7, "P-one cap 40: status is " + strofreal(r40[2]) + ", wanted 7")
    nf = nf + rcap_chk(r40[3] == 40, "P-one cap 40: evaluations " + strofreal(r40[3]) + ", wanted 40")
    nf = nf + rcap_chk(rows(t40) == r40[3], "P-one cap 40: evaluator calls differ from the returned count")
    if (rows(t40) >= 40) {
        nf = nf + rcap_chk(all(t40[|6 \ 40|] :< root), "P-one cap 40: a Fisher point is not below the root (one-sided expected)")
    }
    rb = rcap_run(prm, rcap_nb(), tb)
    nf = nf + rcap_chk(rb[2] == 1, "P-one budget: status is " + strofreal(rb[2]) + ", wanted 1")
    nf = nf + rcap_chk(rb[3] >= 46 & rb[3] <= 50, "P-one budget: evaluations " + strofreal(rb[3]) + " outside 46 .. 50")
    nf = nf + rcap_chk(rb[3] == 48, "P-one budget: evaluations " + strofreal(rb[3]) + ", derived 48 (40 Fisher steps, probe, check)")
    nf = nf + rcap_chk(rows(tb) == rb[3], "P-one budget: evaluator calls differ from the returned count")
    nf = nf + rcap_chk(abs(rb[1] - root) <= tol * (1 + 1e-6), "P-one budget: |rho - root| exceeds tau_S / h")
    if (rows(tb) >= 46) {
        nf = nf + rcap_chk(all(tb[|6 \ 45|] :< root), "P-one budget: a Fisher point is not below the root")
        nf = nf + rcap_chk(tb[46] > root & tb[46] - tb[45] <= trr, "P-one budget: call 46 is not the far-side probe above the root")
    }
    return(nf)
}

// ---- P-near: bracket inside tau_rho, |S| still above tau_S (case 9) -----------------
// prm = (root 0.15, h 1.45, c 1, tau_S 7.5e-15). Walk: e0 (S(0) > 0, d = +1), then
// 0.5 (S < 0): k = 2, L = 0, R = 0.5, W0 = 0.5, x = L. r = 1.45: a Fisher step lands
// on the other side at (1 - r) d = -0.45 d. First step 1.45 * 0.15 = 0.2175 <= W0/2 = 0.25
// (accepted, inside (0, 0.5)); then stp_j / stp_(j-1) = 0.45 <= 0.5: every step accepted,
// alternating sides, bracket = the last two points. Loop evaluation j is call k = 2 + j,
// d_j = 0.15 * 0.45^j, S_j = 1.45 d_j: d_38 = 9.96e-15, S_38 = 1.44e-14 = 1.93 tau_S;
// W_38 = d_38 (1 + 1/0.45) = 3.2e-14 <= tau_rho = 1.2e-13. At the top with k = 40: stop
// fails on |S| alone, cap check returns (x, 7, 40). S_39 = 6.5e-15 = 0.87 tau_S: with the
// larger budget the top at k = 41 stops; endpoint check +2: (rho, 1, 43).
real scalar rcap_l1_near()
{
    real rowvector prm, r40, rb
    real colvector t40, tb
    real scalar nf, root, h, tauS, W, ratio

    nf   = 0
    prm  = (0.15, 1.45, 1, 7.5e-15)
    root = prm[1]
    h    = prm[2]
    tauS = prm[4]
    r40  = rcap_run(prm, 40, t40)
    nf = nf + rcap_chk(r40[2] == 7, "P-near cap 40: status is " + strofreal(r40[2]) + ", wanted 7")
    nf = nf + rcap_chk(r40[3] == 40, "P-near cap 40: evaluations " + strofreal(r40[3]) + ", wanted 40")
    nf = nf + rcap_chk(rows(t40) == r40[3], "P-near cap 40: evaluator calls differ from the returned count")
    if (rows(t40) >= 40) {
        // the state at call 40: the last two points straddle the root, the bracket is
        // inside tau_rho (16 tau_S / c; the 8 ulp term is 2e-16 and left out), |S| in (1, 3) tau_S
        W     = abs(t40[40] - t40[39])
        ratio = h * abs(t40[40] - root) / tauS
        nf = nf + rcap_chk((t40[40] - root) * (t40[39] - root) < 0, "P-near cap 40: the last two points do not straddle the root")
        nf = nf + rcap_chk(W <= 16 * tauS / prm[3], "P-near cap 40: bracket " + strofreal(W) + " is not inside tau_rho")
        nf = nf + rcap_chk(ratio > 1 & ratio < 3, "P-near cap 40: |S| / tau_S is " + strofreal(ratio) + ", wanted in (1, 3)")
        nf = nf + rcap_chk(r40[1] == t40[40], "P-near cap 40: rho is not the last point")
    }
    rb = rcap_run(prm, rcap_nb(), tb)
    nf = nf + rcap_chk(rb[2] == 1, "P-near budget: status is " + strofreal(rb[2]) + ", wanted 1")
    nf = nf + rcap_chk(rb[3] >= 42 & rb[3] <= 44, "P-near budget: evaluations " + strofreal(rb[3]) + " outside 42 .. 44")
    nf = nf + rcap_chk(rb[3] == 43, "P-near budget: evaluations " + strofreal(rb[3]) + ", derived 43")
    nf = nf + rcap_chk(rows(tb) == rb[3], "P-near budget: evaluator calls differ from the returned count")
    nf = nf + rcap_chk(abs(rb[1] - root) <= tauS / h * (1 + 1e-6), "P-near budget: |rho - root| exceeds tau_S / h")
    return(nf)
}

// ---- P-tiny: tiny score, pure bisection (case 1) -----------------------------------
// prm = (root -0.98434653, h 1e-16, c 1e-29, tau_S 1e-30): H/c = 1e13, S near the root
// at the 1e-30 scale. Walk (S(0) < 0, d = -1): e0, then -0.5, -0.75, -0.875, -0.9375,
// -0.96875 (S < 0), -0.984375 (below the root, S > 0: sign change): k = 7, L = -0.984375,
// R = -0.96875, W0 = 2^-6. tau_rho = 16 tau_S / c = 1.6 >> W: the stop is decided by
// |S| <= tau_S, i.e. the near end within delta = tau_S / h = 1e-14 of the root. A Fisher
// step is r d = 1e13 d > 0.1 as long as d > delta, outside the bracket or above half the
// last step: always rejected, so every loop evaluation is the midpoint (L + R) / 2 (dyadic,
// exact). The prediction is simulated below with the same expressions. Total
// k = 7 + m* + 2. m* is about 40 (W0 / 2^m* near 2 delta); the test requires m* >= 34 so
// that cap 40 is reached (top of loop at k = 7 + m* > 40 would otherwise stop first).
real scalar rcap_l1_tiny()
{
    real rowvector prm, r40, rb
    real colvector t40, tb, pred
    real scalar nf, root, h, tauS, L, R, xm, m

    nf   = 0
    prm  = (-0.98434653, 1e-16, 1e-29, 1e-30)
    root = prm[1]
    h    = prm[2]
    tauS = prm[4]
    L    = -0.984375
    R    = -0.96875
    pred = J(0, 1, .)
    m    = 0
    while (min((h * abs(L - root), h * abs(R - root))) > tauS & m < 60) {
        xm   = (L + R) / 2
        pred = pred \ xm
        if (-h * (xm - root) > 0) {
            L = xm
        }
        else {
            R = xm
        }
        m = m + 1
    }
    if (m < 34 | m >= 60) {
        return(rcap_chk(0, "P-tiny: simulated m* = " + strofreal(m) + " is outside [34, 60): pick another root"))
    }
    r40 = rcap_run(prm, 40, t40)
    nf = nf + rcap_chk(r40[2] == 7, "P-tiny cap 40: status is " + strofreal(r40[2]) + ", wanted 7")
    nf = nf + rcap_chk(r40[3] == 40, "P-tiny cap 40: evaluations " + strofreal(r40[3]) + ", wanted 40")
    nf = nf + rcap_chk(rows(t40) == r40[3], "P-tiny cap 40: evaluator calls differ from the returned count")
    if (rows(t40) >= 40) {
        nf = nf + rcap_chk(t40[|8 \ 40|] == pred[|1 \ 33|], "P-tiny cap 40: the points are not the pure bisection midpoints")
    }
    rb = rcap_run(prm, rcap_nb(), tb)
    nf = nf + rcap_chk(rb[2] == 1, "P-tiny budget: status is " + strofreal(rb[2]) + ", wanted 1")
    nf = nf + rcap_chk(rb[3] == 7 + m + 2, "P-tiny budget: evaluations " + strofreal(rb[3]) + ", derived " + strofreal(7 + m + 2))
    nf = nf + rcap_chk(rows(tb) == rb[3], "P-tiny budget: evaluator calls differ from the returned count")
    if (rows(tb) >= 7 + m) {
        nf = nf + rcap_chk(tb[|8 \ 7 + m|] == pred, "P-tiny budget: the points are not the pure bisection midpoints")
    }
    nf = nf + rcap_chk(abs(rb[1] - root) <= tauS / h * (1 + 1e-9), "P-tiny budget: |rho - root| exceeds tau_S / h")
    return(nf)
}

// ---- control: well conditioned, h/c = 1 ----------------------------------------------
// prm = (0.9, 1, 1, 1e-14). Walk as P-alt: k = 5, L = 0.875 (S = 0.025), R = 0.9375.
// Fisher step stp = 0.025 (exact: 0.9 - 0.875 is exact), accepted (<= W0/2 = 0.03125,
// inside), xn = 0.875 + 0.025 = 0.9 exactly, S = -0 (not > 0) so R = e, x = R,
// |S| = 0 <= tau_S, k = 6. W = 0.025 > tau_rho: far-side probe, x != L so
// xn = 0.9 - tau_rho / 2 (S > 0: L = e), k = 7, W = tau_rho / 2. Top: stop. Endpoint
// check +2, status 1: (0.9, 1, 9). The cap is never reached, so cap 40 and the budget
// give identical results: a larger budget changes nothing for a search that converges.
real scalar rcap_l1_ctl()
{
    real rowvector prm, r40, rb
    real colvector t40, tb
    real scalar nf

    nf  = 0
    prm = (0.9, 1, 1, 1e-14)
    r40 = rcap_run(prm, 40, t40)
    rb  = rcap_run(prm, rcap_nb(), tb)
    nf = nf + rcap_chk(r40 == rb & t40 == tb, "control: cap 40 and the budget differ (result or points visited)")
    nf = nf + rcap_chk(rb[1] == 0.9 & rb[2] == 1 & rb[3] == 9,
                       "control: result (" + strofreal(rb[1]) + ", " + strofreal(rb[2]) + ", " + strofreal(rb[3]) + "), derived (0.9, 1, 9)")
    nf = nf + rcap_chk(rows(tb) == rb[3], "control: evaluator calls differ from the returned count")
    return(nf)
}

// ==== LEVEL 2: production on the nine recorded pairs =============================

// The nine tables of cases.txt (rows: item i's category; x = 0 then x = 1).
void rcap_tab(real scalar c, real matrix n0, real matrix n1)
{
    if (c == 1) {
        n0 = (0, 2 \ 0, 0 \ 35, 19)
        n1 = (0, 1 \ 0, 1 \ 267, 105)
    }
    else if (c == 2) {
        n0 = (1, 0, 23 \ 0, 0, 12 \ 0, 0, 11 \ 0, 0, 9)
        n1 = (3, 3, 66 \ 0, 0, 12 \ 0, 0, 3 \ 0, 0, 0)
    }
    else if (c == 3) {
        n0 = (5, 0, 0 \ 9, 0, 0 \ 38, 0, 0 \ 7, 1, 0 \ 0, 0, 0)
        n1 = (1, 0, 0 \ 3, 0, 0 \ 98, 0, 0 \ 89, 0, 0 \ 14, 7, 1)
    }
    else if (c == 4) {
        n0 = (4, 0, 0 \ 0, 0, 0 \ 7, 1, 3)
        n1 = (1, 17, 5 \ 0, 0, 8 \ 0, 0, 0)
    }
    else if (c == 5) {
        n0 = (7, 0, 0, 15, 0 \ 35, 0, 0, 0, 0)
        n1 = (0, 0, 0, 0, 3 \ 4, 1, 5, 0, 0)
    }
    else if (c == 6) {
        n0 = (1, 0 \ 0, 0 \ 2, 0 \ 0, 0 \ 2, 12)
        n1 = (0, 0 \ 6, 0 \ 7, 0 \ 1, 3 \ 0, 0)
    }
    else if (c == 7) {
        n0 = (0, 18, 0, 0, 0 \ 0, 73, 0, 6, 13)
        n1 = (1, 0, 0, 0, 0 \ 4, 0, 1, 0, 0)
    }
    else if (c == 8) {
        n0 = (0, 0, 0, 21, 0 \ 0, 4, 0, 25, 5)
        n1 = (0, 3, 56, 0, 44 \ 6, 0, 5, 0, 0)
    }
    else {
        n0 = (2, 7 \ 6, 25 \ 0, 7)
        n1 = (23, 0 \ 0, 0 \ 0, 15)
    }
}

// Recorded at cap 1000 (mech/trace_cap1000.log, lines TR|result|rho|status|evaluations):
// (rho, status, evaluations). Case 2 is status 6 (endpoint check), the rest status 1.
real rowvector rcap_rec(real scalar c)
{
    if (c == 1) {
        return((-0.9843465284752710, 1, 96))
    }
    else if (c == 2) {
        return((0.8655750793410987, 6, 139))
    }
    else if (c == 3) {
        return((0.9924512747566461, 1, 135))
    }
    else if (c == 4) {
        return((0.9713999082126308, 1, 139))
    }
    else if (c == 5) {
        return((-0.9919302210346441, 1, 142))
    }
    else if (c == 6) {
        return((0.9483529237199626, 1, 96))
    }
    else if (c == 7) {
        return((0.8336534213472614, 1, 135))
    }
    else if (c == 8) {
        return((-0.6950955978491219, 1, 120))
    }
    return((0.7851204790679860, 1, 145))
}

// Raw data from two tables (cap_replay.do's mc_data).
void rcap_data(real matrix n0, real matrix n1, real colvector yi,
               real colvector yj, real colvector x)
{
    real scalar g, r, c
    real matrix n

    yi = J(0, 1, .)
    yj = J(0, 1, .)
    x  = J(0, 1, .)
    for (g = 0; g <= 1; g++) {
        n = (g == 0 ? n0 : n1)
        for (r = 1; r <= rows(n); r++) {
            for (c = 1; c <= cols(n); c++) {
                yi = yi \ J(n[r, c], 1, r - 1)
                yj = yj \ J(n[r, c], 1, c - 1)
                x  = x  \ J(n[r, c], 1, g)
            }
        }
    }
}

// Item fits, the stacked tables and the shifted bounds exactly as dwlsmv_cx_pair_fit()
// builds them (mech/trace.do's tr_case). Outputs by reference.
void rcap_setup(real matrix n0, real matrix n1, real matrix n,
                real colvector ai, real colvector bj)
{
    real colvector yi, yj, x, thi, thj, ti, tj
    real scalar it, Ki, Kj

    rcap_data(n0, n1, yi, yj, x)
    Ki  = rows(n0)
    Kj  = cols(n0)
    it  = .
    thi = dwlsmv_cx_item_fit(yi, x, Ki, 1, it)
    thj = dwlsmv_cx_item_fit(yj, x, Kj, 2, it)
    ti  = thi[|1 \ Ki - 1|]
    tj  = thj[|1 \ Kj - 1|]
    n   = (dwlsmv_crosstab(select(yi, x :== 0), select(yj, x :== 0), Ki, Kj) \
           dwlsmv_crosstab(select(yi, x :== 1), select(yj, x :== 1), Ki, Kj))
    ai  = (ti \ (ti :- thi[Ki]))
    bj  = (tj \ (tj :- thj[Kj]))
}

// (i) and (ii) for case c: the root search through dwlsmv_rootfind() with the
// production evaluator. The cap is passed explicitly, so this does not depend on
// what dwlsmv_cx_pair_fit() passes.
real scalar rcap_l2(real scalar c)
{
    real matrix n0, n1, n
    real colvector ai, bj
    real rowvector r40, rb, rec, r40b
    real scalar nf
    string scalar tag

    nf  = 0
    tag = "L2 case " + strofreal(c)
    rcap_tab(c, n0, n1)
    rcap_setup(n0, n1, n, ai, bj)
    rec = rcap_rec(c)

    // (i) today's behaviour at the old cap
    r40 = dwlsmv_rootfind(&dwlsmv_cx_rf_eval(), n, ai, bj, (., .), 0.999, 40, 54)
    nf = nf + rcap_chk(r40[2] == 7, tag + " cap 40: status is " + strofreal(r40[2]) + ", wanted 7")
    nf = nf + rcap_chk(r40[3] == 40, tag + " cap 40: evaluations " + strofreal(r40[3]) + ", wanted 40")
    r40b = dwlsmv_rootfind(&dwlsmv_cx_rf_eval(), n, ai, bj, (., .), 0.999, 40, 54)
    nf = nf + rcap_chk(r40 == r40b, tag + " cap 40: two calls are not bitwise identical")

    // (ii) the budget: what the cap-1000 run recorded (the path is the same as long
    // as the budget is not reached, so the evaluation count is the same too)
    rb = dwlsmv_rootfind(&dwlsmv_cx_rf_eval(), n, ai, bj, (., .), 0.999, rcap_nb(), 54)
    nf = nf + rcap_chk(rb[2] == rec[2], tag + " budget: status " + strofreal(rb[2]) + ", recorded " + strofreal(rec[2]))
    nf = nf + rcap_chk(abs(rb[1] - rec[1]) <= 1e-9, tag + " budget: rho " + strofreal(rb[1], "%21.16f") + ", recorded " + strofreal(rec[1], "%21.16f"))
    nf = nf + rcap_chk(rb[3] == rec[3], tag + " budget: evaluations " + strofreal(rb[3]) + ", recorded " + strofreal(rec[3]))
    return(nf)
}

// (iii) the entry point: returns rho-hat of the one pair; aborts (3498) if the pair
// cannot enter.
real scalar rcap_entry(real scalar c)
{
    real matrix n0, n1
    real colvector yi, yj, x, K
    struct dwlsmv_stage1 scalar st

    rcap_tab(c, n0, n1)
    rcap_data(n0, n1, yi, yj, x)
    K  = (rows(n0) \ cols(n0))
    st = dwls_sample_stats_cx((yi, yj), K, x)
    return(st.s[st.map.D])
}

// Does the log at f contain the text has? Returns the number of failures.
real scalar rcap_msg(string scalar f, string scalar has, string scalar tag)
{
    string scalar txt

    txt = invtokens(cat(f)', char(10))
    if (!strpos(txt, has)) {
        printf("{err}test_condx_rootcap: case %s: the log lacks: %s\n", tag, has)
        return(1)
    }
    return(0)
}

// The status-6 note: present with all three of its parts when want = 1, absent when
// want = 0. Returns the number of failures.
real scalar rcap_note(string scalar f, real scalar want, string scalar tag)
{
    string scalar    txt
    string rowvector has
    real scalar      j, nf

    txt = invtokens(cat(f)', char(10))
    has = ("the pair y1~~y2 enters at status 6",
           "could not be certified at an end of the range",
           "not ordinary convergence")
    nf  = 0
    for (j = 1; j <= cols(has); j++) {
        if (want & !strpos(txt, has[j])) {
            printf("{err}test_condx_rootcap: case %s: the log lacks: %s\n", tag, has[j])
            nf = nf + 1
        }
    }
    if (!want & strpos(txt, "status 6")) {
        printf("{err}test_condx_rootcap: case %s: a status-6 note on a status-1 pair\n", tag)
        nf = nf + 1
    }
    return(nf)
}

end

local nfail 0
tempfile L

* ---- Level 1: the algorithm on synthetic evaluators ------------------------------
foreach c in rcap_l1_alt rcap_l1_one rcap_l1_near rcap_l1_tiny rcap_l1_ctl {
    capture noisily mata: st_local("nok", strofreal(`c'()))
    if _rc {
        display as error "test_condx_rootcap: `c'() aborted, rc " _rc
        local nfail = `nfail' + 1
    }
    else {
        display as text "  `c'(): `nok' failure(s)"
        local nfail = `nfail' + `nok'
    }
}

* ---- Level 2 (i), (ii): production root search on the nine pairs -----------------
forvalues c = 1/9 {
    capture noisily mata: st_local("nok", strofreal(rcap_l2(`c')))
    if _rc {
        display as error "test_condx_rootcap: rcap_l2(`c') aborted, rc " _rc
        local nfail = `nfail' + 1
    }
    else {
        display as text "  rcap_l2(`c'): `nok' failure(s)"
        local nfail = `nfail' + `nok'
    }
}

* ---- Level 2 (iii): the entry point, before or after the remedy -------------------
* remedy_applied 0: rc 3498 and "status 7" for the pair y1~~y2 in the message.
* remedy_applied 1: rc 0 and rho-hat within 1e-9 of the recorded value. Case 2 is
*   status 6 from the endpoint check at 139 evaluations; v1's entry rule (status 1, or
*   status 6 with evaluations other than 3) lets it enter, so it completes too. That is
*   what the code does, and the user agreed it enters (2026-11-13): in Mata the score
*   is negative and above its error bound up to rho = 0.9988 and an exact zero from
*   0.99889 to B. Not independently verified: the true sign between 0.9988 and B.
forvalues c = 1/9 {
    capture scalar drop rcap_rho
    quietly log using "`L'_e`c'.log", text replace name(rcap)
    capture noisily mata: st_numscalar("rcap_rho", rcap_entry(`c'))
    local rc = _rc
    quietly log close rcap
    if `remedy_applied' {
        if `rc' != 0 {
            display as error "test_condx_rootcap: entry point case `c' gave rc `rc', wanted 0 (remedy applied)"
            local nfail = `nfail' + 1
        }
        else {
            * Status 6 is never silent: case 2 (status 6 at the endpoint check) must be
            * named in a note that says it is not ordinary convergence; the eight
            * status-1 pairs must print no such note.
            local tmpL "`L'_e`c'.log"
            capture noisily mata: st_local("nm", strofreal(rcap_note(st_local("tmpL"), `c' == 2, "e`c'")))
            if _rc {
                display as error "test_condx_rootcap: reading the log of entry case `c' failed, rc " _rc
                local nfail = `nfail' + 1
            }
            else local nfail = `nfail' + `nm'
            mata: st_local("rec", strofreal(rcap_rec(`c')[1], "%21.16f"))
            if abs(scalar(rcap_rho) - `rec') > 1e-9 {
                display as error "test_condx_rootcap: entry point case `c': rho-hat " %21.16f scalar(rcap_rho) ", recorded `rec'"
                local nfail = `nfail' + 1
            }
        }
    }
    else {
        if `rc' != 3498 {
            display as error "test_condx_rootcap: entry point case `c' gave rc `rc', wanted 3498 (remedy not applied)"
            local nfail = `nfail' + 1
        }
        else {
            local tmpL "`L'_e`c'.log"
            local want "status 7"
            capture noisily mata: st_local("nm", strofreal(rcap_msg(st_local("tmpL"), st_local("want"), "e`c'")))
            if _rc {
                display as error "test_condx_rootcap: reading the log of entry case `c' failed, rc " _rc
                local nfail = `nfail' + 1
            }
            else local nfail = `nfail' + `nm'
        }
    }
}

if `nfail' {
    display as error _n "test_condx_rootcap: `nfail' failure(s)"
    exit 9
}

display as result "test_condx_rootcap: the four cap patterns and the controls hold on synthetic evaluators; the nine recorded pairs end at status 7 after 40 evaluations and at the recorded root with the budget; entry point as expected for remedy_applied = `remedy_applied'"
exit 0
