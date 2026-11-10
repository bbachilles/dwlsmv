*! Module 2: polychoric correlations and the Gate 3 information quantities.
*!
*! Executed under StataMP 18 on Windows, 2026-09-19, and it passes. Its first
*! run failed L6 at 1.29e-08 against a 1e-08 tier, which was the finding that
*! produced the Fisher-scoring polish in dwls_polychoric. Since 2026-09-30 the
*! polish and optimize() are replaced by option B's root-finder (derivation
*! @sec-rootfind; tests/unit/test_rootfind.do).
*!
*! Levelled the way the derivation file levels its own verification, so a
*! failure localises rather than just reporting that something is wrong. Each
*! level is built only from things the levels below it have already checked.
*!
*!   L0  Phi2 itself, against the closed form Phi2(0,0;rho) = 1/4 +
*!       asin(rho)/(2 pi). Mata's binormal() is the production bivariate normal
*!       and everything else here rests on it, so it is checked first and
*!       against an identity rather than against another implementation.
*!
*!   L1  The derivatives of pi, against central differences of pi itself.
*!       Everything downstream is built from these. Also the two structural
*!       identities that cost nothing: cells sum to 1, and d/drho sums to 0 --
*!       the latter being why the score has mean zero at the truth.
*!
*!   L2  The information quantities. I_rho_rho against numerical
*!       differentiation of E[g], and I_rho_a BOTH ways -- against numerical
*!       differentiation and against the w-weighted compact form. The second
*!       comparison is the one that catches the two coupled rows being taken
*!       in the wrong order, which is the most likely index error in Module 2
*!       and which the numerical check alone would also catch but less
*!       specifically.
*!
*!   L3  The influence function's exact moments, by summing over the table
*!       rather than simulating. Reproduces the derivation file's finding that
*!       the NAIVE influence function -- propagation dropped -- is also exactly
*!       centred, and understates N*Var(rho-hat) by about 8%. This is the trap
*!       the project's notes record as invisible, and pinning the number here means a
*!       future Module 4 that silently drops the term has a figure to fail
*!       against.
*!
*!   L5  Reduction to the binary oracle of spec 4.5.3, over a 60-point grid in
*!       (q1, q2, rho). The general machinery must reproduce the independently
*!       derived 2x2 closed form, which was obtained by a DIFFERENT route --
*!       the implicit function theorem on the saturated inversion. Two routes
*!       agreeing is the corroboration Gate 3 asked for.
*!
*!   L6  rho-hat itself against lavaan on case_D, all six pairs, at the spec
*!       7.2 tier (1e-8).
*!   L6b rho-hat against Mplus on case D (tier 3), all six pairs.
*!       fixtures/expected/mplus_case_D_correlations.csv holds Mplus's
*!       correlations to the three decimals it prints, so a printed value is
*!       within 5e-4 of Mplus's own estimate; 1e-5 is added for Mplus's
*!       estimate differing from ours, the order by which its theta-hat differs
*!       from lavaan's (test_acov.do check 6). Tolerance 5.1e-4, chosen BEFORE
*!       any run and not measured for rho-hat. Added 2026-10-05 after a review
*!       found the header promised this comparison and no code asserted it.
*!       First run, StataMP 18: worst gap 3.80166e-4 over the six pairs (y2,y3
*!       against Mplus's .476), inside the tolerance. That is also what
*!       lavaan's own rho-hat gives against the same table.
*!
*!   L7  The analytic score at rho-hat, which must vanish. Read this one
*!       carefully: it was an INDEPENDENT check when the optimizer was purely
*!       derivative-free, and it is not any more: since 2026-09-30 (option B)
*!       dwls_polychoric() finds the root of that same score. It would vanish
*!       even for a derivative wrong by a constant factor. It still confirms
*!       the root-finder stopped at a root rather than stalled, but it is not
*!       evidence the derivative is right.
*!
*!   L7b The check that replaces what L7 used to provide. The analytic score
*!       against a central difference of the log-likelihood AWAY from the
*!       optimum, where both are large and the comparison has content. Immune
*!       to the root-finder, which never evaluates those points. This is
*!       what would now catch a sign or scale error in
*!       dwlsmv_dcellprob_rho(), which I_rho_rho cannot notice because it
*!       squares that quantity.
*!
*!   L8  The saturated 2x2 identity, on case_A. A 2x2 table has 3 free cells
*!       and 3 parameters, so the MLE reproduces the observed table EXACTLY
*!       (spec 4.5.3) and Phi2(tau1,tau2;rho-hat) must equal the observed p00
*!       to machine precision. No external reference, no tolerance argument,
*!       and it holds for binary and binary only -- the same over-identification
*!       boundary that stops 4.5.3's closed form extending to K > 2.
*!
*!   L9  Near-duplicate pairs (added 2026-09-26, test-first: before the fix one
*!       table aborted with rc 3498 and the other stalled at rho-hat .9731 with
*!       |S|/N = .55). Two tables whose empty corner cells underflow near the
*!       root (derivation @sec-emptycells), from fixtures/expected/
*!       highrho_polychoric.csv: no abort, converged, the score vanishes,
*!       rho-hat within spec 7.2's 1e-8 of the tier-1 reference, and the
*!       shipped I_rho_rho still the per-observation crossproduct. L9b pins the
*!       rule on a cell that is zero under any rounding: an EMPTY zero cell
*!       leaves the log-likelihood finite and pairinfo zeroes its u; an
*!       OBSERVED one is infeasible, and pairinfo aborts on it (L9c, rc 3498).
*!       Added 2026-09-28 from the Module 2 adversarial review: L9 asserts u is
*!       exactly 0 on empty underflowed cells (mutant M30); L9 and L9b count
*!       missing values term by term, since sum() skips them; L9d holds a boundary 2x2
*!       pair to the flagged-pair contract (converged 0, |rho-hat| <= .999:
*!       rho-hat is exactly 0.999 there); and cellprob_raw's own rho bound is
*!       an error path (mutant M33).
*!
*! LEVEL 4 IS ABSENT ON PURPOSE. The derivation file's decisive test is the
*! pointwise expansion check -- that rho-hat minus rho equals the mean of psi
*! plus o_p(N^-1/2) -- and it is the only test that separates the correct
*! influence function from the naive one at feasible replication counts. It
*! belongs to Module 4, which is where psi is actually assembled; Module 2
*! returns ingredients and never forms psi. Putting it here would mean building
*! psi in the test, which is precisely the assembly the API keeps in one place.
*! Recorded so its absence reads as a decision rather than an oversight.
*!
*! The body is a function because it touches a struct.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in tests/dwlsmv_assert.mata                  ///
             src/mata/dwlsmv_types.mata                ///
             src/mata/dwlsmv_statmap_build.mata        ///
             src/mata/dwlsmv_pos_rho.mata              ///
             src/mata/dwlsmv_catcounts.mata            ///
             src/mata/dwlsmv_catdegen.mata             ///
             src/mata/dwlsmv_cumprop.mata              ///
             src/mata/dwlsmv_thresholds.mata           ///
             src/mata/dwlsmv_psi_tau_table.mata        ///
             src/mata/dwlsmv_cellprob.mata             ///
             src/mata/dwlsmv_dcellprob_rho.mata        ///
             src/mata/dwlsmv_wmat.mata                 ///
             src/mata/dwlsmv_crosstab.mata             ///
             src/mata/dwlsmv_pairinfo.mata             ///
             src/mata/dwls_polychoric.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_polychoric: `f' failed to load"
        exit 9
    }
}

* Reference values. Read as text and converted with strtoreal: import
* delimited types numeric columns as float, which would cap every comparison
* below at about 4e-08. The float-import trap.
import delimited using "fixtures/expected/lavaan_case_D_wls_obs.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: LAB_LAV = st_sdata(., "label")
mata: VAL_LAV = strtoreal(st_sdata(., "value"))

* Mplus's polychoric correlations, three decimals, lower triangle (NA above).
* stringcols(_all) so nothing is typed as float; strtoreal turns NA into missing.
import delimited using "fixtures/expected/mplus_case_D_correlations.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: MPL_LAB = st_sdata(., 1)
mata: MPL = strtoreal(st_sdata(., (2, 3, 4, 5)))

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

import delimited using "fixtures/data/case_A.csv", clear varnames(1) delimiters(",")
mata: YA = st_data(., .)

mata:

real scalar dwlsmv_t_find2(string colvector pool, string scalar want)
{
    real scalar i

    for (i = 1; i <= rows(pool); i++) {
        if (pool[i] == want) {
            return(i)
        }
    }
    return(0)
}

// Expected score E[g] with the TRUE cell probabilities held fixed, as a
// function of the arguments -- the quantity L2 differentiates numerically.
real scalar dwlsmv_t_Eg(real matrix Ptrue, real colvector a, real colvector b,
                        real scalar rho)
{
    return(sum(Ptrue :* (dwlsmv_dcellprob_rho(a, b, rho) :/
                         dwlsmv_cellprob(a, b, rho))))
}

// psi_rho as a K1 x K2 lookup table, assembled here ONLY to check moments.
// Module 4 owns the real assembly; this is a local reconstruction so L3 can
// evaluate exact moments without simulating.
real matrix dwlsmv_t_psirho(real colvector a, real colvector b,
                            real scalar rho, real scalar drop_prop)
{
    struct dwlsmv_pairfit scalar f
    real matrix P, Dr, u, psi, pa, pb, Wa, Wb, uT
    real colvector qa, qb, Ira, Irb
    real scalar K1, K2, r, s, c, d, Irr

    K1 = rows(a) + 1
    K2 = rows(b) + 1
    P  = dwlsmv_cellprob(a, b, rho)
    Dr = dwlsmv_dcellprob_rho(a, b, rho)
    u  = Dr :/ P

    Irr = sum((Dr :* Dr) :/ P)
    Wa  = dwlsmv_wmat(a, b, rho)
    Ira = J(K1 - 1, 1, 0)
    for (c = 1; c <= K1 - 1; c++) {
        Ira[c] = sum(Wa[c, .] :* (u[c, .] - u[c + 1, .]))
    }
    Wb  = dwlsmv_wmat(b, a, rho)
    uT  = u'
    Irb = J(K2 - 1, 1, 0)
    for (d = 1; d <= K2 - 1; d++) {
        Irb[d] = sum(Wb[d, .] :* (uT[d, .] - uT[d + 1, .]))
    }

    // Stage 1 staircases at the POPULATION q implied by these thresholds --
    // this is a theoretical-moment calculation, not a sample one.
    qa = normal(a)
    qb = normal(b)
    pa = dwlsmv_psi_tau_table(qa)
    pb = dwlsmv_psi_tau_table(qb)

    psi = J(K1, K2, 0)
    for (r = 1; r <= K1; r++) {
        for (s = 1; s <= K2; s++) {
            if (drop_prop) {
                psi[r, s] = u[r, s] / Irr
            }
            else {
                psi[r, s] = (u[r, s] - sum(Ira :* pa[r, .]')
                                     - sum(Irb :* pb[s, .]')) / Irr
            }
        }
    }
    return(psi)
}

real scalar dwlsmv_t_polychoric(real matrix Y, string colvector LAB,
                                real colvector VAL, real matrix YA,
                                string colvector MPL_LAB, real matrix MPL)
{
    struct dwlsmv_statmap scalar m
    struct dwlsmv_pairfit scalar f
    string colvector FAIL
    real colvector a, b, K, ti, tj, uvec, vvec
    real matrix P, Dr, u, Wa, n, psiG, psiN, Pm, Wm, Da, V
    real scalar e3, e4, e5, Nn, rI
    real scalar h, rho, g1, g2, g3, gs, gz, i, j, c, d, s, hit
    real scalar Irr, num, alt, worst, worstM, okM, q1, q2, rr, A, B, sg, p0, d2, cf
    real scalar e1, e2, nEg, r1, r2, vg, vn, pct, sc, K1, K2

    FAIL = J(0, 1, "")
    h    = 1e-5

    // Design B of the derivation file: K = 4 and K = 3, rho = 0.5.
    a   = invnormal((0.15 \ 0.50 \ 0.80))
    b   = invnormal((0.40 \ 0.75))
    rho = 0.50

    // ---- L0: Phi2 against a closed form ----------------------------------
    // Phi2(0,0;rho) = 1/4 + asin(rho)/(2 pi). No second implementation
    // involved, which is the point.
    worst = 0
    for (i = -9; i <= 9; i++) {
        rr    = i / 10
        worst = max((worst, abs(binormal(0, 0, rr) - (0.25 + asin(rr) / (2 * pi())))))
    }
    printf("{txt}  L0 Phi2 identity        : worst gap %g over 19 rho values\n", worst)
    if (worst > 1e-14) {
        FAIL = FAIL \ ("L0: binormal disagrees with the closed form by " + strofreal(worst))
    }

    // ---- L1: derivatives of pi against central differences ---------------
    Dr = dwlsmv_dcellprob_rho(a, b, rho)
    g1 = max(abs(Dr - (dwlsmv_cellprob(a, b, rho + h) -
                       dwlsmv_cellprob(a, b, rho - h)) :/ (2 * h)))

    g2 = 0
    for (c = 1; c <= rows(a); c++) {
        ti = a; ti[c] = ti[c] + h
        tj = a; tj[c] = tj[c] - h
        // d pi / d a_c from the w-matrix: +w in row c, -w in row c+1.
        Wa = dwlsmv_wmat(a, b, rho)
        P  = J(rows(a) + 1, rows(b) + 1, 0)
        P[c, .]     =  Wa[c, .]
        P[c + 1, .] = -Wa[c, .]
        g2 = max((g2, max(abs(P - (dwlsmv_cellprob(ti, b, rho) -
                                   dwlsmv_cellprob(tj, b, rho)) :/ (2 * h)))))
    }

    g3 = 0
    for (d = 1; d <= rows(b); d++) {
        ti = b; ti[d] = ti[d] + h
        tj = b; tj[d] = tj[d] - h
        Wa = dwlsmv_wmat(b, a, rho)
        P  = J(rows(b) + 1, rows(a) + 1, 0)
        P[d, .]     =  Wa[d, .]
        P[d + 1, .] = -Wa[d, .]
        g3 = max((g3, max(abs(P' - (dwlsmv_cellprob(a, ti, rho) -
                                    dwlsmv_cellprob(a, tj, rho)) :/ (2 * h)))))
    }

    gs = abs(sum(dwlsmv_cellprob(a, b, rho)) - 1)
    gz = abs(sum(Dr))

    printf("{txt}  L1 dpi/drho vs numeric  : %g\n", g1)
    printf("{txt}  L1 dpi/da_c vs numeric  : %g\n", g2)
    printf("{txt}  L1 dpi/db_d vs numeric  : %g\n", g3)
    printf("{txt}  L1 cells sum to 1       : %g      d/drho sums to 0 : %g\n", gs, gz)

    if (g1 > 1e-8)  FAIL = FAIL \ ("L1: dpi/drho off by "  + strofreal(g1))
    if (g2 > 1e-8)  FAIL = FAIL \ ("L1: dpi/da_c off by "  + strofreal(g2))
    if (g3 > 1e-8)  FAIL = FAIL \ ("L1: dpi/db_d off by "  + strofreal(g3))
    if (gs > 1e-12) FAIL = FAIL \ ("L1: cells do not sum to 1, off by " + strofreal(gs))
    if (gz > 1e-12) FAIL = FAIL \ ("L1: d/drho does not sum to 0, off by " + strofreal(gz))

    // ---- L2: the information quantities ----------------------------------
    P   = dwlsmv_cellprob(a, b, rho)
    u   = Dr :/ P
    Irr = sum((Dr :* Dr) :/ P)

    nEg = abs(dwlsmv_t_Eg(P, a, b, rho))
    e1  = abs(Irr + (dwlsmv_t_Eg(P, a, b, rho + h) -
                     dwlsmv_t_Eg(P, a, b, rho - h)) / (2 * h))

    Wa = dwlsmv_wmat(a, b, rho)
    e2 = 0
    for (c = 1; c <= rows(a); c++) {
        ti  = a; ti[c] = ti[c] + h
        tj  = a; tj[c] = tj[c] - h
        num = -(dwlsmv_t_Eg(P, ti, b, rho) - dwlsmv_t_Eg(P, tj, b, rho)) / (2 * h)
        alt = sum(Wa[c, .] :* (u[c, .] - u[c + 1, .]))
        e2  = max((e2, abs(alt - num)))
    }

    printf("{txt}  L2 E[g] at the truth    : %g   (must be 0)\n", nEg)
    printf("{txt}  L2 I_rho_rho vs numeric : %g\n", e1)
    printf("{txt}  L2 I_rho_a  vs numeric  : %g   (w-weighted form vs central difference)\n", e2)

    if (nEg > 1e-12) FAIL = FAIL \ ("L2: E[g] at the truth is " + strofreal(nEg))
    if (e1  > 1e-6)  FAIL = FAIL \ ("L2: I_rho_rho off by "     + strofreal(e1))
    if (e2  > 1e-6)  FAIL = FAIL \ ("L2: I_rho_a off by "       + strofreal(e2))

    // ---- L2b: dwlsmv_pairinfo at W = pi IS the form L2 just validated ----
    // L2 checks the w-weighted compact expression against numerical
    // differentiation, but it checks a copy written out in this file. What
    // ships is dwlsmv_pairinfo(), which generalises that expression to an
    // arbitrary weight matrix so the empirical and expected normalisers are
    // both reachable. If the generalisation is not EXACT at W = pi then L2 is
    // validating something the production routine no longer computes.
    dwlsmv_pairinfo(f, a, b, rho, P)
    e3 = abs(f.I_rr - Irr)
    for (c = 1; c <= rows(a); c++) {
        e3 = max((e3, abs(f.I_ra[c] - sum(Wa[c, .] :* (u[c, .] - u[c + 1, .])))))
    }
    Wm = dwlsmv_wmat(b, a, rho)
    V  = u'
    for (d = 1; d <= rows(b); d++) {
        e3 = max((e3, abs(f.I_rb[d] - sum(Wm[d, .] :* (V[d, .] - V[d + 1, .])))))
    }
    printf("{txt}  L2b pairinfo at W = pi  : %g   (must reproduce the compact form exactly)\n", e3)
    if (e3 > 1e-14) {
        FAIL = FAIL \ ("L2b: dwlsmv_pairinfo at W = pi does not reproduce the expected " +
                       "informations, off by " + strofreal(e3))
    }

    // ---- L2c: at W = pihat it IS the per-observation crossproduct --------
    // The empirical normaliser is defined as an outer product of observed
    // scores, so that definition is asserted rather than assumed: gather u_i
    // and (d_a pi / pi)_i per observation by lookup and cross the two
    // N-vectors. This is an independent route -- it never forms a table sum --
    // and it is the only check that would catch an error in the weighting that
    // ships. quadcross(), never cross(), as everywhere an O(N) accumulation is
    // asserted tightly.
    ti = dwlsmv_thresholds(Y[., 1], 4)
    tj = dwlsmv_thresholds(Y[., 2], 4)
    f  = dwls_polychoric(Y[., 1], Y[., 2], ti, tj)
    Nn = rows(Y)
    Pm = dwlsmv_cellprob(ti, tj, f.rho)
    Wm = dwlsmv_wmat(ti, tj, f.rho)
    Da = J(4, 4, 0)
    Da[1, .] =  Wm[1, .]
    Da[2, .] = -Wm[1, .]
    V  = Da :/ Pm

    uvec = J(Nn, 1, .)
    vvec = J(Nn, 1, .)
    for (i = 1; i <= Nn; i++) {
        uvec[i] = f.u[Y[i, 1] + 1, Y[i, 2] + 1]
        vvec[i] = V[Y[i, 1] + 1, Y[i, 2] + 1]
    }
    e4 = abs(f.I_rr    - quadcross(uvec, uvec) / Nn)
    e5 = abs(f.I_ra[1] - quadcross(uvec, vvec) / Nn)
    printf("{txt}  L2c pairinfo at W = pihat: I_rr %g, I_ra1 %g   (vs the literal crossproduct)\n",
           e4, e5)
    if (e4 > 1e-12 | e5 > 1e-12) {
        FAIL = FAIL \ ("L2c: the shipped information quantities are not the " +
                       "per-observation crossproduct they are defined to be (" +
                       strofreal(e4) + ", " + strofreal(e5) + ")")
    }

    // The two normalisers must DIFFER on this over-identified 4x4 table --
    // otherwise L2b and L2c are checking the same thing twice and the choice
    // recorded in the development record has quietly stopped having an effect.
    dwlsmv_pairinfo(f, ti, tj, f.rho, Pm)
    rI = abs(f.I_rr / (quadcross(uvec, uvec) / Nn) - 1)
    printf("{txt}  L2d the two normalisers  : relative gap %g on a 4x4 (expect ~1e-2)\n", rI)
    if (rI < 1e-4) {
        FAIL = FAIL \ ("L2d: expected and empirical information agree to " + strofreal(rI) +
                       " on an over-identified table -- the weight argument has stopped mattering")
    }

    // ---- L3: exact moments, and the naive form's 8% understatement -------
    psiG = dwlsmv_t_psirho(a, b, rho, 0)
    psiN = dwlsmv_t_psirho(a, b, rho, 1)
    r1   = abs(sum(P :* psiG))
    r2   = abs(sum(P :* psiN))
    vg   = sum(P :* psiG :* psiG)
    vn   = sum(P :* psiN :* psiN)
    pct  = 100 * (1 - vn / vg)

    printf("{txt}  L3 E[psi] correct %g, naive %g   -- BOTH are centred\n", r1, r2)
    printf("{txt}  L3 N*Var correct %9.6f, naive %9.6f -- naive understates by %4.1f%%\n",
           vg, vn, pct)

    if (r1 > 1e-12) FAIL = FAIL \ ("L3: correct psi is not centred, E[psi] = " + strofreal(r1))
    if (r2 > 1e-12) FAIL = FAIL \ ("L3: naive psi is not centred either -- it should be")
    if (pct < 5 | pct > 12) {
        FAIL = FAIL \ ("L3: naive understates by " + strofreal(pct) +
                       "%, expected about 8% -- the propagation term has changed size")
    }

    // ---- L5: reduction to the binary oracle (spec 4.5.3) -----------------
    worst = 0
    for (i = 1; i <= 4; i++) {
        q1 = (0.15, 0.30, 0.50, 0.75)[i]
        for (j = 1; j <= 3; j++) {
            q2 = (0.25, 0.50, 0.60)[j]
            for (c = 1; c <= 5; c++) {
                rr = (-0.6, -0.2, 0.3, 0.5, 0.8)[c]
                A  = invnormal(q1)
                B  = invnormal(q2)
                sg = sqrt(1 - rr * rr)
                p0 = binormal(A, B, rr)
                d2 = exp(-(A * A - 2 * rr * A * B + B * B) /
                         (2 * (1 - rr * rr))) / (2 * pi() * sg)
                psiG = dwlsmv_t_psirho(J(1, 1, A), J(1, 1, B), rr, 0)
                for (d = 1; d <= 2; d++) {
                    for (s = 1; s <= 2; s++) {
                        cf = (((d == 1) * (s == 1) - p0)
                              - normal((B - rr * A) / sg) * ((d == 1) - q1)
                              - normal((A - rr * B) / sg) * ((s == 1) - q2)) / d2
                        worst = max((worst, abs(psiG[d, s] - cf)))
                    }
                }
            }
        }
    }
    printf("{txt}  L5 binary oracle        : worst gap %g over 60 (q1,q2,rho) combinations\n", worst)
    if (worst > 1e-9) {
        FAIL = FAIL \ ("L5: general machinery does not reduce to the 2x2 closed form, gap " +
                       strofreal(worst))
    }

    // ---- L6/L7: rho-hat against lavaan, and the score cross-check --------
    K = (4 \ 4 \ 3 \ 3)
    m = dwlsmv_statmap_build(K)
    worst  = 0
    worstM = 0
    sc     = 0
    // Mata's & and | do not short-circuit, so the label comparison, which
    // needs four rows, is nested rather than chained.
    okM = (rows(MPL_LAB) == 4 & rows(MPL) == 4 & cols(MPL) == 4)
    if (okM) {
        okM = !any(MPL_LAB :!= ("Y1" \ "Y2" \ "Y3" \ "Y4"))
    }
    if (!okM) {
        FAIL = FAIL \ "L6b: the Mplus correlation fixture is not the 4 x 4 Y1..Y4 table"
    }
    for (i = 1; i <= 3; i++) {
        for (j = i + 1; j <= 4; j++) {
            ti = dwlsmv_thresholds(Y[., i], K[i])
            tj = dwlsmv_thresholds(Y[., j], K[j])
            f  = dwls_polychoric(Y[., i], Y[., j], ti, tj)

            if (!f.converged) {
                FAIL = FAIL \ ("L6: pair (" + strofreal(i) + "," + strofreal(j) +
                               ") did not converge")
            }

            hit = dwlsmv_t_find2(LAB, m.label[dwlsmv_pos_rho(m, i, j)])
            if (hit == 0) {
                FAIL = FAIL \ ("L6: lavaan has no entry for " +
                               m.label[dwlsmv_pos_rho(m, i, j)])
            }
            else {
                worst = max((worst, abs(f.rho - VAL[hit])))
            }

            // L6b. Mplus prints the lower triangle: pair (i,j), i < j, sits
            // in row j, column i. A missing entry is a failure, and max()
            // would skip it silently, so it is counted here.
            if (okM) {
                if (missing(MPL[j, i])) {
                    FAIL = FAIL \ ("L6b: Mplus has no correlation for pair (" +
                                   strofreal(i) + "," + strofreal(j) + ")")
                }
                else {
                    worstM = max((worstM, abs(f.rho - MPL[j, i])))
                }
            }

            // L7. The analytic score vanishes at rho-hat. Since option B the
            // root-finder uses this score, so this confirms it stopped at a
            // root; L7b, not this, checks the derivative itself.
            n  = dwlsmv_crosstab(Y[., i], Y[., j], K[i], K[j])
            sc = max((sc, abs(sum(n :* f.u) / rows(Y))))

            // Shapes, since a transposed I_rb would still be a vector.
            if (rows(f.I_ra) != K[i] - 1 | rows(f.I_rb) != K[j] - 1) {
                FAIL = FAIL \ ("L6: information vectors have the wrong shape for pair (" +
                               strofreal(i) + "," + strofreal(j) + ")")
            }
            if (f.I_rr <= 0) {
                FAIL = FAIL \ ("L6: I_rho_rho is not positive for pair (" +
                               strofreal(i) + "," + strofreal(j) + ")")
            }
        }
    }
    printf("{txt}  L6 rho-hat vs lavaan    : worst gap %g over 6 pairs\n", worst)
    printf("{txt}  L6b rho-hat vs Mplus    : worst gap %g over 6 pairs (tol 5.1e-4)\n", worstM)
    printf("{txt}  L7 analytic score at rho-hat : worst |S|/N %g\n", sc)

    // ---- L7b: the derivative check that the root-finder cannot flatter ---
    // L7 alone is no longer independent evidence that the Plackett derivative
    // is right, because dwls_polychoric() finds the root of that very score --
    // it would drive S to zero even if the derivative were wrong by a
    // constant factor. So compare the analytic score against a central
    // difference of the log-likelihood AWAY from the optimum, where both are
    // large and the comparison has content.
    ti = dwlsmv_thresholds(Y[., 1], K[1])
    tj = dwlsmv_thresholds(Y[., 2], K[2])
    n  = dwlsmv_crosstab(Y[., 1], Y[., 2], K[1], K[2])
    e1 = 0
    for (c = 1; c <= 3; c++) {
        rr  = (-0.30, 0.10, 0.70)[c]
        alt = sum(n :* (dwlsmv_dcellprob_rho(ti, tj, rr) :/
                        dwlsmv_cellprob(ti, tj, rr)))
        num = (sum(n :* ln(dwlsmv_cellprob(ti, tj, rr + h))) -
               sum(n :* ln(dwlsmv_cellprob(ti, tj, rr - h)))) / (2 * h)
        e1  = max((e1, abs(alt - num) / abs(num)))
    }
    printf("{txt}  L7b analytic score vs numeric, off-optimum : worst relative gap %g\n", e1)
    if (e1 > 1e-6) {
        FAIL = FAIL \ ("L7b: analytic score disagrees with the numerical derivative by " +
                       strofreal(e1) + " relative")
    }

    // Spec 7.2: bivariate correlations, eps_abs = eps_rel = 1e-8. Two
    // optimizers on the same data agreeing is a tighter comparison than a
    // cross-tool one, but the spec's tier is what is asserted.
    if (worst > 1e-8) {
        FAIL = FAIL \ ("L6: rho-hat differs from lavaan by " + strofreal(worst))
    }
    if (worstM > 5.1e-4) {
        FAIL = FAIL \ ("L6b: rho-hat differs from Mplus by " + strofreal(worstM))
    }
    if (sc > 1e-7) {
        FAIL = FAIL \ ("L7: analytic score does not vanish at rho-hat, |S|/N = " +
                       strofreal(sc))
    }

    // ---- L8: the saturated 2x2 identity, binary only ---------------------
    // A 2x2 table is saturated -- 3 free cells, 3 parameters -- so the MLE
    // reproduces the observed table EXACTLY (spec 4.5.3). Hence
    // Phi2(tau1, tau2; rho-hat) must equal the observed p00 to machine
    // precision, not approximately. This needs no external reference and no
    // tolerance argument, and it holds for binary and binary only: it is the
    // same over-identification boundary that stops 4.5.3's closed form
    // extending to K > 2.
    if (cols(YA) == 2) {
        ti = dwlsmv_thresholds(YA[., 1], 2)
        tj = dwlsmv_thresholds(YA[., 2], 2)
        f  = dwls_polychoric(YA[., 1], YA[., 2], ti, tj)
        n  = dwlsmv_crosstab(YA[., 1], YA[., 2], 2, 2)
        e2 = abs(binormal(ti[1], tj[1], f.rho) - n[1, 1] / rows(YA))
        printf("{txt}  L8 saturated 2x2 identity : |Phi2(tau,tau;rho-hat) - p00| = %g\n", e2)
        if (e2 > 1e-12) {
            FAIL = FAIL \ ("L8: the 2x2 MLE does not reproduce the observed table, off by " +
                           strofreal(e2))
        }
        if (!f.converged) {
            FAIL = FAIL \ "L8: the binary pair did not converge"
        }
    }

    printf("{txt}test_polychoric: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// A K1 x K2 count table back into two coded columns, row-major. The pair's
// estimate depends on the data only through its table, so this is the data.
real matrix dwlsmv_t_untable(real matrix n)
{
    real matrix Y
    real scalar r, s, k, i

    Y = J(sum(n), 2, .)
    i = 0
    for (r = 1; r <= rows(n); r++) {
        for (s = 1; s <= cols(n); s++) {
            for (k = 1; k <= n[r, s]; k++) {
                i++
                Y[i, .] = (r - 1, s - 1)
            }
        }
    }
    return(Y)
}

// ---- L9: near-duplicate pairs, rho-hat ~ .98 ------------------------------
// Pair (y1, y2) of two Monte Carlo replicates of condition mc_c_98, where
// Module 2 aborted (0007) or stalled with the centring flag set (0002) until
// the fix of @sec-emptycells in the derivation. The tables are banded: their
// off-diagonal corners are EMPTY, and their model probabilities fall below
// binormal()'s absolute resolution from rho ~ .97. Those cells contribute
// 0 * log pi = 0, so nothing production uses may depend on them.
real scalar dwlsmv_t_highrho(string scalar tag, real matrix n, real scalar rref)
{
    struct dwlsmv_pairfit scalar f
    string colvector FAIL
    real matrix Y, P0, mask
    real colvector ti, tj, uvec
    real scalar N, i, sc, e4, nbad, nmiss

    FAIL = J(0, 1, "")
    Y  = dwlsmv_t_untable(n)
    N  = rows(Y)
    ti = dwlsmv_thresholds(Y[., 1], 4)
    tj = dwlsmv_thresholds(Y[., 2], 4)
    f  = dwls_polychoric(Y[., 1], Y[., 2], ti, tj)

    sc = abs(sum(n :* f.u) / N)
    printf("{txt}  L9 %s: rho-hat %20.17f  converged %g  |S|/N %g\n", tag, f.rho, f.converged, sc)
    if (!f.converged) FAIL = FAIL \ ("L9 " + tag + ": the pair is flagged not converged")
    if (abs(f.rho) >= 0.9989) FAIL = FAIL \ ("L9 " + tag + ": rho-hat sits at the bound")
    // L7's budget: the same root-finder on the same score.
    if (sc > 1e-7 | sc >= .) {
        FAIL = FAIL \ ("L9 " + tag + ": the score does not vanish at rho-hat, |S|/N = " + strofreal(sc))
    }
    // An explicit missing count on the terms, as in L9b: sum() treats missing
    // as 0, so I_* alone come back finite even where a u is missing.
    nmiss = missing(f.u) + missing((n :/ N) :* f.u :* f.u) +
            missing(f.I_rr) + missing(f.I_ra) + missing(f.I_rb)
    printf("{txt}  L9 %s: missing values in u, the I_rr terms and I_*: %g\n", tag, nmiss)
    if (nmiss) {
        FAIL = FAIL \ ("L9 " + tag + ": a missing value in u or an information quantity")
    }
    // The empty-cell rule on the cells it exists for: EMPTY and underflowed at
    // rho-hat. There u must be exactly 0, not Dr/1 -- the sums cannot see the
    // difference (those cells carry no weight), only the shipped u can. Not
    // every table has such a cell at rho-hat, so their count is accumulated in
    // L9_nunder and asserted positive over the tables, below.
    P0   = dwlsmv_cellprob_raw(ti, tj, f.rho)
    mask = (n :== 0) :& (P0 :<= 0)
    nbad = sum(mask :& (f.u :!= 0))
    st_numscalar("L9_nunder", st_numscalar("L9_nunder") + sum(mask))
    printf("{txt}  L9 %s: empty underflowed cells %g, u nonzero on %g of them\n", tag, sum(mask), nbad)
    if (nbad) FAIL = FAIL \ ("L9 " + tag + ": u is not exactly 0 on an empty underflowed cell")
    // Spec 7.2: bivariate correlations, eps_abs = eps_rel = 1e-8, as L6.
    // The reference is tier 1, base R's root of the observed-cells score.
    printf("{txt}  L9 %s: rho-hat vs tier-1 reference %g (1e-8 tier)\n", tag, abs(f.rho - rref))
    if (abs(f.rho - rref) > 1e-8 | rref >= .) {
        FAIL = FAIL \ ("L9 " + tag + ": rho-hat differs from the tier-1 reference by " +
                       strofreal(abs(f.rho - rref)))
    }

    // L2c's check where it matters: the shipped informations are still the
    // per-observation crossproducts when empty cells underflow.
    uvec = J(N, 1, .)
    for (i = 1; i <= N; i++) {
        uvec[i] = f.u[Y[i, 1] + 1, Y[i, 2] + 1]
    }
    e4 = abs(f.I_rr - quadcross(uvec, uvec) / N)
    printf("{txt}  L9 %s: I_rr vs the literal crossproduct %g\n", tag, e4)
    if (e4 > 1e-12 | e4 >= .) {
        FAIL = FAIL \ ("L9 " + tag + ": I_rr is not the per-observation crossproduct (" + strofreal(e4) + ")")
    }

    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// ---- L9b: the empty-cell rule, on a cell that is ZERO under any rounding --
// Two equal thresholds make the middle row of cells exactly 0, which is the
// deterministic stand-in for an empty cell that underflows. With that row
// empty in the data, the point is usable (option B's rule R_P, over observed
// cells only) and pairinfo zeroes u there; with a count in it, the point is
// unusable, and pairinfo still aborts -- the L9c error path below. Until
// 2026-09-30 this checked optimize()'s log-likelihood evaluator, which option
// B removed; the rule it pinned now lives in dwlsmv_rf_eval()'s usability.
real scalar dwlsmv_t_emptyrule()
{
    struct dwlsmv_pairfit scalar f
    string colvector FAIL
    real matrix n, W
    real rowvector e1, e2
    real scalar i, nmiss

    FAIL = J(0, 1, "")
    n = (10, 5 \ 0, 0 \ 7, 12)
    e1 = dwlsmv_rf_eval(0.3, n, (0 \ 0), (0.5))
    n[2, 1] = 1
    e2 = dwlsmv_rf_eval(0.3, n, (0 \ 0), (0.5))
    printf("{txt}  L9b usable at rho .3, zero cell empty %g, zero cell observed %g (must be 0)\n", e1[2], e2[2])
    if (e1[2] != 1 | missing(e1)) FAIL = FAIL \ "L9b: an EMPTY zero-probability cell makes the point unusable"
    if (e2[2] != 0) FAIL = FAIL \ "L9b: an OBSERVED zero-probability cell does not make the point unusable"

    W = (10, 5 \ 0, 0 \ 7, 12) :/ 34
    dwlsmv_pairinfo(f, (0 \ 0), (0.5), 0.3, W)
    printf("{txt}  L9b pairinfo on the empty zero row: u there %g %g, I_rr %g\n", f.u[2, 1], f.u[2, 2], f.I_rr)
    if (f.u[2, 1] != 0 | f.u[2, 2] != 0) FAIL = FAIL \ "L9b: u is not zeroed on the empty zero-probability row"
    // An explicit missing count, taken on the terms rather than the totals:
    // sum() treats missing as 0, so I_rr = sum(W :* u :* u) comes back finite
    // even where a u is missing, and a check on I_* alone cannot fire. A
    // missing u makes its term 0 * . = . and is counted here.
    nmiss = missing(f.u) + missing(W :* f.u :* f.u) +
            missing(f.I_rr) + missing(f.I_ra) + missing(f.I_rb)
    printf("{txt}  L9b missing values in u, the I_rr terms and I_*: %g\n", nmiss)
    if (nmiss) {
        FAIL = FAIL \ "L9b: a missing value in u or an information, although the zero row carries no weight"
    }
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// ---- L9d: a boundary 2x2 pair is flagged -------------------------------------
// An empty cell puts the MLE at |rho| = 1 (derivation @sec-emptycells). Under
// option B (since 2026-09-30) the score keeps its sign up to B, so the pair is
// status 2, no sign change up to the bound, with rho-hat = 0.999 exactly. What
// Stage 1 does with it (abort; empty-cell pairs are never kept) is decided 2026-10-03 and
// tested in test_sample_stats.do, which also pins (0.999, 2, 11). The pi
// and missing counts restate pairinfo's own abort, which the capture below
// would count; they are belts, not independent evidence.
real scalar dwlsmv_t_boundary()
{
    struct dwlsmv_pairfit scalar f
    string colvector FAIL
    real matrix n, Y, P
    real colvector ti, tj
    real scalar nneg, nmiss, i

    FAIL = J(0, 1, "")
    n  = (30, 2 \ 0, 8)
    Y  = dwlsmv_t_untable(n)
    ti = dwlsmv_thresholds(Y[., 1], 2)
    tj = dwlsmv_thresholds(Y[., 2], 2)
    f  = dwls_polychoric(Y[., 1], Y[., 2], ti, tj)

    P     = dwlsmv_cellprob_raw(ti, tj, f.rho)
    nneg  = sum((n :> 0) :& (P :<= 0))
    nmiss = missing(f.u) + missing(f.I_rr) + missing(f.I_ra) + missing(f.I_rb)
    printf("{txt}  L9d rho-hat %21x  converged %g  observed cells with pi <= 0: %g  missing: %g\n",
           f.rho, f.converged, nneg, nmiss)
    if (f.converged != 0) FAIL = FAIL \ "L9d: a boundary pair is not flagged"
    if (f.status != 2) FAIL = FAIL \ "L9d: a boundary pair is not status 2"
    if (abs(f.rho) > 0.999 | f.rho >= .) FAIL = FAIL \ "L9d: rho-hat lies beyond the 0.999 bound"
    if (nneg) FAIL = FAIL \ "L9d: an observed cell is not positive at rho-hat"
    if (nmiss) FAIL = FAIL \ "L9d: u or an information has a missing value"
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// L9c's abort, as a function because it touches a struct.
void dwlsmv_t_pairinfo_obs0()
{
    struct dwlsmv_pairfit scalar f

    dwlsmv_pairinfo(f, (0 \ 0), (0.5), 0.3, (10, 5 \ 1, 0 \ 7, 12) :/ 35)
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_polychoric(YD, LAB_LAV, VAL_LAV, YA, MPL_LAB, MPL)))

* L9 runs under capture: until the fix, 0007 ABORTS inside Module 2, and an
* abort must be counted as a failure rather than end the test unreported.
* The tables and their tier-1 rho-hat come from one fixture, written by
* fixtures/generate/10_fixtures.R (derivation @sec-emptycells, Level E7).
import delimited using "fixtures/expected/highrho_polychoric.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: HR_TAG = st_sdata(., "table")
mata: HR_N   = strtoreal(st_sdata(., ("n_0_0", "n_0_1", "n_0_2", "n_0_3",   ///
                                      "n_1_0", "n_1_1", "n_1_2", "n_1_3",   ///
                                      "n_2_0", "n_2_1", "n_2_2", "n_2_3",   ///
                                      "n_3_0", "n_3_1", "n_3_2", "n_3_3")))
mata: HR_RHO = strtoreal(st_sdata(., "rho_hat"))
local nhigh 0
scalar L9_nunder = 0
local nhr = _N
if `nhr' != 2 {
    display as error "test_polychoric: L9 expects 2 tables in highrho_polychoric.csv, found `nhr'"
    local nhigh = `nhigh' + 1
}
forvalues r = 1/`nhr' {
    capture noisily mata: st_local("k", strofreal(dwlsmv_t_highrho(HR_TAG[`r'],     ///
        rowshape(HR_N[`r', .], 4), HR_RHO[`r'])))
    if _rc {
        display as error "test_polychoric: L9 table `r' aborted with rc " _rc
        local nhigh = `nhigh' + 1
    }
    else local nhigh = `nhigh' + `k'
}
* So the u rule cannot pass on no cells at all.
if scalar(L9_nunder) == 0 {
    display as error "test_polychoric: L9 found no empty cell underflowing at rho-hat; the u rule is untested"
    local nhigh = `nhigh' + 1
}
capture noisily mata: st_local("k", strofreal(dwlsmv_t_emptyrule()))
if _rc {
    display as error "test_polychoric: L9b aborted with rc " _rc
    local nhigh = `nhigh' + 1
}
else local nhigh = `nhigh' + `k'
capture noisily mata: st_local("k", strofreal(dwlsmv_t_boundary()))
if _rc {
    display as error "test_polychoric: L9d aborted with rc " _rc
    local nhigh = `nhigh' + 1
}
else local nhigh = `nhigh' + `k'
local nfail = `nfail' + `nhigh'

* ---- error paths -----------------------------------------------------------
* Each asserts rc 3498, not merely a nonzero rc. An _error() text over 100
* characters raises rc 3300 and discards the message, and a `!_rc` check
* passes on that as happily as on the intended abort. That is how the
* non-positive-cell path below shipped message-less until 2026-09-24.
local nerr 0

capture mata: dwlsmv_cellprob((0 \ 1), (0.5), 1)
if _rc != 3498 {
    display as error "test_polychoric: cellprob on rho = 1 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_cellprob((0 \ 1), (0.5), -1)
if _rc != 3498 {
    display as error "test_polychoric: cellprob on rho = -1 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
* The same bound on dwlsmv_cellprob_raw() itself. Through cellprob() above, a
* raw guard weakened to rho > 1 is masked by cellprob's own min(P) <= 0 check;
* the raw routine has no such backstop, and the root-finder and pairinfo call it.
capture mata: dwlsmv_cellprob_raw((0 \ 1), (0.5), 1)
if _rc != 3498 {
    display as error "test_polychoric: cellprob_raw on rho = 1 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_cellprob_raw((0 \ 1), (0.5), -1)
if _rc != 3498 {
    display as error "test_polychoric: cellprob_raw on rho = -1 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
* Two equal thresholds: the middle row of cells is 0 exactly, because its two
* corner rows are the same numbers -- so this reaches min(P) <= 0 under any
* rounding, not by luck.
capture mata: dwlsmv_cellprob((0 \ 0), (0.5), 0.3)
if _rc != 3498 {
    display as error "test_polychoric: cellprob on a zero-probability cell gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
* L9c. The empty-cell rule must not reach an OBSERVED cell: a zero-probability
* cell that carries weight still aborts in pairinfo.
capture mata: dwlsmv_t_pairinfo_obs0()
if _rc != 3498 {
    display as error "test_polychoric: pairinfo on an observed zero-probability cell gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_wmat((0 \ 1), J(0, 1, .), 0.3)
if _rc != 3498 {
    display as error "test_polychoric: wmat on an empty threshold vector gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwls_polychoric((0 \ 1 \ 1), (0 \ 1), (0), (0))
if _rc != 3498 {
    display as error "test_polychoric: polychoric on mismatched lengths gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_polychoric: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_polychoric: Module 2 verified -- Phi2 identity, derivatives, informations, moments, binary oracle, lavaan, score"
exit 0
