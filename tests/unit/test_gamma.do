*! Module 4: the full Gamma -- structural invariants, the binary oracle, exact
*! moment identities, the centring diagnostic, and lavaan.
*!
*! Executed under StataMP 18 on 2026-09-20 and it passes. Checks 5c and 5d
*! (case E) were added and executed 2026-09-25.
*!
*! tests/unit/test_gamma_tau.do covers the threshold block, which is the third
*! of Gamma that needs no Phi2 and no rho. This file covers everything else.
*!
*! Six checks, ordered cheapest and most decisive first, following spec Section
*! 5's validation ladder. The lavaan comparison is LAST and is the least
*! diagnostic of them, because it can agree for the wrong reason.
*!
*!   1. STRUCTURAL INVARIANTS. Symmetry and positive definiteness, on Case D
*!      (p > 2, unequal category counts). Free, and they catch a large share of
*!      assembly errors before any reference is consulted.
*!
*!   2. THE BINARY ORACLE (spec 4.5.3), sub-phase 4a's stated acceptance test.
*!      The general machinery must reproduce the independently derived 2x2
*!      closed form across the parameter space, not at one point. This runs it
*!      through the PRODUCTION routine, dwlsmv_psi_rho_table(), rather than
*!      through a reconstruction -- test_polychoric.do already checks the
*!      mathematics with a test-local copy, and what is unverified until here
*!      is the assembly itself.
*!
*!   3. EXACT MOMENT IDENTITIES, the correlation-block counterpart of the
*!      threshold block's identity. psi_rho is a cell lookup, so wherever two
*!      statistics depend on the SAME pair of variables their Gamma entry is
*!      exactly a sum over that pair's K1 x K2 table -- no sampling error, no
*!      tolerance to argue about. That covers every correlation variance and
*!      every threshold-by-correlation entry whose threshold belongs to the
*!      pair. It is the only check here that would catch an indexing error in
*!      dwlsmv_psi_matrix()'s correlation gather. The entries it CANNOT cover
*!      are the genuinely joint ones -- two pairs sharing one variable need a
*!      trivariate normal, sharing none a quadrivariate -- and those have no
*!      table route by construction. That is exactly why production builds
*!      Gamma as an empirical outer product (spec 4.5.4).
*!
*!   4. CENTRING (spec 4.5.5), SPLIT BY BLOCK, plus a route agreement. The
*!      threshold entries are exactly zero because q-hat is literally the
*!      sample proportion; the correlation entries are zero only to Module 2's
*!      convergence, so a single blanket tolerance would be set by the looser
*!      half and would hide a real bug in the tighter one. Separately,
*!      dwlsmv_psibar() computes from cell COUNTS while Gamma is built from the
*!      GATHERED matrix, so the two are independent routes and are compared as
*!      such -- that is what makes psibar able to say something about the
*!      gather, which on its own it cannot.
*!
*!      What none of this can detect: a MISSING PROPAGATION TERM. The naive
*!      influence function is also exactly centred. Check 2 and check 6 are
*!      what stand between us and that.
*!
*!   5. LAVAAN, DECOMPOSED INTO TWO COMPARISONS, because they have different
*!      achievable agreement and merging them would set the tolerance by the
*!      looser one:
*!
*!        5a. at LAVAAN'S OWN rho-hat, asserted at 1e-12. This is the question
*!            "is our Gamma arithmetic lavaan's arithmetic", and the answer is
*!            yes to 1.4e-13 -- the same agreement the threshold block reaches.
*!        5b. end to end, at OUR rho-hat, asserted at spec 7.2's BIVARIATE
*!            CORRELATION tier of 1e-8 rather than its Gamma tier of 1e-10.
*!            That is not a loosened tolerance, it is the right one: the
*!            comparison is capped by how closely the two rho-hats agree, and
*!            rho-hat has its own tier. Measured 2026-09-20: the rho-hats
*!            differ by 2.6e-10 and Gamma by 9.1e-10, and substituting lavaan's
*!            rho collapses the latter to 1.4e-13, which is what identifies the
*!            cause rather than assuming it. The direction is worth recording:
*!            our score vanishes at rho-hat to 4.8e-16 and lavaan's to 2.6e-10,
*!            so lavaan is the less converged of the two and tightening our
*!            optimizer cannot close the gap.
*!        5c. CASE E at LAVAAN'S OWN rho-hat, at 5a's 1e-12 abs and rel, from
*!            the same spec 7.2 row ("at a matched rho-hat"). It asks 5a's
*!            question on six variables, fifteen pairs and two factors, so the
*!            gather's indexing is tested beyond case D's four. lavaan's rho-hat
*!            comes from lavaan_case_E_sampstat.csv (added 2026-09-25). Stata
*!            reads ratio 0.069 (StataMP 18, 2026-09-25). R read 0.11 in the
*!            same construction (lavaan 0.7.2), set by the threshold block
*!            (2.8e-13 on y5|t1's variance), which rho-hat cannot touch. The
*!            two differ in last-digit arithmetic, and both sit well inside 1.
*!        5d. CASE E end to end, PRINTED, NOT ASSERTED (decided 2026-09-25).
*!            It reads about 3.0 at the rho tier (6.0e-8 on y5~~y6's variance).
*!            5c identifies the cause: lavaan's y5~~y6 rho-hat is 1.67e-8 from
*!            the root, and a first-order propagation of the rho gaps reproduces
*!            the Gamma gap to 1e-13 (R, one-off). lavaan's case E rho-hat
*!            itself misses spec 7.2's bivariate-correlation tier (ratio 1.12),
*!            so no Gamma tier inherited from it can hold. Asserting it would
*!            mean loosening a tolerance to make it pass.
*!
*!   6. THE NORMALISER, as a regression guard on a decision rather than a
*!      numerical claim. dwlsmv_pairinfo() takes the weight matrix explicitly,
*!      so both the shipped empirical form and the expected form of the Gate 3
*!      result are reachable. Rebuilding Gamma the other way must move it FAR
*!      from lavaan. If a later change quietly switches the normaliser back --
*!      the kind of thing that looks like a simplification -- this is what
*!      fails.
*!
*! Reference CSVs are read as text and converted with strtoreal(), and with an
*! EXPLICIT delimiters(","). Both matter and both were found the hard way. A
*! float-typed column loses the low 29 mantissa bits, a relative error near
*! 4e-08 against the 1e-12 asserted here. And `import delimited` AUTO-DETECTS
*! the delimiter: lavaan_case_D_gamma.csv has ten "|" characters in its header,
*! and on a fixture where pipes are as frequent as commas Stata splits on the
*! pipe, silently turning the label "y1|t1" into "y1". Nothing warns; the
*! import succeeds and every label lookup then misses.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_catdegen.mata          ///
             src/mata/dwlsmv_cumprop.mata           ///
             src/mata/dwlsmv_thresholds.mata        ///
             src/mata/dwlsmv_psi_tau_table.mata     ///
             src/mata/dwlsmv_cellprob.mata          ///
             src/mata/dwlsmv_dcellprob_rho.mata     ///
             src/mata/dwlsmv_wmat.mata              ///
             src/mata/dwlsmv_crosstab.mata          ///
             src/mata/dwlsmv_pairinfo.mata          ///
             src/mata/dwls_polychoric.mata          ///
             src/mata/dwlsmv_koff.mata              ///
             src/mata/dwlsmv_psi_tau_stack.mata     ///
             src/mata/dwlsmv_psi_rho_table.mata     ///
             src/mata/dwlsmv_psi_matrix.mata        ///
             src/mata/dwlsmv_psibar.mata            ///
             src/mata/dwls_gamma.mata               ///
             src/mata/dwlsmv_emptycell1.mata      ///
             src/mata/dwlsmv_emptycell_repair.mata ///
             src/mata/dwlsmv_emptycell_pairs.mata ///
             src/mata/dwls_sample_stats.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_gamma: `f' failed to load"
        exit 9
    }
}

* Module 4 exists, so this test never legitimately skips: it passes or fails.

import delimited using "fixtures/expected/lavaan_case_D_gamma.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local lcols `r(varlist)'
mata: G_LAV   = strtoreal(st_sdata(., "`lcols'"))
mata: LAB_LAV = st_sdata(., "label")

import delimited using "fixtures/expected/lavaan_case_D_sampstat.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local scols `r(varlist)'
mata: R_LAV = strtoreal(st_sdata(., "`scols'"))

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

import delimited using "fixtures/data/case_A.csv", clear varnames(1) delimiters(",")
mata: YA = st_data(., .)

import delimited using "fixtures/expected/lavaan_case_E_gamma.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local ecols `r(varlist)'
mata: G_LAVE   = strtoreal(st_sdata(., "`ecols'"))
mata: LAB_LAVE = st_sdata(., "label")

import delimited using "fixtures/expected/lavaan_case_E_sampstat.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local escols `r(varlist)'
mata: R_LAVE = strtoreal(st_sdata(., "`escols'"))

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

// Gamma rebuilt with a chosen normaliser and, optionally, chosen correlations.
// Used by checks 5a and 6. It duplicates dwls_sample_stats()'s loop rather
// than calling it, because the whole point is to vary something that routine
// deliberately does not expose.
real matrix dwlsmv_g_rebuild(real matrix Y, struct dwlsmv_statmap scalar m,
                             real colvector taus, real matrix psi_tau,
                             real scalar use_empirical, real matrix RHO)
{
    struct dwlsmv_pairfit scalar f
    struct dwlsmv_pairfit colvector pairs
    real matrix nab
    real colvector ti, tj
    real scalar k, pos, vi, vj

    pairs = dwlsmv_pairfit(m.P, 1)

    for (k = 1; k <= m.P; k++) {
        pos = m.T + k
        vi  = m.var_i[pos]
        vj  = m.var_j[pos]
        ti  = taus[|m.Coff[vi] + 1 \ m.Coff[vi] + m.C[vi]|]
        tj  = taus[|m.Coff[vj] + 1 \ m.Coff[vj] + m.C[vj]|]
        nab = dwlsmv_crosstab(Y[., vi], Y[., vj], m.K[vi], m.K[vj])

        f = dwls_polychoric(Y[., vi], Y[., vj], ti, tj)
        if (rows(RHO) > 1) {
            f.rho = RHO[vi, vj]
        }
        if (use_empirical) {
            dwlsmv_pairinfo(f, ti, tj, f.rho, nab :/ rows(Y))
        }
        else {
            dwlsmv_pairinfo(f, ti, tj, f.rho, dwlsmv_cellprob(ti, tj, f.rho))
        }
        f.i      = vi
        f.j      = vj
        pairs[k] = f
    }

    return(dwls_gamma(Y, m, pairs, psi_tau))
}

real scalar dwlsmv_g_find3(string colvector pool, string scalar want)
{
    real scalar i

    for (i = 1; i <= rows(pool); i++) {
        if (pool[i] == want) {
            return(i)
        }
    }
    return(0)
}

real scalar dwlsmv_g_full(real matrix Y, real colvector K, real matrix YA,
                          real matrix G_LAV, string colvector LAB_LAV,
                          real matrix R_LAV)
{
    struct dwlsmv_stage1 scalar st, sa
    struct dwlsmv_statmap scalar m
    struct dwlsmv_pairfit scalar f
    struct dwlsmv_pairfit colvector pairs
    string colvector FAIL
    real matrix psi_tau, Psi, tbl, nab, pa, pb, G2, G3, Gex
    real colvector taus, tau, ti, tj, bar, barg
    real scalar N, T, D, j, c, k, i, s, r, pos, vi, vj
    real scalar sym, eig, worst, q1, q2, rr, A, B, sg, p0, d2, cf
    real scalar r_diag, r_cross, cenT, cenR, r_route, dRho, sOur, sLav
    real scalar r_lavA, r_lavB, gap_emp, gap_exp

    FAIL = J(0, 1, "")
    N    = rows(Y)

    st = dwls_sample_stats(Y, K)
    m  = st.map
    T  = m.T
    D  = m.D

    // ---- 1. structural invariants -----------------------------------------
    sym = dwlsmv_symmetry(st.Gamma)
    eig = dwlsmv_mineig(st.Gamma)
    printf("{txt}  1. structural            : symmetry %g, min eigenvalue %g\n", sym, eig)
    if (sym > 1e-15) {
        FAIL = FAIL \ ("structural: Gamma is not symmetric, worst |G - G'| = " +
                       strofreal(sym))
    }
    if (eig <= 0) {
        FAIL = FAIL \ ("structural: Gamma is not positive definite, min eigenvalue " +
                       strofreal(eig))
    }

    // ---- 2. the binary oracle (spec 4.5.3), through the production routine -
    // 60 combinations, as Gate 3's acceptance test specifies. Binary tables are
    // saturated, so the normaliser is irrelevant here and the expected form is
    // passed -- there is no sample at these hypothetical parameter values.
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

                dwlsmv_pairinfo(f, J(1, 1, A), J(1, 1, B), rr,
                                dwlsmv_cellprob(J(1, 1, A), J(1, 1, B), rr))
                pa  = dwlsmv_psi_tau_table(J(1, 1, q1))
                pb  = dwlsmv_psi_tau_table(J(1, 1, q2))
                tbl = dwlsmv_psi_rho_table(f, pa, pb)

                for (r = 1; r <= 2; r++) {
                    for (s = 1; s <= 2; s++) {
                        cf = (((r == 1) * (s == 1) - p0)
                              - normal((B - rr * A) / sg) * ((r == 1) - q1)
                              - normal((A - rr * B) / sg) * ((s == 1) - q2)) / d2
                        worst = max((worst, abs(tbl[r, s] - cf)))
                    }
                }
            }
        }
    }
    printf("{txt}  2. binary oracle (4.5.3) : worst gap %g over 60 (q1,q2,rho) combinations\n",
           worst)
    if (worst > 1e-9) {
        FAIL = FAIL \ ("binary oracle: dwlsmv_psi_rho_table does not reduce to the 2x2 " +
                       "closed form, gap " + strofreal(worst))
    }

    // ---- rebuild the ingredients, for checks 3 onward ----------------------
    taus = J(T, 1, .)
    for (j = 1; j <= m.p; j++) {
        tau = dwlsmv_thresholds(Y[., j], m.K[j])
        for (c = 1; c <= m.C[j]; c++) {
            taus[dwlsmv_pos_tau(m, j, c)] = tau[c]
        }
    }
    psi_tau = dwlsmv_psi_tau_stack(Y, m)
    pairs   = dwlsmv_pairfit(m.P, 1)
    for (k = 1; k <= m.P; k++) {
        pos = m.T + k
        vi  = m.var_i[pos]
        vj  = m.var_j[pos]
        ti  = taus[|m.Coff[vi] + 1 \ m.Coff[vi] + m.C[vi]|]
        tj  = taus[|m.Coff[vj] + 1 \ m.Coff[vj] + m.C[vj]|]
        f   = dwls_polychoric(Y[., vi], Y[., vj], ti, tj)
        f.i = vi
        f.j = vj
        pairs[k] = f
    }
    Psi = dwlsmv_psi_matrix(Y, m, pairs, psi_tau)

    // ---- 3. exact moment identities in the correlation block ---------------
    r_diag  = 0
    r_cross = 0
    for (k = 1; k <= m.P; k++) {
        pos = m.T + k
        vi  = m.var_i[pos]
        vj  = m.var_j[pos]
        pa  = psi_tau[|dwlsmv_koff(m)[vi] + 1, m.Coff[vi] + 1 \
                       dwlsmv_koff(m)[vi] + m.K[vi], m.Coff[vi] + m.C[vi]|]
        pb  = psi_tau[|dwlsmv_koff(m)[vj] + 1, m.Coff[vj] + 1 \
                       dwlsmv_koff(m)[vj] + m.K[vj], m.Coff[vj] + m.C[vj]|]
        tbl = dwlsmv_psi_rho_table(pairs[k], pa, pb)
        nab = dwlsmv_crosstab(Y[., vi], Y[., vj], m.K[vi], m.K[vj]) :/ N

        // the correlation variance
        r_diag = max((r_diag, dwlsmv_gapratio(st.Gamma[pos, pos],
                                              sum(nab :* tbl :* tbl), 1e-12, 0)))

        // threshold-by-correlation, for thresholds belonging to this pair
        for (c = 1; c <= m.C[vi]; c++) {
            s       = dwlsmv_pos_tau(m, vi, c)
            r_cross = max((r_cross,
                           dwlsmv_gapratio(st.Gamma[s, pos],
                                           sum(nab :* (pa[., c] * J(1, m.K[vj], 1)) :* tbl),
                                           1e-12, 0)))
        }
        for (c = 1; c <= m.C[vj]; c++) {
            s       = dwlsmv_pos_tau(m, vj, c)
            r_cross = max((r_cross,
                           dwlsmv_gapratio(st.Gamma[s, pos],
                                           sum(nab :* (J(m.K[vi], 1, 1) * pb[., c]') :* tbl),
                                           1e-12, 0)))
        }
    }
    printf("{txt}  3. exact moments         : correlation variances %g, cross entries %g (1e-12 abs)\n",
           r_diag, r_cross)
    if (r_diag > 1) {
        FAIL = FAIL \ ("moments: a correlation variance disagrees with its table sum, ratio " +
                       strofreal(r_diag))
    }
    if (r_cross > 1) {
        FAIL = FAIL \ ("moments: a threshold-by-correlation entry disagrees with its table " +
                       "sum, ratio " + strofreal(r_cross))
    }

    // ---- 4. centring, split by block, and the two routes -------------------
    bar  = st.psibar
    barg = J(D, 1, 0)
    for (s = 1; s <= D; s++) {
        barg[s] = sum(Psi[., s]) / N
    }
    cenT    = max(abs(bar[|1 \ T|]))
    cenR    = max(abs(bar[|T + 1 \ D|]))
    r_route = max(abs(bar - barg))
    printf("{txt}  4. centring              : thresholds %g, correlations %g, count-vs-gather %g\n",
           cenT, cenR, r_route)
    if (cenT > 1e-14) {
        FAIL = FAIL \ ("centring: threshold entries are not exactly centred, worst " +
                       strofreal(cenT) + " -- q-hat is the sample proportion, so this is exact")
    }
    if (cenR > 1e-10) {
        FAIL = FAIL \ ("centring: correlation entries are off by " + strofreal(cenR) +
                       ", beyond Module 2's convergence tolerance")
    }
    if (r_route > 1e-12) {
        FAIL = FAIL \ ("centring: the count route and the gather route disagree by " +
                       strofreal(r_route) + " -- an indexing error in the gather")
    }

    // ---- 5. lavaan, decomposed --------------------------------------------
    dRho = 0
    sOur = 0
    sLav = 0
    for (k = 1; k <= m.P; k++) {
        pos = m.T + k
        vi  = m.var_i[pos]
        vj  = m.var_j[pos]
        ti  = taus[|m.Coff[vi] + 1 \ m.Coff[vi] + m.C[vi]|]
        tj  = taus[|m.Coff[vj] + 1 \ m.Coff[vj] + m.C[vj]|]
        nab = dwlsmv_crosstab(Y[., vi], Y[., vj], m.K[vi], m.K[vj])
        dRho = max((dRho, abs(pairs[k].rho - R_LAV[vi, vj])))
        sOur = max((sOur, abs(sum(nab :* pairs[k].u) / N)))
        f    = pairs[k]
        dwlsmv_pairinfo(f, ti, tj, R_LAV[vi, vj], nab :/ N)
        sLav = max((sLav, abs(sum(nab :* f.u) / N)))
    }
    printf("{txt}  5. rho-hat vs lavaan     : gap %g; |S/N| ours %g, lavaan's %g\n",
           dRho, sOur, sLav)

    G2     = dwlsmv_g_rebuild(Y, m, taus, psi_tau, 1, R_LAV)
    r_lavA = 0
    r_lavB = 0
    for (s = 1; s <= D; s++) {
        i = dwlsmv_g_find3(LAB_LAV, m.label[s])
        if (i == 0) {
            FAIL = FAIL \ ("lavaan: no entry labelled " + m.label[s])
            continue
        }
        for (r = 1; r <= D; r++) {
            j = dwlsmv_g_find3(LAB_LAV, m.label[r])
            if (j == 0) {
                continue
            }
            r_lavA = max((r_lavA, dwlsmv_gapratio(G2[s, r], G_LAV[i, j], 1e-12, 1e-12)))
            r_lavB = max((r_lavB,
                          dwlsmv_gapratio(st.Gamma[s, r], G_LAV[i, j], 1e-8, 1e-8)))
        }
    }
    printf("{txt}  5a. Gamma at lavaan's rho: ratio %g (1e-12 abs and rel)\n", r_lavA)
    printf("{txt}  5b. Gamma at our rho     : ratio %g (1e-8 abs and rel, the rho tier)\n",
           r_lavB)
    if (r_lavA > 1) {
        FAIL = FAIL \ ("lavaan 5a: at lavaan's own rho the two constructions still " +
                       "disagree, ratio " + strofreal(r_lavA) +
                       " -- this is an arithmetic difference, not a convergence one")
    }
    if (r_lavB > 1) {
        FAIL = FAIL \ ("lavaan 5b: end-to-end ratio " + strofreal(r_lavB))
    }

    // ---- 6. the normaliser, as a regression guard on a decision ------------
    Gex     = dwlsmv_g_rebuild(Y, m, taus, psi_tau, 0, J(0, 0, .))
    gap_emp = max(abs(st.Gamma - G_LAV))
    gap_exp = max(abs(Gex - G_LAV))
    printf("{txt}  6. normaliser guard      : empirical gap %g, expected gap %g (ratio %g)\n",
           gap_emp, gap_exp, gap_exp / gap_emp)
    if (gap_exp <= 100 * gap_emp) {
        FAIL = FAIL \ ("normaliser: the expected-information Gamma is not clearly further " +
                       "from lavaan than the empirical one -- either the weight argument " +
                       "stopped mattering or the shipped normaliser changed")
    }

    printf("{txt}test_gamma: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }

    return(rows(FAIL))
}

// Case E: 5c asserts Gamma at lavaan's own rho-hat at 5a's tier; 5d prints
// the end-to-end ratio at the rho tier and does not assert it (see header).
real scalar dwlsmv_g_caseE(real matrix Y, real colvector K,
                           real matrix G_LAV, string colvector LAB_LAV,
                           real matrix R_LAV)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_statmap scalar m
    struct dwlsmv_pairfit scalar f
    string colvector FAIL
    real matrix psi_tau, nab, G2
    real colvector taus, tau, ti, tj
    real scalar N, j, c, k, s, r, i, pos, vi, vj
    real scalar dRho, sOur, sLav, r_lavC, r_lavD

    FAIL = J(0, 1, "")
    N    = rows(Y)
    st   = dwls_sample_stats(Y, K)
    m    = st.map

    // Mata's max() skips missing values, so a failed strtoreal() or a short
    // fixture would otherwise pass silently.
    if (rows(G_LAV) != m.D | cols(G_LAV) != m.D | missing(G_LAV)) {
        FAIL = FAIL \ "case E: lavaan Gamma fixture is not a complete D x D matrix"
    }
    if (rows(R_LAV) != m.p | cols(R_LAV) != m.p | missing(R_LAV)) {
        FAIL = FAIL \ "case E: lavaan sampstat fixture is not a complete p x p matrix"
    }
    if (rows(FAIL)) {
        for (j = 1; j <= rows(FAIL); j++) {
            printf("{err}    %s\n", FAIL[j])
        }
        return(rows(FAIL))
    }

    taus = J(m.T, 1, .)
    for (j = 1; j <= m.p; j++) {
        tau = dwlsmv_thresholds(Y[., j], m.K[j])
        for (c = 1; c <= m.C[j]; c++) {
            taus[dwlsmv_pos_tau(m, j, c)] = tau[c]
        }
    }
    psi_tau = dwlsmv_psi_tau_stack(Y, m)

    dRho = 0
    sOur = 0
    sLav = 0
    for (k = 1; k <= m.P; k++) {
        pos  = m.T + k
        vi   = m.var_i[pos]
        vj   = m.var_j[pos]
        ti   = taus[|m.Coff[vi] + 1 \ m.Coff[vi] + m.C[vi]|]
        tj   = taus[|m.Coff[vj] + 1 \ m.Coff[vj] + m.C[vj]|]
        nab  = dwlsmv_crosstab(Y[., vi], Y[., vj], m.K[vi], m.K[vj])
        f    = dwls_polychoric(Y[., vi], Y[., vj], ti, tj)
        dRho = max((dRho, abs(f.rho - R_LAV[vi, vj])))
        sOur = max((sOur, abs(sum(nab :* f.u) / N)))
        dwlsmv_pairinfo(f, ti, tj, R_LAV[vi, vj], nab :/ N)
        sLav = max((sLav, abs(sum(nab :* f.u) / N)))
    }
    printf("{txt}  5c. case E rho vs lavaan : gap %g; |S/N| ours %g, lavaan's %g\n",
           dRho, sOur, sLav)

    G2     = dwlsmv_g_rebuild(Y, m, taus, psi_tau, 1, R_LAV)
    r_lavC = 0
    r_lavD = 0
    for (s = 1; s <= m.D; s++) {
        i = dwlsmv_g_find3(LAB_LAV, m.label[s])
        if (i == 0) {
            FAIL = FAIL \ ("case E lavaan: no entry labelled " + m.label[s])
            continue
        }
        for (r = 1; r <= m.D; r++) {
            j = dwlsmv_g_find3(LAB_LAV, m.label[r])
            if (j == 0) {
                continue
            }
            r_lavC = max((r_lavC, dwlsmv_gapratio(G2[s, r], G_LAV[i, j], 1e-12, 1e-12)))
            r_lavD = max((r_lavD,
                          dwlsmv_gapratio(st.Gamma[s, r], G_LAV[i, j], 1e-8, 1e-8)))
        }
    }
    printf("{txt}  5c. case E Gamma, lav rho: ratio %g (1e-12 abs and rel)\n", r_lavC)
    printf("{txt}  5d. case E Gamma, our rho: ratio %g (1e-8 abs and rel; printed, not asserted)\n",
           r_lavD)
    if (r_lavC > 1) {
        FAIL = FAIL \ ("lavaan 5c: at lavaan's own case E rho the two constructions " +
                       "disagree, ratio " + strofreal(r_lavC) +
                       " -- an arithmetic difference, not a convergence one")
    }

    printf("{txt}test_gamma case E: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }
    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_g_full(YD, (4 \ 4 \ 3 \ 3), YA, ///
                  G_LAV, LAB_LAV, R_LAV)))
mata: st_local("nfailE", strofreal(dwlsmv_g_caseE(YE, J(6, 1, 3), ///
                  G_LAVE, LAB_LAVE, R_LAVE)))
local nfail = `nfail' + `nfailE'

* ---------------------------------------------------------------------------
* ERROR PATHS. Each of these returns a plausible-looking Gamma if it is not
* caught, which is the failure mode this project is organised against.
* ---------------------------------------------------------------------------
local nerr 0

* Raw counts where a weight matrix is required: every information comes back
* scaled by N and nothing else reports it.
capture mata: dwlsmv_pairinfo(dwlsmv_pairfit(), (0 \ 1), J(1, 1, 0.2), 0.3, ///
                              J(3, 2, 100))
if _rc != 3498 {
    display as error "test_gamma: pairinfo accepted a weight matrix that does not sum to 1 (rc `=_rc')"
    local nerr = `nerr' + 1
}

* A weight matrix of the wrong shape.
capture mata: dwlsmv_pairinfo(dwlsmv_pairfit(), (0 \ 1), J(1, 1, 0.2), 0.3, ///
                              J(2, 2, 0.25))
if _rc != 3498 {
    display as error "test_gamma: pairinfo accepted a non-conformable weight matrix (rc `=_rc')"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_gamma: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_gamma: Module 4 verified -- structure, binary oracle, exact moments, centring, lavaan at matched and at own rho (case D; case E matched), normaliser guard"
exit 0
