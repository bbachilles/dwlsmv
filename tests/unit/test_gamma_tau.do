*! Module 4: the threshold block of Gamma.
*!
*! Executed under StataMP 18 on 2026-09-20 and it passes.
*!
*! THIS FILE REPLACES A TEMPLATE THAT WOULD HAVE FAILED FOR A REASON THAT HAD
*! NOTHING TO DO WITH MODULE 4. The template compared the empirical Gamma-hat
*! against fixtures/expected/case_D_gamma_tau.csv at 1e-12 and 1e-10. That
*! fixture is TIER 1 and holds the POPULATION threshold block -- built by
*! 10_fixtures.R from the declared category probabilities and the population
*! correlation lambda^2, not from the 5000 rows. Comparing a sample quantity
*! to it as though they were the same object is exactly the trap the project's notes
*! record for the two threshold files, reproduced one level up; it would have
*! failed by about 1e-2, which is sampling error and not a bug. Tier 1 is a
*! BIAS check here and never an equality target. The generator's own
*! self-check makes the same distinction -- it compares against
*! gamma_tau(qnorm(qh), qh, ...) at the SAMPLE proportions, not at the
*! population ones.
*!
*! Four checks. The first two need no external reference at all and are worth
*! more than the two that do, because they cannot agree for the wrong reason.
*!
*!   1. THE WITHIN-VARIABLE IDENTITY (spec 4.5.1), at machine precision.
*!      Z_{j,c} Z_{j,d} = Z_{j,min(c,d)} collapses, so the empirical outer
*!      product and the closed form are the same arithmetic for any sample at
*!      any N -- not asymptotically. A failure is a coding error, full stop.
*!
*!   2. THE WHOLE THRESHOLD BLOCK, cross-variable entries included, against the
*!      closed form built from the OBSERVED joint proportions. This is also an
*!      exact identity and it is the check this file exists for: both psi_tau
*!      are category lookups, so (1/N) sum_i psi_a psi_b is exactly
*!      (phat_ab - qhat_a qhat_b) / (phi_a phi_b) with phat the observed joint
*!      proportion. It covers T^2 of Gamma at machine precision and it is the
*!      only check here that exercises dwlsmv_psi_matrix()'s gather across
*!      DIFFERENT variables -- check 1 cannot, since it never leaves a block.
*!      Note what it does NOT touch: no Phi2, no rho. The stochastic part of
*!      the threshold block is phat estimating Phi2, and that is check 3.
*!
*!   3. BIAS against tier 1, as a z statistic and nothing tighter. Judged
*!      against each entry's own sampling error, because a covariance matrix's
*!      small off-diagonal entries carry large RELATIVE Monte Carlo error and
*!      the worst percentage in the matrix is noise on an entry nobody cares
*!      about (spec 7.2).
*!
*!   4. AGAINST LAVAAN, the external oracle, at 1e-12. The threshold block is
*!      immune to the normaliser question settled on 2026-09-20 -- a univariate
*!      model is saturated, so expected and empirical information coincide
*!      exactly -- which is why this block alone can be asserted this tightly
*!      while the correlation block cannot. See tests/unit/test_gamma.do.
*!
*! Reference CSVs are read as text and converted with strtoreal(). `import
*! delimited` types a numeric column as FLOAT without warning, throwing away
*! the low 29 mantissa bits -- a relative error near 4e-08 against the 1e-12
*! asserted here. A project rule.
*!
*! Accumulation is quadcross(), never cross(): cross() is core-count dependent
*! by up to ~126 ulps, so check 1 would pass on this 6-core machine and fail on
*! the 4-core milestone machine for reasons unrelated to Stata 16.
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
        display as error "test_gamma_tau: `f' failed to load"
        exit 9
    }
}

* Module 4 exists, so this test never legitimately skips: it passes or fails.

local case D

* ---- tier 1: the POPULATION threshold block -- a bias check only -----------
import delimited using "fixtures/expected/case_`case'_gamma_tau.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local gcols `r(varlist)'
mata: G_POP   = strtoreal(st_sdata(., "`gcols'"))
mata: LAB_POP = st_sdata(., "label")

* ---- tier 1: population thresholds, for the standard errors ----------------
import delimited using "fixtures/expected/case_`case'_thresholds.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TAU_POP = strtoreal(st_sdata(., "tau"))
mata: Q_POP   = strtoreal(st_sdata(., "q"))

* ---- tier 2: lavaan's Gamma, the external oracle ---------------------------
import delimited using "fixtures/expected/lavaan_case_`case'_gamma.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
qui ds label, not
local lcols `r(varlist)'
mata: G_LAV   = strtoreal(st_sdata(., "`lcols'"))
mata: LAB_LAV = st_sdata(., "label")

* ---- the data --------------------------------------------------------------
import delimited using "fixtures/data/case_`case'.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

real scalar dwlsmv_g_find(string colvector pool, string scalar want)
{
    real scalar i

    for (i = 1; i <= rows(pool); i++) {
        if (pool[i] == want) {
            return(i)
        }
    }
    return(0)
}

real scalar dwlsmv_g_tau(real matrix Y, real colvector K,
                         real matrix G_POP, string colvector LAB_POP,
                         real colvector TAU_POP, real colvector Q_POP,
                         real matrix G_LAV, string colvector LAB_LAV)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_statmap scalar m
    string colvector FAIL
    real matrix Zc, G_cf
    real colvector q, n, phi, qh
    real scalar p, N, T, j, c, d, s, t, a, b, hit, hjt
    real scalar r_within, nw, r_full, nf, r_lav, maxz, z, se, sa, sb, sab

    FAIL = J(0, 1, "")
    p    = rows(K)
    N    = rows(Y)

    st = dwls_sample_stats(Y, K)
    m  = st.map
    T  = m.T

    // ---- the observed cumulative indicators, Z[i,s] = 1{Y_i < category s} --
    // Built here rather than taken from Module 4, because checks 1 and 2 need
    // an INDEPENDENT construction of the same quantity. Reusing Module 4's
    // gather would make both checks tautological.
    Zc = J(N, T, 0)
    qh = J(T, 1, .)
    phi = J(T, 1, .)
    for (j = 1; j <= p; j++) {
        n = dwlsmv_catcounts(Y[., j], K[j])
        q = dwlsmv_cumprop(n)
        for (c = 1; c <= m.C[j]; c++) {
            s        = dwlsmv_pos_tau(m, j, c)
            Zc[., s] = (Y[., j] :< c)
            qh[s]    = q[c]
            phi[s]   = normalden(invnormal(q[c]))
        }
    }

    // ---- the closed form at the OBSERVED proportions -----------------------
    // Within a variable the joint is min(q_c, q_d) exactly; across variables it
    // is the observed joint proportion. Both are exact identities against the
    // empirical outer product -- neither involves Phi2.
    G_cf = J(T, T, .)
    for (s = 1; s <= T; s++) {
        for (t = 1; t <= T; t++) {
            G_cf[s, t] = (quadcross(Zc[., s], Zc[., t]) / N - qh[s] * qh[t]) /
                         (phi[s] * phi[t])
        }
    }

    // ---- 1. the within-variable identity (spec 7.2: 1e-12 abs, no rel) -----
    r_within = 0
    nw       = 0
    for (j = 1; j <= p; j++) {
        for (c = 1; c <= m.C[j]; c++) {
            for (d = 1; d <= m.C[j]; d++) {
                s = dwlsmv_pos_tau(m, j, c)
                t = dwlsmv_pos_tau(m, j, d)
                r_within = max((r_within,
                                dwlsmv_gapratio(st.Gamma[s, t],
                                                (min((qh[s], qh[t])) - qh[s] * qh[t]) /
                                                (phi[s] * phi[t]), 1e-12, 0)))
                nw = nw + 1
            }
        }
    }
    printf("{txt}  1. within-variable identity : ratio %g over %f entries (1e-12 abs)\n",
           r_within, nw)
    if (r_within > 1) {
        FAIL = FAIL \ ("within-variable identity: ratio " + strofreal(r_within) +
                       " -- a coding error, never sampling noise")
    }

    // ---- 2. the whole threshold block, cross-variable entries included -----
    r_full = 0
    nf     = 0
    for (s = 1; s <= T; s++) {
        for (t = 1; t <= T; t++) {
            r_full = max((r_full, dwlsmv_gapratio(st.Gamma[s, t], G_cf[s, t], 1e-12, 0)))
            nf     = nf + 1
        }
    }
    printf("{txt}  2. full threshold block     : ratio %g over %f entries (1e-12 abs)\n",
           r_full, nf)
    if (r_full > 1) {
        FAIL = FAIL \ ("threshold block: ratio " + strofreal(r_full) +
                       " -- the gather or the crossproduct is wrong, not the mathematics")
    }

    // ---- 3. BIAS against tier 1, as a z statistic --------------------------
    // NOT an equality check: G_POP is the POPULATION block, built from the
    // declared category probabilities and lambda^2, and Gamma-hat differs from
    // it by sampling error by construction. The standard error is the usual
    // one for a covariance entry, on the Z scale, divided through by the two
    // densities -- approximate, because it ignores the variability of the
    // phi(tau-hat) denominators themselves. It is a screen for a systematic
    // error, not a precise test.
    maxz = 0
    for (s = 1; s <= T; s++) {
        hit = dwlsmv_g_find(LAB_POP, m.label[s])
        if (hit == 0) {
            FAIL = FAIL \ ("bias: tier 1 has no entry labelled " + m.label[s])
            continue
        }
        for (t = 1; t <= T; t++) {
            hjt = dwlsmv_g_find(LAB_POP, m.label[t])
            if (hjt == 0) {
                continue
            }
            a   = Q_POP[hit]
            b   = Q_POP[hjt]
            sa  = a * (1 - a)
            sb  = b * (1 - b)
            sab = G_POP[hit, hjt] * normalden(TAU_POP[hit]) * normalden(TAU_POP[hjt])
            se  = sqrt((sa * sb + sab * sab) / N) /
                  (normalden(TAU_POP[hit]) * normalden(TAU_POP[hjt]))
            z    = (st.Gamma[s, t] - G_POP[hit, hjt]) / se
            maxz = max((maxz, abs(z)))
        }
    }
    printf("{txt}  3. tier-1 bias              : worst |z| %g against each entry's own SE\n",
           maxz)
    if (maxz > 5) {
        FAIL = FAIL \ ("bias: worst |z| is " + strofreal(maxz) +
                       " -- further from the population block than sampling error explains")
    }

    // ---- 4. against lavaan, at 1e-12 ---------------------------------------
    // Aligned by label, never by position. The threshold block is saturated in
    // each variable, so this is two implementations doing identical arithmetic
    // and it is asserted as such.
    r_lav = 0
    for (s = 1; s <= T; s++) {
        hit = dwlsmv_g_find(LAB_LAV, m.label[s])
        if (hit == 0) {
            FAIL = FAIL \ ("lavaan: no entry labelled " + m.label[s])
            continue
        }
        for (t = 1; t <= T; t++) {
            hjt = dwlsmv_g_find(LAB_LAV, m.label[t])
            if (hjt == 0) {
                continue
            }
            r_lav = max((r_lav,
                         dwlsmv_gapratio(st.Gamma[s, t], G_LAV[hit, hjt], 1e-12, 1e-12)))
        }
    }
    printf("{txt}  4. lavaan threshold block   : ratio %g (1e-12 abs and rel)\n", r_lav)
    if (r_lav > 1) {
        FAIL = FAIL \ ("lavaan: threshold block ratio " + strofreal(r_lav))
    }

    printf("{txt}test_gamma_tau: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_g_tau(YD, (4 \ 4 \ 3 \ 3), ///
                  G_POP, LAB_POP, TAU_POP, Q_POP, G_LAV, LAB_LAV)))

if `nfail' {
    display as error _n "test_gamma_tau: `nfail' failure(s)"
    exit 9
}

display as result "test_gamma_tau: threshold block verified -- within-variable identity, full block against the observed joint proportions, tier-1 bias, lavaan oracle"
exit 0
