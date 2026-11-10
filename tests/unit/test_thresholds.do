*! Module 1: Stage 1 thresholds and their influence-function staircase.
*!
*! Executed under StataMP 18 on Windows, 2026-09-19, and it passes. Its first
*! run found a real defect -- not in Module 1, but in how this suite reads
*! committed CSV references. See the import block below.
*!
*! Six checks, and they are not redundant -- each fails on something the
*! others do not.
*!
*!   1. ORACLE (external). tau-hat from case_D.csv against lavaan's `th`,
*!      aligned BY LABEL through the Section 4.2 index map, never by position.
*!      This is the only check here with an external source, and it is the one
*!      that pins "invnormal of the cumulative sample proportion, and nothing
*!      more" to a second implementation.
*!
*!   2. BIAS (tier 1, NOT an equality target). The tier-1 fixture holds
*!      POPULATION thresholds and the tier-2 fixture holds SAMPLE estimates.
*!      They sit side by side, look alike, and differ by sampling error --
*!      about 1.8 SE on y1's first threshold. Diffing them as if they were the
*!      same quantity is a convincing false alarm, so this compares them as a z
*!      statistic against each threshold's own standard error and nothing
*!      tighter. A project rule.
*!
*!   3. THE WITHIN-VARIABLE IDENTITY (spec 4.5.1), at machine precision, with
*!      no external dependency at all. N^-1 sum_i psi_a psi_b equals
*!      (min(q_a,q_b) - q_a q_b) / (phi_a phi_b) EXACTLY, for any sample at any
*!      N -- not asymptotically. A failure is a coding error, never sampling
*!      noise. This is the check that is available on day one and is worth more
*!      than the lavaan comparison, because it cannot agree for the wrong
*!      reason.
*!
*!   4. THE TRIDIAGONAL ROUTE. The derivation file's second, independent route
*!      to the same psi: the univariate ordered-probit score solved through its
*!      own information matrix, A^-1 s. It agrees with the staircase at every
*!      category value, and A^-1 agrees with the within-variable block. This
*!      catches a staircase that is self-consistently wrong -- the outer
*!      product in check 3 is comparing psi against a formula derived from the
*!      same q, so a sign or indexing error that moved both would survive it.
*!
*!   5. CENTRING (spec 4.5.5), as a diagnostic and nothing more. It is checked
*!      because a non-zero mean is a real derivative error, and it is not
*!      trusted because the NAIVE influence function is also exactly centred.
*!      It catches derivative errors, not omissions.
*!
*!   6. ERROR PATHS. Every degenerate-category and bad-coding case aborts, and
*!      dwlsmv_catdegen() tells the boundary mode from the interior one, which
*!      spec 4.5.6 requires it to do.
*!
*! This test deliberately makes NO claim about the cross-variable entries of
*! Gamma. Those need Phi2 and the polychoric correlations, which are Modules 2
*! and 4. The N x T psi matrix built here is the gather-by-lookup step Module 4
*! will do for the full N x D matrix, exercised early on the third of Gamma
*! that Stage 1 can already produce.
*!
*! Accumulation is quadcross(), not cross(). cross() is core-count dependent by
*! up to ~126 ulps, so check 3's machine-precision assertion would pass on this
*! 6-core machine and fail on the 4-core milestone machine for reasons having
*! nothing to do with Stata 16. A project rule.
*!
*! The body is a function rather than interactive statements because it touches
*! a struct, and at Mata's interactive level every variable is transmorphic.
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
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_catdegen.mata          ///
             src/mata/dwlsmv_cumprop.mata           ///
             src/mata/dwlsmv_thresholds.mata        ///
             src/mata/dwlsmv_psi_tau_table.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_thresholds: `f' failed to load"
        exit 9
    }
}

* Module 1 exists, so this test never legitimately skips: it passes or fails.

* ---------------------------------------------------------------------------
* READ EVERY REFERENCE CSV AS TEXT AND CONVERT WITH strtoreal().
*
* `import delimited` types a numeric column as FLOAT. It does not warn, the
* column still displays to the digits you expect, and the low 29 bits of every
* mantissa are gone. Measured here on the first run of this test: lavaan's
* -1.05549280152646 arrives as -1.132f800000000X+000 where the file says
* -1.132f8038008b6X+000, a relative error of 4e-08.
*
* That is not a rounding nuisance, it is a ceiling on every comparison in this
* suite that reads a committed CSV. Section 7.2 asks for 1e-10 on Gamma and
* 1e-12 on the within-variable identity; neither is reachable through a float,
* and the natural response to the failure -- loosen the tolerance until it
* passes -- converts a precision bug into a permanent silent one.
*
* tests/unit/test_dta_roundtrip.do already used this idiom for its probe. It
* is the required one for any fixture whose digits matter.
* ---------------------------------------------------------------------------

* ---- tier 2: lavaan's sample thresholds, the external oracle ----------------
import delimited using "fixtures/expected/lavaan_case_D_th.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: LAB_LAV = st_sdata(., "label")
mata: TAU_LAV = strtoreal(st_sdata(., "value"))

* ---- tier 1: POPULATION thresholds, a bias check and never an equality one --
import delimited using "fixtures/expected/case_D_thresholds.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: LAB_POP = st_sdata(., "label")
mata: TAU_POP = strtoreal(st_sdata(., "tau"))
mata: Q_POP   = strtoreal(st_sdata(., "q"))

* ---- the data -------------------------------------------------------------
import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

// Look a label up in a string colvector; 0 if absent. Alignment is by name,
// so a future lavaan that reorders its output fails loudly rather than
// comparing the wrong pair of numbers.
real scalar dwlsmv_t_find(string colvector pool, string scalar want)
{
    real scalar i

    for (i = 1; i <= rows(pool); i++) {
        if (pool[i] == want) {
            return(i)
        }
    }
    return(0)
}

real scalar dwlsmv_t_thresholds(real matrix Y,
                                real colvector K,
                                string colvector LAB_LAV, real colvector TAU_LAV,
                                string colvector LAB_POP, real colvector TAU_POP,
                                real colvector Q_POP)
{
    struct dwlsmv_statmap scalar m
    string colvector FAIL
    real matrix psi, tbl, A, G_emp, G_cf, Ainv
    real colvector n, q, tau, pmf, phi, sc, idx, cm
    real scalar p, N, T, j, c, d, s, hit, y, ratio, se, z, maxz
    real scalar r_within, r_tri, r_ainv, maxcen, nv

    FAIL = J(0, 1, "")
    p    = rows(K)
    N    = rows(Y)
    m    = dwlsmv_statmap_build(K)
    T    = m.T

    if (cols(Y) != p) {
        FAIL = FAIL \ ("data has " + strofreal(cols(Y)) +
                       " columns, K names " + strofreal(p))
        return(rows(FAIL))
    }

    // ---- gather tau and the N x T influence matrix by table lookup --------
    tau = J(T, 1, .)
    psi = J(N, T, 0)

    for (j = 1; j <= p; j++) {
        n   = dwlsmv_catcounts(Y[., j], K[j])
        q   = dwlsmv_cumprop(n)
        tbl = dwlsmv_psi_tau_table(q)
        idx = Y[., j] :+ 1

        if (rows(tbl) != K[j] | cols(tbl) != m.C[j]) {
            FAIL = FAIL \ ("y" + strofreal(j) + ": psi table is " +
                           strofreal(rows(tbl)) + " x " + strofreal(cols(tbl)) +
                           ", expected " + strofreal(K[j]) + " x " +
                           strofreal(m.C[j]))
            continue
        }

        for (c = 1; c <= m.C[j]; c++) {
            s        = dwlsmv_pos_tau(m, j, c)
            tau[s]   = invnormal(q[c])
            psi[., s] = tbl[idx, c]
        }

        // The estimator's own output must agree with the pieces it is built
        // from -- a cheap guard against dwlsmv_thresholds() drifting from
        // invnormal(dwlsmv_cumprop(.)).
        if (max(abs(dwlsmv_thresholds(Y[., j], K[j]) :- invnormal(q))) != 0) {
            FAIL = FAIL \ ("y" + strofreal(j) +
                           ": dwlsmv_thresholds disagrees with invnormal(cumprop)")
        }

        // Strict ordering, which the estimator gets for free and must show.
        for (c = 2; c <= m.C[j]; c++) {
            if (invnormal(q[c]) <= invnormal(q[c - 1])) {
                FAIL = FAIL \ ("y" + strofreal(j) + ": thresholds not strictly increasing at c = " +
                               strofreal(c))
            }
        }
    }

    // ---- 1. ORACLE: against lavaan, by label ------------------------------
    for (s = 1; s <= T; s++) {
        hit = dwlsmv_t_find(LAB_LAV, m.label[s])
        if (hit == 0) {
            FAIL = FAIL \ ("oracle: lavaan has no threshold labelled " + m.label[s])
            continue
        }
        // Spec 7.2, thresholds row: eps_abs = eps_rel = 1e-8.
        ratio = dwlsmv_gapratio(tau[s], TAU_LAV[hit], 1e-8, 1e-8)
        if (ratio > 1) {
            FAIL = FAIL \ ("oracle: " + m.label[s] + " is " + strofreal(tau[s]) +
                           ", lavaan says " + strofreal(TAU_LAV[hit]) +
                           " (tolerance ratio " + strofreal(ratio) + ")")
        }
    }
    ratio = 0
    for (s = 1; s <= T; s++) {
        hit = dwlsmv_t_find(LAB_LAV, m.label[s])
        if (hit > 0) {
            ratio = max((ratio, abs(tau[s] - TAU_LAV[hit])))
        }
    }
    printf("{txt}  1. lavaan oracle      : worst absolute gap %g over %f thresholds\n",
           ratio, T)

    // ---- 2. BIAS against tier 1, as a z statistic -------------------------
    // NOT an equality check. The tier-1 fixture is the POPULATION threshold;
    // the estimate differs from it by sampling error by construction.
    maxz = 0
    for (s = 1; s <= T; s++) {
        hit = dwlsmv_t_find(LAB_POP, m.label[s])
        if (hit == 0) {
            FAIL = FAIL \ ("bias: tier 1 has no threshold labelled " + m.label[s])
            continue
        }
        se   = sqrt(Q_POP[hit] * (1 - Q_POP[hit]) / N) / normalden(TAU_POP[hit])
        z    = (tau[s] - TAU_POP[hit]) / se
        maxz = max((maxz, abs(z)))
    }
    printf("{txt}  2. tier-1 bias        : worst |z| %g against each threshold's own SE\n",
           maxz)
    if (maxz > 4) {
        FAIL = FAIL \ ("bias: worst |z| is " + strofreal(maxz) +
                       " -- the estimator is off the population value by more than sampling error explains")
    }

    // ---- 3. THE WITHIN-VARIABLE IDENTITY, at machine precision ------------
    // quadcross(), never cross(): see the header.
    G_emp = quadcross(psi, psi) :/ N

    G_cf = J(T, T, .)
    for (j = 1; j <= p; j++) {
        n = dwlsmv_catcounts(Y[., j], K[j])
        q = dwlsmv_cumprop(n)
        for (c = 1; c <= m.C[j]; c++) {
            for (d = 1; d <= m.C[j]; d++) {
                G_cf[dwlsmv_pos_tau(m, j, c), dwlsmv_pos_tau(m, j, d)] =
                    (min((q[c], q[d])) - q[c] * q[d]) /
                    (normalden(invnormal(q[c])) * normalden(invnormal(q[d])))
            }
        }
    }

    // Compare only where the closed form is defined: the within-variable
    // blocks. The cross-variable entries need Phi2 and rho, which are Modules
    // 2 and 4, and this test claims nothing about them.
    r_within = 0
    nv       = 0
    for (j = 1; j <= p; j++) {
        for (c = 1; c <= m.C[j]; c++) {
            for (d = 1; d <= m.C[j]; d++) {
                s        = dwlsmv_pos_tau(m, j, c)
                hit      = dwlsmv_pos_tau(m, j, d)
                // Spec 7.2: eps_abs 1e-12, NO relative term.
                r_within = max((r_within,
                                dwlsmv_gapratio(G_emp[s, hit], G_cf[s, hit], 1e-12, 0)))
                nv       = nv + 1
            }
        }
    }
    printf("{txt}  3. within-variable id : tolerance ratio %g over %f entries (1e-12 abs)\n",
           r_within, nv)
    if (r_within > 1) {
        FAIL = FAIL \ ("within-variable identity: tolerance ratio " +
                       strofreal(r_within) + " -- this is a coding error, never noise")
    }

    // ---- 4. THE TRIDIAGONAL ROUTE -----------------------------------------
    // Independent of check 3: the score-and-information route of the
    // derivation file, which never mentions a cumulative proportion.
    r_tri  = 0
    r_ainv = 0
    for (j = 1; j <= p; j++) {
        n   = dwlsmv_catcounts(Y[., j], K[j])
        q   = dwlsmv_cumprop(n)
        tbl = dwlsmv_psi_tau_table(q)
        pmf = n :/ N
        phi = normalden(invnormal(q))

        A = J(m.C[j], m.C[j], 0)
        for (c = 1; c <= m.C[j]; c++) {
            A[c, c] = phi[c] ^ 2 * (1 / pmf[c] + 1 / pmf[c + 1])
            if (c < m.C[j]) {
                A[c, c + 1] = -phi[c] * phi[c + 1] / pmf[c + 1]
                A[c + 1, c] = A[c, c + 1]
            }
        }

        // A^-1 s equals psi, at every category value.
        for (y = 0; y <= m.C[j]; y++) {
            sc = J(m.C[j], 1, 0)
            for (c = 1; c <= m.C[j]; c++) {
                sc[c] = phi[c] * ((y == c - 1) / pmf[c] - (y == c) / pmf[c + 1])
            }
            r_tri = max((r_tri,
                         dwlsmv_gapratio(lusolve(A, sc), tbl[y + 1, .]', 1e-12, 0)))
        }

        // A^-1 equals the within-variable block.
        Ainv = lusolve(A, I(m.C[j]))
        for (c = 1; c <= m.C[j]; c++) {
            for (d = 1; d <= m.C[j]; d++) {
                r_ainv = max((r_ainv,
                              dwlsmv_gapratio(Ainv[c, d],
                                              G_cf[dwlsmv_pos_tau(m, j, c),
                                                   dwlsmv_pos_tau(m, j, d)],
                                              1e-12, 0)))
            }
        }
    }
    printf("{txt}  4. tridiagonal route  : A^-1 s vs psi ratio %g, A^-1 vs block ratio %g\n",
           r_tri, r_ainv)
    if (r_tri > 1) {
        FAIL = FAIL \ ("tridiagonal: A^-1 s does not reproduce the staircase, ratio " +
                       strofreal(r_tri))
    }
    if (r_ainv > 1) {
        FAIL = FAIL \ ("tridiagonal: A^-1 does not equal the within-variable block, ratio " +
                       strofreal(r_ainv))
    }

    // ---- 5. CENTRING, as a diagnostic (spec 4.5.5) ------------------------
    // Reported and loosely asserted. It does NOT establish that a propagation
    // term is present -- the naive influence function is also exactly centred.
    cm     = (quadcolsum(psi) :/ N)'
    maxcen = max(abs(cm))
    printf("{txt}  5. centring diagnostic: worst |mean psi| %g (derivative check only)\n",
           maxcen)
    if (maxcen > 1e-12) {
        FAIL = FAIL \ ("centring: worst |mean psi| is " + strofreal(maxcen))
    }

    printf("{txt}test_thresholds: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_thresholds(YD, (4 \ 4 \ 3 \ 3), ///
                  LAB_LAV, TAU_LAV, LAB_POP, TAU_POP, Q_POP)))

* ---------------------------------------------------------------------------
* 6. ERROR PATHS.
*
* Spec 4.5.6 requires Module 1 to distinguish an empty INTERIOR category from
* an empty BOUNDARY one, because their consequences differ: the boundary mode
* is loud and kills a threshold outright, the interior mode is silent and makes
* Gamma rank-deficient while DWLS point estimates stay correct. Checking that
* both abort is not enough -- the classification itself is asserted.
* ---------------------------------------------------------------------------
local nerr 0

mata: st_local("d_none", strofreal(dwlsmv_catdegen((3 \ 4 \ 5))))
mata: st_local("d_bnd",  strofreal(dwlsmv_catdegen((0 \ 4 \ 5))))
mata: st_local("d_int",  strofreal(dwlsmv_catdegen((3 \ 0 \ 5))))
mata: st_local("d_both", strofreal(dwlsmv_catdegen((0 \ 2 \ 0 \ 5))))

if `d_none' != 0 {
    display as error "test_thresholds: catdegen called a complete table degenerate (`d_none')"
    local nerr = `nerr' + 1
}
if `d_bnd' != 1 {
    display as error "test_thresholds: catdegen did not report an empty boundary category as 1 (`d_bnd')"
    local nerr = `nerr' + 1
}
if `d_int' != 2 {
    display as error "test_thresholds: catdegen did not report an empty interior category as 2 (`d_int')"
    local nerr = `nerr' + 1
}
if `d_both' != 3 {
    display as error "test_thresholds: catdegen did not report both modes as 3 (`d_both')"
    local nerr = `nerr' + 1
}

* Each of these returns a plausible-looking number if it is not caught here.
* Asserted as rc 3498, the _error() code, not merely nonzero: the Rank-deficient
* Gamma decision (E-038, 2026-09-23) rests on the interior abort
* being this one, and a crash for an unrelated reason would also be nonzero.
capture mata: dwlsmv_thresholds((1 \ 1 \ 1 \ 1), 3)
if _rc != 3498 {
    display as error "test_thresholds: empty boundary category gave rc `=_rc', expected 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_thresholds((0 \ 0 \ 2 \ 2), 3)
if _rc != 3498 {
    display as error "test_thresholds: empty interior category gave rc `=_rc', expected 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_catcounts((0 \ 1 \ .), 2)
if _rc != 3498 {
    display as error "test_thresholds: catcounts on a missing value gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_catcounts((0 \ 1 \ 1.5), 2)
if _rc != 3498 {
    display as error "test_thresholds: catcounts on a non-integer code gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_catcounts((0 \ 1 \ 2), 2)
if _rc != 3498 {
    display as error "test_thresholds: catcounts on a code above K-1 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_catcounts((0 \ 1 \ -1), 2)
if _rc != 3498 {
    display as error "test_thresholds: catcounts on a negative code gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_psi_tau_table((0.3 \ 0.3))
if _rc != 3498 {
    display as error "test_thresholds: psi_tau_table on a non-increasing q gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_psi_tau_table((0 \ 0.4))
if _rc != 3498 {
    display as error "test_thresholds: psi_tau_table on q = 0 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_thresholds: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_thresholds: Module 1 verified -- lavaan oracle, tier-1 bias, within-variable identity, tridiagonal route, centring, error paths"
exit 0
