*! dwlsmv -- build the sample-statistic index map (Module 0-1, spec Section 4.2)
*!
*! Keyed on K alone. No data, no estimates, no model -- that is what lets this
*! run before Stage 1 and lets every later module take the map as an input
*! rather than rebuild it. Two reconstructions in two places is how the copies
*! drift apart.
*!
*! The map does NOT shrink when a category happens to be unobserved. C_j as
*! declared governs it; the resulting degenerate row of Gamma is Section
*! 4.5.6's problem, caught by Module 10's rank check. A map that resized would
*! make Gamma non-conformable with Delta for reasons no error message would
*! explain.
*!
*! Ordering is lavaan's and is verified against it, not adopted by convention:
*! all thresholds first, grouped by variable and ascending within variable,
*! then all pairs in column-major lower-triangle order. See the oracle in
*! tests/unit/test_statmap.do, which pins the exact label sequence that
*! lavInspect(fit, "gamma") produces for the four-indicator design.
*!
*! src/mata/dwlsmv_types.mata must be loaded before this file, because the
*! struct it returns is declared there. Order beyond that does not matter: no
*! definition file issues `mata clear`, so loading one cannot wipe another.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

struct dwlsmv_statmap scalar dwlsmv_statmap_build(real colvector K)
{
    struct dwlsmv_statmap scalar m
    real scalar i, j, c, pos, run

    if (rows(K) < 1) {
        _error("dwlsmv_statmap_build: K must name at least one variable")
    }
    if (cols(K) != 1) {
        _error("dwlsmv_statmap_build: K must be a column vector")
    }
    for (j = 1; j <= rows(K); j++) {
        if (K[j] < 2 | K[j] != floor(K[j])) {
            _error("dwlsmv_statmap_build: every K must be an integer >= 2")
        }
    }

    m.p = rows(K)
    m.K = K
    m.C = K :- 1
    m.T = sum(m.C)
    m.P = m.p * (m.p - 1) / 2
    m.D = m.T + m.P

    // Coff[j] = sum_{m<j} C[m]. Built by running total rather than by a
    // library call, so the file leans on as little as possible.
    m.Coff = J(m.p, 1, 0)
    run    = 0
    for (j = 1; j <= m.p; j++) {
        m.Coff[j] = run
        run       = run + m.C[j]
    }

    m.stat_type = J(m.D, 1, 0)
    m.var_i     = J(m.D, 1, 0)
    m.var_j     = J(m.D, 1, 0)
    m.thr_k     = J(m.D, 1, 0)
    m.label     = J(m.D, 1, "")

    // Threshold block.
    for (j = 1; j <= m.p; j++) {
        for (c = 1; c <= m.C[j]; c++) {
            pos              = m.Coff[j] + c
            m.stat_type[pos] = 1
            m.var_i[pos]     = j
            m.var_j[pos]     = 0
            m.thr_k[pos]     = c
            m.label[pos]     = "y" + strofreal(j) + "|t" + strofreal(c)
        }
    }

    // Correlation block. p == 1 leaves this empty, which is correct.
    for (i = 1; i <= m.p - 1; i++) {
        for (j = i + 1; j <= m.p; j++) {
            pos              = m.T + (i - 1) * m.p - i * (i - 1) / 2 + (j - i)
            m.stat_type[pos] = 2
            m.var_i[pos]     = i
            m.var_j[pos]     = j
            m.thr_k[pos]     = 0
            m.label[pos]     = "y" + strofreal(i) + "~~y" + strofreal(j)
        }
    }

    return(m)
}

end
