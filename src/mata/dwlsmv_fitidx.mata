*! dwlsmv -- fit indices, option A: scaled RMSEA, CFI and TLI (Module 12)
*!
*! Written 2026-10-21 from derivation @sec-fit-indices, Result FI1, with the
*! conventions F1-F4 frozen by the user 2026-10-21 (F4 as amended).
*!
*!     real rowvector dwlsmv_fitidx(real scalar T3, real scalar D,
*!                                  real scalar n, real scalar G,
*!                                  real scalar T3B, real scalar DB)
*!
*! Returns (RMSEA, CFI, TLI), family S (lavaan 0.7.2's *.scaled). T3, D the
*! model's statistic and df; n = N - G, the divisor of T = (N - G) F; G the
*! number of groups; T3B, DB the c1 baseline's (dwlsmv_fitidx_base).
*!   F3  D <= 0, or T3 or D missing: all three missing.
*!   F1  RMSEA = sqrt(G) sqrt(max(E_L, 0) / (n D)), E_L = T3 - D.
*!   F4  two INDEPENDENT conditions:
*!       model side E_L decides only the model's clamp (RMSEA 0 when
*!         E_L <= 0; max(E_L, 0) in CFI); never whether CFI/TLI exist;
*!       baseline side E_B = T3B - DB decides only whether CFI and TLI
*!         exist (missing when E_B <= 0, T3B missing or DB <= 0); never
*!         RMSEA.
*!       E_B > 0:  CFI = 1 - max(E_L, 0) / max(E_L, E_B)
*!                 TLI = 1 - E_L DB / (E_B D)   (raw E_L: unclipped, F2)
*! Exact comparisons throughout. n <= 0 or missing, or G not a positive
*! integer, aborts rc 3498.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real rowvector dwlsmv_fitidx(real scalar T3, real scalar D, real scalar n,
    real scalar G, real scalar T3B, real scalar DB)
{
    real scalar EL, EB, rmsea, cfi, tli

    if (missing(n) | n <= 0) {
        errprintf("dwlsmv_fitidx: the divisor n = %g must be positive.\n", n)
        _error("dwlsmv_fitidx: n must be positive")
    }
    if (missing(G) | G < 1 | G != trunc(G)) {
        errprintf("dwlsmv_fitidx: G = %g must be a positive integer.\n", G)
        _error("dwlsmv_fitidx: G must be a positive integer")
    }
    if (missing(T3) | missing(D)) {
        return((., ., .))
    }
    if (D <= 0) {
        return((., ., .))
    }
    // model side: the clamp only
    EL    = T3 - D
    rmsea = sqrt(G) * sqrt((EL > 0 ? EL : 0) / (n * D))
    // baseline side: whether CFI and TLI exist, from the baseline alone
    if (missing(T3B) | missing(DB)) {
        return((rmsea, ., .))
    }
    if (DB <= 0) {
        return((rmsea, ., .))
    }
    EB = T3B - DB
    if (EB <= 0) {
        return((rmsea, ., .))
    }
    cfi = 1 - (EL > 0 ? EL : 0) / (EL > EB ? EL : EB)
    tli = 1 - (EL * DB) / (EB * D)
    return((rmsea, cfi, tli))
}

end
