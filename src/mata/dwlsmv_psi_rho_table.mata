*! dwlsmv -- one pair's psi_rho lookup table (Module 4, spec 4.5.2 / Gate 3)
*!
*! Executed under StataMP 18 on 2026-09-20. tests/unit/test_gamma.do runs it
*! through spec 4.5.3's binary oracle at 1.1e-14 over 60 (q1,q2,rho) combinations,
*! and tests/unit/test_gamma_expansion.do shows the psi it builds satisfies the
*! pointwise expansion while the naive form does not.
*!
*! THIS IS WHERE THE TWO-STAGE PROPAGATION TERM LIVES, and it is the only
*! place in the project where it appears. That is the whole point of the
*! Module 2 / Module 4 split in spec 5.1: Module 2 ESTIMATES and returns
*! ingredients, Module 4 ASSEMBLES. A propagation term that appeared once per
*! pair inside Module 2 would be a term that has to be audited once per pair.
*!
*!     psi_rho[r,s] = ( u[r,s] - sum_c I_ra[c] psi_a[r,c]
*!                             - sum_d I_rb[d] psi_b[s,d] ) / I_rr
*!
*! The two sums are the propagation: they carry the uncertainty in Module 1's
*! ESTIMATED thresholds through into rho-hat, which is a function of them.
*! Dropping them leaves `u / I_rr`, which is a different (wrong) estimator --
*! spec 4.5.2 -- and understates N*Var(rho-hat) by about 8% on a realistic 4x3
*! design.
*!
*! THE CENTRING DIAGNOSTIC CANNOT DETECT THAT OMISSION. The naive table is
*! also exactly centred, which tests/unit/test_polychoric.do measures rather
*! than asserts on faith (both E[psi] come back 0.0). So do not read a small
*! psibar as evidence this term is present. The test that separates them is
*! the derivation file's Level 4 pointwise expansion, in
*! tests/unit/test_gamma_expansion.do, which is Module 4's to own for exactly
*! the reason this routine exists: Level 4 needs an assembled psi, and this is
*! the assembly.
*!
*! psi_a and psi_b are the per-variable staircases of spec 4.5.1 --
*! dwlsmv_psi_tau_table() output, or equivalently a variable's row/column block
*! of the dwlsmv_psi_tau_stack() container. They are indexed by CATEGORY, which
*! is why the stacked tables rather than a gathered N x T matrix are what
*! crosses the module boundary.
*!
*! The result is a K1 x K2 table and nothing larger: psi_rho depends on
*! observation i only through the cell (r_i, s_i) it lands in, because the
*! threshold influence functions are themselves category lookups. Built once
*! per pair, then indexed. The O(N) work in Gamma is the gather and the
*! crossproduct, not this.
*!
*! Requires dwlsmv_types.mata.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real matrix dwlsmv_psi_rho_table(struct dwlsmv_pairfit scalar f,
                                 real matrix psi_a, real matrix psi_b)
{
    real scalar K1, K2

    K1 = rows(f.u)
    K2 = cols(f.u)

    if (rows(psi_a) != K1 | cols(psi_a) != rows(f.I_ra)) {
        _error("dwlsmv_psi_rho_table: psi_a does not match u and I_ra")
    }
    if (rows(psi_b) != K2 | cols(psi_b) != rows(f.I_rb)) {
        _error("dwlsmv_psi_rho_table: psi_b does not match u and I_rb")
    }
    if (f.I_rr <= 0) {
        errprintf("dwlsmv_psi_rho_table: I_rho_rho is %s, which cannot be a" +
                  " Fisher information.\n", strofreal(f.I_rr))
        errprintf("  This pair's information quantities are not usable and" +
                  " Gamma cannot be formed\n")
        errprintf("  from them. Inspect the pair's contingency table -- spec" +
                  " 4.5.6.\n")
        _error("dwlsmv_psi_rho_table: I_rho_rho is not positive")
    }

    // psi_a * I_ra is K1 x 1 and psi_b * I_rb is K2 x 1; the row vector
    // transpose broadcasts across columns. Written as matrix products rather
    // than as the two literal sums so the index order cannot be written
    // backwards without a conformability error.
    return((f.u :- (psi_a * f.I_ra) :- (psi_b * f.I_rb)') :/ f.I_rr)
}

end
