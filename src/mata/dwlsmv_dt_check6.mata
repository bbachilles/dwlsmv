*! dwlsmv -- difftest's ONE entry to check 6 (Module 12)
*!
*! Written 2026-10-23. Check 6's contract is open (E-023
*! plan.md: target B, reported-fit consistency including diagonal SEs and
*! T3 + df, is the leading candidate; A, e(b) and e(T3), the fallback). The
*! front end calls check 6 only through this function, so the quantities,
*! the measure and the threshold can be changed here without touching it.
*!
*!     void dwlsmv_dt_check6(real rowvector b0, real matrix V0, real scalar T0,
*!                           real matrix bound0, real scalar df0,
*!                           real rowvector b1, real matrix V1, real scalar T1,
*!                           real matrix bound1, real scalar df1,
*!                           string scalar name)
*!
*! (b0, V0, T0, bound0, df0) the stored e(b), e(V), e(T3), e(bound),
*! e(df_T3); (b1, ...) the refit's, with V1 posted as e(V) would be. TODAY:
*! dwlsmv_dtrefit on e(b), e(T3), e(bound) and e(df_T3) at dwlsmv_dt_reftol()
*! (1e-10, G = 1); V0 and V1 are carried for the contract and NOT compared.
*! Returns silently or aborts rc 3498 (dwlsmv_dtrefit's messages).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_dt_check6(real rowvector b0, real matrix V0, real scalar T0,
                      real matrix bound0, real scalar df0,
                      real rowvector b1, real matrix V1, real scalar T1,
                      real matrix bound1, real scalar df1,
                      string scalar name)
{
    if (rows(V0) != cols(b0) | rows(V1) != cols(b1)) {
        _error("dwlsmv check 6: e(V) does not match e(b)")
    }
    dwlsmv_dtrefit(b0, T0, bound0, df0, b1, T1, bound1, df1, dwlsmv_dt_reftol(), name)
}

end
