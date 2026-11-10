*! dwlsmv -- Module 11: Stage 1 per group (invariance() flow)
*!
*!     void dwlsmv_m11_stage1_g(real scalar keepbound)
*!
*! Build step 3 of the group surface (2026-10-12). dwls_sample_stats on each
*! group's recoded Y (dwlsmv_m11_gsupport()), with that group's row of K, in
*! O13 order; the vector is the external dwlsmv_m11_st (O14: one struct per
*! group). Called under capture noisily (D8a), as dwlsmv_m11_stage1().
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_invariance_flow.do.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_stage1_g(real scalar keepbound)
{
    struct dwlsmv_stage1 vector sts
    pointer(real matrix) colvector Yg
    real matrix K
    real scalar G, g
    pointer() scalar p

    Yg  = *findexternal("dwlsmv_m11_Yg")
    K   = *findexternal("dwlsmv_m11_K")
    G   = rows(K)
    sts = dwlsmv_stage1(G, 1)
    for (g = 1; g <= G; g++) {
        sts[g] = dwls_sample_stats(*Yg[g], K[g, .]', keepbound)
    }
    p = crexternal("dwlsmv_m11_st")
    if (p == NULL) {
        p = findexternal("dwlsmv_m11_st")
    }
    *p = sts
}

end
