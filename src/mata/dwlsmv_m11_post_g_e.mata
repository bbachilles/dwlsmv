*! dwlsmv -- Module 11: the canonical texts into e() (group surface, step 4)
*!
*!     void dwlsmv_m11_post_g_e()
*!
*! Build step 4 (2026-10-12). Called by dwlsmv.ado after ereturn post. Sets,
*! from Mata because the texts hold char(10) (spec 4.8, what each fit posts
*! for the exemption, decided 2026-10-12):
*!   e(canon)             the fit's own canonical text, dwlsmv_levelgen's
*!                        output (dwlsmv_m11_gtext), byte for byte
*!   e(canon_configural)  the generator's configural output for the same
*!                        configural text and K (dwlsmv_m11_gconf); equal to
*!                        e(canon) for a configural fit
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_invariance_post.do.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_post_g_e()
{
    st_global("e(canon)",            *findexternal("dwlsmv_m11_gtext"))
    st_global("e(canon_configural)", *findexternal("dwlsmv_m11_gconf"))
}

end
