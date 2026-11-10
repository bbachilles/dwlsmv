*! dwlsmv -- Module 11: removes the driver's Mata externals (three, plus D11's dwlsmv_m11_rec)
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*! Called by dwlsmv.ado on every exit; see dwlsmv_m11_fit.mata for the steps.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_clear()
{
    rmexternal("dwlsmv_m11_st")
    rmexternal("dwlsmv_m11_Y")
    rmexternal("dwlsmv_m11_txt")
    rmexternal("dwlsmv_m11_rec")
    // the invariance() flow (build step 3, 2026-10-12)
    rmexternal("dwlsmv_m11_Yg")
    rmexternal("dwlsmv_m11_K")
    rmexternal("dwlsmv_m11_sup")
    rmexternal("dwlsmv_m11_gtext")
    rmexternal("dwlsmv_m11_gpar")
    rmexternal("dwlsmv_m11_gfit")
    rmexternal("dwlsmv_m11_gd0")
    rmexternal("dwlsmv_m11_gconf")
}

end
