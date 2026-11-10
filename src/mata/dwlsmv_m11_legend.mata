*! dwlsmv -- Module 11: D8b's legend after a Stage 1 abort
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*! One line per indicator, yj = name, in e(oyvars) order.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_legend(string scalar indlist)
{
    string rowvector ind
    real scalar j

    ind = tokens(indlist)
    errprintf("Positional labels above, in the model's indicator order:\n")
    for (j = 1; j <= cols(ind); j++) {
        errprintf("    y%f = %s\n", j, ind[j])
    }
}

// The same legend after a first stage that succeeded but printed positional
// labels: a pair kept at the bound or given the one-empty-cell repair, in any
// group. Not an error, so printf.
void dwlsmv_m11_legend_note(string scalar indlist, real scalar grouped)
{
    struct dwlsmv_stage1 vector sts
    pointer(real matrix) colvector Yg
    real matrix Y
    string rowvector ind
    real scalar g, j, show

    sts = *findexternal("dwlsmv_m11_st")
    if (grouped) {
        Yg = *findexternal("dwlsmv_m11_Yg")
    }
    else {
        Y  = *findexternal("dwlsmv_m11_Y")
        Yg = &Y
    }
    show = 0
    for (g = 1; g <= length(sts); g++) {
        if (any(sts[g].bound)) {
            show = 1
        }
        else if (any(dwlsmv_emptycell_pairs(*Yg[g], sts[g].map))) {
            show = 1
        }
    }
    if (show) {
        ind = tokens(indlist)
        printf("Positional labels above, in the model's indicator order:\n")
        for (j = 1; j <= cols(ind); j++) {
            printf("    y%f = %s\n", j, ind[j])
        }
    }
}

end
