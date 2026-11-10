*! dwlsmv -- row offsets into the stacked psi_tau block (Module 4)
*!
*! Executed under StataMP 18 on 2026-09-20. Exercised by tests/unit/test_loading.do,
*! which checks the offsets by inspection on four binary variables, and by every
*! Module 4 routine.
*!
*! Koff[j] = sum_{m<j} K[m], the number of stacked rows that precede variable
*! j's lookup table in the psi_tau container described in
*! src/mata/dwlsmv_psi_tau_stack.mata.
*!
*! This is deliberately a function and NOT a field on struct dwlsmv_statmap.
*! It is derivable from map.K in three lines, and a derived quantity stored in
*! a struct is a second copy that can disagree with the first -- the project's
*! rule for the locked API. It is a function rather than three inline running
*! sums because the stacked layout is read in three places (the stack builder,
*! the psi gather and the centring diagnostic) and an offset convention that
*! lives in three copies is an offset convention that will eventually differ in
*! one of them.
*!
*! Note the asymmetry with map.Coff, which IS a field: Coff indexes s itself,
*! so every module needs it and the map is built from C. Koff indexes a
*! container that only Module 4 uses.
*!
*! src/mata/dwlsmv_types.mata must be loaded first -- the argument is a struct
*! declared there.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_koff(struct dwlsmv_statmap scalar map)
{
    real colvector off
    real scalar j, run

    off = J(map.p, 1, 0)
    run = 0
    for (j = 1; j <= map.p; j++) {
        off[j] = run
        run    = run + map.K[j]
    }

    return(off)
}

end
