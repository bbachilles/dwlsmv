*! dwlsmv -- DIFFTEST validation check 8: H0 is nested in H1, on the parameter
*! tables (Module 12)
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9; rule 0 (factor order)
*! added by the main session (user). Executed 2026-10-08: test_dtnest.do green,
*! mutants N1-N14 and C1-C7 each turned it red.
*!
*!     real scalar dwlsmv_dtnest(struct dwlsmv_parmap scalar par0,
*!                               struct dwlsmv_parmap scalar par1,
*!                               real scalar df0, real scalar df1,
*!                               string scalar name0, string scalar name1
*!                               | cc0, can0, inv0, K0, cc1, can1, inv1, K1,
*!                               real scalar exempt)
*!
*! GENERATED CHAINS (built 2026-10-13, exemption gate step 2b; spec check 8).
*! The eight optional inputs are each model's e(canon_configural), e(canon),
*! e(invariance) and e(K); exempt (optional output) is set to 1 if the pair
*! passed through the exemption. Order: the beta guard, then
*! dwlsmv_dtgenerated, then (if not exempt) the tie refusal and rules 0-4. An
*! exempt pair skips rules 0-2, may have d = 0, and must have d = df0 - df1.
*! Called with 6 arguments it behaves exactly as before (never exempt).
*!
*! GUARD (user, 2026-10-10, spec 4.8 O2, option (e)): a pair in which either
*! table has groups and a beta row (free or fixed nonzero) is refused first, rc
*! 3498: beta DIFFTEST was not validated (follow-up NOT MET). Executed
*! 2026-10-10: tests/unit/test_dtnest_mgbeta.do green, mutants G1-G5 each red.
*!
*! Spec Module 12, Validation contract, check 8 (decided 2026-10-08). par0 is
*! H0 (the more restricted model), par1 is H1. df0 and df1 are the two fits'
*! e(df_T3) (the fit struct's D_stat). name0 and name1 only label messages.
*! Returns d = q1 - q0 on success; aborts with rc 3498 otherwise.
*!
*! G > 1 (user, 2026-10-12, option (b)): after the beta guard, a pair in which
*! either table has a TIE (two free rows sharing one fidx) is refused, rc 3498,
*! fail-closed: the rules below never read ties, so H0 tying loadings and H1
*! tying thresholds would otherwise pass though not nested. A tie rule is a
*! deferred gate. Executed 2026-10-12: tests/unit/test_dtnest_groups.do green;
*! mutants listed in its header.
*!
*! THE RULE. Rows are identified by (ptype, grp, row, col) (spec 4.8 G6; an
*! empty or missing grp reads 1). A row ABSENT from a table counts as fixed at
*! its type's absent-row default (spec 4.8 G4: omega 1, every other type 0).
*! Then:
*!   0. every row present in both carries the same name (user, 2026-10-08).
*!      Names hold the factor and indicator names, so a pair declaring its
*!      factors in a different order, index-identical otherwise, is refused
*!      here instead of being compared position by position;
*!   1. every row free in H0 is present and free in H1;
*!   2. every row fixed in both, present or absent, is fixed at the same value
*!      (exact comparison; both come from the same text through strtoreal());
*!   3. d = q1 - q0 >= 1, q = par.npar;
*!   4. d equals df0 - df1.
*! A row fixed in H0 and free in H1 is the restriction and is accepted whatever
*! the values. A free row's value is its START and is never read.
*!
*! A pair nested only through a reparameterization (a different marker, or
*! marker against unit variance) fails rule 1 and is refused; the message says
*! why. Nothing is reordered or repaired.
*!
*! WHAT IT DOES NOT DO. No data, no estimate, no theta. It does not compare
*! indicator sets, indicator order, keepbound or metric (checks 4-5, which read
*! e()), does not look at T0 - T1 (check 10) or at the boundary (check 11), and
*! does not run check 6 (the exemption's condition (4) is the front end's). It
*! does not check nesting under ties (refused, above) outside the exemption. It does not modify either table (Mata passes by
*! reference; nothing here assigns to an argument).
*!
*! Load order: dwlsmv_types, dwlsmv_canon_lines, dwlsmv_levelgen,
*! dwlsmv_dtgenerated.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Group of row i: grp[i], or 1 when grp is empty or missing there (a
// single-group table built without it).
real scalar dwlsmv_dtn_grp(struct dwlsmv_parmap scalar par, real scalar i)
{
    if (rows(par.grp) < i) {
        return(1)
    }
    if (missing(par.grp[i])) {
        return(1)
    }
    return(par.grp[i])
}

// Absent-row default (spec 4.8 G4, decided 2026-10-09): omega (ptype 6) reads 1,
// every other ptype 0.
real scalar dwlsmv_dtn_absent(real scalar t)
{
    return(t == 6 ? 1 : 0)
}

// Position in par of the row (t, g, r, c), or 0 if the table has no such row.
// The key includes the group (spec 4.8 G6).
real scalar dwlsmv_dtn_find(struct dwlsmv_parmap scalar par, real scalar t,
                            real scalar g, real scalar r, real scalar c)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == t) {
            if (dwlsmv_dtn_grp(par, i) == g) {
                if (par.row[i] == r) {
                    if (par.col[i] == c) {
                        return(i)
                    }
                }
            }
        }
    }
    return(0)
}

// 1 if two free rows share one fidx (a tie, spec 4.8 G2).
real scalar dwlsmv_dtn_hastie(struct dwlsmv_parmap scalar par)
{
    real colvector f

    f = select(par.fidx, par.free :!= 0)
    if (rows(f) < 2) {
        return(0)
    }
    return(rows(uniqrows(f)) < rows(f))
}

// 1 if the table has rows in a group above 1 AND a beta row free or fixed nonzero.
real scalar dwlsmv_dtn_mgbeta(struct dwlsmv_parmap scalar par)
{
    real scalar i, mg, hasb

    mg   = 0
    hasb = 0
    for (i = 1; i <= rows(par.grp); i++) {
        if (!missing(par.grp[i])) {
            if (par.grp[i] > 1) {
                mg = 1
            }
        }
    }
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 3) {
            if (par.free[i] | par.value[i] != 0) {
                hasb = 1
            }
        }
    }
    return(mg & hasb)
}

real scalar dwlsmv_dtnest(struct dwlsmv_parmap scalar par0,
                          struct dwlsmv_parmap scalar par1,
                          real scalar df0, real scalar df1,
                          string scalar name0, string scalar name1,
                          | string scalar cc0, string scalar can0,
                          string scalar inv0, real matrix K0,
                          string scalar cc1, string scalar can1,
                          string scalar inv1, real matrix K1,
                          real scalar exempt)
{
    real scalar i, j, k, v0, v1, d, nbad, fixed1, ex
    string colvector bad

    if (args() != 6 & args() < 14) {
        errprintf("dwlsmv check 8 (%s in %s): the generated-chain inputs are incomplete.\n",
                  name0, name1)
        _error("dwlsmv check 8: incomplete generated-chain inputs")
    }
    if (args() == 15) {
        exempt = 0
    }

    if (rows(par0.ptype) != par0.nall | rows(par1.ptype) != par1.nall) {
        errprintf("dwlsmv check 8 (%s in %s): a parameter table is malformed.\n",
                  name0, name1)
        _error("dwlsmv check 8: a parameter table is malformed")
    }
    if (missing(df0) | missing(df1)) {
        errprintf("dwlsmv check 8 (%s in %s): a model's e(df_T3) is missing.\n",
                  name0, name1)
        _error("dwlsmv check 8: e(df_T3) is missing")
    }

    // Guard (user, 2026-10-10, spec 4.8 O2, option (e)): DIFFTEST for a pair with
    // groups and a structural path is refused until its own validation gate is met
    // (the beta DIFFTEST follow-up was NOT MET). Groups are counted row by row:
    // max() of an empty grp would be missing, and missing > 1 is true.
    if (dwlsmv_dtn_mgbeta(par0) | dwlsmv_dtn_mgbeta(par1)) {
        errprintf("dwlsmv check 8 (%s in %s): DIFFTEST with groups and a structural path (beta) is\n",
                  name0, name1)
        errprintf("not available in this version: its validation was not met (spec 4.8 O2). Compare\n")
        errprintf("untied slopes with a Wald test on e(b) and e(V) instead.\n")
        exit(3498)
    }

    // The generated-chain exemption (spec check 8, built 2026-10-13): after the
    // beta guard, before the tie refusal. An exempt pair is nested by
    // construction (any two levels from one base); rules 0-2 are not applied,
    // d may be 0 (zero-df pairs, option (c): the caller reports "no test"), and
    // d must still equal the difference in e(df_T3).
    ex = 0
    if (args() >= 14) {
        ex = dwlsmv_dtgenerated(cc0, can0, inv0, K0, cc1, can1, inv1, K1, name0, name1)
    }
    if (ex) {
        d = par1.npar - par0.npar
        if (d < 0) {
            errprintf("dwlsmv check 8 (%s nested in %s): generated levels, but q1 - q0 = %g - %g < 0;\n",
                      name0, name1, par1.npar, par0.npar)
            errprintf("%s is the less restricted model.\n", name0)
            exit(3498)
        }
        if (df0 - df1 != d) {
            errprintf("dwlsmv check 8 (%s nested in %s): the tables give d = %g but e(df_T3)\n",
                      name0, name1, d)
            errprintf("differs by %g (%g in %s, %g in %s).\n", df0 - df1, df0, name0, df1, name1)
            _error("dwlsmv check 8: d disagrees with the difference in e(df_T3)")
        }
        if (args() == 15) {
            exempt = 1
        }
        return(d)
    }

    // Ties (user, 2026-10-12, option (b), fail-closed): the rules below read
    // free/fixed status only, never ties, so a pair tying loadings in H0 and
    // thresholds in H1 would pass though not nested. Refused until a tie rule is
    // validated; generated invariance chains are to come through the exemption.
    if (dwlsmv_dtn_hastie(par0) | dwlsmv_dtn_hastie(par1)) {
        errprintf("dwlsmv check 8 (%s in %s): a model has equality constraints (tied parameters).\n",
                  name0, name1)
        errprintf("Nesting under ties is not checked in this version, so DIFFTEST for such a pair\n")
        errprintf("is not available (spec Module 12 check 8).\n")
        exit(3498)
    }

    // Rule 0: rows matched by index must be the same parameter by name.
    bad = J(0, 1, "")
    for (i = 1; i <= par0.nall; i++) {
        j = dwlsmv_dtn_find(par1, par0.ptype[i], dwlsmv_dtn_grp(par0, i),
                                par0.row[i], par0.col[i])
        if (j > 0) {
            if (par0.name[i] != par1.name[j]) {
                bad = bad \ ("position of " + par0.name[i] + " in " + name0 +
                             " holds " + par1.name[j] + " in " + name1)
            }
        }
    }
    nbad = rows(bad)
    if (nbad > 0) {
        errprintf("dwlsmv check 8 (%s nested in %s): %g matched row(s) differ by name:\n",
                  name0, name1, nbad)
        for (k = 1; k <= min((nbad, 8)); k++) {
            errprintf("    %s\n", bad[k])
        }
        if (nbad > 8) {
            errprintf("    ... and %g more\n", nbad - 8)
        }
        errprintf("The models must declare the same factors and indicators in the same order;\n")
        errprintf("nothing is reordered.\n")
        exit(3498)
    }

    bad = J(0, 1, "")

    // Every row of H0.
    for (i = 1; i <= par0.nall; i++) {
        j = dwlsmv_dtn_find(par1, par0.ptype[i], dwlsmv_dtn_grp(par0, i),
                                par0.row[i], par0.col[i])
        if (par0.free[i]) {
            if (j == 0) {
                bad = bad \ (par0.name[i] + " is free in " + name0 +
                             " and absent from " + name1 + " (fixed at " +
                             strofreal(dwlsmv_dtn_absent(par0.ptype[i])) + ")")
            }
            else {
                if (par1.free[j] == 0) {
                    bad = bad \ (par0.name[i] + " is free in " + name0 +
                                 " and fixed at " + strofreal(par1.value[j]) +
                                 " in " + name1)
                }
            }
        }
        else {
            // Fixed in H0: free in H1 is the restriction; otherwise the values
            // must agree, an absent H1 row reading its type's default (G4).
            v0     = par0.value[i]
            v1     = dwlsmv_dtn_absent(par0.ptype[i])
            fixed1 = 1
            if (j > 0) {
                if (par1.free[j]) {
                    fixed1 = 0
                }
                else {
                    v1 = par1.value[j]
                }
            }
            if (fixed1) {
                if (missing(v0) | missing(v1)) {
                    errprintf("dwlsmv check 8 (%s in %s): %s has a missing fixed value.\n",
                              name0, name1, par0.name[i])
                    _error("dwlsmv check 8: a fixed parameter value is missing")
                }
                if (v0 != v1) {
                    bad = bad \ (par0.name[i] + " is fixed at " + strofreal(v0) +
                                 " in " + name0 + " and at " + strofreal(v1) +
                                 " in " + name1)
                }
            }
        }
    }

    // Rows of H1 that H0 does not have: absent from H0 reads its type's default (G4).
    for (j = 1; j <= par1.nall; j++) {
        i = dwlsmv_dtn_find(par0, par1.ptype[j], dwlsmv_dtn_grp(par1, j),
                                par1.row[j], par1.col[j])
        if (i == 0) {
            if (par1.free[j] == 0) {
                v1 = par1.value[j]
                if (missing(v1)) {
                    errprintf("dwlsmv check 8 (%s in %s): %s has a missing fixed value.\n",
                              name0, name1, par1.name[j])
                    _error("dwlsmv check 8: a fixed parameter value is missing")
                }
                if (v1 != dwlsmv_dtn_absent(par1.ptype[j])) {
                    bad = bad \ (par1.name[j] + " is absent from " + name0 +
                                 " (fixed at " + strofreal(dwlsmv_dtn_absent(par1.ptype[j])) +
                                 ") and fixed at " + strofreal(v1) +
                                 " in " + name1)
                }
            }
        }
    }

    nbad = rows(bad)
    if (nbad > 0) {
        errprintf("dwlsmv check 8 (%s nested in %s): the models are not nested; %g row(s):\n",
                  name0, name1, nbad)
        for (k = 1; k <= min((nbad, 8)); k++) {
            errprintf("    %s\n", bad[k])
        }
        if (nbad > 8) {
            errprintf("    ... and %g more\n", nbad - 8)
        }
        errprintf("Nesting needs every parameter free in %s to be free in %s, and every\n",
                  name0, name1)
        errprintf("parameter fixed in both to be fixed at the same value (a row with no\n")
        errprintf("directive is fixed at 0, a scale (omega) at 1). A pair nested only through a reparameterization,\n")
        errprintf("a different marker for example, is refused: as tabled, neither model is a\n")
        errprintf("restriction of the other.\n")
        exit(3498)
    }

    d = par1.npar - par0.npar
    if (d < 1) {
        errprintf("dwlsmv check 8 (%s nested in %s): q1 - q0 = %g - %g = %g; the less\n",
                  name0, name1, par1.npar, par0.npar, d)
        errprintf("restricted model needs at least one more free parameter.\n")
        exit(3498)
    }
    if (df0 - df1 != d) {
        errprintf("dwlsmv check 8 (%s nested in %s): the tables give d = %g but e(df_T3)\n",
                  name0, name1, d)
        errprintf("differs by %g (%g in %s, %g in %s).\n", df0 - df1, df0, name0, df1, name1)
        exit(3498)
    }

    return(d)
}

end
