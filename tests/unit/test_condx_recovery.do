*! Executed 2026-11-15 by the main session (StataMP 18): passes after three fixes to the draft (a Mata reserved word, a helper redefined after tools/build_mlib.do clears Mata, the var(F) stripe of an exogenous factor). Drafted 2026-11-15 by a read-only agent.
*! D11's marker-metric recovery for a conditional fit (one binary covariate, one group): the owner's decision of 2026-11-15
*! (E-107); derivation draft
*! E-097 (Results C1 to C6).
*!
*! Needs the files of recovery_build/ (dwlsmv_m11_estimate, _rectrigger, _receligible, _rectransform, _recstarts,
*! dwlsmv_cx_admissnote .mata; dwlsmv.ado, dwlsmv.sthlp). Simulated data, fixed seeds, no fixture file. Tools: tools/build_mlib.do.
*!
*! DESIGNS. tm_data is test_condx_marker.do's: N 500, loadings .7 .7 -.7 -.7 (the reverse-keyed trap), factor shifted .5 by x,
*! direct effect +.3 of x on y4*, 3 categories (thresholds -.5, .5). Its seeds 20262415 20262516 20262617 20262718 are the four on
*! which the conditional marker fit did not converge before this change (probe_inadm.log) and unit variance did.
*!
*! Part M, Mata level, no fit: the map against sigma(theta) (R1), recstarts round trip (R2), eligibility (R3), trigger (R4).
*!   R1  a two-factor model with beta, gamma on both factors, a free and a fixed kappa (the verification model of
*!       recovery_verify.R): for 200 random unit-variance points, dwlsmv_m11_rectransform() gives a marker point with the
*!       SAME sigma(theta) as the unit-variance point, max |gap| <= 1e-12 (the floating-point tier of the derivation's
*!       Result C1; chosen before measuring). This is the deciding criterion for the map; the independent comparator is
*!       dwls_model() on two different tables, so the code under test (the map) shares nothing with it.
*!   R2  dwlsmv_m11_recstarts() then dwls_identify(): every free non-threshold row of the restarted table equals the transformed
*!       value EXACTLY (gamma and kappa lines carry their start; no second-start substitution); fixed rows are untouched.
*!   R3  dwlsmv_m11_receligible(): fixed gamma 0.5 and fixed gamma -0.2 ineligible, each naming its row; fixed gamma 0 eligible;
*!       fixed kappa of any value eligible; the old rule (fixed beta 0.3) still ineligible.
*!   R4  dwlsmv_m11_rectrigger(): with Bd = J(0,0,.) the trigger is non-convergence or an admissibility flag, and does not abort;
*!       a nonempty Bd with row 3 flagged still fires (an ordinary fit is unchanged).
*! Part C, through the command:
*!   C1  the four trap seeds, marker metric: e(recovery) = recovered, e(converged) = 1, e(V) 15 x 15, T3 and p posted, the recovery
*!       note printed, no "second identification start", no first-fit retry note; replay agrees.
*!   C2  the recovered e(b) satisfies the map against a SEPARATE metric(unitvariance) fit of the same data (lambda/c, c*gamma,
*!       kappa and tau unchanged, psi = c^2), tolerance 1e-5 relative on e(b) and 1e-6 relative on T3, df equal. These two
*!       are chosen before measuring (2026-11-15); if they fail, the observed gap is reported, they are not loosened to pass.
*!       The same T3 also says the recovery reused Stage 1 unchanged (a different s or Gamma gives a different T3).
*!   C3  a fixed NONZERO gamma, (F <- x@.3), first fit triggered: ineligible, with the note naming gamma[F,x]; not a refit.
*!   C4  a fixed gamma 0, (F <- x@0), and a fixed kappa, (y4 <- x@.3): never ineligible; recovered or failed with its note.
*!   C5  metric(unitvariance) on the trap data: e(recovery) none, no marker-metric note.
*!   C6  a converged admissible conditional marker fit (test_condx_output.do's data): e(recovery) none, no note, as before.
*!   C7  an ORDINARY fit (no covariate) on the same data: recovery unchanged. If it is recovered, its e(b) satisfies the same map
*!       (without gamma and kappa) against a unit-variance fit; no conditional text is printed.
*!   C8  the header of a nonconverged conditional fit says which recovery outcome applies (ineligible, failed, other): hand-set
*!       e(recovery) on a posted conditional fit (no real fit is forced to fail; spec D8 has no iteration-limit option).
*! Part T, Mata level, the bread cause of the conditional trigger (drafted 2026-11-15 by a read-only agent; executed by the main
*! session the same day, StataMP 18: passes on its first run; T1 recovered, T2 ineligible, twin ratio 1.9e-13; mutants T-M1..T-M5 killed):
*!   T0  What cannot be built. The bread cause fires only on a first fit whose W^1/2 Delta is numerically rank deficient at
*!       theta-hat. For this class (one factor, one covariate) a scale-unidentified table (marker loading, gamma and psi all
*!       free) is refused by dwls_identify() at theta0 and at theta2 (dwlsmv_id_rank, rc 3498) before any fit, and Module 3's
*!       marker line fixes the loading, so no sem-style text produces it; a deficiency that is generic is refused at the
*!       start, and one that is point-specific needs gamma-hat and psi-hat both within 1e-8 of 0 (see the report). So the
*!       first fit is real and only the bread handed to the recovery is the scale twin's (probe_cx_delta.do, ratio ~4e-13).
*!   T1  Population theta0 (lambda 1 .8 .6, psi .5, gamma .3), marker, gamma free. CONTROL through dwlsmv_m11_estimate: rec none,
*!       converged, no admissibility flag, Bd row 3 clear, ratio > 1e-4. The twin's Bd alone fires the trigger with exactly the bread
*!       cause; dwlsmv_m11_recover() with it returns recovered, the note names the bread cause and the unit-variance restart, not
*!       the other two causes, and the recovered theta is theta0 (1e-5).
*!   T2  The same population and Bd, gamma[F,x] fixed at .3: ineligible, the note names the bread cause and "gamma[F,x] is fixed at",
*!       and no restart is mentioned.
*!   T1/T2 certify the composition trigger(bread) -> eligibility -> recovery on a real first fit and a real Bd flag; they do not
*!   certify that production reaches a deficient bread on a conditional fit (nothing here shows that it does).
*!   Mutants (each must turn this file red):
*!   T-M1  rectrigger: `Bd[3, 1] == 1` -> `Bd[3, 1] == 0`.           Killed by T1/T2 (control fires; the twin's bread does not) and R4.
*!   T-M2  receligible E4: delete the `for` loop over ptype 8.       Killed by T2 (recovered, not ineligible) and R3, C3.
*!   T-M3  estimate: `bread_cx = 1` -> `bread_cx = 0`.               Killed by T1/T2 control (Bd is 0 x 0: rows(Bd) < 4) and diag check 6.
*!   T-M4  rectrigger: the bread phrase "could not determine all parameters reliably (the information in the fit was insufficient)"
*!         edited (any word).                                         Killed by the exact-text assertion in T1/T2 and R4.
*!   T-M5  recover: pass `Bd` unused, i.e. `trig = dwlsmv_m11_rectrigger(fit, A, J(0, 0, .))`. Killed by T1/T2 (rec "none", no note).
*!
*! WHAT THIS CERTIFIES AND DOES NOT. It certifies the map's algebra on random points (R1), that recovery is attempted, succeeds on
*! the four named seeds and reproduces the unit-variance fit through the map (C1, C2), the eligibility rule and the trigger's
*! composition. It does NOT certify: a rate of recovery or of the trap in conditional mode (four seeds, one design); recovery at
*! G > 1 or with several covariates; Module 10 diagnostics (their own file, test_condx_diag.do). The bread-rank term is in the conditional
*! trigger since 2026-11-15 (dwlsmv_m11_estimate.mata, bread_cx = 1; dwlsmv_m11_rectrigger.mata).
*!
*! MUTANTS (edit under src/, run, revert; each must turn this file red). Changes are to the proposed files.
*!   R-M1  rectransform, ptype 8: `u * c[r]` -> `u`.                Killed by R1 (sigma gap ~1) and C2 (F:x column).
*!   R-M2  rectransform, ptype 8: `u * c[r]` -> `u / c[r]`.        Killed by R1 and C2.
*!   R-M3  rectransform, ptype 9: `u` -> `u * c[r]`.               Killed by R1 and C2 (y4:x column).
*!   R-M4  rectransform: delete the ptype 8 and ptype 9 branches (the unknown-ptype error fires). Killed by R1 (abort) and C1.
*!   R-M5  rectransform, ptype 2: `u * c[r] * c[cc]` -> `u * c[r]`. Killed by R1 and C2 (var(e.F)); the ordinary path by C7.
*!   R-M6  recstarts: delete the `gamma` branch. Killed by R2 (restarted value is the default 0, not the transformed value) and by
*!         C1 (the note says the second identification start was used).
*!   R-M7  recstarts: delete the `kappa` branch. Killed by R2 and C1 likewise.
*!   R-M8  recstarts: read the covariate list from the `ind` line (cvn = ind). Killed by R2.
*!   R-M9  receligible: delete the E4 loop. Killed by R3 (fixed gamma 0.5 eligible) and C3.
*!   R-M10 receligible, E4: `par.value[i] != 0` -> `par.value[i] == 0`. Killed by R3 (gamma fixed 0 ineligible) and C4.
*!   R-M11 receligible, E4: `par.ptype[i] == 8` -> `par.ptype[i] == 8 | par.ptype[i] == 9`. Killed by R3 (fixed kappa) and C4.
*!   R-M12 estimate: `if (metric == "marker") { fit = dwlsmv_fit_noretry` -> `fit = dwls_fit(st, par)` (retry back for every fit).
*!         Killed by C1 (retry note; and the trap seeds' first fit may converge, so recovery is none).
*!   R-M13 estimate: guard the recovery with `& !cx`. Killed by C1.
*!   R-M14 rectrigger: delete the `rows(Bd) >= 3` guard. Killed by R4 (Bd empty aborts).
*!   R-M15 rectrigger: `Bd[3, 1] == 1` -> `Bd[3, 1] == 0`. Killed by R4 (a nonempty Bd with row 3 flagged must fire) and by
*!         test_marker_recovery.do.
*!   R-M16 admissnote: restore the `if (metric == "marker") {` paragraph "No recovery refit was attempted for a model with a covariate".
*!         Killed by C1, C3 (hasnt).
*!   R-M17 dwlsmv.ado header: swap the `ineligible` and `failed` branches' text. Killed by C8 and C3 (when its fit is nonconverged).
*!   R-M18 dwlsmv.ado header: restore the old line "No recovery refit is attempted". Killed by C8, C3.
*!   NOT KILLED BY THIS FILE (stated, not hidden): (a) a mutant in dwlsmv_m11_recrun.mata (not changed). A mutant that turns the
*!   bread-rank term OFF for a conditional fit (`bread_cx = 0`) is now killed by Part T's control (T-M3) and test_condx_diag.do.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

global TR_NFAIL 0

* =====================================================================================================
* Part M: Mata level. Sources loaded, types first (as test_condx_marker.do Part A).
* =====================================================================================================
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_condx_recovery: `f' failed to load"
        exit 9
    }
}

mata:

real scalar tr_fail(string scalar what)
{
    printf("{err}test_condx_recovery: %s\n", what)
    return(1)
}

// the model of recovery_verify.R, as canonical text lines (marker metric). Items 1..6, K = (2,3,2,3,2,2): 8 thresholds.
// F: y1 (marker), y2, y3; G: y3, y4 (marker), y5, y6. beta G <- F free, psi G F fixed 0, gamma free on both,
// kappa y2 fixed .15, kappa y5 free.
string rowvector tr_base()
{
    return(("factor F G", "ind y1 y2 y3 y4 y5 y6",
            "lambda F y2 free", "lambda F y3 free", "lambda G y3 free", "lambda G y5 free", "lambda G y6 free",
            "beta G F free", "psi F F free", "psi G G free", "psi G F fixed 0",
            "metric marker F y1", "metric marker G y4",
            "covariate x", "gamma F x free", "gamma G x free",
            "kappa y2 x fixed 0.15", "kappa y5 x free"))
}

// replace the line that starts with `pre' by `new' (new == "" deletes it); a line that is not there is an error of the test
string rowvector tr_swap(string rowvector base, string scalar pre, string scalar repl)
{
    string rowvector out
    real scalar i, hit

    out = J(1, 0, "")
    hit = 0
    for (i = 1; i <= cols(base); i++) {
        if (substr(base[i], 1, strlen(pre)) == pre) {
            hit = hit + 1
            if (repl != "") {
                out = out, repl
            }
        }
        else {
            out = out, base[i]
        }
    }
    if (hit != 1) {
        _error("test_condx_recovery: tr_swap found " + strofreal(hit) + " lines starting with " + pre)
    }
    return(out)
}

string scalar tr_join(string rowvector lines)
{
    string scalar canon
    real scalar i

    canon = ""
    for (i = 1; i <= cols(lines); i++) {
        if (i > 1) {
            canon = canon + ";"
        }
        canon = canon + lines[i]
    }
    return(canon)
}

// a unit-variance point: thresholds from thr (in table order), every other free row drawn from the stream
real colvector tr_uvpoint(struct dwlsmv_parmap scalar paru, real colvector thr)
{
    real colvector th
    real scalar i, t, u

    th = J(paru.npar, 1, .)
    t  = 0
    for (i = 1; i <= paru.nall; i++) {
        if (paru.free[i]) {
            u = runiform(1, 1)
            if (paru.ptype[i] == 4) {
                t = t + 1
                th[paru.fidx[i]] = thr[t]
            }
            else if (paru.ptype[i] == 1) {
                th[paru.fidx[i]] = (runiform(1, 1) < .35 ? -1 : 1) * (.25 + .35 * u)
            }
            else if (paru.ptype[i] == 3) {
                th[paru.fidx[i]] = -.3 + .6 * u
            }
            else if (paru.ptype[i] == 8) {
                th[paru.fidx[i]] = -.8 + 1.6 * u
            }
            else if (paru.ptype[i] == 9) {
                th[paru.fidx[i]] = -.4 + .8 * u
            }
            else {
                _error("test_condx_recovery: tr_uvpoint met a free row of an unexpected ptype")
            }
        }
    }
    return(th)
}

// R1 and R2
real scalar tr_maprt()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar parm, paru, park
    struct dwlsmv_statmap vector mv
    real colvector K, thr, thu, thm, c, marker, su, sm
    string colvector fac, reasons
    string scalar canon, uvtext, kt
    real scalar nf, rep, rc, gap, worst, i, ok, nneg, nrt, worst2

    nf    = 0
    K     = (2 \ 3 \ 2 \ 3 \ 2 \ 2)
    thr   = (0.2 \ -0.5 \ 0.4 \ -0.1 \ -0.3 \ 0.6 \ 0.0 \ 0.5)
    canon = tr_join(tr_base())
    st.map = dwlsmv_statmap_build(K, 1)
    mv     = dwlsmv_statmap(1, 1)
    mv[1]  = st.map
    parm   = dwls_identify(subinstr(canon, ";", char(10)), st.map)
    uvtext = dwlsmv_m11_recuvtext(canon)
    paru   = dwls_identify(subinstr(uvtext, ";", char(10)), st.map)
    fac    = ("F" \ "G")
    if (!dwlsmv_m11_receligible(parm, fac, marker, reasons)) {
        return(tr_fail("R1: the verification model is not eligible: " + invtokens(reasons')))
    }
    if (parm.npar != paru.npar) {
        return(tr_fail("R1: the two tables have " + strofreal(parm.npar) + " and " + strofreal(paru.npar) + " free rows"))
    }

    rseed(20261115)
    worst = 0
    worst2 = 0
    nneg  = 0
    nrt   = 0
    for (rep = 1; rep <= 200; rep++) {
        thu = tr_uvpoint(paru, thr)
        rc  = dwlsmv_m11_rectransform(parm, paru, thu, marker, thm, c)
        if (rc != 0) {
            return(nf + tr_fail("R1: rectransform returned " + strofreal(rc) + " at draw " + strofreal(rep)))
        }
        if (any(c :< 0)) {
            nneg = nneg + 1
        }
        su  = dwls_model(thu, paru, mv)
        sm  = dwls_model(thm, parm, mv)
        gap = max(abs(su - sm))
        if (gap > worst) {
            worst = gap
        }
        // R2 on the first 20 draws: the restarted table carries the transformed value exactly
        if (rep <= 20) {
            kt   = dwlsmv_m11_recstarts(canon, parm, thm)
            park = dwls_identify(subinstr(kt, ";", char(10)), st.map)
            ok   = 1
            if (park.nall != parm.nall) {
                ok = 0
            }
            else {
                for (i = 1; i <= park.nall; i++) {
                    if (park.free[i] != parm.free[i]) {
                        ok = 0
                    }
                    else if (park.free[i]) {
                        if (park.ptype[i] != 4) {
                            if (park.value[i] != thm[park.fidx[i]]) {
                                ok = 0
                            }
                        }
                    }
                    else {
                        if (park.value[i] != parm.value[i]) {
                            ok = 0
                        }
                    }
                }
            }
            nrt = nrt + 1
            if (!ok) {
                nf = nf + tr_fail("R2: draw " + strofreal(rep) + ": the restarted table does not carry the transformed values (a gamma or kappa start missing, or dwls_identify used its second point)")
            }
        }
    }
    if (worst > 1e-12) {
        nf = nf + tr_fail("R1: max |sigma_M(h(theta_U)) - sigma_U(theta_U)| = " + strofreal(worst, "%9.2e") + " over 200 draws, tolerance 1e-12")
    }
    if (nneg < 20) {
        nf = nf + tr_fail("R1: only " + strofreal(nneg) + " of 200 draws had a negative c_f; the sign path is not exercised")
    }
    printf("{txt}  R1: max sigma gap %s over 200 draws (%f with a negative c_f); R2: %f restarted tables checked\n",
           strofreal(worst, "%9.2e"), nneg, nrt)
    return(nf)
}

// R3
real scalar tr_elig()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar parm
    real colvector K, marker
    string colvector fac, reasons
    string scalar canon, why
    real scalar nf, ok

    nf     = 0
    K      = (2 \ 3 \ 2 \ 3 \ 2 \ 2)
    st.map = dwlsmv_statmap_build(K, 1)
    fac    = ("F" \ "G")

    // the base model: eligible (a fixed kappa .15 is in it)
    parm = dwls_identify(subinstr(tr_join(tr_base()), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (!ok) {
        nf = nf + tr_fail("R3 base: not eligible, reasons: " + invtokens(reasons'))
    }

    // fixed gamma 0.5 on F: ineligible, exactly one reason, naming gamma[F,x]; marker zeroed
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_base(), "gamma F x", "gamma F x fixed 0.5")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (ok | rows(reasons) != 1) {
        nf = nf + tr_fail("R3 gamma F fixed 0.5: eligible " + strofreal(ok) + ", " + strofreal(rows(reasons)) + " reasons, wanted ineligible with 1")
    }
    else {
        why = reasons[1]
        if (strpos(why, "gamma[F,x] is fixed at") == 0 | strpos(why, ", not 0") == 0) {
            nf = nf + tr_fail("R3 gamma F fixed 0.5: the reason is not 'gamma[F,x] is fixed at ..., not 0': " + why)
        }
        if (any(marker :!= 0)) {
            nf = nf + tr_fail("R3 gamma F fixed 0.5: marker is not zeroed on an ineligible model")
        }
    }

    // fixed gamma -0.2 on G: ineligible, names gamma[G,x]
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_base(), "gamma G x", "gamma G x fixed -0.2")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (ok | rows(reasons) != 1) {
        nf = nf + tr_fail("R3 gamma G fixed -0.2: eligible " + strofreal(ok) + ", " + strofreal(rows(reasons)) + " reasons, wanted ineligible with 1")
    }
    else if (strpos(reasons[1], "gamma[G,x] is fixed at") == 0) {
        nf = nf + tr_fail("R3 gamma G fixed -0.2: the reason does not name gamma[G,x]: " + reasons[1])
    }

    // both fixed: two reasons (all reasons are collected)
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_swap(tr_base(), "gamma G x", "gamma G x fixed -0.2"), "gamma F x", "gamma F x fixed 0.5")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (ok | rows(reasons) != 2) {
        nf = nf + tr_fail("R3 both gammas fixed: eligible " + strofreal(ok) + ", " + strofreal(rows(reasons)) + " reasons, wanted ineligible with 2")
    }

    // fixed gamma 0 on F (the other free): eligible
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_base(), "gamma F x", "gamma F x fixed 0")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (!ok) {
        nf = nf + tr_fail("R3 gamma F fixed 0: not eligible, reasons: " + invtokens(reasons'))
    }

    // fixed kappa of any value: eligible (the free one fixed at .4, the fixed one at -.7 and at 0)
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_base(), "kappa y5 x", "kappa y5 x fixed 0.4")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (!ok) {
        nf = nf + tr_fail("R3 kappa y5 fixed 0.4: not eligible, reasons: " + invtokens(reasons'))
    }
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_base(), "kappa y2 x", "kappa y2 x fixed -0.7")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (!ok) {
        nf = nf + tr_fail("R3 kappa y2 fixed -0.7: not eligible, reasons: " + invtokens(reasons'))
    }
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_base(), "kappa y2 x", "kappa y2 x fixed 0")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (!ok) {
        nf = nf + tr_fail("R3 kappa y2 fixed 0: not eligible, reasons: " + invtokens(reasons'))
    }

    // the old rule is intact: a fixed nonzero beta is still a reason, and only that one
    parm = dwls_identify(subinstr(tr_join(tr_swap(tr_base(), "beta G F", "beta G F fixed 0.3")), ";", char(10)), st.map)
    ok   = dwlsmv_m11_receligible(parm, fac, marker, reasons)
    if (ok | rows(reasons) != 1) {
        nf = nf + tr_fail("R3 beta fixed 0.3: eligible " + strofreal(ok) + ", " + strofreal(rows(reasons)) + " reasons, wanted ineligible with 1")
    }
    return(nf)
}

// R4
real scalar tr_trig()
{
    struct dwlsmv_fit scalar f
    real matrix A, Bd
    string scalar r
    real scalar nf

    nf = 0
    A  = J(6, 4, 0)
    Bd = J(0, 0, .)

    f.converged = 1
    r = dwlsmv_m11_rectrigger(f, A, Bd)
    if (r != "") {
        nf = nf + tr_fail("R4: converged, admissible, empty Bd fires: " + r)
    }
    f.converged = 0
    r = dwlsmv_m11_rectrigger(f, A, Bd)
    if (r != "did not converge") {
        nf = nf + tr_fail("R4: not converged, empty Bd gives: " + r)
    }
    f.converged = 1
    A[2, 1] = 1
    r = dwlsmv_m11_rectrigger(f, A, Bd)
    if (r != "was not admissible") {
        nf = nf + tr_fail("R4: an admissibility flag, empty Bd gives: " + r)
    }
    f.converged = 0
    r = dwlsmv_m11_rectrigger(f, A, Bd)
    if (r != "did not converge and was not admissible") {
        nf = nf + tr_fail("R4: both, empty Bd gives: " + r)
    }
    // a missing admissibility flag does not fire
    A = J(6, 4, 0)
    A[2, 1] = .
    f.converged = 1
    r = dwlsmv_m11_rectrigger(f, A, Bd)
    if (r != "") {
        nf = nf + tr_fail("R4: a missing flag fires: " + r)
    }
    // an ordinary Bd (nonempty): row 3 flagged still fires; unflagged does not
    Bd = J(4, 3, 0)
    r = dwlsmv_m11_rectrigger(f, A, Bd)
    if (r != "") {
        nf = nf + tr_fail("R4: a nonempty Bd with row 3 clear fires: " + r)
    }
    Bd[3, 1] = 1
    r = dwlsmv_m11_rectrigger(f, A, Bd)
    if (strpos(r, "could not determine all parameters reliably") == 0) {
        nf = nf + tr_fail("R4: a nonempty Bd with row 3 flagged does not fire: " + r)
    }
    return(nf)
}

// ---------------------------------------------------------------------------------------------------------------
// Part T (executed 2026-11-15; see the header): the bread-rank cause of the conditional D11 trigger, and eligibility E4 when it fires.
// The command cannot reach a conditional marker fit that is bread-deficient at its first fit (see the header, T0), so
// the first fit is REAL (dwlsmv_m11_estimate on an exact conditional population, Gamma = I, as tm_bcase) and only the
// bread diagnostic handed to dwlsmv_m11_recover() is the one a scale twin gives (dwls_delta() of the same model with
// the marker loading freed, dwlsmv_breaddiag() of it: a genuine row 3 flag, ratio about 1e-13). The call is the
// one dwlsmv_m11_estimate makes, `rec = dwlsmv_m11_recover(canon, par, fit, A, Bd, parr, fitr, rnote)`, with Bd replaced.
// ---------------------------------------------------------------------------------------------------------------

// the twin's table: the marker loading (first fixed ptype 1 row with value 1) freed, fidx renumbered in table order;
// thnew = theta in the new numbering (that loading at its fixed value 1). As probe_cx_delta.do's pb_twin.
struct dwlsmv_parmap scalar trx_twin(struct dwlsmv_parmap scalar par, real colvector th, real colvector thnew)
{
    struct dwlsmv_parmap scalar p2
    real colvector vals
    real scalar i, r1, k

    p2   = par
    vals = J(par.nall, 1, 0)
    r1   = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            vals[i] = th[par.fidx[i]]
        }
        else {
            vals[i] = par.value[i]
            if (r1 == 0) {
                if (par.ptype[i] == 1) {
                    if (par.value[i] == 1) {
                        r1 = i
                    }
                }
            }
        }
    }
    if (r1 == 0) {
        _error("test_condx_recovery: trx_twin found no fixed marker loading")
    }
    p2.free[r1] = 1
    k = 0
    for (i = 1; i <= par.nall; i++) {
        p2.fidx[i] = 0
        if (p2.free[i]) {
            k = k + 1
            p2.fidx[i] = k
        }
    }
    p2.npar = k
    thnew   = J(k, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        if (p2.free[i]) {
            thnew[p2.fidx[i]] = vals[i]
        }
    }
    return(p2)
}

// theta of a table by ptype: thresholds from tau, free loadings (ptype 1) from lam in table order, gamma (8), psi (2)
real colvector trx_theta(struct dwlsmv_parmap scalar par, real colvector tau, real colvector lam, real scalar gam, real scalar psi)
{
    real colvector th
    real scalar i, t, l

    th = J(par.npar, 1, .)
    t  = 0
    l  = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            if (par.ptype[i] == 4) {
                t = t + 1
                th[par.fidx[i]] = tau[t]
            }
            else if (par.ptype[i] == 1) {
                l = l + 1
                th[par.fidx[i]] = lam[l]
            }
            else if (par.ptype[i] == 8) {
                th[par.fidx[i]] = gam
            }
            else if (par.ptype[i] == 2) {
                th[par.fidx[i]] = psi
            }
            else {
                _error("test_condx_recovery: trx_theta met a free row of an unexpected ptype")
            }
        }
    }
    if (t != rows(tau) | l != rows(lam)) {
        _error("test_condx_recovery: trx_theta: the table does not have the expected threshold and loading rows")
    }
    return(th)
}

real scalar trx_ck(string scalar txt, string scalar ph, real scalar want, string scalar tag)
{
    if ((strpos(txt, ph) > 0) != want) {
        return(tr_fail(tag + ": the note " + (want ? "lacks" : "contains") + ": " + ph))
    }
    return(0)
}

// one case. gfixed 0: gamma[F,x] free ("T1"); gfixed 1: gamma[F,x] fixed at .3, its population value ("T2").
// The population is the SAME in both (s = sigma(theta0), theta0: lambda = (1, .8, .6), psi = .5, gamma = .3, so every
// residual variance is positive and the first fit is admissible), so the only difference between the cases is the table.
real scalar trx_case(real scalar gfixed, string scalar tag)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar parF, par, parr, p2
    struct dwlsmv_fit    scalar fit, fitr
    struct dwlsmv_statmap vector mv
    real matrix A, Bd, BdT, D2
    real colvector K, tau, lam, thF, th2
    string rowvector linesF, lines
    string scalar canon, canonF, rec, rnote, rec2, note2, trig, want
    real scalar nf, D
    pointer() scalar pst

    nf     = 0
    K      = (3 \ 3 \ 2)
    tau    = (-0.5 \ 0.5 \ -0.3 \ 0.6 \ 0.2)
    lam    = (0.8 \ 0.6)
    linesF = ("factor F", "ind y1 y2 y3", "lambda F y2 free", "lambda F y3 free", "psi F F free",
              "metric marker F y1", "covariate x", "gamma F x free")
    lines  = linesF
    if (gfixed) {
        lines = tr_swap(linesF, "gamma F x", "gamma F x fixed 0.3")
    }
    canonF = tr_join(linesF)
    canon  = tr_join(lines)

    st.map = dwlsmv_statmap_build(K, 1)
    D      = st.map.D
    mv     = dwlsmv_statmap(1, 1)
    mv[1]  = st.map
    parF   = dwls_identify(subinstr(canonF, ";", char(10)), st.map)
    thF    = trx_theta(parF, tau, lam, 0.3, 0.5)
    st.s      = dwls_model(thF, parF, mv)
    st.Gamma  = I(D)
    st.N      = 1000
    st.psibar = J(D, 1, 0)
    st.bound  = J(st.map.P, 1, 0)

    // the bread diagnostic of the scale twin of the free-gamma model at theta0 (with gamma fixed nonzero the model has no
    // scale null direction at all, so the same Bd is used in both cases: it is the diagnostic of the free-gamma analogue)
    th2 = J(0, 1, 0)
    p2  = trx_twin(parF, thF, th2)
    D2  = dwls_delta(th2, p2, mv)
    BdT = dwlsmv_breaddiag(D2, st.Gamma)
    // preconditions of the construction, not of the code under test
    if (BdT[3, 1] != 1 | BdT[2, 1] != 0 | !(BdT[4, 2] < 1e-10)) {
        return(nf + tr_fail(tag + ": the scale twin's bread is not flagged by row 3 alone (row 3 " + strofreal(BdT[3, 1]) +
                            ", row 2 " + strofreal(BdT[2, 1]) + ", ratio " + strofreal(BdT[4, 2], "%9.2e") + "): the construction failed"))
    }

    // the external D11's recovery reads (dwlsmv_m11_recrun), as dwlsmv_m11_stage1() leaves it; as tm_bcase
    rmexternal("dwlsmv_m11_st")
    rmexternal("dwlsmv_m11_rec")
    pst = crexternal("dwlsmv_m11_st")
    *pst = st

    // CONTROL: the real first fit through dwlsmv_m11_estimate: converged, admissible, healthy bread, no recovery
    dwlsmv_m11_estimate(canon, "marker", st, par, fit, A, Bd, rec, rnote)
    if (fit.converged != 1) {
        nf = nf + tr_fail(tag + " control: the exact-population first fit did not converge")
    }
    if (rec != "none" | rnote != "") {
        nf = nf + tr_fail(tag + " control: e(recovery) is " + rec + " on a converged, admissible fit with a healthy bread; note: " + rnote)
    }
    if (rows(Bd) < 4) {
        rmexternal("dwlsmv_m11_st")
        rmexternal("dwlsmv_m11_rec")
        return(nf + tr_fail(tag + " control: Bd has " + strofreal(rows(Bd)) + " rows (the conditional bread diagnostic is not computed)"))
    }
    if (Bd[3, 1] != 0 | !(Bd[4, 2] > 1e-4)) {
        nf = nf + tr_fail(tag + " control: the healthy bread is flagged or its ratio " + strofreal(Bd[4, 2], "%9.2e") + " is not above 1e-4")
    }
    if (sum(A[., 1] :== 1) != 0) {
        nf = nf + tr_fail(tag + " control: the first fit is flagged inadmissible, so the bread cause is not the only one")
    }

    // THE TRIGGER with the twin's Bd: the bread cause and no other (exact text)
    want = "could not determine all parameters reliably (the information in the fit was insufficient)"
    trig = dwlsmv_m11_rectrigger(fit, A, BdT)
    if (trig != want) {
        nf = nf + tr_fail(tag + ": the trigger with the twin's bread says '" + trig + "', wanted exactly '" + want + "'")
    }

    // THE RECOVERY CALL, as dwlsmv_m11_estimate makes it, with the twin's Bd
    rec2 = dwlsmv_m11_recover(canon, par, fit, A, BdT, parr, fitr, note2)
    rmexternal("dwlsmv_m11_st")
    rmexternal("dwlsmv_m11_rec")
    nf = nf + trx_ck(note2, "could not determine all parameters reliably", 1, tag)
    nf = nf + trx_ck(note2, "(the information in the fit was insufficient)", 1, tag)
    nf = nf + trx_ck(note2, "did not converge", 0, tag)
    nf = nf + trx_ck(note2, "was not admissible", 0, tag)
    if (!gfixed) {
        // T1: eligible; the unit-variance table is identified (only the twin was not), so the restart succeeds
        if (rec2 != "recovered") {
            nf = nf + tr_fail(tag + ": e(recovery) is " + rec2 + ", wanted recovered; note: " + note2)
        }
        else {
            nf = nf + trx_ck(note2, "restarted from the unit-variance solution", 1, tag)
            nf = nf + trx_ck(note2, "second identification start", 0, tag)
            nf = nf + trx_ck(note2, "not eligible", 0, tag)
            if (fitr.converged != 1) {
                nf = nf + tr_fail(tag + ": recovered, but the recovered fit is not converged")
            }
            // exact population: the recovered marker solution is theta0 (1e-5, chosen before measuring; a diagnostic of the
            // map through the command's own recovery, not the deciding criterion, which is the status and the note)
            if (rows(fitr.theta) != rows(thF)) {
                nf = nf + tr_fail(tag + ": recovered theta has " + strofreal(rows(fitr.theta)) + " rows, wanted " + strofreal(rows(thF)))
            }
            else if (max(abs(fitr.theta - thF)) > 1e-5) {
                nf = nf + tr_fail(tag + ": recovered theta is " + strofreal(max(abs(fitr.theta - thF)), "%9.2e") + " from the population theta0 (1e-5)")
            }
        }
    }
    else {
        // T2: gamma[F,x] fixed at .3 is ineligible (E4); the bread cause names the trigger, the note names the reason
        if (rec2 != "ineligible") {
            nf = nf + tr_fail(tag + ": e(recovery) is " + rec2 + ", wanted ineligible; note: " + note2)
        }
        nf = nf + trx_ck(note2, "the model is not eligible for it:", 1, tag)
        nf = nf + trx_ck(note2, "gamma[F,x] is fixed at", 1, tag)
        nf = nf + trx_ck(note2, "restarted from the unit-variance solution", 0, tag)
        nf = nf + trx_ck(note2, "Recovery attempt failed", 0, tag)
    }
    printf("{txt}  %s: control %s (ratio %s); twin ratio %s; recovery %s\n", tag, rec, strofreal(Bd[4, 2], "%9.2e"),
           strofreal(BdT[4, 2], "%9.2e"), rec2)
    return(nf)
}

// Part T entry: T1 gamma free, T2 gamma fixed .3
real scalar tr_cxrec()
{
    return(trx_case(0, "T1 gamma free") + trx_case(1, "T2 gamma fixed .3"))
}

end

foreach t in maprt elig trig cxrec {
    capture noisily mata: st_local("n", strofreal(tr_`t'()))
    if _rc {
        display as error "test_condx_recovery: Part M `t' aborted, rc " _rc
        global TR_NFAIL = ${TR_NFAIL} + 1
    }
    else {
        display as text "  Part M `t': `n' failure(s)"
        global TR_NFAIL = ${TR_NFAIL} + `n'
    }
}

* =====================================================================================================
* Part C: through the command
* =====================================================================================================
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

mata:

// tools/build_mlib.do clears Mata, so Part A's tr_fail() is defined again here (main session, 2026-11-15)
real scalar tr_fail(string scalar what)
{
    printf("{err}test_condx_recovery: %s\n", what)
    return(1)
}

string scalar tr_txt(string scalar fn)
{
    return(invtokens(cat(fn)', " "))
}

real scalar tr_has(string scalar fn, string scalar ph)
{
    return(strpos(tr_txt(fn), ph) > 0)
}

// the map of the derivation, checked on e(b): bm = the marker fit's e(b), bu = the unit-variance fit's, same stripe.
// One factor F, indicators y1..y4, covariate x. c = bu[y1:F]. Returns the number of columns off by more than tol
// (relative to max(1, |want|)), and prints the worst. The rule by column: thresholds and kappa (eq y?, name x) unchanged;
// gamma (eq F, name x) times c; loadings (name F) over c; /var(e.F) = c^2 times the unit-variance value.
real scalar tr_map(string scalar bm, string scalar bu, real scalar tol, string scalar tag)
{
    real rowvector m, u
    string matrix sm, su
    real scalar j, c, want, nbad, worst, rel, k

    m  = st_matrix(bm)
    u  = st_matrix(bu)
    sm = st_matrixcolstripe(bm)
    su = st_matrixcolstripe(bu)
    if (cols(m) != cols(u) | rows(sm) != rows(su) | any(sm :!= su)) {
        return(tr_fail(tag + ": the two fits do not have the same e(b) stripe"))
    }
    c = .
    for (j = 1; j <= cols(m); j++) {
        if (sm[j, 1] == "y1" & sm[j, 2] == "F") {
            c = u[j]
        }
    }
    if (missing(c) | c == 0) {
        return(tr_fail(tag + ": c_F, the marker's unit-variance loading, is " + strofreal(c)))
    }
    nbad  = 0
    worst = 0
    for (j = 1; j <= cols(m); j++) {
        k = 1
        if (substr(sm[j, 1], 1, 1) == "/" & substr(sm[j, 2], 1, 3) == "cut") {
            want = u[j]
        }
        else if (sm[j, 1] == "/" & (sm[j, 2] == "var(e.F)" | sm[j, 2] == "var(F)")) {
            // var(e.F) for a factor x predicts, var(F) for C7's exogenous factor (main session, 2026-11-15)
            want = u[j] * c^2
        }
        else if (sm[j, 1] == "F" & sm[j, 2] == "x") {
            want = u[j] * c
        }
        else if (sm[j, 2] == "x") {
            want = u[j]
        }
        else if (sm[j, 2] == "F") {
            want = u[j] / c
        }
        else {
            return(tr_fail(tag + ": no map rule for the column " + sm[j, 1] + ":" + sm[j, 2]))
        }
        rel = abs(m[j] - want) / max((1, abs(want)))
        if (rel > worst) {
            worst = rel
        }
        if (rel > tol) {
            nbad = nbad + 1
            printf("{err}test_condx_recovery: %s: column %s:%s is %s, the map gives %s\n", tag, sm[j, 1], sm[j, 2],
                   strofreal(m[j], "%12.8g"), strofreal(want, "%12.8g"))
        }
    }
    printf("{txt}    %s: map vs the unit-variance fit, worst relative gap %s, c_F = %s\n", tag, strofreal(worst, "%9.2e"), strofreal(c, "%9.5f"))
    return(nbad)
}

end

capture program drop tr_run
program tr_run, rclass
    gettoken fn 0 : 0
    quietly log using "`fn'", text replace name(trlog)
    capture noisily `0'
    local rc = _rc
    quietly log close trlog
    return scalar rc = `rc'
end

capture program drop tr_hit
program tr_hit, rclass
    args fn phrase
    mata: st_local("hit", strofreal(tr_has(st_local("fn"), st_local("phrase"))))
    if "`phrase'" == "Std. err." {
        mata: st_local("hit2", strofreal(tr_has(st_local("fn"), "Std. Err.")))
        local hit = `hit' | `hit2'
    }
    return scalar hit = `hit'
end

capture program drop tr_expect
program tr_expect
    args fn kind phrase tag
    tr_hit "`fn'" "`phrase'"
    local hit = r(hit)
    if ("`kind'" == "has" & `hit' == 0) | ("`kind'" == "hasnt" & `hit' == 1) {
        display as error "test_condx_recovery: `tag': the output " cond("`kind'" == "has", "lacks", "contains") ": `phrase'"
        global TR_NFAIL = ${TR_NFAIL} + 1
    }
end

capture program drop tr_cond
program tr_cond
    args ok tag
    if `ok' != 1 {
        display as error "test_condx_recovery: `tag'"
        global TR_NFAIL = ${TR_NFAIL} + 1
    }
end

* the trap data (test_condx_marker.do's tm_data)
capture program drop tr_data
program tr_data
    args n l1 l2 l3 l4 seed
    clear
    set seed `seed'
    quietly set obs `n'
    quietly generate byte x = runiform() < .5
    quietly generate double eta = .5 * x + rnormal()
    local j 0
    foreach l in `l1' `l2' `l3' `l4' {
        local ++j
        local sd = sqrt(1 - (`l')^2)
        quietly generate double z`j' = `l' * eta + `sd' * rnormal() + (`j' == 4) * .3 * x
        quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
    }
    drop eta z1 z2 z3 z4
end

* test_condx_output.do's data (converged, admissible)
capture program drop tr_data2
program tr_data2
    clear
    set seed 20261114
    quietly set obs 1500
    quietly generate byte x = runiform() < .5
    quietly generate double eta = .4 * x + rnormal()
    forvalues j = 1/4 {
        quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
        quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
    }
    drop eta z1 z2 z3 z4
end

capture program drop tr_sethand
program tr_sethand, eclass
    args conv rec
    ereturn scalar converged = `conv'
    ereturn local recovery "`rec'"
end

tempfile LG LG2
local mB "(F -> y1 y2 y3 y4) (F <- x) (y4 <- x)"
local title "DWLS estimation, WLSMV test"
local seeds 20262415 20262516 20262617 20262718
local nexC3 0
local nexC4 0

* ---- C1, C2, C5: the four trap seeds ----
foreach seed of local seeds {
    display as text _n "  C1, C2, C5: seed `seed'"
    tr_data 500 .70 .70 -.70 -.70 `seed'

    tr_run "`LG'" dwlsmv `mB'
    local rc = r(rc)
    tr_cond `=`rc' == 0' "C1 seed `seed': the marker fit gave rc `rc', wanted 0"
    if `rc' continue
    local rec "`e(recovery)'"
    if "`rec'" != "recovered" {
        display as error "test_condx_recovery: C1 seed `seed': e(recovery) is `rec', wanted recovered" ///
            cond("`rec'" == "failed", " (the unit-variance first fit without the retry, or the restart, did not converge on this seed: the premise that it does is measured here, not assumed)", "")
        global TR_NFAIL = ${TR_NFAIL} + 1
        continue
    }
    tr_cond `=e(converged) == 1' "C1 seed `seed': recovered but e(converged) = `=e(converged)'"
    mata: st_local("nv", strofreal(rows(st_matrix("e(V)"))))
    tr_cond `=`nv' == 15' "C1 seed `seed': e(V) has `nv' rows, wanted 15"
    tr_cond `=("`e(T3)'" != "" & "`e(p_T3)'" != "" & e(df_T3) == 4)' "C1 seed `seed': e(T3), e(p_T3) or e(df_T3) not as for a converged fit"
    tr_cond `=("`e(covariates)'" == "x")' "C1 seed `seed': e(covariates) is `e(covariates)'"
    display as text "    iterations of the reported (restarted) fit: `=e(iterations)' (a diagnostic, not asserted)"
    tr_expect "`LG'" has "The marker-metric fit did not converge" "C1 seed `seed'"
    tr_expect "`LG'" has "it was restarted from the unit-variance solution, transformed to the" "C1 seed `seed'"
    tr_expect "`LG'" hasnt "second identification start" "C1 seed `seed'"
    tr_expect "`LG'" hasnt "Recovery attempt failed" "C1 seed `seed'"
    tr_expect "`LG'" hasnt "convergence not achieved" "C1 seed `seed'"
    tr_expect "`LG'" hasnt "last iterate" "C1 seed `seed'"
    tr_expect "`LG'" hasnt "No recovery refit" "C1 seed `seed'"
    tr_expect "`LG'" hasnt "did not converge from the ULS starting values" "C1 seed `seed': the first conditional marker fit takes no retry"
    tr_expect "`LG'" has "Std. err." "C1 seed `seed'"
    tr_expect "`LG'" has "Prob > chi2" "C1 seed `seed'"
    tr_expect "`LG'" has "Not computed for a model with a covariate in this version: the check for thresholds" "C1 seed `seed'"

    * the reported fit is the restarted one: its estimates and T3 are kept for C2
    tempname bM
    matrix `bM' = e(b)
    local t3M = e(T3)
    local dfM = e(df_T3)

    * replay: the same
    tr_run "`LG2'" dwlsmv
    tr_cond `=(r(rc) == 0 & "`e(recovery)'" == "recovered")' "C1 seed `seed': replay rc `=r(rc)' or recovery `e(recovery)'"
    tr_expect "`LG2'" hasnt "convergence not achieved" "C1 replay, seed `seed'"
    tr_expect "`LG2'" has "Std. err." "C1 replay, seed `seed'"

    * C5 first: the same data under metric(unitvariance); then the map against it (C2)
    tr_run "`LG'" dwlsmv `mB', metric(unitvariance)
    local rcu = r(rc)
    tr_cond `=(`rcu' == 0)' "C5 seed `seed': the unit-variance fit gave rc `rcu'"
    if `rcu' continue
    tr_cond `=(e(converged) == 1)' "C5 seed `seed': the unit-variance fit did not converge"
    if e(converged) != 1 continue
    tr_cond `=("`e(recovery)'" == "none")' "C5 seed `seed': metric(unitvariance) posts e(recovery) `e(recovery)', wanted none"
    tr_expect "`LG'" hasnt "marker-metric fit" "C5 seed `seed'"
    tempname bU
    matrix `bU' = e(b)
    local t3U = e(T3)
    local dfU = e(df_T3)
    mata: st_local("nbad", strofreal(tr_map("`bM'", "`bU'", 1e-5, "C2 seed `seed'")))
    tr_cond `=(`nbad' == 0)' "C2 seed `seed': `nbad' column(s) of the recovered e(b) are off the map by more than 1e-5 relative"
    tr_cond `=(abs(`t3M' - `t3U') <= 1e-6 * max(1, abs(`t3U')) & `dfM' == `dfU')' "C2 seed `seed': T3 `t3M' vs `t3U' (relative 1e-6), df `dfM' vs `dfU'"
}

* ---- C3: a fixed NONZERO gamma: ineligible when the first fit triggers ----
foreach seed of local seeds {
    display as text _n "  C3: (F <- x@.3), seed `seed'"
    tr_data 500 .70 .70 -.70 -.70 `seed'
    tr_run "`LG'" dwlsmv (F -> y1 y2 y3 y4) (F <- x@.3) (y4 <- x)
    local rc = r(rc)
    tr_cond `=`rc' == 0' "C3 seed `seed': rc `rc', wanted 0"
    if `rc' continue
    if "`e(recovery)'" == "none" {
        display as text "    first fit not triggered on this seed (converged and admissible): not exercised"
        continue
    }
    local ++nexC3
    tr_cond `=("`e(recovery)'" == "ineligible")' "C3 seed `seed': e(recovery) is `e(recovery)', wanted ineligible"
    tr_expect "`LG'" has "the model is not eligible for it:" "C3 seed `seed'"
    tr_expect "`LG'" has "gamma[F,x] is fixed at" "C3 seed `seed'"
    tr_expect "`LG'" hasnt "restarted from the unit-variance solution" "C3 seed `seed'"
    tr_expect "`LG'" hasnt "No recovery refit is attempted" "C3 seed `seed'"
    tr_expect "`LG'" hasnt "No recovery refit was attempted for a model with a covariate" "C3 seed `seed'"
    if e(converged) == 0 {
        tr_expect "`LG'" has "No recovery refit was attempted:" "C3 seed `seed' (header)"
        tr_expect "`LG'" has "the model is not eligible for it (see the notes)." "C3 seed `seed' (header)"
        tr_expect "`LG'" hasnt "Std. err." "C3 seed `seed'"
    }
}

* ---- C4: a fixed gamma 0 and a fixed kappa stay eligible ----
foreach spec in "(F <- x@0) (y4 <- x)" "(F <- x) (y4 <- x@.3)" {
    foreach seed of local seeds {
        display as text _n "  C4: `spec', seed `seed'"
        tr_data 500 .70 .70 -.70 -.70 `seed'
        tr_run "`LG'" dwlsmv (F -> y1 y2 y3 y4) `spec'
        local rc = r(rc)
        tr_cond `=`rc' == 0' "C4 `spec' seed `seed': rc `rc', wanted 0"
        if `rc' continue
        if "`e(recovery)'" == "none" {
            display as text "    first fit not triggered: not exercised"
            continue
        }
        local ++nexC4
        tr_cond `=("`e(recovery)'" != "ineligible")' "C4 `spec' seed `seed': ineligible, but nothing fixed here is a nonzero gamma"
        tr_expect "`LG'" hasnt "the model is not eligible for it" "C4 `spec' seed `seed'"
        if "`e(recovery)'" == "recovered" {
            tr_cond `=(e(converged) == 1)' "C4 `spec' seed `seed': recovered but not converged"
            tr_expect "`LG'" hasnt "second identification start" "C4 `spec' seed `seed'"
            if strpos("`spec'", "x@.3") {
                tr_cond `=(abs(_b[y4:x] - .3) < 1e-14)' "C4 `spec' seed `seed': the fixed kappa moved: _b[y4:x] = `=_b[y4:x]'"
            }
            else {
                // a path fixed at 0 is not posted in e(b) (as sem); absent or exactly 0 both pass (main session, 2026-11-15)
                capture local gF = _b[F:x]
                local absent = (_rc != 0)
                if !`absent' {
                    tr_cond `=(`gF' == 0)' "C4 `spec' seed `seed': the fixed gamma 0 moved: _b[F:x] = `gF'"
                }
            }
        }
        display as text "    e(recovery) = `e(recovery)'"
    }
}

* ---- C6: a converged, admissible conditional marker fit: none, as before ----
display as text _n "  C6: converged, admissible"
tr_data2
tr_run "`LG'" dwlsmv (F -> y1 y2 y3 y4) (F <- x)
tr_cond `=(r(rc) == 0 & e(converged) == 1 & "`e(recovery)'" == "none")' "C6: rc `=r(rc)', converged `=e(converged)', recovery `e(recovery)'"
tr_expect "`LG'" hasnt "marker-metric fit" "C6"
tr_expect "`LG'" hasnt "restarted" "C6"
tr_expect "`LG'" hasnt "No recovery refit" "C6"
tr_expect "`LG'" has "Std. err." "C6"

* ---- C7: an ordinary fit on the trap data: recovery unchanged ----
display as text _n "  C7: ordinary fit (no covariate), trap data"
foreach seed of local seeds {
    tr_data 500 .70 .70 -.70 -.70 `seed'
    tr_run "`LG'" dwlsmv (F -> y1 y2 y3 y4)
    tr_cond `=(r(rc) == 0)' "C7 seed `seed': rc `=r(rc)'"
    if r(rc) continue
    tr_cond `=inlist("`e(recovery)'", "none", "recovered")' "C7 seed `seed': e(recovery) is `e(recovery)'"
    tr_expect "`LG'" hasnt "Not computed for a model with a covariate" "C7 seed `seed'"
    tr_expect "`LG'" hasnt "No recovery refit" "C7 seed `seed'"
    if "`e(recovery)'" == "recovered" {
        tempname oM oU
        matrix `oM' = e(b)
        local t3M = e(T3)
        quietly dwlsmv (F -> y1 y2 y3 y4), metric(unitvariance)
        matrix `oU' = e(b)
        * the one-factor map without gamma and kappa: the same rule table (no x columns occur)
        mata: st_local("nbad", strofreal(tr_map("`oM'", "`oU'", 1e-5, "C7 seed `seed'")))
        tr_cond `=(`nbad' == 0)' "C7 seed `seed': the ordinary recovered fit is off the map in `nbad' column(s)"
        tr_cond `=(abs(`t3M' - e(T3)) <= 1e-6 * max(1, abs(e(T3))))' "C7 seed `seed': T3 differs from the unit-variance fit"
    }
    display as text "    seed `seed': e(recovery) = `e(recovery)'"
}

* ---- C8: the header of a nonconverged conditional fit, by recovery outcome (hand-set; spec D8 has no way to force it) ----
display as text _n "  C8: header by e(recovery), hand-set on a posted conditional fit"
tr_data2
quietly dwlsmv (F -> y1 y2 y3 y4) (F <- x)
foreach rec in ineligible failed none {
    tr_sethand 0 `rec'
    tr_run "`LG'" dwlsmv
    tr_cond `=(r(rc) == 0)' "C8 `rec': replay gave rc `=r(rc)' (the hand-set e() may need adjusting, not necessarily the code)"
    tr_expect "`LG'" has "convergence not achieved: the optimizer stopped after" "C8 `rec'"
    tr_expect "`LG'" has "Refit with metric(unitvariance)." "C8 `rec'"
    tr_expect "`LG'" hasnt "No recovery refit is attempted" "C8 `rec'"
    if "`rec'" == "ineligible" {
        tr_expect "`LG'" has "No recovery refit was attempted:" "C8 `rec'"
        tr_expect "`LG'" has "the model is not eligible for it (see the notes)." "C8 `rec'"
        tr_expect "`LG'" hasnt "and failed (see the notes)" "C8 `rec'"
    }
    else if "`rec'" == "failed" {
        tr_expect "`LG'" has "A recovery refit was attempted" "C8 `rec'"
        tr_expect "`LG'" has "and failed (see the notes)." "C8 `rec'"
        tr_expect "`LG'" hasnt "not eligible for it" "C8 `rec'"
    }
    else {
        tr_expect "`LG'" hasnt "recovery refit was attempted" "C8 `rec'"
        tr_expect "`LG'" hasnt "not eligible for it" "C8 `rec'"
    }
}

adopath - "`here'/src/ado"
adopath - "`here'/build"

if ${TR_NFAIL} {
    display as error _n "test_condx_recovery: ${TR_NFAIL} failure(s)"
    global TR_NFAIL
    exit 9
}
global TR_NFAIL
local ex3 = cond(`nexC3' > 0, "a fixed nonzero gamma was ineligible on `nexC3' triggered fit(s)", "C3 was NOT exercised (no seed's first fit triggered with gamma fixed at .3); the Mata-level eligibility check R3 stands alone")
local ex4 = cond(`nexC4' > 0, "fixed gamma 0 and fixed kappa were never ineligible on `nexC4' triggered fit(s)", "C4 was NOT exercised at the command level (R3 stands alone)")
display as result "test_condx_recovery: the map reproduces sigma(theta) on 200 random points and restarts exactly; the four trap seeds are recovered and agree with a separate unit-variance fit through the map; `ex3'; `ex4'; an ordinary fit is unchanged"
exit 0
