*! test_s2_guard.do -- Stage 2 refuses a conditional index map (spec 5.1).
*!
*! Executed 2026-11-12. Drafted by a read-only agent from tests/unit/test_o14.do
*! (loader, fixtures) and test_condx_stage1.do (log read-back); run by the main
*! session, which changed two expected messages to match the guard as built.
*! The tolerances 1e-14 and 1e-12 in the success checks were set by the agent
*! without measuring; they are closed-form identities, not tier comparisons.
*!
*! The guard is dwlsmv_s2_guard() / dwlsmv_s2_guard_st() (src/mata/
*! dwlsmv_s2_guard.mata), called first by dwls_identify, dwls_model, dwls_delta,
*! dwlsmv_fit_stack (so dwls_fit), dwls_inference, dwlsmv_sens_hess,
*! dwlsmv_admissible, dwlsmv_mg_admissible, dwlsmv_fitidx_base.
*! 2026-11-14: dwlsmv_admissible() now accepts one conditional map (test_condx_admissible.do); its
*! refusal case was removed from this file.
*!
*! dwls_identify() alone accepts a conditional map when the model text declares the covariate
*! (qok, 3rd argument of the guard); that acceptance and the two disagreements are checked in
*! tests/unit/test_condx_identify.do. Here every call carries a text with NO covariate.
*!
*! Checks:
*!   i    each entry point, handed a conditional map (q = 1) or Stage 1 result,
*!        aborts with rc == 3498 AND its own message on the log ("<who>: the
*!        index map of group <g> is conditional"). The message is checked
*!        because several of these calls would abort at 3498 anyway without the
*!        guard (a sigma never written, a group count that does not match), so
*!        the return code alone cannot tell the guard from the old refusals.
*!   iii  G = 2: dwls_identify with only the SECOND map conditional; and with
*!        both maps unconditional but only the second Stage 1 result's map
*!        conditional (built so that every older check passes: same K, s sized
*!        by the unconditional D, so without the guard the call succeeds). The
*!        same with dwls_model and dwls_delta (a G = 2 maps vector against a
*!        one-group table: the guard fires before the count check).
*!   ii   the same calls with the unconditional map of the same K succeed
*!        (rc 0) and return what they returned before the guard: table size,
*!        sigma's threshold block and one closed-form correlation, Delta's
*!        threshold block (identity, exactly), a converged fit and inference
*!        (D - npar df), the admissibility table's shape, the c1 baseline
*!        against its closed form. The G = 2 unconditional dwls_identify gives
*!        npar 7 (the value test_identify_o7_policy.do pins for this text).
*!        dwlsmv_sens_hess and dwlsmv_mg_admissible have no success case here:
*!        a valid call needs a G = 2 table or a fitted theta; tests/unit/
*!        test_sens_hess.do and test_mg_admissible.do cover them unchanged.
*!   iv   the guard itself: a q = 2 map is refused; an unfilled struct slot and
*!        a zero-length vector are not (left to the callers' own checks).
*!
*! Mutants for the main session (edit, run, revert); each must turn the file
*! red: delete the guard call from any one of the nine entry points (that
*! entry's case i); change "q > 0" to "q > 1" in dwlsmv_s2_guard (every q = 1
*! case); start the loop at g = 2 (group-1 cases); drop the args() test around
*! the st guard in dwls_identify (case 3 of iii: the call would then succeed);
*! make the _error literal 101 or more characters (rc 3300, every case).
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

* Load every source, types first, as test_o14.do does; the new file is picked
* up by the directory listing, so it needs no line here.
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
        display as error "test_s2_guard: `f' failed to load"
        exit 9
    }
}

capture mata: dwlsmv_s2_guard(dwlsmv_statmap(0), "load check")
if _rc {
    display as error "test_s2_guard: dwlsmv_s2_guard() is not defined (src/mata/dwlsmv_s2_guard.mata missing?)"
    exit 9
}

mata:

// ---- fixtures (after test_o14.do) -------------------------------------------

string scalar gt_spec()
{
    return("factor F" + char(10) + "ind y1 y2 y3 y4" + char(10) +
           "lambda F y1 free" + char(10) + "lambda F y2 free" + char(10) +
           "lambda F y3 free" + char(10) + "lambda F y4 free" + char(10) +
           "metric unitvariance F")
}

struct dwlsmv_statmap scalar gt_map()
{
    return(dwlsmv_statmap_build((3 \ 3 \ 2 \ 2)))
}

struct dwlsmv_statmap scalar gt_cmap()
{
    return(dwlsmv_statmap_build((3 \ 3 \ 2 \ 2), 1))
}

// theta for the table of gt_spec(): the T thresholds, then lambda y1..y4.
real colvector gt_theta(struct dwlsmv_parmap scalar par)
{
    real colvector th

    th = J(par.npar, 1, 0)
    th[1::6]                   = (-0.5 \ 0.5 \ -0.4 \ 0.4 \ -0.3 \ 0.3)
    th[(par.npar - 3)::par.npar] = (0.6 \ 0.7 \ 0.5 \ 0.8)
    return(th)
}

// An unconditional Stage 1 result that dwls_fit() takes to convergence
// (test_o14.do's to_st(): the thresholds at 0 and a small fixed perturbation).
struct dwlsmv_stage1 scalar gt_st()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    real colvector th

    st.map = gt_map()
    par    = dwls_identify(gt_spec(), st.map)
    th     = J(par.npar, 1, 0)
    th[(par.npar - 3)::par.npar] = (0.6 \ 0.7 \ 0.5 \ 0.8)
    st.s      = dwls_model(th, par, st.map) + 0.01 * (mod(1::st.map.D, 3) :- 1)
    st.Gamma  = I(st.map.D)
    st.N      = 500
    st.psibar = J(st.map.D, 1, 0)
    st.bound  = J(st.map.P, 1, 0)
    return(st)
}

// A conditional Stage 1 result, filled consistently with its own map
// (D = 16 at q = 1): nothing but the guard objects to it.
struct dwlsmv_stage1 scalar gt_cst()
{
    struct dwlsmv_stage1 scalar st

    st.map    = gt_cmap()
    st.s      = J(st.map.D, 1, 0.3)
    st.Gamma  = I(st.map.D)
    st.N      = 500
    st.psibar = J(st.map.D, 1, 0)
    st.bound  = J(st.map.P, 1, 0)
    return(st)
}

// The one-group table and theta for gt_spec() on the unconditional map.
void gt_table(struct dwlsmv_parmap scalar par, real colvector th)
{
    par = dwls_identify(gt_spec(), gt_map())
    th  = gt_theta(par)
}

// The G = 2 canonical text of test_identify_o7_policy.do (top_id), K = (3, 3)
// in both groups. Group 2's tau rows are tied to group 1's.
string scalar gt_gspec(string scalar al)
{
    string rowvector a
    string scalar s
    real scalar i

    a = ("groups 2", "group 1", "factor F", "ind y1 y2",
         "lambda F y1 fixed 0.7", "lambda F y2 fixed 0.7", "metric unitvariance F",
         "group 2", "factor F", "ind y1 y2", "lambda F y1 free", "lambda F y2 free",
         "psi F F fixed 1", "alpha F " + al,
         "tie tau 1 y1 1 tau 2 y1 1", "tie tau 1 y1 2 tau 2 y1 2",
         "tie tau 1 y2 1 tau 2 y2 1", "tie tau 1 y2 2 tau 2 y2 2")
    s = ""
    for (i = 1; i <= cols(a); i++) {
        s = s + a[i] + char(10)
    }
    return(s)
}

// Two groups' maps (K = (3, 3)) and Stage 1 results: group g's map is
// conditional (q = 1) where cm[g] is 1, its Stage 1 result's map where cs[g]
// is 1. Each s is sized by its own map's D (5, or 7 when conditional), with the
// first four entries (the thresholds) set.
void gt_two(real rowvector cm, real rowvector cs,
            struct dwlsmv_statmap vector m, struct dwlsmv_stage1 vector st)
{
    real scalar g

    m  = dwlsmv_statmap(1, 2)
    st = dwlsmv_stage1(1, 2)
    for (g = 1; g <= 2; g++) {
        m[g]      = dwlsmv_statmap_build(J(2, 1, 3), cm[g])
        st[g].map = dwlsmv_statmap_build(J(2, 1, 3), cs[g])
        st[g].s   = J(st[g].map.D, 1, 0.3)
        st[g].s[1..4] = (-0.4 \ 0.4 \ -0.3 \ 0.5)
    }
}

// ---- (i), (iii): each call must abort; the caller reads the code and the log ---

real scalar gt_c_id1()
{
    (void) dwls_identify(gt_spec(), gt_cmap())
    return(0)
}

real scalar gt_c_id2()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st

    gt_two((0, 1), (0, 1), m, st)
    (void) dwls_identify(gt_gspec("free"), m, st)
    return(0)
}

// maps unconditional, only the second Stage 1 result's map conditional, s sized
// by the unconditional D: every older check passes, so only the guard refuses.
real scalar gt_c_id3()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st

    gt_two((0, 0), (0, 1), m, st)
    st[2].s = J(5, 1, 0.3)
    st[2].s[1..4] = (-0.4 \ 0.4 \ -0.3 \ 0.5)
    (void) dwls_identify(gt_gspec("free"), m, st)
    return(0)
}

real scalar gt_c_model()
{
    struct dwlsmv_parmap scalar par
    real colvector th

    gt_table(par, th)
    (void) dwls_model(th, par, gt_cmap())
    return(0)
}

// G = 2 maps, second conditional, against a one-group table: the count check
// would refuse this at 3498 too; the guard must speak first.
real scalar gt_c_model2()
{
    struct dwlsmv_statmap vector mv
    struct dwlsmv_parmap scalar par
    real colvector th

    gt_table(par, th)
    mv = dwlsmv_statmap(2, 1)
    mv[1] = gt_map()
    mv[2] = gt_cmap()
    (void) dwls_model(th, par, mv)
    return(0)
}

real scalar gt_c_delta()
{
    struct dwlsmv_parmap scalar par
    real colvector th

    gt_table(par, th)
    (void) dwls_delta(th, par, gt_cmap())
    return(0)
}

real scalar gt_c_delta2()
{
    struct dwlsmv_statmap vector mv
    struct dwlsmv_parmap scalar par
    real colvector th

    gt_table(par, th)
    mv = dwlsmv_statmap(2, 1)
    mv[1] = gt_map()
    mv[2] = gt_cmap()
    (void) dwls_delta(th, par, mv)
    return(0)
}

// dwls_fit()'s own stack call (flag 1) with the SECOND group conditional: dwls_fit() itself would
// reach dwlsmv_mg_check() first on a one-group table, so the stack function is called directly.
real scalar gt_c_fit()
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_statmap vector maps
    real colvector s, w
    real scalar N

    sv = dwlsmv_stage1(2, 1)
    sv[1] = gt_st()
    sv[2] = gt_cst()
    dwlsmv_fit_stack(sv, s, w, N, maps, 1)
    return(0)
}

// The diagnostics' call: no flag, so one conditional group is still refused.
real scalar gt_c_fit0()
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_statmap vector maps
    real colvector s, w
    real scalar N

    sv = dwlsmv_stage1(1, 1)
    sv[1] = gt_cst()
    dwlsmv_fit_stack(sv, s, w, N, maps)
    return(0)
}

// A fit assembled by hand, as test_o14.do's to_inf2() does, carrying a
// conditional Stage 1 result.
real scalar gt_c_inf()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f

    st  = gt_st()
    par = dwls_identify(gt_spec(), st.map)
    f   = dwls_fit(st, par)
    sv  = dwlsmv_stage1(2, 1)
    sv[1] = gt_st()
    sv[2] = gt_cst()
    f.stage1 = sv
    dwls_inference(f)
    return(0)
}

real scalar gt_c_hess()
{
    struct dwlsmv_statmap vector mv
    struct dwlsmv_parmap scalar par
    real colvector th

    gt_table(par, th)
    mv = dwlsmv_statmap(1, 1)
    mv[1] = gt_cmap()
    (void) dwlsmv_sens_hess(th, par, mv)
    return(0)
}

real scalar gt_c_mgadm()
{
    struct dwlsmv_parmap scalar par
    real colvector th

    gt_table(par, th)
    (void) dwlsmv_mg_admissible(th, par, 1, gt_cmap())
    return(0)
}

real scalar gt_c_fidx()
{
    struct dwlsmv_stage1 vector sv

    sv = dwlsmv_stage1(1, 1)
    sv[1] = gt_cst()
    (void) dwlsmv_fitidx_base(sv)
    return(0)
}

// the guard itself, q = 2
real scalar gt_c_q2()
{
    dwlsmv_s2_guard(dwlsmv_statmap_build((3 \ 3 \ 2 \ 2), 2), "direct")
    return(0)
}

// ---- (ii), (iv): each returns its number of failures --------------------------

real scalar gt_fail(string scalar what)
{
    printf("{err}test_s2_guard: %s\n", what)
    return(1)
}

real scalar gt_ok_id()
{
    struct dwlsmv_parmap scalar p
    struct dwlsmv_statmap scalar m
    real scalar nf

    nf = 0
    m  = gt_map()
    p  = dwls_identify(gt_spec(), m)
    // T free thresholds and four free lambdas; unit variance fixes psi
    if (p.npar != m.T + 4) {
        nf = nf + gt_fail("unconditional dwls_identify: npar is " + strofreal(p.npar) +
                          ", wanted T + 4 = " + strofreal(m.T + 4))
    }
    return(nf)
}

real scalar gt_ok_id2()
{
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st
    struct dwlsmv_parmap scalar p
    real scalar nf

    nf = 0
    gt_two((0, 0), (0, 0), m, st)
    p = dwls_identify(gt_gspec("free 0.5"), m, st)
    if (p.npar != 7) {
        nf = nf + gt_fail("unconditional G = 2 dwls_identify: npar is " + strofreal(p.npar) + ", wanted 7")
    }
    return(nf)
}

real scalar gt_ok_model_delta()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap scalar m
    real colvector th, sg
    real matrix Dm
    real scalar nf, pos, want

    nf = 0
    m  = gt_map()
    gt_table(par, th)
    sg = dwls_model(th, par, m)
    if (rows(sg) != m.D) {
        nf = nf + gt_fail("dwls_model: sigma has " + strofreal(rows(sg)) + " rows, wanted D")
    }
    else {
        if (sg[1::m.T] != th[1::m.T]) {
            nf = nf + gt_fail("dwls_model: sigma's threshold block is not theta's")
        }
        // one factor at unit variance: R_12 = lambda_1 lambda_2 (th[T+1], th[T+2])
        pos  = dwlsmv_pos_rho(m, 1, 2)
        want = th[m.T + 1] * th[m.T + 2]
        if (reldif(sg[pos], want) > 1e-14) {
            nf = nf + gt_fail("dwls_model: R_12 is not lambda_1 lambda_2 to 1e-14")
        }
    }
    Dm = dwls_delta(th, par, m)
    if (rows(Dm) != m.D | cols(Dm) != par.npar) {
        nf = nf + gt_fail("dwls_delta: Delta is " + strofreal(rows(Dm)) + " x " + strofreal(cols(Dm)))
    }
    else if (Dm[1::m.T, 1::m.T] != I(m.T)) {
        // single group, no invariance constraint: the threshold block is the identity
        nf = nf + gt_fail("dwls_delta: the threshold block of Delta is not the identity")
    }
    return(nf)
}

real scalar gt_ok_fit_inf_adm()
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar f
    real matrix A
    real scalar nf

    nf  = 0
    st  = gt_st()
    par = dwls_identify(gt_spec(), st.map)
    f   = dwls_fit(st, par)
    if (!f.converged) {
        nf = nf + gt_fail("dwls_fit: the unconditional fit did not converge")
    }
    if (rows(f.theta) != par.npar | rows(f.Delta) != st.map.D) {
        nf = nf + gt_fail("dwls_fit: theta or Delta has the wrong size")
    }
    dwls_inference(f)
    if (f.T3_df != st.map.D - par.npar) {
        nf = nf + gt_fail("dwls_inference: T3_df is not D - npar")
    }
    A = dwlsmv_admissible(f.theta, par, st.map)
    if (rows(A) != 6 | cols(A) != 4) {
        nf = nf + gt_fail("dwlsmv_admissible: the table is " + strofreal(rows(A)) + " x " + strofreal(cols(A)))
    }
    return(nf)
}

// The c1 baseline in closed form on this fixture: Gamma = I, so every
// polychoric variance is 1, R_g = I, a_B = 1, and T3_B = T_B = (N - 1) sum rho^2.
real scalar gt_ok_fidx()
{
    struct dwlsmv_stage1 vector sv
    real rowvector base
    real colvector r
    real scalar nf, want

    nf = 0
    sv = dwlsmv_stage1(1, 1)
    sv[1] = gt_st()
    base  = dwlsmv_fitidx_base(sv)
    r     = sv[1].s[(sv[1].map.T + 1)::sv[1].map.D]
    want  = (sv[1].N - 1) * quadsum(r :^ 2)
    if (base[2] != sv[1].map.P) {
        nf = nf + gt_fail("dwlsmv_fitidx_base: D_B is not P")
    }
    if (reldif(base[3], want) > 1e-12 | reldif(base[1], want) > 1e-12) {
        nf = nf + gt_fail("dwlsmv_fitidx_base: T_B or T3_B is not (N - 1) sum rho^2 to 1e-12")
    }
    return(nf)
}

// An unfilled struct slot and a zero-length vector are not the guard's business.
real scalar gt_ok_unfilled()
{
    struct dwlsmv_statmap vector mv
    struct dwlsmv_stage1 vector sv

    mv = dwlsmv_statmap(2, 1)
    mv[1] = gt_map()
    dwlsmv_s2_guard(mv, "unfilled")
    sv = dwlsmv_stage1(1, 1)
    dwlsmv_s2_guard_st(sv, "unfilled")
    dwlsmv_s2_guard(dwlsmv_statmap(0), "empty")
    dwlsmv_s2_guard_st(dwlsmv_stage1(0), "empty")
    return(0)
}

// Does the log at f contain the text has? Returns the number of failures.
real scalar gt_msg(string scalar f, string scalar has, string scalar tag)
{
    string scalar txt

    txt = invtokens(cat(f)', char(10))
    if (!strpos(txt, has)) {
        printf("{err}test_s2_guard: case %s: the log lacks: %s\n", tag, has)
        return(1)
    }
    return(0)
}

end

local nfail 0
tempfile L

* ---- (i), (iii): rc 3498 and the guard's own message ---------------------------
local n 0
local ++n
local c`n' gt_c_id1()
local w`n' dwls_identify: the index map of group 1 is conditional
local ++n
local c`n' gt_c_id2()
local w`n' dwls_identify: the index map of group 2 is conditional
local ++n
local c`n' gt_c_id3()
local w`n' dwls_identify: the index map of group 2 is conditional
local ++n
local c`n' gt_c_model2()
local w`n' dwls_model: the index map of group 2 is conditional
local ++n
local c`n' gt_c_delta2()
local w`n' dwls_delta: the index map of group 2 is conditional
local ++n
local c`n' gt_c_fit()
local w`n' dwls_fit: the index map of group 2 is conditional
local ++n
local c`n' gt_c_fit0()
local w`n' dwls_fit: the index map of group 1 is conditional
local ++n
local c`n' gt_c_inf()
local w`n' dwls_inference: the index map of group 2 is conditional
local ++n
local c`n' gt_c_hess()
local w`n' dwlsmv_sens_hess: the index map of group 1 is conditional
local ++n
local c`n' gt_c_mgadm()
local w`n' dwlsmv_mg_admissible: the index map of group 1 is conditional
local ++n
local c`n' gt_c_fidx()
local w`n' dwlsmv_fitidx_base: the index map of group 1 is conditional
local ++n
local c`n' gt_c_q2()
local w`n' direct: the index map of group 1 is conditional (q = 2)

forvalues i = 1/`n' {
    quietly log using "`L'_`i'.log", text replace name(gt)
    capture noisily mata: (void) `c`i''
    local rc = _rc
    quietly log close gt
    display as text "  `c`i'': rc `rc'"
    if `rc' != 3498 {
        display as error "test_s2_guard: `c`i'' gave rc `rc', wanted 3498"
        local nfail = `nfail' + 1
    }
    else {
        local tmpL "`L'_`i'.log"
        local want "`w`i''"
        capture noisily mata: st_local("nm", strofreal(gt_msg(st_local("tmpL"), st_local("want"), "`i'")))
        if _rc {
            display as error "test_s2_guard: reading the log of `c`i'' failed, rc " _rc
            local nfail = `nfail' + 1
        }
        else local nfail = `nfail' + `nm'
    }
}

* ---- (ii), (iv): unconditional maps still work, results as before ---------------
* gt_c_model and gt_c_delta: since 2026-11-13 dwls_model() and dwls_delta() accept ONE group
* with a conditional map (slope rows, tests/unit/test_condx_model_delta.do); they must now
* complete. Two groups with a conditional map stay refused (gt_c_model2, gt_c_delta2 above).
foreach c in gt_c_model gt_c_delta gt_ok_id gt_ok_id2 gt_ok_model_delta gt_ok_fit_inf_adm gt_ok_fidx gt_ok_unfilled {
    capture noisily mata: st_local("nok", strofreal(`c'()))
    if _rc {
        display as error "test_s2_guard: `c'() aborted, rc " _rc ", wanted 0"
        local nfail = `nfail' + 1
    }
    else {
        display as text "  `c'(): `nok' failure(s)"
        local nfail = `nfail' + `nok'
    }
}

if `nfail' {
    display as error _n "test_s2_guard: `nfail' failure(s)"
    exit 9
}

display as result "test_s2_guard: every Stage 2 entry point refuses a conditional map or Stage 1 result at rc 3498 with its own message, G = 2 included, except dwls_model, dwls_delta, dwls_fit and dwls_inference at one group, which accept one covariate (tests/unit/test_condx_model_delta.do, test_condx_fit.do); the unconditional calls are unchanged"
exit 0
