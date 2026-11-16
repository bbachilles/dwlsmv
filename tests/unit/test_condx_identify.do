*! test_condx_identify.do -- dwls_identify() and the conditional-covariate directives.
*!
*! Executed by the main session (Stata 18 MP) and mutant-tested, 2026-11-14 at the latest. Drafted by a read-only agent 2026-11-13 (revision 2, to the owner's
*! decisions of 2026-11-13) against identify_changes.md (dwlsmv_id_cx.mata, the guard's optional
*! third argument, the gamma and kappa loops in dwls_identify, the optional line-number argument
*! of dwlsmv_id_parse). The main session loads, runs and fixes it; remove this banner only then.
*!
*! DEPENDS ON THE PARALLEL DRAFT: dwls_identify() ends in the rank check, which calls dwls_delta().
*! Checks a1 to a4 therefore pass only once dwls_delta() accepts the conditional map and ptypes 8
*! and 9 (the slope rows). The refusal checks r and the line-number checks L do not reach the rank
*! check and need only this draft.
*!
*! Fixture: one factor F, four indicators y1..y4 with K = (3, 3, 2, 2) (T = 6), one covariate x,
*! unconditional map D = 11 and conditional map (q = 1) D = 16. No sample statistic enters: the
*! table is a function of the text and the map's shape alone, and no number below is a fixture
*! value, only the table's own structure.
*!
*! Checks:
*!   a1  theta order and ptypes on a model with gamma and kappa: the table is tau (6), lambda (4),
*!       gamma (ptype 8), kappa (ptype 9, ordered by indicator), psi; names, rows, columns, free
*!       and fixed values, starts (gamma 0, a kappa start 0.1, a fixed kappa -0.3), fidx 1..12
*!       over the free rows with 0 on the fixed ones, grp all 1, nall 14, npar 12. A Stage 1
*!       vector whose map is conditional is accepted too.
*!   a2  unitvariance on a predicted factor: with gamma on F2 only, psi[F2,F2] is a fixed ptype-2
*!       row at 1 (the disturbance variance); no gamma row for F1; gamma's row and column.
*!   a3  no-covariate invariance: the unconditional model's table has the old shape (nall 11,
*!       ptypes, fidx), and the covariate model's tau, lambda and psi rows equal it in every
*!       column but fidx (compared within this test from the two calls; no stored numbers).
*!   a4  a kappa fixed at 0 is not free (owner, 2026-11-13): kappa on y1 to y3 free and y4
*!       `fixed 0` is accepted, y4 being the anchor; the row is carried (free 0, value 0).
*!   a5  screen N accepts (each model also passes the rank check; @sec-condx-anchor): (a) a factor
*!       with no gamma whose indicators all carry a free kappa; (b) a factor whose gamma is fixed at
*!       a nonzero value, all its kappas free; (c) a free-gamma F1 with all its own kappas free,
*!       anchored through `beta F2 F1 free` to F2, whose indicators carry no kappa; (d) a kappa
*!       FIXED at a nonzero value is no free kappa, so its indicator is the only anchor of a
*!       free-gamma factor.
*!   r   refusals, each rc == 3498 exactly: every indicator of a free-gamma F carries a free
*!       kappa and no beta leaves F; the same F1 with its only beta `fixed 0` (not present), and
*!       with the beta REVERSED (`beta F1 F2`: F2 predicts F1, so F2 is no descendant of F1); the
*!       only kappa-free indicator has its loading fixed at 0 (a loading fixed at 0 is not
*!       present); a conditional map and no
*!       covariate; a covariate and an unconditional map; a covariate and a conditional map of the
*!       wrong q; two covariates on one line; two 'covariate' lines; a covariate under 'groups';
*!       gamma without a covariate line; kappa on an undeclared indicator; gamma on an undeclared
*!       factor; gamma naming an undeclared covariate; a duplicate gamma; a covariate named like
*!       an indicator and like a factor (decision 4); a kappa line with too few tokens, a bad
*!       status word, and a non-numeric start; a conditional map but an unconditional Stage 1
*!       result.
*!   L   line numbers (decision 5), read back from the log: every message names the line as
*!       written, blank lines and comments counted. L1 an error after covariate lines (in
*!       dwlsmv_id_parse); L2 an error in a gamma line; L3 in a kappa line; L4 the second
*!       'covariate' directive; L5 a text with NO covariate and blank/comment lines above the error;
*!       L6 CRLF line ends. Each also asserts that the number the earlier numbering would have
*!       printed is absent.
*!
*! The rc is exact (3498) wherever a message is read back.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

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
        display as error "test_condx_identify: `f' failed to load"
        exit 9
    }
}

capture mata: (void) dwlsmv_id_cxsplit("", J(0, 1, ""), J(0, 1, .), J(0, 1, ""), J(0, 1, .), J(0, 1, ""), .)
if _rc {
    display as error "test_condx_identify: dwlsmv_id_cxsplit() is not defined (src/mata/dwlsmv_id_cx.mata missing?)"
    exit 9
}

mata:

// ---- fixtures -----------------------------------------------------------------

string scalar cx_join(string rowvector a)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(a); i++) {
        s = s + a[i] + char(10)
    }
    return(s)
}

// One factor, four indicators, unit variance; extra lines are appended.
string scalar cx_spec(string rowvector extra)
{
    string rowvector base

    base = ("factor F", "ind y1 y2 y3 y4",
            "lambda F y1 free", "lambda F y2 free", "lambda F y3 free", "lambda F y4 free",
            "metric unitvariance F")
    if (cols(extra) == 0) {
        return(cx_join(base))
    }
    return(cx_join((base, extra)))
}

// The covariate model of check a1.
string scalar cx_spec1()
{
    return(cx_spec(("covariate x", "gamma F x free", "kappa y3 x free 0.1",
                    "kappa y4 x fixed -0.3")))
}

struct dwlsmv_statmap scalar cx_umap()
{
    return(dwlsmv_statmap_build((3 \ 3 \ 2 \ 2)))
}

struct dwlsmv_statmap scalar cx_cmap(real scalar q)
{
    return(dwlsmv_statmap_build((3 \ 3 \ 2 \ 2), q))
}

// Two factors, two indicators each; F2 is predicted by x, F1 is not.
string scalar cx_spec2()
{
    return(cx_join(("factor F1 F2", "ind y1 y2 y3 y4", "covariate x",
                    "lambda F1 y1 free", "lambda F1 y2 free",
                    "lambda F2 y3 free", "lambda F2 y4 free",
                    "gamma F2 x free", "psi F1 F2 free",
                    "metric unitvariance F1", "metric unitvariance F2")))
}

// Six indicators, three categories each: the conditional map for the two-factor screen models.
struct dwlsmv_statmap scalar cx_cmap6(real scalar q)
{
    return(dwlsmv_statmap_build((3 \ 3 \ 3 \ 3 \ 3 \ 3), q))
}

// Two factors, three indicators each (F1: y1-y3, F2: y4-y6), unit variance on both; the
// covariate line and the rest are appended by the caller.
string scalar cx_spec6(string rowvector extra)
{
    string rowvector base

    base = ("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
            "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
            "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
            "metric unitvariance F1", "metric unitvariance F2")
    return(cx_join((base, extra)))
}

// Row of the table with the given name, or 0.
real scalar cx_prow(struct dwlsmv_parmap scalar par, string scalar nm)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.name[i] == nm) {
            return(i)
        }
    }
    return(0)
}

real scalar cx_chk(real scalar ok, string scalar msg)
{
    if (!ok) {
        printf("{err}test_condx_identify: %s\n", msg)
        return(1)
    }
    return(0)
}

// Equal, with the dimensions checked first (Mata's & evaluates both operands).
real scalar cx_eqr(real colvector a, real colvector b)
{
    if (rows(a) != rows(b)) {
        return(0)
    }
    return(all(a :== b))
}

real scalar cx_eqs(string colvector a, string colvector b)
{
    if (rows(a) != rows(b)) {
        return(0)
    }
    return(all(a :== b))
}

// ---- a1: theta order, ptypes, rows, values --------------------------------------

real scalar cx_a1()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_stage1 vector st
    real colvector want_pt, want_fidx
    real scalar nf

    nf  = 0
    par = dwls_identify(cx_spec1(), cx_cmap(1))

    nf = nf + cx_chk(par.nall == 14, "a1: nall is " + strofreal(par.nall) + ", wanted 14")
    nf = nf + cx_chk(par.npar == 12, "a1: npar is " + strofreal(par.npar) + ", wanted 12")

    want_pt = (J(6, 1, 4) \ J(4, 1, 1) \ 8 \ 9 \ 9 \ 2)
    nf = nf + cx_chk(cx_eqr(par.ptype, want_pt), "a1: the ptype sequence is not tau, lambda, gamma, kappa, kappa, psi")

    // Rows 11..14: gamma[F,x], kappa[y3,x], kappa[y4,x], psi[F,F].
    nf = nf + cx_chk(par.name[11] == "gamma[F,x]", "a1: row 11 is named " + par.name[11])
    nf = nf + cx_chk(par.name[12] == "kappa[y3,x]", "a1: row 12 is named " + par.name[12])
    nf = nf + cx_chk(par.name[13] == "kappa[y4,x]", "a1: row 13 is named " + par.name[13])
    nf = nf + cx_chk(par.name[14] == "psi[F,F]", "a1: row 14 is named " + par.name[14])
    nf = nf + cx_chk(par.name[7] == "lambda[y1,F]" & par.name[10] == "lambda[y4,F]",
                     "a1: the lambda rows moved")

    // row = factor / indicator, col = covariate.
    nf = nf + cx_chk(par.row[11] == 1 & par.col[11] == 1, "a1: gamma (row, col) is not (1, 1)")
    nf = nf + cx_chk(par.row[12] == 3 & par.col[12] == 1, "a1: kappa[y3,x] (row, col) is not (3, 1)")
    nf = nf + cx_chk(par.row[13] == 4 & par.col[13] == 1, "a1: kappa[y4,x] (row, col) is not (4, 1)")

    // Free and fixed, start and fixed values. Gamma starts at 0 (as beta does).
    nf = nf + cx_chk(par.free[11] == 1 & par.value[11] == 0, "a1: gamma is not free with start 0")
    nf = nf + cx_chk(par.free[12] == 1 & par.value[12] == 0.1, "a1: kappa[y3,x] is not free with start 0.1")
    nf = nf + cx_chk(par.free[13] == 0 & par.value[13] == -0.3, "a1: kappa[y4,x] is not fixed at -0.3")
    nf = nf + cx_chk(par.free[14] == 0 & par.value[14] == 1, "a1: psi[F,F] is not fixed at 1")

    // fidx numbers the free rows 1..12 in table order; fixed rows carry 0.
    want_fidx = (1::12 \ 0 \ 0)
    nf = nf + cx_chk(cx_eqr(par.fidx, want_fidx), "a1: fidx is not 1..12 over the free rows with 0 on the fixed ones")
    nf = nf + cx_chk(all(par.grp :== 1), "a1: grp is not all 1")

    // A Stage 1 result carrying the conditional map is accepted (it is not read).
    st = dwlsmv_stage1(1, 1)
    st[1].map = cx_cmap(1)
    par = dwls_identify(cx_spec1(), cx_cmap(1), st)
    nf = nf + cx_chk(par.nall == 14, "a1: the call with a conditional Stage 1 result gave nall " + strofreal(par.nall))
    return(nf)
}

// ---- a2: unitvariance fixes the disturbance variance of a predicted factor -------

real scalar cx_a2()
{
    struct dwlsmv_parmap scalar par
    real scalar nf, i, ipsi2, ngam

    nf = 0
    par = dwls_identify(cx_spec2(), cx_cmap(1))

    ipsi2 = 0
    ngam  = 0
    for (i = 1; i <= par.nall; i++) {
        if (par.name[i] == "psi[F2,F2]") {
            ipsi2 = i
        }
        if (par.ptype[i] == 8) {
            ngam = ngam + 1
            nf = nf + cx_chk(par.name[i] == "gamma[F2,x]" & par.row[i] == 2 & par.col[i] == 1,
                             "a2: the gamma row is " + par.name[i])
        }
    }
    nf = nf + cx_chk(ngam == 1, "a2: " + strofreal(ngam) + " gamma rows, wanted 1 (F1 is not predicted)")
    nf = nf + cx_chk(ipsi2 > 0, "a2: no psi[F2,F2] row")
    if (ipsi2 > 0) {
        nf = nf + cx_chk(par.ptype[ipsi2] == 2 & par.free[ipsi2] == 0 & par.value[ipsi2] == 1,
                         "a2: psi[F2,F2], F2 being predicted, is not a fixed ptype-2 row at 1")
    }
    return(nf)
}

// ---- a3: a model with no covariate is what it was ---------------------------------

real scalar cx_a3()
{
    struct dwlsmv_parmap scalar p0, p1
    real colvector keep
    real scalar nf

    nf = 0
    p0 = dwls_identify(cx_spec(J(1, 0, "")), cx_umap())

    nf = nf + cx_chk(p0.nall == 11 & p0.npar == 10, "a3: the unconditional model has nall " +
                     strofreal(p0.nall) + ", npar " + strofreal(p0.npar) + ", wanted 11 and 10")
    nf = nf + cx_chk(cx_eqr(p0.ptype, (J(6, 1, 4) \ J(4, 1, 1) \ 2)),
                     "a3: the unconditional ptype sequence changed")
    nf = nf + cx_chk(cx_eqr(p0.fidx, (1::10 \ 0)), "a3: the unconditional fidx changed")
    nf = nf + cx_chk(!any(p0.ptype :== 8) & !any(p0.ptype :== 9), "a3: the unconditional table has gamma or kappa rows")
    nf = nf + cx_chk(all(p0.grp :== 1), "a3: the unconditional grp changed")

    // The covariate model, less its gamma and kappa rows, is the same table.
    p1   = dwls_identify(cx_spec1(), cx_cmap(1))
    keep = selectindex(p1.ptype :!= 8 :& p1.ptype :!= 9)
    nf = nf + cx_chk(rows(keep) == p0.nall, "a3: the covariate model less gamma and kappa has " +
                     strofreal(rows(keep)) + " rows, the unconditional one " + strofreal(p0.nall))
    if (rows(keep) == p0.nall) {
        nf = nf + cx_chk(cx_eqs(p1.name[keep], p0.name), "a3: names differ")
        nf = nf + cx_chk(cx_eqr(p1.ptype[keep], p0.ptype), "a3: ptypes differ")
        nf = nf + cx_chk(cx_eqr(p1.row[keep], p0.row), "a3: rows differ")
        nf = nf + cx_chk(cx_eqr(p1.col[keep], p0.col), "a3: columns differ")
        nf = nf + cx_chk(cx_eqr(p1.free[keep], p0.free), "a3: free flags differ")
        nf = nf + cx_chk(cx_eqr(p1.value[keep], p0.value), "a3: values differ")
        nf = nf + cx_chk(cx_eqr(p1.grp[keep], p0.grp), "a3: grp differs")
    }
    return(nf)
}

// ---- a4: a kappa fixed at 0 is not present -----------------------------------------

real scalar cx_a4()
{
    struct dwlsmv_parmap scalar par
    real scalar nf

    nf  = 0
    par = dwls_identify(cx_spec(("covariate x", "gamma F x free", "kappa y1 x free",
                                 "kappa y2 x free", "kappa y3 x free", "kappa y4 x fixed 0")),
                        cx_cmap(1))

    // tau 6, lambda 4, gamma, kappa x 4, psi = 16 rows; free: 6 + 4 + 1 + 3 = 14.
    nf = nf + cx_chk(par.nall == 16, "a4: nall is " + strofreal(par.nall) + ", wanted 16")
    nf = nf + cx_chk(par.npar == 14, "a4: npar is " + strofreal(par.npar) + ", wanted 14")
    nf = nf + cx_chk(par.name[15] == "kappa[y4,x]", "a4: row 15 is named " + par.name[15])
    nf = nf + cx_chk(par.ptype[15] == 9 & par.free[15] == 0 & par.value[15] == 0,
                     "a4: kappa[y4,x] is not a carried fixed row at 0")
    return(nf)
}

// ---- a5: screen N accepts what the rank check accepts -----------------------------------
// Each call also runs the rank check of Delta inside dwls_identify(); an abort here is a failure.

real scalar cx_a5()
{
    struct dwlsmv_parmap scalar par
    real scalar nf, r

    nf = 0

    // (a) F2 has no gamma and a free kappa on every indicator: no requirement on F2.
    par = dwls_identify(cx_spec6(("covariate x", "gamma F1 x free", "psi F1 F2 free",
                                  "kappa y4 x free", "kappa y5 x free", "kappa y6 x free")),
                        cx_cmap6(1))
    nf = nf + cx_chk(sum(par.ptype :== 8) == 1 & sum(par.ptype :== 9) == 3,
                     "a5a: wanted one gamma row and three kappa rows")

    // (b) gamma[F2,x] fixed at 0.5 (not free), kappa free on all of y4-y6.
    par = dwls_identify(cx_spec6(("covariate x", "gamma F1 x free", "gamma F2 x fixed 0.5",
                                  "psi F1 F2 free",
                                  "kappa y4 x free", "kappa y5 x free", "kappa y6 x free")),
                        cx_cmap6(1))
    r = cx_prow(par, "gamma[F2,x]")
    nf = nf + cx_chk(r > 0, "a5b: no gamma[F2,x] row")
    if (r > 0) {
        nf = nf + cx_chk(par.free[r] == 0 & par.value[r] == 0.5, "a5b: gamma[F2,x] is not fixed at 0.5")
    }
    nf = nf + cx_chk(sum(par.ptype :== 9) == 3, "a5b: wanted three kappa rows")

    // (c) F1: free gamma, free kappa on y1-y3; anchored through the free beta to F2 (no kappa).
    par = dwls_identify(cx_spec6(("covariate x", "gamma F1 x free", "kappa y1 x free",
                                  "kappa y2 x free", "kappa y3 x free", "beta F2 F1 free")),
                        cx_cmap6(1))
    nf = nf + cx_chk(sum(par.ptype :== 3) == 1 & sum(par.ptype :== 9) == 3,
                     "a5c: wanted one beta row and three kappa rows")

    // (d) kappa[y4,x] fixed at 0.2: not free, so y4 is the anchor of F.
    par = dwls_identify(cx_spec(("covariate x", "gamma F x free", "kappa y1 x free",
                                 "kappa y2 x free", "kappa y3 x free", "kappa y4 x fixed 0.2")),
                        cx_cmap(1))
    r = cx_prow(par, "kappa[y4,x]")
    nf = nf + cx_chk(r > 0, "a5d: no kappa[y4,x] row")
    if (r > 0) {
        nf = nf + cx_chk(par.free[r] == 0 & par.value[r] == 0.2, "a5d: kappa[y4,x] is not fixed at 0.2")
    }
    return(nf)
}

// ---- r: refusals. Each function must abort; the caller reads the code. -------------

real scalar cx_r_anchor()
{
    // every indicator of F carries a kappa
    (void) dwls_identify(cx_spec(("covariate x", "gamma F x free", "kappa y1 x free",
                                  "kappa y2 x free", "kappa y3 x free", "kappa y4 x free")),
                         cx_cmap(1))
    return(0)
}

real scalar cx_r_beta0()
{
    // F1 has a free gamma and a free kappa on every indicator; its only beta is fixed at 0, so
    // it is not present and F2 (kappa-free) is not reachable from F1
    (void) dwls_identify(cx_spec6(("covariate x", "gamma F1 x free", "kappa y1 x free",
                                   "kappa y2 x free", "kappa y3 x free", "psi F1 F2 free",
                                   "beta F2 F1 fixed 0")),
                         cx_cmap6(1))
    return(0)
}

real scalar cx_r_beta_rev()
{
    // `beta F1 F2` makes F2 predict F1: F2 is not a descendant of F1, so F1 has no anchor
    (void) dwls_identify(cx_spec6(("covariate x", "gamma F1 x free", "kappa y1 x free",
                                   "kappa y2 x free", "kappa y3 x free", "beta F1 F2 free")),
                         cx_cmap6(1))
    return(0)
}

real scalar cx_r_anchor_lam0()
{
    // y4 has no kappa but its loading is fixed at 0: not a loading, so not an anchor
    string scalar s

    s = cx_join(("factor F", "ind y1 y2 y3 y4",
                 "lambda F y1 free", "lambda F y2 free", "lambda F y3 free", "lambda F y4 fixed 0",
                 "metric unitvariance F",
                 "covariate x", "gamma F x free", "kappa y1 x free", "kappa y2 x free",
                 "kappa y3 x free"))
    (void) dwls_identify(s, cx_cmap(1))
    return(0)
}

real scalar cx_r_cmap_nocov()
{
    (void) dwls_identify(cx_spec(J(1, 0, "")), cx_cmap(1))
    return(0)
}

real scalar cx_r_cov_umap()
{
    (void) dwls_identify(cx_spec1(), cx_umap())
    return(0)
}

real scalar cx_r_wrongq()
{
    // one covariate declared, a q = 2 map
    (void) dwls_identify(cx_spec1(), cx_cmap(2))
    return(0)
}

real scalar cx_r_two_one_line()
{
    (void) dwls_identify(cx_spec(("covariate x z", "gamma F x free")), cx_cmap(2))
    return(0)
}

real scalar cx_r_two_lines()
{
    (void) dwls_identify(cx_spec(("covariate x", "covariate z", "gamma F x free")), cx_cmap(1))
    return(0)
}

real scalar cx_r_groups()
{
    struct dwlsmv_statmap vector m
    string scalar s
    real scalar g

    m = dwlsmv_statmap(1, 2)
    for (g = 1; g <= 2; g++) {
        m[g] = dwlsmv_statmap_build(J(2, 1, 3), 1)
    }
    s = cx_join(("groups 2", "group 1", "covariate x", "factor F", "ind y1 y2",
                 "lambda F y1 fixed 0.7", "lambda F y2 fixed 0.7", "metric unitvariance F",
                 "group 2", "factor F", "ind y1 y2", "lambda F y1 free", "lambda F y2 free",
                 "psi F F fixed 1", "alpha F free",
                 "tie tau 1 y1 1 tau 2 y1 1", "tie tau 1 y1 2 tau 2 y1 2",
                 "tie tau 1 y2 1 tau 2 y2 1", "tie tau 1 y2 2 tau 2 y2 2"))
    (void) dwls_identify(s, m)
    return(0)
}

real scalar cx_r_gamma_nocov()
{
    (void) dwls_identify(cx_spec(("gamma F x free")), cx_cmap(1))
    return(0)
}

real scalar cx_r_kappa_badind()
{
    (void) dwls_identify(cx_spec(("covariate x", "gamma F x free", "kappa y9 x free")), cx_cmap(1))
    return(0)
}

real scalar cx_r_gamma_badfac()
{
    (void) dwls_identify(cx_spec(("covariate x", "gamma G x free")), cx_cmap(1))
    return(0)
}

real scalar cx_r_gamma_badcov()
{
    (void) dwls_identify(cx_spec(("covariate x", "gamma F z free")), cx_cmap(1))
    return(0)
}

real scalar cx_r_gamma_dup()
{
    (void) dwls_identify(cx_spec(("covariate x", "gamma F x free", "gamma F x fixed 0.2")), cx_cmap(1))
    return(0)
}

real scalar cx_r_cov_name()
{
    // the covariate has an indicator's name
    (void) dwls_identify(cx_spec(("covariate y1", "gamma F y1 free")), cx_cmap(1))
    return(0)
}

real scalar cx_r_cov_name_fac()
{
    // the covariate has a factor's name
    (void) dwls_identify(cx_spec(("covariate F", "gamma F F free")), cx_cmap(1))
    return(0)
}

real scalar cx_r_kappa_short()
{
    (void) dwls_identify(cx_spec(("covariate x", "gamma F x free", "kappa y1 x")), cx_cmap(1))
    return(0)
}

real scalar cx_r_kappa_status()
{
    (void) dwls_identify(cx_spec(("covariate x", "gamma F x free", "kappa y1 x maybe")), cx_cmap(1))
    return(0)
}

real scalar cx_r_kappa_nonnum()
{
    (void) dwls_identify(cx_spec(("covariate x", "gamma F x free", "kappa y1 x free abc")), cx_cmap(1))
    return(0)
}

real scalar cx_r_st_umap()
{
    // conditional maps and a covariate, but an unconditional Stage 1 result
    struct dwlsmv_stage1 vector st

    st = dwlsmv_stage1(1, 1)
    st[1].map = cx_umap()
    (void) dwls_identify(cx_spec1(), cx_cmap(1), st)
    return(0)
}

// ---- L: line numbers. Each function aborts; the caller reads the log. ----------------

// L1. Original lines: 1 comment, 2 factor, 3 ind, 4 blank, 5 covariate, 6 gamma, 7 kappa,
// 8 lambda with too few tokens. dwlsmv_id_parse() must say 8. (Numbering the lines that reach it
// would say 3; numbering the lines with covariate lines kept but blanks dropped would say 6.)
real scalar cx_l1()
{
    (void) dwls_identify(cx_join(("# model", "factor F", "ind y1 y2 y3 y4", "",
                                  "covariate x", "gamma F x free", "kappa y3 x free 0.1",
                                  "lambda F y1")),
                         cx_cmap(1))
    return(0)
}

// L2. Line 6 is a gamma line with a name missing.
real scalar cx_l2()
{
    (void) dwls_identify(cx_join(("# model", "factor F", "ind y1 y2 y3 y4", "",
                                  "covariate x", "gamma F x",
                                  "lambda F y1 free", "lambda F y2 free", "lambda F y3 free",
                                  "lambda F y4 free", "metric unitvariance F")),
                         cx_cmap(1))
    return(0)
}

// L3. Line 8 is a kappa line naming an undeclared indicator; a comment sits on line 7.
real scalar cx_l3()
{
    (void) dwls_identify(cx_join(("# model", "factor F", "ind y1 y2 y3 y4", "",
                                  "covariate x", "gamma F x free", "# direct effect",
                                  "kappa y9 x free",
                                  "lambda F y1 free", "lambda F y2 free", "lambda F y3 free",
                                  "lambda F y4 free", "metric unitvariance F")),
                         cx_cmap(1))
    return(0)
}

// L4. The second covariate directive is on line 7; the first on line 5.
real scalar cx_l4()
{
    (void) dwls_identify(cx_join(("# model", "factor F", "ind y1 y2 y3 y4", "",
                                  "covariate x", "gamma F x free", "covariate z")),
                         cx_cmap(1))
    return(0)
}

// L5. No covariate. A comment, a blank and a comment sit above the error on line 6.
real scalar cx_l5()
{
    (void) dwls_identify(cx_join(("# c", "factor F", "", "ind y1 y2 y3 y4", "# note",
                                  "lambda F y1")),
                         cx_umap())
    return(0)
}

// L6. CRLF line ends, error on line 3.
real scalar cx_l6()
{
    string scalar s, crlf

    crlf = char(13) + char(10)
    s    = "factor F" + crlf + "ind y1 y2 y3 y4" + crlf + "lambda F y1" + crlf
    (void) dwls_identify(s, cx_umap())
    return(0)
}

// ---- a6: the screen follows a path of two betas (F1 -> F2 -> F3) --------------------------
// Three factors, two indicators each (F1: y1 y2, F2: y3 y4, F3: y5 y6), unit variance on all
// three, conditional map cx_cmap6(1). facline is the factor declaration, so the caller sets the
// order of the factors; the rest of the text is appended by the caller.
string scalar cx_spec3(string scalar facline, string rowvector extra)
{
    string rowvector base

    base = (facline, "ind y1 y2 y3 y4 y5 y6",
            "lambda F1 y1 free", "lambda F1 y2 free",
            "lambda F2 y3 free", "lambda F2 y4 free",
            "lambda F3 y5 free", "lambda F3 y6 free",
            "metric unitvariance F1", "metric unitvariance F2", "metric unitvariance F3")
    return(cx_join((base, extra)))
}

// gamma on F1 only; kappa free on every indicator of F1 and F2; none on y5 y6.
string rowvector cx_x3()
{
    return(("covariate x", "gamma F1 x free", "kappa y1 x free", "kappa y2 x free",
            "kappa y3 x free", "kappa y4 x free"))
}

real scalar cx_a6()
{
    struct dwlsmv_parmap scalar par
    real scalar nf, r, s, j
    string rowvector facs

    nf   = 0
    facs = ("factor F1 F2 F3", "factor F1 F3 F2")
    for (j = 1; j <= 2; j++) {
        // (a) and (a2, F2 declared last): F1 reaches F3 through F2 only, by two present betas.
        // The screen accepts. The rank check is expected to accept too: at theta0 (beta = 0)
        // the column of gamma vanishes (derivation, free-zero case), theta2 passes.
        par = dwls_identify(cx_spec3(facs[j], (cx_x3(), "beta F2 F1 free", "beta F3 F2 free")),
                            cx_cmap6(1))
        nf = nf + cx_chk(sum(par.ptype :== 3) == 2 & sum(par.ptype :== 8) == 1 &
                         sum(par.ptype :== 9) == 4,
                         "a6" + strofreal(j) + ": wanted two beta rows, one gamma row, four kappa rows")
        r = cx_prow(par, "beta[F3,F2]")
        nf = nf + cx_chk(r > 0, "a6" + strofreal(j) + ": no beta[F3,F2] row")
        if (r > 0) {
            nf = nf + cx_chk(par.free[r] == 1, "a6" + strofreal(j) + ": beta[F3,F2] is not free")
            // theta0 is deficient, so the check passed at theta2, which becomes the start.
            nf = nf + cx_chk(par.value[r] != 0,
                             "a6" + strofreal(j) + ": beta[F3,F2] still starts at 0 (theta2 not used?)")
        }
        s = cx_prow(par, "beta[F2,F1]")
        nf = nf + cx_chk(s > 0, "a6" + strofreal(j) + ": no beta[F2,F1] row")
    }
    return(nf)
}

// (b) the second step fixed at 0 is absent: R(F1) = {F1, F2}, both with a kappa on every item.
real scalar cx_l8()
{
    (void) dwls_identify(cx_spec3("factor F1 F2 F3",
                                  (cx_x3(), "beta F2 F1 free", "beta F3 F2 fixed 0")),
                         cx_cmap6(1))
    return(0)
}

// (c) the path runs F3 -> F2 -> F1: F1 reaches nothing, F3 is its ancestor.
real scalar cx_l9()
{
    (void) dwls_identify(cx_spec3("factor F1 F2 F3",
                                  (cx_x3(), "beta F2 F3 free", "beta F1 F2 free")),
                         cx_cmap6(1))
    return(0)
}

// L7. A covariate under 'groups': the refusal cites the covariate on line 5 and the groups text
// on line 3 (blank and comment lines count). Raised by dwlsmv_id_cxsplit() before any map is read.
real scalar cx_l7()
{
    (void) dwls_identify(cx_join(("# model", "", "groups 2", "group 1", "covariate x")),
                         cx_cmap(1))
    return(0)
}

// Does the log at f contain the text has? Returns the number of failures.
real scalar cx_msg(string scalar f, string scalar has, string scalar tag)
{
    string scalar txt

    txt = invtokens(cat(f)', char(10))
    if (!strpos(txt, has)) {
        printf("{err}test_condx_identify: case %s: the log lacks: %s\n", tag, has)
        return(1)
    }
    return(0)
}

// Does the log at f NOT contain the text hasnot? Returns the number of failures.
real scalar cx_nomsg(string scalar f, string scalar hasnot, string scalar tag)
{
    string scalar txt

    if (hasnot == "") {
        return(0)
    }
    txt = invtokens(cat(f)', char(10))
    if (strpos(txt, hasnot)) {
        printf("{err}test_condx_identify: case %s: the log has: %s\n", tag, hasnot)
        return(1)
    }
    return(0)
}

end

local nfail 0
tempfile L

* ---- a1 to a4: the table ---------------------------------------------------------
foreach c in cx_a1 cx_a2 cx_a3 cx_a4 cx_a5 cx_a6 {
    capture noisily mata: st_local("nok", strofreal(`c'()))
    if _rc {
        display as error "test_condx_identify: `c'() aborted, rc " _rc ", wanted 0"
        local nfail = `nfail' + 1
    }
    else {
        display as text "  `c'(): `nok' failure(s)"
        local nfail = `nfail' + `nok'
    }
}

* ---- r: every refusal is rc 3498 exactly --------------------------------------------
local calls cx_r_anchor cx_r_anchor_lam0 cx_r_beta0 cx_r_beta_rev cx_r_cmap_nocov cx_r_cov_umap
local calls `calls' cx_r_wrongq cx_r_two_one_line cx_r_two_lines cx_r_groups cx_r_gamma_nocov
local calls `calls' cx_r_kappa_badind cx_r_gamma_badfac cx_r_gamma_badcov cx_r_gamma_dup
local calls `calls' cx_r_cov_name cx_r_cov_name_fac cx_r_kappa_short cx_r_kappa_status
local calls `calls' cx_r_kappa_nonnum cx_r_st_umap
foreach c of local calls {
    capture noisily mata: (void) `c'()
    local rc = _rc
    display as text "  `c'(): rc `rc'"
    if `rc' != 3498 {
        display as error "test_condx_identify: `c'() gave rc `rc', wanted 3498"
        local nfail = `nfail' + 1
    }
}

* ---- L: line numbers, read back from the log ------------------------------------------
* w: the text the log must contain; x: the text it must not (the number an earlier numbering gave).
local n 0
local ++n
local c`n' cx_l1()
local w`n' dwls_identify: line 8:
local x`n' dwls_identify: line 3:
local ++n
local c`n' cx_l1()
local w`n' needs two names then
local x`n' dwls_identify: line 6:
local ++n
local c`n' cx_l2()
local w`n' dwls_identify: line 6:
local x`n' dwls_identify: line 4:
local ++n
local c`n' cx_l3()
local w`n' dwls_identify: line 8:
local x`n' dwls_identify: line 5:
local ++n
local c`n' cx_l4()
local w`n' dwls_identify: line 7: a second
local x`n' dwls_identify: line 4:
local ++n
local c`n' cx_l5()
local w`n' dwls_identify: line 6:
local x`n' dwls_identify: line 3:
local ++n
local c`n' cx_l6()
local w`n' dwls_identify: line 3:
local x`n'
local ++n
local c`n' cx_r_anchor()
local w`n' dwls_identify: factor F has a free gamma on x
local x`n'
local ++n
local c`n' cx_r_anchor_lam0()
local w`n' dwls_identify: factor F has a free gamma on x
local x`n'
local ++n
local c`n' cx_r_beta0()
local w`n' dwls_identify: factor F1 has a free gamma on x
local x`n'
local ++n
local c`n' cx_r_beta_rev()
local w`n' dwls_identify: factor F1 has a free gamma on x
local x`n'
local ++n
local c`n' cx_r_cov_name()
local w`n' has the name of a factor or an indicator
local x`n'
local ++n
local c`n' cx_r_cov_name_fac()
local w`n' has the name of a factor or an indicator
local x`n'
local ++n
local c`n' cx_l8()
local w`n' dwls_identify: factor F1 has a free gamma on x
local x`n'
local ++n
local c`n' cx_l9()
local w`n' dwls_identify: factor F1 has a free gamma on x
local x`n'
local ++n
local c`n' cx_l7()
local w`n' dwls_identify: line 5: a covariate together with a
local x`n' (line 2)
local ++n
local c`n' cx_l7()
local w`n' text (line 3)
local x`n' (line 4)

forvalues i = 1/`n' {
    quietly log using "`L'_`i'.log", text replace name(gt)
    capture noisily mata: (void) `c`i''
    local rc = _rc
    quietly log close gt
    display as text "  `c`i'': rc `rc'"
    if `rc' != 3498 {
        display as error "test_condx_identify: `c`i'' gave rc `rc', wanted 3498"
        local nfail = `nfail' + 1
    }
    else {
        local tmpL "`L'_`i'.log"
        local want "`w`i''"
        local nothave "`x`i''"
        capture noisily mata: st_local("nm", strofreal(cx_msg(st_local("tmpL"), st_local("want"), "L`i'") + cx_nomsg(st_local("tmpL"), st_local("nothave"), "L`i'")))
        if _rc {
            display as error "test_condx_identify: reading the log of `c`i'' failed, rc " _rc
            local nfail = `nfail' + 1
        }
        else local nfail = `nfail' + `nm'
    }
}

if `nfail' {
    display as error _n "test_condx_identify: `nfail' failure(s)"
    exit 9
}

display as result "test_condx_identify: dwls_identify accepts a declared covariate on a conditional map, orders gamma and kappa between beta and psi, refuses the unspecified cases at rc 3498 naming the line as written, and a model without covariates is unchanged"
exit 0
