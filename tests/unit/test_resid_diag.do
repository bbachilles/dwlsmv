*! Residual diagnostics (SRMR, normalized residuals, their ranking) against the lavaan
*! fixtures of fixtures/generate/99t_lavaan_srmr.R.
*!
*! Drafted 2026-11-22 by a read-only agent; executed and fixed by the main session (StataMP 18,
*! 2026-11-22): passes. Main-session changes: the fits file's lowercased n; the printed-part split.
*!
*! Derivation @sec-resid-diag (SR1-SR7), owner decisions E-117,
*! package 3b-3d (O1-O7). Subject: src/mata/dwlsmv_resid.mata, dwlsmv_resid_rank.mata.
*!
*!   A  Arithmetic at 1e-12. A Stage 1 struct is hand-built from the FIXTURE's own r,
*!      gamma_kk, N and the fixture rho, so only the arithmetic is tested, not Stage 1
*!      or a fit: E, Z (where eligible), SRMR, the Z missing pattern.
*!   B  End to end, every exported fit: dwlsmv's own Stage 1, own fit, own repaired-pair
*!      set, against lavaan. Tolerances are the existing tiers of the generator header
*!      (rho 1e-5 abs + 1e-4 rel; r and gamma_kk 1e-8 + 1e-8; e within the fixture's e_tol;
*!      srmr within srmr_tol). z is NOT compared directly (owner 2(a)): its parts are, and
*!      z is checked against its own formula at 1e-12. End-to-end r and gamma_kk are
*!      PRINTED, NOT ASSERTED where lavaan's own rho-hat misses the tier: case_E both
*!      (test_gamma.do 5d, 2026-09-25), case_F both (owner 2026-11-22, after the cause
*!      check in E-117), emptycell gamma_kk only
*!      (test_lavaan_emptycell.do check 3; its r passes and stays asserted). Where gamma_kk
*!      is printed, gamma_kk at lavaan's own rho-hat is asserted at 1e-12 (as 5c). Where r is
*!      printed, e and srmr are asserted against tolerances REBUILT from the observed r gap
*!      (rho term stays at its tier); wherever either is printed, the ranking intervals are
*!      rebuilt from this run's observed gaps for the printed parts. The ranking is asserted
*!      only where the (fixture or rebuilt) intervals determine it.
*!   C  Hand-built pair at the bound (no fixture has one; lavaan has no keepbound), repaired
*!      and bound pairs together, and ties in the ranking.
*!   D  Every guard, _rc == 3498 exactly, and the phrase of its errprintf.
*!
*! Fit routes. EFA fits: dwlsmv ..., efa() and rho from e(Lambda) e(Phi) (rotation-invariant,
*! SR6). CFA fits: e(Lambda) and e(Phi) are posted only for efa(); e(b)'s labels are not used.
*! The CFA implied rho comes from the Mata layer on dwlsmv's own Stage 1: dwls_identify on a
*! canonical specification (all loadings free, metric unitvariance, f1 ~~ f2 for m = 2),
*! dwls_fit, dwls_model, as test_lavaan_emptycell.do check 4.
*!
*! MUTATION LIST (one-line diffs; the part that must kill each):
*!   M1  dwlsmv_resid: divisor (P + p) -> P                         A srmr, B srmr, C srmr
*!   M2  dwlsmv_resid: gk :/ N -> gk :/ (N - 1)                      A z, B z formula, C z
*!   M3  dwlsmv_resid: gk[k] = Gamma[k, k] (row k, not T + k)        A, B (thresholds shift it)
*!   M4  dwlsmv_resid: drop "& st.bound[k] == 0"                     C (Z present at bound)
*!   M5  dwlsmv_resid: drop "repaired[k] == 0 &"                     A and B on R_cfa1, C
*!   M6  dwlsmv_resid_rank: order idx (1,2,3) -> (1,-2,3)-style, ties
*!       by larger i (e.g. order((-az, -ii, jj), ...))               C ties (and B where tied)
*!   M7  dwlsmv_resid: e = rho - r (sign flipped)                    A e, B e
*!   M8  dwlsmv_resid: Zo = J(p, p, 0) (diagonal 0 not missing)      A, C (diagonal missing)
*!   M9  dwlsmv_resid_rank: abs(Z[j,i]) -> Z[j,i] (signed order)     C (mixed signs), B
*!   M10 dwlsmv_resid: add the lavaan 1e-5 rule, z = e where |e|<1e-5   B z formula (n_small_e)
*!       (only bites where a fit has |e| < 1e-5: the generator's n_small_e says where)
*!   M11 dwlsmv_resid: Eo[vj,vi] not written (E lower triangle lost) A E symmetry
*!   M12 dwlsmv_resid: srmr over eligible pairs only (sum over z present) C (bound pair
*!       included in SRMR), A R_cfa1
*!   M13 dwlsmv_resid: drop the "gk :<= 0" guard                     D case 9, 10
*!   M14 dwlsmv_resid_rank: ranking over |e| instead of |z|          B where |e| and |z| order differ, C
*! End-of-module review mutants (2026-11-22, run on the whole branch): the pair-index guard
*! without "vi >= vj" is killed by D cases 21-22; the Stage-1-length guard disabled, by D 23-24.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

* builds build/ldwlsmv.mlib and leaves every function and struct in memory
quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index
do "tests/dwlsmv_assert.mata"

* run one command, capturing what it prints; r(rc), and r(found) for a phrase
capture program drop tr_run
program tr_run, rclass
    gettoken phrase 0 : 0
    tempfile lg
    quietly log using "`lg'", text replace name(trlog)
    capture noisily `0'
    local rc = _rc
    quietly log close trlog
    mata: st_local("found", strofreal(strpos(invtokens(cat(st_local("lg"))', " "), st_local("phrase")) > 0))
    return scalar rc    = `rc'
    return scalar found = `found'
end

* ---- fixtures: stringcols(_all), delimiters(","), strtoreal() (a project rule) --------
import delimited using "fixtures/expected/lavaan_srmr.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: TR_F = st_sdata(., "fit")
mata: TR_Q = st_sdata(., "quantity")
mata: TR_I = strtoreal(st_sdata(., "i"))
mata: TR_J = strtoreal(st_sdata(., "j"))
mata: TR_V = strtoreal(st_sdata(., "value"))
import delimited using "fixtures/expected/lavaan_srmr_fits.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
* columns: fit data items model m converged df n attempted (import delimited lowercases N)
mata: TR_FT = st_sdata(., ("fit", "data", "items", "model", "m", "converged", "df", "n", "attempted"))
mata: TR_FAIL = J(0, 1, "")

mata:

void tr_fail(string scalar msg)
{
    external string colvector TR_FAIL

    TR_FAIL = TR_FAIL \ msg
}

void tr_check(real scalar ok, string scalar msg)
{
    if (!ok) {
        tr_fail(msg)
    }
}

// the fixture's scalar quantity, or missing if there is not exactly one row
real scalar tr_sc(string scalar fid, string scalar q)
{
    external string colvector TR_F, TR_Q
    external real colvector TR_I, TR_V
    real colvector ix

    ix = selectindex(TR_F :== fid :& TR_Q :== q :& TR_I :>= .)
    if (rows(ix) != 1) {
        return(.)
    }
    return(TR_V[ix])
}

// the fixture's pair quantity as a P x 1 vector in the spec's pair order (i outer)
real colvector tr_pr(string scalar fid, string scalar q, real scalar p)
{
    external string colvector TR_F, TR_Q
    external real colvector TR_I, TR_J, TR_V
    real colvector out, ix
    real scalar a, i, j

    out = J(p * (p - 1) / 2, 1, .)
    ix  = selectindex(TR_F :== fid :& TR_Q :== q :& TR_I :< .)
    for (a = 1; a <= rows(ix); a++) {
        i = TR_I[ix[a]]
        j = TR_J[ix[a]]
        out[(i - 1) * p - i * (i - 1) / 2 + (j - i)] = TR_V[ix[a]]
    }
    return(out)
}

string scalar tr_ft(string scalar fid, real scalar col)
{
    external string matrix TR_FT
    real colvector ix

    ix = selectindex(TR_FT[., 1] :== fid)
    if (rows(ix) != 1) {
        return("")
    }
    return(TR_FT[ix, col])
}

// "y1-y6" -> "y1 y2 y3 y4 y5 y6"
string scalar tr_items(string scalar rng)
{
    real scalar pos, a, b, k
    string scalar out

    pos = strpos(rng, "-")
    if (substr(rng, 1, 1) != "y" | pos < 3) {
        _error("tr_items: expects y<a>-y<b>")
    }
    a = strtoreal(substr(rng, 2, pos - 2))
    b = strtoreal(substr(rng, pos + 2, .))
    out = ""
    for (k = a; k <= b; k++) {
        out = out + (k > a ? " " : "") + "y" + strofreal(k)
    }
    return(out)
}

// the off-diagonal of a p x p matrix at each pair of the map, P x 1, pair order
real colvector tr_lt(real matrix M, struct dwlsmv_statmap scalar map)
{
    real colvector out
    real scalar k, pos

    out = J(map.P, 1, .)
    for (k = 1; k <= map.P; k++) {
        pos    = map.T + k
        out[k] = M[map.var_i[pos], map.var_j[pos]]
    }
    return(out)
}

// a hand-built Stage 1: only what dwlsmv_resid() reads is meaningful
struct dwlsmv_stage1 scalar tr_mkst(real colvector r, real colvector gk, real scalar N, real scalar p)
{
    struct dwlsmv_stage1 scalar st
    real scalar T, P, D

    st.map    = dwlsmv_statmap_build(J(p, 1, 3))
    T         = st.map.T
    P         = st.map.P
    D         = st.map.D
    st.s      = J(T, 1, 0.1) \ r
    st.Gamma  = diag(J(T, 1, 1) \ gk)
    st.N      = N
    st.psibar = J(D, 1, 0)
    st.bound  = J(P, 1, 0)
    return(st)
}

// pair index (spec order, i outer) of (i, j), i < j
real scalar tr_pk(real scalar i, real scalar j, real scalar p)
{
    return((i - 1) * p - i * (i - 1) / 2 + (j - i))
}

// ===================== A. arithmetic from the fixture's own values ======================
void tr_a(string scalar fid)
{
    struct dwlsmv_stage1 scalar st
    real scalar p, P, N, srmr, k, el
    real colvector r, gk, rho, e_fx, z_fx, rep, elig, eo, zo
    real matrix E, Z

    p    = tr_sc(fid, "p")
    P    = tr_sc(fid, "P")
    N    = tr_sc(fid, "N")
    r    = tr_pr(fid, "r", p)
    rho  = tr_pr(fid, "rho", p)
    gk   = tr_pr(fid, "gamma_kk", p)
    e_fx = tr_pr(fid, "e", p)
    z_fx = tr_pr(fid, "z_hand", p)
    rep  = tr_pr(fid, "repaired", p)
    elig = tr_pr(fid, "eligible", p)
    if (missing(r) | missing(rho) | missing(gk) | missing(e_fx) | missing(z_fx) | missing(rep) | missing(elig) | rows(r) != P) {
        tr_fail("A " + fid + ": a fixture column is missing")
        return
    }
    st   = tr_mkst(r, gk, N, p)
    srmr = dwlsmv_resid(st, rho, rep, E, Z)
    eo   = tr_lt(E, st.map)
    zo   = tr_lt(Z, st.map)

    tr_check(dwlsmv_gapratio(eo, e_fx, 1e-12, 1e-12) <= 1, "A " + fid + ": e differs from the fixture at 1e-12")
    tr_check(dwlsmv_gapratio(srmr, tr_sc(fid, "srmr_hand"), 1e-12, 1e-12) <= 1, "A " + fid + ": srmr differs from srmr_hand at 1e-12")
    tr_check(max(abs(E - E')) == 0, "A " + fid + ": E is not symmetric")
    tr_check(max(abs(diagonal(E))) == 0, "A " + fid + ": the diagonal of E is not 0")
    tr_check(missing(diagonal(Z)) == p, "A " + fid + ": the diagonal of Z is not missing")
    tr_check(missing(Z) == p + 2 * sum(elig :== 0), "A " + fid + ": Z missing pattern: count")
    tr_check(all((zo :>= .) :== (elig :== 0)), "A " + fid + ": Z is missing at a different set of pairs than the fixture's ineligible")
    el = selectindex(elig :== 1)
    if (rows(el) > 0) {
        tr_check(dwlsmv_gapratio(zo[el], z_fx[el], 1e-12, 1e-12) <= 1, "A " + fid + ": z differs from z_hand at 1e-12")
    }
    // Z symmetric (a missing equals a missing in Mata)
    tr_check(Z == Z', "A " + fid + ": Z is not symmetric")
    printf("{txt}  A %-7s srmr %12.10f (fixture %12.10f)  pairs %g, repaired %g\n", fid, srmr, tr_sc(fid, "srmr_hand"), P, sum(rep))
}

// ===================== B. end to end against lavaan =====================================
// Gamma from production's dwls_gamma() with each pair's rho and information taken at the given
// rho (lavaan's rho-hat): test_lavaan_emptycell.do's tle_gamma() at rep_plus = 0, which is
// production's construction (a) for a repaired pair too.
real matrix tr_gamma_at(real matrix Y, struct dwlsmv_statmap scalar m,
                        real colvector s, real colvector rho)
{
    struct dwlsmv_pairfit colvector pairs
    struct dwlsmv_pairfit scalar f
    real colvector ti, tj
    real matrix n
    real scalar k, pos, vi, vj, N

    N     = rows(Y)
    pairs = dwlsmv_pairfit(m.P, 1)
    for (k = 1; k <= m.P; k++) {
        pos = m.T + k
        vi  = m.var_i[pos]
        vj  = m.var_j[pos]
        ti  = s[|m.Coff[vi] + 1 \ m.Coff[vi] + m.C[vi]|]
        tj  = s[|m.Coff[vj] + 1 \ m.Coff[vj] + m.C[vj]|]
        n   = dwlsmv_crosstab(Y[., vi], Y[., vj], m.K[vi], m.K[vj])
        // a full dwls_polychoric() fit first, as tle_gamma() (an empty element crashed Stata 18)
        f     = dwls_polychoric(Y[., vi], Y[., vj], ti, tj)
        f.rho = rho[k]
        dwlsmv_pairinfo(f, ti, tj, rho[k], n :/ N)
        f.i = vi
        f.j = vj
        pairs[k] = f
    }
    return(dwls_gamma(Y, m, pairs, dwlsmv_psi_tau_stack(Y, m)))
}

string scalar tr_cfa_spec(string scalar items, real scalar m)
{
    string rowvector it
    string scalar s
    real scalar p, n1, j

    it = tokens(items)
    p  = cols(it)
    n1 = floor(p / 2)
    if (m == 1) {
        s = "factor f" + char(10) + "ind " + items + char(10)
        for (j = 1; j <= p; j++) {
            s = s + "lambda f " + it[j] + " free" + char(10)
        }
        return(s + "metric unitvariance f" + char(10))
    }
    s = "factor f1 f2" + char(10) + "ind " + items + char(10)
    for (j = 1; j <= p; j++) {
        s = s + "lambda " + (j <= n1 ? "f1 " : "f2 ") + it[j] + " free" + char(10)
    }
    return(s + "psi f1 f2 free" + char(10) + "metric unitvariance f1" + char(10) + "metric unitvariance f2" + char(10))
}

void tr_b(string scalar fid, string scalar items, string scalar model, real scalar m,
          string scalar touse, string scalar dat)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    pointer() scalar pst, pY
    string rowvector it
    string matrix rn
    real matrix L, Ph, S, E, Z, Rk
    real colvector r, gk, rho, rep, r_fx, rho_fx, gk_fx, e_fx, etol_fx, se_fx, elig, rk_fx
    real colvector sig, eo, zo, tol, zlo, zhi, DE, EPSG, den, el, ord_fx, k3, rest, ours, a
    real scalar p, P, N, k, srmr, srmr_fx, srmr_tol, df, set_det, ord_det
    real scalar rg_rho, rg_r, rg_g, rg_e, rg_s, rg_g2, pr_r, pr_g, rebuilt
    real matrix G2

    p      = tr_sc(fid, "p")
    P      = tr_sc(fid, "P")
    df     = tr_sc(fid, "df")
    // End-to-end r and gamma_kk printed, not asserted, where lavaan's own rho-hat misses the tier:
    //   case_E    r and gamma_kk  test_gamma.do 5d (E-060)
    //   case_F    r and gamma_kk  owner 2026-11-22, after the cause check of
    //                             E-117 (comparator convergence:
    //                             dwlsmv's r within 8.2e-14 of the root, lavaan's up to 3.3e-8 off)
    //   emptycell gamma_kk only   test_lavaan_emptycell.do check 3; r passes and stays asserted
    // Where gamma_kk is printed, gamma_kk at lavaan's own rho-hat is asserted at 1e-12 (5c).
    pr_r   = (dat == "case_E" | dat == "case_F")
    pr_g   = (dat == "case_E" | dat == "case_F" | dat == "emptycell")
    rebuilt = (pr_r | pr_g)
    it     = tokens(items)

    dwlsmv_m11_stage1(touse, items, 0)
    pst = findexternal("dwlsmv_m11_st")
    pY  = findexternal("dwlsmv_m11_Y")
    st  = *pst
    rep = dwlsmv_emptycell_pairs(*pY, st.map)

    // ---- model-implied rho at theta-hat, in the map's pair order
    if (model == "efa") {
        L  = st_matrix("e(Lambda)")
        rn = st_matrixrowstripe("e(Lambda)")
        if (rows(L) != p | cols(L) != m) {
            tr_fail("B " + fid + ": e(Lambda) is not p x m")
            dwlsmv_m11_clear()
            return
        }
        for (k = 1; k <= p; k++) {
            if (rn[k, 2] != it[k]) {
                tr_fail("B " + fid + ": e(Lambda) row order differs from the item order")
                dwlsmv_m11_clear()
                return
            }
        }
        if (m == 1) {
            Ph = I(1)
        }
        else {
            Ph = st_matrix("e(Phi)")
        }
        S   = L * Ph * L'
        rho = tr_lt(S, st.map)
    }
    else {
        par = dwls_identify(tr_cfa_spec(items, m), st.map)
        fit = dwls_fit(st, par)
        tr_check(fit.converged == 1, "B " + fid + ": the CFA fit did not converge")
        sig = dwls_model(fit.theta, par, st.map)
        rho = sig[|st.map.T + 1 \ st.map.D|]
    }
    srmr = dwlsmv_resid(st, rho, rep, E, Z)

    // ---- fixture
    N       = tr_sc(fid, "N")
    r_fx    = tr_pr(fid, "r", p)
    rho_fx  = tr_pr(fid, "rho", p)
    gk_fx   = tr_pr(fid, "gamma_kk", p)
    e_fx    = tr_pr(fid, "e", p)
    etol_fx = tr_pr(fid, "e_tol", p)
    se_fx   = tr_pr(fid, "se_r", p)
    elig    = tr_pr(fid, "eligible", p)
    srmr_fx = tr_sc(fid, "srmr_hand")
    srmr_tol = tr_sc(fid, "srmr_tol")
    r       = st.s[|st.map.T + 1 \ st.map.D|]
    gk      = J(P, 1, .)
    for (k = 1; k <= P; k++) {
        gk[k] = st.Gamma[st.map.T + k, st.map.T + k]
    }

    // exact: N, p, P, the repaired set, the eligible set
    tr_check(st.N == N, "B " + fid + ": N differs")
    tr_check(st.map.p == p & st.map.P == P, "B " + fid + ": p or P differs")
    tr_check(all(rep :== tr_pr(fid, "repaired", p)), "B " + fid + ": the repaired pairs differ from lavaan's")
    tr_check(all(elig :== (1 :- rep)), "B " + fid + ": eligible is not 1 - repaired (no pair is at the bound here)")
    tr_check(sum(st.bound) == 0, "B " + fid + ": a pair is at the bound; no fixture has one")

    // tiers
    rg_rho = dwlsmv_gapratio(rho, rho_fx, 1e-5, 1e-4)
    rg_r   = dwlsmv_gapratio(r, r_fx, 1e-8, 1e-8)
    rg_g   = dwlsmv_gapratio(gk, gk_fx, 1e-8, 1e-8)
    tr_check(rg_rho <= 1, "B " + fid + ": rho ratio " + strofreal(rg_rho))
    if (!pr_r) {
        tr_check(rg_r <= 1, "B " + fid + ": r ratio " + strofreal(rg_r))
    }
    if (!pr_g) {
        tr_check(rg_g <= 1, "B " + fid + ": gamma_kk ratio " + strofreal(rg_g))
    }
    rg_g2 = .
    if (pr_g) {
        G2 = tr_gamma_at(*pY, st.map, st.s, r_fx)
        rg_g2 = 0
        for (k = 1; k <= P; k++) {
            rg_g2 = max((rg_g2, dwlsmv_gapratio(G2[st.map.T + k, st.map.T + k], gk_fx[k], 1e-12, 1e-12)))
        }
        tr_check(rg_g2 <= 1, "B " + fid + ": gamma_kk at lavaan's rho-hat, ratio " + strofreal(rg_g2))
    }
    eo = tr_lt(E, st.map)
    zo = tr_lt(Z, st.map)
    // e and srmr: the fixture's tolerances; where r is printed, its tier term is replaced by the
    // observed r gap (a measured quantity, not a tolerance); the rho term stays at its tier
    if (pr_r) {
        tol      = (1e-5 :+ 1e-4 :* abs(rho_fx)) :+ abs(r - r_fx)
        srmr_tol = sqrt(sum(tol :^ 2) / (P + p))
    }
    else {
        tol = etol_fx
    }
    rg_e = max(abs(eo - e_fx) :/ tol)
    rg_s = abs(srmr - srmr_fx) / srmr_tol
    tr_check(rg_e <= 1, "B " + fid + ": e outside its tolerance, ratio " + strofreal(rg_e))
    tr_check(rg_s <= 1, "B " + fid + ": srmr outside its tolerance, ratio " + strofreal(rg_s))
    // z: its own formula, 1e-12, at the eligible pairs; missing exactly at the ineligible
    tr_check(all((zo :>= .) :== (elig :== 0)), "B " + fid + ": Z is missing at a different set of pairs than the fixture's ineligible")
    el = selectindex(elig :== 1)
    if (rows(el) > 0) {
        tr_check(dwlsmv_gapratio(zo[el], ((r - rho) :/ sqrt(gk :/ st.N))[el], 1e-12, 1e-12) <= 1, "B " + fid + ": z differs from e / sqrt(Gamma_kk / N)")
    }
    printf("{txt}  B %-7s srmr %10.8f vs %10.8f (ratio %5.3f); ratios rho %5.3f, r %5.3f%s, gamma_kk %5.3f%s, e %5.3f; repaired %g\n",
           fid, srmr, srmr_fx, rg_s, rg_rho, rg_r, (pr_r ? " (printed)" : ""), rg_g, (pr_g ? " (printed)" : ""), rg_e, sum(rep))
    if (pr_g) {
        printf("{txt}    gamma_kk at lavaan's rho-hat: ratio %g (1e-12 abs + rel, asserted)\n", rg_g2)
    }

    // ---- ranking (df > 0)
    if (df > 0) {
        rk_fx = tr_pr(fid, "rank", p)
        el    = selectindex(rk_fx :< .)
        ord_fx = el[order(rk_fx[el], 1)]
        k3    = ord_fx[|1 \ min((3, rows(ord_fx)))|]
        if (rows(ord_fx) >= 4) {
            rest = ord_fx[|4 \ rows(ord_fx)|]
        }
        else {
            rest = J(0, 1, .)
        }
        if (rebuilt) {
            // each printed part's tier term replaced by this run's observed gap
            if (pr_r) {
                DE = (1e-5 :+ 1e-4 :* abs(rho_fx)) :+ abs(r - r_fx)
            }
            else {
                DE = etol_fx
            }
            if (pr_g) {
                EPSG = abs(gk - gk_fx) :/ gk_fx
            }
            else {
                EPSG = (1e-8 :+ 1e-8 :* gk_fx) :/ gk_fx
            }
            a    = abs(e_fx) :- DE
            zlo  = (a :* (a :> 0)) :/ (se_fx :* sqrt(1 :+ EPSG))
            den  = 1 :- EPSG
            zhi  = J(P, 1, 1e300)
            for (k = 1; k <= P; k++) {
                if (den[k] > 0) {
                    zhi[k] = (abs(e_fx[k]) + DE[k]) / (se_fx[k] * sqrt(den[k]))
                }
            }
            set_det = 0
            if (rows(rest) == 0) {
                set_det = 1
            }
            else {
                set_det = (min(zlo[k3]) > max(zhi[rest]))
            }
            ord_det = 0
            if (set_det == 1) {
                ord_det = 1
                for (k = 1; k < rows(k3); k++) {
                    if (!(zlo[k3[k]] > zhi[k3[k + 1]])) {
                        ord_det = 0
                    }
                }
            }
        }
        else {
            set_det = tr_sc(fid, "top3_set_determined")
            ord_det = tr_sc(fid, "top3_order_determined")
        }
        Rk = dwlsmv_resid_rank(Z, 3)
        ours = J(rows(Rk), 1, .)
        for (k = 1; k <= rows(Rk); k++) {
            ours[k] = tr_pk(Rk[k, 1], Rk[k, 2], p)
        }
        tr_check(rows(ours) == rows(k3), "B " + fid + ": the ranking has " + strofreal(rows(ours)) + " rows, the fixture " + strofreal(rows(k3)))
        if (rows(ours) == rows(k3)) {
            if (set_det == 1) {
                tr_check(all(sort(ours, 1) :== sort(k3, 1)), "B " + fid + ": the top 3 set differs although the intervals determine it")
            }
            if (ord_det == 1) {
                tr_check(all(ours :== k3), "B " + fid + ": the top 3 order differs although the intervals determine it")
            }
        }
        printf("{txt}    ranking %-7s: set %s, order %s%s; ours %s, fixture %s\n", fid,
               (set_det == 1 ? "determined, asserted" : "UNDETERMINED, printed only"),
               (ord_det == 1 ? "determined, asserted" : "UNDETERMINED, printed only"),
               (rebuilt ? " (intervals rebuilt from this run's gaps)" : ""),
               invtokens(strofreal(ours')), invtokens(strofreal(k3')))
    }
    dwlsmv_m11_clear()
}

// ===================== C. hand-built: bound pair, repaired pair, ties ====================
// p = 4, pairs in order (1,2) (1,3) (1,4) (2,3) (2,4) (3,4)
void tr_c()
{
    struct dwlsmv_stage1 scalar st
    real colvector r, rho, gk, rep, z, ci, cj
    real matrix E, Z, Rk
    real scalar srmr, N, p, k

    p   = 4
    N   = 100
    r   = (.3 \ .4 \ .999 \ .2 \ .25 \ .35)
    rho = (.28 \ .41 \ .8 \ .22 \ .2 \ .33)
    gk  = (1 \ 2 \ 3 \ 4 \ 5 \ 6)
    rep = J(6, 1, 0)
    st  = tr_mkst(r, gk, N, p)

    // C1: nothing excluded
    srmr = dwlsmv_resid(st, rho, rep, E, Z)
    tr_check(abs(srmr - sqrt(sum((r - rho) :^ 2) / 10)) <= 1e-14, "C1: srmr")
    Rk = dwlsmv_resid_rank(Z, 3)
    tr_check(rows(Rk) == 3 & Rk == (1, 4 \ 2, 4 \ 1, 2), "C1: with no exclusion the bound-sized residual (1,4) ranks first, then (2,4), (1,2)")

    // C2: pair (1,4) at the bound
    st.bound[3] = 1
    srmr = dwlsmv_resid(st, rho, rep, E, Z)
    tr_check(missing(Z[1, 4]) & missing(Z[4, 1]), "C2: Z is not missing at the pair at the bound")
    tr_check(abs(E[1, 4] - (.999 - .8)) <= 1e-15 & E[4, 1] == E[1, 4], "C2: E is not present at the pair at the bound")
    tr_check(abs(srmr - sqrt(sum((r - rho) :^ 2) / 10)) <= 1e-14, "C2: SRMR must include the pair at the bound")
    tr_check(missing(Z) == 4 + 2, "C2: exactly the diagonal and one pair are missing in Z")
    Rk = dwlsmv_resid_rank(Z, 3)
    tr_check(rows(Rk) == 3 & Rk == (2, 4 \ 1, 2 \ 2, 3), "C2: ranking must skip the pair at the bound")
    Rk = dwlsmv_resid_rank(Z, 10)
    tr_check(rows(Rk) == 5 & cols(Rk) == 2, "C2: k larger than the eligible count returns the eligible count")
    tr_check(all((Rk[., 1] :== 1 :& Rk[., 2] :== 4) :== 0), "C2: (1,4) appears in the ranking")

    // C3: pair (2,4) repaired, (1,4) still at the bound
    rep[5] = 1
    srmr = dwlsmv_resid(st, rho, rep, E, Z)
    tr_check(missing(Z[2, 4]) & missing(Z[4, 2]) & missing(Z[1, 4]), "C3: Z must be missing at the repaired and the bound pair")
    tr_check(abs(srmr - sqrt(sum((r - rho) :^ 2) / 10)) <= 1e-14, "C3: SRMR must include the repaired pair")
    Rk = dwlsmv_resid_rank(Z, 3)
    tr_check(rows(Rk) == 3 & Rk == (1, 2 \ 2, 3 \ 3, 4), "C3: ranking must skip repaired and bound pairs")

    // C4: ties, and mixed signs: |z| equal at (1,2), (1,3), (2,3), (3,4); signs +,-,+,-
    z  = (.5 \ -.5 \ .1 \ .5 \ .2 \ -.5)
    ci = (1 \ 1 \ 1 \ 2 \ 2 \ 3)
    cj = (2 \ 3 \ 4 \ 3 \ 4 \ 4)
    Z  = J(4, 4, .)
    for (k = 1; k <= 6; k++) {
        Z[ci[k], cj[k]] = z[k]
        Z[cj[k], ci[k]] = z[k]
    }
    Rk = dwlsmv_resid_rank(Z, 6)
    tr_check(rows(Rk) == 6 & Rk[|1, 1 \ 4, 2|] == (1, 2 \ 1, 3 \ 2, 3 \ 3, 4), "C4: ties by smaller i then smaller j, signs ignored")
    tr_check(Rk[5, .] == (2, 4) & Rk[6, .] == (1, 4), "C4: the rest by descending |z|")

    // C5: no eligible pair
    Z = J(4, 4, .)
    Rk = dwlsmv_resid_rank(Z, 3)
    tr_check(rows(Rk) == 0 & cols(Rk) == 2, "C5: no eligible pair gives a 0 x 2 matrix")

    // C6: inputs are not modified
    rho = (.28 \ .41 \ .8 \ .22 \ .2 \ .33)
    rep = J(6, 1, 0)
    st  = tr_mkst(r, gk, N, p)
    srmr = dwlsmv_resid(st, rho, rep, E, Z)
    tr_check(rho == (.28 \ .41 \ .8 \ .22 \ .2 \ .33) & rep == J(6, 1, 0) & st.s[|st.map.T + 1 \ st.map.D|] == r & st.bound == J(6, 1, 0) & st.N == N, "C6: an input was modified")
    printf("{txt}  C hand-built cases done\n")
}

// ===================== D. every guard, one case per call =================================
void tr_d(real scalar c)
{
    struct dwlsmv_stage1 scalar st
    real colvector r, rho, gk, rep
    real matrix E, Z, Rk
    real scalar srmr, p, N, T

    p   = 4
    N   = 100
    r   = (.3 \ .4 \ .7 \ .2 \ .25 \ .35)
    rho = (.28 \ .41 \ .6 \ .22 \ .2 \ .33)
    gk  = (1 \ 2 \ 3 \ 4 \ 5 \ 6)
    rep = J(6, 1, 0)
    st  = tr_mkst(r, gk, N, p)
    T   = st.map.T

    if (c == 1)  st.map = dwlsmv_statmap_build(J(4, 1, 3), 1)
    if (c == 2)  st = tr_mkst(J(0, 1, .), J(0, 1, .), N, 1)
    if (c == 2)  rho = J(0, 1, .)
    if (c == 2)  rep = J(0, 1, 0)
    if (c == 3)  rho = rho[|1 \ 5|]
    if (c == 4)  rep = rep[|1 \ 5|]
    if (c == 5)  st.map.stat_type[T + 2] = 1
    if (c == 6)  st.s[T + 3] = .
    if (c == 7)  rho[2] = .
    if (c == 8)  st.Gamma[T + 4, T + 4] = .
    if (c == 9)  st.Gamma[T + 4, T + 4] = 0
    if (c == 10) st.Gamma[T + 4, T + 4] = -1
    if (c == 11) rep[2] = 2
    if (c == 12) rep[1] = .
    if (c == 13) st.bound[1] = 2
    if (c == 14) st.bound = J(5, 1, 0)
    if (c == 15) st.N = .
    if (c == 16) st.N = 0
    // 21-24: the pair-index and length guards of dwlsmv_resid() (review 2026-11-22; p = 4,
    // pair 2 is (1, 3))
    if (c == 21) st.map.var_i[T + 2] = 3
    if (c == 21) st.map.var_j[T + 2] = 1
    if (c == 22) st.map.var_j[T + 2] = 1
    if (c == 23) st.s = st.s[|1 \ T + 5|]
    if (c == 24) st.Gamma = st.Gamma[|1, 1 \ T + 5, T + 5|]

    // 17-20 are the rank function's guards, on a valid Z
    if (c >= 17 & c <= 20) {
        srmr = dwlsmv_resid(st, rho, rep, E, Z)
        if (c == 17) Rk = dwlsmv_resid_rank(J(3, 4, 0), 3)
        if (c == 18) Rk = dwlsmv_resid_rank(Z, 0)
        if (c == 19) Rk = dwlsmv_resid_rank(Z, .)
        if (c == 20) {
            Z[1, 2] = .
            Rk = dwlsmv_resid_rank(Z, 3)
        }
        return
    }
    srmr = dwlsmv_resid(st, rho, rep, E, Z)
}

// ===================== driver pieces ======================================================
end

* ---- A. arithmetic ------------------------------------------------------------------------
display as text _n "A. arithmetic from the fixture's own values (1e-12)"
foreach f in E_cfa1 F_efa3 R_cfa1 E3_df0 G_efa2 {
    mata: tr_a("`f'")
}

* ---- B. end to end -------------------------------------------------------------------------
display as text _n "B. end to end against lavaan, every exported fit"
mata: st_local("nfit", strofreal(rows(TR_FT)))
forvalues f = 1/`nfit' {
    mata: st_local("fid", TR_FT[`f', 1])
    mata: st_local("conv", TR_FT[`f', 6])
    if "`conv'" != "1" {
        display as text "  B `fid': not converged in lavaan, nothing exported"
        continue
    }
    mata: st_local("dat", TR_FT[`f', 2])
    mata: st_local("rng", TR_FT[`f', 3])
    mata: st_local("model", TR_FT[`f', 4])
    mata: st_local("mm", TR_FT[`f', 5])
    mata: st_local("items", tr_items("`rng'"))
    quietly import delimited using "fixtures/data/`dat'.csv", clear varnames(1) delimiters(",") stringcols(_all)
    quietly destring y*, replace
    tempvar touse
    if "`model'" == "efa" {
        local fl
        forvalues k = 1/`mm' {
            local fl `fl' F`k'
        }
        quietly dwlsmv (`fl' -> `items'), efa(`fl')
        if e(converged) != 1 {
            mata: tr_fail("B `fid': dwlsmv EFA did not converge")
            continue
        }
        quietly generate byte `touse' = e(sample)
    }
    else {
        quietly generate byte `touse' = 1
    }
    mata: tr_b("`fid'", "`items'", "`model'", `mm', "`touse'", "`dat'")
    drop `touse'
}

* ---- C. hand-built -------------------------------------------------------------------------
display as text _n "C. hand-built bound, repaired and tie cases"
mata: tr_c()

* ---- D. guards -----------------------------------------------------------------------------
display as text _n "D. every guard aborts with rc 3498 and says why (the errprintf lines below are expected)"
local nph 24
local ph1  "conditional"
local ph2  "at least two indicators"
local ph3  "rho has"
local ph4  "repaired has"
local ph5  "not a polychoric row"
local ph6  "missing polychoric"
local ph7  "missing implied correlation"
local ph8  "missing Gamma entries"
local ph9  "non-positive polychoric variance"
local ph10 "non-positive polychoric variance"
local ph11 "repaired must be 0 or 1"
local ph12 "repaired must be 0 or 1"
local ph13 "bound must be 0 or 1"
local ph14 "bound has"
local ph15 "N must be positive"
local ph16 "N must be positive"
local ph17 "not square"
local ph18 "k must be at least 1"
local ph19 "k must be at least 1"
local ph20 "not symmetric"
local ph21 "outside 1 <= i < j"
local ph22 "outside 1 <= i < j"
local ph23 "shorter than"
local ph24 "shorter than"

* the control: an unmodified valid call must run cleanly
quietly tr_run "" mata: tr_d(0)
if r(rc) {
    display as error "D0: the valid control call failed, rc " r(rc)
    mata: tr_fail("D0: the valid control call failed")
}
forvalues c = 1/`nph' {
    quietly tr_run "`ph`c''" mata: tr_d(`c')
    if r(rc) != 3498 {
        display as error "D`c': rc " r(rc) ", expected 3498"
        mata: tr_fail("D`c': rc " + strofreal(`=r(rc)') + ", expected 3498")
    }
    else if r(found) != 1 {
        display as error "D`c': rc 3498 but the message does not say: `ph`c''"
        mata: tr_fail("D`c': rc 3498 but the errprintf phrase is absent: `ph`c''")
    }
}

* ---- verdict -------------------------------------------------------------------------------
mata: st_local("nfail", strofreal(rows(TR_FAIL)))
mata: for (i = 1; i <= rows(TR_FAIL); i++) printf("{err}%s\n", TR_FAIL[i])
adopath - "`here'/src/ado"
adopath - "`here'/build"
if `nfail' {
    display as error _n "test_resid_diag: `nfail' failure(s)"
    exit 9
}
display as result "test_resid_diag: SRMR, E, Z and the ranking agree with lavaan; excluded pairs, ties and every guard behave"
exit 0
