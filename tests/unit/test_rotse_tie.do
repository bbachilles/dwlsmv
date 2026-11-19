*! Release A, unit A4 near-tie: does dwlsmv post the same rotated solution as lavaan at a tie in
*! the column-ordering statistic, and if the column order differs, is it the same solution up to
*! the signed permutation that RS5(c) carries through the covariance?
*!
*! Drafted by a read-only agent 2026-11-17; executed by the main session 2026-11-18 (6 and 4 processors).
*!
*! Design and acceptance criterion: E-109 (A4 near-tie criterion) (declared before any
*! data; implemented here as written). Inputs: fixtures/data/efatie_T.csv, efatie_B.csv
*! (92_efa_tie_data.R). Comparator: fixtures/expected/lavaan_efa_tie.csv, lavaan_efa_tie_vcov.csv
*! (99r_lavaan_efa_tie.R; lavaan 0.7-3, pinned tier-2 settings, as test_rotse_post.do).
*!
*! WHAT THIS CERTIFIES. On Case T (all eight primaries .70, both column statistics 4.5 by
*! construction) and Case B (y8 on f1 at .56, exactly the population 0.8 cutoff), dwlsmv (F1 F2 ->
*! y1-y8), efa(F1 F2) posts a rotated solution that equals lavaan's up to a signed column
*! permutation P, and the rotated e(b)/e(V) mapped through that P agree with lavaan's SEs and
*! covariance at the spec 7.2 / 7.2.1 tolerances.
*!
*! DECIDING CRITERION (one per case, the same for T and B; fixed before any data):
*!   Step A  of the 8 signed column permutations of m = 2 (2 orders x 4 sign patterns), EXACTLY ONE
*!           must put every rotated estimate (loadings, factor correlation, thresholds) within
*!           1e-5 abs + 1e-4 rel of lavaan's, and every other must miss. Else FAIL.
*!   Step B  with e(b), e(V) mapped through that P (RS5(c): V -> Pi V Pi'):
*!             SE          |a - b| <= 1e-8 + 1e-5 |b|
*!             covariance  |a - b| <= 1e-8 + 1e-5 sqrt(V_L,ii V_L,jj), lower triangle, V_L lavaan's
*!           (spec 7.2.1's attribution rule is NOT applied: a miss is a miss.)
*!   Outcome per case, one greppable line "TIEOUTCOME case X: ...":
*!             AGREE          Steps A and B pass and P is the identity
*!             ORDER DIFFERS  Steps A and B pass and P swaps the columns, both signs +
*!             FAIL           Step A or B fails, or P has a negative sign (a sign is a rule
*!                            mismatch, positive column sums, not a tie effect)
*!   The test passes iff neither case is FAIL. It does NOT assert AGREE or ORDER DIFFERS: which of
*!   them occurs is the finding (A3.7's revisit trigger), not a tolerance.
*!
*! DIAGNOSTICS (printed, never deciding, never counted):
*!   - dwlsmv's per-column index statistic from its posted e(Lambda) (mean row index of
*!     |l| >= 0.8 max|l| in the column), and whether both are 4.5 (the Case T precondition; if
*!     not, "tie not reached" is PRINTED and reported, it is not a failure, and nothing is reseeded)
*!   - |l| / (0.8 max) of every loading within 10% of the cutoff (the y8 ratio on Case B)
*!   - lavaan's "order" and "tie" rows and its criterion value Q, from the fixture
*!   - e(rot_cond) and the minima information dwlsmv posts: e(rot_starts), e(rot_converged),
*!     e(rot_minima), e(rot_atbest), e(rot_Q), e(rot_gpnorm), e(rot_gap)
*!   - the worst estimate ratio of each of the 8 signed permutations (the Step A table)
*!
*! WHAT IT DOES NOT CERTIFY. Ties at m >= 3 or among more than two columns; Mplus's order (tier 3
*! is not built); calibration (A6); any tie rule (A3.7 N1 is unchanged unless the owner changes
*! it); a near tie within estimation error (about 1e-8) of the 0.8 cutoff.
*!
*! DETERMINISM. Run at the default processor count and at set processors 4 (a project rule: wrapper
*! do-file in %TEMP%). The lines "TPDET ..." print trace(V), the largest SE and one element to
*! 15 digits per case; the lines "TIEOUTCOME" should be identical between the two logs (if dwlsmv's
*! order differs between processor counts, that is a finding). Nothing is compared across runs
*! here.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 120

* ---- a missing fixture is a deliberate skip (exit 8), never a pass -----------------------
foreach f in fixtures/expected/lavaan_efa_tie.csv fixtures/expected/lavaan_efa_tie_vcov.csv ///
             fixtures/data/efatie_T.csv fixtures/data/efatie_B.csv {
    capture confirm file "`f'"
    if _rc {
        display as text "test_rotse_tie: SKIP, fixture `f' not found (generators 92_ and 99r_ not yet run)"
        exit 8
    }
}

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

display "TPENV stata_version " c(stata_version) " processors " c(processors)

local nfail 0

* ---- the tier-2 fixture, as text, converted with strtoreal() (the float-import trap) -----------------
import delimited using "fixtures/expected/lavaan_efa_tie.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: XC = st_sdata(., "case")
mata: XQ = st_sdata(., "quantity")
mata: XL = st_sdata(., "label")
mata: XV = strtoreal(st_sdata(., "value"))
import delimited using "fixtures/expected/lavaan_efa_tie_vcov.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
mata: WC = st_sdata(., "case")
mata: WR = st_sdata(., "row")
mata: WK = st_sdata(., "col")
mata: WV = strtoreal(st_sdata(., "value"))

mata:

// the Stata column of e(b) for lavaan's label (as test_rotse_post.do): fK=~yJ -> yJ : FK,
// yJ|tK -> /yJ : cutK, fA~~fB -> /: cov(FA,FB) or cov(FB,FA). Missing if absent.
real scalar tpt_col(string matrix S, string scalar lab)
{
    string scalar l, r, a, b, f
    real colvector k

    if (strpos(lab, "=~") > 0) {
        l = substr(lab, 1, strpos(lab, "=~") - 1)
        r = substr(lab, strpos(lab, "=~") + 2, .)
        f = strupper(substr(l, 1, 1)) + substr(l, 2, .)
        k = selectindex((S[., 1] :== r) :& (S[., 2] :== f))
    }
    else if (strpos(lab, "|") > 0) {
        l = substr(lab, 1, strpos(lab, "|") - 1)
        r = substr(lab, strpos(lab, "|") + 1, .)
        k = selectindex((S[., 1] :== "/" + l) :& (S[., 2] :== "cut" + substr(r, 2, .)))
    }
    else if (strpos(lab, "~~") > 0) {
        l = substr(lab, 1, strpos(lab, "~~") - 1)
        r = substr(lab, strpos(lab, "~~") + 2, .)
        a = strupper(substr(l, 1, 1)) + substr(l, 2, .)
        b = strupper(substr(r, 1, 1)) + substr(r, 2, .)
        k = selectindex((S[., 1] :== "/") :& ((S[., 2] :== "cov(" + a + "," + b + ")") :|
                                              (S[., 2] :== "cov(" + b + "," + a + ")")))
    }
    else {
        return(.)
    }
    return(length(k) == 1 ? k : .)
}

// lavaan's label under the permutation q (lavaan column k is dwlsmv column q[k]): the label of
// the dwlsmv parameter that plays its role. Only loadings move; thresholds and the correlation
// keep their label (m = 2 has one correlation).
string scalar tpt_srclabel(string scalar lab, real rowvector q)
{
    string scalar l, r
    real scalar k

    if (strpos(lab, "=~") > 0) {
        l = substr(lab, 1, strpos(lab, "=~") - 1)
        r = substr(lab, strpos(lab, "=~") + 2, .)
        k = strtoreal(substr(l, 2, .))
        return("f" + strofreal(q[k]) + "=~" + r)
    }
    return(lab)
}

// the sign with which dwlsmv's source parameter maps to lavaan's: a loading in lavaan column k
// carries s[k]; the factor correlation carries s[1]*s[2] (one-column reflection flips it);
// thresholds are untouched.
real scalar tpt_sign(string scalar lab, real rowvector s)
{
    real scalar k

    if (strpos(lab, "=~") > 0) {
        k = strtoreal(substr(lab, 2, strpos(lab, "=~") - 2))
        return(s[k])
    }
    if (strpos(lab, "~~") > 0) {
        return(s[1] * s[2])
    }
    return(1)
}

// Steps A and B for one case against the posted e(b), e(V). Prints the Step A table, the
// outcome line, and returns the outcome code in out[1]: 0 AGREE, 1 ORDER DIFFERS, 2 FAIL.
real rowvector tpt_case(string scalar cs, real scalar hasV)
{
    external string colvector XC, XQ, XL, WC, WR, WK
    external real colvector XV, WV
    string matrix S
    string colvector lab
    string scalar reason, name
    real rowvector b, qv, sv, out
    real matrix V, PERM, SGN, SRC, SG
    real colvector ix, ixs, est, sev, vd, worst, unm, pass, src, sg, w
    real scalar n, ip, i, c, x, r, np, win, nse, ncov, wse, wcov, tol, unmB, a, d
    real colvector ni, nj, jj
    string scalar wlse, wlcov

    b = st_matrix("e(b)")
    S = st_matrixcolstripe("e(b)")
    ix  = selectindex((XC :== cs) :& (XQ :== "est"))
    lab = XL[ix]
    est = XV[ix]
    n   = length(ix)

    PERM = (1, 2 \ 1, 2 \ 1, 2 \ 1, 2 \ 2, 1 \ 2, 1 \ 2, 1 \ 2, 1)
    SGN  = (1, 1 \ 1, -1 \ -1, 1 \ -1, -1 \ 1, 1 \ 1, -1 \ -1, 1 \ -1, -1)
    SRC  = J(n, 8, .)
    SG   = J(n, 8, .)
    worst = J(8, 1, 0)
    unm   = J(8, 1, 0)

    // ---- Step A: all 8 signed permutations --------------------------------------------
    for (ip = 1; ip <= 8; ip++) {
        qv = PERM[ip, .]
        sv = SGN[ip, .]
        for (i = 1; i <= n; i++) {
            c = tpt_col(S, tpt_srclabel(lab[i], qv))
            if (c >= .) {
                unm[ip] = unm[ip] + 1
                continue
            }
            SRC[i, ip] = c
            SG[i, ip]  = tpt_sign(lab[i], sv)
            x = SG[i, ip] * b[c]
            r = abs(x - est[i]) / (1e-5 + 1e-4 * abs(est[i]))
            if (r > worst[ip]) {
                worst[ip] = r
            }
        }
        if (unm[ip] > 0) {
            worst[ip] = .
        }
        // Mata's printf has no %+g; the signs are printed as strings
        printf("  STEPA case %s: order (%g,%g) signs (%s,%s): worst estimate ratio %9.3e, unmatched %g\n",
               cs, qv[1], qv[2], (sv[1] > 0 ? "+1" : "-1"), (sv[2] > 0 ? "+1" : "-1"), worst[ip], unm[ip])
    }
    pass = (worst :<= 1)
    np   = sum(pass)
    win  = 0
    if (np == 1) {
        win = selectindex(pass)
    }
    out = (2, np, win, ., ., ., ., ., ., .)

    if (np != 1) {
        reason = "Step A: " + strofreal(np) + " of 8 signed permutations within the theta-hat tier (exactly 1 required)"
        printf("TIEOUTCOME case %s: FAIL, %s\n", cs, reason)
        return(out)
    }

    qv = PERM[win, .]
    sv = SGN[win, .]
    printf("  STEPA case %s: unique aligning P = order (%g,%g), signs (%s,%s)\n", cs, qv[1], qv[2],
           (sv[1] > 0 ? "+1" : "-1"), (sv[2] > 0 ? "+1" : "-1"))

    // ---- Step B: SEs and covariance through that P --------------------------------------
    if (!hasV) {
        printf("TIEOUTCOME case %s: FAIL, Step B: no e(V) posted (e(rot_nose) set)\n", cs)
        return(out)
    }
    V   = st_matrix("e(V)")
    src = SRC[., win]
    sg  = SG[., win]

    // lavaan's SEs and variances, indexed by position in the est label list
    sev = J(n, 1, .)
    ixs = selectindex((XC :== cs) :& (XQ :== "se"))
    for (i = 1; i <= length(ixs); i++) {
        jj = selectindex(lab :== XL[ixs[i]])
        if (length(jj) == 1) {
            sev[jj] = XV[ixs[i]]
        }
    }
    vd = J(n, 1, .)
    w  = selectindex((WC :== cs) :& (WR :== WK))
    for (i = 1; i <= length(w); i++) {
        jj = selectindex(lab :== WR[w[i]])
        if (length(jj) == 1) {
            vd[jj] = WV[w[i]]
        }
    }

    unmB = 0
    nse = 0
    wse = 0
    wlse = ""
    for (i = 1; i <= n; i++) {
        if (sev[i] >= . | vd[i] >= .) {
            unmB = unmB + 1
        }
        if (sev[i] >= .) {
            continue
        }
        x = sqrt(V[src[i], src[i]])
        r = abs(x - sev[i]) / (1e-8 + 1e-5 * abs(sev[i]))
        nse = nse + 1
        if (r > wse | r >= .) {
            wse  = r
            wlse = lab[i]
        }
    }
    ncov = 0
    wcov = 0
    wlcov = ""
    w = selectindex(WC :== cs)
    for (i = 1; i <= length(w); i++) {
        ni = selectindex(lab :== WR[w[i]])
        nj = selectindex(lab :== WK[w[i]])
        if (length(ni) != 1 | length(nj) != 1) {
            unmB = unmB + 1
            continue
        }
        a = ni[1]
        d = nj[1]
        if (vd[a] >= . | vd[d] >= .) {
            continue
        }
        // RS5(c): V_aligned = Pi V Pi'
        x   = sg[a] * sg[d] * V[src[a], src[d]]
        tol = 1e-8 + 1e-5 * sqrt(vd[a] * vd[d])
        r   = abs(x - WV[w[i]]) / tol
        ncov = ncov + 1
        if (r > wcov | r >= .) {
            wcov  = r
            wlcov = WR[w[i]] + " , " + WK[w[i]]
        }
    }
    printf("  STEPB case %s: SEs %g, worst ratio %9.2e (%s); covariances %g (expected %g), worst %9.2e (%s); unmatched %g\n",
           cs, nse, wse, wlse, ncov, n * (n + 1) / 2, wcov, wlcov, unmB)
    out[4..10] = (nse, wse, ncov, wcov, unmB, n, .)

    reason = ""
    if (unmB != 0 | nse != n | ncov != n * (n + 1) / 2) {
        reason = "Step B: fixture rows not all matched"
    }
    else if (!(wse <= 1)) {
        reason = "Step B: a standard error misses its tolerance, worst ratio " + strofreal(wse, "%9.2e")
    }
    else if (!(wcov <= 1)) {
        reason = "Step B: a covariance misses its tolerance, worst ratio " + strofreal(wcov, "%9.2e")
    }
    else if (sv[1] < 0 | sv[2] < 0) {
        reason = "P has a negative sign (a rule mismatch, positive column sums, not a tie effect)"
    }
    if (reason != "") {
        printf("TIEOUTCOME case %s: FAIL, %s\n", cs, reason)
        return(out)
    }
    if (qv[1] == 1) {
        out[1] = 0
        printf("TIEOUTCOME case %s: AGREE (P is the identity; Steps A and B pass)\n", cs)
    }
    else {
        out[1] = 1
        printf("TIEOUTCOME case %s: ORDER DIFFERS (P swaps the columns, both signs +; Steps A and B pass; A3.7 revisit trigger)\n", cs)
    }
    return(out)
}

// diagnostics: never deciding. dwlsmv's index statistic from e(Lambda); lavaan's rows.
void tpt_diag(string scalar cs)
{
    external string colvector XC, XQ, XL
    external real colvector XV
    real matrix L
    real colvector ix, hi
    real rowvector idx, cmax
    real scalar k, j, p, m, both, rt
    string scalar nm

    L = st_matrix("e(Lambda)")
    p = rows(L)
    m = cols(L)
    idx  = J(1, m, .)
    cmax = J(1, m, .)
    for (k = 1; k <= m; k++) {
        cmax[k] = max(abs(L[., k]))
        hi = selectindex(abs(L[., k]) :>= 0.8 * cmax[k])
        idx[k] = mean(hi)
    }
    both = (max(abs(idx :- 4.5)) < 1e-12)
    printf("  DIAG case %s: dwlsmv index statistic (%9.6f, %9.6f); both 4.5: %g\n", cs, idx[1], idx[2], both)
    if (cs == "T" & !both) {
        printf("  DIAG case T: TIE NOT REACHED in dwlsmv's posted solution (reported, not a failure, not reseeded)\n")
    }
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            rt = abs(L[j, k]) / (0.8 * cmax[k])
            if (rt >= 0.9 & rt <= 1.1) {
                printf("  DIAG case %s: dwlsmv |l|/(0.8 max) for y%g on F%g = %9.6f (loading %9.6f)\n", cs, j, k, rt, L[j, k])
            }
        }
    }
    ix = selectindex((XC :== cs) :& ((XQ :== "order") :| (XQ :== "tie")))
    for (k = 1; k <= length(ix); k++) {
        printf("  DIAG case %s: lavaan %s %s = %15.10g\n", cs, XQ[ix[k]], XL[ix[k]], XV[ix[k]])
    }
    ix = selectindex((XC :== cs) :& (XQ :== "rot") :& (XL :== "Q"))
    if (length(ix) == 1) {
        printf("  DIAG case %s: lavaan criterion Q = %15.10g; dwlsmv e(rot_Q) = %15.10g\n", cs, XV[ix], st_numscalar("e(rot_Q)"))
    }
}

// the determinism line, diffed between the 6- and 4-processor logs by the main session
void tpt_det(string scalar cs)
{
    real matrix V

    V = st_matrix("e(V)")
    printf("TPDET %s trace %22.15e maxdiag %22.15e V21 %22.15e Vlast %22.15e\n", cs, trace(V),
           max(diagonal(V)), V[2, 1], V[rows(V), cols(V) - 1])
}
end

* ---- per case ----------------------------------------------------------------------------
foreach c in T B {
    import delimited using "fixtures/data/efatie_`c'.csv", clear varnames(1) delimiters(",")
    capture noisily dwlsmv (F1 F2 -> y1-y8), efa(F1 F2)
    if _rc {
        display as error "TIEOUTCOME case `c': FAIL, dwlsmv rc " _rc
        local ++nfail
        continue
    }
    capture confirm matrix e(V)
    local hasv = (_rc == 0)
    display as text "  case `c': n* " colsof(e(b)) ", e(rank) " e(rank) ", e(rot_cond) " %9.3e e(rot_cond) ///
        ", e(rot_nose) '`e(rot_nose)''"
    display as text "  case `c': rotation starts " e(rot_starts) ", converged " e(rot_converged) ///
        ", distinct minima " e(rot_minima) ", starts at the reported minimum " e(rot_atbest) ///
        ", e(rot_Q) " %12.8f e(rot_Q) ", e(rot_gpnorm) " %9.2e e(rot_gpnorm) ", e(rot_gap) " %9.2e e(rot_gap)

    mata: tpt_diag("`c'")
    mata: TPO = tpt_case("`c'", `hasv')
    mata: st_local("code", strofreal(TPO[1], "%21.0g"))
    if `hasv' {
        mata: tpt_det("`c'")
    }
    if `code' == 2 {
        display as error "test_rotse_tie: case `c' is FAIL"
        local ++nfail
    }
}

if `nfail' {
    display as error _n "test_rotse_tie: `nfail' failing case(s)"
    exit 9
}
display as result "test_rotse_tie: neither case is FAIL (see the TIEOUTCOME lines for AGREE or ORDER DIFFERS, which this test does not assert)"
exit 0
