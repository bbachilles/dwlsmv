*! Conditional-covariate Stage 1 (one binary x, one group): s, Gamma in the
*! Hessian form, and the blocks H_j, U_rho, U_theta, against the analytic R
*! reference fixtures/expected/condx_ref_{A,D}_*.csv (99n_condx_ref.R;
*! MANIFEST_condx_ref.json), with the refusals and the unconditional path.
*!
*! Executed 2026-11-10 (suite green). Drafted by a read-only agent together
*! with the Mata it tests; run, and extended with refusal cases 10 and 11, by
*! the main session. Source mutants M1, M2, M4 and two on the P1 rule were
*! run by hand that day and all failed this file.
*!
*! THIS IS A DEVELOPMENT TEST. It is not link 2 of the conditional-mode
*! validation (E-098, FROZEN 2026-11-10).
*! The recorded link 2 run needs its own do-file,
*! tests/unit/test_condx_link2.do, covering the whole table of prereg section 3
*! (s, Gamma, theta-hat, SEs, T3, a, b, df), whose SHA-256 goes into a dated
*! prereg addendum BEFORE the run that is recorded as the gate result ([V1]).
*! Passing this file is not that result.
*!
*! Section 8 (2026-11-12): the refusals of the user's package of that day,
*! each at rc 3498 with the name the decision requires (item; pair and cell;
*! pair and group), and dwlsmv_gammadiag()'s refusal of a conditional result.
*! Several are reached only by direct calls with crafted arguments; three item
*! refusals and pair statuses other than 2 have no input and no test
*! (E-104).
*!
*! Checks, per dataset (condx_A, condx_D; K = (3, 3, 2); N = 1500 each):
*!   1. s, every element, at 1e-8 abs + 1e-8 rel (prereg section 3: the
*!      threshold tier, which the slopes take by [V5]); Gamma, every entry of
*!      the symmetric matrix, at 1e-8 abs + 1e-8 rel ([V6] (i)). Matched by
*!      LABEL, never by position (spec 4.2); the reference's "y<j>~x" is the
*!      map's "y<j>~x1".
*!   2. Exact structure, no tolerance: the phi2 copy guard (the four-corner
*!      difference of dwlsmv_cx_phi2grid() equals dwlsmv_dcellprob_rho() bit
*!      for bit, every pair and group: the prereg's phi2 clause requires one
*!      density form); U_theta's slope columns are exactly 0 in every
*!      observed x = 0 cell (@sec-condx-hess); Gamma exactly symmetric; bound
*!      all zero; Gamma from dwlsmv_cx_psi_matrix() at ctl 0 equals
*!      st.Gamma exactly.
*!   2b. The blocks H_j, U_rho, U_theta against the reference: DIAGNOSTIC
*!      ONLY, printed, never failing (prereg section 3a item 2 makes them
*!      diagnostics, and prereg section 3 sets no tier for them; this file
*!      does not invent one).
*!   5. Negative controls, in-file: Gamma rebuilt with (ctl 1) the slope part
*!      of U_theta dropped, (2) all of U_theta dropped, (3) H_j's (tau, b)
*!      cross entries dropped, (4) the quotient term -u dP/P dropped. Each
*!      MUST FAIL the Gamma tier (gap ratio > 1); a control that passes
*!      means the comparison cannot see that term, and fails this test.
*! Once:
*!   3. The unconditional path is untouched: dwlsmv_statmap_build(K) and
*!      dwlsmv_statmap_build(K, 0) agree field by field (==) on case_D's K;
*!      pos_rho equals the pre-change closed form T + (i-1)p - i(i-1)/2 +
*!      (j-i) on every pair; pos_inv round-trips every position; q = 0. The
*!      conditional map's labels, slope positions and round trip are pinned.
*!      dwls_sample_stats.mata itself is not edited by this build; its
*!      bitwise identity is checked by the scratch tool snapshot_v1.do (save
*!      before installing, compare after; not part of the suite), and
*!      by the existing suite (test_gamma.do, test_sample_stats.do,
*!      test_statmap.do).
*!   6. The item fit's stopping rule (@sec-condx-itemstop): the stop helper
*!      on small, large, signed and missing steps and scores; the fit takes
*!      2 to 49 iterations; through the tolc hook, the iterate-50 branch is
*!      reached (it == 50) and its point is within 1.1e-12 of production's;
*!      a tolerance nothing can meet is refused at rc 3498.
*!   7. Targeted (user, 2026-11-11): dwlsmv_cx_item_ok() on ordered, equal,
*!      reversed and missing thresholds; on three constructed tables (chosen
*!      after measuring) every returned point is strictly ordered; omitting
*!      tolc equals (1e-12, 1e-10) exactly and a looser d1 stops earlier,
*!      which pins the two production constants; a malformed tolc is
*!      refused at rc 3498.
*!   4. Refusals, each rc == 3498: x coded 1/2; x coded -1/+1; constant x;
*!      a missing x; X with two columns; X of the wrong length; x with three
*!      values; an item separated by x (y3 := x, Result P1); a duplicated item
*!      (y2 := y1: a boundary residual rho, P2, no repair); pos_slope on an
*!      unconditional map. The messages of four of them are read back from a
*!      log.
*!
*! MUTANTS for the main session to run by hand (edit, run this file, revert):
*!   M1 dwlsmv_cx_rf_eval(): return group 0 alone (w0 = 1, w1 = 0). Check 1
*!      (s, the residual correlations) must fail.
*!   M2 dwlsmv_cx_pair_tab(): shift the bounds by + thi[Ki] * xg instead of
*!      - (both a and b). Check 1 (Gamma) must fail.
*!   M3 dwlsmv_cx_uniP(): drop the upper-tail branch. EXPECTED TO SURVIVE on
*!      these data (prereg [V4]: the difference is far below the tier at
*!      this domain's probabilities); record it as a known survivor, not a
*!      defect of the test.
*!   M4 dwlsmv_cx_check_x(): accept two values other than 0/1. Check 4 must
*!      fail (cases 1 and 2).
*!
*! Reference CSVs are read as text with an explicit delimiters(",") and
*! converted with strtoreal() (a project rule: float columns and the "|" in labels).
*! Every reference vector is checked for missing values explicitly, because
*! Mata's max() skips them and dwlsmv_gapratio() would then pass vacuously.
*!
*! Return codes: 0 pass, 8 skip (the conditional Stage 1 is not installed),
*! 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

capture confirm file "src/mata/dwls_sample_stats_cx.mata"
if _rc {
    display as text "test_condx_stage1: src/mata/dwls_sample_stats_cx.mata is absent -- skip"
    exit 8
}

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_statmap_q.mata         ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
             src/mata/dwlsmv_pos_slope.mata         ///
             src/mata/dwlsmv_pos_inv.mata           ///
             src/mata/dwlsmv_catcounts.mata         ///
             src/mata/dwlsmv_catdegen.mata          ///
             src/mata/dwlsmv_cumprop.mata           ///
             src/mata/dwlsmv_thresholds.mata        ///
             src/mata/dwlsmv_psi_tau_table.mata     ///
             src/mata/dwlsmv_cellprob.mata          ///
             src/mata/dwlsmv_dcellprob_rho.mata     ///
             src/mata/dwlsmv_wmat.mata              ///
             src/mata/dwlsmv_crosstab.mata          ///
             src/mata/dwlsmv_pairinfo.mata          ///
             src/mata/dwls_polychoric.mata          ///
             src/mata/dwlsmv_koff.mata              ///
             src/mata/dwlsmv_psi_tau_stack.mata     ///
             src/mata/dwlsmv_psi_rho_table.mata     ///
             src/mata/dwlsmv_psi_matrix.mata        ///
             src/mata/dwlsmv_psibar.mata            ///
             src/mata/dwls_gamma.mata               ///
             src/mata/dwlsmv_emptycell1.mata        ///
             src/mata/dwlsmv_emptycell_repair.mata  ///
             src/mata/dwlsmv_emptycell_pairs.mata   ///
             src/mata/dwls_sample_stats.mata        ///
             src/mata/dwlsmv_cx_check.mata          ///
             src/mata/dwlsmv_cx_item_fit.mata           ///
             src/mata/dwlsmv_cx_pair_fit.mata           ///
             src/mata/dwlsmv_cx_psi_matrix.mata     ///
             src/mata/dwlsmv_cx_psibar.mata         ///
             src/mata/dwls_gamma_cx.mata            ///
             src/mata/dwls_sample_stats_cx.mata     ///
             src/mata/dwlsmv_asym.mata              ///
             src/mata/dwlsmv_npatterns.mata         ///
             src/mata/dwlsmv_thralias.mata          ///
             src/mata/dwlsmv_dalias.mata            ///
             src/mata/dwlsmv_gammadiag.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_condx_stage1: `f' failed to load"
        exit 9
    }
}

* ---- fixtures --------------------------------------------------------------
foreach t in A D {
    import delimited using "fixtures/expected/condx_ref_`t'_s.csv", ///
        clear varnames(1) delimiters(",") stringcols(_all)
    mata: SL_`t' = st_sdata(., "label")
    mata: SV_`t' = strtoreal(st_sdata(., "value"))

    import delimited using "fixtures/expected/condx_ref_`t'_gamma.csv", ///
        clear varnames(1) delimiters(",") stringcols(_all)
    mata: GR_`t' = st_sdata(., "row")
    mata: GC_`t' = st_sdata(., "col")
    mata: GV_`t' = strtoreal(st_sdata(., "value"))

    import delimited using "fixtures/expected/condx_ref_`t'_blocks.csv", ///
        clear varnames(1) delimiters(",") stringcols(_all)
    mata: BB_`t' = st_sdata(., "block")
    mata: BR_`t' = st_sdata(., "row")
    mata: BC_`t' = st_sdata(., "col")
    mata: BV_`t' = strtoreal(st_sdata(., "value"))

    import delimited using "fixtures/data/condx_`t'.csv", ///
        clear varnames(1) delimiters(",") stringcols(_all)
    mata: Y_`t' = strtoreal(st_sdata(., ("y1", "y2", "y3")))
    mata: X_`t' = strtoreal(st_sdata(., "x"))
}

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

// The reference's slope label "y<j>~x" is the map's "y<j>~x1".
string scalar dwlsmv_t_cx_lab(string scalar L)
{
    if (strlen(L) >= 2) {
        if (substr(L, strlen(L) - 1, 2) == "~x") {
            return(L + "1")
        }
    }
    return(L)
}

// Position of a reference label in the map; 0 if absent, -1 if ambiguous.
real scalar dwlsmv_t_cx_pos(struct dwlsmv_statmap scalar m, string scalar L)
{
    string scalar M
    real scalar k, hit

    M   = dwlsmv_t_cx_lab(L)
    hit = 0
    for (k = 1; k <= m.D; k++) {
        if (m.label[k] == M) {
            if (hit) {
                return(-1)
            }
            hit = k
        }
    }
    return(hit)
}

// A theta label "y<j>|t<c>" or "y<j>~x" -> (j, index within theta_j).
real rowvector dwlsmv_t_cx_th(string scalar L, real colvector K)
{
    real scalar p, j

    p = strpos(L, "|t")
    if (p) {
        j = strtoreal(substr(L, 2, p - 2))
        return((j, strtoreal(substr(L, p + 2, .))))
    }
    p = strpos(L, "~x")
    if (p) {
        j = strtoreal(substr(L, 2, p - 2))
        return((j, K[j]))
    }
    return((., .))
}

// theta_j (K_j x 1) out of s.
real colvector dwlsmv_t_cx_theta(struct dwlsmv_statmap scalar m, real colvector s,
                                 real scalar j)
{
    real colvector th
    real scalar c

    th = J(m.K[j], 1, .)
    for (c = 1; c <= m.C[j]; c++) {
        th[c] = s[dwlsmv_pos_tau(m, j, c)]
    }
    th[m.K[j]] = s[dwlsmv_pos_slope(m, j, 1)]
    return(th)
}

real scalar dwlsmv_t_cx_one(real matrix Y, real colvector X,
                            string colvector SL, real colvector SV,
                            string colvector GR, string colvector GC,
                            real colvector GV, string colvector BB,
                            string colvector BR, string colvector BC,
                            real colvector BV, string scalar tag)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_statmap scalar m
    string colvector FAIL
    real colvector K, R, n, P, gbar, thi, thj, th, a, b, av, bv, Urv
    real matrix    Gr, Gx, Psi, G, H, Hst, tab, F, D4, Prm, Utm
    real rowvector Ut, ti, tc, hoff
    real scalar    k, i, j, pa, pb, rs, rg, rr, ctl, ll, Ur, xg, r, Ki, Kj
    real scalar    pk, D, nexp, blk, gabs, grel, nmiss

    K    = (3 \ 3 \ 2)
    FAIL = J(0, 1, "")

    st = dwls_sample_stats_cx(Y, K, X)
    m  = st.map
    D  = m.D

    // ---- structure ---------------------------------------------------------
    if (st.N != rows(Y)) {
        FAIL = FAIL \ (tag + ": N is not rows(Y)")
    }
    if (D != 11 | dwlsmv_statmap_q(m) != 1) {
        FAIL = FAIL \ (tag + ": the map is not T + p + P = 11 with q = 1")
    }
    if (rows(st.bound) != m.P | cols(st.bound) != 1) {
        FAIL = FAIL \ (tag + ": bound is not P x 1")
    }
    else if (any(st.bound :!= 0)) {
        FAIL = FAIL \ (tag + ": bound is not all zero")
    }
    if (missing(st.s) | missing(st.Gamma) | missing(st.psibar)) {
        FAIL = FAIL \ (tag + ": s, Gamma or psibar has missing values")
        for (k = 1; k <= rows(FAIL); k++) {
            printf("{err}%s\n", FAIL[k])
        }
        return(rows(FAIL))
    }
    if (any(st.Gamma :!= st.Gamma')) {
        FAIL = FAIL \ (tag + ": Gamma is not exactly symmetric")
    }
    printf("%s: max |psibar| %9.2e (diagnostic, never asserted)\n", tag, max(abs(st.psibar)))

    // ---- 1a. s, by label ---------------------------------------------------
    nmiss = missing(SV)
    if (rows(SL) != D | nmiss) {
        FAIL = FAIL \ (tag + ": the s fixture is not D complete rows")
    }
    else {
        R = J(D, 1, .)
        for (k = 1; k <= D; k++) {
            pa = dwlsmv_t_cx_pos(m, SL[k])
            if (pa <= 0) {
                FAIL = FAIL \ (tag + ": s label not in the map: " + SL[k])
                continue
            }
            if (!missing(R[pa])) {
                FAIL = FAIL \ (tag + ": s label repeated: " + SL[k])
            }
            R[pa] = SV[k]
        }
        if (missing(R)) {
            FAIL = FAIL \ (tag + ": the s fixture does not cover every position")
        }
        else {
            rs = dwlsmv_gapratio(st.s, R, 1e-8, 1e-8)
            printf("%s: s      gap ratio %9.3e (tier 1e-8 abs + 1e-8 rel)\n", tag, rs)
            if (missing(rs) | rs > 1) {
                FAIL = FAIL \ (tag + ": s outside its tier, ratio " + strofreal(rs, "%9.3e"))
            }
        }
    }

    // ---- 1b. Gamma, every entry, by label ----------------------------------
    Gr   = J(D, D, .)
    nexp = D * (D + 1) / 2
    if (rows(GR) != nexp | missing(GV)) {
        FAIL = FAIL \ (tag + ": the Gamma fixture is not D(D+1)/2 complete rows")
    }
    else {
        for (k = 1; k <= nexp; k++) {
            pa = dwlsmv_t_cx_pos(m, GR[k])
            pb = dwlsmv_t_cx_pos(m, GC[k])
            if (pa <= 0 | pb <= 0) {
                FAIL = FAIL \ (tag + ": Gamma label not in the map: " + GR[k] + ", " + GC[k])
                continue
            }
            Gr[pa, pb] = GV[k]
            Gr[pb, pa] = GV[k]
        }
    }
    if (missing(Gr)) {
        FAIL = FAIL \ (tag + ": the Gamma fixture does not fill the D x D matrix")
    }
    else {
        rg = dwlsmv_gapratio(st.Gamma, Gr, 1e-8, 1e-8)
        printf("%s: Gamma  gap ratio %9.3e (tier 1e-8 abs + 1e-8 rel, every entry)\n", tag, rg)
        if (missing(rg) | rg > 1) {
            FAIL = FAIL \ (tag + ": Gamma outside its tier, ratio " + strofreal(rg, "%9.3e"))
        }
    }

    // ---- 2. exact structure ------------------------------------------------
    Psi = dwlsmv_cx_psi_matrix(Y, X, m, st.s, 0)
    if (any((quadcross(Psi, Psi) :/ rows(Y)) :!= st.Gamma)) {
        FAIL = FAIL \ (tag + ": psi_matrix at ctl 0 does not reproduce st.Gamma exactly")
    }

    // The blocks at s-hat, kept for 2b.
    hoff = J(1, m.p, 0)
    for (j = 2; j <= m.p; j++) {
        hoff[j] = hoff[j - 1] + m.K[j - 1]
    }
    Hst = J(hoff[m.p] + m.K[m.p], max(K), .)
    for (j = 1; j <= m.p; j++) {
        th = dwlsmv_t_cx_theta(m, st.s, j)
        n  = dwlsmv_cx_item_counts(Y[., j], X, m.K[j])
        dwlsmv_cx_item_eval(th, m.K[j], n, 0, P, G, H, gbar, ll)
        Hst[|hoff[j] + 1, 1 \ hoff[j] + m.K[j], m.K[j]|] = H
    }
    Urv = J(m.P, 1, .)
    Utm = J(m.P, 2 * max(K), .)
    pk  = 0
    for (i = 1; i <= m.p - 1; i++) {
        for (j = i + 1; j <= m.p; j++) {
            pk  = pk + 1
            Ki  = m.K[i]
            Kj  = m.K[j]
            thi = dwlsmv_t_cx_theta(m, st.s, i)
            thj = dwlsmv_t_cx_theta(m, st.s, j)
            dwlsmv_cx_pair_blocks(Y[., i], Y[., j], X, thi, thj,
                                  st.s[dwlsmv_pos_rho(m, i, j)], 0, Ur, Ut, tab, n)
            Urv[pk] = Ur
            Utm[|pk, 1 \ pk, Ki + Kj|] = Ut

            // phi2 copy guard and the x = 0 slope zeros
            for (xg = 0; xg <= 1; xg++) {
                a   = thi[|1 \ Ki - 1|] :- thi[Ki] * xg
                b   = thj[|1 \ Kj - 1|] :- thj[Kj] * xg
                F   = dwlsmv_cx_phi2grid(a, b, st.s[dwlsmv_pos_rho(m, i, j)])
                D4  = F[|2, 2 \ Ki + 1, Kj + 1|] :- F[|1, 2 \ Ki, Kj + 1|] :-
                      F[|2, 1 \ Ki + 1, Kj|] :+ F[|1, 1 \ Ki, Kj|]
                Prm = dwlsmv_dcellprob_rho(a, b, st.s[dwlsmv_pos_rho(m, i, j)])
                if (any(D4 :!= Prm)) {
                    FAIL = FAIL \ (tag + ": phi2 copy guard: grid and dwlsmv_dcellprob_rho differ, y" +
                                   strofreal(i) + "~~y" + strofreal(j) + ", x = " + strofreal(xg))
                }
            }
            for (r = 1; r <= Ki * Kj; r++) {
                if (n[r] == 0) {
                    continue
                }
                if (tab[r, 3 + Ki] != 0 | tab[r, 3 + Ki + Kj] != 0) {
                    FAIL = FAIL \ (tag + ": a slope column of du/dtheta is not exactly 0 at x = 0, y" +
                                   strofreal(i) + "~~y" + strofreal(j))
                    break
                }
            }
        }
    }

    // ---- 2b. blocks against the reference: DIAGNOSTIC ONLY -----------------
    if (missing(BV)) {
        printf("{txt}%s: blocks fixture has missing values; block diagnostic skipped\n", tag)
    }
    else {
        av = J(rows(BV), 1, .)
        for (k = 1; k <= rows(BV); k++) {
            if (substr(BB[k], 1, 3) == "H_y") {
                j  = strtoreal(substr(BB[k], 4, .))
                ti = dwlsmv_t_cx_th(BR[k], K)
                tc = dwlsmv_t_cx_th(BC[k], K)
                if (ti[1] == j & tc[1] == j) {
                    av[k] = Hst[hoff[j] + ti[2], tc[2]]
                }
            }
            else if (substr(BB[k], 1, 5) == "Urho_") {
                pa = strpos(BR[k], "~~")
                i  = strtoreal(substr(BR[k], 2, pa - 2))
                j  = strtoreal(substr(BR[k], pa + 3, .))
                pk = (i - 1) * m.p - i * (i - 1) / 2 + (j - i)
                av[k] = Urv[pk]
            }
            else if (substr(BB[k], 1, 7) == "Utheta_") {
                pa = strpos(BR[k], "~~")
                i  = strtoreal(substr(BR[k], 2, pa - 2))
                j  = strtoreal(substr(BR[k], pa + 3, .))
                pk = (i - 1) * m.p - i * (i - 1) / 2 + (j - i)
                tc = dwlsmv_t_cx_th(BC[k], K)
                if (tc[1] == i) {
                    av[k] = Utm[pk, tc[2]]
                }
                else if (tc[1] == j) {
                    av[k] = Utm[pk, K[i] + tc[2]]
                }
            }
        }
        blk = missing(av)
        if (blk) {
            printf("{txt}%s: %s block rows could not be matched (diagnostic)\n", tag, strofreal(blk))
        }
        else {
            gabs = max(abs(av :- BV))
            grel = max(abs(av :- BV) :/ abs(BV))
            printf("%s: blocks max |gap| %9.3e, max rel gap %9.3e, over %s entries (diagnostic only)\n",
                   tag, gabs, grel, strofreal(rows(BV)))
        }
    }

    // ---- 5. negative controls: each must FAIL the Gamma tier ---------------
    if (!missing(Gr)) {
        for (ctl = 1; ctl <= 4; ctl++) {
            Psi = dwlsmv_cx_psi_matrix(Y, X, m, st.s, ctl)
            Gx  = quadcross(Psi, Psi) :/ rows(Y)
            rr  = dwlsmv_gapratio(Gx, Gr, 1e-8, 1e-8)
            printf("%s: control %s, Gamma gap ratio %9.3e (must exceed 1)\n",
                   tag, strofreal(ctl), rr)
            if (missing(rr) | rr <= 1) {
                FAIL = FAIL \ (tag + ": negative control " + strofreal(ctl) +
                               " passes the Gamma tier; the check cannot see it")
            }
        }
    }

    for (k = 1; k <= rows(FAIL); k++) {
        printf("{err}%s\n", FAIL[k])
    }
    return(rows(FAIL))
}

// Check 3: the unconditional map is untouched; the conditional map's layout.
real scalar dwlsmv_t_cx_map()
{
    struct dwlsmv_statmap scalar m0, m1, m2
    string colvector FAIL, want
    real colvector K
    real rowvector v
    real scalar i, j, s

    FAIL = J(0, 1, "")
    K    = (4 \ 4 \ 3 \ 3)
    m0   = dwlsmv_statmap_build(K)
    m1   = dwlsmv_statmap_build(K, 0)

    if (m0.p != m1.p | m0.T != m1.T | m0.P != m1.P | m0.D != m1.D) {
        FAIL = FAIL \ "3: the scalars of build(K) and build(K, 0) differ"
    }
    if (m0.K != m1.K | m0.C != m1.C | m0.Coff != m1.Coff) {
        FAIL = FAIL \ "3: K, C or Coff of build(K) and build(K, 0) differ"
    }
    if (m0.stat_type != m1.stat_type | m0.var_i != m1.var_i |
        m0.var_j != m1.var_j | m0.thr_k != m1.thr_k | m0.label != m1.label) {
        FAIL = FAIL \ "3: the rows of build(K) and build(K, 0) differ"
    }
    if (m0.D != m0.T + m0.P | dwlsmv_statmap_q(m0) != 0) {
        FAIL = FAIL \ "3: an unconditional map does not have D = T + P and q = 0"
    }
    if (any(m0.stat_type :== 3)) {
        FAIL = FAIL \ "3: an unconditional map has a slope row"
    }
    for (i = 1; i <= m0.p - 1; i++) {
        for (j = i + 1; j <= m0.p; j++) {
            if (dwlsmv_pos_rho(m0, i, j) != m0.T + (i - 1) * m0.p - i * (i - 1) / 2 + (j - i)) {
                FAIL = FAIL \ ("3: pos_rho moved on an unconditional map at (" +
                               strofreal(i) + ", " + strofreal(j) + ")")
            }
        }
    }
    for (s = 1; s <= m0.D; s++) {
        v = dwlsmv_pos_inv(m0, s)
        if (v != (m0.stat_type[s], m0.var_i[s], m0.var_j[s], m0.thr_k[s])) {
            FAIL = FAIL \ ("3: pos_inv does not round-trip position " + strofreal(s) + " (q = 0)")
        }
    }

    m2   = dwlsmv_statmap_build((3 \ 3 \ 2), 1)
    want = ("y1|t1" \ "y1|t2" \ "y2|t1" \ "y2|t2" \ "y3|t1" \
            "y1~x1" \ "y2~x1" \ "y3~x1" \ "y1~~y2" \ "y1~~y3" \ "y2~~y3")
    if (m2.D != 11 | rows(m2.label) != 11) {
        FAIL = FAIL \ "3: the conditional map is not 11 long"
    }
    else if (m2.label != want) {
        FAIL = FAIL \ "3: the conditional map's labels are not spec 4.2's order"
    }
    for (j = 1; j <= 3; j++) {
        if (dwlsmv_pos_slope(m2, j, 1) != 5 + j) {
            FAIL = FAIL \ ("3: pos_slope(" + strofreal(j) + ", 1) is not T + j")
        }
    }
    if (dwlsmv_pos_rho(m2, 1, 2) != 9 | dwlsmv_pos_rho(m2, 2, 3) != 11) {
        FAIL = FAIL \ "3: pos_rho is not shifted by pq on the conditional map"
    }
    for (s = 1; s <= m2.D; s++) {
        v = dwlsmv_pos_inv(m2, s)
        if (v != (m2.stat_type[s], m2.var_i[s], m2.var_j[s], m2.thr_k[s])) {
            FAIL = FAIL \ ("3: pos_inv does not round-trip position " + strofreal(s) + " (q = 1)")
        }
    }

    for (s = 1; s <= rows(FAIL); s++) {
        printf("{err}%s\n", FAIL[s])
    }
    return(rows(FAIL))
}

// Check 4: one refusal per call; the result is assigned so that nothing is
// printed on the (wrong) success path.
void dwlsmv_t_cx_err(real matrix Y, real colvector X, real scalar c)
{
    struct dwlsmv_stage1 scalar st
    real matrix Yb, Xb

    Yb = Y
    Xb = X
    if (c == 1) {
        Xb = X :+ 1
    }
    else if (c == 2) {
        Xb = 2 :* X :- 1
    }
    else if (c == 3) {
        Xb = J(rows(X), 1, 0)
    }
    else if (c == 4) {
        Xb[5] = .
    }
    else if (c == 5) {
        Xb = (X, X)
    }
    else if (c == 6) {
        Xb = X[|1 \ rows(X) - 1|]
    }
    else if (c == 7) {
        Xb[1] = 2
    }
    else if (c == 8) {
        Yb[., 3] = X
    }
    else if (c == 9) {
        Yb[., 2] = Y[., 1]
    }
    else if (c == 10) {
        // P1 boundary, refused: the groups share exactly one category
        // (max0 = min1), the example of @sec-condx-policy.
        Yb[., 1] = mod((1::rows(Y)), 2) :+ X
    }
    else if (c == 11) {
        // P1 boundary, finite: one group in a single middle category, the
        // other in the two outer ones (max0 > min1 and max1 > min0).
        Yb[., 1] = (1 :- X) :+ X :* (2 :* mod((1::rows(Y)), 2))
    }
    st = dwls_sample_stats_cx(Yb, (3 \ 3 \ 2), Xb)
}

void dwlsmv_t_cx_err_slope()
{
    struct dwlsmv_statmap scalar m
    real scalar p

    m = dwlsmv_statmap_build((3 \ 3 \ 2))
    p = dwlsmv_pos_slope(m, 1, 1)
}

end

* An abort inside a check (a refusal on the reference data, a load-order
* problem) must count as a failure, not leave the local empty.
foreach t in A D {
    local n`t' 1
    capture noisily mata: st_local("n`t'", strofreal(dwlsmv_t_cx_one(Y_`t', X_`t', ///
        SL_`t', SV_`t', GR_`t', GC_`t', GV_`t', BB_`t', BR_`t', BC_`t', BV_`t', "condx_`t'")))
    if _rc {
        display as error "test_condx_stage1: the condx_`t' checks aborted, rc " _rc
        local n`t' 1
    }
}
local nM 1
capture noisily mata: st_local("nM", strofreal(dwlsmv_t_cx_map()))
if _rc {
    display as error "test_condx_stage1: the map checks aborted, rc " _rc
    local nM 1
}

* The v1 entry point still runs on case_D after the shared map files changed.
* Its values are pinned by test_gamma.do and test_sample_stats.do, and its
* bitwise identity by the snapshot pair (header, check 3).
* A SUCCESS check, not an error path: any nonzero rc is a failure here.
local nV 0
capture noisily mata: (void) dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
if _rc {
    display as error "test_condx_stage1: dwls_sample_stats() on case_D failed, rc " _rc
    local nV 1
}

* ---- 4. refusals, each at rc 3498 ------------------------------------------
local nerr 0
forvalues c = 1/10 {
    capture mata: dwlsmv_t_cx_err(Y_A, X_A, `c')
    if _rc != 3498 {
        display as error "test_condx_stage1: refusal case `c' gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}
* Case 11 is a SUCCESS check, the other side of the P1 boundary: the item is
* finite (@sec-condx-policy), so Stage 1 must not refuse it. The R reference
* and the frozen cxs_rep refuse it; that difference is recorded in
* E-103.
capture noisily mata: dwlsmv_t_cx_err(Y_A, X_A, 11)
if _rc {
    display as error "test_condx_stage1: the finite P1 boundary case was refused, rc " _rc
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_t_cx_err_slope()
if _rc != 3498 {
    display as error "test_condx_stage1: pos_slope on an unconditional map gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

* The messages of cases 1, 4, 8 and 9, read back from a log.
tempfile L
foreach c in 1 4 8 9 {
    quietly log using "`L'_`c'.log", text replace name(cx)
    capture noisily mata: dwlsmv_t_cx_err(Y_A, X_A, `c')
    quietly log close cx
}

mata:
real scalar dwlsmv_t_cx_msg(string scalar f, string rowvector has, string scalar tag)
{
    string scalar txt
    real scalar i, bad

    txt = invtokens(cat(f)', char(10))
    bad = 0
    for (i = 1; i <= cols(has); i++) {
        if (!strpos(txt, has[i])) {
            printf("{err}test_condx_stage1: case %s message lacks: %s\n", tag, has[i])
            bad = bad + 1
        }
    }
    if (strpos(txt, "\" + "n")) {
        printf("{err}test_condx_stage1: case %s message has a literal backslash-n\n", tag)
        bad = bad + 1
    }
    return(bad)
}
end

local tmpL "`L'"
local nmsg 1
capture noisily mata: st_local("nmsg", strofreal(                                                          ///
    dwlsmv_t_cx_msg(st_local("tmpL") + "_1.log", ("takes the values 1 and 2, not 0 and 1", ///
        "coded exactly 0/1"), "1") +                                                        ///
    dwlsmv_t_cx_msg(st_local("tmpL") + "_4.log", ("1 of 1500 values of the covariate are missing", ///
        "Delete listwise across the covariate AND all indicators"), "4") +                  ///
    dwlsmv_t_cx_msg(st_local("tmpL") + "_8.log", ("separation by the covariate",            ///
        "y3: categories observed at x = 0: 0; at x = 1: 1"), "8") +                         ///
    dwlsmv_t_cx_msg(st_local("tmpL") + "_9.log", ("y1~~y2", "status 2",                    ///
        "no cell repair is applied"), "9")))
if _rc {
    display as error "test_condx_stage1: the message checks aborted, rc " _rc
    local nmsg 1
}

* ---- 6. item stopping rule: dwlsmv_cx_item_stop and dwlsmv_cx_item_fit -----
mata:

// One helper case: 1 on a mismatch against want, else 0.
real scalar dwlsmv_t_stop_one(real colvector g, real colvector st,
                              real scalar tol, real scalar want,
                              string scalar tag)
{
    real scalar got

    got = dwlsmv_cx_item_stop(g, st, tol)
    if (got != want) {
        printf("{err}test_condx_stage1: item stop %s gave %s, not %s\n",
               tag, strofreal(got), strofreal(want))
        return(1)
    }
    return(0)
}

// Returns the number of failed checks in section 6.
real scalar dwlsmv_t_stop_all(real matrix Y, real colvector X)
{
    real colvector g, st, K, th
    real scalar j, it, nf

    nf = 0

    // (a) the helper, every case against its expected value
    // a1: small non-missing step, both entries under tol
    g  = (1e-3 \ 1e-3)
    st = (1e-13 \ -5e-13)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 1, "a1 small step")

    // a2: a step entry over tol
    st = (2e-12 \ 0)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 0, "a2 step over tol")

    // a3: each entry under tol, the sum over tol (kills sum() for max())
    st = (6e-13 \ 6e-13)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 1, "a3 sum over tol, max under")

    // a4: a large negative entry (kills max(step) for max(abs(step)))
    st = (-2e-12 \ 1e-13)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 0, "a4 negative entry over tol")

    // a5: small step, large score (kills gbar in place of step)
    g  = (1e-3 \ 0)
    st = (1e-13 \ 0)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 1, "a5 large score, small step")

    // a6: tiny score, large step (kills gbar in place of step)
    g  = (1e-15 \ 0)
    st = (2e-12 \ 0)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 0, "a6 tiny score, large step")

    // a7: one missing step entry, the other tiny (kills missing(step))
    g  = (1e-3 \ 1e-3)
    st = (1e-15 \ .)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 0, "a7 one step entry missing")

    // a8: one missing score entry, the step tiny (kills missing(gbar))
    g  = (1e-15 \ .)
    st = (1e-13 \ 1e-13)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 0, "a8 one score entry missing")

    // a9: the whole step missing
    st = (. \ .)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 0, "a9 step all missing")

    // a10: exactly at tol
    g  = (1e-3 \ 1e-3)
    st = (1e-12 \ -1e-12)
    nf = nf + dwlsmv_t_stop_one(g, st, 1e-12, 1, "a10 exactly at tol")

    // (b) the fit on condx_A's three items: converges, 2 <= it <= 49.
    // A fit that refuses aborts this whole call; the Stata side then counts
    // section 6 as failed (see the block below).
    K = (3 \ 3 \ 2)
    for (j = 1; j <= 3; j++) {
        it = .
        th = dwlsmv_cx_item_fit(Y[., j], X, K[j], j, it)
        if (rows(th) != K[j]) {
            printf("{err}test_condx_stage1: item fit y%s returned %s rows, not %s\n",
                   strofreal(j), strofreal(rows(th)), strofreal(K[j]))
            nf = nf + 1
        }
        else if (missing(th)) {
            printf("{err}test_condx_stage1: item fit y%s has missing thresholds\n",
                   strofreal(j))
            nf = nf + 1
        }
        else if (it < 2 | it > 49) {
            printf("{err}test_condx_stage1: item fit y%s took %s iterations, not 2 to 49\n",
                   strofreal(j), strofreal(it))
            nf = nf + 1
        }
    }
    return(nf)
}


// (c) the iterate-50 branch and the refusal, through the tolc hook.
// what 1: tolc = (-1, 1e-10), so the loop cannot stop and iterate 50 decides;
//         it must be 50 and the point within 1.1e-12 of production's (the
//         d1 certificate; both are within it of the maximiser, and the
//         iterate-50 point sits at the arithmetic's floor).
// what 2: tolc = (-1, -1): nothing can pass; must abort (rc 3498, checked by
//         the caller). Returns the failure count for what 1.
real scalar dwlsmv_t_stop_fb(real matrix Y, real colvector X, real scalar what)
{
    real colvector K, th, th50
    real scalar j, it, nf, gap

    nf = 0
    K  = (3 \ 3 \ 2)
    if (what == 2) {
        it = .
        th = dwlsmv_cx_item_fit(Y[., 1], X, K[1], 1, it, (-1, -1))
        return(1)
    }
    for (j = 1; j <= 3; j++) {
        it   = .
        th   = dwlsmv_cx_item_fit(Y[., j], X, K[j], j, it)
        it   = .
        th50 = dwlsmv_cx_item_fit(Y[., j], X, K[j], j, it, (-1, 1e-10))
        if (it != 50) {
            printf("{err}test_condx_stage1: forced fallback y%s stopped at iterate %s, not 50\n",
                   strofreal(j), strofreal(it))
            nf = nf + 1
        }
        gap = max(abs(th :- th50))
        printf("{txt}  item stop: y%s production vs iterate 50, max gap %9.2e (limit 1.1e-12)\n",
               strofreal(j), gap)
        if (missing(th) + missing(th50) > 0) {
            printf("{err}test_condx_stage1: forced fallback y%s has missing estimates\n",
                   strofreal(j))
            nf = nf + 1
        }
        else if (gap > 1.1e-12) {
            printf("{err}test_condx_stage1: y%s production is %9.2e from the iterate-50 point\n",
                   strofreal(j), gap)
            nf = nf + 1
        }
    }
    return(nf)
}

end

local nstop 1
capture noisily mata: st_local("nstop", strofreal(dwlsmv_t_stop_all(Y_A, X_A) + dwlsmv_t_stop_fb(Y_A, X_A, 1) + dwlsmv_t_stop_fb(Y_D, X_D, 1)))
if _rc {
    display as error "test_condx_stage1: the item stop checks aborted, rc " _rc
    local nstop 1
}
capture mata: dwlsmv_t_stop_fb(Y_A, X_A, 2)
if _rc != 3498 {
    display as error "test_condx_stage1: an unreachable step tolerance gave rc " _rc ", not 3498"
    local nstop = `nstop' + 1
}
capture noisily mata: dwlsmv_t_stop_fb(Y_A, X_A, 2)

* ---- 7. item fit: threshold ordering and the tolc hook (targeted) -----------
* The three tables are CONSTRUCTED, AND CHOSEN AFTER MEASURING: T2 (a narrow
* middle category, thresholds 0.027 apart) and T6, T10 were picked from ten
* tried because each has an iterate whose step lies between 1e-12 and 1e-10,
* which is what lets a test see the two production constants. They are not
* data of the project.
mata:

// y and x from a 2 x K table of counts (row 1: x = 0).
void dwlsmv_t_tab2yx(real matrix n, real colvector y, real colvector x)
{
    real scalar g, c

    y = J(0, 1, .)
    x = J(0, 1, .)
    for (g = 1; g <= 2; g++) {
        for (c = 1; c <= cols(n); c++) {
            y = y \ J(n[g, c], 1, c - 1)
            x = x \ J(n[g, c], 1, g - 1)
        }
    }
}

// 1 unless the thresholds of th (its first K - 1 rows) increase strictly and
// nothing is missing. Written out here, not through dwlsmv_cx_item_ok().
real scalar dwlsmv_t_unordered(real colvector th, real scalar K)
{
    real scalar c

    if (missing(th) > 0) {
        return(1)
    }
    for (c = 2; c <= K - 1; c++) {
        if (th[c] <= th[c - 1]) {
            return(1)
        }
    }
    return(0)
}

real scalar dwlsmv_t_ok_one(real colvector th, real scalar K, real scalar want,
                            string scalar tag)
{
    if (dwlsmv_cx_item_ok(th, K) != want) {
        printf("{err}test_condx_stage1: item_ok %s did not give %s\n", tag, strofreal(want))
        return(1)
    }
    return(0)
}

// what 1: the checks; returns the failure count.
// what 2..5: a malformed tolc; each must abort at rc 3498 (caller checks).
real scalar dwlsmv_t_stop_tgt(real scalar what)
{
    pointer(real matrix) rowvector T
    real colvector y, x, th, tha, thb, th50
    real scalar t, K, it, ita, itb, it50, nf

    T = (&(700, 5, 95 \ 400, 8, 392), &(300, 200 \ 150, 350),
         &(120, 260, 120 \ 60, 240, 200))
    if (what > 1) {
        dwlsmv_t_tab2yx(*T[1], y, x)
        it = .
        if (what == 2) {
            th = dwlsmv_cx_item_fit(y, x, 3, 1, it, 1e-12)
        }
        if (what == 3) {
            th = dwlsmv_cx_item_fit(y, x, 3, 1, it, (1e-12, 1e-10, 1e-10))
        }
        if (what == 4) {
            th = dwlsmv_cx_item_fit(y, x, 3, 1, it, (1e-12, .))
        }
        if (what == 5) {
            th = dwlsmv_cx_item_fit(y, x, 3, 1, it, (1e-12 \ 1e-10))
        }
        return(1)
    }

    nf = 0
    // (a) dwlsmv_cx_item_ok: the last row is the slope and is not ordered
    nf = nf + dwlsmv_t_ok_one((-1 \ 1 \ 0), 3, 1, "increasing")
    nf = nf + dwlsmv_t_ok_one((-1 \ 1 \ -5), 3, 1, "increasing, slope below both")
    nf = nf + dwlsmv_t_ok_one((1 \ 1 \ 0), 3, 0, "equal thresholds")
    nf = nf + dwlsmv_t_ok_one((1 \ -1 \ 0), 3, 0, "reversed thresholds")
    nf = nf + dwlsmv_t_ok_one((-1 \ . \ 0), 3, 0, "a missing threshold")
    nf = nf + dwlsmv_t_ok_one((-1 \ 1 \ .), 3, 0, "a missing slope")
    nf = nf + dwlsmv_t_ok_one((-1 \ 0 \ 1 \ 0.5), 4, 1, "K = 4 increasing")
    nf = nf + dwlsmv_t_ok_one((-1 \ 0 \ 0 \ 0.5), 4, 0, "K = 4, last pair equal")
    nf = nf + dwlsmv_t_ok_one((-1 \ 0 \ -0.5 \ 0.5), 4, 0, "K = 4, last pair reversed")
    nf = nf + dwlsmv_t_ok_one((0.3 \ -2), 2, 1, "K = 2, one threshold")

    for (t = 1; t <= 3; t++) {
        dwlsmv_t_tab2yx(*T[t], y, x)
        K    = cols(*T[t])
        it   = .
        th   = dwlsmv_cx_item_fit(y, x, K, 1, it)
        ita  = .
        tha  = dwlsmv_cx_item_fit(y, x, K, 1, ita, (1e-12, 1e-10))
        itb  = .
        thb  = dwlsmv_cx_item_fit(y, x, K, 1, itb, (1e-10, 1e-10))
        it50 = .
        th50 = dwlsmv_cx_item_fit(y, x, K, 1, it50, (-1, 1e-10))
        printf("{txt}  item fit T%g: it %g; with (1e-12, 1e-10) %g; with (1e-10, 1e-10) %g; forced %g\n",
               t, it, ita, itb, it50)
        // (b) the production constants: omitting tolc is (1e-12, 1e-10)
        // exactly, and a looser d1 stops earlier on these tables
        if (th != tha) {
            printf("{err}test_condx_stage1: T%g: omitting tolc is not (1e-12, 1e-10)\n", t)
            nf = nf + 1
        }
        if (it != ita) {
            printf("{err}test_condx_stage1: T%g: iterations differ from tolc (1e-12, 1e-10)\n", t)
            nf = nf + 1
        }
        if (!(itb < ita)) {
            printf("{err}test_condx_stage1: T%g: d1 = 1e-10 did not stop earlier\n", t)
            nf = nf + 1
        }
        if (it50 != 50) {
            printf("{err}test_condx_stage1: T%g: the forced fit stopped at %g\n", t, it50)
            nf = nf + 1
        }
        // (c) every returned point has strictly increasing thresholds
        if (dwlsmv_t_unordered(th, K) + dwlsmv_t_unordered(thb, K) +
            dwlsmv_t_unordered(th50, K) > 0) {
            printf("{err}test_condx_stage1: T%g: a returned point is unordered or missing\n", t)
            nf = nf + 1
        }
        if (max(abs(th :- th50)) > 1.1e-12) {
            printf("{err}test_condx_stage1: T%g: production is %9.2e from the iterate-50 point\n",
                   t, max(abs(th :- th50)))
            nf = nf + 1
        }
    }
    return(nf)
}

end

local ntgt 1
capture noisily mata: st_local("ntgt", strofreal(dwlsmv_t_stop_tgt(1)))
if _rc {
    display as error "test_condx_stage1: the targeted item-fit checks aborted, rc " _rc
    local ntgt 1
}
forvalues w = 2/5 {
    capture mata: dwlsmv_t_stop_tgt(`w')
    if _rc != 3498 {
        display as error "test_condx_stage1: malformed tolc case `w' gave rc " _rc ", not 3498"
        local ntgt = `ntgt' + 1
    }
}

* ---- 8. refusals named by the decision of 2026-11-12 (R4 - R11) ----------------
* Drafted by a read-only agent 2026-11-12; run by the main session. Case ids are numbers on purpose: no macro appears inside a Mata string.
mata:

// Counts the phrases in has missing from the log and the phrases in hasnot present in it.
real scalar dwlsmv_t_s8_msg(string scalar f, string rowvector has,
                            string rowvector hasnot, string scalar tag)
{
    string scalar txt
    real scalar i, bad

    txt = invtokens(cat(f)', char(10))
    bad = 0
    for (i = 1; i <= cols(has); i++) {
        if (!strpos(txt, has[i])) {
            printf("{err}test_condx_stage1: case %s message lacks: %s\n", tag, has[i])
            bad = bad + 1
        }
    }
    for (i = 1; i <= cols(hasnot); i++) {
        if (strpos(txt, hasnot[i])) {
            printf("{err}test_condx_stage1: case %s message has, and must not: %s\n", tag, hasnot[i])
            bad = bad + 1
        }
    }
    if (strpos(txt, "\" + "n")) {
        printf("{err}test_condx_stage1: case %s message has a literal backslash-n\n", tag)
        bad = bad + 1
    }
    return(bad)
}

// The phrases of each case (see the header). Outputs by reference.
void dwlsmv_t_s8_exp(real scalar c, string rowvector has, string rowvector hasnot)
{
    has    = J(1, 0, "")
    hasnot = J(1, 0, "")
    if (c == 1) {
        has    = ("item y2 has an empty category in the pooled sample",
                  "no observations in category 2 of 3")
        hasnot = ("item y1 has")
    }
    else if (c == 2) {
        has    = ("item y1 has an empty category in the pooled sample",
                  "no observations in category 1 of 3", "INTERIOR category")
        hasnot = ("item y2 has")
    }
    else if (c == 3) {
        has    = ("y2~~y3",
                  "no positive probability at the observed cell (x = 1, y2 = 1, y3 = 2)")
        hasnot = ("a missing value")
    }
    else if (c == 4) {
        has    = ("y2~~y3",
                  "a missing value at the observed cell (x = 0, y2 = 0, y3 = 0)")
        hasnot = ("no positive probability at")
    }
    else if (c == 5) {
        has    = ("y2~~y4", "rho = 0 is not usable in the group x = 1",
                  "dwls_polychoric: rho = 0 is not usable")
        hasnot = ("group x = 0")
    }
    else if (c == 6) {
        has    = ("y2~~y4", "rho = 0 is not usable in the group x = 0",
                  "dwls_polychoric: rho = 0 is not usable")
        hasnot = ("group x = 1")
    }
    else if (c == 7) {
        has    = ("y2~~y4", "in the group x = 0", "in the group x = 1",
                  "dwls_polychoric: rho = 0 is not usable")
    }
    else if (c == 8) {
        has    = ("item y2 did not converge", "after 50 iterations")
        hasnot = ("item y1 did not")
    }
    else if (c == 9) {
        has    = ("item y2's H_j is not positive definite")
        hasnot = ("item y1's H_j")
    }
    else if (c == 10) {
        has    = ("y1~~y2", "not negative")
    }
    else if (c == 11) {
        has    = ("1 of 3 pairs cannot enter the sample statistics",
                  "y1~~y2: residual rho-hat", "status 2")
    }
    else if (c == 12) {
        has    = ("dwlsmv_gammadiag", "conditional")
    }
    else if (c == 13) {
        hasnot = ("empty category in the pooled sample", "has no starting point",
                  "at the observed cell")
    }
}

// R7: a pooled empty category. c = 1: item 2 loses its top category (boundary);
// c = 2: item 1 is recoded to {0, 2} (interior). Both are refused by the pooled-category
// loop before the separation rule, which is why neither can be a P1 refusal.
void dwlsmv_t_s8_r7(real matrix Y, real colvector X, real scalar c)
{
    struct dwlsmv_stage1 scalar st
    real matrix Yb

    Yb = Y
    if (c == 1) {
        Yb[., 2] = (Y[., 2] :< 2) :* Y[., 2] :+ (Y[., 2] :>= 2)
    }
    else {
        Yb[., 1] = 2 :* (Y[., 1] :>= 1)
    }
    st = dwls_sample_stats_cx(Yb, (3 \ 3 \ 2), X)
}

// R6: dwlsmv_cx_pair_blocks called directly with crafted thresholds; (i, k) = (2, 3).
// Observed rows (x Ki Kj + c_i Kj + c_j + 1): 1 = (x 0, 0, 0), 9 = (x 0, 2, 2),
// 15 = (x 1, 1, 2), 17 = (x 1, 2, 1). c = 1: item i's thresholds crossed, so category 1 has
// a negative probability (row 15). c = 2: item i's first threshold missing, so row 1 is
// missing. The missing case rests on binormal() returning missing for a missing argument
// (not verified).
void dwlsmv_t_s8_r6(real scalar c)
{
    real colvector yi, yj, x, thi, thj, n
    real rowvector Ut
    real matrix    tab
    real scalar    Ur

    yi  = (0 \ 2 \ 1 \ 2)
    yj  = (0 \ 1 \ 2 \ 2)
    x   = (0 \ 1 \ 1 \ 0)
    thj = (-0.5 \ 0.5 \ 0)
    if (c == 1) {
        thi = (0.5 \ -0.5 \ 0)
    }
    else {
        thi = (. \ 0.5 \ 0)
    }
    dwlsmv_cx_pair_blocks(yi, yj, x, thi, thj, 0.3, 0, Ur, Ut, tab, n, (2, 3))
}

// R8: dwlsmv_cx_pair_fit called directly; (i, k) = (2, 4). c = 1: x = 1 bounds (-13, -11)
// so group 1 is unusable; c = 2: x = 0 bounds (-20, -19), x = 1 bounds (0, 1), so only
// group 0; c = 3: both groups' bounds (-20, -19).
void dwlsmv_t_s8_r8(real scalar c)
{
    real colvector yi, yj, x, thi, thj
    real rowvector r

    yi  = (0 \ 1 \ 2 \ 0 \ 1 \ 2)
    yj  = (0 \ 1 \ 2 \ 1 \ 0 \ 2)
    x   = (0 \ 0 \ 0 \ 1 \ 1 \ 1)
    thj = (-1 \ 1 \ 0)
    if (c == 1) {
        thi = (-1 \ 1 \ 12)
    }
    else if (c == 2) {
        thi = (-20 \ -19 \ -20)
    }
    else {
        thi = (-20 \ -19 \ 0)
    }
    r = dwlsmv_cx_pair_fit(yi, yj, x, thi, thj, (2, 4))
}

// R9: the iterate-50 refusal through the tolc hook. j = 2 is passed on item 1's data (j is for
// messages only), so the message must carry j, not a constant.
void dwlsmv_t_s8_r9(real matrix Y, real colvector X)
{
    real colvector th
    real scalar it

    it = .
    th = dwlsmv_cx_item_fit(Y[., 1], X, 3, 2, it, (-1, -1))
}

// R4: item 2's two thresholds made equal in s: its middle category has probability exactly 0
// and is observed, so H_j is missing. Relies on cholinv() returning missing for it.
void dwlsmv_t_s8_r4(real matrix Y, real colvector X)
{
    struct dwlsmv_stage1 scalar st
    real colvector s
    real matrix Psi

    st = dwls_sample_stats_cx(Y, (3 \ 3 \ 2), X)
    s  = st.s
    s[dwlsmv_pos_tau(st.map, 2, 2)] = s[dwlsmv_pos_tau(st.map, 2, 1)]
    Psi = dwlsmv_cx_psi_matrix(Y, X, st.map, s, 0)
}

// R5: item 2 := item 1, so only diagonal cells are observed; s from the items' own fits and
// rho(1, 2) = 0.9, where U_rho is clearly positive (R finite differences on condx_A: 6.3).
void dwlsmv_t_s8_r5(real matrix Y, real colvector X)
{
    struct dwlsmv_statmap scalar m
    real matrix Yd, Psi
    real colvector K, s, th
    real scalar j, c, it

    K  = (3 \ 3 \ 2)
    m  = dwlsmv_statmap_build(K, 1)
    Yd = Y
    Yd[., 2] = Y[., 1]
    s  = J(m.D, 1, 0)
    for (j = 1; j <= 3; j++) {
        it = .
        th = dwlsmv_cx_item_fit(Yd[., j], X, K[j], j, it)
        for (c = 1; c <= K[j] - 1; c++) {
            s[dwlsmv_pos_tau(m, j, c)] = th[c]
        }
        s[dwlsmv_pos_slope(m, j, 1)] = th[K[j]]
    }
    s[dwlsmv_pos_rho(m, 1, 2)] = 0.9
    Psi = dwlsmv_cx_psi_matrix(Yd, X, m, s, 0)
}

// R11: the diagnostic handed a conditional Stage 1 result.
void dwlsmv_t_s8_r11(real matrix Y, real colvector X)
{
    struct dwlsmv_stage1 scalar st
    real matrix A

    st = dwls_sample_stats_cx(Y, (3 \ 3 \ 2), X)
    A  = dwlsmv_gammadiag(st, Y)
}

// The unconditional diagnostic is untouched: 1 on a failed check, else 0.
real scalar dwlsmv_t_s8_v1gd(real matrix YD)
{
    struct dwlsmv_stage1 scalar st
    real matrix A

    st = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    A  = dwlsmv_gammadiag(st, YD)
    if (rows(A) != 7) {
        return(1)
    }
    if (cols(A) != 4) {
        return(1)
    }
    return(0)
}

// Dispatcher. Cases 1..12 must abort; 13 is a success (the quiet-success check).
void dwlsmv_t_s8(real matrix Y, real colvector X, real scalar c)
{
    struct dwlsmv_stage1 scalar st

    if (c <= 2) {
        dwlsmv_t_s8_r7(Y, X, c)
    }
    else if (c <= 4) {
        dwlsmv_t_s8_r6(c - 2)
    }
    else if (c <= 7) {
        dwlsmv_t_s8_r8(c - 4)
    }
    else if (c == 8) {
        dwlsmv_t_s8_r9(Y, X)
    }
    else if (c == 9) {
        dwlsmv_t_s8_r4(Y, X)
    }
    else if (c == 10) {
        dwlsmv_t_s8_r5(Y, X)
    }
    else if (c == 11) {
        dwlsmv_t_cx_err(Y, X, 9)
    }
    else if (c == 12) {
        dwlsmv_t_s8_r11(Y, X)
    }
    else {
        st = dwls_sample_stats_cx(Y, (3 \ 3 \ 2), X)
    }
}

// All logs 1..14 (14 = case 13 on condx_D): message content. stem = the tempfile path.
real scalar dwlsmv_t_s8_msgs(string scalar stem)
{
    string rowvector has, hasnot
    real scalar c, e, bad

    bad = 0
    for (c = 1; c <= 14; c++) {
        e = c
        if (c == 14) {
            e = 13
        }
        dwlsmv_t_s8_exp(e, has, hasnot)
        bad = bad + dwlsmv_t_s8_msg(stem + "_s8_" + strofreal(c) + ".log", has, hasnot,
                                    "s8." + strofreal(c))
    }
    return(bad)
}

end

local n8 0
forvalues c = 1/13 {
    quietly log using "`L'_s8_`c'.log", text replace name(cx)
    capture noisily mata: dwlsmv_t_s8(Y_A, X_A, `c')
    local rc8 = _rc
    quietly log close cx
    if `c' < 13 & `rc8' != 3498 {
        display as error "test_condx_stage1: section 8 case `c' gave rc `rc8', not 3498"
        local n8 = `n8' + 1
    }
    if `c' == 13 & `rc8' != 0 {
        display as error "test_condx_stage1: section 8 case 13 (Stage 1 on condx_A) gave rc `rc8'"
        local n8 = `n8' + 1
    }
}
quietly log using "`L'_s8_14.log", text replace name(cx)
capture noisily mata: dwlsmv_t_s8(Y_D, X_D, 13)
local rc8 = _rc
quietly log close cx
if `rc8' != 0 {
    display as error "test_condx_stage1: section 8 case 14 (Stage 1 on condx_D) gave rc `rc8'"
    local n8 = `n8' + 1
}

local nmsg8 1
capture noisily mata: st_local("nmsg8", strofreal(dwlsmv_t_s8_msgs(st_local("tmpL"))))
if _rc {
    display as error "test_condx_stage1: the section 8 message checks aborted, rc " _rc
    local nmsg8 1
}
local n8 = `n8' + `nmsg8'

* A SUCCESS check, not an error path: any nonzero rc is a failure here.
local n8v 1
capture noisily mata: st_local("n8v", strofreal(dwlsmv_t_s8_v1gd(YD)))
if _rc {
    display as error "test_condx_stage1: dwlsmv_gammadiag() on case_D failed, rc " _rc
    local n8v 1
}
local n8 = `n8' + `n8v'

* rho = 0 unusable, reached FROM DATA (user, 2026-11-12: the abort is the guard
* for K >= 3). Item 1 has counts (0, m, 1) at x = 0 and (1, 0, m) at x = 1 with
* m = 1000: it passes P1 and fits, but its one observation in category 0 at
* x = 1 has a fitted probability near 9e-17, below delta_P, so rho = 0 is not
* usable in that group for any pair of the item. The table was chosen after
* measuring (m = 500 gives 2e-14, above delta_P). This is also the only case
* that reaches the refusal through dwls_sample_stats_cx(), so it is what sees
* the pair being passed at that call.
mata:

void dwlsmv_t_s8_rho0(real scalar m)
{
    struct dwlsmv_stage1 scalar st
    real colvector y1, y2, x

    y1 = J(m, 1, 1) \ 2 \ 0 \ J(m, 1, 2)
    x  = J(m + 1, 1, 0) \ J(m + 1, 1, 1)
    y2 = mod(1::rows(y1), 2)
    st = dwls_sample_stats_cx((y1, y2), (3 \ 2), x)
}

real scalar dwlsmv_t_s8_rho0msg(string scalar f)
{
    string rowvector has, hasnot

    has    = ("the pair y1~~y2 has no starting point",
              "rho = 0 is not usable in the group x = 1",
              "dwls_polychoric: rho = 0 is not usable")
    hasnot = J(1, 1, "group x = 0")
    return(dwlsmv_t_s8_msg(f, has, hasnot, "s8.rho0"))
}

end

local f15 "`L'_s8_15.log"
quietly log using "`f15'", text replace name(cx)
capture noisily mata: dwlsmv_t_s8_rho0(1000)
local rc8 = _rc
quietly log close cx
if `rc8' != 3498 {
    display as error "test_condx_stage1: the rho = 0 case from data gave rc `rc8', not 3498"
    local n8 = `n8' + 1
}
local n15 1
capture noisily mata: st_local("n15", strofreal(dwlsmv_t_s8_rho0msg(st_local("f15"))))
if _rc {
    display as error "test_condx_stage1: the rho = 0 message check aborted, rc " _rc
    local n15 1
}
local n8 = `n8' + `n15'

if `nA' | `nD' | `nM' | `nV' | `nerr' | `nmsg' | `nstop' | `ntgt' | `n8' {
    display as error _n "test_condx_stage1: condx_A `nA', condx_D `nD', map `nM', v1 `nV', refusals `nerr', messages `nmsg', item stop `nstop', targeted `ntgt', named refusals `n8' failure(s)"
    exit 9
}

display as result "test_condx_stage1: conditional Stage 1 (one binary x) -- s and Gamma at the prereg tiers on condx_A and condx_D, exact structure, four controls fail, map untouched, refusals at rc 3498"
exit 0
