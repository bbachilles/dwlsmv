*! dwlsmv -- A2b probe 2: RS7 route 1, dwlsmv_delta_cols() at the rotated point
*!
*! Executed under StataMP 18 (6 processors) on 2026-11-17, first run unchanged;
*! record in tools/probe_a2b_delta_route_stata18.txt. Drafted by a read-only agent.
*!
*! A probe, not a test (tools/); the suite runner ignores it. Roadmap unit A2b
*! (E-108), derivation @sec-rot-se RS7 route 1. Question: does
*! dwlsmv_delta_cols(), called on the EFA echelon parameter table with FIXED rows
*! selected (every loading, the echelon anchor zeros included, and the fixed factor
*! covariance rows), at the ROTATED point (Lambda and Phi in the values vector in
*! place of A and I), give the rotated Delta* = d sigma / d theta*'? It reports what
*! happened and chooses no convention. No src/ or tests/ change.
*!
*! Case E (fixtures/data/case_E.csv, p = 6 items of 3 categories, m = 2), the echelon
*! table written by hand as tests/unit/test_geomin.do does (lambda f2 y1 fixed 0, the
*! rest free, unitvariance on both factors) plus one line psi f1 f2 fixed 0, the row
*! dwlsmv_m11_efa_canon() emits, so the table carries the fixed covariance row.
*! n_u = 12 thresholds + 11 loadings = 23, n* = 12 + 12 + 1 = 25, D = 12 + 15 = 27.
*!
*! Construction:
*!   fit the unrotated echelon model (dwls_sample_stats, dwls_identify, dwls_fit);
*!   A from the table; (Lambda, Phi) = dwlsmv_geomin(A, 1e-4, 30, Phi, info);
*!   vals* = dwlsmv_parvals(theta-hat) with each loading row set to Lambda[row,col]
*!     (fixed zero row included) and each psi off-diagonal row to Phi[row,col];
*!   sel = threshold rows (table order), then loading rows in vec(Lambda) order
*!     (column-major), then psi off-diagonal rows (strict lower triangle of Phi,
*!     column-major);
*!   D1 = dwlsmv_delta_cols(vals*, par, map, sel).
*!
*! Comparators, written here and independent of dwlsmv_corrvec/dwlsmv_rystar:
*!   sigma(z), z = (tau in sel order, vec(Lambda), phi): thresholds placed at
*!     map.Coff[row] + col, then the correlations as the off-diagonal elements of
*!     Lambda Phi Lambda' (unit-variance indicators), pairs (i,j), i < j, i the outer
*!     loop. That is the map's order: dwlsmv_pos_rho() = D - P + (i-1)p - i(i-1)/2
*!     + (j-i). The agreement with dwlsmv_corrvec() at vals* is a diagnostic.
*!   D2  central finite differences of that sigma, step 1e-4 (delta_cols uses 1e-5).
*!   D3  analytic: d(L Phi L') = dL Phi L' + L Phi dL' + L dPhi L', entrywise
*!       dR_ij/dL_ka = [k=i] (L Phi)[j,a] + [k=j] (L Phi)[i,a],
*!       dR_ij/dphi_rc = L_ir L_jc + L_ic L_jr (r > c); threshold columns unit vectors.
*!
*! DECIDING CRITERION (fixed before running): max |D1 - D3| <= 1e-8, absolute.
*! Everything else is a diagnostic: max|D1-D2|, max|D2-D3|, dimensions, the numerical
*! rank of D1 (singular values, 1e-8 relative cut; RS2.1 predicts n_u = n* - m(m-1)
*! = 23), the threshold block of D1 against the identity, and D1 at the unrotated
*! point (Lambda = A, Phi = I) against the free columns of dwls_delta().
*!
*! Lines: A2BD <check> <values>, and one A2BD VERDICT line. Run from the repository
*! root (PowerShell):
*!     Start-Process -Wait -FilePath "C:\Program Files\Stata18\StataMP-64.exe" -ArgumentList "/e","do","tools/probe_a2b_delta_route.do"
*! (Git Bash: StataMP-64.exe //e do tools/probe_a2b_delta_route.do), then
*!     grep -E "^A2BD" probe_a2b_delta_route.log | diff - tools/probe_a2b_delta_route_stata18.txt
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 200
display "A2BDENV stata_version " c(stata_version) " processors " c(processors)

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_mg_check.mata          ///
             src/mata/dwlsmv_mg_sub.mata            ///
             src/mata/dwlsmv_mg_parts.mata          ///
             src/mata/dwlsmv_mg_sigma.mata          ///
             src/mata/dwlsmv_mg_dcols.mata          ///
             src/mata/dwlsmv_id_parse.mata          ///
             src/mata/dwlsmv_id_cx.mata             ///
             src/mata/dwlsmv_id_rank.mata           ///
             src/mata/dwlsmv_id_tie.mata            ///
             src/mata/dwlsmv_id_groups.mata         ///
             src/mata/dwlsmv_id_opname.mata         ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_statmap_q.mata         ///
             src/mata/dwlsmv_s2_guard.mata          ///
             src/mata/dwlsmv_pos_slope.mata         ///
             src/mata/dwlsmv_cx1_check.mata         ///
             src/mata/dwlsmv_cx1_fit_check.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
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
             src/mata/dwlsmv_nfac.mata              ///
             src/mata/dwlsmv_parvals.mata           ///
             src/mata/dwlsmv_matrices.mata          ///
             src/mata/dwlsmv_rystar.mata            ///
             src/mata/dwlsmv_corrvec.mata           ///
             src/mata/dwls_model.mata               ///
             src/mata/dwlsmv_delta_cols.mata        ///
             src/mata/dwls_delta.mata               ///
             src/mata/dwlsmv_canon_lines.mata       ///
             src/mata/dwlsmv_canon_find.mata        ///
             src/mata/dwlsmv_canon_status.mata      ///
             src/mata/dwlsmv_numrank.mata           ///
             src/mata/dwls_identify.mata            ///
             src/mata/dwlsmv_fdwls.mata             ///
             src/mata/dwlsmv_minimize.mata          ///
             src/mata/dwlsmv_signcanon.mata         ///
             src/mata/dwls_fit.mata                 ///
             src/mata/dwlsmv_asym.mata              ///
             src/mata/dwlsmv_acov.mata              ///
             src/mata/dwlsmv_traces.mata            ///
             src/mata/dwlsmv_scaleshift.mata        ///
             src/mata/dwls_inference.mata           ///
             src/mata/dwlsmv_geomin.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "probe_a2b_delta_route: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
mata: YE = st_data(., .)

mata:

// the echelon text of tests/unit/test_geomin.do's tg_spec(), plus the fixed psi row
string scalar a2bd_spec(real scalar p, real scalar m)
{
    string scalar s
    real scalar j, k, r

    s = "factor"
    for (k = 1; k <= m; k++) {
        s = s + " f" + strofreal(k)
    }
    s = s + char(10) + "ind"
    for (j = 1; j <= p; j++) {
        s = s + " y" + strofreal(j)
    }
    s = s + char(10)
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            s = s + "lambda f" + strofreal(k) + " y" + strofreal(j) +
                (j < k ? " fixed 0" : " free") + char(10)
        }
    }
    for (k = 1; k <= m - 1; k++) {
        for (r = k + 1; r <= m; r++) {
            s = s + "psi f" + strofreal(k) + " f" + strofreal(r) + " fixed 0" + char(10)
        }
    }
    for (k = 1; k <= m; k++) {
        s = s + "metric unitvariance f" + strofreal(k) + char(10)
    }
    return(s)
}

// Independent sigma. z = (tau [nT], vec(L) [p*m], phi [m(m-1)/2, strict lower, by
// column]); tpos = position of each tau in sigma; D = rows of sigma = nT + p(p-1)/2.
real colvector a2bd_sig(real colvector z, real colvector tpos, real scalar p,
                        real scalar m, real scalar D)
{
    real matrix L, Phi
    real colvector s
    real scalar nT, i, j, k, c, r, idx

    nT = rows(tpos)
    L  = J(p, m, .)
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            L[j, k] = z[nT + (k - 1) * p + j]
        }
    }
    Phi = I(m)
    idx = nT + p * m
    for (c = 1; c <= m - 1; c++) {
        for (r = c + 1; r <= m; r++) {
            idx++
            Phi[r, c] = z[idx]
            Phi[c, r] = z[idx]
        }
    }
    s = J(D, 1, .)
    for (k = 1; k <= nT; k++) {
        s[tpos[k]] = z[k]
    }
    idx = 0
    for (i = 1; i <= p - 1; i++) {
        for (j = i + 1; j <= p; j++) {
            idx++
            s[nT + idx] = L[i, .] * Phi * L[j, .]'
        }
    }
    return(s)
}

// Analytic Jacobian of a2bd_sig, columns in z's order.
real matrix a2bd_jac(real colvector z, real colvector tpos, real scalar p,
                     real scalar m, real scalar D)
{
    real matrix L, Phi, LP, Jm
    real scalar nT, i, j, k, c, r, a, idx, q, row

    nT = rows(tpos)
    L  = J(p, m, .)
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            L[j, k] = z[nT + (k - 1) * p + j]
        }
    }
    Phi = I(m)
    idx = nT + p * m
    for (c = 1; c <= m - 1; c++) {
        for (r = c + 1; r <= m; r++) {
            idx++
            Phi[r, c] = z[idx]
            Phi[c, r] = z[idx]
        }
    }
    LP = L * Phi
    Jm = J(D, rows(z), 0)
    for (k = 1; k <= nT; k++) {
        Jm[tpos[k], k] = 1
    }
    idx = 0
    for (i = 1; i <= p - 1; i++) {
        for (j = i + 1; j <= p; j++) {
            idx++
            row = nT + idx
            for (a = 1; a <= m; a++) {
                Jm[row, nT + (a - 1) * p + i] = LP[j, a]
                Jm[row, nT + (a - 1) * p + j] = LP[i, a]
            }
            q = 0
            for (c = 1; c <= m - 1; c++) {
                for (r = c + 1; r <= m; r++) {
                    q++
                    Jm[row, nT + p * m + q] = L[i, r] * L[j, c] + L[i, c] * L[j, r]
                }
            }
        }
    }
    return(Jm)
}

// numerical rank, singular values against a relative cut
real scalar a2bd_nrank(real matrix X, real scalar cut)
{
    real colvector sv
    sv = svdsv(X)
    return(sum(sv :> cut * max(sv)))
}

real scalar a2bd_run(real matrix Y)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit    scalar fit
    struct dwlsmv_statmap scalar map
    real matrix A, Psi, B, L, Phi, D1, D2, D3, D1u, Dfree, Dpl
    real rowvector info
    real colvector K, vals, vstar, sel, tauSel, tpos, tvals, z, zp, zm, svs, idx, sg1, sg2
    real scalar p, m, j, k, i, c, r, nT, nL, nC, n, Dd, h, eps, crit, dmax13, dmax12, dmax23
    real scalar dthr, dstr, dcons, ncons, nfix, dsig, pl, rk, dtau, nsel

    p = cols(Y)
    m = 2
    K = J(p, 1, .)
    for (j = 1; j <= p; j++) {
        K[j] = max(Y[., j]) + 1
    }
    st  = dwls_sample_stats(Y, K)
    map = st.map
    par = dwls_identify(a2bd_spec(p, m), map)
    fit = dwls_fit(st, par)
    printf("A2BD fit converged %g npar %g nall %g D %g T %g P %g\n", fit.converged,
           par.npar, par.nall, map.D, map.T, map.P)
    if (fit.converged != 1) {
        errprintf("A2BD: the unrotated fit did not converge\n")
        _error("A2BD: the unrotated fit did not converge")
    }
    if (map.D != map.T + map.P) {
        _error("A2BD: the map has slope rows")
    }

    vals = dwlsmv_parvals(fit.theta, par)
    dwlsmv_matrices(vals, par, p, m, A, Psi, B)
    eps = 1e-4
    L   = dwlsmv_geomin(A, eps, 30, Phi, info)
    printf("A2BD rotation converged %g of %g starts, Q %12.8f, gpnorm %9.2e\n",
           info[2], info[1], info[5], info[6])

    // the selection, table order for thresholds, vec(Lambda) order, strict lower of Phi
    tauSel = selectindex(par.ptype :== 4)
    nT = rows(tauSel)
    sel = tauSel
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            idx = selectindex((par.ptype :== 1) :& (par.row :== j) :& (par.col :== k))
            if (rows(idx) != 1) {
                _error("A2BD: a loading row is not unique in the table")
            }
            sel = sel \ idx
        }
    }
    nL = p * m
    nC = 0
    for (c = 1; c <= m - 1; c++) {
        for (r = c + 1; r <= m; r++) {
            idx = selectindex((par.ptype :== 2) :& (((par.row :== r) :& (par.col :== c)) :| ((par.row :== c) :& (par.col :== r))))
            if (rows(idx) != 1) {
                _error("A2BD: a factor covariance row is not unique in the table")
            }
            sel = sel \ idx
            nC++
        }
    }
    n = rows(sel)
    printf("A2BD dims sel %g (thresholds %g, loadings %g, covariances %g); D rows %g\n",
           n, nT, nL, nC, map.D)
    nfix = 0
    for (k = 1; k <= n; k++) {
        nfix = nfix + (par.free[sel[k]] == 0)
    }
    printf("A2BD sel fixed rows %g free rows %g (expect fixed = m(m-1) = %g)\n",
           nfix, n - nfix, m * (m - 1))

    // vals*: the rotated point
    vstar = vals
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1) {
            vstar[i] = L[par.row[i], par.col[i]]
        }
        else if (par.ptype[i] == 2 & par.row[i] != par.col[i]) {
            vstar[i] = Phi[par.row[i], par.col[i]]
        }
    }

    tpos = J(nT, 1, .)
    dtau = 0
    for (k = 1; k <= nT; k++) {
        tpos[k] = map.Coff[par.row[tauSel[k]]] + par.col[tauSel[k]]
        dtau = max((dtau, abs(tpos[k] - dwlsmv_pos_tau(map, par.row[tauSel[k]], par.col[tauSel[k]]))))
    }
    printf("A2BD tau positions Coff+col vs dwlsmv_pos_tau max diff %g\n", dtau)

    // the point z* for the comparators
    z = vstar[tauSel] \ vec(L)
    for (c = 1; c <= m - 1; c++) {
        for (r = c + 1; r <= m; r++) {
            z = z \ Phi[r, c]
        }
    }
    Dd = map.D

    // diagnostic: independent sigma against dwlsmv_corrvec at vals*
    sg1 = a2bd_sig(z, tpos, p, m, Dd)
    sg2 = dwlsmv_corrvec(vstar, par, map)
    dsig = max(abs(sg1[(nT + 1)::Dd] - sg2))
    printf("A2BD sigma independent vs dwlsmv_corrvec at vals*, max abs diff %9.2e (diagnostic)\n", dsig)

    // D1
    D1 = dwlsmv_delta_cols(vstar, par, map, sel)
    printf("A2BD D1 dims %g x %g\n", rows(D1), cols(D1))

    // D2 central finite differences, step 1e-4
    h  = 1e-4
    D2 = J(Dd, n, .)
    for (k = 1; k <= n; k++) {
        zp = z
        zm = z
        zp[k] = z[k] + h
        zm[k] = z[k] - h
        D2[., k] = (a2bd_sig(zp, tpos, p, m, Dd) - a2bd_sig(zm, tpos, p, m, Dd)) / (2 * h)
    }
    // D3 analytic
    D3 = a2bd_jac(z, tpos, p, m, Dd)

    dmax13 = max(abs(D1 - D3))
    dmax12 = max(abs(D1 - D2))
    dmax23 = max(abs(D2 - D3))
    printf("A2BD max|D1-D3| %9.3e (deciding)\n", dmax13)
    printf("A2BD max|D1-D2| %9.3e (diagnostic, FD step %g)\n", dmax12, h)
    printf("A2BD max|D2-D3| %9.3e (diagnostic)\n", dmax23)
    printf("A2BD max|D1|    %9.3e\n", max(abs(D1)))

    // by block
    dthr = max(abs(D1[tpos, 1..nT] - I(nT)))
    dstr = max(abs(D1[tpos, (nT + 1)..n]))
    printf("A2BD D1 threshold block vs identity max abs %9.2e; threshold rows of the structural columns max abs %9.2e\n",
           dthr, dstr)
    printf("A2BD D1 non-threshold rows of the threshold columns max abs %9.2e\n",
           max(abs(D1[(nT + 1)::Dd, 1..nT])))
    printf("A2BD D1-D3 on the loading columns %9.3e, on the covariance columns %9.3e\n",
           max(abs(D1[., (nT + 1)..(nT + nL)] - D3[., (nT + 1)..(nT + nL)])),
           max(abs(D1[., (nT + nL + 1)..n] - D3[., (nT + nL + 1)..n])))

    // rank
    rk = a2bd_nrank(D1, 1e-8)
    svs = svdsv(D1)
    printf("A2BD rank D1 %g of %g columns (RS2.1 predicts n_u = %g); sigma_min/sigma_max %9.3e\n",
           rk, n, n - m * (m - 1), min(svs) / max(svs))
    printf("A2BD rank D3 %g, rank D2 %g\n", a2bd_nrank(D3, 1e-8), a2bd_nrank(D2, 1e-8))

    // D1 at the unrotated point against dwls_delta's free columns
    D1u   = dwlsmv_delta_cols(vals, par, map, sel)
    Dfree = dwls_delta(fit.theta, par, map)
    dcons = 0
    ncons = 0
    for (k = 1; k <= n; k++) {
        if (par.free[sel[k]]) {
            ncons++
            dcons = max((dcons, max(abs(D1u[., k] - Dfree[., par.fidx[sel[k]]]))))
        }
    }
    printf("A2BD unrotated point: %g free shared columns, D1u vs dwls_delta max abs diff %9.2e; dwls_delta dims %g x %g\n",
           ncons, dcons, rows(Dfree), cols(Dfree))
    printf("A2BD unrotated point: fixed columns of D1u (anchor zero and psi) max abs %9.3e\n",
           max(abs(select(D1u, (par.free[sel] :== 0)'))))
    printf("A2BD rank D1u (sel columns) %g; rank Dfree %g of %g\n", a2bd_nrank(D1u, 1e-8),
           a2bd_nrank(Dfree, 1e-8), cols(Dfree))

    crit = 1e-8
    printf("A2BD VERDICT %s: max|D1-D3| %9.3e against criterion %g\n",
           (dmax13 <= crit ? "pass" : "fail"), dmax13, crit)
    return(dmax13 <= crit)
}

end

capture noisily mata: (void) a2bd_run(YE)
display as result "A2BD run rc " _rc
