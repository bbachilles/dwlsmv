*! sens_t3_setup.mata -- the fits and inputs shared by test_sens_t3.do and
*! tools/export_sens_t3_in.do (the reference values' provenance)
*!
*! Executed under StataMP 18 on 2026-10-16: tests/unit/test_sens_t3.do green at six processors.
*!
*! Design E (fixtures/data/mg_E.csv, G = 2; consistency-only data for the
*! check 6 investigation) at the metric and scalar levels, built as
*! test_sens.do builds its fit, then dwls_inference(). The step is
*! theta_K - theta_{K-1}, with theta_{K-1} from maxiter K - 2 (optimize() takes
*! m + 1 steps at maxiter m; test_sens.do s1).
*!
*! The synthetic measured inputs (the "full" case), fixed before any run:
*!     ds_k    = 1e-13 sin(k)
*!     eps_k   = 3e-16 cos(k)
*!     dOm_ij  = 4u |Omega_ij| sin(i + j)        (symmetric)
*! sized like the grid's measured Stage 1 differences (ulp scale); nothing
*! here is a measurement.
*!
*! Needs e(K) and e(canon_configural) from a configural dwlsmv fit in memory,
*! and the data loaded. Loaded with `do`; defines no test of its own.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

mata:

void st3_put(string scalar nm, transmorphic x)
{
    pointer() scalar p

    p = crexternal(nm)
    if (p == NULL) p = findexternal(nm)
    *p = x
}

void st3_setup()
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_statmap vector maps
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    string scalar cc, lev
    string rowvector levs
    real matrix K, Y
    real colvector gv, s, w, start, thu, thK1
    real scalar G, g, N, i, fu, iu, fK1, iK1, L

    K    = st_matrix("e(K)")
    cc   = st_global("e(canon_configural)")
    G    = rows(K)
    gv   = uniqrows(st_data(., "g"))
    levs = ("metric", "scalar")
    for (L = 1; L <= cols(levs); L++) {
        lev  = levs[L]
        sts  = dwlsmv_stage1(G, 1)
        maps = dwlsmv_statmap(G, 1)
        for (g = 1; g <= G; g++) {
            Y       = select(st_data(., ("y1", "y2", "y3", "y4", "y5")), st_data(., "g") :== gv[g])
            sts[g]  = dwls_sample_stats(Y, K[g, .]')
            maps[g] = sts[g].map
        }
        par = dwls_identify(dwlsmv_levelgen(cc, K, lev), maps, sts)
        fit = dwlsmv_tsu_fit(sts, par)          // unpolished: the diagnostics model optimize()'s stop (sens_unpolished.mata)
        dwls_inference(fit)
        dwlsmv_fit_stack(sts, s, w, N, maps)

        // dwls_fit's sequence again: ULS from the table's starts, then DWLS to K - 1
        start = J(par.npar, 1, .)
        for (i = 1; i <= par.nall; i++) {
            if (par.free[i]) start[par.fidx[i]] = par.value[i]
        }
        (void) dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, start, "d2", "nr", 300, thu, fu, iu)
        (void) dwlsmv_minimize(s, w, par, maps, thu, "d2", "nr", fit.iterations - 2, thK1, fK1, iK1)

        st3_put("ST3_F_" + lev, fit)
        st3_put("ST3_STEP_" + lev, fit.theta - thK1)
    }
}

void st3_pert(struct dwlsmv_fit scalar f, real colvector ds, real colvector eps, real matrix dOm)
{
    real matrix Om
    real scalar Dn, i, j

    Om  = dwlsmv_omega(f.stage1)
    Dn  = rows(Om)
    ds  = 1e-13 * sin(1::Dn)
    eps = 3e-16 * cos(1::Dn)
    dOm = J(Dn, Dn, .)
    for (i = 1; i <= Dn; i++) {
        for (j = 1; j <= Dn; j++) {
            dOm[i, j] = 4 * (epsilon(1) / 2) * abs(Om[i, j]) * sin(i + j)
        }
    }
}

// A case's fit, step and measured inputs; returns 1 if the case gives (ds, eps, dOm).
// Cases: "metric_full", "scalar_full", "metric_none", "metric_bigeps" (full's ds and dOm,
// with eps_k = 1e-6 cos(k): large enough that ||eps||_inf's term in beta_2's gbar shows
// at the reference's 1e-12; added after mutant M9 survived the first three cases)
real scalar st3_inputs(string scalar cs, struct dwlsmv_fit scalar f, real colvector step,
                       real colvector ds, real colvector eps, real matrix dOm)
{
    string rowvector t

    t    = tokens(subinstr(cs, "_", " "))
    f    = *findexternal("ST3_F_" + t[1])
    step = *findexternal("ST3_STEP_" + t[1])
    if (t[2] == "none") {
        return(0)
    }
    st3_pert(f, ds, eps, dOm)
    if (t[2] == "bigeps") {
        eps = 1e-6 * cos(1::rows(eps))
    }
    return(1)
}

real rowvector st3_case(string scalar cs, real scalar epsV, real scalar rhoT)
{
    struct dwlsmv_fit scalar f
    real colvector step, ds, eps
    real matrix dOm

    if (st3_inputs(cs, f, step, ds, eps, dOm)) {
        return(dwlsmv_sens_t3(f, step, epsV, rhoT, ds, eps, dOm))
    }
    return(dwlsmv_sens_t3(f, step, epsV, rhoT))
}

// The same case on the exact path (dwlsmv_sens_t3_x, Result MS6), the path
// test_sens_t3's r1x and r2x compare (user, plan step 5, 2026-10-19).
real rowvector st3_case_x(string scalar cs, real scalar epsV, real scalar rhoT)
{
    struct dwlsmv_fit scalar f
    real colvector step, ds, eps
    real matrix dOm

    if (st3_inputs(cs, f, step, ds, eps, dOm)) {
        return(dwlsmv_sens_t3_x(f, step, epsV, rhoT, ds, eps, dOm))
    }
    return(dwlsmv_sens_t3_x(f, step, epsV, rhoT))
}

end
