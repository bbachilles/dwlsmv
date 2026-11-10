*! Test-side helper: the fit as dwls_fit() built it BEFORE the post-convergence polish (adopted
*! 2026-10-23), i.e. ULS, then DWLS by dwlsmv_minimize(), then dwlsmv_signcanon(), with no polish.
*!
*! The sensitivity diagnostics (dwlsmv_sens, dwlsmv_sens_x, dwlsmv_sens_t3; diagnostic only, no production
*! caller) model optimize()'s stopping point: S3 takes K = f.iterations and optimize()'s measured last
*! step. Under plan P8 (a) they keep describing optimize(), so their tests (test_sens.do, test_sens_hess.do,
*! test_sens_t3.do via sens_t3_setup.mata) build their fits with this helper (user, 2026-10-23, option (a);
*! E-014). Their pinned references were exported from exactly this fit and are
*! unchanged.
*!
*!     struct dwlsmv_fit scalar dwlsmv_tsu_fit(struct dwlsmv_stage1 vector sts,
*!                                             struct dwlsmv_parmap scalar par)
*!
*! A verbatim copy of dwls_fit() and dwlsmv_fit_from() at 08b5b99^ (src/mata/dwls_fit.mata), the two
*! merged, the error checks kept. NOT PRODUCTION CODE.

mata:

struct dwlsmv_fit scalar dwlsmv_tsu_fit(struct dwlsmv_stage1 vector sts,
                                        struct dwlsmv_parmap  scalar par)
{
    struct dwlsmv_fit scalar f
    struct dwlsmv_statmap vector maps
    real colvector start, thu, s, w, th
    real rowvector gr
    real matrix H
    real scalar i, fu, iu, N, G, fm, it, conv

    if (length(sts) < 1) {
        _error("dwlsmv_tsu_fit: no Stage 1 result was passed")
    }
    if (length(sts) > 1) {
        dwlsmv_mg_check(par, length(sts), "dwlsmv_tsu_fit")
    }
    G = length(sts)
    dwlsmv_fit_stack(sts, s, w, N, maps)

    start = J(par.npar, 1, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) {
            start[par.fidx[i]] = par.value[i]
        }
    }
    if (missing(start)) {
        _error("dwlsmv_tsu_fit: fidx does not number the free rows 1..npar")
    }

    (void) dwlsmv_minimize(s, J(rows(s), 1, 1), par, maps, start,
                           "d2", "nr", 300, thu, fu, iu)
    conv = dwlsmv_minimize(s, w, par, maps, thu,
                           "d2", "nr", 300, th, fm, it)

    th = dwlsmv_signcanon(th, par)
    dwlsmv_fdwls(0, th', s, w, par, maps, fm, gr, H)

    f.stage1     = sts
    f.par        = par
    f.theta      = th
    f.Delta      = dwls_delta(th, par, maps)
    f.fmin       = fm
    f.converged  = conv
    f.iterations = it
    f.T_raw      = (N - G) * fm
    f.D_stat     = rows(s) - par.npar

    return(f)
}

end
