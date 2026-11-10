*! dwlsmv -- Gamma-hat rank study, design B: threshold aliasing under a
*! CORRECTLY SPECIFIED model (Module 10; E-054 "Rank-deficient
*! Gamma"). Companion to tools/gamma_rank_study.do, whose aliased conditions
*! were misspecified. A diagnostic tool, not a pass/fail test.
*!
*! Design (chosen 2026-10-05 after seeing design A's result, before any
*! design-B data): one factor, p = 4, K = 3, marker on y1, N = 500. y1 and y2
*! load .97 (rho12 = .94) and share thresholds (-2, .5); y3, y4 load .7 with
*! (-.5, .7). Conditions keep a draw only when exactly m rows are discordant
*! on the lowest category, (y1 == 0) != (y2 == 0): ctrl (no selection),
*! alias0 (m = 0, exact aliasing), alias1, alias5. The model is true for every
*! draw; the aliasing is selected, not built. CAVEAT: the Monte Carlo SD is
*! then conditional on the selection event, which by itself may differ from
*! the sandwich SE; read alias* against that, and ctrl as the unconditional
*! check. The est/se columns are relative to the marker metric (lambda_y1 = 1).
*!
*! FAILURES ARE RECORDED, NEVER EXCLUDED: each stage runs through _stata(),
*! and a failure writes status "error" with its stage and rc.
*!
*! Usage, from the repository root through a one-line wrapper in %TEMP%:
*!     do tools/gamma_rank_study2.do [R]      (default R = 500)
*! Output: c(tmpdir)/gamma_rank_study2.csv, one row per replicate; the draws
*! taken per condition are displayed.
*!
*! RESULT, 2026-10-05, StataMP 18 (serial [dev], 6 processors), R = 500:
*! 2000 of 2000 replicates ok, 171 s. Reading in E-054,
*! "Rank-deficient Gamma" (design B).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
args nrep
if "`nrep'" == "" local nrep = 500
version 16
clear all
set more off
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
quietly do src/mata/dwlsmv_types.mata
foreach f of local srcs {
    quietly do "src/mata/`f'"
}
local out "`c(tmpdir)'/gamma_rank_study2.csv"

mata:
string scalar gr_spec()
{
    return("factor f" + char(10) + "ind y1 y2 y3 y4" + char(10) +
        "lambda f y2 free" + char(10) + "lambda f y3 free" + char(10) +
        "lambda f y4 free" + char(10) + "psi f f free" + char(10) +
        "metric marker f y1" + char(10))
}

string scalar gr_n(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

// one draw from the correctly specified one-factor model: y1, y2 load .97
// with thresholds (-2, .5); y3, y4 load .7 with thresholds (-.5, .7)
real matrix gr_draw(real scalar N)
{
    real colvector eta, z
    real rowvector lam, t1, t2
    real matrix Y
    real scalar j

    lam = (.97, .97, .7, .7)
    t1 = (-2, -2, -.5, -.5)
    t2 = (.5, .5, .7, .7)
    eta = rnormal(N, 1, 0, 1)
    Y = J(N, 4, .)
    for (j = 1; j <= 4; j++) {
        z = lam[j] * eta + sqrt(1 - lam[j]^2) * rnormal(N, 1, 0, 1)
        Y[., j] = (z :> t1[j]) + (z :> t2[j])
    }
    return(Y)
}

// m = . : the first draw (control); m >= 0 : redraw until exactly m rows
// have (y1 == 0) != (y2 == 0). GR_DRAWS counts the draws taken.
real matrix gr_data(real scalar N, real scalar m)
{
    external real scalar GR_DRAWS
    real matrix Y

    Y = gr_draw(N)
    GR_DRAWS = GR_DRAWS + 1
    if (m < .) {
        while (sum((Y[., 1] :== 0) :!= (Y[., 2] :== 0)) != m) {
            Y = gr_draw(N)
            GR_DRAWS = GR_DRAWS + 1
        }
    }
    return(Y)
}

real scalar _gr_try1(real matrix Y, real colvector K,
                     struct dwlsmv_stage1 scalar st)
{
    external real matrix GR_Y
    external real colvector GR_K
    external struct dwlsmv_stage1 scalar GR_ST
    real scalar rc

    GR_Y = Y
    GR_K = K
    rc = _stata("mata: GR_ST = dwls_sample_stats(GR_Y, GR_K)", 1)
    if (!rc) {
        st = GR_ST
    }
    return(rc)
}

real scalar _gr_try2(struct dwlsmv_stage1 scalar st,
                     struct dwlsmv_parmap scalar par,
                     struct dwlsmv_fit scalar fit)
{
    external struct dwlsmv_stage1 scalar GR_ST
    external struct dwlsmv_parmap scalar GR_PAR
    external struct dwlsmv_fit scalar GR_FIT
    real scalar rc

    GR_ST = st
    rc = _stata("mata: GR_PAR = dwls_identify(gr_spec(), GR_ST.map)", 1)
    if (!rc) {
        rc = _stata("mata: GR_FIT = dwls_fit(GR_ST, GR_PAR)", 1)
    }
    if (!rc) {
        rc = _stata("mata: dwls_inference(GR_FIT)", 1)
    }
    if (!rc) {
        par = GR_PAR
        fit = GR_FIT
    }
    return(rc)
}

void gr_run(real scalar nrep, string scalar outf)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    string rowvector cn, row
    string scalar hdr, stage
    real rowvector mv
    real matrix Y, G, B, eV
    real colvector K, sv, svb, s, se, WD
    real scalar fh, c, r, i, rc, N, D, nsv, nest

    external real scalar GR_DRAWS

    cn = ("ctrl", "alias0", "alias1", "alias5")
    mv = (., 0, 1, 5)
    K = J(4, 1, 3)
    N = 500
    nsv = 14
    nest = 12
    if (fileexists(outf)) {
        unlink(outf)
    }
    fh = fopen(outf, "w")
    hdr = "condition,replicate,seed,status,stage,rc,converged,T3,T3_df"
    for (i = 1; i <= 7; i++) {
        hdr = hdr + ",g" + strofreal(i) + "f,g" + strofreal(i) + "v"
    }
    hdr = hdr + ",g_rstar,centring_z"
    for (i = 1; i <= nsv; i++) {
        hdr = hdr + ",sv" + strofreal(i)
    }
    for (i = 1; i <= 4; i++) {
        hdr = hdr + ",b" + strofreal(i) + "v"
    }
    hdr = hdr + ",bread_cond"
    for (i = 1; i <= nest; i++) {
        hdr = hdr + ",est" + strofreal(i) + ",se" + strofreal(i)
    }
    fput(fh, hdr)
    for (c = 1; c <= cols(cn); c++) {
        GR_DRAWS = 0
        for (r = 1; r <= nrep; r++) {
            rseed(20261005 + 100000 * c + r)
            row = J(1, 9 + 16 + nsv + 5 + 2 * nest, ".")
            row[1..3] = (cn[c], strofreal(r),
                strofreal(20261005 + 100000 * c + r, "%12.0f"))
            Y = gr_data(N, mv[c])
            stage = "stage1"
            rc = _gr_try1(Y, K, st)
            if (!rc) {
                G = dwlsmv_gammadiag(st, Y)
                for (i = 1; i <= 7; i++) {
                    row[8 + 2 * i] = gr_n(G[i, 1])
                    row[9 + 2 * i] = gr_n(G[i, 2])
                }
                row[24] = gr_n(G[7, 3])
                row[25] = gr_n(G[5, 2])
                D = st.map.D
                s = sqrt(diagonal(st.Gamma))
                sv = sort(abs(symeigenvalues(st.Gamma :/ (s * s'))'), -1)
                sv = sv / sv[1]
                for (i = 1; i <= min((nsv, D)); i++) {
                    row[25 + i] = gr_n(sv[i])
                }
                stage = "identify/fit/inference"
                rc = _gr_try2(st, par, fit)
            }
            if (!rc) {
                B = dwlsmv_breaddiag(fit.Delta, st.Gamma)
                for (i = 1; i <= 4; i++) {
                    row[25 + nsv + i] = gr_n(B[i, 2])
                }
                WD = 1 :/ diagonal(st.Gamma)
                svb = svdsv(sqrt(WD) :* fit.Delta)
                row[30 + nsv] = gr_n((max(svb) / min(svb))^2)
                row[7..9] = (gr_n(fit.converged), gr_n(fit.T3), gr_n(fit.T3_df))
                eV = dwlsmv_ev(fit.acov, st.N)
                se = sqrt(diagonal(eV))
                for (i = 1; i <= min((nest, rows(fit.theta))); i++) {
                    row[29 + nsv + 2 * i] = gr_n(fit.theta[i])
                    row[30 + nsv + 2 * i] = gr_n(se[i])
                }
                row[4] = (fit.converged == 1 ? "ok" : "nonconverged")
            }
            else {
                row[4..6] = ("error", stage, strofreal(rc))
            }
            hdr = row[1]
            for (i = 2; i <= cols(row); i++) {
                hdr = hdr + "," + row[i]
            }
            fput(fh, hdr)
        }
        printf("gamma_rank_study2: %s, %g draws for %g replicates\n",
            cn[c], GR_DRAWS, nrep)
    }
    fclose(fh)
}
end

timer clear 1
timer on 1
mata: gr_run(`nrep', "`out'")
timer off 1
timer list 1
display "gamma_rank_study2: wrote `out'"
