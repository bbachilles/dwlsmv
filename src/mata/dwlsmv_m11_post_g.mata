*! dwlsmv -- Module 11: what an invariance() fit posts (group surface, step 4)
*!
*!     void dwlsmv_m11_post_g(canon, indlist, gvar, gvalslist, semnames,
*!                            bn, vn, cn, kn)
*!
*! Build step 4 (2026-10-12). Reads the step-3 externals (dwlsmv_m11_gpar,
*! _gfit, _gtext, _st, _K) and fills the Stata matrices named bn, vn, cn, kn
*! and locals for dwlsmv.ado to post.
*!
*! e(b) (spec D5, multi-group posting convention S4-1 (a), order S4-2 (a)):
*! the full group-expanded table at every level, in blocks: the parameter
*! map's tau, lambda, beta, psi rows in its order; then nu, omega, alpha, each
*! over every group (group 1 first; indicators, or factors, in order), a map
*! row where the level has one, else the absent default posted fixed (nu 0,
*! omega 1, alpha 0, G4). A posted fixed default is not a table row: e(rank)
*! is the map's free count, e(canon) the table.
*!
*! Names (spec 4.8, decided 2026-10-10), tag = v + "bn." + gvar for group 1,
*! v + "." + gvar after, v the group's value:
*!   tau    /y1:tag#c.cut1      lambda y1:tag#c.F       beta   F2:tag#c.F1
*!   psi    /var(F)#tag, /cov(F,G)#tag (e. prefix if endogenous)
*!   nu     y1:tag              omega  /scale(y1)#tag   alpha  /mean(F)#tag
*! D6 at G > 1: every built name other than tau and omega is in sem's stripe
*! (ginvariant(none)) and every sem name is built, sem's indicator residual
*! variances /var(e.y)#tag set aside (Theta is derived, as at G = 1); a
*! mismatch is an internal error.
*!
*! e(V) = ACOV / (N - G) embedded, exact zeros in fixed rows and columns; tied
*! members share their theta index, so their rows repeat (spec 4.8 G2).
*! e(Cns): one row per posted fixed value. e(K): G x p, columns named by
*! indicator, rows by group value, asserted against the thresholds in e(b)
*! (K - 1 per group and indicator; an internal error if not).
*! e(bread_ratio) at G > 1 (O-T10 (a)): Module 10's bread on the stacked
*! Delta with W = blockdiag(pi_g / diag(Gamma_g)), pi_g = N_g / N.
*! Locals: N, rank, nfix, T3, df, converged, iterations, bratio. The configural text (the generator's own
*! configural output) goes to the external dwlsmv_m11_gconf for
*! dwlsmv_m11_post_g_e().
*!
*! Executed under StataMP 18 on 2026-10-12: tests/unit/test_invariance_post.do.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_post_g(string scalar canon, string scalar indlist, string scalar gvar,
                       string scalar gvalslist, string scalar semnames,
                       string scalar bn, string scalar vn, string scalar cn,
                       string scalar kn)
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    struct dwlsmv_stage1 vector sts
    string rowvector ind, fac, gv, sem, built, w
    string colvector L, tag, nm
    string matrix stp
    real matrix K, P, V, Cns, Bd
    real colvector endo, bvec, gd
    real scalar G, i, k, g, j, f, n, N, nfix, c
    string scalar a, z, conf
    pointer() scalar p

    par = *findexternal("dwlsmv_m11_gpar")
    fit = *findexternal("dwlsmv_m11_gfit")
    sts = *findexternal("dwlsmv_m11_st")
    K   = *findexternal("dwlsmv_m11_K")
    ind = tokens(indlist)
    gv  = tokens(gvalslist)
    G   = rows(sts)
    L   = dwlsmv_canon_lines(subinstr(canon, ";", char(10)))
    fac = J(1, 0, "")
    for (i = 1; i <= rows(L); i++) {
        w = tokens(L[i])
        if (w[1] == "factor") {
            fac = w[2..cols(w)]
            break
        }
    }
    tag = J(G, 1, "")
    for (g = 1; g <= G; g++) {
        tag[g] = gv[g] + (g == 1 ? "bn." : ".") + gvar
    }

    // ---- the posted rows: (ptype, grp, row, col, free, value, fidx) --------
    N = 0
    for (g = 1; g <= G; g++) {
        N = N + sts[g].N
    }
    dwlsmv_m11_gbv(par, fit, N, G, cols(ind), cols(fac), P, bvec, V)
    n = rows(P)

    // ---- names --------------------------------------------------------------
    endo = J(cols(fac), 1, 0)
    for (i = 1; i <= n; i++) {
        if (P[i, 1] == 3) {
            endo[P[i, 3]] = 1
        }
    }
    stp = J(n, 2, "")
    for (i = 1; i <= n; i++) {
        g = P[i, 2]
        if (P[i, 1] == 4) {
            stp[i, .] = ("/" + ind[P[i, 3]], tag[g] + "#c.cut" + strofreal(P[i, 4]))
        }
        else if (P[i, 1] == 1) {
            stp[i, .] = (ind[P[i, 3]], tag[g] + "#c." + fac[P[i, 4]])
        }
        else if (P[i, 1] == 3) {
            stp[i, .] = (fac[P[i, 3]], tag[g] + "#c." + fac[P[i, 4]])
        }
        else if (P[i, 1] == 2) {
            a = (endo[P[i, 4]] ? "e." : "") + fac[P[i, 4]]
            z = (endo[P[i, 3]] ? "e." : "") + fac[P[i, 3]]
            stp[i, .] = ("/", (P[i, 3] == P[i, 4] ? "var(" + z + ")" : "cov(" + a + "," + z + ")") + "#" + tag[g])
        }
        else if (P[i, 1] == 5) {
            stp[i, .] = (ind[P[i, 3]], tag[g])
        }
        else if (P[i, 1] == 6) {
            stp[i, .] = ("/", "scale(" + ind[P[i, 3]] + ")#" + tag[g])
        }
        else {
            stp[i, .] = ("/", "mean(" + fac[P[i, 3]] + ")#" + tag[g])
        }
    }
    nm = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        nm[i] = (stp[i, 1] == "/" ? "/" + stp[i, 2] : stp[i, 1] + ":" + stp[i, 2])
    }

    // ---- D6 at G > 1 ---------------------------------------------------------
    built = select(nm, P[., 1] :!= 4 :& P[., 1] :!= 6)'
    w     = tokens(semnames)
    sem   = J(1, 0, "")
    for (i = 1; i <= cols(w); i++) {
        c = 0
        if (substr(w[i], 1, 7) == "/var(e.") {
            for (j = 1; j <= cols(ind); j++) {
                if (substr(w[i], 8, strlen(ind[j]) + 2) == ind[j] + ")#") {
                    c = 1
                }
            }
        }
        if (!c) {
            sem = sem, w[i]
        }
    }
    for (i = 1; i <= cols(built); i++) {
        if (!anyof(sem, built[i])) {
            errprintf("dwlsmv: internal error: built name %s is not in sem's stripe.\n", built[i])
            errprintf("    built: %s\n    sem:   %s\n", invtokens(built), invtokens(sem))
            _error("dwlsmv: internal error, parameter names differ from sem's")
        }
    }
    for (i = 1; i <= cols(sem); i++) {
        if (!anyof(built, sem[i])) {
            errprintf("dwlsmv: internal error: sem's name %s was not built.\n", sem[i])
            errprintf("    built: %s\n    sem:   %s\n", invtokens(built), invtokens(sem))
            _error("dwlsmv: internal error, parameter names differ from sem's")
        }
    }

    // ---- e(K), asserted against e(b)'s thresholds ----------------------------
    for (g = 1; g <= G; g++) {
        for (j = 1; j <= cols(ind); j++) {
            if (sum(P[., 1] :== 4 :& P[., 2] :== g :& P[., 3] :== j) != K[g, j] - 1) {
                errprintf("dwlsmv: internal error: %s in group %f has K = %f but %f thresholds in e(b).\n",
                          ind[j], g, K[g, j], sum(P[., 1] :== 4 :& P[., 2] :== g :& P[., 3] :== j))
                _error("dwlsmv: internal error, e(K) disagrees with e(b)")
            }
        }
    }
    st_matrix(kn, K)
    st_matrixcolstripe(kn, (J(cols(ind), 1, ""), ind'))
    st_matrixrowstripe(kn, (J(G, 1, ""), gv'))

    // ---- e(b), e(V) (built above by dwlsmv_m11_gbv), e(Cns) ----------------------
    nfix = n - sum(P[., 5])
    Cns  = J(nfix, n + 1, 0)
    k    = 0
    for (i = 1; i <= n; i++) {
        if (!P[i, 5]) {
            k++
            Cns[k, i]     = 1
            Cns[k, n + 1] = P[i, 6]
        }
    }
    st_matrix(bn, bvec')
    st_matrixcolstripe(bn, stp)
    st_matrix(vn, V)
    st_matrixcolstripe(vn, stp)
    st_matrixrowstripe(vn, stp)
    if (nfix) {
        st_matrix(cn, Cns)
        st_matrixcolstripe(cn, stp \ ("", "_r"))
    }

    // ---- e(bread_ratio), weighting (a) ----------------------------------------
    gd = J(0, 1, .)
    for (g = 1; g <= G; g++) {
        gd = gd \ (diagonal(sts[g].Gamma) :/ (sts[g].N / N))
    }
    Bd = dwlsmv_breaddiag(fit.Delta, diag(gd))

    conf = subinstr(canon, ";", char(10)) + char(10)
    p = crexternal("dwlsmv_m11_gconf")
    if (p == NULL) {
        p = findexternal("dwlsmv_m11_gconf")
    }
    *p = dwlsmv_levelgen(conf, K, "configural")
    st_local("N",          strofreal(N, "%21.0g"))
    st_local("rank",       strofreal(par.npar))
    st_local("nfix",       strofreal(nfix))
    st_local("T3",         strofreal(fit.T3, "%21.0g"))
    st_local("df",         strofreal(fit.T3_df, "%21.0g"))
    dwlsmv_m11_fitidx(sts, fit.T3, fit.T3_df)
    st_local("converged",  strofreal(fit.converged))
    st_local("iterations", strofreal(fit.iterations))
    st_local("bratio",     strofreal(Bd[4, 1] == . ? . : Bd[4, 2], "%21.0g"))
}

// The posted rows P (ptype, grp, row, col, free, value, fidx) in e(b)'s order,
// e(b) as a column and e(V) = ACOV / (N - G), from a G > 1 table and its fit.
// Shared by dwlsmv_m11_post_g and difftest's refit (dwlsmv_dtfront), so the
// refit's e(b) is built by the code that posted the stored one.
void dwlsmv_m11_gbv(struct dwlsmv_parmap scalar par, struct dwlsmv_fit scalar fit,
                    real scalar N, real scalar G, real scalar nind, real scalar nfac,
                    real matrix P, real colvector bvec, real matrix V)
{
    real matrix Ev
    real scalar i, k, g, j, t, n, pos

    P = J(0, 7, .)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] <= 4) {
            P = P \ (par.ptype[i], par.grp[i], par.row[i], par.col[i],
                     par.free[i], par.value[i], par.fidx[i])
        }
    }
    for (t = 5; t <= 7; t++) {
        n = (t == 7 ? nfac : nind)
        for (g = 1; g <= G; g++) {
            for (j = 1; j <= n; j++) {
                pos = 0
                for (i = 1; i <= par.nall; i++) {
                    if (par.ptype[i] == t) {
                        if (par.grp[i] == g & par.row[i] == j) {
                            pos = i
                        }
                    }
                }
                if (pos) {
                    P = P \ (t, g, j, 1, par.free[pos], par.value[pos], par.fidx[pos])
                }
                else {
                    P = P \ (t, g, j, 1, 0, (t == 6 ? 1 : 0), 0)
                }
            }
        }
    }
    n = rows(P)
    bvec = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        bvec[i] = (P[i, 5] ? fit.theta[P[i, 7]] : P[i, 6])
    }
    Ev = dwlsmv_ev(fit.acov, N, G)
    V  = J(n, n, 0)
    for (i = 1; i <= n; i++) {
        if (P[i, 5]) {
            for (k = 1; k <= n; k++) {
                if (P[k, 5]) {
                    V[i, k] = Ev[P[i, 7], P[k, 7]]
                }
            }
        }
    }
}

end
