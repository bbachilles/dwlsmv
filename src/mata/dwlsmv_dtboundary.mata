*! dwlsmv -- DIFFTEST validation check 11: boundary restrictions forced by
*! H0's fixed values (Module 12)
*!
*! Executed under StataMP 18 on 2026-10-08: tests/unit/test_dtboundary.do passes,
*! 20 cases, and twelve mutants each turn it red (list in that file's header).
*! Drafted by a read-only agent; reviewed, run and landed by the main session.
*! Made group-aware 2026-10-13 (exemption gate step 1): test_dtboundary.do
*! unchanged and green; tests/unit/test_dtboundary_groups.do green, mutants
*! B1-B8 and P red (list in that file's header).
*!
*!     void dwlsmv_dtboundary(struct dwlsmv_parmap scalar par,
*!                            string scalar model, string scalar consequence)
*!
*! Spec Module 12, Validation contract, check 11 (decided 2026-10-08; see
*! E-091, "DIFFTEST checks 10 and 11 refined"). Runs on H0's
*! parameter table alone: no data, no theta-hat, no start value of a free
*! parameter is ever read. Returns silently or aborts with rc 3498.
*!
*! MESSAGES (decided 2026-10-13, user, mechanism (a); spec check 11). The
*! routine states the condition, the caller the consequence: consequence is
*! printed after the condition, then the routine aborts where it always has.
*! An empty consequence is refused (rc 3498), so a caller cannot fail open by
*! forgetting it. Callers: the configural preflight (dwlsmv_m11_fit_g); the
*! user-written comparison and the generated-level backstop (dwlsmv_dtstep).
*!
*! WHAT "FIXED" MEANS HERE. A row with free == 0 is fixed at its value. A row
*! ABSENT from the table is fixed at its type's absent-row default (spec 4.8
*! G4: omega 1, every other type 0; an element with no directive is not a
*! parameter). So the fixed-status matrices start at 1 everywhere and only free
*! rows clear them. A tied row is free (spec 4.8 G2), so it is never "fixed".
*!
*! G GROUPS (derivation @sec-inv-check11, C11-1 to C11-3; spec Module 12 check
*! 11, 2026-10-13). The rules run group by group on each group's own rows, the
*! key being (ptype, grp, row, col) (spec 4.8 G6; an empty or missing grp reads
*! 1, so a single-group table runs exactly as before). p and nfac come from the
*! whole table, so a row missing from one group reads its default there. Rule
*! (b) becomes Theta_gjj = 1 - omega_gj^2 (Lambda_g Phi_g Lambda_g')_jj: omega
*! multiplies at every level (Q-2, decided 2026-10-09); a free omega_gj leaves
*! the row undetermined (a boundary forced through a free parameter, not
*! refused); a fixed one enters as its value; an absent one reads 1. Messages
*! name the group when the table has more than one. The configural preflight
*! (dwlsmv_m11_fit_g) and the generated-level backstop wording (step 2 of the
*! exemption gate) are the callers' business, not this routine's.
*!
*! (a) PSI. Index sets S whose every Psi entry (diagonal and off-diagonal) is
*! fixed are the cliques of the "fixed pair" graph on the factors with a fixed
*! diagonal. Only the MAXIMAL ones are examined, and the argument is:
*!   - if a maximal block is positive definite, every principal submatrix of it
*!     is positive definite, hence nonsingular, so no subset needs a look;
*!   - if it is not positive definite it is refused.
*! The second line also refuses an indefinite block with no singular principal
*! submatrix (variances 1, a covariance fixed at 1.2): outside the parameter
*! space rather than on its boundary, refused by the user's decision of
*! 2026-10-08 (spec Module 12 check 11). What is NOT an option is checking the
*! maximal blocks for singularity alone: [[0,1],[1,0]] is nonsingular while its
*! 1x1 block (a variance fixed at 0, which the spec names) is singular. For a
*! PSD matrix, interlacing makes "maximal block singular" equivalent to "some
*! principal submatrix singular"; for an indefinite one it does not.
*!
*! The positive-definiteness test is an LDL' factorization with a strict
*! pivot > 0 test: no square root, no ratio, no cutoff. On the values users fix
*! (0, +-1, 0.5, ...) the pivots are exact, so a correlation fixed at 1 gives a
*! pivot of exactly 0. A mathematically singular block whose entries are not
*! exact in binary can leave a pivot of +-1 ulp; "exact" is literal (user, 2026-10-08).
*!
*! (b) THETA. Theta_jj = 1 - (Lambda A Psi A' Lambda')_jj with A = (I - B)^-1.
*! It is "a function of fixed values only" when, writing K_j for the factors
*! with a nonzero loading in row j (free, or fixed nonzero) and Anc for K_j and
*! all their ancestors through B's paths (free, or fixed nonzero, betas):
*!   - every loading in row j is fixed (absent = fixed 0);
*!   - every beta whose TO-factor is in Anc is fixed (a free beta into Anc
*!     enters A's rows for K_j; a free beta elsewhere cannot);
*!   - every Psi entry with both indices in Anc is fixed.
*! These are syntactic (structural) conditions. Outside Anc, (Lambda A)_jm is 0,
*! so free entries there are irrelevant and are zeroed before computing. A
*! free parameter can still drop out by numerical cancellation of fixed values;
*! Theta_jj is then "not determined" and is not refused (structural, as spec
*! check 11 reads). B is acyclic (dwls_identify() refuses cycles), so B is nilpotent and A
*! is formed as I + B + ... + B^(nfac-1), with no inversion.
*!
*! (b) uses PSI, never Phi: an endogenous factor's total variance includes its
*! explained part (spec 4.8, ptype 2). The test's endogenous case refuses only
*! if this is right.
*!
*! Load order: dwlsmv_types, dwlsmv_nfac.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Bron-Kerbosch without pivoting: every maximal clique of G (0/1, zero
// diagonal) within candidate set P, appended to out as an indicator row.
// R, P and X are indicator rowvectors. Mata passes by reference, so P and X
// are copied before the loop alters them; the recursive calls receive
// temporaries.
void dwlsmv_dtb_cliques(real rowvector R, real rowvector P, real rowvector X,
                        real matrix G, real matrix out)
{
    real rowvector Pl, Xl, Rv
    real scalar v

    if (sum(P) == 0 & sum(X) == 0) {
        if (sum(R) > 0) {
            out = out \ R
        }
        return
    }

    Pl = P
    Xl = X
    for (v = 1; v <= cols(Pl); v++) {
        if (Pl[v]) {
            Rv    = R
            Rv[v] = 1
            dwlsmv_dtb_cliques(Rv, Pl :* G[v, .], Xl :* G[v, .], G, out)
            Pl[v] = 0
            Xl[v] = 1
        }
    }
}

// 1 if the symmetric M is positive definite, 0 otherwise. LDL' with a strict
// pivot test: a pivot that is 0 (singular) or negative (indefinite) or
// missing returns 0. No square roots, so dyadic inputs give exact pivots.
real scalar dwlsmv_dtb_pd(real matrix M)
{
    real scalar n, k, i, t, d
    real matrix L
    real colvector dv

    n  = rows(M)
    L  = I(n)
    dv = J(n, 1, 0)

    for (k = 1; k <= n; k++) {
        d = M[k, k]
        for (t = 1; t < k; t++) {
            d = d - L[k, t] * L[k, t] * dv[t]
        }
        if (missing(d)) {
            return(0)
        }
        if (d <= 0) {
            return(0)
        }
        dv[k] = d
        for (i = k + 1; i <= n; i++) {
            L[i, k] = M[i, k]
            for (t = 1; t < k; t++) {
                L[i, k] = L[i, k] - L[i, t] * L[k, t] * dv[t]
            }
            L[i, k] = L[i, k] / d
        }
    }

    return(1)
}

// Group of row i: grp[i], or 1 when grp is empty or missing there (a
// single-group table built without it), as dwlsmv_dtn_grp reads it.
real scalar dwlsmv_dtb_grp(struct dwlsmv_parmap scalar par, real scalar i)
{
    if (rows(par.grp) < i) {
        return(1)
    }
    if (missing(par.grp[i])) {
        return(1)
    }
    return(par.grp[i])
}

// Item of loading row j in group g, read from a lambda row's name
// ("lambda[y1,F1]", "g2:" prefixed at G > 1); "indicator j" if none parses.
string scalar dwlsmv_dtb_item(struct dwlsmv_parmap scalar par, real scalar g,
                              real scalar j)
{
    real scalar i, a, b
    string scalar s

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1 & par.row[i] == j) {
            if (dwlsmv_dtb_grp(par, i) == g) {
                s = par.name[i]
                a = strpos(s, "lambda[")
                b = strpos(s, ",")
                if (a > 0 & b > a + 7) {
                    return(substr(s, a + 7, b - a - 7))
                }
            }
        }
    }
    return("indicator " + strofreal(j))
}

// Rules (a) and (b) on the rows of par in group g. p and nfac are the whole
// table's; glab labels messages ("" at G = 1, ", group g" otherwise);
// consequence is the caller's sentence, printed before each abort.
void dwlsmv_dtb_group(struct dwlsmv_parmap scalar par, real scalar g,
                      real scalar p, real scalar nfac,
                      string scalar model, string scalar glab,
                      string scalar consequence)
{
    real scalar i, j, k, m, r, c, q, isdet, th
    real matrix Lam, Psi, B, Lfx, Pfx, Bfx, G, cl, E, A, Bk
    real rowvector idx, nz, anc, a
    real colvector om, omfx, ing

    ing = J(par.nall, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        ing[i] = (dwlsmv_dtb_grp(par, i) == g)
    }

    // Values of FIXED rows only; free rows are zeroed and flagged. Absent
    // rows keep the defaults: 0 for lambda, psi, beta; 1 for omega (G4).
    Lam  = J(p, nfac, 0)
    Psi  = J(nfac, nfac, 0)
    B    = J(nfac, nfac, 0)
    om   = J(p, 1, 1)
    Lfx  = J(p, nfac, 1)
    Pfx  = J(nfac, nfac, 1)
    Bfx  = J(nfac, nfac, 1)
    omfx = J(p, 1, 1)
    for (i = 1; i <= par.nall; i++) {
        if (ing[i] == 0) {
            continue
        }
        r = par.row[i]
        c = par.col[i]
        if (par.ptype[i] == 1) {
            if (par.free[i]) {
                Lfx[r, c] = 0
            }
            else {
                Lam[r, c] = par.value[i]
            }
        }
        else if (par.ptype[i] == 2) {
            if (par.free[i]) {
                Pfx[r, c] = 0
                Pfx[c, r] = 0
            }
            else {
                Psi[r, c] = par.value[i]
                Psi[c, r] = par.value[i]
            }
        }
        else if (par.ptype[i] == 3) {
            if (par.free[i]) {
                Bfx[r, c] = 0
            }
            else {
                B[r, c] = par.value[i]
            }
        }
        else if (par.ptype[i] == 6) {
            if (par.free[i]) {
                omfx[r] = 0
            }
            else {
                om[r] = par.value[i]
            }
        }
    }
    if (missing(Lam) > 0 | missing(Psi) > 0 | missing(B) > 0 | missing(om) > 0) {
        errprintf("dwlsmv check 11 (%s%s): a fixed parameter value is missing.\n", model, glab)
        _error("dwlsmv check 11: a fixed parameter value is missing")
    }

    // =====================================================================
    // (a) Every all-fixed principal block of Psi_g must be positive definite.
    // =====================================================================
    G = Pfx
    for (m = 1; m <= nfac; m++) {
        G[m, m] = 0
    }
    cl = J(0, nfac, 0)
    dwlsmv_dtb_cliques(J(1, nfac, 0), diagonal(Pfx)', J(1, nfac, 0), G, cl)

    for (q = 1; q <= rows(cl); q++) {
        idx = selectindex(cl[q, .])
        if (dwlsmv_dtb_pd(Psi[idx, idx]) == 0) {
            errprintf("dwlsmv check 11 (%s%s): the fixed values in Psi force a boundary. Every entry\n",
                      model, glab)
            errprintf("of this block is fixed and the block is singular or not positive semidefinite:\n")
            for (i = 1; i <= par.nall; i++) {
                if (ing[i] & par.ptype[i] == 2) {
                    if (cl[q, par.row[i]] & cl[q, par.col[i]]) {
                        errprintf("    %s = %g\n", par.name[i], par.value[i])
                    }
                }
            }
            errprintf("(a covariance with no directive is fixed at 0; a factor correlation fixed at 1,\n")
            errprintf("two factors made one, is such a block).\n")
            errprintf("%s\n", consequence)
            exit(3498)
        }
    }

    // =====================================================================
    // (b) Every Theta_gjj determined by fixed values alone must be > 0.
    // =====================================================================
    // E[to, from] = 1 for every beta path, free or fixed nonzero; its closure
    // E[k, m] = 1 iff m is an ancestor of k.
    E = (Bfx :== 0) :| (B :!= 0)
    for (k = 1; k <= nfac; k++) {
        for (i = 1; i <= nfac; i++) {
            for (j = 1; j <= nfac; j++) {
                if (E[i, k]) {
                    if (E[k, j]) {
                        E[i, j] = 1
                    }
                }
            }
        }
    }

    // A = (I - B)^-1 over the FIXED betas, as a finite series (B nilpotent).
    A  = I(nfac)
    Bk = I(nfac)
    for (k = 1; k < nfac; k++) {
        Bk = Bk * B
        A  = A + Bk
    }

    for (j = 1; j <= p; j++) {
        if (sum(Lfx[j, .] :== 0) > 0) {
            continue
        }
        if (omfx[j] == 0) {
            continue                        // a free omega multiplies the row
        }
        nz = (Lam[j, .] :!= 0)
        if (sum(nz) == 0) {
            continue                        // Theta_jj = 1
        }
        anc = nz :| ((nz * E) :> 0)

        isdet = 1
        for (r = 1; r <= nfac; r++) {
            if (anc[r]) {
                if (sum(Bfx[r, .] :== 0) > 0) {
                    isdet = 0
                }
                for (c = 1; c <= nfac; c++) {
                    if (anc[c]) {
                        if (Pfx[r, c] == 0) {
                            isdet = 0
                        }
                    }
                }
            }
        }
        if (isdet == 0) {
            continue
        }

        a  = Lam[j, .] * A
        th = 1 - om[j] * om[j] * (a * Psi * a')
        if (missing(th)) {
            errprintf("dwlsmv check 11 (%s%s): indicator %f's residual variance is missing.\n",
                      model, glab, j)
            _error("dwlsmv check 11: a fixed residual variance evaluated to missing")
        }
        if (th <= 0) {
            errprintf("dwlsmv check 11 (%s%s): the fixed values force the residual variance of %s\n",
                      model, glab, dwlsmv_dtb_item(par, g, j))
            errprintf("to 1 - omega^2 (Lambda Phi Lambda')_jj = %g, which is zero or negative (omega = %g).\n",
                      th, om[j])
            errprintf("Its loadings:\n")
            for (i = 1; i <= par.nall; i++) {
                if (ing[i] & par.ptype[i] == 1) {
                    if (par.row[i] == j) {
                        errprintf("    %s = %g (fixed)\n", par.name[i], par.value[i])
                    }
                }
            }
            errprintf("(Phi from the fixed beta and psi values; omega reads 1 when the table has no\n")
            errprintf("omega row for the item.)\n")
            errprintf("%s\n", consequence)
            exit(3498)
        }
    }
}

void dwlsmv_dtboundary(struct dwlsmv_parmap scalar par, string scalar model,
                       string scalar consequence)
{
    real scalar nfac, p, i, g, ng

    if (strtrim(consequence) == "") {
        errprintf("dwlsmv check 11 (%s): the caller gave no consequence sentence.\n", model)
        _error("dwlsmv check 11: no consequence sentence was given")
    }

    nfac = dwlsmv_nfac(par)
    p    = 0
    ng   = 1
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 1) {
            p = max((p, par.row[i]))
        }
        ng = max((ng, dwlsmv_dtb_grp(par, i)))
    }
    // Backstop: dwls_identify() requires a metric and a loading per factor.
    if (nfac == 0 | p == 0) {
        errprintf("dwlsmv check 11 (%s): the parameter table has no factors or no loadings.\n",
                  model)
        _error("dwlsmv check 11: the parameter table has no factors or loadings")
    }

    for (g = 1; g <= ng; g++) {
        dwlsmv_dtb_group(par, g, p, nfac, model, (ng > 1 ? ", group " + strofreal(g) : ""),
                         consequence)
    }
}

end
