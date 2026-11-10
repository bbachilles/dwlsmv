*! dwlsmv -- one pair's information quantities, at an explicit weight matrix
*! (Module 2, Gate 3; spec 4.5.2 and the derivation file's @sec-normaliser)
*!
*! Executed under StataMP 18 on 2026-09-20. tests/unit/test_polychoric.do checks
*! it three ways: at W = pi it reproduces the compact expected form to 2.8e-17, at
*! W = pihat it equals the literal per-observation crossproduct to 2.2e-16, and the
*! two differ by 3.5e-02 on an over-identified 4x4 so the argument demonstrably
*! matters.
*!
*! Fills f.u, f.I_rr, f.I_ra and f.I_rb from (a, b, rho) and a weight matrix:
*!
*!     I_rho_rho(W) = sum_rs  W_rs u_rs^2
*!     I_rho_a_c(W) = sum_rs  W_rs u_rs (d_{a_c} pi_rs / pi_rs)
*!                  = sum_s   w_cs [ W_cs u_cs / pi_cs
*!                                 - W_{c+1,s} u_{c+1,s} / pi_{c+1,s} ]
*!
*! with u = d_rho pi / pi. The second line is the compact form: a_c appears in
*! exactly the two rows it separates, so the sum runs over one variable's
*! categories rather than the whole table.
*!
*! WHY THE WEIGHT IS AN ARGUMENT. The derivation's two-step expansion needs
*! A = -E[d g_i / d rho] and G_a = E[d g_i / d a_c], and ANY consistent
*! estimator of those expectations gives a valid influence function. Two are
*! natural and the choice is not a correctness question:
*!
*!     W = pi     the EXPECTED information -- the model-implied cell
*!                probabilities, which is how the Gate 3 result is written
*!     W = pihat  the EMPIRICAL (outer-product / BHHH) information -- the
*!                observed cell proportions n_rs/N, which makes both quantities
*!                literal per-observation crossproducts
*!
*! They differ at O(N^-1/2) and converge to the same Gamma; dropping the
*! propagation term, by contrast, converges to a different one. The derivation
*! file measures exactly that: the empirical-vs-expected gap halves as N
*! quadruples while the naive-vs-expected gap sits at 8.2% and does not move.
*!
*! v1 SHIPS W = pihat, because lavaan does. This is the same compatibility
*! argument spec 4.5.4 already makes one level up for Gamma itself: an
*! O(N^-1/2) disagreement with the primary reference is reported as a
*! standard-error bug, not read as a modelling convention. The size is not
*! academic -- on case_D at n = 5000 the expected form differs from
*! lavInspect(fit, "gamma") by 0.051 absolute and 6.56 relative, against
*! 1.2e-09 relative for the empirical form. lavaan's own construction was read
*! rather than inferred: lav_m84_a22 is `sum(sc_cor[,i]^2)` and lav_m84_a21 is
*! `crossprod(sc_cor[,pstar], pair$dx_th_i)`.
*!
*! THE EXPECTED FORM IS NOT DISCARDED, which is why this takes an argument
*! instead of hard-coding pihat. W = pi is the unit-test oracle: it must
*! reproduce the w-weighted compact form that Module 2's test already validates
*! against numerical differentiation, so the shipped path and the verified path
*! are the same code. It is also the candidate fallback if the empirical form
*! proves the worse behaved under sparse cells -- pihat gives an empty cell zero
*! weight and drops it, where pi retains it. That is open (E-001).
*!
*! WHERE THE TWO COINCIDE EXACTLY, which is worth knowing before debugging a
*! difference that cannot exist: whenever the table is SATURATED, pihat = pi by
*! construction. That covers every univariate threshold -- so the whole
*! threshold block of Gamma is untouched by this choice -- and every binary
*! pair. Cases A and B of spec Section 5 therefore cannot detect the difference
*! at all. Only K > 2 separates them.
*!
*! Requires dwlsmv_types, dwlsmv_cellprob, dwlsmv_dcellprob_rho, dwlsmv_wmat.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_pairinfo(struct dwlsmv_pairfit scalar f,
                     real colvector a, real colvector b,
                     real scalar rho, real matrix Wt)
{
    real matrix P, Dr, u, Wa, Wb, uT, PT, WT, bad
    real scalar K1, K2, c, d

    K1 = rows(a) + 1
    K2 = rows(b) + 1

    if (rows(Wt) != K1 | cols(Wt) != K2) {
        errprintf("dwlsmv_pairinfo: the weight matrix is %s x %s, expected %s x %s.\n",
                  strofreal(rows(Wt)), strofreal(cols(Wt)),
                  strofreal(K1), strofreal(K2))
        _error("dwlsmv_pairinfo: weight matrix is not conformable with the table")
    }
    // A weight matrix is a distribution over the cells -- either the model
    // probabilities or the observed proportions. Handing this routine raw
    // COUNTS is the easy mistake and it is silent: every information comes
    // back scaled by N, psi by 1/N, and Gamma by 1/N^2, all of which look like
    // numbers. So it is refused rather than trusted.
    if (abs(sum(Wt) - 1) > 1e-8) {
        errprintf("dwlsmv_pairinfo: the weight matrix sums to %s, not 1.\n",
                  strofreal(sum(Wt)))
        errprintf("  Pass cell PROPORTIONS (n:/N) or model cell probabilities," +
                  " never raw counts --\n")
        errprintf("  counts would scale every information by N and nothing" +
                  " would report it.\n")
        _error("dwlsmv_pairinfo: weight matrix does not sum to 1")
    }

    // Positivity is required only where the weight is. A cell with W = 0
    // enters no sum below and is never looked up by Module 4, and near
    // rho = 1 an EMPTY corner cell's probability is below binormal()'s
    // resolution and comes back 0 or -1 ulp (derivation @sec-emptycells). On
    // exactly those cells -- W = 0 and P <= 0 -- u is zeroed, as the
    // derivation's psi_tab0 does. The sums below would not show the difference:
    // sum() treats missing as 0, so a 0 * (x/0) term would silently drop out
    // rather than poison them, and a cell with W = 0 adds nothing anyway. What
    // the zeroing protects is f.u itself: dwlsmv_psi_rho_table() builds psi
    // over every cell from it, and without the zeroing u there is missing
    // (P = 0) or Dr/(-1 ulp) (tests L9, L9b). A cell that CARRIES weight and is
    // not positive still aborts: that is a degenerate configuration, not an
    // underflow nobody uses.
    P  = dwlsmv_cellprob_raw(a, b, rho)
    Dr = dwlsmv_dcellprob_rho(a, b, rho)
    if (min(P :+ (Wt :== 0)) <= 0) {
        errprintf("dwlsmv_pairinfo: a cell with weight %s has model probability %s at rho = %s.\n",
                  strofreal(max(Wt :* (P :<= 0)), "%10.3e"),
                  strofreal(min(P :+ (Wt :== 0)), "%10.3e"), strofreal(rho, "%12.9f"))
        errprintf("  Only a cell carrying no weight may underflow; see dwlsmv_cellprob().\n")
        _error("dwlsmv_pairinfo: a weighted cell has a non-positive probability")
    }
    bad = (Wt :== 0) :& (P :<= 0)
    P   = P :* (1 :- bad) :+ bad
    Dr  = Dr :* (1 :- bad)
    u   = Dr :/ P

    f.u    = u
    f.I_rr = sum(Wt :* u :* u)

    Wa     = dwlsmv_wmat(a, b, rho)
    f.I_ra = J(K1 - 1, 1, 0)
    for (c = 1; c <= K1 - 1; c++) {
        f.I_ra[c] = sum(Wa[c, .] :* (Wt[c,     .] :* u[c,     .] :/ P[c,     .] :-
                                     Wt[c + 1, .] :* u[c + 1, .] :/ P[c + 1, .]))
    }

    // The same expression with the roles swapped -- arguments reversed and the
    // table transposed -- rather than a mirrored twin of the loop above. The
    // derivation writes the identity d_b pi = [d_a pi(b,a)]', and a second
    // routine would be a second place for the same index error to live.
    Wb     = dwlsmv_wmat(b, a, rho)
    uT     = u'
    PT     = P'
    WT     = Wt'
    f.I_rb = J(K2 - 1, 1, 0)
    for (d = 1; d <= K2 - 1; d++) {
        f.I_rb[d] = sum(Wb[d, .] :* (WT[d,     .] :* uT[d,     .] :/ PT[d,     .] :-
                                     WT[d + 1, .] :* uT[d + 1, .] :/ PT[d + 1, .]))
    }
}

end
