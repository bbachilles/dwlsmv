*! dwlsmv -- Stage 1 as a whole: s, Gamma and the map (Modules 1-4, spec 5.1)
*!
*! Executed under StataMP 18 on 2026-09-20 and it passes. On case_D its thresholds
*! agree with lavaan to 3.3e-15, its Gamma is symmetric to 0.0 with minimum
*! eigenvalue 0.442, and its centring diagnostic is 2.3e-16 on both blocks.
*!
*! TAKES (Y, K) AND NOTHING MODEL-IMPLIED, and that is the constraint the
*! locked API exists to protect. No model specification, no starting values, no
*! theta. Stage 1 produces s and Gamma with no knowledge of what structural
*! model will eventually be fitted to them, which is exactly what makes spec
*! 7.1's paired validation meaningful: this whole routine is checkable against
*! lavaan with no model specified at all. A signature change that let anything
*! model-implied reach here would quietly remove that.
*!
*! The order of operations is the two-stage estimator of Muthen (1984), and the
*! order is the point rather than an implementation detail:
*!
*!     1. thresholds, per variable, from the cumulative sample proportions
*!     2. the staircase lookup tables those thresholds imply (spec 4.5.1)
*!     3. each pairwise polychoric, with the thresholds HELD FIXED (Module 2)
*!     4. Gamma, assembled from 2 and 3 together (Module 4)
*!
*! Step 3 holding a random quantity fixed is what creates the propagation term
*! that step 4 carries. Steps 1-3 are not independent M-estimators and Gamma is
*! wrong if they are treated as such.
*!
*! WHAT IS NOT DONE HERE.
*!
*! Listwise deletion. dwlsmv_listwise() is deliberately not called: deletion
*! changes N, and N is what every influence function is divided by, so a
*! deletion buried inside Stage 1 would make the N that Gamma was built with
*! invisible to the caller reporting it. The caller deletes, sees rows() move,
*! and hands this routine a matrix whose row count is the N it will report.
*! Module 11 does the authoritative deletion at the dataset level, because
*! e(sample) cannot be marked from a matrix whose rows are already gone.
*!
*! Which pairs enter s and Gamma (spec Module 2, decided 2026-10-03). Status 1
*! and status 6 reached by the endpoint check enter. Every other pair aborts
*! Stage 1 (rc 3498) after ALL pairs are classified and reported: status 5,
*! zero-probe status 6 (status 6 at three evaluations), a perfect-association
*! pair (a 2x2 table with two empty cells on a diagonal; @sec-kindE: no
*! interior maximum), the other status-2 pairs and statuses 3, 4, 7 and 8.
*! The optional keep_bound = 1 instead includes the status-2 pairs other than
*! empty-cell 2x2 pairs at +-B, as lavaan does, and marks them in st.bound.
*! Degenerate CATEGORIES abort upstream in dwlsmv_thresholds().
*!
*! The one-empty-cell repair (O6, decided 2026-10-06; spec Module 2,
*! *One-empty-cell repair*). A binary pair whose 2x2 table has exactly one
*! empty cell (dwlsmv_emptycell_pairs(), read from the table, never from a
*! status) is refitted on the repaired table n+ (dwlsmv_emptycell_repair(),
*! c = 0.5, margins kept) with the same root finder, B and caps. If that root
*! is accepted by the rule above, rho-hat+ enters s with bound = 0, and its
*! information quantities are taken at rho-hat+ with the RAW weights n/N:
*! construction (a), lavaan's and Mplus's, so Gamma stays (1/N) sum psi psi'
*! over the real observations and psi-bar_rho != 0 by construction. If it is
*! not accepted (|rho-hat+| >= B reaches status 2), the pair aborts (O3):
*! nothing is truncated. Repaired pairs are printed on every call. No other
*! table is repaired, collapsed or corrected.
*!
*! The centring diagnostic is carried in the returned struct, not asserted on.
*! See dwlsmv_psibar() for what it can and cannot detect -- in particular it
*! cannot detect a missing propagation term, because the naive influence
*! function is also exactly centred.
*!
*! Requires the whole of Modules 1, 2 and 4.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The status in words (spec Module 2's table). A zero-probe status 6 is
// worded apart: it is the one that never reached a bracket.
string scalar dwlsmv_ss_statusword(real scalar st, real scalar it)
{
    if (st == 2) return("no sign change up to the bound")
    if (st == 3) return("no sign change up to the usable end")
    if (st == 4) return("gap in the usable region")
    if (st == 5) return("more than one sign change")
    if (st == 6 & it == 3) return("ambiguous sign at an exact zero of the score at 0")
    if (st == 6) return("ambiguous endpoint sign")
    if (st == 7) return("evaluation cap")
    if (st == 8) return("unresolvable point inside the bracket")
    return("root")
}

// The readings of a correlation at or near the bound, stated once. Module 10's
// rule: the fact first, a reading second, never a diagnosis.
void dwlsmv_ss_readings()
{
    errprintf("  A correlation at or near the bound has several possible readings:\n")
    errprintf("  near-duplicate (redundant) items, strong skew that collapses the\n")
    errprintf("  table, or too few observations. The data cannot tell these apart,\n")
    errprintf("  and this note is not a diagnosis.\n")
}

struct dwlsmv_stage1 scalar dwls_sample_stats(real matrix Y, real colvector K,
                                              | real scalar keep_bound)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_pairfit scalar f
    struct dwlsmv_pairfit colvector pairs
    real matrix psi_tau, n
    real colvector tau_i, tau_j, cls, ecell, rep
    real rowvector r
    real scalar j, c, k, pos, vi, vj, kb, nbad, nkeep, nbad2, allkeepable, nrep
    string scalar line

    // keep_bound: 0 (default) aborts on every pair outside the v1 policy; 1
    // includes the status-2 pairs other than empty-cell 2x2 pairs at +-B (spec 5.1).
    if (args() < 3) {
        kb = 0
    }
    else {
        kb = keep_bound
    }
    if (!(kb == 0 | kb == 1)) {
        errprintf("dwls_sample_stats: keep_bound is %s; it must be 0 or 1.\n",
                  strofreal(kb))
        _error("dwls_sample_stats: keep_bound must be 0 or 1")
    }

    st.map = dwlsmv_statmap_build(K)

    if (cols(Y) != st.map.p) {
        errprintf("dwls_sample_stats: Y has %s columns but K names %s variables.\n",
                  strofreal(cols(Y)), strofreal(st.map.p))
        _error("dwls_sample_stats: Y and K disagree on the number of variables")
    }
    st.N = rows(Y)
    if (st.N < 2) {
        _error("dwls_sample_stats: Y needs at least two observations")
    }

    st.s = J(st.map.D, 1, .)

    // ---- 1. thresholds ------------------------------------------------------
    // invnormal() of a cumulative sample proportion, and nothing more. No
    // optimizer, no starting value, no iteration; a Stage 1 that needs one of
    // those is a signal that something upstream is wrong.
    for (j = 1; j <= st.map.p; j++) {
        tau_i = dwlsmv_thresholds(Y[., j], st.map.K[j])
        for (c = 1; c <= st.map.C[j]; c++) {
            st.s[dwlsmv_pos_tau(st.map, j, c)] = tau_i[c]
        }
    }

    // ---- 2. the staircase tables -------------------------------------------
    psi_tau = dwlsmv_psi_tau_stack(Y, st.map)

    // ---- 3. the pairwise polychorics ---------------------------------------
    pairs = dwlsmv_pairfit(st.map.P, 1)
    // cls: 0 enters, 1 included at the bound under keep_bound, 2 aborts.
    cls   = J(st.map.P, 1, 0)
    ecell = J(st.map.P, 1, 0)
    // rep: the one-empty-cell pairs, from the table (dwlsmv_emptycell_pairs.mata).
    rep   = dwlsmv_emptycell_pairs(Y, st.map)

    for (k = 1; k <= st.map.P; k++) {
        pos = st.map.T + k
        vi  = st.map.var_i[pos]
        vj  = st.map.var_j[pos]

        // Thresholds are sliced out of s rather than recomputed. The map makes
        // a variable's thresholds contiguous, so this is a slice and not a
        // gather -- and re-running the estimator would be a second route to a
        // quantity already in hand.
        tau_i = st.s[|st.map.Coff[vi] + 1 \ st.map.Coff[vi] + st.map.C[vi]|]
        tau_j = st.s[|st.map.Coff[vj] + 1 \ st.map.Coff[vj] + st.map.C[vj]|]

        f = dwls_polychoric(Y[., vi], Y[., vj], tau_i, tau_j)

        // dwls_polychoric() is not told which variables it was handed, which
        // is deliberate -- a pair estimator has no business knowing its
        // position in the model. The assembler fills them, and this is the
        // assembler.
        f.i = vi
        f.j = vj

        // The one-empty-cell repair: refit on n+, information at the RAW
        // weights (construction (a)). The raw fit above is kept only as the
        // classification it would have had; it is replaced here.
        if (rep[k]) {
            n = dwlsmv_crosstab(Y[., vi], Y[., vj], 2, 2)
            r = dwlsmv_rootfind(&dwlsmv_rf_eval(), dwlsmv_emptycell_repair(n),
                                tau_i, tau_j, (., .), 0.999, 40, 54)
            f.rho        = r[1]
            f.status     = r[2]
            f.converged  = (f.status == 1)
            f.iterations = r[3]
            dwlsmv_pairinfo(f, tau_i, tau_j, f.rho, n :/ st.N)
        }

        pairs[k]  = f
        st.s[pos] = f.rho

        // Status 1 enters. Status 6 enters only from the endpoint check, which
        // evaluates at least once per side, so it has five or more
        // evaluations; the zero-probe branch returns at exactly three.
        if (f.status == 1 | (f.status == 6 & f.iterations != 3)) {
            continue
        }
        cls[k] = 2
        if (f.status == 2 & !rep[k]) {
            // An empty-cell 2x2 pair is read from the observed table, never inferred from
            // "2x2 at status 2" (decided 2026-10-02). A 2x2 table reaching
            // Module 2 has both margins positive, so any empty cell is one
            // cell or a diagonal pair; neither has an interior maximum. The
            // one-cell case was repaired above, so here it is the diagonal
            // pair: the perfect-association pair, never kept.
            if (st.map.K[vi] == 2 & st.map.K[vj] == 2) {
                n = dwlsmv_crosstab(Y[., vi], Y[., vj], 2, 2)
                ecell[k] = any(n :== 0)
            }
            if (!ecell[k] & kb) {
                cls[k] = 1
            }
        }
    }

    st.bound = (cls :== 1)
    nbad  = sum(cls :== 2)
    nkeep = sum(cls :== 1)

    if (nbad) {
        nbad2 = 0
        allkeepable = 1
        errprintf("dwls_sample_stats: %s of %s pairs cannot enter the sample statistics:\n",
                  strofreal(nbad), strofreal(st.map.P))
        for (k = 1; k <= st.map.P; k++) {
            if (cls[k] != 2) continue
            f    = pairs[k]
            line = "  " + st.map.label[st.map.T + k] + ": rho-hat " +
                   strofreal(f.rho, "%9.4f") + ", status " + strofreal(f.status) +
                   " (" + dwlsmv_ss_statusword(f.status, f.iterations) + ")"
            if (ecell[k]) {
                line = line + "; 2x2 table with two empty cells: no interior maximum, not repaired"
            }
            if (rep[k]) {
                line = line + "; after the one-empty-cell repair"
                if (f.status == 2) {
                    line = line + ", the repaired correlation reaches the bound"
                }
            }
            errprintf("%s\n", line)
            if (f.status == 2) nbad2 = nbad2 + 1
            if (f.status != 2 | ecell[k] | rep[k]) allkeepable = 0
        }
        if (nbad2) {
            dwlsmv_ss_readings()
        }
        if (allkeepable) {
            errprintf("  These pairs can be included at the bound, as lavaan does, with\n")
            errprintf("  keep_bound = 1; their effect on standard errors and T3 is not\n")
            errprintf("  established.\n")
        }
        errprintf("  The statuses are defined in spec Module 2.\n")
        exit(3498)
    }

    if (nkeep) {
        errprintf("dwls_sample_stats: %s of %s pairs are included at the bound (keep_bound):\n",
                  strofreal(nkeep), strofreal(st.map.P))
        for (k = 1; k <= st.map.P; k++) {
            if (cls[k] != 1) continue
            errprintf("  %s: rho-hat %s, status 2 (no sign change up to the bound)\n",
                      st.map.label[st.map.T + k], strofreal(pairs[k].rho, "%9.4f"))
        }
        errprintf("  They enter s and Gamma at the bound, as lavaan does. Their effect on\n")
        errprintf("  standard errors and T3 is not established: a pre-registered experiment\n")
        errprintf("  measured it as small in one design and did not support inclusion.\n")
        dwlsmv_ss_readings()
    }

    // An intentional repair is reported, not raised as an error: printf, not
    // errprintf. Positional labels; Module 11 prints the names.
    nrep = sum(rep)
    if (nrep) {
        printf("dwls_sample_stats: one-empty-cell repair applied to %s of %s pairs:\n",
               strofreal(nrep), strofreal(st.map.P))
        for (k = 1; k <= st.map.P; k++) {
            if (!rep[k]) continue
            n = dwlsmv_crosstab(Y[., st.map.var_i[st.map.T + k]],
                                Y[., st.map.var_j[st.map.T + k]], 2, 2)
            printf("  %s: rho-hat %s, table (%s, %s \\ %s, %s)\n",
                   st.map.label[st.map.T + k], strtrim(strofreal(pairs[k].rho, "%9.4f")),
                   strofreal(n[1, 1]), strofreal(n[1, 2]), strofreal(n[2, 1]),
                   strofreal(n[2, 2]))
        }
        printf("  0.5 is added to the empty cell and its diagonal partner and taken from\n")
        printf("  the other two; margins and thresholds are unchanged. lavaan and Mplus\n")
        printf("  apply the same repair and the same Gamma. Its effect on standard errors\n")
        printf("  and T3 at small N is not established.\n")
    }

    // ---- 4. Gamma, and the diagnostic --------------------------------------
    st.Gamma  = dwls_gamma(Y, st.map, pairs, psi_tau)
    st.psibar = dwlsmv_psibar(Y, st.map, pairs, psi_tau)

    return(st)
}

end
