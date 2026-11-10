*! Module 4: the pointwise expansion test -- Level 4 of the derivation file.
*!
*! Executed under StataMP 18 on 2026-09-20 and it passes.
*!
*! THIS IS THE ONE TEST THAT SEPARATES THE CORRECT INFLUENCE FUNCTION FROM THE
*! NAIVE ONE. Everything else in the suite can pass with the two-stage
*! threshold-propagation term missing:
*!
*!   - the centring diagnostic of spec 4.5.5 cannot see it, because the naive
*!     influence function is ALSO exactly centred (measured in
*!     test_polychoric.do, both means come back 0.0);
*!   - the binary oracle of spec 4.5.3 cannot see it, because a 2x2 table is
*!     saturated and the two forms coincide there;
*!   - the threshold block of Gamma cannot see it, because it contains no
*!     correlation entries at all;
*!   - and a Monte Carlo variance comparison cannot see it at any feasible
*!     replication count -- the naive form understates N*Var(rho-hat) by about
*!     8%, while the Monte Carlo standard error of that quantity at R = 150 is
*!     around 11%. The development record gives that as the reason this test
*!     exists in this form rather than as a variance check.
*!
*! WHY IT LIVES HERE AND NOT IN MODULE 2'S TEST. An influence function makes a
*! stronger claim than a variance:
*!
*!     rho-hat - rho = (1/N) sum_i psi_rho,i + o_p(N^-1/2)
*!
*! Testing it pointwise needs an ASSEMBLED psi, and assembly is Module 4's job
*! -- Module 2 returns ingredients and never forms psi. Putting this in Module
*! 2's test would have meant building psi inside a test, which is the one thing
*! the locked API keeps in a single place. The project's notes record its absence there
*! as a decision rather than a gap; this file is the other half of it.
*!
*! HOW IT SEPARATES THEM. If the expansion holds, the residual
*! (rho-hat - rho) - (1/N) sum_i psi_i is o_p(N^-1/2), so its standard
*! deviation times sqrt(N) SHRINKS as N grows. If a term is missing, the
*! residual contains a surviving O_p(N^-1/2) piece and that product does not
*! move. A difference between "vanishes at the right rate" and "does not vanish
*! at all" is qualitative, which is why it is decisive where a percentage
*! comparison is not.
*!
*! psi IS EVALUATED AT THE POPULATION PARAMETERS, and that is essential rather
*! than incidental. Evaluated at its own estimates, any correct psi sums to
*! identically zero -- the pairwise score vanishes at rho-hat by construction
*! and the threshold influence functions vanish because q-hat is the sample
*! proportion -- so the comparison would have no content and would "pass" for
*! both forms. The weight matrix passed to dwlsmv_pairinfo() is therefore the
*! POPULATION cell probability matrix, which is also the expected-information
*! form; at the population parameters the empirical and expected normalisers
*! agree to O(N^-1/2) and neither is favoured by this test. This test is about
*! the propagation term, not about the normaliser -- see test_gamma.do check 6
*! for that.
*!
*! The naive form is built here, inline, as `u / I_rr` rather than through a
*! flag on the production routine. A production routine carrying a switch that
*! makes it deliberately wrong is a switch someone can leave flipped.
*!
*! Runtime is a few seconds: 2 sample sizes x 400 replicates, one polychoric
*! fit each. 400 rather than the derivation file's 150 because the SHRINK RATIO
*! is a ratio of two standard deviations each estimated from R replicates, and
*! at 150 it landed within 5% of its own threshold -- a tolerance sitting on
*! its boundary, which the project's notes name as the most likely and least interesting
*! cause of a milestone failure. `set seed` is bitwise portable between this
*! machine and the Stata
*! 16 milestone machine -- the rng block of tools/probe_machine.do is identical
*! on both -- so this test is reproducible there rather than merely re-run.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
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
             src/mata/dwlsmv_psi_rho_table.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_gamma_expansion: `f' failed to load"
        exit 9
    }
}

* Module 4 exists, so this test never legitimately skips: it passes or fails.

set seed 20260920

mata:

// Categorise a latent column against a threshold vector: 0 .. C.
real colvector dwlsmv_x_cut(real colvector ystar, real colvector tau)
{
    real colvector y
    real scalar c

    y = J(rows(ystar), 1, 0)
    for (c = 1; c <= rows(tau); c++) {
        y = y :+ (ystar :> tau[c])
    }
    return(y)
}

real scalar dwlsmv_x_sd(real colvector x)
{
    real scalar mu, R

    R  = rows(x)
    mu = sum(x) / R
    return(sqrt(sum((x :- mu) :^ 2) / (R - 1)))
}

// One sample size: returns (resid_correct, resid_naive, cor_correct, cor_naive).
real rowvector dwlsmv_x_run(real scalar N, real scalar R,
                            real colvector t1, real colvector t2,
                            real scalar rho,
                            real matrix psiG, real matrix psiN)
{
    struct dwlsmv_pairfit scalar f
    real colvector dev, Tc, Tn, z1, z2, y1, y2
    real matrix n
    real scalar r, sig

    dev = J(R, 1, 0)
    Tc  = J(R, 1, 0)
    Tn  = J(R, 1, 0)
    sig = sqrt(1 - rho * rho)

    for (r = 1; r <= R; r++) {
        z1 = rnormal(N, 1, 0, 1)
        z2 = rnormal(N, 1, 0, 1)
        y1 = dwlsmv_x_cut(z1, t1)
        y2 = dwlsmv_x_cut(rho :* z1 :+ sig :* z2, t2)

        // Module 2's actual estimator, on the sample's own thresholds. Spec
        // 7.3 is explicit that a Monte Carlo check only has teeth if it runs
        // the literal code path being validated, not a simpler stand-in.
        f = dwls_polychoric(y1, y2,
                            dwlsmv_thresholds(y1, rows(t1) + 1),
                            dwlsmv_thresholds(y2, rows(t2) + 1))

        n      = dwlsmv_crosstab(y1, y2, rows(t1) + 1, rows(t2) + 1)
        dev[r] = f.rho - rho
        Tc[r]  = sum(n :* psiG) / N
        Tn[r]  = sum(n :* psiN) / N
    }

    return((dwlsmv_x_sd(dev - Tc) * sqrt(N),
            dwlsmv_x_sd(dev - Tn) * sqrt(N),
            dwlsmv_x_sd(dev) * sqrt(N)))
}

real scalar dwlsmv_x_expansion()
{
    struct dwlsmv_pairfit scalar f
    string colvector FAIL
    real colvector t1, t2, q1, q2
    real matrix P, psiG, psiN, pa, pb
    real rowvector o1, o2
    real scalar rho, R, j, shrink, hold, sep

    FAIL = J(0, 1, "")

    // The derivation file's Design B: K = 4 and K = 3, rho = 0.5.
    q1  = (0.15 \ 0.50 \ 0.80)
    q2  = (0.40 \ 0.75)
    t1  = invnormal(q1)
    t2  = invnormal(q2)
    rho = 0.50
    R   = 400

    // psi at the POPULATION parameters. The weight is the population cell
    // probability matrix -- see the header for why this must not be the
    // sample's.
    P = dwlsmv_cellprob(t1, t2, rho)
    dwlsmv_pairinfo(f, t1, t2, rho, P)

    pa   = dwlsmv_psi_tau_table(q1)
    pb   = dwlsmv_psi_tau_table(q2)
    psiG = dwlsmv_psi_rho_table(f, pa, pb)

    // The naive form: the direct score with the propagation term dropped.
    // Built inline rather than through a flag on the production routine.
    psiN = f.u :/ f.I_rr

    // Both are exactly centred -- the property that makes the centring
    // diagnostic useless here. Asserted so the claim is measured, not quoted.
    if (abs(sum(P :* psiG)) > 1e-12 | abs(sum(P :* psiN)) > 1e-12) {
        FAIL = FAIL \ ("setup: one of the two psi tables is not centred (" +
                       strofreal(sum(P :* psiG)) + ", " +
                       strofreal(sum(P :* psiN)) + ")")
    }
    printf("{txt}  E[psi] correct %g, naive %g -- both centred, which is the point\n",
           sum(P :* psiG), sum(P :* psiN))
    printf("{txt}  N*Var correct %9.6f, naive %9.6f -- naive understates by %4.1f%%\n",
           sum(P :* psiG :* psiG), sum(P :* psiN :* psiN),
           100 * (1 - sum(P :* psiN :* psiN) / sum(P :* psiG :* psiG)))

    o1 = dwlsmv_x_run(2000, R, t1, t2, rho, psiG, psiN)
    o2 = dwlsmv_x_run(8000, R, t1, t2, rho, psiG, psiN)

    printf("{txt}  N = 2000 : residual*sqrt(N) correct %8.5f, naive %8.5f  (sd(dev)*sqrt(N) %8.5f)\n",
           o1[1], o1[2], o1[3])
    printf("{txt}  N = 8000 : residual*sqrt(N) correct %8.5f, naive %8.5f  (sd(dev)*sqrt(N) %8.5f)\n",
           o2[1], o2[2], o2[3])

    shrink = o2[1] / o1[1]
    hold   = o2[2] / o1[2]
    sep    = o2[2] / o2[1]
    printf("{txt}  correct shrinks by %5.3f (want <= 0.85), naive holds at %5.3f (want >= 0.80), separation %5.1fx\n",
           shrink, hold, sep)
    printf("{txt}  correct residual as a share of sd(dev): %5.3f%% at N=2000, %5.3f%% at N=8000 (want < 10%%)\n",
           100 * o1[1] / o1[3], 100 * o2[1] / o2[3])

    // THE PRIMARY ASSERTION is the share, not the ratio. "The residual is
    // negligible beside the deviation it is supposed to reproduce" is the same
    // claim as the expansion holding, and it carries a 60x margin here where
    // the shrink RATIO is itself a noisy statistic -- it is a ratio of two
    // standard deviations each estimated from R replicates, and at R = 150 it
    // moved enough to sit within 5% of its own threshold. That is precisely
    // the near-boundary tolerance the project's notes name as the most likely and least
    // interesting cause of a milestone failure, so it is not what the test
    // turns on.
    if (100 * o2[1] / o2[3] > 10) {
        FAIL = FAIL \ ("expansion: at N = 8000 the correct residual is " +
                       strofreal(100 * o2[1] / o2[3]) + "% of sd(rho-hat - rho), " +
                       "so psi does not reproduce the deviation it should")
    }
    // The rate is still checked, with the threshold moved away from its
    // boundary. Quadrupling N should roughly halve the residual; 0.85 excludes
    // "does not move", which is what the naive form does at 1.14, while
    // leaving room for Monte Carlo noise.
    if (shrink > 0.85) {
        FAIL = FAIL \ ("expansion: the correct residual did not shrink -- ratio " +
                       strofreal(shrink) + " over a fourfold N. The pointwise " +
                       "expansion does not hold, so psi_rho is wrong")
    }
    // The naive residual must NOT shrink. If it does, the two forms are not
    // being distinguished and this test has stopped doing its job.
    if (hold < 0.80) {
        FAIL = FAIL \ ("expansion: the NAIVE residual shrank too -- ratio " +
                       strofreal(hold) + ". Either the naive table is not " +
                       "actually naive or this test has lost its power")
    }
    if (sep < 3) {
        FAIL = FAIL \ ("expansion: correct and naive are only " + strofreal(sep) +
                       "x apart at N = 8000 -- not a decisive separation")
    }

    printf("{txt}test_gamma_expansion: %f failure(s)\n", rows(FAIL))
    for (j = 1; j <= rows(FAIL); j++) {
        printf("{err}    %s\n", FAIL[j])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_x_expansion()))

if `nfail' {
    display as error _n "test_gamma_expansion: `nfail' failure(s)"
    exit 9
}

display as result "test_gamma_expansion: the pointwise expansion holds for the assembled psi_rho and fails for the naive form -- the propagation term is present"
exit 0
