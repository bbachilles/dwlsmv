*! Machine-comparison probe. NOT a test -- it has no pass/fail and the suite
*! runner ignores it (the runner globs tests/unit for test_*.do only).
*!
*! Purpose: decide whether two machines agree BITWISE on the operations this
*! estimator depends on. Run it on each machine, diff the PROBE| lines:
*!
*!     cd <repo root>
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/probe_machine.do
*!     grep "^PROBE|" probe_machine.log > probe_<machine-label>.txt
*!
*! Then: diff probe_devbox.txt probe_laptop.txt
*! Identical output means the two machines are bitwise equivalent on every
*! operation exercised here. That is the result we want and it is decisive.
*!
*! Processor count is pinned inside, because MP varies reduction order with it
*! and that difference is already known and measured -- see E-040.
*! Pinning it removes core count as a variable so that anything left is
*! architecture or version.
*!
*! Input is fixtures/data/case_D.csv, committed and therefore byte-identical
*! across machines. Do not substitute generated data.
*!
*! Everything prints in %21x -- exact hexadecimal float. Do not compare the
*! decimal displays; they round and will hide exactly the differences this is
*! looking for.
*!
*! WHAT COUNTS AS A DIFFERENCE
*!
*! Compare the p1 block across machines. That is the clean comparison: one
*! processor, so reduction order is fixed and core count cannot contribute.
*!
*! The env| lines are EXPECTED to differ -- different Stata version, core count
*! and machine type is the whole point. They are recorded, not compared.
*!
*! The p1 and p4 blocks are EXPECTED to differ from each other WITHIN a single
*! machine, on cross11, crosssum, invsym*, chol11 and eig*, and by about 1 ulp
*! on qcross*. That is the measured core-count effect of E-040, not
*! a finding. Comparing p1 against p4 tells you nothing new.
*!
*! A real result is the p1 block differing BETWEEN machines. The rng| block
*! differing is more serious still: it means the two machines do not even agree
*! on generated data, and nothing downstream is comparable.
*!
*! RESULT ON RECORD, re-measured 2026-09-20
*!
*! Development machine (StataMP 18, 6-core, native x86-64, licence serial
*! [dev]) against the milestone machine (StataMP 16.1, 4-core, Windows
*! ARM running the x86 build under emulation, serial [milestone]):
*! 143 of 146 computational lines identical bitwise.
*!
*! The three that differ are all symeigenvalues(), all by exactly 1 ulp:
*! p1 eig1, p4 eig1, p4 eigsum. Everything else matches exactly, including all
*! 35 binormal() values, normal(), lnnormal(), invnormal(), cross(),
*! quadcross(), invsym(), cholesky(), exp(), ln() and sqrt().
*!
*! The 2026-09-18 version of this note said "138 of 141". The routines it named
*! were right and the count was not; 146 is the number of computational PROBE|
*! lines this file emits. Recorded rather than quietly corrected, because a
*! count nobody can reproduce is how a reader decides the rest is approximate.
*!
*! cholesky() deserves its own line. Stata 17 moved Mata's decompositions onto
*! LAPACK/Intel MKL and names cholesky() explicitly in that whatsnew -- and it
*! agrees BITWISE across 16 and 18 here. So the LAPACK migration is not a
*! general source of cross-version difference in this project, and the 1-ulp
*! gap is specific to the eigen path rather than evidence of a broad shift.
*! symeigenvalues() is used only for the rank checks of spec 4.8 and 4.5.6,
*! which discriminate on a sixteen-order-of-magnitude singular-value gap.
*!
*! The milestone machine also shows ZERO drift between its 2026-09-18 and
*! 2026-09-20 runs: all 146 computational lines byte-identical, two days and a
*! full Stata suite apart. Only the env| block changed, because this file was
*! fixed to record c(MP).
*!
*! The rng| block is identical too, so `set seed` is bitwise portable between
*! the machines and the Section 7.3 Monte Carlo work can be run on either.
*!
*! EXTENDED 2026-09-23: svd() and luinv() were appended to the end of the
*! battery, 12 lines per block, so the file now emits 170 computational PROBE|
*! lines. The counts above describe the 146 that existed before.
*!
*! RESULT, 2026-09-23 Modules 3/5/6 milestone: the 146 are bitwise unchanged on
*! the milestone machine. In the new p1 lines against StataMP 18, sv1, svsum,
*! luinv11 and luinvsum match bitwise; svmin moves 3 ulp and one null-vector
*! entry 1 ulp; svnoise and the two structurally-zero null-vector entries differ
*! at rounding-residue level, as they must; svnull_lt_1e12 is 1 on both. So
*! svd() sits with symeigenvalues() on the eigen path, and luinv() with
*! cholesky() and invsym(). The committed baseline now includes these lines.
*! Detail in E-090.
*!
*! EXTENDED 2026-09-23 for Module 8: cholinv() appended after luinv(), 6 lines
*! per block (revised 2026-09-24), so the file now emits 182 computational
*! PROBE| lines.
*!
*! RESULT, 2026-09-24 Module 8 milestone: the 170 are bitwise unchanged on the
*! milestone machine (sixth zero-drift run). The six new p1 lines match
*! StataMP 18 bitwise, so cholinv() sits with cholesky(), invsym() and luinv();
*! zerocol 1, sing14 1, sing2 0, nearsing 0 on both. At p4 cholinvsum differs
*! by 1 ulp between the machines. The committed baseline now includes these
*! lines. Detail in E-090.
*!
*! The milestone machine's full output is committed as
*! tools/probe_baseline_stata16.txt. Diff a future run against it rather than
*! against a fresh run of this machine -- that file is the reference, and it
*! records the environment it came from (Stata 16.1, 4 processors).
*!
*! Trap worth knowing: c(machine_type) reads "PC (64-bit x86-64)" on the ARM
*! machine, because the emulated x86 build cannot see through the emulation.
*! That field does NOT identify the architecture. stata_version together with
*! processors_max is what distinguishes the two machines.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

display "PROBE|env|stata_version|" c(stata_version)
display "PROBE|env|version|" c(version)
* c(flavor) does NOT identify MP. It reads "IC" on StataMP -- verified on
* StataMP 18, where c(flavor) is IC while c(MP) is 1. Reading the flavor line
* alone led to a confident and wrong conclusion that the milestone machine was
* Stata/IC. c(MP) is the field that answers the question, so both are recorded
* and the misleading one is labelled. This is the second field in this probe
* that does not mean what its name suggests; see c(machine_type) below.
display "PROBE|env|flavor_MISLEADING|" c(flavor)
display "PROBE|env|MP|" c(MP)
display "PROBE|env|SE|" c(SE)
display "PROBE|env|processors_max|" c(processors_max)
display "PROBE|env|os|" c(os)
display "PROBE|env|machine_type|" c(machine_type)
display "PROBE|env|bit|" c(bit)

* ---------------------------------------------------------------------------
* RNG determinism. Stata's generator is integer-based and should be portable;
* if this block differs, nothing downstream is comparable and the rest of the
* output is uninterpretable.
* ---------------------------------------------------------------------------
set seed 20260918
quietly set obs 5
quietly generate double rr = rnormal()
mata:
rv = st_data(., "rr")
for (i = 1; i <= 5; i++) printf("PROBE|rng|%f|%21x\n", i, rv[i])
end

* ---------------------------------------------------------------------------
* The battery, defined once so the echoed source does not pollute the grep.
* ---------------------------------------------------------------------------
mata:
void battery(real matrix Y, string scalar tag)
{
    real matrix X, S, Q, Si, C, A, U, V
    real colvector L, sv
    real rowvector z, p, rh
    real scalar i, j

    X = Y :- mean(Y)
    S = cross(X, X) / rows(X)
    Q = quadcross(X, X) / rows(X)

    printf("PROBE|%s|cross11|%21x\n",    tag, S[1,1])
    printf("PROBE|%s|crosssum|%21x\n",   tag, sum(S))
    printf("PROBE|%s|qcross11|%21x\n",   tag, Q[1,1])
    printf("PROBE|%s|qcrosssum|%21x\n",  tag, sum(Q))

    Si = invsym(Q)
    printf("PROBE|%s|invsym11|%21x\n",   tag, Si[1,1])
    printf("PROBE|%s|invsymsum|%21x\n",  tag, sum(Si))

    C = cholesky(Q)
    printf("PROBE|%s|chol11|%21x\n",     tag, C[1,1])
    printf("PROBE|%s|cholsum|%21x\n",    tag, sum(C))

    L = symeigenvalues(Q)'
    printf("PROBE|%s|eig1|%21x\n",       tag, L[1])
    printf("PROBE|%s|eigsum|%21x\n",     tag, sum(L))

    /* Special functions. These matter most: the polychoric step is built on
       normal() and binormal(), and a library difference would surface here
       long before it surfaced as a wrong correlation. */
    z = (-2.5, -1, -0.25, 0, 0.33, 1.5, 2.75)
    p = (0.001, 0.025, 0.1, 0.5, 0.75, 0.975, 0.999)
    rh = (-0.8, -0.3, 0, 0.45, 0.9)

    for (i = 1; i <= cols(z); i++) printf("PROBE|%s|normal|%f|%21x\n",    tag, i, normal(z[i]))
    for (i = 1; i <= cols(z); i++) printf("PROBE|%s|lnnormal|%f|%21x\n",  tag, i, lnnormal(z[i]))
    for (i = 1; i <= cols(p); i++) printf("PROBE|%s|invnormal|%f|%21x\n", tag, i, invnormal(p[i]))

    for (i = 1; i <= cols(rh); i++) {
        for (j = 1; j <= cols(z); j++) {
            printf("PROBE|%s|binormal|%f|%f|%21x\n", tag, i, j, binormal(z[j], z[j]*0.5 - 0.2, rh[i]))
        }
    }

    printf("PROBE|%s|exp|%21x\n",   tag, exp(1.234567))
    printf("PROBE|%s|ln|%21x\n",    tag, ln(3.14159))
    printf("PROBE|%s|sqrt2|%21x\n", tag, sqrt(2))

    /* svd() and luinv(), added 2026-09-23 for the Modules 3/5/6 milestone.
       Appended LAST so every line above keeps its position and still diffs
       against the pre-2026-09-23 baseline. svd() is what the section 4.8 rank
       check reads, so A has the rank check's shape: tall, with one column
       that is a sum of two others, so one singular value is rounding noise.

       svnoise is that noise value and is NOT expected to match across
       versions -- it is whatever the decomposition's rounding leaves. The
       project relies only on it sitting far below 1e-8 relative, which is
       what svnull_lt_1e12 records. sv1, svsum and svmin (the smallest GENUINE
       singular value) are the real comparison. Null-vector entries print in
       absolute value because a singular vector's sign is arbitrary. */
    A = X[|1, 1 \ 20, cols(X)|]
    A = A, A[., 1] + A[., 2]
    U  = J(0, 0, .)
    V  = J(0, 0, .)
    sv = J(0, 1, .)
    svd(A, U, sv, V)
    printf("PROBE|%s|sv1|%21x\n",            tag, sv[1])
    printf("PROBE|%s|svsum|%21x\n",          tag, sum(sv))
    printf("PROBE|%s|svmin|%21x\n",          tag, sv[rows(sv) - 1])
    printf("PROBE|%s|svnull_lt_1e12|%f\n",   tag, sv[rows(sv)] / sv[1] < 1e-12)
    printf("PROBE|%s|svnoise|%21x\n",        tag, sv[rows(sv)])
    for (i = 1; i <= cols(V); i++) printf("PROBE|%s|svnullvec|%f|%21x\n", tag, i, abs(V[rows(V), i]))

    Si = luinv(Q)
    printf("PROBE|%s|luinv11|%21x\n",    tag, Si[1,1])
    printf("PROBE|%s|luinvsum|%21x\n",   tag, sum(Si))

    /* cholinv(), added 2026-09-23 for Module 8, appended last for the same
       reason as svd(). dwlsmv_acov() inverts the sandwich bread with it.
       Besides the value, four behaviours are recorded. cholinv() returns
       missing only when a FLOATING-POINT pivot is <= 0, so whether a matrix
       that is singular in exact arithmetic aborts depends on rounding:

       cholinv_zerocol_missing -- a bread with an identically zero column,
         the one case whose pivot is 0 under any rounding. This is the only
         abort dwlsmv_acov() can rely on. Must read 1.
       cholinv_sing14_missing  -- [14 14; 14 14]: its second pivot rounds to
         exactly 0 on StataMP 18 (reads 1). Recorded because the reviewer
         argued a fused multiply-add path could leave it +4.2e-16; a 0 on
         another machine would be that, not a version regression.
       cholinv_sing2_missing   -- [2 2; 2 2]: pivot +4.4e-16, a finite
         inverse near 2.3e15. Reads 0 on StataMP 18.
       cholinv_nearsing_missing -- A'A with A's summed column, rounding
         noise: a finite inverse, entries ~1.3e14. Reads 0 on StataMP 18.

       History, 2026-09-23/24: the first version had one line, expecting a
       singular matrix to return missing; its 0 corrected dwlsmv_acov()'s
       header. The next version used [14 14; 14 14] as "the exactly singular
       case" -- which worked only because 14 rounds to a zero pivot. */
    Si = cholinv(Q)
    printf("PROBE|%s|cholinv11|%21x\n",  tag, Si[1,1])
    printf("PROBE|%s|cholinvsum|%21x\n", tag, sum(Si))
    printf("PROBE|%s|cholinv_zerocol_missing|%f\n",  tag, missing(cholinv((14, 0 \ 0, 0))) > 0)
    printf("PROBE|%s|cholinv_sing14_missing|%f\n",   tag, missing(cholinv((14, 14 \ 14, 14))) > 0)
    printf("PROBE|%s|cholinv_sing2_missing|%f\n",    tag, missing(cholinv((2, 2 \ 2, 2))) > 0)
    printf("PROBE|%s|cholinv_nearsing_missing|%f\n", tag, missing(cholinv(quadcross(A, A))) > 0)
}
end

quietly import delimited using "fixtures/data/case_D.csv", clear varnames(1)
mata: YD = st_data(., .)
mata: printf("PROBE|data|rows|%f\n", rows(YD))
mata: printf("PROBE|data|cols|%f\n", cols(YD))
mata: printf("PROBE|data|sum|%21x\n", sum(YD))

* p1 is the cross-machine comparison: one processor, so reduction order is
* fixed and core count cannot contribute. p4 is available on both machines and
* is the bonus check.
foreach np in 1 4 {
    capture set processors `np'
    if _rc {
        display "PROBE|skip|processors|`np'"
        continue
    }
    mata: battery(YD, "p`=c(processors)'")
}
