*! EXECUTED 2026-10-30 on StataMP 18 (dev machine); results: E-085.
*!
*! Usage: do tools/probe_binormal.do outdir
*!
*! delta_Phi certification, phase 1: characterize Stata's binormal() over v1's
*! delta_Phi domain D (|rho| <= 0.999, h, k in [-2.5, 2.5]). Design:
*! E-085. NOT a test: no pass/fail. It prints
*! PB| lines and writes pb_points.csv (every probed point, its binormal()
*! value and the P4 symmetry flags) for tools/probe_binormal_ref.R, which adds
*! the MPFR references. Every regime this finds is OBSERVED, NOT DOCUMENTED.
*!
*! Probes (design note section 2):
*!   P1  flush boundary along rays (h, k) = t (cos th, sin th), th 180..270 by
*!       5 degrees, at rho in {-0.999, -0.95, ..., 0.95, 0.999}: a coarse scan
*!       of t in [0, 40] by 0.01, then bisection of every 0/nonzero change down
*!       to adjacent t. More than one change on a ray is reported, not assumed
*!       away. Fixed before any run (user, 2026-10-30, option b); refinement
*!       points added later are CHOSEN AFTER MEASURING and are not in this file.
*!   P2  small-output structure: 50 points stepping in from each first
*!       boundary (t_pos - 0.002 j); per value, the exponent e of its lowest
*!       set bit (v = m 2^-e, m odd) and its significant bits.
*!   P3  rho scan at nine fixed (h, k): uniform rho = -0.999 + 0.001 i, plus
*!       bands |rho| in {0.3, 0.75, 0.925} +- 0.002 by 1e-5, both signs (user,
*!       2026-10-30, option b; the bands test a HYPOTHESIS about the method).
*!       Pairs 7-9 of the nine are Arm A's in-domain over-budget corners:
*!       CHOSEN AFTER MEASURING (E-063).
*!   P4  on every P1-P3 point: swap and the three reflection identities, bit
*!       for bit.
*!   P5  anchors without the integrand: rho = 0 against normal(h)*normal(k)
*!       (A-L certifies normal() on (beta_h, 0] only; anything else is a
*!       diagnostic), and the margin binormal(h, 8e307, rho) against normal(h).
*!
*! Run (repo root, Git Bash):
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/probe_binormal.do <outdir>
*!     grep "^PB|" <last token of outdir>.log   (Stata names the batch log so)

version 16
clear all
args outdir
if `"`outdir'"' == "" {
    display as error "usage: do tools/probe_binormal.do outdir"
    exit 198
}

mata:
mata set matastrict on

real scalar pb_lowbit(real scalar v)
{
    real scalar x, e

    if (v == 0) return(.)
    x = abs(v)
    e = 0
    while (x != floor(x) & e < 1100) {
        x = x * 2
        e++
    }
    return(e)
}

string scalar pb_g(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

// one output row: probe, tag, h, k, rho, value, then the P4 flags
void pb_put(real scalar fh, string scalar probe, string scalar tag,
            real scalar h, real scalar k, real scalar r, real scalar want_ref)
{
    real scalar v, sw, r1, r2, r3, inD

    v  = binormal(h, k, r)
    sw = (binormal(k, h, r) == v)
    r1 = ((normal(h) - binormal(h, -k, -r)) == v)
    r2 = ((normal(k) - binormal(-h, k, -r)) == v)
    r3 = ((normal(h) + normal(k) - 1 + binormal(-h, -k, r)) == v)
    inD = (abs(h) <= 2.5 & abs(k) <= 2.5 & abs(r) <= 0.999)
    fput(fh, probe + "," + tag + "," + pb_g(h) + "," + pb_g(k) + "," + pb_g(r)
        + "," + pb_g(v) + "," + strofreal(pb_lowbit(v)) + "," + strofreal(inD)
        + "," + strofreal(sw) + "," + strofreal(r1) + "," + strofreal(r2)
        + "," + strofreal(r3) + "," + strofreal(want_ref))
}

void pb_run(string scalar outdir)
{
    real rowvector rhos, ths
    real matrix hk
    real scalar fh, i, j, q, r, th, c, s, t, ta, tb, tm, z0, z1, nchg
    real scalar nray, nmulti, ninD, v, e, emax, sbmin, h, k, nP2
    real scalar d, dmax0, dmaxo, nflush0, pflush0
    real scalar nmarg, neqm
    string scalar tag

    rhos = (-0.999, (-19..19) :* 0.05, 0.999)
    ths  = (36..54) :* 5

    fh = fopen(outdir + "/pb_points.csv", "w")
    fput(fh, "probe,tag,h,k,rho,stata,lowbit,inD,swap_eq,refl1_eq,refl2_eq,refl3_eq,want_ref")

    // ---- P1 and P2 ----------------------------------------------------------
    nray = 0; nmulti = 0; ninD = 0; emax = 0; sbmin = 53; nP2 = 0
    for (q = 1; q <= cols(rhos); q++) {
        r = rhos[q]
        for (i = 1; i <= cols(ths); i++) {
            th = ths[i]
            // exact axes at 180 and 270 degrees (cos/sin of pi are not exact)
            c = cos(th * pi() / 180)
            s = sin(th * pi() / 180)
            if (th == 180) {
                c = -1
                s = 0
            }
            if (th == 270) {
                c = 0
                s = -1
            }
            nray++
            nchg = 0
            z0 = (binormal(0, 0, r) == 0)
            for (j = 1; j <= 4000; j++) {
                t = j * 0.01
                z1 = (binormal(t * c, t * s, r) == 0)
                if (z1 != z0) {
                    nchg++
                    ta = (j - 1) * 0.01; tb = t
                    while (1) {
                        tm = (ta + tb) / 2
                        if (tm <= ta | tm >= tb) break
                        if ((binormal(tm * c, tm * s, r) == 0) == z0) ta = tm
                        else tb = tm
                    }
                    tag = "r" + strofreal(q) + "_th" + strofreal(th) + "_c" + strofreal(nchg)
                    pb_put(fh, "P1a", tag, ta * c, ta * s, r, 1)
                    pb_put(fh, "P1b", tag, tb * c, tb * s, r, 1)
                    if (abs(ta * c) <= 2.5 & abs(ta * s) <= 2.5 & abs(r) <= 0.999) ninD++
                    if (nchg == 1 & z0 == 0) {
                        // P2: step in from the last nonzero point of the first change
                        for (d = 0; d < 50; d++) {
                            tm = ta - 0.002 * d
                            if (tm < 0) break
                            h = tm * c; k = tm * s
                            v = binormal(h, k, r)
                            if (v == 0) continue
                            e = pb_lowbit(v)
                            emax = max((emax, e))
                            sbmin = min((sbmin, e + floor(ln(v) / ln(2)) + 1))
                            nP2++
                            pb_put(fh, "P2", tag + "_d" + strofreal(d), h, k, r, mod(d, 10) == 0)
                        }
                    }
                    z0 = z1
                }
            }
            if (nchg > 1) {
                nmulti++
                printf("PB|P1|multi|rho=%s th=%g changes=%g\n", pb_g(r), th, nchg)
            }
            if (nchg == 0) printf("PB|P1|nochange|rho=%s th=%g (no flush up to t = 40)\n", pb_g(r), th)
        }
    }
    printf("PB|P1|rays=%g rays_with_more_than_one_change=%g boundaries_inside_D=%g\n", nray, nmulti, ninD)
    printf("PB|P2|values=%g max_lowbit_exponent=%g min_significant_bits=%g\n", nP2, emax, sbmin)

    // ---- P3 ---------------------------------------------------------------------
    hk = (-2.5, -2.5 \ -2.5, 0 \ -1, -1 \ 0, 0 \ -2.5, 2.5 \ -2, 1 \
          -2.5, -2.25 \ -1.75, 1.25 \ -2, 0.75)
    for (i = 1; i <= rows(hk); i++) {
        tag = "hk" + strofreal(i)
        for (j = 0; j <= 1998; j++) pb_put(fh, "P3u", tag, hk[i, 1], hk[i, 2], -0.999 + 0.001 * j, 1)
        for (q = 1; q <= 3; q++) {
            c = (0.3, 0.75, 0.925)[q]
            for (j = -200; j <= 200; j++) {
                pb_put(fh, "P3b", tag, hk[i, 1], hk[i, 2],  c + 1e-5 * j, 1)
                pb_put(fh, "P3b", tag, hk[i, 1], hk[i, 2], -c + 1e-5 * j, 1)
            }
        }
    }
    printf("PB|P3|pairs=%g (pairs 7-9 chosen after measuring)\n", rows(hk))

    // ---- P5 -----------------------------------------------------------------------
    dmax0 = 0; dmaxo = 0; nflush0 = 0; pflush0 = 0
    for (i = -80; i <= 50; i++) {
        h = (i < -50 ? -2.5 + 0.5 * (i + 50) : i * 0.05)
        for (j = i; j <= 50; j++) {
            k = (j < -50 ? -2.5 + 0.5 * (j + 50) : j * 0.05)
            v = binormal(h, k, 0)
            d = abs(v - normal(h) * normal(k)) / 2^-52
            if (h <= 0 & k <= 0 & h > -8.3743889230674569 & k > -8.3743889230674569) dmax0 = max((dmax0, d))
            else dmaxo = max((dmaxo, d))
            if (v == 0 & normal(h) * normal(k) > 0) {
                nflush0++
                pflush0 = max((pflush0, normal(h) * normal(k)))
            }
            pb_put(fh, "P5", "rho0", h, k, 0, 0)
        }
    }
    printf("PB|P5|rho0 max|bn - normal*normal| / 2^-52: A-L domain=%s other (diagnostic)=%s\n", pb_g(dmax0), pb_g(dmaxo))
    printf("PB|P5|rho0 flushed to 0 with normal*normal > 0: %g, largest such product=%s\n", nflush0, pb_g(pflush0))
    nmarg = 0; neqm = 0
    for (q = 1; q <= cols(rhos); q++) {
        for (i = -50; i <= 50; i++) {
            nmarg++
            neqm = neqm + (binormal(i * 0.05, 8e307, rhos[q]) == normal(i * 0.05))
        }
    }
    printf("PB|P5|margin binormal(h, 8e307, rho) == normal(h) bitwise: %g of %g\n", neqm, nmarg)
    fclose(fh)

    printf("PB|done|P4 flags are in pb_points.csv; summary by tools/probe_binormal_ref.R\n")
}
end

mata: pb_run(st_local("outdir"))
display "PB|stata|" c(stata_version) "|" c(processors) " processors|" c(born_date)
