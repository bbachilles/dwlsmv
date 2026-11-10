*! EXECUTED 2026-10-30 on StataMP 18 (dev machine); results: E-085, section 2b.
*!
*! Usage: do tools/probe_binormal_p1r.do outdir
*!
*! delta_Phi phase 1, the P1 refinement (design note section 6, point 2). The
*! base run (tools/probe_binormal.do) found the flush boundary inside D only at
*! rho <= -0.80. This file re-runs P1 alone, finer, where that happens:
*!   rho in {-0.999} and -0.995, -0.990, ..., -0.700 (step 0.005)
*!   theta 180..270 degrees by 1 degree
*! The whole grid is CHOSEN AFTER MEASURING (E-085).
*! Same scan and bisection as P1 of the base file, and the same output columns,
*! written to outdir/pb_points.csv so tools/probe_binormal_ref.R reads it
*! unchanged. want_ref = 1 only where |h|, |k| <= 2.6 (D plus a margin, also
*! chosen after measuring), to keep the MPFR step to the points phase 2 needs.
*!
*! Run (repo root, Git Bash):
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/probe_binormal_p1r.do <outdir>
*!     grep "^PB|" <last token of outdir>.log   (Stata names the batch log so)

version 16
clear all
args outdir
if `"`outdir'"' == "" {
    display as error "usage: do tools/probe_binormal_p1r.do outdir"
    exit 198
}

mata:
mata set matastrict on

real scalar pr_lowbit(real scalar v)
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

string scalar pr_g(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

void pr_put(real scalar fh, string scalar probe, string scalar tag,
            real scalar h, real scalar k, real scalar r)
{
    real scalar v, sw, r1, r2, r3, inD, want

    v  = binormal(h, k, r)
    sw = (binormal(k, h, r) == v)
    r1 = ((normal(h) - binormal(h, -k, -r)) == v)
    r2 = ((normal(k) - binormal(-h, k, -r)) == v)
    r3 = ((normal(h) + normal(k) - 1 + binormal(-h, -k, r)) == v)
    inD = (abs(h) <= 2.5 & abs(k) <= 2.5 & abs(r) <= 0.999)
    want = (abs(h) <= 2.6 & abs(k) <= 2.6)
    fput(fh, probe + "," + tag + "," + pr_g(h) + "," + pr_g(k) + "," + pr_g(r)
        + "," + pr_g(v) + "," + strofreal(pr_lowbit(v)) + "," + strofreal(inD)
        + "," + strofreal(sw) + "," + strofreal(r1) + "," + strofreal(r2)
        + "," + strofreal(r3) + "," + strofreal(want))
}

void pr_run(string scalar outdir)
{
    real rowvector rhos, ths
    real scalar fh, i, j, q, r, th, c, s, t, ta, tb, tm, z0, z1, nchg
    real scalar nray, nmulti, ninD
    string scalar tag

    rhos = (-0.999, (-199..-140) :* 0.005)
    ths  = 180..270

    fh = fopen(outdir + "/pb_points.csv", "w")
    fput(fh, "probe,tag,h,k,rho,stata,lowbit,inD,swap_eq,refl1_eq,refl2_eq,refl3_eq,want_ref")

    nray = 0; nmulti = 0; ninD = 0
    for (q = 1; q <= cols(rhos); q++) {
        r = rhos[q]
        for (i = 1; i <= cols(ths); i++) {
            th = ths[i]
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
                    ta = (j - 1) * 0.01
                    tb = t
                    while (1) {
                        tm = (ta + tb) / 2
                        if (tm <= ta | tm >= tb) break
                        if ((binormal(tm * c, tm * s, r) == 0) == z0) ta = tm
                        else tb = tm
                    }
                    tag = "q" + strofreal(q) + "_th" + strofreal(th) + "_c" + strofreal(nchg)
                    pr_put(fh, "P1a", tag, ta * c, ta * s, r)
                    pr_put(fh, "P1b", tag, tb * c, tb * s, r)
                    if (abs(ta * c) <= 2.5 & abs(ta * s) <= 2.5 & abs(r) <= 0.999) ninD++
                    z0 = z1
                }
            }
            if (nchg > 1) nmulti++
        }
    }
    fclose(fh)
    printf("PB|P1r|rays=%g rays_with_more_than_one_change=%g boundaries_inside_D=%g\n", nray, nmulti, ninD)
}
end

mata: pr_run(st_local("outdir"))
display "PB|stata|" c(stata_version) "|" c(processors) " processors|" c(born_date)
