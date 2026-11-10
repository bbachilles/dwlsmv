*! EXECUTED 2026-10-31 on StataMP 18 (dev machine); log: E-086.
*!
*! Usage: do tools/probe_binormal_pilot.do outdir
*!
*! delta_Phi phase 2 G-in ACCEPTANCE PILOT (user, 2026-10-31). A design input,
*! never evidence: placement only, binormal() alone, no reference read. It
*! applies phase 2's exact G-in acceptance rule to random rays drawn as phase 2
*! will draw them, so that the per-band counts are binomial:
*!   per band b (23 bands: [-0.999,-0.99), [-0.99,-0.98), ..., [-0.78,-0.77)),
*!   2,000 rays with rho ~ U[band) and theta ~ U[90, 360] degrees;
*!   rng mt64, seed 7310 (pilot only; distinct from rate.do and the study).
*! The rule (prereg section 8.4a, 8.4d):
*!   1. scan t_j = 0.001 j, j = 1..3536; reject the ray if no in-D point has
*!      binormal() = 0;
*!   2. scan t_j = 0.0001 j, j = 1..35356; flushed in-D t values, zones split
*!      at gaps > 0.05; the inner end is the smallest t of the first zone;
*!   3. refine the 0.01 window below the inner end at step 1e-6; the
*!      designated point is the farthest flushed in-D point found.
*!   The ray is ACCEPTED when step 1 passes and step 2 finds a flushed in-D
*!   point (step 3 can then only move the point, staying in D).
*! Writes outdir/pilot_rays.csv (band, rho, theta, acc1, acc, t_desig) and
*! prints PB|PILOT| lines.

version 16
clear all
args outdir
if `"`outdir'"' == "" {
    display as error "usage: do tools/probe_binormal_pilot.do outdir"
    exit 198
}
set rng mt64
set seed 7310

mata:
mata set matastrict on

string scalar pl_g(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

real scalar pl_flushD(real scalar t, real scalar c, real scalar s, real scalar r)
{
    real scalar h, k

    h = t * c
    k = t * s
    if (abs(h) > 2.5 | abs(k) > 2.5) return(0)
    return(binormal(h, k, r) == 0)
}

// designated t, or missing if the ray is rejected; acc1 returns step 1
real scalar pl_ray(real scalar c, real scalar s, real scalar r, real scalar acc1)
{
    real scalar j, t, lo, best

    acc1 = 0
    for (j = 1; j <= 3536; j++) {
        if (pl_flushD(j * 0.001, c, s, r)) {
            acc1 = 1
            break
        }
    }
    if (!acc1) return(.)
    lo = .
    for (j = 1; j <= 35356; j++) {
        t = j * 0.0001
        if (pl_flushD(t, c, s, r)) {
            lo = t
            break
        }
    }
    if (lo == .) return(.)
    best = lo
    for (j = 1; j <= 10000; j++) {
        t = lo - j * 1e-6
        if (t <= 0) break
        if (pl_flushD(t, c, s, r)) best = t
    }
    return(best)
}

void pl_run(string scalar outdir)
{
    real rowvector lo, hi
    real scalar fh, b, i, r, th, c, s, td, a1, nacc, n1

    lo = (-0.999, (-99..-78) :* 0.01)
    hi = ((-99..-77) :* 0.01)
    fh = fopen(outdir + "/pilot_rays.csv", "w")
    fput(fh, "band,rho,theta,acc1,acc,t_desig")
    for (b = 1; b <= cols(lo); b++) {
        nacc = 0
        n1 = 0
        for (i = 1; i <= 2000; i++) {
            r = lo[b] + (hi[b] - lo[b]) * runiform(1, 1)
            th = 90 + 270 * runiform(1, 1)
            c = cos(th * pi() / 180)
            s = sin(th * pi() / 180)
            a1 = 0
            td = pl_ray(c, s, r, a1)
            n1 = n1 + a1
            nacc = nacc + (td < .)
            fput(fh, strofreal(b) + "," + pl_g(r) + "," + pl_g(th) + "," + strofreal(a1)
                + "," + strofreal(td < .) + "," + pl_g(td))
        }
        printf("PB|PILOT|band=[%s,%s) drawn=2000 step1=%g accepted=%g\n", pl_g(lo[b]), pl_g(hi[b]), n1, nacc)
    }
    fclose(fh)
}
end

mata: pl_run(st_local("outdir"))
display "PB|stata|" c(stata_version) "|" c(processors) " processors|" c(born_date)
