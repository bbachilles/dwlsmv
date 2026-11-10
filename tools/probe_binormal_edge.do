*! EXECUTED 2026-10-31 on StataMP 18 (dev machine); M1 summary: E-086.
*!
*! Usage: do tools/probe_binormal_edge.do outdir
*!
*! delta_Phi phase 2 DIAGNOSTIC, not a criterion (user, 2026-10-31): does the
*! edge of each flicker zone move at finer resolution? Same rays as
*! tools/probe_binormal_res.do (eight rho, theta 90..360 by 1 degree). Per ray:
*!   1. scan t_j = 0.0001 j (j = 1..35356) and keep every in-D t where
*!      binormal() = 0; group them into zones, a new zone after a gap > 0.05 in
*!      t (CHOSEN AFTER MEASURING: zones are at most 0.013 wide);
*!   2. at each zone end (lo = smallest t, hi = largest t), scan a window of
*!      0.01 in t outward from the zone at steps 1e-5 and 1e-6, and record the
*!      flushed point farthest from the zone found at each step.
*! Writes outdir/pe_edges.csv: rho, theta, zone, end, t and (h, k) of the
*! extreme flushed point at 1e-4, 1e-5 and 1e-6. Placement only: no
*! reference, no error read. M1 values are added by
*! E-086.
*!
*! Run (repo root, Git Bash):
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/probe_binormal_edge.do <outdir>

version 16
clear all
args outdir
if `"`outdir'"' == "" {
    display as error "usage: do tools/probe_binormal_edge.do outdir"
    exit 198
}

mata:
mata set matastrict on

string scalar pe_g(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

real scalar pe_flushD(real scalar t, real scalar c, real scalar s, real scalar r)
{
    real scalar h, k

    h = t * c
    k = t * s
    if (abs(h) > 2.5 | abs(k) > 2.5) return(0)
    return(binormal(h, k, r) == 0)
}

// extreme flushed t in the window beyond t0, direction dir (-1 or +1), step dt
real scalar pe_extreme(real scalar t0, real scalar dir, real scalar dt,
                       real scalar c, real scalar s, real scalar r)
{
    real scalar j, n, t, best

    best = t0
    n = round(0.01 / dt)
    for (j = 1; j <= n; j++) {
        t = t0 + dir * j * dt
        if (t <= 0) break
        if (pe_flushD(t, c, s, r)) best = t
    }
    return(best)
}

void pe_put(real scalar fh, real scalar r, real scalar th, real scalar z,
            string scalar e, real scalar t4, real scalar t5, real scalar t6,
            real scalar c, real scalar s)
{
    fput(fh, pe_g(r) + "," + strofreal(th) + "," + strofreal(z) + "," + e
        + "," + pe_g(t4) + "," + pe_g(t4 * c) + "," + pe_g(t4 * s)
        + "," + pe_g(t5) + "," + pe_g(t5 * c) + "," + pe_g(t5 * s)
        + "," + pe_g(t6) + "," + pe_g(t6 * c) + "," + pe_g(t6 * s))
}

void pe_run(string scalar outdir)
{
    real rowvector rhos, ths
    real colvector tz
    real scalar fh, q, i, j, r, th, c, s, t, nz, a, b, z, nedge, n5, n6
    real scalar lo5, lo6, hi5, hi6

    rhos = (-995, -950, -900, -850, -820, -800, -790, -785) :* 0.001
    ths  = 90..360
    fh = fopen(outdir + "/pe_edges.csv", "w")
    fput(fh, "rho,theta,zone,end,t4,h4,k4,t5,h5,k5,t6,h6,k6")
    nedge = 0; n5 = 0; n6 = 0
    for (q = 1; q <= cols(rhos); q++) {
        r = rhos[q]
        for (i = 1; i <= cols(ths); i++) {
            th = ths[i]
            c = cos(th * pi() / 180)
            s = sin(th * pi() / 180)
            if (th == 90) {
                c = 0
                s = 1
            }
            if (th == 180) {
                c = -1
                s = 0
            }
            if (th == 270) {
                c = 0
                s = -1
            }
            if (th == 360) {
                c = 1
                s = 0
            }
            tz = J(35356, 1, .)
            nz = 0
            for (j = 1; j <= 35356; j++) {
                t = j * 0.0001
                if (pe_flushD(t, c, s, r)) {
                    nz++
                    tz[nz] = t
                }
            }
            if (nz == 0) continue
            z = 0
            a = 1
            for (j = 1; j <= nz; j++) {
                if (j == nz) b = 1
                else b = (tz[j + 1] - tz[j] > 0.05)
                if (b) {
                    z++
                    lo5 = pe_extreme(tz[a], -1, 1e-5, c, s, r)
                    lo6 = pe_extreme(tz[a], -1, 1e-6, c, s, r)
                    hi5 = pe_extreme(tz[j], 1, 1e-5, c, s, r)
                    hi6 = pe_extreme(tz[j], 1, 1e-6, c, s, r)
                    pe_put(fh, r, th, z, "lo", tz[a], lo5, lo6, c, s)
                    pe_put(fh, r, th, z, "hi", tz[j], hi5, hi6, c, s)
                    nedge = nedge + 2
                    n5 = n5 + (lo5 != tz[a]) + (hi5 != tz[j])
                    n6 = n6 + (lo6 != tz[a]) + (hi6 != tz[j])
                    a = j + 1
                }
            }
        }
    }
    fclose(fh)
    printf("PB|EDGE|edges=%g moved_at_1e-5=%g moved_at_1e-6=%g\n", nedge, n5, n6)
}
end

mata: pe_run(st_local("outdir"))
display "PB|stata|" c(stata_version) "|" c(processors) " processors|" c(born_date)
