*! EXECUTED 2026-10-31 on StataMP 18 (dev machine); results: E-085, section 2c.
*!
*! Usage: do tools/probe_binormal_p1m.do outdir
*!
*! delta_Phi phase 1, SUPPLEMENT over the mixed-sign quadrants (correction of
*! 2026-10-31, user option a): P1 and its refinement had theta in [180, 270]
*! only, but F meets D elsewhere (binormal(-2.5, 0.5, -0.999) = 0). This file
*! scans the mixed quadrants on the phase 2 grid convention:
*!   rho in {-0.999, -0.99, -0.98, ..., 0.99, 0.999} (step 0.01)
*!   theta 90..179 and 271..360 degrees by 1 degree (exact axes at 90, 360)
*!   t_j = 0.001 j, j = 1..3536 (t <= 3.536 >= 2.5 sqrt 2: all of D on the ray)
*! The grid is CHOSEN AFTER MEASURING. Same bisection and output columns as
*! P1, written to outdir/pb_points.csv for tools/probe_binormal_ref.R.
*! want_ref = 1 where |h|, |k| <= 2.6.
*!
*! Run (repo root, Git Bash):
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/probe_binormal_p1m.do <outdir>
*!     grep "^PB|" <last token of outdir>.log   (Stata names the batch log so)

version 16
clear all
args outdir
if `"`outdir'"' == "" {
    display as error "usage: do tools/probe_binormal_p1m.do outdir"
    exit 198
}

mata:
mata set matastrict on

real scalar pm_lowbit(real scalar v)
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

string scalar pm_g(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

void pm_put(real scalar fh, string scalar probe, string scalar tag,
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
    fput(fh, probe + "," + tag + "," + pm_g(h) + "," + pm_g(k) + "," + pm_g(r)
        + "," + pm_g(v) + "," + strofreal(pm_lowbit(v)) + "," + strofreal(inD)
        + "," + strofreal(sw) + "," + strofreal(r1) + "," + strofreal(r2)
        + "," + strofreal(r3) + "," + strofreal(want))
}

void pm_run(string scalar outdir)
{
    real rowvector rhos, ths
    real scalar fh, i, j, q, r, th, c, s, t, ta, tb, tm, z0, z1, nchg
    real scalar nray, nmulti, ninD
    string scalar tag

    rhos = (-0.999, (-99..99) :* 0.01, 0.999)
    ths  = (90..179, 271..360)

    fh = fopen(outdir + "/pb_points.csv", "w")
    fput(fh, "probe,tag,h,k,rho,stata,lowbit,inD,swap_eq,refl1_eq,refl2_eq,refl3_eq,want_ref")

    nray = 0; nmulti = 0; ninD = 0
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
            if (th == 360) {
                c = 1
                s = 0
            }
            nray++
            nchg = 0
            z0 = (binormal(0, 0, r) == 0)
            for (j = 1; j <= 3536; j++) {
                t = j * 0.001
                z1 = (binormal(t * c, t * s, r) == 0)
                if (z1 != z0) {
                    nchg++
                    ta = (j - 1) * 0.001
                    tb = t
                    while (1) {
                        tm = (ta + tb) / 2
                        if (tm <= ta | tm >= tb) break
                        if ((binormal(tm * c, tm * s, r) == 0) == z0) ta = tm
                        else tb = tm
                    }
                    tag = "q" + strofreal(q) + "_th" + strofreal(th) + "_c" + strofreal(nchg)
                    pm_put(fh, "P1a", tag, ta * c, ta * s, r)
                    pm_put(fh, "P1b", tag, tb * c, tb * s, r)
                    if (abs(ta * c) <= 2.5 & abs(ta * s) <= 2.5 & abs(r) <= 0.999) ninD++
                    z0 = z1
                }
            }
            if (nchg > 1) nmulti++
        }
    }
    fclose(fh)
    printf("PB|P1m|rays=%g rays_with_more_than_one_change=%g boundaries_inside_D=%g\n", nray, nmulti, ninD)
}
end

mata: pm_run(st_local("outdir"))
display "PB|stata|" c(stata_version) "|" c(processors) " processors|" c(born_date)
