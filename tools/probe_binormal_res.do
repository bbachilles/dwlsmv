*! EXECUTED 2026-10-31 on StataMP 18 (dev machine); comparison: E-086.
*!
*! Usage: do tools/probe_binormal_res.do outdir
*!
*! delta_Phi phase 2 DIAGNOSTIC, not a criterion (user, 2026-10-31): does a
*! 0.0001 step in t find crossings the 0.001 acceptance scan missed? Same
*! rays as tools/probe_binormal_acc.do at eight rho values (chosen after
*! measuring: band interiors and the edge), every 1-degree theta in [90, 360]:
*!   rho in {-0.995, -0.95, -0.90, -0.85, -0.82, -0.80, -0.79, -0.785}
*!   t_j = 0.0001 j, j = 1..35356, each change bisected to adjacent doubles
*! Placement only: no reference, no error read.
*!
*! Run (repo root, Git Bash):
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/probe_binormal_res.do <outdir>
*!     grep "^PB|" <last token of outdir>.log   (Stata names the batch log so)

version 16
clear all
args outdir
if `"`outdir'"' == "" {
    display as error "usage: do tools/probe_binormal_res.do outdir"
    exit 198
}

mata:
mata set matastrict on

real scalar pq_lowbit(real scalar v)
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

string scalar pq_g(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

void pq_put(real scalar fh, string scalar probe, string scalar tag,
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
    fput(fh, probe + "," + tag + "," + pq_g(h) + "," + pq_g(k) + "," + pq_g(r)
        + "," + pq_g(v) + "," + strofreal(pq_lowbit(v)) + "," + strofreal(inD)
        + "," + strofreal(sw) + "," + strofreal(r1) + "," + strofreal(r2)
        + "," + strofreal(r3) + "," + strofreal(want))
}

void pq_run(string scalar outdir)
{
    real rowvector rhos, ths
    real scalar fh, i, j, q, r, th, c, s, t, ta, tb, tm, z0, z1, nchg
    real scalar nray, nmulti, ninD
    string scalar tag

    rhos = (-995, -950, -900, -850, -820, -800, -790, -785) :* 0.001  // as in the acc scan
    ths  = 90..360

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
            nray++
            nchg = 0
            z0 = (binormal(0, 0, r) == 0)
            for (j = 1; j <= 35356; j++) {
                t = j * 0.0001
                z1 = (binormal(t * c, t * s, r) == 0)
                if (z1 != z0) {
                    nchg++
                    ta = (j - 1) * 0.0001
                    tb = t
                    while (1) {
                        tm = (ta + tb) / 2
                        if (tm <= ta | tm >= tb) break
                        if ((binormal(tm * c, tm * s, r) == 0) == z0) ta = tm
                        else tb = tm
                    }
                    tag = "q" + strofreal(q) + "_th" + strofreal(th) + "_c" + strofreal(nchg)
                    pq_put(fh, "P1a", tag, ta * c, ta * s, r)
                    pq_put(fh, "P1b", tag, tb * c, tb * s, r)
                    if (abs(ta * c) <= 2.5 & abs(ta * s) <= 2.5 & abs(r) <= 0.999) ninD++
                    z0 = z1
                }
            }
            if (nchg > 1) nmulti++
        }
    }
    fclose(fh)
    printf("PB|RES|rays=%g rays_with_more_than_one_change=%g boundaries_inside_D=%g\n", nray, nmulti, ninD)
}
end

mata: pq_run(st_local("outdir"))
display "PB|stata|" c(stata_version) "|" c(processors) " processors|" c(born_date)
