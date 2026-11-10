*! STUDY RUN NOT YET EXECUTED. Drafted 2026-10-31; smoke-run 2026-10-31 (docs/notes/deltaphi-phase2/smoke/).
*!
*! Usage: do tools/phase2_binormal.do outdir mode      (mode: study | smoke)
*!
*! delta_Phi phase 2 (docs/notes/deltaphi-phase2/prereg.md): point placement,
*! binormal() values and Claim 1. Placement reads binormal() only, never a
*! reference. M4/M1 and Claims 2-3: tools/phase2_binormal_ref.R.
*!
*! RNG (prereg 4): set rng mt64s, seed 20261101; each design draws from its own
*! stream, reset before it starts: 1 G-in, 2 Claim 1, 3 G-out.
*!
*! G-in (8.4d): 23 bands [-0.999,-0.99), [-0.99,-0.98), ..., [-0.78,-0.77).
*!   Rays rho ~ U[band), theta ~ U[90, 360]. Acceptance rule as the pilot
*!   (tools/probe_binormal_pilot.do): 0.001 scan for an in-D flush; 1e-4 scan,
*!   first flushed in-D t is the inner end; 1e-6 refinement over the 0.01 window
*!   below it; the designated point is the farthest flushed in-D point found.
*!   STOP RULE: a band stops at its 150th accepted ray or at its cap, whichever
*!   comes first; no further ray is drawn. Caps from the 8.4d table; the last
*!   band is terminal (reduced coverage, never inconclusive).
*! Claim 1 (layout B): 100,000 rays, same region; t_j = 0.001 j, j = 1..3536,
*!   z0 at t = 0, each change bisected to adjacent doubles; the nonzero side of
*!   each in-D change is a Claim 1 point. Claim 1 is also read on every G-out
*!   point. Claim 1: binormal() != 0 implies binormal() >= 1e-15.
*! G-out: 4,000 uniform in D; 1,000 with |rho| in [0.99, 0.999]; 1,000 with
*!   0 < binormal() < 1e-14 (cap 5e6 draws). Draws with binormal() = 0 are
*!   rejected (they are in F). Safety cap 1e6 draws on the first two strata.
*!
*! Writes outdir/p2_points.csv (G-in designated points and G-out points, for the
*! R script), outdir/p2_claim1.csv, outdir/p2_counts.csv (per band and per
*! stratum: drawn, accepted, quota, cap, terminal). Prints P2| lines.
*!
*! mode smoke: seed 1 (never the study seed), quota 3, caps 1,000, 200 Claim 1
*!   rays, G-out quotas 20/5/5 (small-output cap 1e5), and
*!   - a forced-acceptance band (every ray accepted, quota 150, cap 1,000):
*!     must write exactly 150 rows (stop rule);
*!   - a stream check: each stream, reset, must repeat its first draw;
*!   Claim 1 summaries are WITHHELD in smoke mode (structure only).
*!
*! Run (repo root, Git Bash):
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/phase2_binormal.do <outdir> study

version 16
clear all
args outdir mode
if `"`outdir'"' == "" | !inlist(`"`mode'"', "study", "smoke") {
    display as error "usage: do tools/phase2_binormal.do outdir study|smoke"
    exit 198
}
set rng mt64s

mata:
mata set matastrict on

string scalar p2_g(real scalar x)
{
    return(strtrim(strofreal(x, "%21.17g")))
}

void p2_stream(real scalar s, real scalar seed)
{
    stata("set rngstream " + strofreal(s))
    stata("set seed " + strofreal(seed))
}

real scalar p2_inD(real scalar h, real scalar k)
{
    return(abs(h) <= 2.5 & abs(k) <= 2.5)
}

real scalar p2_flushD(real scalar t, real scalar c, real scalar s, real scalar r)
{
    real scalar h, k

    h = t * c
    k = t * s
    if (!p2_inD(h, k)) return(0)
    return(binormal(h, k, r) == 0)
}

// G-in designated t, or missing if the ray is rejected (pilot rule)
real scalar p2_ray(real scalar c, real scalar s, real scalar r)
{
    real scalar j, t, lo, best, acc1

    acc1 = 0
    for (j = 1; j <= 3536; j++) {
        if (p2_flushD(j * 0.001, c, s, r)) {
            acc1 = 1
            break
        }
    }
    if (!acc1) return(.)
    lo = .
    for (j = 1; j <= 35356; j++) {
        t = j * 0.0001
        if (p2_flushD(t, c, s, r)) {
            lo = t
            break
        }
    }
    if (lo == .) return(.)
    best = lo
    for (j = 1; j <= 10000; j++) {
        t = lo - j * 1e-6
        if (t <= 0) break
        if (p2_flushD(t, c, s, r)) best = t
    }
    return(best)
}

void p2_putpt(real scalar fh, string scalar gate, string scalar stratum,
              real scalar band, real scalar id, real scalar h, real scalar k,
              real scalar r, real scalar v)
{
    fput(fh, gate + "," + stratum + "," + strofreal(band) + "," + strofreal(id)
        + "," + p2_g(h) + "," + p2_g(k) + "," + p2_g(r) + "," + p2_g(v))
}

void p2_putcount(real scalar fc, string scalar gate, string scalar stratum,
                 real scalar drawn, real scalar acc, real scalar quota,
                 real scalar cap, real scalar term)
{
    fput(fc, gate + "," + stratum + "," + strofreal(drawn) + "," + strofreal(acc)
        + "," + strofreal(quota) + "," + strofreal(cap) + "," + strofreal(term))
}

// one G-in band; force = 1 accepts every ray at t = 0.1 (smoke stop-rule check)
real scalar p2_band(real scalar fh, real scalar fc, real scalar b, real scalar lo,
                    real scalar hi, real scalar quota, real scalar cap,
                    real scalar term, real scalar force)
{
    real scalar n, nacc, r, th, c, s, td, v

    n = 0
    nacc = 0
    while (nacc < quota & n < cap) {
        n++
        r = lo + (hi - lo) * runiform(1, 1)
        th = 90 + 270 * runiform(1, 1)
        c = cos(th * pi() / 180)
        s = sin(th * pi() / 180)
        if (force) td = 0.1
        else td = p2_ray(c, s, r)
        if (td < .) {
            nacc++
            v = binormal(td * c, td * s, r)
            p2_putpt(fh, "Gin", (force ? "forced" : "band"), b, n, td * c, td * s, r, v)
        }
    }
    p2_putcount(fc, "Gin", (force ? "forced" : strofreal(b)), n, nacc, quota, cap, term)
    printf("P2|GIN|band=%g [%s,%s) drawn=%g accepted=%g quota=%g cap=%g%s%s\n", b,
        p2_g(lo), p2_g(hi), n, nacc, quota, cap, (term ? " TERMINAL" : ""),
        (force ? " FORCED" : ""))
    return(nacc)
}

// Claim 1 rays; returns (points, violations, min nonzero, count == 1e-15)
real rowvector p2_claim1(real scalar f1, real scalar nray)
{
    real scalar i, j, r, th, c, s, t, ta, tb, tm, z0, z1, h, k, v
    real scalar npt, nbad, vmin, neq

    npt = 0
    nbad = 0
    vmin = .
    neq = 0
    for (i = 1; i <= nray; i++) {
        r = -0.999 + 0.229 * runiform(1, 1)
        th = 90 + 270 * runiform(1, 1)
        c = cos(th * pi() / 180)
        s = sin(th * pi() / 180)
        z0 = (binormal(0, 0, r) == 0)
        for (j = 1; j <= 3536; j++) {
            t = j * 0.001
            z1 = (binormal(t * c, t * s, r) == 0)
            if (z1 != z0) {
                ta = (j - 1) * 0.001
                tb = t
                while (1) {
                    tm = (ta + tb) / 2
                    if (tm <= ta | tm >= tb) break
                    if ((binormal(tm * c, tm * s, r) == 0) == z0) ta = tm
                    else tb = tm
                }
                // the nonzero side of the change
                if (z0) t = tb
                else t = ta
                h = t * c
                k = t * s
                if (p2_inD(h, k)) {
                    v = binormal(h, k, r)
                    npt++
                    nbad = nbad + (v > 0 & v < 1e-15)
                    neq = neq + (v == 1e-15)
                    if (v > 0 & (vmin == . | v < vmin)) vmin = v
                    p2_putpt(f1, "C1", "ray", 0, i, h, k, r, v)
                }
                z0 = z1
            }
        }
    }
    return((npt, nbad, vmin, neq))
}

// one G-out stratum; kind 1 uniform, 2 |rho| in [0.99, 0.999], 3 small output
real rowvector p2_gout(real scalar fh, real scalar fc, real scalar kind,
                       real scalar quota, real scalar cap)
{
    real scalar n, nacc, h, k, r, v, nbad, ok
    string scalar nm

    nm = ("uniform", "rho_edge", "small")[kind]
    n = 0
    nacc = 0
    nbad = 0
    while (nacc < quota & n < cap) {
        n++
        h = 5 * runiform(1, 1) - 2.5
        k = 5 * runiform(1, 1) - 2.5
        if (kind == 2) {
            r = 0.99 + 0.009 * runiform(1, 1)
            if (runiform(1, 1) < 0.5) r = -r
        }
        else r = 1.998 * runiform(1, 1) - 0.999
        v = binormal(h, k, r)
        if (kind == 3) ok = (v > 0 & v < 1e-14)
        else ok = (v != 0)
        if (ok) {
            nacc++
            nbad = nbad + (v < 1e-15)
            p2_putpt(fh, "Gout", nm, 0, n, h, k, r, v)
        }
    }
    p2_putcount(fc, "Gout", nm, n, nacc, quota, cap, 0)
    printf("P2|GOUT|%s drawn=%g accepted=%g quota=%g cap=%g%s\n", nm, n, nacc, quota,
        cap, (nacc < quota ? " SHORT" : ""))
    return((nacc, nbad))
}

void p2_run(string scalar outdir, string scalar mode)
{
    real rowvector lo, hi, cap, c1, g
    real scalar smoke, seed, quota, fh, fc, f1, b, x1, x2, s, nc1, nbadg, ng
    real colvector q

    smoke = (mode == "smoke")
    seed = (smoke ? 1 : 20261101)
    quota = (smoke ? 3 : 150)
    lo = (-0.999, (-99..-78) :* 0.01)
    hi = ((-99..-77) :* 0.01)
    cap = (J(1, 9, 1000), J(1, 4, 2000), 3000, 3000, 4000, 4000, 5000, 8000,
           11000, 15000, 69000, 100000)
    if (smoke) cap = J(1, 23, 1000)

    fh = fopen(outdir + "/p2_points.csv", "w")
    fput(fh, "gate,stratum,band,id,h,k,rho,stata")
    fc = fopen(outdir + "/p2_counts.csv", "w")
    fput(fc, "gate,stratum,drawn,accepted,quota,cap,terminal")
    f1 = fopen(outdir + "/p2_claim1.csv", "w")
    fput(f1, "gate,stratum,band,id,h,k,rho,stata")

    // ---- G-in ----
    p2_stream(1, seed)
    for (b = 1; b <= cols(lo); b++) {
        (void) p2_band(fh, fc, b, lo[b], hi[b], quota, cap[b], (b == cols(lo)), 0)
    }

    // ---- Claim 1 ----
    p2_stream(2, seed)
    nc1 = (smoke ? 200 : 100000)
    c1 = p2_claim1(f1, nc1)

    // ---- G-out ----
    p2_stream(3, seed)
    q = (smoke ? (20, 5, 5) : (4000, 1000, 1000))'
    nbadg = 0
    ng = 0
    for (s = 1; s <= 3; s++) {
        g = p2_gout(fh, fc, s, q[s], (s == 3 ? (smoke ? 1e5 : 5e6) : 1e6))
        ng = ng + g[1]
        nbadg = nbadg + g[2]
    }

    if (smoke) {
        // stop rule on a forced-acceptance band, then stream reproducibility
        p2_stream(1, seed)
        b = p2_band(fh, fc, 99, lo[1], hi[1], 150, 1000, 0, 1)
        printf("P2|SMOKE|stop rule: forced band wrote %g rows (must be 150)\n", b)
        for (s = 1; s <= 3; s++) {
            p2_stream(s, seed)
            x1 = runiform(1, 1)
            p2_stream(s, seed)
            x2 = runiform(1, 1)
            printf("P2|SMOKE|stream %g first draw repeats: %g\n", s, (x1 == x2))
        }
        p2_stream(1, seed)
        x1 = runiform(1, 1)
        p2_stream(3, seed)
        x2 = runiform(1, 1)
        printf("P2|SMOKE|streams 1 and 3 differ: %g\n", (x1 != x2))
        printf("P2|SMOKE|Claim 1 rays=%g points=%g (summaries withheld)\n", nc1, c1[1])
    }
    else {
        printf("P2|CLAIM1|rays=%g near-cutoff points=%g violations (0<out<1e-15)=%g min nonzero=%s out==1e-15: %g\n",
            nc1, c1[1], c1[2], p2_g(c1[3]), c1[4])
        printf("P2|CLAIM1|G-out points=%g violations=%g\n", ng, nbadg)
        printf("P2|CLAIM1|%s\n", (c1[2] + nbadg == 0 ? "holds on the test points" : "CONTRADICTED"))
    }
    fclose(fh)
    fclose(fc)
    fclose(f1)
}
end

mata: p2_run(st_local("outdir"), st_local("mode"))
display "P2|stata|" c(stata_version) "|" c(processors) " processors|" c(born_date)
