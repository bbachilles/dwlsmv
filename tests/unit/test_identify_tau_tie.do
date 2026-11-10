*! dwls_identify: a tied threshold starts at the sample threshold of its
*! lowest-numbered member group (spec 4.8 O5a, amended 2026-11-01, user).
*!
*! Expected starts are literals, recomputed outside Mata in R (see
*! E-074, 2026-11-01). s(g, j, k, K) = qnorm(k/K) + .02 j + .05 (g-1).
*!
*! Cases:
*!   T        K 2,3,3; tie g2,g3, k = 2: accepted (aborted in dwlsmv_pos_tau before)
*!   silent   K 3,4,4; tie g2,g3, k = 2: starts at g2's 2/4 split, not g1's 2/3
*!   control  K 3,3,3; tie g2,g3, k = 2: starts at g2's, not g1's
*!   reversed K 3,3,3,3; tie written g4 then g3: starts at g3's
*!   incl     K 3,3,3; tie g1,g2,g3: starts at g1's, as before
*!   mismatch K 3,3,4; tie g2,g3 is refused by G7's K rule, rc 3498
*!   fit      y1 K 3,4,4 on simulated data, all y1 thresholds tied in g2,g3, nu and
*!            omega free in g3 (without them tau separates from the rest of the
*!            objective and its start cannot matter): the fit from the rule's start
*!            and from the former start reach one solution. The tolerance 1e-5
*!            was chosen after measuring.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    quietly do "`f'"
}

mata:
string scalar ttt_text(real scalar G, string scalar tie, | string scalar x3)
{
    string scalar s
    real scalar g

    s = "groups " + strofreal(G) + char(10)
    for (g = 1; g <= G; g++) {
        s = s + "group " + strofreal(g) + char(10) + "factor F" + char(10) +
            "ind y1 y2 y3 y4" + char(10) + "lambda F y1 free" + char(10) +
            "lambda F y2 free" + char(10) + "lambda F y3 free" + char(10) +
            "lambda F y4 free" + char(10)
        if (g == 1) {
            s = s + "metric unitvariance F" + char(10)
        }
        else {
            s = s + "psi F F fixed 1" + char(10)
            if (g == 3 & x3 != "") {
                s = s + x3
            }
        }
    }
    return(s + tie + char(10))
}

struct dwlsmv_stage1 vector ttt_st(real rowvector K1)
{
    struct dwlsmv_stage1 vector st
    struct dwlsmv_statmap scalar m
    real scalar g, i

    st = dwlsmv_stage1(1, cols(K1))
    for (g = 1; g <= cols(K1); g++) {
        m         = dwlsmv_statmap_build(K1[g] \ J(3, 1, 3))
        st[g].map = m
        st[g].s   = J(m.D, 1, 0.3)
        for (i = 1; i <= m.T; i++) {
            st[g].s[i] = invnormal(m.thr_k[i] / m.K[m.var_i[i]]) + 0.02 * m.var_i[i] + 0.05 * (g - 1)
        }
    }
    return(st)
}

struct dwlsmv_parmap scalar ttt_id(real rowvector K1, string scalar tie)
{
    struct dwlsmv_stage1 vector st
    struct dwlsmv_statmap vector m
    real scalar g

    st = ttt_st(K1)
    m  = dwlsmv_statmap(1, cols(K1))
    for (g = 1; g <= cols(K1); g++) {
        m[g] = st[g].map
    }
    return(dwls_identify(ttt_text(cols(K1), tie), m, st))
}

// every member (g in gs) of y1's threshold k: start == want, one shared fidx
real scalar ttt_check(string scalar lab, struct dwlsmv_parmap scalar p,
                      real rowvector gs, real scalar k, real scalar want)
{
    real scalar i, n, f, bad

    n   = 0
    f   = .
    bad = 0
    for (i = 1; i <= p.nall; i++) {
        if (p.ptype[i] == 4 & p.row[i] == 1 & p.col[i] == k & anyof(gs, p.grp[i])) {
            n++
            if (f == .) {
                f = p.fidx[i]
            }
            if (p.fidx[i] != f | abs(p.value[i] - want) > 1e-14) {
                bad = 1
            }
            printf("  %s g%f start %21.15f fidx %f\n", lab, p.grp[i], p.value[i], p.fidx[i])
        }
    }
    if (n != cols(gs)) {
        bad = 1
    }
    printf("%s: want %21.15f %s\n", lab, want, (bad ? "FAIL" : "ok"))
    return(bad)
}

void ttt_starts()
{
    real scalar bad

    bad = 0
    bad = bad + ttt_check("T", ttt_id((2, 3, 3), "tie tau 2 y1 2 tau 3 y1 2"), (2, 3), 2, 0.500727299295457)
    bad = bad + ttt_check("silent", ttt_id((3, 4, 4), "tie tau 2 y1 2 tau 3 y1 2"), (2, 3), 2, 0.07)
    bad = bad + ttt_check("control", ttt_id((3, 3, 3), "tie tau 2 y1 2 tau 3 y1 2"), (2, 3), 2, 0.500727299295457)
    bad = bad + ttt_check("reversed", ttt_id((3, 3, 3, 3), "tie tau 4 y1 2 tau 3 y1 2"), (3, 4), 2, 0.550727299295457)
    bad = bad + ttt_check("incl", ttt_id((3, 3, 3), "tie tau 1 y1 2 tau 2 y1 2 tau 3 y1 2"), (1, 2, 3), 2, 0.450727299295457)
    if (bad) {
        _error(3498, "ttt_starts: a start differs")
    }
}

void ttt_mismatch()
{
    struct dwlsmv_parmap scalar p

    p = ttt_id((3, 3, 4), "tie tau 2 y1 2 tau 3 y1 2")
}

real matrix ttt_data(real scalar n, real rowvector cuts)
{
    real matrix Y
    real colvector f, ys
    real scalar j, c

    f = rnormal(n, 1, 0, 1)
    Y = J(n, 4, 0)
    for (j = 1; j <= 4; j++) {
        ys = 0.7 * f + sqrt(0.51) * rnormal(n, 1, 0, 1)
        if (j == 1) {
            for (c = 1; c <= cols(cuts); c++) {
                Y[., j] = Y[., j] + (ys :> cuts[c])
            }
        }
        else {
            Y[., j] = (ys :> -0.5) + (ys :> 0.5)
        }
    }
    return(Y)
}

void ttt_fit()
{
    struct dwlsmv_stage1 vector st
    struct dwlsmv_statmap vector m
    struct dwlsmv_parmap scalar p, q
    struct dwlsmv_fit scalar a, b
    real rowvector Kv
    real scalar g, i, d
    string scalar tie

    rseed(20261101)
    Kv = (3, 4, 4)
    st = dwlsmv_stage1(3, 1)
    m  = dwlsmv_statmap(3, 1)
    st[1] = dwls_sample_stats(ttt_data(1500, (-0.5, 0.5)), 3 \ 3 \ 3 \ 3)
    st[2] = dwls_sample_stats(ttt_data(1500, (-0.8, 0, 0.8)), 4 \ 3 \ 3 \ 3)
    st[3] = dwls_sample_stats(ttt_data(1500, (-0.8, 0, 0.8)), 4 \ 3 \ 3 \ 3)
    for (g = 1; g <= 3; g++) {
        m[g] = st[g].map
    }
    tie = "tie tau 2 y1 1 tau 3 y1 1" + char(10) + "tie tau 2 y1 2 tau 3 y1 2" + char(10) + "tie tau 2 y1 3 tau 3 y1 3"
    p = dwls_identify(ttt_text(3, tie, "nu y1 free" + char(10) + "omega y1 free" + char(10)), m, st)
    q = p
    for (i = 1; i <= q.nall; i++) {
        if (q.ptype[i] == 4 & q.row[i] == 1 & q.col[i] <= 2 & q.grp[i] >= 2) {
            printf("  rule start g%f %12.8f; former start %12.8f\n", q.grp[i], q.value[i], st[1].s[q.col[i]])
            q.value[i] = st[1].s[q.col[i]]
        }
    }
    a = dwls_fit(st, p)
    b = dwls_fit(st, q)
    d = max(abs(a.theta - b.theta))
    printf("fit: converged %f %f; iterations %f %f; fmin %21.15g %21.15g; max theta diff %9.2e\n",
           a.converged, b.converged, a.iterations, b.iterations, a.fmin, b.fmin, d)
    if (a.converged != 1 | b.converged != 1 | d > 1e-5) {
        _error(3498, "ttt_fit: the two starts disagree")
    }
}
end

local nfail 0
capture noisily mata: ttt_starts()
display "starts rc = " _rc
if _rc {
    local ++nfail
}
capture noisily mata: ttt_mismatch()
display "mismatch rc = " _rc
if _rc != 3498 {
    local ++nfail
}
capture noisily mata: ttt_fit()
display "fit rc = " _rc
if _rc {
    local ++nfail
}
if `nfail' {
    display as error "test_identify_tau_tie: `nfail' failure(s)"
    exit 9
}
display "test_identify_tau_tie: pass"
