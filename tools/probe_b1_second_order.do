*! dwlsmv -- release B unit B1: what sem returns for a second-order model
*!
*! A probe, not a test (tools/); the suite runner ignores it. Roadmap
*! (E-108), B1; questions from
*! E-114 (architecture, section 2; recovery, 1.7).
*! Record: E-114 (the B1 finding).
*!
*! Runs sem, noestimate on a benign dummy dataset (12 four-category items
*! from one normal factor, N = 200, as Module 11's dummy data) for:
*!   P1  (G -> F1 F2 F3) under sem's default anchoring
*!   P2  the same with var(G@1)
*!   P3  forcenoanchor, no constraint (what Module 11's D2 pass sees)
*!   P4  forcenoanchor with var(G@1) and var(e.Fk@1) on all three Fk
*!   P5  the same constraints without forcenoanchor (does sem still anchor?)
*!   P6  mixed: G by var(G@1), the Fk by item markers (forcenoanchor)
*!   P7  mixed the other way: G by the path G -> F1 @1, the Fk by var(e.Fk@1)
*!   P8  var(F1@1) (total variance) on an endogenous F1: expected refused
*!   P9  X outside the hierarchy covarying with e.F1 e.F2 e.F3 (sem default
*!       adds cov(G*X)); P9b the same with cov(G*X@0); P9c X with sem's
*!       defaults only
*!   P10 an item also loading on G, written after the F's (P10) and before
*!       them (P10b): which path anchors G?
*! For each: rc, e(lxvars), e(lyvars), e(oyvars), and every e(Cns) row as
*! the stripe name it fixes and the value. Then Module 11 itself (the
*! built library): D2's injection on the forcenoanchor pass for P1's paths
*! and for P7's user-anchored G (does metric(unitvariance) express a mixed
*! metric?), and _dwlsmv_translate's final rc for P1, P2, P7 in each metric.
*!
*! Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_b1_second_order.do   (log: probe_b1_second_order.log)
*!     grep -E "^B1|^r\(" probe_b1_second_order.log | diff - tools/probe_b1_second_order_stata18.txt
*! Only the B1ENV line is expected to differ.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 220
display "B1ENV stata_version " c(stata_version) " processors " c(processors)

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

mata:
void b1_show(string scalar lab)
{
    real matrix C
    string matrix st
    real scalar r, j, nz, cj

    printf("B1LV %s lx=[%s] ly=[%s] oy=[%s]\n", lab, st_global("e(lxvars)"),
        st_global("e(lyvars)"), st_global("e(oyvars)"))
    st = st_matrixcolstripe("e(b)")
    C  = st_matrix("e(Cns)")
    printf("B1NB %s %f columns, %f Cns rows\n", lab, rows(st), rows(C))
    for (r = 1; r <= rows(C); r++) {
        nz = 0
        for (j = 1; j < cols(C); j++) {
            if (C[r, j] != 0) {
                nz++
                cj = j
            }
        }
        if (nz == 1) {
            printf("B1CNS %s %s:%s = %g\n", lab, st[cj, 1], st[cj, 2], C[r, cols(C)] / C[r, cj])
        }
        else {
            printf("B1CNS %s multi-column row (%f entries)\n", lab, nz)
        }
    }
}
void b1_stripe(string scalar lab)
{
    string matrix st
    real scalar j
    string scalar s

    st = st_matrixcolstripe("e(b)")
    for (j = 1; j <= rows(st); j++) {
        if (st[j, 2] != "_cons") printf("B1STR %s %f %s:%s\n", lab, j, st[j, 1], st[j, 2])
    }
}
void b1_printcanon(string scalar lab)
{
    string rowvector c
    real scalar i

    c = tokens(st_global("r(canon)"), ";")
    for (i = 1; i <= cols(c); i++) {
        if (c[i] != ";") printf("B1CAN %s %s\n", lab, c[i])
    }
}
void b1_canon(string scalar lab)
{
    string colvector P
    real scalar i

    P = dwlsmv_translate(0)
    for (i = 1; i <= rows(P); i++) printf("B1TXT %s %s\n", lab, P[i])
    printf("B1INJ %s [%s]\n", lab, dwlsmv_inject(P))
}
end

* benign dummy data, as Module 11's
set seed 20261119
quietly set obs 200
quietly generate double f = rnormal()
forvalues j = 1/12 {
    quietly generate double z = 0.6 * f + 0.8 * rnormal()
    quietly generate double y`j' = (z > -0.6) + (z > 0) + (z > 0.6)
    drop z
}
drop f

local M  (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9)
local SO (G -> F1 F2 F3) `M'
local FV var(e.F1@1) var(e.F2@1) var(e.F3@1)

capture program drop b1_run
program b1_run
    args lab paths opts
    capture noisily sem `paths', `opts' noestimate
    display "B1RC `lab' " _rc
    if _rc exit
    mata: b1_show("`lab'")
end

b1_run P1 "`SO'" ""
mata: b1_stripe("P1")
b1_run P2 "`SO'" "var(G@1)"
mata: b1_stripe("P2")
b1_run P3 "`SO'" "forcenoanchor"
b1_run P4 "`SO'" "var(G@1) `FV' forcenoanchor"
b1_run P5 "`SO'" "var(G@1) `FV'"
b1_run P6 "(G -> F1 F2 F3) (F1 -> y1@1 y2 y3) (F2 -> y4@1 y5 y6) (F3 -> y7@1 y8 y9)" "var(G@1) forcenoanchor"
b1_run P7 "(G -> F1@1 F2 F3) `M'" "`FV' forcenoanchor"
b1_run P8 "`SO'" "var(F1@1) forcenoanchor"
b1_run P9 "`SO' (X -> y10 y11 y12)" "cov(X*e.F1) cov(X*e.F2) cov(X*e.F3)"
b1_run P9b "`SO' (X -> y10 y11 y12)" "cov(X*e.F1) cov(X*e.F2) cov(X*e.F3) cov(G*X@0)"
* P9c: the outside factor's covariance with G alone, sem's default
b1_run P9c "`SO' (X -> y10 y11 y12)" ""
* P10: an item loading directly on G; which path does sem anchor G on?
b1_run P10 "(G -> F1 F2 F3 y10) `M'" ""
b1_run P10b "(G -> y10 F1 F2 F3) `M'" ""

* Module 11's D2 pass (the built library): what metric(unitvariance) injects
quietly sem `SO', forcenoanchor noestimate
mata: b1_canon("D2a")
quietly sem (G -> F1@1 F2 F3) `M', forcenoanchor noestimate
mata: b1_canon("D2b")
quietly sem (G -> F1 F2 F3) (F1 -> y1@1 y2 y3) (F2 -> y4@1 y5 y6) (F3 -> y7@1 y8 y9), var(G@1) forcenoanchor noestimate
mata: b1_canon("D2c")
* D2b's injection applied: is the result over-scaled (G fixed twice)?
capture noisily sem (G -> F1@1 F2 F3) `M', var(G@1) `FV' noestimate
display "B1RC D2b_applied " _rc
if !_rc mata: b1_show("D2b_applied")

* Module 11 end to end: final rc in each metric
foreach met in marker unitvariance {
    capture noisily _dwlsmv_translate , paths(`SO') metric(`met')
    display "B1M11 P1 `met' rc " _rc
    capture noisily _dwlsmv_translate , paths(`SO') semopts(var(G@1)) metric(`met')
    display "B1M11 P2 `met' rc " _rc
    capture noisily _dwlsmv_translate , paths((G -> F1@1 F2 F3) `M') metric(`met')
    display "B1M11 P7 `met' rc " _rc
}

* P10 through Module 11 and the estimator: does an item on G open the door?
foreach met in marker unitvariance {
    capture noisily _dwlsmv_translate , paths((G -> F1 F2 F3 y10) `M') metric(`met')
    display "B1M11 P10 `met' rc " _rc
    if !_rc {
        mata: b1_printcanon("P10 `met'")
    }
    capture noisily dwlsmv (G -> F1 F2 F3 y10) `M', metric(`met')
    display "B1FIT P10 `met' rc " _rc
    if !_rc {
        display "B1FIT P10 `met' converged " e(converged) " df " e(df_T3)
        matrix b = e(b)
        local nb = colsof(b)
        local cn : colfullnames b
        forvalues j = 1/`nb' {
            local n : word `j' of `cn'
            if strpos("`n'", "G") display "B1FITB P10 `met' `n' " %9.4f b[1, `j']
        }
    }
}

adopath - "`here'/build"
adopath - "`here'/src/ado"
display "B1END"
