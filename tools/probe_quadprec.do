*! dwlsmv -- probe: in what precision do quadsum() and quadcross() accumulate?
*!
*! Executed under StataMP 18 on 2026-10-14 (tools/probe_quadprec_stata18.txt).
*!
*! @sec-mg-sensitivity Result S3 bounds the processor-count dependence of
*! dwlsmv_fdwls's quadcross() reductions with the accumulator's unit roundoff
*! u_q. This probe sums (1, x, -1) with x = 2^-e for e = 50, 60, ..., 220,
*! in that order, and prints whether each routine returns x exactly. Plain sum()
*! (double) loses x for e >= 54. An 80-bit extended accumulator would lose it for
*! e >= 64 and IEEE binary128 for e >= 113; a double-double (pair) accumulator
*! keeps it at every e. The probe does not establish u_q itself: it rules out
*! the first two, which is all the section's bound (u_q <= 2^-104) uses.
*!
*! Run from anywhere:
*!     StataMP-64.exe /e do tools/probe_quadprec.do
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set linesize 200
display "QP|env|stata_version|`c(stata_version)'"
display "QP|env|processors|`c(processors)'"
mata:
for (e = 50; e <= 220; e = e + 10) {
    x = 2^-e
    v = (1 \ x \ -1)
    printf("QP|e=%3.0f|quadsum_exact|%f|quadcross_exact|%f|sum_exact|%f\n", e,
           quadsum(v) == x, quadcross(v, J(3, 1, 1)) == x, sum(v) == x)
}
end
