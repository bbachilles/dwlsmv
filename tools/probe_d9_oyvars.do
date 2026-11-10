*! dwlsmv -- Module 11 D9 probe: how sem orders e(oyvars)
*!
*! A probe, not a test (tools/). D9 (2026-10-04, the user's): the indicator order
*! is sem's e(oyvars) (D1), which fixes Stage 1's (Y, K) column order, the
*! section 4.2 index map and the /y:cutk block of e(b). This probe records, on the
*! machine it runs on, e(oyvars) after sem, noestimate for models whose path
*! order differs from the dataset order (y1..y6 created in that order):
*!   rev    (F1 -> y3 y1 y2)
*!   inter  (F1 -> y4 y2 y6) (F2 -> y1 y5 y3)
*!   struct (F1 -> y5 y6 y4) (F2 -> y3 y1 y2) (F2 <- F1)
*!   cross  (F1 -> y2 y1 y3 y4) (F2 -> y6 y4 y5)   y4 loads on both
*!   split  (F1 -> y6 y2) (F2 -> y4 y5) (F1 -> y1)  F1's paths in two pieces
*! Output lines prefixed D9 are the record. Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_d9_oyvars.do      (log: probe_d9_oyvars.log)
*!     grep -E "^(====|D9)" probe_d9_oyvars.log | diff - tools/probe_d9_oyvars_stata18.txt
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
clear all
set more off
set linesize 200
display "D9ENV stata_version " c(stata_version) " processors " c(processors)
set seed 20261004
set obs 200
generate double f = rnormal()
forvalues j = 1/6 {
    generate double y`j' = 0.7*f + rnormal()
}

local m_rev    "(F1 -> y3 y1 y2)"
local m_inter  "(F1 -> y4 y2 y6) (F2 -> y1 y5 y3)"
local m_struct "(F1 -> y5 y6 y4) (F2 -> y3 y1 y2) (F2 <- F1)"
local m_cross  "(F1 -> y2 y1 y3 y4) (F2 -> y6 y4 y5)"
local m_split  "(F1 -> y6 y2) (F2 -> y4 y5) (F1 -> y1)"
foreach m in rev inter struct cross split {
    display "==== `m'"
    capture noisily quietly sem `m_`m'', noestimate
    display "D9 `m' rc " _rc
    if _rc == 0 {
        display "D9 `m' oyvars " e(oyvars)
        display "D9 `m' lyvars " e(lyvars)
    }
}
display "==== done"
