* Probe (2026-10-12, S4-4): why ereturn display labels a posted, fixed ν = 0
* coefficient (y1:1bn.g) "(omitted)" while /mean(F)#1bn.g = 0 shows
* "(constrained)", and whether any posting or display option shows the
* former as constrained WITHOUT renaming (the user's option (c); fallback (a)).
* Also how esttab renders the same columns, with and without noomitted.
* Diagnostic, not a pass/fail test. Executed under StataMP 18 on 2026-10-12
* (output alongside, *_stata18.txt); Stata 16 NOT YET EXECUTED.
*
* SPDX-License-Identifier: GPL-3.0-or-later
* SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
set more off
set linesize 255
set obs 10
generate byte g = 1 + (_n > 5)

* columns: nu g1 (fixed 0), nu g2 (fixed 0), lambda g1 (free), alpha g1
* (fixed 0), omega g1 (fixed 1), alpha g2 (free)
capture program drop mk
program mk, eclass
    args cnsopt
    matrix b = (0, 0, .7, 0, 1, .2)
    matrix colnames b = y1:1bn.g y1:2.g y1:1bn.g#c.F /:mean(F)#1bn.g /:scale(y1)#1bn.g /:mean(F)#2.g
    matrix V = J(6, 6, 0)
    matrix V[3, 3] = .01
    matrix V[6, 6] = .04
    matrix rownames V = y1:1bn.g y1:2.g y1:1bn.g#c.F /:mean(F)#1bn.g /:scale(y1)#1bn.g /:mean(F)#2.g
    matrix colnames V = y1:1bn.g y1:2.g y1:1bn.g#c.F /:mean(F)#1bn.g /:scale(y1)#1bn.g /:mean(F)#2.g
    matrix C = (1,0,0,0,0,0,0 \ 0,1,0,0,0,0,0 \ 0,0,0,1,0,0,0 \ 0,0,0,0,1,0,1)
    matrix colnames C = y1:1bn.g y1:2.g y1:1bn.g#c.F /:mean(F)#1bn.g /:scale(y1)#1bn.g /:mean(F)#2.g _Cns:_r
    if "`cnsopt'" == "nocns" ereturn post b V
    else if "`cnsopt'" == "fvinfo" ereturn post b V C, buildfvinfo
    else ereturn post b V C
    ereturn local cmd "probe"
end

display _n "==== V1 ereturn post b V C"
mk
ereturn display
_ms_omit_info e(b)
matrix list r(omit)

display _n "==== V2 ereturn post b V (no Cns)"
mk nocns
ereturn display

display _n "==== V3 ereturn post b V C, buildfvinfo"
capture noisily mk fvinfo
if !_rc ereturn display

display _n "==== V4 display options: noomitted, baselevels, allbaselevels"
mk
capture noisily ereturn display, noomitted
capture noisily ereturn display, allbaselevels

display _n "==== V5 b stripe with the o. operator removed explicitly (repost)"
mk
matrix bb = e(b)
local cn : colfullnames bb
display "colfullnames: `cn'"

display _n "==== E1 esttab default"
mk
estimates store P
capture noisily esttab P, se nostar
display _n "==== E2 esttab noomitted"
capture noisily esttab P, se nostar noomitted
display _n "==== E3 esttab omitted (explicit)"
capture noisily esttab P, se nostar omitted
display _n "==== E4 esttab nobaselevels"
capture noisily esttab P, se nostar nobaselevels
