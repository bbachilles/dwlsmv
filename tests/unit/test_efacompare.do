*! Drafted 2026-11-18 by a read-only agent; executed by the main session the same day, StataMP 18.
*! What broke on the first run: checks 7 and 12 (80 columns) failed only because the file set
*! linesize 120, so {p} paragraphs wrapped at 120; it now sets linesize 80. Nothing else changed.
*! Renamed from test_dwlsmv_nf.do with the command (dwlsmv_nf -> efacompare, owner 2026-11-18,
*! E-113); the checks are unchanged.
*!
*! efacompare (owner 2026-11-17 / 2026-11-18; E-111, review section 3; spec Module 12,
*! Surface): one row per number of factors m, each row the stored results of the standalone fit
*!     dwlsmv (F1..Fm -> items), efa(F1..Fm) [epsilon() rstarts()]
*! on the same sample. The deciding check is equality with that standalone fit, exact (== on the
*! columns that are copies, == on the two derived columns recomputed here with different code from
*! the ado's: ascending sort and a double loop, against the ado's descending sort and a diagonal
*! mask). No tolerance anywhere. Error paths assert the exact rc (198, 2000, 3498).
*!
*! Rows by status: 1 fitted; 2 the unrotated fit did not converge (rc 0, e(converged) 0): m and df
*! only, conv 0; 3 the geomin rotation did not converge from any start (rc 3498 AND the Mata external
*! dwlsmv_rotfail, owner 2026-11-18 option A): m and df only. Any other rc stops the command.
*!
*! End to end through build/ldwlsmv.mlib and src/ado, as test_efa_surface.do and test_unconv_post.do.
*! Private subprograms are reached the way test_efa_display.do reaches _dwlsmv_efaci: drop the
*! autoloaded program, `run' the file (every test below then runs that file's code).
*!
*! Checks:
*!   1  Case E (p = 6, ordinal), m = 1..3 (maxfactors(3) reaches the df = 0 row): every row against its
*!      standalone fit; r(Ev), r(N), r(p), r(maxfactors); the default range ends at m = 2; the
*!      caller's e() (a CFA fit made before) is untouched after the call
*!   2  Case G (p = 6, binary), m = 1..3: the same
*!   3  an if and an in subset on Case E: rows equal the standalone fits on the same subset; r(N)
*!   4  epsilon() and rstarts() reach every fit: rows equal the standalone fits with the same options
*!   5  the range at p = 4 (df(4,1) = 2, df(4,2) < 0): default m = 1 only; maxfactors(2) is rc 198
*!   6  Case F, N = 500 (keep in 1/500), m = 1..5: every row against its standalone fit; the m = 5
*!      row is a not-converged row (premise, asserted); df = 27 19 12 6 1; the display's notes
*!   7  Case F, N = 5000, full default run (m = 1..5), output logged: rc 0 past a failed m = 4 and
*!      m = 5; statuses 1 1 1 2 for m = 1..4, and m = 5 blank (status 2 or 3: rotation failure at 6
*!      processors, not converged at 4, measured 2026-11-18); the notes; no output line wider than
*!      80 columns; rows 1..3 against their standalone fits
*!   8  (b) Case F, N = 5000, m = 4: PRECONDITION asserted first, standalone: rc 0, e(converged) 0,
*!      e(rot_nose) "unconverged", no e(V), e(T3) missing, AND e(rot_minima), e(Lambda), e(Phi) posted
*!      (so that a blank row proves the command blanked them). The row at m = 4 is built and checked in
*!      9 (status 2, conv 0, df 6, every other column missing)
*!   9  9a: check 7's m = 5 row against its standalone fit, at whichever status it has. Then (c) the
*!      mixed sequence through the private row program, on Case F and Case E side by side: F m = 3,
*!      F m = 4, E m = 3, F m = 3. PRECONDITION asserted: standalone E m = 3 is rc 3498 and sets
*!      dwlsmv_rotfail, and the failed call leaves the F m = 4 fit in e() (the stale-result trap). The
*!      E m = 3 row is the rotation-failure row with nothing carried from F m = 4; the two F m = 3 rows
*!      are identical; the F rows equal the full run's rows of check 7. (Revised 2026-11-18, after
*!      measuring: F m = 5 was the rotation failure until a four-processor run made it non-converged.)
*!  10  the mark: dropped on entry to an estimation call (set by hand, then a converging fit), kept by
*!      a replay; NOT set by a different rc 3498 (Stage 1 refusal of a constant item); the row program
*!      at m = 2 on that constant item exits 3498 even with the mark set by hand beforehand (a mark
*!      from an earlier call is not read), and efacompare on it exits 3498 with the caller's e() restored
*!  11  error paths, exact rc: 198 (p < 4, an item named F1, a duplicated item, a string item,
*!      maxfactors() above the cap and below 1, a bad epsilon() and rstarts(), group(), standardized),
*!      2000 (empty sample), 3498 (a Stage 1 abort at m = 1 passes through)
*!  12  the display of a default run on Case E: the two blocks, no note for rows that do not exist, the
*!      one note that applies, and width
*!
*! Mutants planned (review.md section 3 (e)); each should turn a check red, run by the main session:
*!   M1  the default range off by one (ends at the largest m with df >= 0, or one early): 1, 5, 7
*!   M2  the df cap at > 1 (m = 1 only dropped, or df = 1 rows dropped): 5, 7
*!   M3  the second- or the fourth-largest |lambda| taken instead of the third: every full row (1, 2, 6)
*!   M4  max instead of min over the factors for lam3: 1, 2, 6
*!   M5  |phi| without abs (a negative factor correlation can be the largest in absolute value): 1, 2, 6
*!   M6  the diagonal included in max |phi|: 1, 2, 6
*!   M7  minima and gap swapped: 1, 2, 6
*!   M8  the admissibility flag inverted: 1, 2, 6 (comm1 is compared on every converged row)
*!   M9  `ereturn clear' removed: may SURVIVE, since a row is decided by rc and the mark, not by e();
*!       it is defence in depth. A survivor needs a demonstrated contract violation (memory rule)
*!  M10  the caller's e() not restored (the `restore' of _estimates hold removed): 1 and 10
*!  M11  `if' not passed to the inner call: 3, and every case whose N differs
*!  M12  the not-converged row filled from e(): 7, 8
*!  M13  the rotation-failure row filled from the previous row (stale e()): 7, 9
*!  M14  any rc 3498 at m >= 2 treated as a rotation failure (the mark not read): 10 (the constant
*!       item at m = 2 must exit 3498)
*!  M15  the mark not dropped on entry (dwlsmv, or the row program): 10
*!  M16  epsilon() / rstarts() not passed on: 4
*!  M17  the non-converged row shows T3 (or any fit column): 7, 8
*!  M18  the rotation-failure row reports conv 0 instead of missing: 7, 9
*!
*! Not certified: that any column identifies the true number of factors.
*! Not measured: whether the Case F failures at N = 5000 hold at 4 processors. The suite runs at both
*! counts, and checks 8 and 9 assert the premises, so a change shows as a failure, not a skip.
*! Time: Case F at N = 5000 fits m = 4 and m = 5 in about 11 and 14 seconds each (measured, one
*! machine); the file fits each of them a handful of times only (checks 7, 8, 9).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
* 80, the width the display checks (7, 12) assert: {p} paragraphs wrap at the line size
set linesize 80

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

global TNF_NF 0
* a tempfile, not a fixed name in c(tmpdir): two runs at once (mutants, parallel suites) would lock
* each other's log (r(608), seen 2026-11-18)
tempfile tnflg
global TNF_LG "`tnflg'"

* the private subprograms: drop the autoloaded command and run the file
capture program drop efacompare
capture program drop _efacompare_row
capture program drop _efacompare_show
quietly run src/ado/efacompare.ado

* independent recomputation of the two derived columns (different code from the ado's)
mata:
real scalar tnf_lam3(real matrix L)
{
    real scalar k, best
    real colvector x, v

    v = J(cols(L), 1, .)
    for (k = 1; k <= cols(L); k++) {
        x = sort(abs(L[., k]), 1)
        v[k] = x[rows(x) - 2]
    }
    best = min(v)
    return(best)
}

real scalar tnf_phimax(real matrix P)
{
    real scalar i, j, best

    if (rows(P) < 2) {
        return(.)
    }
    best = 0
    for (i = 1; i <= rows(P); i++) {
        for (j = 1; j <= cols(P); j++) {
            if (i != j) {
                best = max((best, abs(P[i, j])))
            }
        }
    }
    return(best)
}
end

capture program drop tnf_ok
program tnf_ok
    args cond label
    if !(`cond') {
        display as error "test_efacompare: `label'"
        global TNF_NF = $TNF_NF + 1
    }
end

* 1 if the mark exists
capture program drop tnf_mark
program tnf_mark, rclass
    mata: st_local("has", strofreal(findexternal("dwlsmv_rotfail") != NULL))
    return scalar has = `has'
end

* want (0/1): is the phrase in the log $TNF_LG (lines joined by one space)
capture program drop tnf_chk
program tnf_chk
    args want phrase label
    mata: st_local("got", strofreal(strpos(invtokens(cat(st_global("TNF_LG"))', " "), st_local("phrase")) > 0))
    if `got' != `want' {
        display as error "test_efacompare: `label' (phrase [`phrase'] found `got', expected `want')"
        global TNF_NF = $TNF_NF + 1
    }
end

* the longest line of the log
capture program drop tnf_width
program tnf_width, rclass
    mata: st_local("w", strofreal(max(strlen(cat(st_global("TNF_LG"))))))
    return scalar w = `w'
end

* a full row (name R) against the posted e() of the standalone fit; e() is a converged fit
capture program drop tnf_full
program tnf_full
    args R lab
    tempname E
    local m = el(`R', 1, 1)
    matrix `E' = e(Ev)
    tnf_ok "el(`R',1,2) == e(df_T3)"  "`lab': df"
    tnf_ok "el(`R',1,3) == 1"         "`lab': status"
    tnf_ok "el(`R',1,4) == e(converged)" "`lab': conv"
    tnf_ok "el(`R',1,5) == e(T3)"     "`lab': T3"
    tnf_ok "el(`R',1,6) == e(p_T3)"   "`lab': p"
    tnf_ok "el(`R',1,7) == e(rmsea_s)" "`lab': RMSEA"
    tnf_ok "el(`R',1,8) == e(cfi_s)"  "`lab': CFI"
    tnf_ok "el(`R',1,9) == e(tli_s)"  "`lab': TLI"
    tnf_ok "el(`R',1,10) == e(rot_minima)" "`lab': minima"
    tnf_ok "el(`R',1,11) == e(rot_gap)" "`lab': gap"
    mata: st_numscalar("tnf_l3", tnf_lam3(st_matrix("e(Lambda)"))); st_numscalar("tnf_ph", tnf_phimax(st_matrix("e(Phi)")))
    tnf_ok "el(`R',1,12) == scalar(tnf_l3)" "`lab': lam3"
    tnf_ok "el(`R',1,13) == scalar(tnf_ph)" "`lab': phimax"
    local c1 = ("`e(rot_nose)'" == "inadmissible")
    tnf_ok "el(`R',1,14) == `c1'"     "`lab': comm1"
    tnf_ok "el(`R',1,15) == el(`E',1,`m')" "`lab': ev_m"
end

* a blank row: status st (2 not converged, 3 rotation failed); only m and df are set
capture program drop tnf_blank
program tnf_blank
    args R lab st df
    tnf_ok "el(`R',1,2) == `df'" "`lab': df"
    tnf_ok "el(`R',1,3) == `st'" "`lab': status"
    if `st' == 2 {
        tnf_ok "el(`R',1,4) == 0" "`lab': conv is 0"
    }
    else {
        tnf_ok "missing(el(`R',1,4))" "`lab': conv is missing"
    }
    forvalues c = 5/15 {
        tnf_ok "missing(el(`R',1,`c'))" "`lab': column `c' is missing"
    }
end

* two rows equal element by element (missing equals missing)
capture program drop tnf_same
program tnf_same
    args A B lab
    forvalues c = 1/15 {
        tnf_ok "el(`A',1,`c') == el(`B',1,`c')" "`lab': column `c'"
    }
end

* every m = 1..mmax of the table TAB against its standalone fit. The standalone call is classified by
* its own rc: rc 0 and converged, rc 0 and not converged, rc 3498 with the mark; any other rc fails.
* e() is cleared first, so it is read only after rc 0.
capture program drop tnf_vs
program tnf_vs
    syntax varlist [if] [in], TAB(name) MMAX(integer) LAB(string) EV(name) N(integer) [OPTS(string)]
    local p : word count `varlist'
    tempname R E1
    forvalues m = 1/`mmax' {
        local fl ""
        forvalues k = 1/`m' {
            local fl `fl' F`k'
        }
        local fl : list retokenize fl
        matrix `R' = `tab'[`m', 1..15]
        ereturn clear
        capture mata: rmexternal("dwlsmv_rotfail")
        capture quietly dwlsmv (`fl' -> `varlist') `if' `in', efa(`fl') `opts'
        local rc = _rc
        tnf_mark
        local mk = r(has)
        local dfm = `p'*(`p'-1)/2 - (`p'*`m' - `m'*(`m'-1)/2)
        tnf_ok "el(`R',1,1) == `m'" "`lab' m=`m': the m column"
        if `rc' == 0 {
            tnf_ok "e(N) == `n'" "`lab' m=`m': N"
            matrix `E1' = e(Ev)
            tnf_ok "mreldif(`E1', `ev') == 0" "`lab' m=`m': r(Ev) equals e(Ev)"
            tnf_ok "e(df_T3) == `dfm'" "`lab' m=`m': e(df_T3) is the formula"
        }
        if `rc' == 0 & e(converged) == 1 {
            tnf_full `R' "`lab' m=`m' (fitted)"
        }
        else if `rc' == 0 {
            tnf_blank `R' "`lab' m=`m' (not converged)" 2 `dfm'
        }
        else if `rc' == 3498 & `mk' {
            tnf_blank `R' "`lab' m=`m' (rotation failed)" 3 `dfm'
        }
        else {
            tnf_ok "0" "`lab' m=`m': the standalone fit has rc `rc' and no mark, which efacompare would not have shown as a row"
        }
        display as text "  `lab' m=`m': status " el(`R',1,3) " (standalone rc `rc')"
    }
end


* ---- data ----
tempname TAB EV0 B0 TAB5 EV5 R3a R3b R4 R5 RX

* ---- 1  Case E, with a caller's e() ----
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), metric(unitvariance)
matrix `B0' = e(b)
local cmd0 `"`e(cmdline)'"'
local N0 = e(N)
quietly log using "$TNF_LG", text replace name(tnflog)
capture noisily efacompare y1-y6, maxfactors(3)
local rc = _rc
matrix `TAB' = r(table)
matrix `EV0' = r(Ev)
local NR = r(N)
local RP = r(p)
local RM = r(maxfactors)
quietly log close tnflog
tnf_ok "`rc' == 0" "1: the run is rc 0 (rc `rc')"
tnf_ok "`RP' == 6 & `RM' == 3" "1: r(p) and r(maxfactors)"
tnf_ok "rowsof(`TAB') == 3 & colsof(`TAB') == 15" "1: the table is 3 x 15"
tnf_ok "`NR' == 5000" "1: r(N)"
local same = ("`e(cmd)'" == "dwlsmv" & "`e(efa)'" == "" & `"`e(cmdline)'"' == `"`cmd0'"')
tnf_ok "`same' & e(N) == `N0'" "1: the caller's e() is back"
tnf_ok "mreldif(e(b), `B0') == 0" "1: the caller's e(b) is untouched"
tnf_chk 1 "df = 0 (m = 3)" "1: the df = 0 note names m = 3"
tnf_chk 1 "The default geomin epsilon changes with m" "1: the epsilon note (m up to 3)"
tnf_vs y1-y6, tab(`TAB') mmax(3) lab("E") ev(`EV0') n(`NR')
* the default range ends at the largest m with df >= 1: df(6,2) = 4, df(6,3) = 0
quietly efacompare y1-y6
local RM = r(maxfactors)
local nrow = rowsof(r(table))
tnf_ok "`RM' == 2 & `nrow' == 2" "1: the default range is m = 1, 2"

* ---- 2  Case G (binary) ----
import delimited using "fixtures/data/case_G.csv", clear varnames(1) delimiters(",")
quietly efacompare y1-y6, maxfactors(3)
matrix `TAB' = r(table)
matrix `EV0' = r(Ev)
local NR = r(N)
tnf_vs y1-y6, tab(`TAB') mmax(3) lab("G") ev(`EV0') n(`NR')

* ---- 3  if and in on Case E ----
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
generate byte sub = mod(_n, 3) != 0
quietly efacompare y1-y6 if sub in 1/4000
matrix `TAB' = r(table)
matrix `EV0' = r(Ev)
local NR = r(N)
quietly count if sub in 1/4000
local cnt = r(N)
tnf_ok "`NR' == `cnt' & `NR' < 4000" "3: r(N) is the marked sample (`NR' of `cnt')"
tnf_vs y1-y6 if sub in 1/4000, tab(`TAB') mmax(2) lab("E if/in") ev(`EV0') n(`NR')

* ---- 4  epsilon() and rstarts() ----
quietly efacompare y1-y6, maxfactors(2) epsilon(1e-3) rstarts(10)
matrix `TAB' = r(table)
matrix `EV0' = r(Ev)
local NR = r(N)
tnf_vs y1-y6, tab(`TAB') mmax(2) lab("E eps/starts") ev(`EV0') n(`NR') opts(epsilon(1e-3) rstarts(10))

* ---- 5  the range at p = 4 ----
quietly efacompare y1-y4
local RM = r(maxfactors)
local nrow = rowsof(r(table))
tnf_ok "`RM' == 1 & `nrow' == 1" "5: p = 4: the default is m = 1 only"
capture efacompare y1-y4, maxfactors(2)
local rc = _rc
tnf_ok "`rc' == 198" "5: p = 4, maxfactors(2) is rc 198 (rc `rc')"

* ---- 6  Case F, N = 500 ----
import delimited using "fixtures/data/case_F.csv", clear varnames(1) delimiters(",")
keep in 1/500
quietly log using "$TNF_LG", text replace name(tnflog)
capture noisily efacompare y1-y9
local rc = _rc
matrix `TAB' = r(table)
matrix `EV0' = r(Ev)
local NR = r(N)
local RM = r(maxfactors)
quietly log close tnflog
tnf_ok "`rc' == 0" "6: the run is rc 0 (rc `rc')"
tnf_ok "`NR' == 500 & `RM' == 5" "6: r(N), r(maxfactors)"
forvalues m = 1/5 {
    local want : word `m' of 27 19 12 6 1
    tnf_ok "el(`TAB',`m',2) == `want'" "6: df at m = `m'"
}
local st5 = el(`TAB', 5, 3)
tnf_ok "`st5' == 2" "6: the premise: the unrotated fit at m = 5 does not converge on N = 500 (status `st5')"
tnf_chk 1 "Not converged (m = 5)" "6: the not-converged note"
tnf_chk 0 "Rotation failed" "6: no rotation-failure note"
tnf_chk 1 "The default geomin epsilon changes with m" "6: the epsilon note"
tnf_chk 1 "Rotated quantities (minima" "6: the conditional-minimum note"
tnf_chk 0 "df = 0 (m" "6: no df = 0 note"
tnf_vs y1-y9, tab(`TAB') mmax(5) lab("F500") ev(`EV0') n(`NR')

* ---- 7  Case F, N = 5000: the full run ----
import delimited using "fixtures/data/case_F.csv", clear varnames(1) delimiters(",")
tnf_ok "_N == 5000" "7: Case F has 5000 rows"
quietly log using "$TNF_LG", text replace name(tnflog)
capture noisily efacompare y1-y9
local rc = _rc
matrix `TAB5' = r(table)
matrix `EV5' = r(Ev)
local NR = r(N)
local RP = r(p)
local RM = r(maxfactors)
quietly log close tnflog
tnf_ok "`rc' == 0" "7: the run continues past the failed rows and is rc 0 (rc `rc')"
tnf_ok "`NR' == 5000 & `RM' == 5 & `RP' == 9" "7: r(N), r(maxfactors), r(p)"
* m = 1..4 hold at 6 and at 4 processors (measured 2026-11-18). m = 5 depends on the processor
* count: at 6 the rotation fails, at 4 the unrotated fit does not converge (measured the same day,
* after the first four-processor suite run). Either is a blank row; which one it must be is
* decided against the standalone fit in check 9.
forvalues m = 1/4 {
    local want : word `m' of 1 1 1 2
    tnf_ok "el(`TAB5',`m',3) == `want'" "7: status at m = `m' (premise, measured 2026-11-18)"
}
forvalues m = 1/5 {
    local wdf : word `m' of 27 19 12 6 1
    tnf_ok "el(`TAB5',`m',2) == `wdf'" "7: df at m = `m'"
}
local st5 = el(`TAB5', 5, 3)
tnf_ok "`st5' == 2 | `st5' == 3" "7: m = 5 is a blank row, not converged or rotation failed (status `st5')"
tnf_chk 1 "Not converged (m = 4" "7: the not-converged note names m = 4"
if `st5' == 3 {
    tnf_chk 1 "Rotation failed (m = 5)" "7: the rotation-failure note names m = 5"
}
else {
    tnf_chk 0 "Rotation failed" "7: no rotation-failure note when m = 5 is not converged"
}
tnf_chk 0 "df = 0 (m" "7: no df = 0 note"
tnf_width
local W = r(w)
tnf_ok "`W' <= 80" "7: the output fits in 80 columns (longest line `W')"
* rows 1..3 against their standalone fits (cheap); rows 4 and 5 are checked in 8 and 9
tnf_vs y1-y9, tab(`TAB5') mmax(3) lab("F5000") ev(`EV5') n(`NR')

* ---- 8  (b) the real non-converged efa() path: PRECONDITION, standalone m = 4 ----
* If a later change makes this fit converge, the test FAILS here; it does not skip. The fixture
* would then have to be replaced.
ereturn clear
capture quietly dwlsmv (F1 F2 F3 F4 -> y1-y9), efa(F1 F2 F3 F4)
local rc = _rc
tnf_ok "`rc' == 0" "8: standalone m = 4 is rc 0 (rc `rc')"
tnf_ok "e(converged) == 0" "8: standalone m = 4 did not converge"
local nose = ("`e(rot_nose)'" == "unconverged")
tnf_ok "`nose'" "8: e(rot_nose) is unconverged"
capture confirm matrix e(V)
local rcv = _rc
tnf_ok "`rcv' == 111" "8: no e(V) (rc `rcv')"
tnf_ok "missing(e(T3))" "8: e(T3) is missing"
* what the row must blank is really posted, at the last iterate
capture confirm matrix e(Lambda)
local rcl = _rc
capture confirm matrix e(Phi)
local rcp = _rc
tnf_ok "`rcl' == 0 & `rcp' == 0 & e(rot_minima) < . & e(rot_minima) >= 1" "8: e(Lambda), e(Phi), e(rot_minima) are posted (blanking is not vacuous)"
tnf_ok "e(df_T3) == 6" "8: df 6"

* ---- 9  (c) the mixed sequence: success, not converged, rotation failure, success ----
* Revised 2026-11-18 after the first four-processor suite run (chosen after measuring): Case F at
* m = 5 fails to rotate at 6 processors but does not converge at 4, so it cannot carry the
* rotation failure. Case E at m = 3 fails to rotate at both (measured the same day). Case F and
* Case E are put side by side in one data set (both N = 5000), so the sequence runs with no reload:
* F m = 3 (fitted), F m = 4 (not converged), E m = 3 (rotation failed), F m = 3 (fitted).
*
* 9a. the m = 5 row of check 7 against its standalone fit, whichever status it has
ereturn clear
capture mata: rmexternal("dwlsmv_rotfail")
capture quietly dwlsmv (F1 F2 F3 F4 F5 -> y1-y9), efa(F1 F2 F3 F4 F5)
local rc = _rc
tnf_mark
local mk = r(has)
matrix `RX' = `TAB5'[5, 1..15]
if `rc' == 0 & e(converged) == 0 {
    tnf_blank `RX' "9a: F m = 5 (standalone not converged)" 2 1
}
else if `rc' == 3498 & `mk' {
    tnf_blank `RX' "9a: F m = 5 (standalone rotation failed)" 3 1
}
else {
    tnf_ok "0" "9a: F m = 5 standalone rc `rc', mark `mk': not a state check 7 allows"
}
display as text "  9a: F m = 5 standalone rc `rc', converged " e(converged) ", mark `mk'"

* the two cases side by side
tempfile fe
preserve
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
rename (y1 y2 y3 y4 y5 y6) (e1 e2 e3 e4 e5 e6)
quietly save `fe'
restore
quietly merge 1:1 _n using `fe', nogenerate
tnf_ok "_N == 5000" "9: Case F and Case E side by side, 5000 rows"

* PRECONDITION: standalone E m = 3 is rc 3498 and sets the mark, and the failed call leaves the
* previous fit (F m = 4) in e(): the stale-result trap. No ereturn clear here, on purpose.
capture quietly dwlsmv (F1 F2 F3 F4 -> y1-y9), efa(F1 F2 F3 F4)
capture quietly dwlsmv (F1 F2 F3 -> e1-e6), efa(F1 F2 F3)
local rc = _rc
tnf_mark
local mk = r(has)
tnf_ok "`rc' == 3498" "9: standalone E m = 3 is rc 3498 (rc `rc')"
tnf_ok "`mk' == 1" "9: the failed rotation sets dwlsmv_rotfail"
local stale = ("`e(efa)'" == "F1 F2 F3 F4")
tnf_ok "`stale' & e(df_T3) == 6" "9: the failed call leaves the F m = 4 fit in e()"

_efacompare_row y1-y9 if 1, m(3) matname(`R3a')
_efacompare_row y1-y9 if 1, m(4) matname(`R4')
_efacompare_row e1-e6 if 1, m(3) matname(`R5')
_efacompare_row y1-y9 if 1, m(3) matname(`R3b')
tnf_blank `R4' "9: the F m = 4 row" 2 6
tnf_blank `R5' "9: the E m = 3 row (nothing carried from F m = 4)" 3 0
tnf_same `R3a' `R3b' "9: the two F m = 3 rows"
tnf_ok "el(`R3a',1,3) == 1" "9: the F m = 3 row is a fitted row"
* each F row equals the full run's row of check 7 (rows 1..3 equal their standalone fits there)
matrix `RX' = `TAB5'[3, 1..15]
tnf_same `R3a' `RX' "9: F m = 3 against the full run"
matrix `RX' = `TAB5'[4, 1..15]
tnf_same `R4' `RX' "9: F m = 4 against the full run"

* ---- 10  the mark ----
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
* dropped on entry to an estimation call
mata: dwlsmv_rotfail = 1
quietly dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
tnf_mark
tnf_ok "r(has) == 0" "10: an estimation call drops a mark left by an earlier call"
* kept by a replay
mata: dwlsmv_rotfail = 1
quietly dwlsmv
tnf_mark
tnf_ok "r(has) == 1" "10: a replay leaves the mark alone"
capture mata: rmexternal("dwlsmv_rotfail")
* a different rc 3498 (Stage 1: a constant item) does not set the mark
generate byte cst = 1
capture quietly dwlsmv (F1 F2 -> y1 y2 y3 y4 y5 cst), efa(F1 F2)
local rc = _rc
tnf_mark
tnf_ok "`rc' == 3498" "10: a constant item is a Stage 1 refusal, rc 3498 (rc `rc')"
tnf_ok "r(has) == 0" "10: that rc 3498 does not set the mark"
* the row program at m = 2 does not read a mark left by an earlier call
mata: dwlsmv_rotfail = 1
capture quietly _efacompare_row y1 y2 y3 y4 y5 cst if 1, m(2) matname(`RX')
local rc = _rc
tnf_ok "`rc' == 3498" "10: the row program exits 3498 on a non-rotation rc 3498 (rc `rc')"
* efacompare on it: m = 1 stops with the rc, and the caller's e() is back
quietly dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), metric(unitvariance)
matrix `B0' = e(b)
local cmd0 `"`e(cmdline)'"'
capture efacompare y1 y2 y3 y4 y5 cst
local rc = _rc
tnf_ok "`rc' == 3498" "10: efacompare passes a Stage 1 rc 3498 at m = 1 through (rc `rc')"
local same = ("`e(cmd)'" == "dwlsmv" & "`e(efa)'" == "" & `"`e(cmdline)'"' == `"`cmd0'"')
tnf_ok "`same'" "10: the caller's e() is restored after the error exit"
tnf_ok "mreldif(e(b), `B0') == 0" "10: the caller's e(b) is untouched after the error exit"

* ---- 11  error paths, exact rc ----
capture efacompare y1 y2 y3
local rc = _rc
tnf_ok "`rc' == 198" "11: p = 3 is rc 198 (rc `rc')"
generate F1 = y1
capture efacompare F1 y2 y3 y4 y5
local rc = _rc
tnf_ok "`rc' == 198" "11: an item named F1 is rc 198 (rc `rc')"
drop F1
capture efacompare y1 y2 y3 y4 y5 y1
local rc = _rc
tnf_ok "`rc' == 198" "11: a duplicated item is rc 198 (rc `rc')"
generate str1 s = "a"
capture efacompare y1 y2 y3 y4 y5 s
local rc = _rc
tnf_ok "`rc' == 198" "11: a string item is rc 198 (rc `rc')"
drop s
capture efacompare y1-y6, maxfactors(4)
local rc = _rc
tnf_ok "`rc' == 198" "11: maxfactors(4) above the cap of 3 is rc 198 (rc `rc')"
capture efacompare y1-y6, maxfactors(0)
local rc = _rc
tnf_ok "`rc' == 198" "11: maxfactors(0) is rc 198 (rc `rc')"
capture efacompare y1-y6, epsilon(0)
local rc = _rc
tnf_ok "`rc' == 198" "11: epsilon(0) is rc 198 (rc `rc')"
capture efacompare y1-y6, epsilon(abc)
local rc = _rc
tnf_ok "`rc' == 198" "11: epsilon(abc) is rc 198 (rc `rc')"
capture efacompare y1-y6, rstarts(0)
local rc = _rc
tnf_ok "`rc' == 198" "11: rstarts(0) is rc 198 (rc `rc')"
capture efacompare y1-y6, rstarts(2.5)
local rc = _rc
tnf_ok "`rc' == 198" "11: rstarts(2.5) is rc 198 (rc `rc')"
generate byte g = 1 + (mod(_n, 2) == 0)
capture efacompare y1-y6, group(g)
local rc = _rc
tnf_ok "`rc' == 198" "11: group() is not an option, rc 198 (rc `rc')"
capture efacompare y1-y6, standardized
local rc = _rc
tnf_ok "`rc' == 198" "11: standardized is not an option, rc 198 (rc `rc')"
capture efacompare y1-y6 if y1 > 100
local rc = _rc
tnf_ok "`rc' == 2000" "11: an empty sample is rc 2000 (rc `rc')"

* ---- 12  the display of a default run on Case E ----
quietly log using "$TNF_LG", text replace name(tnflog)
capture noisily efacompare y1-y6
local rc = _rc
matrix `TAB' = r(table)
quietly log close tnflog
tnf_ok "`rc' == 0 & el(`TAB',1,3) == 1 & el(`TAB',2,3) == 1" "12: Case E m = 1, 2 are fitted rows (premise)"
tnf_chk 1 "Exploratory factor analysis by number of factors" "12: the title"
tnf_chk 1 "RMSEA" "12: block 1 header"
tnf_chk 1 "minima" "12: block 2 header"
tnf_chk 1 "Rotated quantities (minima" "12: the conditional-minimum note"
tnf_chk 0 "Not converged" "12: no not-converged note"
tnf_chk 0 "Rotation failed" "12: no rotation-failure note"
tnf_chk 0 "df = 0 (m" "12: no df = 0 note"
tnf_chk 0 "The default geomin epsilon changes with m" "12: no epsilon note (m up to 2)"
tnf_width
local W = r(w)
tnf_ok "`W' <= 80" "12: the output fits in 80 columns (longest line `W')"
* with epsilon() given the epsilon note goes
quietly log using "$TNF_LG", text replace name(tnflog)
capture noisily efacompare y1-y6, maxfactors(3) epsilon(1e-3)
quietly log close tnflog
tnf_chk 0 "The default geomin epsilon changes with m" "12: no epsilon note when epsilon() is given"

* ---- 13  a negative factor correlation (added by the main session 2026-11-18) ----
* Mutant M5 (|phi| without abs) survived checks 1-12: every fixture's factor correlation is
* positive. Case E with y4-y6 reverse-coded has phi_12 < 0 (asserted as the premise), so phimax
* must be its absolute value.
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
foreach v of varlist y4 y5 y6 {
    quietly summarize `v'
    quietly replace `v' = r(max) + r(min) - `v'
}
ereturn clear
capture quietly dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)
local rc = _rc
local ph21 = .
if `rc' == 0 local ph21 = el(e(Phi), 2, 1)
tnf_ok "`rc' == 0 & e(converged) == 1 & `ph21' < 0" "13: premise: Case E reversed, m = 2 converges with phi_12 < 0 (`ph21')"
quietly efacompare y1-y6
matrix `TAB' = r(table)
matrix `EV0' = r(Ev)
local NR = r(N)
tnf_ok "el(`TAB',2,13) == abs(`ph21') & el(`TAB',2,13) > 0" "13: phimax is |phi_12|"
tnf_vs y1-y6, tab(`TAB') mmax(2) lab("E reversed") ev(`EV0') n(`NR')

if $TNF_NF {
    display as error _n "test_efacompare: $TNF_NF failure(s)"
    exit 9
}
display as result "test_efacompare: every row of efacompare equals its standalone fit; non-converged and rotation-failed rows, the mark, the range and the error paths"
