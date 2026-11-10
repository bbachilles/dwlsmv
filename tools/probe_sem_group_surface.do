* Probe (2026-10-10): sem's group grammar as dwlsmv's surface (D1). (1) group() + ginvariant(none) +
* noestimate; (2) e(Cns) for a cross-group @name tie; (3) group-specific paths.
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
set seed 4
set obs 400
gen byte grp = 1 + (_n > 200)
gen fz = rnormal()
forvalues j = 1/3 {
    gen y`j' = 0.7*fz + rnormal()
}
display "--- (1) ginvariant(none), noestimate"
capture noisily sem (F -> y1 y2 y3), group(grp) ginvariant(none) noestimate
display "rc = " _rc
if !_rc {
    local names : colfullnames e(b)
    display "`names'"
    capture matrix list e(Cns)
    display "Cns rc = " _rc
    if !_rc matrix list e(Cns)
}
display "--- (2) cross-group tie @a on y2's loading"
capture noisily sem (F -> y1 y2@a y3), group(grp) ginvariant(none) noestimate
display "rc = " _rc
if !_rc {
    local names : colfullnames e(b)
    display "`names'"
    matrix list e(Cns)
}
display "--- (3) group-specific tie: only group-specific paths"
capture noisily sem (1: F -> y1 y2@b y3) (2: F -> y1 y2@b y3), group(grp) ginvariant(none) noestimate
display "rc = " _rc
if !_rc matrix list e(Cns)
display "--- (4) default ginvariant (no option), noestimate"
capture noisily sem (F -> y1 y2 y3), group(grp) noestimate
display "rc = " _rc
if !_rc matrix list e(Cns)

* Extension (2026-10-12, decision "Group-specific paths under invariance()", option (a)):
* can Module 11 read each group's paths from e(b)'s stripe and e(Cns), so that a
* non-mirroring model is detectable from e() and a mirrored one is not flagged?
* Five indicators, two factors, so a cross-loading can differ by group.
* Run on StataMP 18; Stata 16 NOT YET EXECUTED.
capture drop y4 y5
gen y4 = 0.7*fz + rnormal()
gen y5 = 0.7*fz + rnormal()
* indexed locals: ": word of" splits compound-quoted lists (trap, 2026-10-08)
local m1 (F -> y1 y2 y3) (G -> y3 y4 y5)
local m2 (1: F -> y1 y2 y3) (2: F -> y1 y2 y3) (G -> y3 y4 y5)
local m3 (F -> y1 y2 y3) (G -> y4 y5) (2: G -> y3)
local m4 (F -> y1 y2 y3) (G -> y3 y4 y5) (2: F -> y4)
local m5 (F -> y1 y2 y3) (G -> y3 y4 y5) (2: F -> y4@0)
local m6 (F -> y1 y2 y3) (G -> y3 y4 y5) (2: F -> y4@0.3)
local m7 (F -> y1 y2 y3) (G -> y3 y4 y5) (2: y5 <- y4)
forvalues c = 1/7 {
    local m `m`c''
    display "--- (5.`c') `m', ginvariant(none) noestimate"
    capture noisily sem `m', group(grp) ginvariant(none) noestimate
    display "rc = " _rc
    if !_rc {
        local names : colfullnames e(b)
        display "ncols = " colsof(e(b))
        foreach n of local names {
            display "  `n'"
        }
        capture confirm matrix e(Cns)
        if !_rc matrix list e(Cns), noheader
        else display "no e(Cns)"
        matrix list e(b), noheader
    }
}
