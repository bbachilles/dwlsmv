* Probe (2026-10-10): sem's e(b) stripe under group(), for the G > 1 naming decision (spec 4.8 O-list, e() naming); is estout installed?
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
set seed 1
set obs 2000
gen byte grp = 1 + (_n > 1000)
gen fz = rnormal()
forvalues j = 1/3 {
    gen y`j' = 0.7*fz + rnormal()
}
gen gz = 0.5*fz + rnormal()
forvalues j = 4/6 {
    gen y`j' = 0.7*gz + rnormal()
}
capture noisily sem (F -> y1 y2 y3) (G -> y4 y5 y6) (F -> G), group(grp) nolog
matrix b = e(b)
local names : colfullnames b
foreach n of local names {
    display "`n'"
}
display "k = " colsof(b)
display "properties: " e(properties)
capture which esttab
display "esttab rc = " _rc
capture which estout
display "estout rc = " _rc
