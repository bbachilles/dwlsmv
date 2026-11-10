* Probe (2026-10-10): is partial() a sem option? (A dwlsmv option must not collide.)
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
set obs 200
set seed 5
gen fz = rnormal()
forvalues j = 1/3 {
    gen y`j' = 0.7*fz + rnormal()
}
capture noisily sem (F -> y1 y2 y3), partial(y3) noestimate
display "partial rc = " _rc
capture noisily sem (F -> y1 y2 y3), invariance(metric) noestimate
display "invariance rc = " _rc
