* Probe (2026-10-10): does a non-sem keyword name /scale(y1)#1bn.grp post and resolve in _b[] and test,
* beside sem-style names (spec 4.8, names of nu, omega, alpha at G > 1)?
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
matrix bb = (1, 1.2, 0, 0.3)
matrix VV = (0, 0, 0, 0 \ 0, .04, 0, 0 \ 0, 0, 0, 0 \ 0, 0, 0, .01)
local s "/scale(y1)#1bn.grp /scale(y1)#2.grp y1:1bn.grp y1:2.grp"
matrix colnames bb = `s'
matrix colnames VV = `s'
matrix rownames VV = `s'
capture noisily ereturn post bb VV
display "post rc = " _rc
capture noisily display _b[/scale(y1)#2.grp]
display "_b scale rc = " _rc
capture noisily display _b[y1:2.grp]
display "_b nu rc = " _rc
capture noisily test _b[/scale(y1)#2.grp] = 1
display "test rc = " _rc
capture noisily ereturn display
display "display rc = " _rc
local names : colfullnames e(b)
display "`names'"
