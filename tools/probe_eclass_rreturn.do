* Probe (2026-10-10): can an eclass command leave matrices in r() (via an rclass helper called last),
* so that `estadd matrix b_std = r(b_std)` works after `dwlsmv, standardized`?
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
program define stdhelper, rclass
    matrix S = (.6, .8)
    matrix colnames S = y1:F y2:F
    return matrix b_std = S
end
program define fakecmd, eclass
    matrix bb = (1, .5)
    matrix VV = (.01, 0 \ 0, .01)
    matrix colnames bb = y1:F y2:F
    matrix colnames VV = y1:F y2:F
    matrix rownames VV = y1:F y2:F
    ereturn post bb VV
    ereturn local cmd "fakecmd"
    stdhelper
end
fakecmd
capture noisily matrix list r(b_std)
display "r(b_std) rc = " _rc
display "e(cmd) = " e(cmd)
capture which estadd
display "estadd rc = " _rc
capture noisily estadd matrix b_std = r(b_std)
display "estadd rc = " _rc
capture noisily matrix list e(b_std)
display "e(b_std) rc = " _rc
