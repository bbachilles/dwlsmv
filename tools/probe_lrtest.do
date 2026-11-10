* Probe (2026-10-10): does lrtest refuse a dwlsmv-style posting (b, V, N, T3; no e(ll))?
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
program define fakefit, eclass
    args t3
    matrix bb = (1, .5)
    matrix VV = (.01, 0 \ 0, .01)
    matrix colnames bb = y1:F y2:F
    matrix colnames VV = y1:F y2:F
    matrix rownames VV = y1:F y2:F
    ereturn post bb VV
    ereturn scalar N = 500
    ereturn scalar T3 = `t3'
    ereturn scalar df_T3 = 4
    ereturn scalar rank = 2
    ereturn local cmd "dwlsmv"
end
fakefit 10
estimates store m1
fakefit 12
estimates store m0
capture noisily lrtest m1 m0
display "lrtest rc = " _rc
