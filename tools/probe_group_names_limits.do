* Probe (2026-10-10): (1) gsem's cutpoint names under group() for ordinal outcomes;
* (2) the 32-character limit with group tags: which part of a grouped name it applies to.
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
set seed 2
set obs 2000
gen byte grp = 1 + (_n > 1000)
gen fz = rnormal()
forvalues j = 1/3 {
    gen z = 0.8*fz + rnormal()
    gen byte o`j' = 1 + (z > -0.5) + (z > 0.6)
    drop z
}
display "--- (1) gsem, ologit outcomes, group(grp)"
capture noisily gsem (F -> o1 o2 o3, ologit), group(grp) nolog
if !_rc {
    matrix b = e(b)
    local names : colfullnames b
    foreach n of local names {
        display "`n'"
    }
}
display "--- (2) sem, group(), a 32-character indicator name"
rename o1 abcdefghijklmnopqrstuvwxyz012345
capture noisily sem (F -> abcdefghijklmnopqrstuvwxyz012345 o2 o3), group(grp) nolog
if !_rc {
    matrix b = e(b)
    local names : colfullnames b
    foreach n of local names {
        display "`n'  [" strlen("`n'") "]"
    }
}
display "--- (3) ereturn post with long grouped stripes"
foreach s in "/var(e.abcdefghijklmnopqrstuvwxyz012345)#1bn.grp" "abcdefghijklmnopqrstuvwxyz012345:2.grp#c.F" "/abcdefghijklmnopqrstuvwxyz012345:cut1#2.grp" "/abcdefghijklmnopqrstuvwxyz012345:2.grp#cut1" {
    matrix bb = (1)
    matrix VV = (1)
    capture matrix colnames bb = `s'
    local rc1 = _rc
    capture matrix colnames VV = `s'
    capture matrix rownames VV = `s'
    capture ereturn post bb VV
    display "`s'  colnames rc=`rc1'  post rc=" _rc
    if !_rc {
        capture noisily display _b[`s']
        display "  _b[] rc=" _rc
    }
}
