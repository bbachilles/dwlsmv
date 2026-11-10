* Probe (2026-10-10): which group values can appear in sem-style names (n.grp#c.F), and what
* sem's group() accepts: non-integer, negative, zero, large, string.
* Run on StataMP 18 only; Stata 16 NOT YET EXECUTED (the milestone machine is the user's).
version 16
clear all
set seed 3
set obs 2000
gen fz = rnormal()
forvalues j = 1/3 {
    gen y`j' = 0.7*fz + rnormal()
}
gen half = 1 + (_n > 1000)
gen g_nonint = cond(half == 1, 1.5, 2.5)
gen g_neg    = cond(half == 1, -1, 2)
gen g_zero   = cond(half == 1, 0, 1)
gen g_big    = cond(half == 1, 1, 100000)
gen str5 g_str = cond(half == 1, "boys", "girls")
foreach v in g_nonint g_neg g_zero g_big g_str {
    display "--- sem, group(`v')"
    capture noisily quietly sem (F -> y1 y2 y3), group(`v')
    display "sem rc = " _rc
    if !_rc {
        matrix b = e(b)
        local names : colfullnames b
        display "first names: " word("`names'", 1) " " word("`names'", 2)
    }
}
display "--- hand-posted names"
foreach s in "y1:1.5.grp#c.F" "y1:-1.grp#c.F" "y1:0bn.grp#c.F" "y1:100000.grp#c.F" {
    matrix bb = (1)
    capture matrix colnames bb = `s'
    display "`s'  colnames rc=" _rc
}
