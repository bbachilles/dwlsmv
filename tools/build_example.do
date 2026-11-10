*! dwlsmv -- build the help file's example datasets, build/dwlsmv_example.dta and
*! build/dwlsmv_example_mg.dta (the invariance example, 2026-10-12)
*!
*! Decided 2026-10-05 (the user's): the help file's examples use a
*! distributed dataset, dwlsmv_example.dta, not a test fixture. It is the
*! data of fixtures/data/case_E.csv (fixtures/generate/10_fixtures.R),
*! unchanged, with labels: six three-category indicators y1-y6, 5,000
*! observations, generated from two factors with loadings 0.7 and factor
*! correlation 0.5. fixtures/ is only read. The file goes to build/, like the
*! library, and is not committed.
*!
*! Run from the repository root:
*!     StataMP-64.exe /e do tools/build_example.do
*! test_dwlsmv_gaps.do runs it and checks that the fit on the dataset is the
*! fit on case_E.csv, bit for bit.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
preserve
import delimited using "fixtures/data/case_E.csv", ///
    clear varnames(1) delimiters(",") stringcols(_all)
quietly destring y1 y2 y3 y4 y5 y6, replace
confirm numeric variable y1 y2 y3 y4 y5 y6
assert _N == 5000
forvalues j = 1/6 {
    assert inlist(y`j', 0, 1, 2)
    label variable y`j' "Item `j' (0, 1, 2)"
}
label data "dwlsmv example: 6 ordinal items, 2 correlated factors (simulated)"
notes _dta: Simulated: y1-y3 load 0.7 on one factor, y4-y6 0.7 on another; factor correlation 0.5; 3 categories; N = 5,000. The data of the dwlsmv test fixture case_E.csv.
compress
capture mkdir build
quietly save "build/dwlsmv_example.dta", replace
restore
display "build_example: build/dwlsmv_example.dta"

* the invariance example (P5-0 (b), approved in principle by the user
* 2026-10-12, pending the packaging check): the data of the fixture
* mg_E.csv (fixtures/generate/98_lavaan_multigroup.R), unchanged except that
* the group variable g is named grp, with labels. fixtures/ is only read.
preserve
import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
confirm numeric variable g y1 y2 y3 y4 y5
assert _N == 4000
assert inlist(g, 1, 2)
forvalues j = 1/3 {
    assert inlist(y`j', 0, 1, 2, 3)
    label variable y`j' "Item `j' (0-3)"
}
forvalues j = 4/5 {
    assert inlist(y`j', 0, 1)
    label variable y`j' "Item `j' (0, 1)"
}
rename g grp
label variable grp "Group (1 = reference, 2)"
label data "dwlsmv invariance example: 5 ordinal items, 1 factor, 2 groups (simulated)"
notes _dta: Simulated: one factor; y1-y3 four categories, y4-y5 binary; group 1 (N = 1,200) is the reference, group 2 (N = 2,800) has factor mean 0.3, factor variance 1.44 and the same loadings and thresholds. The data of the dwlsmv test fixture mg_E.csv.
compress
quietly save "build/dwlsmv_example_mg.dta", replace
restore
display "build_example: build/dwlsmv_example_mg.dta"
