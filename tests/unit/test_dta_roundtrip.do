*! The R-to-Stata transfer channel (spec Section 7.1 step 2).
*!
*! Section 7.1 writes each replicate to .dta because both tools read it
*! natively, and requires the round-trip be verified lossless BEFORE anything
*! downstream is trusted. A silent precision loss here would corrupt every
*! paired comparison in the validation design while looking exactly like an
*! estimator disagreement -- which is the most expensive kind of bug this
*! project can have.
*!
*! This is the half that counts. fixtures/generate/40_dta_roundtrip.R checks
*! that haven can read back what haven wrote, which proves only that haven
*! agrees with itself. The question is whether STATA reads what R wrote, and
*! only Stata can answer it.
*!
*! Two comparisons:
*!
*!   case_D          the real payload. The .dta is compared against the .csv,
*!                   both read by Stata, so the two routes into Stata are
*!                   independent of each other.
*!   roundtrip_probe doubles chosen to break a lossy channel -- values with no
*!                   finite binary expansion, 2^53 and its neighbour, and both
*!                   ends of Stata's representable range. The reference reaches
*!                   Stata as %.17g TEXT, which round-trips an IEEE double
*!                   exactly, so the .dta is checked against something that did
*!                   not travel by .dta.
*!
*! Known and deliberate: negative zero is compared numerically, not bitwise,
*! because -0 == 0 is true in IEEE arithmetic and Mata's == follows it. The
*! probe carries -0 so the channel is exercised on it, but this test cannot
*! distinguish -0 from 0 and does not claim to.
*!
*! This test has no module dependency and never legitimately skips.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

local D "fixtures/data"

* ---------------------------------------------------------------------------
* 1. the real payload: case_D.dta against case_D.csv, both read by Stata
* ---------------------------------------------------------------------------
capture noisily use "`D'/case_D.dta", clear
if _rc {
    display as error "test_dta_roundtrip: cannot read `D'/case_D.dta"
    display as error "Generate it: Rscript fixtures/generate/40_dta_roundtrip.R"
    exit 9
}
mata: DTA = st_data(., .)

capture noisily import delimited "`D'/case_D.csv", clear varnames(1)
if _rc {
    display as error "test_dta_roundtrip: cannot read `D'/case_D.csv"
    exit 9
}
mata: CSV = st_data(., .)

mata:
nfail = 0
if (rows(DTA) != rows(CSV) | cols(DTA) != cols(CSV)) {
    printf("    case_D: %f x %f from .dta, %f x %f from .csv\n",
           rows(DTA), cols(DTA), rows(CSV), cols(CSV))
    nfail = nfail + 1
}
else {
    if (DTA != CSV) {
        printf("    case_D: %f cell(s) differ between .dta and .csv\n",
               sum(DTA :!= CSV))
        nfail = nfail + 1
    }
}
st_local("nfail1", strofreal(nfail))
end

* ---------------------------------------------------------------------------
* 2. the precision probe: .dta against %.17g text
* ---------------------------------------------------------------------------
capture noisily use "`D'/roundtrip_probe.dta", clear
if _rc {
    display as error "test_dta_roundtrip: cannot read `D'/roundtrip_probe.dta"
    exit 9
}
mata: PDTA = st_data(., .)

capture noisily import delimited "`D'/roundtrip_probe.csv", clear varnames(1) stringcols(_all)
if _rc {
    display as error "test_dta_roundtrip: cannot read `D'/roundtrip_probe.csv"
    exit 9
}

mata:
PTXT = st_sdata(., .)
nfail = 0

if (rows(PDTA) != rows(PTXT) | cols(PDTA) != cols(PTXT)) {
    printf("    probe: shape differs between .dta and text reference\n")
    nfail = nfail + 1
}
else {
    for (i = 1; i <= rows(PDTA); i++) {
        for (j = 1; j <= cols(PDTA); j++) {
            want = strtoreal(PTXT[i, j])
            got  = PDTA[i, j]
            if (got != want) {
                printf("    probe[%f,%f]: .dta %21x  text %21x  (%s)\n",
                       i, j, got, want, PTXT[i, j])
                nfail = nfail + 1
            }
        }
    }
}

printf("test_dta_roundtrip: %f x %f payload, %f x %f probe, %f probe failure(s)\n",
       rows(DTA), cols(DTA), rows(PDTA), cols(PDTA), nfail)
st_local("nfail2", strofreal(nfail))
end

if `nfail1' | `nfail2' {
    display as error _n "test_dta_roundtrip: `nfail1' payload failure(s), " ///
                        "`nfail2' probe failure(s)."
    display as error "The R-to-Stata channel is NOT lossless. Do not trust any"
    display as error "paired comparison built on it until this is understood."
    exit 9
}

display as result "test_dta_roundtrip: .dta channel is lossless into Stata"
exit 0
