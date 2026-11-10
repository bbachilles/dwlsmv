*! The loading contract: every definition file can be loaded into one session,
*! in any order, without any of them wiping the others.
*!
*! This exists because they could not. Both src/mata/dwlsmv_types.mata and
*! tests/dwlsmv_assert.mata used to issue `mata clear`, so whichever loaded
*! second silently destroyed the first. Nothing failed loudly -- the structs
*! or the comparator were simply gone, and the error surfaced later as
*! "function not found" in a test that looked unrelated. Module 4's Gamma test
*! needs the structs AND the tolerance comparator together, so this had to be
*! fixed before it was written rather than discovered inside it.
*!
*! The fix: no definition file clears. Clearing belongs to the caller, and
*! every test opens with `clear all`, which clears Mata too. That was verified
*! rather than assumed -- a Mata function defined before `clear all` is gone
*! after it, raising rc 3499.
*!
*! A passing load is not enough to prove anything, because a wiped definition
*! looks exactly like one never loaded until it is called. So this test CALLS
*! something from every file it loads. If a `mata clear` is ever reintroduced,
*! one of those calls fails and this test fails with it.
*!
*! THE FILE LIST IS READ FROM DISK, not written out. Until 2026-09-26 it was a
*! hand list that stopped at Module 4, so nothing from Modules 3 or 5-10 was
*! covered, and its summary line claimed 20 files while loading 27. Every
*! src/mata file except dwlsmv_types.mata defines a function of its own name,
*! so after everything is loaded, each of those names must still resolve
*! (check 2 below). A new file is covered on the day it is written. The
*! per-module calls stay: a surviving name is not a working function.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

* Deliberately load the comparator FIRST and the types SECOND. That is the
* order that used to break: dwlsmv_types.mata cleared, so it destroyed the
* comparator loaded before it. The reverse order was the one that happened to
* work, which is why the problem went unnoticed. Everything else follows in
* directory order; only struct declarations must precede their use, and those
* are all in dwlsmv_types.mata.
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local srcs : list sort srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load tests/dwlsmv_assert.mata src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
local load `load' tests/dwlsmv_scan.mata
local nload : word count `load'

foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_loading: `f' failed to load"
        exit 9
    }
}

* ---- check 2: every file's own function survived every later load ----------
* The probe must be able to fail, or a capture that always returns 0 would
* pass everything. A name that was never defined is the negative control.
capture quietly mata: mata describe dwlsmv_no_such_function()
if !_rc {
    display as error "test_loading: the existence probe found a function that does not exist"
    exit 9
}
local nsurv 0
foreach f of local srcs {
    local fn = subinstr("`f'", ".mata", "", 1)
    capture quietly mata: mata describe `fn'()
    if _rc {
        display as error "test_loading: `fn'() is gone after loading every file"
        exit 9
    }
    local ++nsurv
}
display as text "test_loading: `nsurv' src/mata function(s) survive a full load"

mata:

real scalar dwlsmv_t_loading()
{
    struct dwlsmv_statmap scalar m
    struct dwlsmv_pairfit scalar pf
    struct dwlsmv_stage1 scalar s1
    real scalar nfail, gap
    real matrix YL, ST
    real rowvector worst, z
    string colvector names

    nfail = 0

    // ---- from tests/dwlsmv_assert.mata ------------------------------------
    // Loaded first, so this is the call that used to come back "not found".
    gap = dwlsmv_gapratio((1, 2), (1, 2), 1e-12, 1e-12)
    if (gap != 0) {
        printf("    dwlsmv_gapratio on identical input returned %21x, expected 0\n", gap)
        nfail = nfail + 1
    }
    worst = dwlsmv_worstcell((1, 2), (1, 2), 1e-12, 1e-12)
    if (cols(worst) != 5) {
        printf("    dwlsmv_worstcell returned %f columns, expected 5\n", cols(worst))
        nfail = nfail + 1
    }
    if (dwlsmv_symmetry((1, 2 \ 2, 1)) != 0) {
        printf("    dwlsmv_symmetry called a symmetric matrix asymmetric\n")
        nfail = nfail + 1
    }
    if (dwlsmv_mineig((2, 0 \ 0, 3)) != 2) {
        printf("    dwlsmv_mineig returned the wrong smallest eigenvalue\n")
        nfail = nfail + 1
    }

    // ---- from src/mata/dwlsmv_types.mata and the Module 0-1 files ---------
    // A struct declared in one file, built in another, read here: if the type
    // declarations had been wiped this would not compile.
    m = dwlsmv_statmap_build((4 \ 4 \ 3 \ 3))
    if (m.D != 16) {
        printf("    statmap D is %f, expected 16\n", m.D)
        nfail = nfail + 1
    }
    if (dwlsmv_pos_tau(m, 2, 1) != 4) {
        printf("    pos_tau(2,1) is %f, expected 4\n", dwlsmv_pos_tau(m, 2, 1))
        nfail = nfail + 1
    }
    if (dwlsmv_pos_rho(m, 1, 2) != 11) {
        printf("    pos_rho(1,2) is %f, expected 11\n", dwlsmv_pos_rho(m, 1, 2))
        nfail = nfail + 1
    }
    z = dwlsmv_pos_inv(m, 11)
    if (z[1] != 2 | z[2] != 1 | z[3] != 2) {
        printf("    pos_inv(11) did not return the (1,2) correlation\n")
        nfail = nfail + 1
    }

    // ---- from the Module 1 files ------------------------------------------
    // Hand-built counts, so the right answer is obvious by inspection: four
    // observations, categories 0 and 1 seen twice each, so q = 0.5 and the
    // threshold is 0 exactly.
    if (dwlsmv_catcounts((0 \ 0 \ 1 \ 1), 2) != (2 \ 2)) {
        printf("    dwlsmv_catcounts did not return (2,2) on a balanced binary column\n")
        nfail = nfail + 1
    }
    if (dwlsmv_catdegen((2 \ 0 \ 2)) != 2) {
        printf("    dwlsmv_catdegen did not classify an empty interior category as 2\n")
        nfail = nfail + 1
    }
    if (dwlsmv_cumprop((2 \ 2)) != 0.5) {
        printf("    dwlsmv_cumprop did not return 0.5 on (2,2)\n")
        nfail = nfail + 1
    }
    if (dwlsmv_thresholds((0 \ 0 \ 1 \ 1), 2) != 0) {
        printf("    dwlsmv_thresholds did not return 0 at q = 0.5\n")
        nfail = nfail + 1
    }
    if (rows(dwlsmv_psi_tau_table((0.5))) != 2) {
        printf("    dwlsmv_psi_tau_table returned the wrong number of rows\n")
        nfail = nfail + 1
    }
    if (rows(dwlsmv_listwise((0, 1 \ ., 1 \ 1, 0))) != 2) {
        printf("    dwlsmv_listwise did not drop the one incomplete row\n")
        nfail = nfail + 1
    }

    // ---- from the Module 2 files ------------------------------------------
    // Independence is the case with an obvious answer: at rho = 0 the cell
    // probabilities are the product of the margins, and with thresholds at 0
    // every cell of a 2x2 is exactly 1/4.
    if (abs(dwlsmv_cellprob((0), (0), 0)[1, 1] - 0.25) > 1e-15) {
        printf("    dwlsmv_cellprob did not give 1/4 for an independent 2x2 at tau = 0\n")
        nfail = nfail + 1
    }
    // At rho = 0 and tau = 0 the density corner is 1/(2 pi), and the rectangle
    // difference puts it in all four cells with alternating sign.
    if (abs(dwlsmv_dcellprob_rho((0), (0), 0)[1, 1] - 1 / (2 * pi())) > 1e-15) {
        printf("    dwlsmv_dcellprob_rho did not give 1/(2 pi) at tau = 0, rho = 0\n")
        nfail = nfail + 1
    }
    if (rows(dwlsmv_wmat((0 \ 1), (0), 0.3)) != 2) {
        printf("    dwlsmv_wmat returned the wrong number of rows\n")
        nfail = nfail + 1
    }
    if (dwlsmv_crosstab((0 \ 0 \ 1), (0 \ 1 \ 1), 2, 2) != (1, 1 \ 0, 1)) {
        printf("    dwlsmv_crosstab did not tabulate a hand-built pair correctly\n")
        nfail = nfail + 1
    }
    // dwls_polychoric returns the struct declared in dwlsmv_types.mata, so
    // this is the call that exercises the two files together. Eight
    // observations, all four cells filled, both margins balanced.
    pf = dwls_polychoric((0 \ 0 \ 1 \ 1 \ 0 \ 1 \ 0 \ 1),
                         (0 \ 1 \ 0 \ 1 \ 0 \ 0 \ 1 \ 1),
                         (0), (0))
    if (rows(pf.u) != 2 | cols(pf.u) != 2 | pf.I_rr <= 0) {
        printf("    dwls_polychoric returned a malformed pairfit\n")
        nfail = nfail + 1
    }

    // ---- from the Module 4 files ------------------------------------------
    // Four variables of two categories each: Koff runs 0, 2, 4, 6 and the
    // stack is 8 x 4, which is checkable by inspection rather than by running
    // anything.
    m = dwlsmv_statmap_build((2 \ 2 \ 2 \ 2))
    if (dwlsmv_koff(m) != (0 \ 2 \ 4 \ 6)) {
        printf("    dwlsmv_koff did not return the stacked row offsets\n")
        nfail = nfail + 1
    }
    YL = (0, 0, 0, 0 \ 0, 1, 1, 0 \ 1, 0, 1, 1 \ 1, 1, 0, 1 \
          0, 1, 0, 1 \ 1, 0, 0, 0 \ 1, 1, 1, 1 \ 0, 0, 1, 1)
    ST = dwlsmv_psi_tau_stack(YL, m)
    if (rows(ST) != 8 | cols(ST) != 4) {
        printf("    dwlsmv_psi_tau_stack returned %f x %f, expected 8 x 4\n",
               rows(ST), cols(ST))
        nfail = nfail + 1
    }
    // dwlsmv_pairinfo fills the struct declared in dwlsmv_types.mata, and
    // dwlsmv_psi_rho_table reads it -- two files plus the type boundary.
    dwlsmv_pairinfo(pf, (0), (0), 0, J(2, 2, 0.25))
    if (abs(pf.I_rr - 1 / (pi() * pi() * 0.25)) > 1e-12) {
        printf("    dwlsmv_pairinfo gave I_rr = %g on the independent 2x2\n", pf.I_rr)
        nfail = nfail + 1
    }
    if (rows(dwlsmv_psi_rho_table(pf, dwlsmv_psi_tau_table((0.5)),
                                  dwlsmv_psi_tau_table((0.5)))) != 2) {
        printf("    dwlsmv_psi_rho_table returned the wrong number of rows\n")
        nfail = nfail + 1
    }
    // dwls_sample_stats reaches every Module 4 file plus all of Modules 1-2,
    // and returns the struct they all share.
    s1 = dwls_sample_stats(YL, (2 \ 2 \ 2 \ 2))
    if (rows(s1.Gamma) != 10 | cols(s1.Gamma) != 10 | rows(s1.psibar) != 10) {
        printf("    dwls_sample_stats returned a malformed stage1 struct\n")
        nfail = nfail + 1
    }
    if (dwlsmv_symmetry(s1.Gamma) > 1e-12) {
        printf("    dwls_gamma returned an asymmetric Gamma\n")
        nfail = nfail + 1
    }

    // ---- from the diagnostic files ----------------------------------------
    if (dwlsmv_cellsparse((50, 60 \ 70, 0), 5)[1] != 2) {
        printf("    dwlsmv_cellsparse did not report an empty cell as code 2\n")
        nfail = nfail + 1
    }

    // ---- from tests/dwlsmv_scan.mata --------------------------------------
    names = dwlsmv_scanfile("src/mata/dwlsmv_pos_tau.mata")
    if (rows(names) == 0) {
        printf("    dwlsmv_scanfile returned nothing for a file with calls in it\n")
        nfail = nfail + 1
    }

    printf("test_loading: %f file(s) co-loaded, %f failure(s)\n",
           strtoreal(st_local("nload")), nfail)
    return(nfail)
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_loading()))

if `nfail' {
    display as error _n "test_loading: `nfail' call(s) into co-loaded files failed."
    display as error "A definition file has probably reintroduced `mata clear'."
    exit 9
}

display as result "test_loading: all definition files coexist and remain callable"
exit 0
