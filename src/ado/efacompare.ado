*! version 0.6.1  efacompare -- fit and over-extraction diagnostics by number of factors, one row per m
*!
*!   efacompare varlist [if] [in] [, maxfactors(#) epsilon(#) rstarts(#) keepbound]
*!
*! Drafted 2026-11-18 by a read-only agent; run by the main session the same day (StataMP 18,
*! tests/unit/test_dwlsmv_nf.do, every help example), with no change needed to this file.
*! Renamed from dwlsmv_nf before publication (owner 2026-11-18; the development record): the command,
*! its private subprograms, its help and its test (now tests/unit/test_efacompare.do). Behaviour
*! is unchanged.
*!
*! Design: the design package and its review in the development record; owner calls of 2026-11-17 and
*! 2026-11-18 (the development record); spec Module 12, Surface, "efacompare". It loops
*!     dwlsmv (F1..Fm -> varlist) if touse, efa(F1..Fm) [epsilon() rstarts() keepbound]
*! for m = 1, 2, ... under quietly and reads one row of results per m. Stage 1 is recomputed at every
*! m (no reuse), so a row equals the standalone fit by construction. No verdict, no cutoffs, no
*! difference test between m. The caller's e() is held and restored, so nothing is left in e().
*!
*! The row of one m is built by the private _efacompare_row, which decides the row from the return
*! code, never from e() alone (dwlsmv does not clear e() on entry, so after a failed call e() still
*! holds the previous fit):
*!   rc 0, e(converged) = 1        the full row (status 1)
*!   rc 0, e(converged) = 0        m, df, conv = 0, all else missing (status 2)
*!   rc 3498, m >= 2, and the Mata external dwlsmv_rotfail exists (dwlsmv sets it on the geomin
*!                                 "no converged start" exit only, owner 2026-11-18, option A):
*!                                 m, df, all else missing (status 3)
*!   any other rc                  the command stops with that rc (one line names m and the
*!                                 standalone command that shows dwlsmv's own message)
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

program efacompare, rclass
    version 16
    syntax varlist [if] [in] [, MAXfactors(integer -1) EPSilon(string) RSTARTS(string) KEEPBOUND]
    local items `varlist'
    local p : word count `items'

    * refusals, rc 198, before anything is run
    if `p' < 4 {
        display as error "efacompare needs at least 4 items: with fewer, even one factor leaves no degrees of freedom."
        exit 198
    }
    foreach v of local items {
        if regexm("`v'", "^F[0-9]+$") {
            display as error "item `v' is named like a factor: efacompare names its factors F1, F2, ..."
            display as error "Rename the item (for example with rename) and run the command again."
            exit 198
        }
    }
    local dups : list dups items
    if "`dups'" != "" {
        display as error "item `: word 1 of `dups'' is listed twice."
        exit 198
    }
    capture confirm numeric variable `items'
    if _rc {
        display as error "the items must be numeric variables."
        exit 198
    }
    if `"`epsilon'"' != "" {
        capture confirm number `epsilon'
        if _rc | real(`"`epsilon'"') <= 0 {
            display as error "epsilon() must be a positive number."
            exit 198
        }
    }
    if `"`rstarts'"' != "" {
        capture confirm integer number `rstarts'
        if _rc | real(`"`rstarts'"') < 1 {
            display as error "rstarts() must be a positive integer."
            exit 198
        }
    }

    * the range. df(p, m) = p(p-1)/2 - (p*m - m(m-1)/2), thresholds saturated (the formula
    * dwlsmv_m11_efa_canon.mata refuses a negative value of). df falls as m rises up to m = p, so the
    * caps are found by stepping m up; cap1 is the default's end, cap0 the most maxfactors() may reach
    local cap1 0
    local cap0 0
    forvalues k = 1/`p' {
        local dfk = `p'*(`p'-1)/2 - (`p'*`k' - `k'*(`k'-1)/2)
        if `dfk' >= 1 local cap1 = `k'
        if `dfk' >= 0 local cap0 = `k'
    }
    if `maxfactors' == -1 local mmax = `cap1'
    else {
        if `maxfactors' < 1 {
            display as error "maxfactors() must be at least 1."
            exit 198
        }
        if `maxfactors' > `cap0' {
            display as error "maxfactors(`maxfactors') is too large for `p' items: the largest number of factors"
            display as error "that leaves df >= 0 is `cap0'; the default range ends at `cap1' (df >= 1)."
            exit 198
        }
        local mmax = `maxfactors'
    }

    * the sample, marked once (items, if, in); each fit is given if `touse', so a row's N is the
    * standalone fit's N on the same if/in
    marksample touse
    quietly count if `touse'
    local N = r(N)
    if `N' == 0 error 2000

    local opts ""
    if `"`epsilon'"' != "" local opts `opts' epsilon(`epsilon')
    if `"`rstarts'"' != "" local opts `opts' rstarts(`rstarts')
    if "`keepbound'" != "" local opts `opts' keepbound

    * the caller's e() is held and restored when this program ends, on every exit path
    tempname eh tab row evm ev
    _estimates hold `eh', restore nullok

    local evset 0
    local rn ""
    display as text "Fitting m = 1, ..., `mmax' (one dot per m): " _continue
    forvalues m = 1/`mmax' {
        capture matrix drop `evm'
        _efacompare_row `items' if `touse', m(`m') matname(`row') evname(`evm') `opts'
        matrix `tab' = nullmat(`tab') \ `row'
        local rn `rn' m`m'
        if !`evset' {
            capture confirm matrix `evm'
            if !_rc {
                matrix `ev' = `evm'
                local evset 1
            }
        }
        display as text "." _continue
    }
    display as text ""
    matrix rownames `tab' = `rn'

    local epsg = (`"`epsilon'"' != "")
    _efacompare_show, matname(`tab') n(`N') p(`p') items(`items') mmax(`mmax') epsgiven(`epsg')

    return scalar N = `N'
    return scalar p = `p'
    return scalar maxfactors = `mmax'
    if `evset' return matrix Ev = `ev'
    return matrix table = `tab'
end

* one row for one m. Columns: m df status conv T3 p RMSEA CFI TLI minima gap lam3 phimax comm1 ev_m
* status 1 fitted, 2 not converged, 3 rotation failed. matname() receives the 1 x 15 row; evname(),
* if given, receives e(Ev) when dwlsmv returned rc 0. eclass because of ereturn clear; tests call it
* directly, in any order, after `program drop efacompare' and `run src/ado/efacompare.ado'
program _efacompare_row, eclass
    version 16
    syntax varlist [if], M(integer) MATname(name) [EVname(name) EPSilon(string) RSTARTS(string) KEEPBOUND]
    local p : word count `varlist'
    local fl ""
    forvalues k = 1/`m' {
        local fl `fl' F`k'
    }
    local fl : list retokenize fl
    local dfm = `p'*(`p'-1)/2 - (`p'*`m' - `m'*(`m'-1)/2)
    local opts ""
    if `"`epsilon'"' != "" local opts `opts' epsilon(`epsilon')
    if `"`rstarts'"' != "" local opts `opts' rstarts(`rstarts')
    if "`keepbound'" != "" local opts `opts' keepbound

    * dwlsmv never clears e() on entry: a failed call leaves the previous fit there. Clear it, so e()
    * is read only after rc 0. The mark is dropped here as well (dwlsmv drops it on entry too).
    ereturn clear
    capture mata: rmexternal("dwlsmv_rotfail")
    capture quietly dwlsmv (`fl' -> `varlist') `if', efa(`fl') `opts'
    local rc = _rc

    local rot 0
    if `rc' == 3498 & `m' >= 2 {
        mata: st_local("rot", strofreal(findexternal("dwlsmv_rotfail") != NULL))
    }

    local cn m df status conv T3 p RMSEA CFI TLI minima gap lam3 phimax comm1 ev_m
    if `rc' == 0 {
        if "`evname'" != "" matrix `evname' = e(Ev)
        if e(converged) == 1 {
            tempname L3 T3k PH EvM
            matrix `EvM' = e(Ev)
            * lam3: per factor column, the 3rd-largest absolute loading; the minimum over columns.
            * One Mata line per column, no subscript on a call result, no brace block.
            scalar `L3' = .
            forvalues k = 1/`m' {
                mata: dnfL = abs(st_matrix("e(Lambda)")); dnfx = sort(dnfL[., `k'], -1); st_numscalar("`T3k'", dnfx[3])
                scalar `L3' = min(`L3', `T3k')
            }
            * phimax: the largest absolute off-diagonal of Phi (missing at m = 1)
            mata: dnfP = abs(st_matrix("e(Phi)")); st_numscalar("`PH'", (cols(dnfP) > 1 ? max(dnfP - diag(diagonal(dnfP))) : .))
            capture mata: mata drop dnfL dnfx dnfP
            local comm1 = ("`e(rot_nose)'" == "inadmissible")
            matrix `matname' = (`m', e(df_T3), 1, 1, e(T3), e(p_T3), e(rmsea_s), e(cfi_s), e(tli_s), ///
                e(rot_minima), e(rot_gap), `L3', `PH', `comm1', el(`EvM', 1, `m'))
        }
        else {
            * the unrotated fit did not converge: the last iterate is not a solution, so nothing but
            * m and df is read, although e(rot_minima), e(Lambda), e(Phi) are posted
            matrix `matname' = (`m', e(df_T3), 2, 0, ., ., ., ., ., ., ., ., ., ., .)
        }
    }
    else if `rot' {
        * the geomin "no converged start" exit (the mark): nothing posted, e() is cleared
        matrix `matname' = (`m', `dfm', 3, ., ., ., ., ., ., ., ., ., ., ., .)
    }
    else {
        local hint ""
        if `rc' == 3498 & `m' >= 2 local hint " (not a rotation failure)"
        display as error _n "efacompare: dwlsmv stopped at m = `m' with rc `rc'`hint'. To see its message, run (with the same if/in): dwlsmv (`fl' -> `varlist'), efa(`fl')`opts'"
        exit `rc'
    }
    matrix colnames `matname' = `cn'
end

* the display of the table: two blocks of at most 68 columns, then the footnotes that apply
program _efacompare_show
    version 16
    syntax, MATname(name) N(integer) P(integer) ITEMS(string) MMAX(integer) EPSGiven(integer)
    local nr = rowsof(`matname')
    local l2 ""
    local l3 ""
    local l0 ""
    forvalues r = 1/`nr' {
        local mm = el(`matname', `r', 1)
        local st = el(`matname', `r', 3)
        if `st' == 2 local l2 `l2' `mm'
        if `st' == 3 local l3 `l3' `mm'
        if el(`matname', `r', 2) == 0 local l0 `l0' `mm'
    }

    display as text _n "Exploratory factor analysis by number of factors m (oblique geomin)"
    display as text "{p 0 6 2}Observations: `n'.  Items (p = `p'): `items'{p_end}"
    display as text "Columns are defined in {help efacompare}."

    display as text _n "Fit"
    display as text "{hline 68}"
    display as text _col(2) %3s "m" _col(7) %4s "df" _col(13) %-8s "status" _col(23) %9s "T3" ///
        _col(34) %7s "p" _col(43) %7s "RMSEA" _col(52) %7s "CFI" _col(61) %7s "TLI"
    display as text "{hline 68}"
    forvalues r = 1/`nr' {
        local st = el(`matname', `r', 3)
        local lab = cond(`st' == 1, "fitted", cond(`st' == 2, "not conv", "rot fail"))
        display as result _col(2) %3.0f el(`matname', `r', 1) _col(7) %4.0f el(`matname', `r', 2) ///
            _col(13) %-8s "`lab'" _col(23) %9.3f el(`matname', `r', 5) _col(34) %7.3f el(`matname', `r', 6) ///
            _col(43) %7.3f el(`matname', `r', 7) _col(52) %7.3f el(`matname', `r', 8) ///
            _col(61) %7.3f el(`matname', `r', 9)
    }
    display as text "{hline 68}"

    display as text _n "Rotated diagnostics and eigenvalue"
    display as text "{hline 68}"
    display as text _col(2) %3s "m" _col(7) %6s "minima" _col(15) %10s "gap" _col(27) %8s "lam3" ///
        _col(37) %8s "phimax" _col(47) %6s "comm1" _col(55) %9s "ev_m"
    display as text "{hline 68}"
    forvalues r = 1/`nr' {
        display as result _col(2) %3.0f el(`matname', `r', 1) _col(7) %6.0f el(`matname', `r', 10) ///
            _col(15) %10.4g el(`matname', `r', 11) _col(27) %8.3f el(`matname', `r', 12) ///
            _col(37) %8.3f el(`matname', `r', 13) _col(47) %6.0f el(`matname', `r', 14) ///
            _col(55) %9.3f el(`matname', `r', 15)
    }
    display as text "{hline 68}"

    * footnotes, only the ones that apply
    if "`l2'" != "" {
        display as text "{p 0 6 2}Not converged (m = `l2'): the unrotated fit did not converge. It is the last iterate, not a solution, so the fit statistics and the rotated diagnostics are not reported.{p_end}"
    }
    if "`l3'" != "" {
        display as text "{p 0 6 2}Rotation failed (m = `l3'): the geomin rotation did not converge from any of its starts, and nothing is reported for these rows. dwlsmv itself suggests fewer factors or more starts with rstarts().{p_end}"
    }
    if "`l0'" != "" {
        display as text "{p 0 6 2}df = 0 (m = `l0'): T3, p and the fit indices are shown as dwlsmv posts them; a value dwlsmv does not post is missing.{p_end}"
    }
    if !`epsgiven' & `mmax' >= 3 {
        display as text "{p 0 6 2}The default geomin epsilon changes with m (1e-4 at m = 2, 1e-3 at m = 3, 1e-2 at m >= 4); give epsilon() to use one value throughout.{p_end}"
    }
    if `mmax' >= 2 {
        display as text "{p 0 6 2}Rotated quantities (minima, gap, lam3, phimax) are conditional on the selected minimum.{p_end}"
    }
end
