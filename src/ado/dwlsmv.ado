*! version 0.7.0  dwlsmv -- SEM/CFA for ordinal indicators by DWLS with the WLSMV test
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do;
*! standardized on 2026-10-05: tests/unit/test_standardized.do.
*!
*!   dwlsmv <sem paths> [if] [in] [, metric(marker|unitvariance) keepbound
*!                                   standardized level(#) <sem options>]
*!   dwlsmv (F1 F2 -> y1-y6) [if] [in], efa(F1 F2) [epsilon(#) rstarts(#)
*!                                   keepbound]           (EFA, 0.3.0)
*!   dwlsmv [, standardized level(#) se ci thresholds nothresholds]   (replay)
*!
*! The group surface (spec principle 4; decisions 2026-10-10 and 2026-10-12;
*! build step 2, test_translate_groups.do): group() and invariance() go
*! together; partial(), ginvariant(), a missing or marker metric with
*! invariance(), and fewer than two groups are refused, rc 198; sem reads the
*! groups (ginvariant(none)) and _dwlsmv_translate returns the configural
*! text. Build step 3 (2026-10-12; test_invariance_flow.do): D9 per group,
*! the equal-categories and binary scalar refusals before any fit, Stage 1 per
*! group, dwlsmv_levelgen and the pooled fit. Build step 4 (2026-10-12;
*! test_invariance_post.do): posting under spec D5's multi-group convention
*! (S4-1, S4-2), e(K), e(canon), e(canon_configural), e(invariance), the
*! pooled e(bread_ratio); groupwise Stage 1 notes (S4-3); standardized
*! at G > 1 built 2026-11-01 (O17), at estimation and on replay.
*!
*! sem's estimation options (vce(), method(), ...) are refused, rc 198, and
*! level() is honoured (spec D10, 2026-10-05; test_dwlsmv_gaps.do G11, G12).
*! Help: dwlsmv.sthlp.
*!
*! standardized (D7 step 3, 2026-10-05; test_standardized.do): the table shows
*! g(theta) with delta-method SEs and the derived Theta rows, from e() alone
*! (dwlsmv_m11_std); nothing is stored.
*!
*! The coefficient table (release A, owner 2026-11-18) is the compact grouped table of
*! dwlsmv_m11_ctab.mata, for every fit but efa(): sections, the marker as "(marker)", user
*! constraints as "(constrained)", metric fixings and (at G > 1) level conventions omitted,
*! thresholds last; e() is never held or reposted; r(table) is returned (replay: the last
*! call). Run 2026-11-18 (tests/unit/test_ctab.do, suite).
*!
*! Spec Module 11, D1-D9. sem reads the model (D1, _dwlsmv_translate); dwlsmv
*! fits it. Return codes (D8a): Module 11's own refusals rc 198; a Stage 1
*! abort is caught under capture noisily, followed by the legend and the
*! pointer line, and rethrown as rc 3498; any later Mata abort is rethrown as
*! it is. A fit that does not converge (e(converged) = 0) is posted for inspection
*! without inference, whatever the fit (owner 2026-11-18; a covariate fit since
*! 2026-11-14): e(b), the last iterate, and e(Cns_b); no e(V), e(T3), e(p_T3) or fit
*! indices; the header and the table say so, the replay says the same, and
*! standardized is refused (rc 198).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
program dwlsmv, eclass
    version 16
    if replay() {
        if "`e(cmd)'" != "dwlsmv" error 301
        * nothresholds first: with THResholds declared in the same call, syntax takes nothresholds
        * as its negation and records neither (verified, StataMP 18, 2026-11-18)
        syntax [, NOTHResholds *]
        local 0 , `options'
        syntax [, STANDardized Level(cilevel) SE CI THResholds]
        local isefa = ("`e(efa)'" != "")
        _dwlsmv_dispchk, isefa(`isefa') `se' `ci' `thresholds' `nothresholds'
        if "`standardized'" != "" & "`e(efa)'" != "" {
            display as error "standardized is not available with efa(): the rotated loadings are already standardized."
            exit 198
        }
        if "`standardized'" != "" & "`e(covariates)'" != "" {
            display as error "standardized is not available with a covariate: the standardized solution of a model"
            display as error "conditional on `e(covariates)' has not been derived."
            exit 198
        }
        if "`standardized'" != "" & e(converged) == 0 {
            display as error "standardized is not available for a fit that did not converge: the estimates are"
            display as error "the last iterate, not a solution, and no inference is reported for them."
            exit 198
        }
        _dwlsmv_header
        if "`e(efa)'" != "" _dwlsmv_efatable, level(`level') `se' `ci' `thresholds' `nothresholds'
        else _dwlsmv_table, `standardized' level(`level') `thresholds' `nothresholds'
        exit
    }
    * the rotation-failure mark (owner 2026-11-18, option A; spec Module 12, Surface): a Mata external
    * dwlsmv_rotfail, set only on the geomin "no converged start" exit below, read by efacompare.
    * Dropped on entry to every estimation call (not on replay), so a mark never outlives the call
    * that set it
    capture mata: rmexternal("dwlsmv_rotfail")
    * nothresholds is parsed on its own first (see the replay branch); `0' is kept for e(cmdline)
    local zero `"`0'"'
    syntax anything(name=paths equalok) [if] [in] [, NOTHResholds *]
    local 0 `"`paths' `if' `in' , `options'"'
    syntax anything(name=paths equalok) [if] [in] ///
        [, METRIC(string) KEEPBOUND STANDardized Level(cilevel) ///
        VCE(string) METHod(string) ITERate(string) TECHnique(string) ///
        TOLerance(string) LTOLerance(string) NRTOLerance(string) NONRTOLerance ///
        FROM(string) STARTValues(string) STARTGrouping(string) NOIVSTart ///
        DIFficult GRoup(varname) INVariance(string) PARTIAL(string) ///
        GINVariant(string) EFA(string) EPSilon(string) RSTARTS(string) ///
        SE CI THResholds *]
    local 0 `"`zero'"'

    * display options (release A, owner 2026-11-17): se and ci belong to efa(); thresholds and
    * nothresholds override the default for hiding or showing the threshold rows. rc 198, nothing run.
    local isefa = (`"`efa'"' != "")
    _dwlsmv_dispchk, isefa(`isefa') `se' `ci' `thresholds' `nothresholds'

    * EFA with oblique geomin (0.3.0; spec section 5 Module 12, Surface,
    * 2026-10-29): one complete block, every refusal rc 198 before anything is
    * run. The block's own shape is checked on the canonical text below.
    if `"`efa'"' == "" & regexm(`"`0'"', "[, ]efa\( *\)") {
        display as error "efa() needs the block's factor names, as in efa(F1 F2)."
        exit 198
    }
    if `"`efa'"' == "" & (`"`epsilon'"' != "" | `"`rstarts'"' != "") {
        display as error "epsilon() and rstarts() need efa(); they set the geomin rotation."
        exit 198
    }
    if `"`efa'"' != "" {
        foreach o in group invariance partial ginvariant metric {
            if `"``o''"' != "" {
                display as error "option `o'() is not available with efa() in this version."
                if "`o'" == "metric" {
                    display as error "The EFA block identifies itself: factor variances are fixed at 1."
                }
                display as error "No model was fitted."
                exit 198
            }
        }
        if "`standardized'" != "" {
            display as error "standardized is not available with efa(): the rotated loadings are"
            display as error "already standardized (indicator and factor variances 1). No model was fitted."
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
        local metric unitvariance
    }

    * sem's estimation options would reach sem noestimate on dummy data and
    * change nothing, so a user asking for vce(cluster id) would get the DWLS
    * sandwich without a word: refused, rc 198 (decided 2026-10-05)
    foreach o in vce method iterate technique tolerance ltolerance nrtolerance ///
        nonrtolerance from startvalues startgrouping noivstart difficult {
        if `"``o''"' != "" {
            display as error "option `o' is not available with dwlsmv."
            if "`o'" == "vce" {
                display as error "Standard errors are always the robust (sandwich) DWLS ones;"
                display as error "clustering and sampling weights are planned for a later version."
            }
            else {
                display as error "It is a sem estimation option; dwlsmv estimates by DWLS and"
                display as error "does not use it."
            }
            exit 198
        }
    }
    * the group surface (spec principle 4; decided 2026-10-10 and 2026-10-12),
    * every refusal rc 198 before anything is run, in this order
    if `"`partial'"' != "" {
        display as error "partial invariance is not available in this version."
        exit 198
    }
    if `"`ginvariant'"' != "" {
        display as error "option ginvariant() is not available with dwlsmv; use"
        display as error "invariance(configural|thresholds|metric|scalar). dwlsmv's invariance"
        display as error "levels are not sem's ginvariant() classes; see help dwlsmv."
        exit 198
    }
    if `"`invariance'"' != "" & "`group'" == "" {
        display as error "invariance() needs group(); it compares groups. No model was fitted."
        exit 198
    }
    if "`group'" != "" & `"`invariance'"' == "" {
        display as error "group() without invariance() is not available in this version."
        display as error "Specify invariance(configural), invariance(thresholds), invariance(metric)"
        display as error "or invariance(scalar). No model was fitted."
        exit 198
    }
    if `"`invariance'"' != "" {
        if !inlist(`"`invariance'"', "configural", "thresholds", "metric", "scalar") {
            display as error "invariance() must be configural, thresholds, metric or scalar."
            exit 198
        }
        if "`metric'" != "unitvariance" {
            display as error "invariance() requires metric(unitvariance) in this version; specify metric(unitvariance)."
            display as error "The invariance levels are defined with each factor's variance fixed at 1 in"
            display as error "group 1; the marker metric is not yet available with them. No model was fitted."
            exit 198
        }
    }
    if "`metric'" == "" local metric marker
    if !inlist("`metric'", "marker", "unitvariance") {
        display as error "metric() must be marker or unitvariance"
        exit 198
    }

    marksample touse, novarlist
    local gvals
    if "`group'" != "" {
        * O13 (decided 2026-10-09, amended 2026-10-10): nonnegative integer
        * values; a string variable rc 109; missing values excluded listwise
        capture confirm numeric variable `group'
        if _rc {
            display as error "group variable `group' is a string variable; group() needs nonnegative"
            display as error "integer values. Create one with, e.g., egen newgrp = group(`group')"
            display as error "or encode `group', generate(newgrp)."
            exit 109
        }
        markout `touse' `group'
        capture assert `group' == floor(`group') & `group' >= 0 if `touse'
        if _rc {
            display as error "group variable `group' must take nonnegative integer values; recode it,"
            display as error "e.g. egen newgrp = group(`group')."
            exit 198
        }
        quietly levelsof `group' if `touse', local(gvals)
        if `: word count `gvals'' == 1 {
            display as error "group() has only one value in the estimation sample (`group' = `gvals')."
            display as error "No model was fitted."
            exit 198
        }
        local gopts group(`group') gvals(`gvals')
    }
    _dwlsmv_translate , paths(`paths') semopts(`options') metric(`metric') `gopts'
    local canon    `"`r(canon)'"'
    local ind      `r(ind)'
    local semnames `r(semnames)'
    local cov      `r(cov)'
    local sofac    `r(sofac)'
    if `"`efa'"' != "" & "`sofac'" != "" {
        * release B, B3.5: a second-order factor with efa() is refused before anything is built
        display as error "efa() is not available with a second-order factor (`sofac') in this version. No model was fitted."
        exit 198
    }
    if "`cov'" != "" {
        * a model conditional on a binary covariate (decisions 2026-11-13; owner's Stage 2 package 2026-11-14): standardized,
        * keepbound and efa() are refused (rc 198). keepbound is unavailable with a covariate: a boundary estimate (status 2)
        * is not treated as ordinary without a derivation of its consequences for estimation, admissibility, covariance and inference.
        foreach o in standardized keepbound efa {
            if `"``o''"' != "" {
                display as error "`o' is not available with a covariate (`cov') in this version. No model was fitted."
                exit 198
            }
        }
    }
    if `"`efa'"' != "" {
        mata: st_local("canon", dwlsmv_m11_efa_canon(`"`canon'"', `"`efa'"'))
        if `"`efamsg'"' != "" {
            display as error `"`efamsg'"'
            display as error "No model was fitted."
            exit 198
        }
    }
    * a non-numeric indicator is Module 11's refusal, rc 198 (D8a; decided
    * 2026-10-05), not confirm's own rc
    capture confirm numeric variable `ind'
    if _rc {
        local crc = _rc
        local bad
        foreach v of local ind {
            capture confirm numeric variable `v'
            if _rc local bad `bad' `v'
        }
        display as error "indicators must be numeric; not numeric: `bad'"
        exit 198
    }
    markout `touse' `ind'
    * PROVISIONAL P5 (the owner has not confirmed): a covariate is numeric, takes only 0 and 1 in the
    * estimation sample, and takes both. The same listwise rule as the indicators (one complete-case
    * N). Stage 1 repeats the coding check at rc 3498 (dwlsmv_cx_check).
    if "`cov'" != "" {
        capture confirm numeric variable `cov'
        if _rc {
            display as error "the covariate `cov' must be numeric with the values 0 and 1. No model was fitted."
            exit 198
        }
        markout `touse' `cov'
        quietly count if `touse'
        local ncx = r(N)
        if `ncx' == 0 error 2000
        quietly count if `touse' & `cov' == 0
        local ncx0 = r(N)
        quietly count if `touse' & `cov' == 1
        local ncx1 = r(N)
        if `ncx0' + `ncx1' != `ncx' {
            display as error "the covariate `cov' must take only the values 0 and 1 (a binary covariate coded 0/1)."
            display as error "Recode it; the coding fixes the unit of its effect and the reported thresholds. No model was fitted."
            exit 198
        }
        if `ncx0' == 0 | `ncx1' == 0 {
            display as error "the covariate `cov' must take both the values 0 and 1 in the estimation sample; it is constant."
            display as error "A constant covariate has no estimable effect. No model was fitted."
            exit 198
        }
    }
    quietly count if `touse'
    if r(N) == 0 error 2000
    if "`group'" != "" {
        * the estimation sample is final only now: a group whose every
        * observation has a missing indicator is gone (O13 maps the values
        * in the estimation sample), so the values are read again
        quietly levelsof `group' if `touse', local(gvals2)
        if `: word count `gvals2'' == 1 {
            display as error "group() has only one value in the estimation sample (`group' = `gvals2')."
            display as error "No model was fitted."
            exit 198
        }
        if "`gvals2'" != "`gvals'" {
            local gvals `gvals2'
            _dwlsmv_translate , paths(`paths') semopts(`options') metric(`metric') ///
                group(`group') gvals(`gvals')
            local canon    `"`r(canon)'"'
            local semnames `r(semnames)'
        }
        * build step 3 (2026-10-12): D9 per group, the table checks before
        * any fit (F4 order), Stage 1 per group, the level and its fit
        local kb = ("`keepbound'" != "")
        mata: dwlsmv_m11_gsupport("`touse'", "`ind'", "`group'", "`gvals'")
        capture noisily mata: dwlsmv_m11_gcheck(`"`canon'"', "`invariance'", ///
            "`ind'", "`group'", "`gvals'")
        local rc = _rc
        if `rc' | "`gref'" == "1" {
            capture mata: dwlsmv_m11_clear()
            if `rc' exit `rc'
            exit 198
        }
        capture noisily mata: dwlsmv_m11_stage1_g(`kb')
        local rc = _rc
        if `rc' {
            mata: dwlsmv_m11_legend("`ind'")
            display as error "dwlsmv's default differs from lavaan's and Mplus's here; see"
            display as error "help dwlsmv, section Comparison with lavaan and Mplus."
            if !`kb' {
                capture quietly mata: dwlsmv_m11_stage1_g(1)
                if !_rc {
                    display as error "keep_bound = 1 above is the dwlsmv option keepbound."
                }
            }
            capture mata: dwlsmv_m11_clear()
            exit `rc'
        }
        mata: dwlsmv_m11_legend_note("`ind'", 1)
        capture noisily mata: dwlsmv_m11_fit_g(`"`canon'"', "`invariance'")
        local rc = _rc
        if `rc' {
            capture mata: dwlsmv_m11_clear()
            exit `rc'
        }
        * build step 4 (2026-10-12): posting (spec D5's multi-group
        * convention S4-1, S4-2; what each fit posts for the exemption; S4-3)
        tempname b V C K
        capture noisily mata: dwlsmv_m11_post_g(`"`canon'"', "`ind'", "`group'", ///
            "`gvals'", "`semnames'", "`b'", "`V'", "`C'", "`K'")
        local rc = _rc
        if `rc' {
            capture mata: dwlsmv_m11_clear()
            exit `rc'
        }
        * non-convergence is an inference boundary (owner 2026-11-18): the last iterate and the fixed
        * rows only (e(Cns_b); e(Cns) cannot be set without e(V)); no e(V), T3, p or fit indices
        if `converged' == 0 {
            ereturn post `b', esample(`touse') obs(`N')
            if `nfix' ereturn matrix Cns_b = `C'
        }
        else if `nfix' ereturn post `b' `V' `C', esample(`touse') obs(`N')
        else           ereturn post `b' `V', esample(`touse') obs(`N')
        ereturn matrix K = `K'
        ereturn scalar rank        = `rank'
        ereturn scalar df_T3       = `df'
        if `converged' != 0 {
            ereturn scalar T3          = `T3'
            ereturn scalar p_T3        = cond(`df' > 0, chi2tail(`df', `T3'), .)
            ereturn scalar rmsea_s     = `rmsea_s'
            ereturn scalar cfi_s       = `cfi_s'
            ereturn scalar tli_s       = `tli_s'
            ereturn scalar T3_b        = `T3_b'
            ereturn scalar df_T3_b     = `df_T3_b'
        }
        ereturn scalar converged   = `converged'
        ereturn scalar iterations  = `iterations'
        ereturn scalar bread_ratio = `bratio'
        ereturn local invariance "`invariance'"
        ereturn local metric     "`metric'"
        ereturn local recovery   "none"
        ereturn local indicators "`ind'"
        ereturn local cmdline    `"dwlsmv `0'"'
        ereturn local version    "0.7.0"
        // what difftest refits from at G > 1 (spec Module 12, 2026-10-31):
        // e(canon) and e(canon_configural) (dwlsmv_m11_post_g_e), the group
        // variable, keepbound, and the data signature over the indicators
        // and the group variable
        ereturn local groupvar   "`group'"
        ereturn local keepbound  "`keepbound'"
        signestimationsample `ind' `group'
        ereturn local cmd        "dwlsmv"
        mata: dwlsmv_m11_post_g_e()

        _dwlsmv_header
        * a standardized refusal (O17) leaves the fit posted; notes and cleanup still run
        capture noisily _dwlsmv_table, `standardized' level(`level') `thresholds' `nothresholds'
        local tabrc = _rc
        mata: dwlsmv_m11_gnotes(`"`canon'"', "`ind'","`group'", "`gvals'")
        capture mata: dwlsmv_m11_clear()
        if `tabrc' exit `tabrc'
        exit
    }

    local kb = ("`keepbound'" != "")
    capture noisily mata: dwlsmv_m11_stage1("`touse'", "`ind'", `kb', "`cov'")
    local rc = _rc
    if `rc' {
        mata: dwlsmv_m11_legend("`ind'")
        if "`cov'" == "" {
            display as error "dwlsmv's default differs from lavaan's and Mplus's here; see"
            display as error "help dwlsmv, section Comparison with lavaan and Mplus."
        }
        * only when keepbound would let the data through: Stage 1 rerun with
        * keep_bound = 1, quietly (decided 2026-10-04); not with a covariate
        if !`kb' & "`cov'" == "" {
            capture quietly mata: dwlsmv_m11_stage1("`touse'", "`ind'", 1)
            if !_rc {
                display as error "keep_bound = 1 above is the dwlsmv option keepbound."
            }
        }
        capture mata: dwlsmv_m11_clear()
        exit `rc'
    }

    if "`cov'" == "" {
        mata: dwlsmv_m11_legend_note("`ind'", 0)
    }

    tempname b V C bnd
    capture noisily mata: dwlsmv_m11_fit(`"`canon'"', "`ind'", "`semnames'", ///
        "`b'", "`V'", "`C'", "`bnd'", "`metric'")
    local rc = _rc
    if `rc' {
        * Module 3's message stays as it is (D8a); the remedy in the
        * command's terms for each factor left without a scale
        mata: st_local("unanc", dwlsmv_m11_unanchored(`"`canon'"'))
        while "`unanc'" != "" {
            gettoken one unanc : unanc, parse(";")
            if "`one'" == ";" continue
            gettoken fac y1 : one
            local y1 = strtrim("`y1'")
            if "`y1'" == "" {
                * a factor with no loading line (second-order): its scale is its variance
                display as error "Factor `fac' has no fixed scale: fix its variance with var(`fac'@1)."
                continue
            }
            display as error "Factor `fac' has no fixed scale: fix one of its loadings at 1,"
            display as error "e.g. (`fac' -> `y1'@1 ...), or specify metric(unitvariance)."
        }
        capture mata: dwlsmv_m11_clear()
        exit `rc'
    }

    if `"`efa'"' != "" {
        * rotated estimates; their covariance e(V) from release A (unit A5) unless withheld (A3.2, A3.4)
        tempname rb lam phi
        capture mata: dwlsmv_m11_efa_rotate("`b'", `"`canon'"', ///
            "`epsilon'", "`rstarts'", "`rb'", "`lam'", "`phi'", `"`efa'"')
        local rc = _rc
        if `rc' & "`rot_stage'" != "rotate" {
            * not the rotation's non-convergence: show the real abort
            * (CLAUDE.md: capture hides it), rethrown as it is
            capture noisily mata: dwlsmv_m11_efa_rotate("`b'", `"`canon'"', ///
                "`epsilon'", "`rstarts'", "`rb'", "`lam'", "`phi'", `"`efa'"')
            capture mata: dwlsmv_m11_clear()
            exit `rc'
        }
        if `rc' {
            * no converged start is an error (spec Module 12), explained in
            * the command's terms (the user's choice, 2026-10-29); nothing posted
            local ns = cond("`rstarts'" == "", 30, real("`rstarts'"))
            display as error "The geomin rotation did not converge from any of its `ns' starts."
            display as error "The unrotated fit " ///
                cond(`converged' == 1, "converged", "did not converge") ///
                cond("`rot_inadm'" == "1", " and is inadmissible (an indicator's communality exceeds 1)", "") "."
            display as error "This often means more factors than the data support: try fewer factors,"
            display as error "or more starts with rstarts(). No results were posted."
            capture mata: dwlsmv_m11_clear()
            * option A (owner 2026-11-18): the mark efacompare reads, set after the clear (the clear
            * removes only the dwlsmv_m11_* externals, so the order is not load-bearing), only here
            capture mata: dwlsmv_rotfail = 1
            exit 3498
        }
        * release A (unit A5, spec Module 12 A3.2-A3.5, A3.13-A3.15): the covariance of the rotated
        * estimates, J ACOV_u J' at the GPA stall point, over every rotated e(b) column. No e(V) (the
        * rotated e(b) alone, A3.4) for an unrotated fit that did not converge or is inadmissible, for
        * the section 4.8 rank cut (A3.2), and for a degenerate critical point or a missing value.
        * e(rot_nose) names the reason; it is absent when e(V) is posted.
        tempname Vr
        local nose ""
        local rot_cond .
        if `converged' == 0 local nose unconverged
        else if "`rot_inadm'" == "1" local nose inadmissible
        if "`nose'" == "" {
            capture noisily mata: st_local("sest", strofreal(dwlsmv_m11_efa_se("`b'", "`V'", ///
                `"`canon'"', "`rb'", "`lam'", "`phi'", "`Vr'", `"`efa'"')))
            local rc = _rc
            if `rc' {
                capture mata: dwlsmv_m11_clear()
                exit `rc'
            }
            if `sest' == 1 local nose rank
            else if `sest' == 2 local nose critical
            else if `sest' == 3 local nose acov
        }
        * e(Ev) (release A, owner 2026-11-17): the eigenvalues of the sample polychoric/tetrachoric
        * correlation matrix, descending, 1 x p; read from Stage 1's s before dwlsmv_m11_clear(); properties(... eigen) is what screeplot checks (has_eprop)
        tempname ev
        capture noisily mata: dwlsmv_m11_efa_ev("`ev'")
        if _rc {
            local rc = _rc
            capture mata: dwlsmv_m11_clear()
            exit `rc'
        }
        if "`nose'" == "" ereturn post `rb' `Vr', esample(`touse') obs(`N') properties(b V eigen)
        else               ereturn post `rb', esample(`touse') obs(`N') properties(b eigen)
        ereturn matrix Lambda = `lam'
        ereturn matrix Phi    = `phi'
        ereturn matrix Ev     = `ev'
        foreach s in starts converged minima atbest Q gpnorm eps gap {
            ereturn scalar rot_`s' = `rot_`s''
        }
        if "`nose'" != "" ereturn local rot_nose "`nose'"
        if `rot_cond' < . ereturn scalar rot_cond = `rot_cond'
        ereturn local efa      "`efa'"
        ereturn local rotation "geomin oblique"
    }
    else {
        if `converged' == 0 {
            // non-convergence is an inference boundary (owner 2026-11-18; covariate fits since 2026-11-14):
            // the last iterate as point estimates only. No e(V) (the sandwich at a non-solution is not
            // inference), as the EFA block posts e(b) only; e(Cns_b) (e(Cns) cannot be set without e(V))
            // is kept so that the table can mark the fixed parameters
            ereturn post `b', esample(`touse') obs(`N')
            if `nfix' ereturn matrix Cns_b = `C'   // e(Cns) is reserved: ereturn matrix refuses the name (rc 198, probed 2026-11-14)
        }
        else {
            if `nfix' ereturn post `b' `V' `C', esample(`touse') obs(`N')
            else      ereturn post `b' `V', esample(`touse') obs(`N')
        }
    }
    if "`cov'" == "" ereturn matrix bound = `bnd'
    ereturn scalar rank       = `rank'
    ereturn scalar df_T3      = `df'
    // T3, p and the fit indices only for a converged fit (owner 2026-11-18), CFA or efa()
    if `converged' != 0 {
        ereturn scalar T3     = `T3'
        ereturn scalar p_T3   = cond(`df' > 0, chi2tail(`df', `T3'), .)
    }
    // fit indices, option A (family S, frozen 2026-10-21); printed by _dwlsmv_header.
    // Omitted, not posted as missing, for a conditional fit (decision 2026-11-13 (a)) and for a
    // fit that did not converge, with the baseline's T3 and df (2026-11-18): the header prints
    // them only when posted
    if "`cov'" == "" & `converged' != 0 {
        ereturn scalar rmsea_s    = `rmsea_s'
        ereturn scalar cfi_s      = `cfi_s'
        ereturn scalar tli_s      = `tli_s'
        ereturn scalar T3_b       = `T3_b'
        ereturn scalar df_T3_b    = `df_T3_b'
        // residual diagnostics, under the same rule (@sec-resid-diag SR5, 2026-11-22): SRMR and the
        // signed and normalized residual matrices; set by dwlsmv_m11_fit for a converged fit
        ereturn scalar srmr       = `srmr'
        ereturn matrix Resid      = `resid_m'
        ereturn matrix Nresid     = `nresid_m'
    }
    ereturn scalar converged  = `converged'
    ereturn scalar iterations = `iterations'
    // Module 10's ratios, stored and not printed (2026-10-07); no cutoff. A conditional fit posts all four,
    // as an ordinary one does (2026-11-15: Gamma's three, and the bread's once measured on production's
    // conditional Delta, "Yes, both"); posted whether or not the fit converged
    ereturn scalar gamma_ratio     = `gratio'
    ereturn scalar gamma_ratio_raw = `gratior'
    ereturn scalar gamma_rank      = `grank'
    ereturn scalar bread_ratio     = `bratio'
    ereturn local metric   "`metric'"
    ereturn local recovery "`recovery'"
    ereturn local indicators "`ind'"
    if "`cov'" != "" ereturn local covariates "`cov'"
    ereturn local cmdline  `"dwlsmv `0'"'
    // the release; must equal line 1 (test_dwlsmv_gaps.do G14, stage_release.ps1)
    ereturn local version  "0.7.0"
    // what difftest refits from (spec Module 12, D8d addition, 2026-10-23):
    // the canonical text as dwlsmv_m11_fit received it, keepbound, and the
    // data signature that check 2 tests with checkestimationsample
    ereturn local canon     `"`canon'"'
    ereturn local keepbound "`keepbound'"
    signestimationsample `ind' `cov'
    ereturn local cmd      "dwlsmv"

    _dwlsmv_header
    local tabrc 0
    if `"`efa'"' != "" _dwlsmv_efatable, level(`level') `se' `ci' `thresholds' `nothresholds'
    else {
        * standardized on a fit that did not converge is refused (rc 198); the fit stays posted and
        * the notes and cleanup still run, as at G > 1 (O17)
        capture noisily _dwlsmv_table, `standardized' level(`level') `thresholds' `nothresholds'
        local tabrc = _rc
    }
    mata: dwlsmv_m11_notes()
    capture mata: dwlsmv_m11_clear()
    if `tabrc' exit `tabrc'
end

* EFA (0.3.0, Unit 6; standard errors from release A, unit A5): rotated loadings and factor
* correlations; from e() alone, so replay works. Release A display (owner 2026-11-17, revising
* A3.6's SE/z/p/CI columns): by default the two matrices only, no standard errors, intervals or
* threshold rows, and one line saying e(V) is posted. Option se adds the standard-error matrices;
* option ci puts a bracketed Wald interval at level() beside each estimate (a grouped long table
* when m is too wide for the line); either prints the caveat of A3.1. thresholds / nothresholds:
* the threshold rows of e(b), hidden by default. Without e(V) (e(rot_nose), A3.2, A3.4): the two
* matrices and no standard errors; the header carries the reason.
* The wording of the caveat, the minima and tie notes and the no-SE reasons is A7's (owner 2026-11-18).
program _dwlsmv_efatable
    syntax [, NOTHResholds *]
    local 0 , `options'
    syntax [, Level(cilevel) SE CI THResholds]
    _dwlsmv_thrshow, `thresholds' `nothresholds'
    local showthr = r(show)
    tempname L P SL SP b V
    matrix `L' = e(Lambda)
    matrix `P' = e(Phi)
    display _n as text "Rotation: geomin, oblique (epsilon = " as result e(rot_eps) ///
        as text ", " as result e(rot_starts) as text " starts, " ///
        as result e(rot_converged) as text " converged)"
    capture confirm matrix e(V)
    * no e(V) to use when the fit did not converge (owner 2026-11-18), also when e(converged) is set by hand
    local hasv = (_rc == 0 & "`e(rot_nose)'" == "" & e(converged) != 0)
    local wantse = ("`se'" != "")
    local wantci = ("`ci'" != "")
    if `hasv' & (`wantse' | `wantci') {
        * plain copies, so the Mata routine reads names it can rely on
        matrix `b' = e(b)
        matrix `V' = e(V)
        mata: dwlsmv_m11_efa_semat("`b'", "`V'", "`L'", "`SL'", "`SP'")
    }
    if `hasv' & `wantci' {
        _dwlsmv_efaci `L', se(`SL') level(`level') title("Rotated loadings")
        _dwlsmv_efaci `P', se(`SP') level(`level') title("Factor correlations") phi
    }
    else {
        matlist `L', title("Rotated loadings") format(%9.3f) border(rows) lines(oneline)
        matlist `P', title("Factor correlations") format(%9.3f) border(rows) lines(oneline)
    }
    if `hasv' & `wantse' {
        matlist `SL', title("Standard errors of the rotated loadings") format(%9.3f) border(rows) lines(oneline)
        matlist `SP', title("Standard errors of the factor correlations") format(%9.3f) border(rows) lines(oneline)
    }
    if `hasv' & (`wantse' | `wantci') {
        * A3.1 and A3.6: one caveat for every rotated estimate (A7 wording, owner 2026-11-18, variant S)
        display as text "Standard errors of the rotated estimates are a local first-order approximation,"
        display as text "conditional on the selected rotation. Their accuracy depends on how far a"
        display as text "loading is from zero relative to sqrt(epsilon), and on N. In simulation,"
        display as text "intervals for some loadings near zero covered well below the nominal level."
        display as text "See {help dwlsmv##efase:help dwlsmv}, Remarks on rotated standard errors."
    }
    else if `hasv' {
        display as text "The covariance of the rotated estimates is posted in e(V), so lincom and test can use it;"
        display as text "options se and ci display standard errors and Wald intervals."
        * the caveat also where e(V) is offered without se or ci (owner 2026-11-18)
        display as text "They are approximate, least accurate near zero; see {help dwlsmv##efase:help dwlsmv}."
    }
    else {
        display as text "No standard errors are reported for the rotated estimates (see above)."
        if `wantse' | `wantci' {
            display as text "Options se and ci need e(V), which is not posted for this fit."
        }
    }
    if `showthr' _dwlsmv_efathr, level(`level')
    * no threshold note on efa() output (owner veto, 2026-11-18): thresholds alone adds the table
    if e(rot_minima) > 1 {
        display as text "The rotation found " as result e(rot_minima) as text " distinct minima; " ///
            as result e(rot_atbest) as text " of the converged starts reached the reported one" ///
            " (criterion gap to the next " as result %9.2e e(rot_gap) as text ")."
        if `hasv' {
            * A3.3 (A7 wording, owner 2026-11-18)
            display as text "Standard errors are conditional on the reported minimum and differ between"
            display as text "minima."
        }
    }
    * A3.7 revisited (owner 2026-11-18): name a tie on the column-ordering statistic and the key
    * that decided it (spec Module 12). Computed from e(Lambda) and e(Phi). A7 wording (owner
    * 2026-11-18); the pointer to the help once, after the last tying pair.
    tempname tk
    local nt = colsof(e(Lambda)) - 1
    if `nt' > 0 {
        mata: st_matrix("`tk'", dwlsmv_geomin_tiekeys(st_matrix("e(Lambda)"), st_matrix("e(Phi)")))
    }
    local fn : colnames e(Lambda)
    local anytie 0
    forvalues j = 1/`nt' {
        local k = `tk'[1, `j']
        if `k' != 1 {
            local anytie 1
            local a : word `j' of `fn'
            local b : word `=`j'+1' of `fn'
            display as text "Factors `a' and `b' tie on the statistic that orders the factors."
            if `k' == 2 display as text "They are ordered by the sum of squared loadings, larger first."
            else if `k' == 3 display as text "They are ordered by their loadings item by item, larger first."
            else if `k' == 4 display as text "They are ordered by their factor correlations; their loadings are identical."
            else display as text "Their loadings and correlations are identical, so the order changes no estimate."
        }
    }
    if `anytie' {
        display as text "See {help dwlsmv##efaorder:help dwlsmv}, Remarks on the order of the factors."
    }
    if e(df_T3) == 0 {
        display as text "With df = 0 the model fits exactly; T3 and the fit indices are not informative."
    }
end

* Display options (release A, owner 2026-11-17). se and ci belong to efa(): an ordinary fit's
* coefficient table already carries standard errors and intervals. thresholds and nothresholds
* contradict each other. Every refusal is rc 198, as Module 11's (D8a).
program _dwlsmv_dispchk
    syntax [, NOTHResholds *]
    local 0 , `options'
    syntax , ISEFA(integer) [SE CI THResholds]
    if "`thresholds'" != "" & "`nothresholds'" != "" {
        display as error "thresholds and nothresholds contradict each other; specify one."
        exit 198
    }
    if `isefa' == 0 & ("`se'" != "" | "`ci'" != "") {
        display as error "options se and ci belong to efa(): the coefficient table already shows"
        display as error "standard errors and intervals for an ordinary fit."
        exit 198
    }
end

* Whether the threshold rows are shown (r(show)): thresholds shows, nothresholds hides; with neither,
* shown only when invariance() constrains them across groups (thresholds, metric, scalar), hidden
* otherwise. An EFA fit never has e(invariance), so its default is hidden.
program _dwlsmv_thrshow, rclass
    syntax [, NOTHResholds *]
    local 0 , `options'
    syntax [, THResholds]
    if "`thresholds'" != "" & "`nothresholds'" != "" {
        display as error "thresholds and nothresholds contradict each other; specify one."
        exit 198
    }
    local show 0
    if "`thresholds'" != "" local show 1
    else if "`nothresholds'" == "" & inlist("`e(invariance)'", "thresholds", "metric", "scalar") local show 1
    return scalar show = `show'
end

program _dwlsmv_thrnote
    display as text "Thresholds are not shown; the option thresholds shows them (they are in e(b))."
end

* The thresholds of an EFA fit, from e(b) and e(V), as a coefficient table (or point estimates only
* without e(V)). Its own program, so that _estimates restore runs before the caller shows anything else.
program _dwlsmv_efathr
    syntax [, Level(cilevel)]
    tempname b V tb tV hold
    matrix `b' = e(b)
    capture confirm matrix e(V)
    local hasv = (_rc == 0 & "`e(rot_nose)'" == "" & e(converged) != 0)
    if `hasv' {
        matrix `V' = e(V)
        mata: st_local("nk", strofreal(dwlsmv_m11_dispsub("`b'", "`V'", "`tb'", "`tV'")))
        if `nk' == 0 exit
        _estimates hold `hold', copy restore
        ereturn post `tb' `tV'
        display _n as text "Thresholds"
        ereturn display, level(`level')
    }
    else {
        mata: st_local("nk", strofreal(dwlsmv_m11_dispsub("`b'", "", "`tb'", "")))
        if `nk' == 0 exit
        matrix `tb' = `tb''
        matlist `tb', title("Thresholds (no standard errors)") format(%9.3f) border(rows) lines(oneline)
    }
end

* One estimate matrix with a bracketed Wald interval beside each entry: est [lo, hi], the interval
* est -/+ z se at level(), z the normal quantile, printed as computed (a correlation's interval
* is not clipped to [-1, 1]; a note says so when it crosses). Wide form when the m columns fit
* c(linesize); otherwise a long table grouped by factor. phi: the matrix is the factor correlations,
* so only the strict lower triangle carries an interval (the diagonal is fixed at 1).
program _dwlsmv_efaci
    syntax namelist(max=1), SE(name) Level(cilevel) TItle(string) [PHI]
    tempname E S
    matrix `E' = `namelist'
    matrix `S' = `se'
    local p = rowsof(`E')
    local m = colsof(`E')
    local rn : rownames `E'
    local cn : colnames `E'
    local z = invnormal(.5 + `level'/200)
    local rw 12
    local cw 28
    local wide = (`rw' + 2 + `m' * `cw' <= c(linesize))
    local clip 0
    display _n as text "`title'" " (estimate and " `level' "% Wald interval)"
    if `wide' {
        display as text %-`rw's "" " {c |}" _continue
        forvalues k = 1/`m' {
            local f : word `k' of `cn'
            local hd = "   " + abbrev("`f'", `=`cw' - 4')
            display as text %-`cw's "`hd'" _continue
        }
        display
        display as text "{hline `=`rw' + 1'}{c +}{hline `=`m' * `cw''}"
        forvalues j = 1/`p' {
            local r : word `j' of `rn'
            display as text %-`rw's abbrev("`r'", `rw') " {c |}" _continue
            forvalues k = 1/`m' {
                local est = `E'[`j', `k']
                local sev = `S'[`j', `k']
                if "`phi'" != "" & `j' < `k' {
                    display as text %-`cw's "" _continue
                }
                else if "`phi'" != "" & `j' == `k' {
                    * the estimate right-aligned in at least 6 characters (owner 2026-11-18), so signs line up; never cut
                    local cell = max(0, 6 - strlen(strofreal(`est', "%6.3f"))) * " " + strofreal(`est', "%6.3f")
                    display as result %-`cw's "  `cell'" _continue
                }
                else {
                    local lo = `est' - `z' * `sev'
                    local hi = `est' + `z' * `sev'
                    if "`phi'" != "" & `sev' < . {
                        if `lo' < -1 | `hi' > 1 local clip 1
                    }
                    local cell = max(0, 6 - strlen(strofreal(`est', "%6.3f"))) * " " + strofreal(`est', "%6.3f") + " [" + strofreal(`lo', "%7.3f") + ", " + strofreal(`hi', "%7.3f") + "]"
                    display as result %-`cw's "  `cell'" _continue
                }
            }
            display
        }
    }
    else {
        forvalues k = 1/`m' {
            local f : word `k' of `cn'
            local j0 = cond("`phi'" != "", `k' + 1, 1)
            if `j0' > `p' continue
            if "`phi'" != "" display _n as text "Correlations with `f'"
            else display _n as text "Loadings on `f'"
            display as text %-`rw's "" " {c |}" %9s "Estimate" "   [" `level' "% Wald interval]"
            display as text "{hline `=`rw' + 1'}{c +}{hline 38}"
            forvalues j = `j0'/`p' {
                local r : word `j' of `rn'
                local est = `E'[`j', `k']
                local sev = `S'[`j', `k']
                local lo = `est' - `z' * `sev'
                local hi = `est' + `z' * `sev'
                if "`phi'" != "" & `sev' < . {
                    if `lo' < -1 | `hi' > 1 local clip 1
                }
                local iv = "[" + strofreal(`lo', "%8.3f") + ", " + strofreal(`hi', "%8.3f") + "]"
                display as text %-`rw's abbrev("`r'", `rw') " {c |} " as result %9.3f `est' ///
                    "   `iv'"
            }
        }
    }
    if `clip' {
        display as text "The Wald interval of a factor correlation is printed as computed; at least one crosses -1 or 1."
        display as text "A transformed interval is not offered in this version."
    }
end

* A model with a covariate that did not converge: the last iterate as a plain table, from e(b) and
* e(Cns_b) alone (replay works). No standard errors, z, p or interval exist to show: e(V) is not posted.
*
* Release A, owner 2026-11-18 (compact table; dwlsmv_m11_ctab.mata): an Estimate column only, the marker
* as "(marker)" and any other fixed parameter as "(constrained)", sections as in the coefficient table.
* Optional level() and rtab(name): the level of r(table)'s crit row, and the Stata matrix that
* receives the full r(table) (one column per e(b) column, hidden rows included).
program _dwlsmv_nctable
    syntax , SHOW(integer) [Level(cilevel) RTAB(name)]
    tempname b C T
    matrix `b' = e(b)
    * a hand-set e(converged) = 0 on a fit posted with e(V) carries e(Cns), not e(Cns_b). Looked up in
    * e(matrices): "matrix C = e(Cns_b)" does not fail when there is none, it gives a 1 x 1 missing
    * (rc 0; found 2026-11-18 by test_ctab.do 5g)
    local emats : e(matrices)
    local cn ""
    if `: list posof "Cns_b" in emats' {
        matrix `C' = e(Cns_b)
        local cn "`C'"
    }
    else if `: list posof "Cns" in emats' {
        matrix `C' = e(Cns)
        local cn "`C'"
    }
    display as text "Last iterate, not a solution; no standard errors"
    * the threshold rows are hidden unless thresholds was given (release A, owner 2026-11-17)
    mata: dwlsmv_m11_ctab("`b'", "", "`cn'", `show', 0, `level', `=c(linesize)', "`T'")
    if "`rtab'" != "" matrix `rtab' = `T'
    _dwlsmv_ctnote, nlev(`ct_nlev')
    if !`show' _dwlsmv_thrnote
end

* The line under the table at G > 1 when level conventions were left out (C1 (a), owner 2026-11-18)
program _dwlsmv_ctnote
    syntax , NLEV(integer)
    if `nlev' > 0 {
        display as text "Fixed by the level, not shown: variances 1, intercepts 0, scales 1, means 0."
    }
end

* D7 step (3): standardized replaces the coefficient table; Theta rows appear
* only here; nothing is stored (D8d). Computed from e() alone, so replay
* works. Release A (owner 2026-11-18): one compact table, printed by dwlsmv_m11_ctab; e() is never
* held, reposted or restored, so lincom and test are untouched. r(table) is returned, one column
* per e(b) column (or per standardized row), rows b se z pvalue ll ul df crit eform.
program _dwlsmv_table, rclass
    syntax [, NOTHResholds *]
    local 0 , `options'
    syntax [, STANDardized Level(cilevel) THResholds]
    * the threshold rows (release A, owner 2026-11-17): hidden by default, shown by default only
    * when invariance() constrains them, thresholds / nothresholds override either default.
    * Display only: e(b) and its order are unchanged. A one-line note says when they are hidden.
    _dwlsmv_thrshow, `thresholds' `nothresholds'
    local show = r(show)
    * a fit that did not converge (owner 2026-11-18; a covariate fit since 2026-11-14): point
    * estimates only, no e(V) to display; standardized needs e(V) for its standard errors
    if e(converged) == 0 {
        if "`standardized'" != "" {
            display as error "standardized is not available for a fit that did not converge: the estimates are"
            display as error "the last iterate, not a solution, and no inference is reported for them."
            exit 198
        }
        tempname T0
        _dwlsmv_nctable, show(`show') level(`level') rtab(`T0')
        return matrix table = `T0'
        exit
    }
    tempname bs Vs Cs T
    local inv "`e(invariance)'"
    local nfixs 0
    if "`standardized'" == "" {
        matrix `bs' = e(b)
        matrix `Vs' = e(V)
        * e(matrices), not capture: "matrix Cs = e(Cns)" gives a 1 x 1 missing, rc 0, when there is none
        local emats : e(matrices)
        local nfixs = (`: list posof "Cns" in emats' > 0)
        if `nfixs' matrix `Cs' = e(Cns)
    }
    else {
        * at G > 1 the group routine fails closed before anything is shown (O17)
        if "`inv'" != "" mata: dwlsmv_m11_std_g("`bs'", "`Vs'", "`Cs'")
        else             mata: dwlsmv_m11_std("`bs'", "`Vs'", "`Cs'")
        display as text "Standardized (std.all); /var(e.y) rows are the derived residual variances"
        * C2 (a), owner 2026-11-18: the constant rows are not shown, at every G, and the title says so
        display as text "Constant rows (fixed by the level or by the standardization) are not shown."
    }
    local cn = cond(`nfixs', "`Cs'", "")
    mata: dwlsmv_m11_ctab("`bs'", "`Vs'", "`cn'", `show', ("`standardized'" != ""), `level', `=c(linesize)', "`T'")
    * the level conventions of a group fit are said once, under the table (C1 (a)); the standardized
    * title has said it
    if "`standardized'" == "" _dwlsmv_ctnote, nlev(`ct_nlev')
    if !`show' _dwlsmv_thrnote
    return matrix table = `T'
    if "`standardized'" == "" exit
    * Q-C (a), user 2026-11-01: one note where parameters are tied
    if !inlist("`inv'", "", "configural") {
        display as text "Note: equal raw parameters can have different standardized values across groups."
    }
end

program _dwlsmv_header
    display _n as text "DWLS estimation, WLSMV test (ordinal indicators)" ///
        _col(51) "Number of obs" _col(67) "= " as result %9.0f e(N)
    * a fit that did not converge (owner 2026-11-18, every fit; a covariate fit since 2026-11-14): said
    * first, in error colour
    if e(converged) == 0 {
        display as error "convergence not achieved: the optimizer stopped after " e(iterations) " iterations."
        display as error "The estimates below are the last iterate, not a solution."
        display as error "No standard errors, T3 or p-value are reported for them."
        if "`e(metric)'" == "marker" {
            * recovery (D11) now runs for a conditional marker fit, so a fit still not converged here is
            * one whose recovery was ineligible or failed (a fit that converged after it is not in this block)
            if "`e(recovery)'" == "ineligible" {
                display as error "Refit with metric(unitvariance). No recovery refit was attempted:"
                display as error "the model is not eligible for it (see the notes)."
            }
            else if "`e(recovery)'" == "failed" {
                display as error "Refit with metric(unitvariance). A recovery refit was attempted"
                display as error "and failed (see the notes)."
            }
            else {
                display as error "Refit with metric(unitvariance)."
            }
        }
        else if "`e(covariates)'" != "" {
            display as error "No further remedy is built for a model with a covariate;"
            display as error "dwlsmv has no option for other starting values."
        }
        else {
            display as error "dwlsmv has no option for other starting values."
        }
    }
    if "`e(invariance)'" != "" {
        * build step 4 (2026-10-12); wording (a), approved by the user the same day
        tempname k
        matrix `k' = e(K)
        display as text "Invariance: `e(invariance)', " rowsof(`k') " groups;" ///
            " factor variances fixed at 1 in the first group (the reference group)"
    }
    else if "`e(metric)'" == "marker" {
        display as text "Metric: marker (first-listed indicator)"
    }
    else {
        * H-1 (a), decided 2026-10-12 (user): plain terms; psi_ff is the
        * disturbance variance for an endogenous factor, never its total
        display as text "Metric: factor variances fixed at 1 (for a factor predicted by other factors, its residual variance)"
    }
    if "`e(covariates)'" != "" {
        display as text "Conditional on covariate: " as result "`e(covariates)'" as text ///
            " (binary 0/1); thresholds are those at " "`e(covariates)'" " = 0, and (y <- " "`e(covariates)'" ") is minus a threshold shift"
    }
    if e(converged) == 0 {
        display as text "Scaled-and-shifted test (T3): not reported (not a solution); df = " as result e(df_T3)
    }
    else {
        display as text "Scaled-and-shifted test (T3):" _col(51) "chi2(" as result e(df_T3) ///
            as text ")" _col(67) "= " as result %9.3f e(T3)
        display as text _col(51) "Prob > chi2" _col(67) "= " as result %9.4f e(p_T3)
    }
    * fit indices, option A (family S; F1-F6 frozen 2026-10-21); shown only
    * when posted, so a replay of results without them prints nothing; a
    * missing value (F3, F4) prints as "."
    if "`e(rmsea_s)'" != "" & e(converged) != 0 {
        display as text _col(51) "RMSEA (scaled)" _col(67) "= " as result %9.4f e(rmsea_s)
        display as text _col(51) "CFI (scaled)"   _col(67) "= " as result %9.4f e(cfi_s)
        display as text _col(51) "TLI (scaled)"   _col(67) "= " as result %9.4f e(tli_s)
    }
    * SRMR (@sec-resid-diag, owner decision 2 of 2026-11-21): under the indices, only when posted
    if "`e(srmr)'" != "" & e(converged) != 0 {
        display as text _col(51) "SRMR"           _col(67) "= " as result %9.4f e(srmr)
    }
    if "`e(covariates)'" != "" {
        display as text "Note: no CFI or TLI is defined for a model with a covariate, and RMSEA is not reported in this version."
    }
    * release A, A3.2 and A3.4: said in error colour, with the reason (A7 wording, owner 2026-11-18;
    * critical is status 2, which also covers a missing value; acov is status 3, dwlsmv_m11_efa_se)
    if "`e(efa)'" != "" & "`e(rot_nose)'" != "" {
        display as error "No standard errors are reported for the rotated estimates; e(V) is not stored."
        if "`e(rot_nose)'" == "unconverged" {
            display as error "Reason: the unrotated fit did not converge."
        }
        else if "`e(rot_nose)'" == "inadmissible" {
            display as error "Reason: the unrotated fit is inadmissible (a communality exceeds 1)."
        }
        else if "`e(rot_nose)'" == "rank" {
            display as error "Reason: the weighted Jacobian of the rotated model is rank deficient"
            display as error "(smallest to largest singular value below 1e-8)."
        }
        else if "`e(rot_nose)'" == "critical" {
            display as error "Reason: the rotated covariance cannot be computed: the geomin criterion has a"
            display as error "degenerate critical point there, or a needed matrix has a missing value."
        }
        else {
            display as error "Reason: the covariance of the unrotated estimates is not available (it holds a"
            display as error "missing value)."
        }
    }
end
