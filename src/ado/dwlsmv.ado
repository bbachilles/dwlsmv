*! version 0.4.1  dwlsmv -- SEM/CFA for ordinal indicators by DWLS with the WLSMV test
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do;
*! standardized on 2026-10-05: tests/unit/test_standardized.do.
*!
*!   dwlsmv <sem paths> [if] [in] [, metric(marker|unitvariance) keepbound
*!                                   standardized level(#) <sem options>]
*!   dwlsmv (F1 F2 -> y1-y6) [if] [in], efa(F1 F2) [epsilon(#) rstarts(#)
*!                                   keepbound]           (EFA, 0.3.0)
*!   dwlsmv [, standardized level(#)]               (replay)
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
*! Spec Module 11, D1-D9. sem reads the model (D1, _dwlsmv_translate); dwlsmv
*! fits it. Return codes (D8a): Module 11's own refusals rc 198; a Stage 1
*! abort is caught under capture noisily, followed by the legend and the
*! pointer line, and rethrown as rc 3498; any later Mata abort is rethrown as
*! it is. A fit that does not converge is posted with e(converged) = 0 and a
*! note, as sem does (decided 2026-10-04). EXCEPTION, a model with a covariate
*! (confirmed by the owner 2026-11-14): e(b) only, no e(V), no e(T3), no e(p_T3); the
*! header and the table say so and the replay says the same.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
program dwlsmv, eclass
    version 16
    if replay() {
        if "`e(cmd)'" != "dwlsmv" error 301
        syntax [, STANDardized Level(cilevel)]
        if "`standardized'" != "" & "`e(efa)'" != "" {
            display as error "standardized is not available with efa(): the rotated loadings are already standardized."
            exit 198
        }
        if "`standardized'" != "" & "`e(covariates)'" != "" {
            display as error "standardized is not available with a covariate: the standardized solution of a model"
            display as error "conditional on `e(covariates)' has not been derived."
            exit 198
        }
        _dwlsmv_header
        if "`e(efa)'" != "" _dwlsmv_efatable
        else _dwlsmv_table, `standardized' level(`level')
        exit
    }
    syntax anything(name=paths equalok) [if] [in] ///
        [, METRIC(string) KEEPBOUND STANDardized Level(cilevel) ///
        VCE(string) METHod(string) ITERate(string) TECHnique(string) ///
        TOLerance(string) LTOLerance(string) NRTOLerance(string) NONRTOLerance ///
        FROM(string) STARTValues(string) STARTGrouping(string) NOIVSTart ///
        DIFficult GRoup(varname) INVariance(string) PARTIAL(string) ///
        GINVariant(string) EFA(string) EPSilon(string) RSTARTS(string) *]

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
        if `nfix' ereturn post `b' `V' `C', esample(`touse') obs(`N')
        else      ereturn post `b' `V', esample(`touse') obs(`N')
        ereturn matrix K = `K'
        ereturn scalar rank        = `rank'
        ereturn scalar T3          = `T3'
        ereturn scalar df_T3       = `df'
        ereturn scalar p_T3        = cond(`df' > 0, chi2tail(`df', `T3'), .)
        ereturn scalar rmsea_s     = `rmsea_s'
        ereturn scalar cfi_s       = `cfi_s'
        ereturn scalar tli_s       = `tli_s'
        ereturn scalar T3_b        = `T3_b'
        ereturn scalar df_T3_b     = `df_T3_b'
        ereturn scalar converged   = `converged'
        ereturn scalar iterations  = `iterations'
        ereturn scalar bread_ratio = `bratio'
        ereturn local invariance "`invariance'"
        ereturn local metric     "`metric'"
        ereturn local recovery   "none"
        ereturn local indicators "`ind'"
        ereturn local cmdline    `"dwlsmv `0'"'
        ereturn local version    "0.4.1"
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
        capture noisily _dwlsmv_table, `standardized' level(`level')
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
            display as error "Factor `fac' has no fixed scale: fix one of its loadings at 1,"
            display as error "e.g. (`fac' -> `y1'@1 ...), or specify metric(unitvariance)."
        }
        capture mata: dwlsmv_m11_clear()
        exit `rc'
    }

    if `"`efa'"' != "" {
        * rotated point estimates only, no e(V) (the user's choice, 2026-10-29)
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
            exit 3498
        }
        ereturn post `rb', esample(`touse') obs(`N')
        ereturn matrix Lambda = `lam'
        ereturn matrix Phi    = `phi'
        foreach s in starts converged minima atbest Q gpnorm eps gap {
            ereturn scalar rot_`s' = `rot_`s''
        }
        ereturn local efa      "`efa'"
        ereturn local rotation "geomin oblique"
    }
    else {
        if "`cov'" != "" & `converged' == 0 {
            // OUTSTANDING DECISION (owner, 2026-11-14): an ordinary fit without a covariate that does not converge is unchanged in this
            // release (it still posts e(V) and T3); its inconsistency with this branch is an explicit outstanding decision, not a policy.
            // a model with a covariate that did not converge: the last iterate as point estimates only. No
            // e(V) (the sandwich at a non-solution is not inference), as the EFA block posts e(b) only;
            // e(Cns_b) (e(Cns) cannot be set without e(V)) is kept so that the table can mark the fixed parameters
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
    if !("`cov'" != "" & `converged' == 0) {
        ereturn scalar T3     = `T3'
        ereturn scalar p_T3   = cond(`df' > 0, chi2tail(`df', `T3'), .)
    }
    // fit indices, option A (family S, frozen 2026-10-21); printed by _dwlsmv_header.
    // Omitted, not posted as missing, for a conditional fit (decision 2026-11-13 (a)): the
    // header prints them only when posted
    if "`cov'" == "" {
        ereturn scalar rmsea_s    = `rmsea_s'
        ereturn scalar cfi_s      = `cfi_s'
        ereturn scalar tli_s      = `tli_s'
        ereturn scalar T3_b       = `T3_b'
        ereturn scalar df_T3_b    = `df_T3_b'
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
    ereturn local version  "0.4.1"
    // what difftest refits from (spec Module 12, D8d addition, 2026-10-23):
    // the canonical text as dwlsmv_m11_fit received it, keepbound, and the
    // data signature that check 2 tests with checkestimationsample
    ereturn local canon     `"`canon'"'
    ereturn local keepbound "`keepbound'"
    signestimationsample `ind' `cov'
    ereturn local cmd      "dwlsmv"

    _dwlsmv_header
    if `"`efa'"' != "" _dwlsmv_efatable
    else _dwlsmv_table, `standardized' level(`level')
    mata: dwlsmv_m11_notes()
    capture mata: dwlsmv_m11_clear()
end

* EFA (0.3.0, Unit 6): rotated loadings and factor correlations, point
* estimates only; from e() alone, so replay works
program _dwlsmv_efatable
    tempname L P
    matrix `L' = e(Lambda)
    matrix `P' = e(Phi)
    display _n as text "Rotation: geomin, oblique (epsilon = " as result e(rot_eps) ///
        as text ", " as result e(rot_starts) as text " starts, " ///
        as result e(rot_converged) as text " converged)"
    matlist `L', title("Rotated loadings") format(%9.3f) border(rows) lines(oneline)
    matlist `P', title("Factor correlations") format(%9.3f) border(rows) lines(oneline)
    display as text "Standard errors of the rotated estimates are not available in this version."
    if e(rot_minima) > 1 {
        display as text "The rotation found " as result e(rot_minima) as text " distinct minima; " ///
            as result e(rot_atbest) as text " of the converged starts reached the reported one" ///
            " (criterion gap to the next " as result %9.2e e(rot_gap) as text ")."
    }
    if e(df_T3) == 0 {
        display as text "With df = 0 the model fits exactly; T3 and the fit indices are not informative."
    }
end

* A model with a covariate that did not converge: the last iterate as a plain table, from e(b) and
* e(Cns_b) alone (replay works). No standard errors, z, p or interval exist to show: e(V) is not posted.
program _dwlsmv_nctable
    tempname b C M F
    matrix `b' = e(b)
    local nb = colsof(`b')
    matrix `F' = J(`nb', 1, 0)
    capture matrix `C' = e(Cns_b)
    if !_rc {
        forvalues k = 1/`=rowsof(`C')' {
            forvalues i = 1/`nb' {
                if `C'[`k', `i'] == 1 matrix `F'[`i', 1] = 1
            }
        }
    }
    matrix `M' = `b'' , `F'
    matrix colnames `M' = Estimate Fixed
    matlist `M', title("Last iterate, not a solution; no standard errors (Fixed = 1: set by the model)") ///
        format(%12.5g) border(rows) lines(oneline)
end

* D7 step (3): standardized replaces the coefficient table; Theta rows appear
* only here; nothing is stored (D8d). Computed from e() alone, so replay
* works; the posted results are held and restored around the display.
program _dwlsmv_table
    syntax [, STANDardized Level(cilevel)]
    * a model with a covariate that did not converge: point estimates only, no e(V) to display
    if "`e(covariates)'" != "" & e(converged) == 0 {
        _dwlsmv_nctable
        exit
    }
    if "`standardized'" == "" {
        ereturn display, level(`level')
        exit
    }
    tempname hold bs Vs Cs
    * at G > 1 the group routine fails closed before anything is shown (O17)
    local inv "`e(invariance)'"
    if "`inv'" != "" mata: dwlsmv_m11_std_g("`bs'", "`Vs'", "`Cs'")
    else             mata: dwlsmv_m11_std("`bs'", "`Vs'", "`Cs'")
    _estimates hold `hold', copy restore
    if `nfixs' ereturn post `bs' `Vs' `Cs'
    else       ereturn post `bs' `Vs'
    display as text "Standardized (std.all); /var(e.y) rows are the derived residual variances"
    ereturn display, level(`level')
    * Q-C (a), user 2026-11-01: one note where parameters are tied
    if !inlist("`inv'", "", "configural") {
        display as text "Note: equal raw parameters can have different standardized values across groups."
    }
end

program _dwlsmv_header
    display _n as text "DWLS estimation, WLSMV test (ordinal indicators)" ///
        _col(51) "Number of obs" _col(67) "= " as result %9.0f e(N)
    * a model with a covariate that did not converge (confirmed by the owner 2026-11-14): said first, in error colour
    if "`e(covariates)'" != "" & e(converged) == 0 {
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
        else {
            display as error "No further remedy is built for a model with a covariate;"
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
    if "`e(covariates)'" != "" & e(converged) == 0 {
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
    if "`e(rmsea_s)'" != "" {
        display as text _col(51) "RMSEA (scaled)" _col(67) "= " as result %9.4f e(rmsea_s)
        display as text _col(51) "CFI (scaled)"   _col(67) "= " as result %9.4f e(cfi_s)
        display as text _col(51) "TLI (scaled)"   _col(67) "= " as result %9.4f e(tli_s)
    }
    if "`e(covariates)'" != "" {
        display as text "Note: no CFI or TLI is defined for a model with a covariate, and RMSEA is not reported in this version."
    }
    if e(converged) == 0 & "`e(covariates)'" == "" {
        display as error "convergence not achieved"
    }
end
