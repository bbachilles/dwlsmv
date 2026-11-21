*! version 0.6.0  dwlsmv -- Module 11: sem's reading of the model, as the canonical form
*!
*! Internal; called by dwlsmv.ado. Spec Module 11 D1, D2, D9.
*!
*!   _dwlsmv_translate , paths(<sem paths>) [semopts(<sem options>) metric(marker|unitvariance)
*!                       group(<varname>) gvals(<sorted group values>)]
*!
*! With group() (invariance() only, 2026-10-12) the dummy data carry the
*! group variable with the user's values, sem runs with ginvariant(none), and
*! r(canon) is the multiple-group configural text (dwlsmv_translate_g), every
*! factor scaled by var(F@1) injected in every group.
*!
*! Runs sem, noestimate on a dummy dataset holding only the model's observed
*! variables (N = 200, one factor, asserted benign), so the canonical form
*! cannot depend on the user's data. Under metric(unitvariance) a forcenoanchor
*! pass first finds the factors the user did not normalise, and the model is
*! run again with var(F@1) / var(e.F@1) injected for them (D2). Under
*! metric(marker) the ordinary pass is read first, and only a factor with no
*! observed indicator and a path out of it (a second-order factor) that the user
*! did not normalise gets var(G@1) injected (release B, B3.2); sem then runs
*! again and one Note line says so. The user's
*! data, sort order, RNG state and estimation results are restored on every
*! exit, error or not. sem's refusals reach the user as sem's message and rc.
*!
*! Returns
*!   r(canon)     the canonical lines, separated by ";"
*!   r(ind)       the indicators in e(oyvars) order (D9)
*!   r(semnames)  sem's structural stripe, eq:name, for D6's check
*!   r(injected)  what D2 injected, if anything
*!   r(sofac)     the factors with no observed indicator and a path out of them
*!                (second-order factors, release B), or empty
*!   r(cov)       the observed exogenous variable (the binary covariate of a conditional
*!                fit), or empty; the canonical text then carries a covariate line
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
program _dwlsmv_translate, rclass
    version 16
    syntax , PATHS(string asis) [SEMopts(string asis) METRIC(string) ///
        GROUP(name) GVALS(numlist integer >=0 sort)]
    if "`metric'" == "" local metric marker
    local G : word count `gvals'

    * the model's observed variables: tokens of the text that name variables
    * in the data. Capitalized names are latent to sem unless nocapslatent.
    local caps = strpos(`"`semopts'"', "nocapslatent") > 0
    local txt `"`paths' `semopts'"'
    local txt : subinstr local txt "->" " ", all
    local txt : subinstr local txt "<-" " ", all
    foreach ch in "(" ")" "@" "*" "," "=" "." {
        local txt : subinstr local txt "`ch'" " ", all
    }
    local obs
    foreach t of local txt {
        if !`caps' & ustrregexm("`t'", "^[A-Z]") continue
        capture unab u : `t'
        if !_rc local obs `obs' `u'
    }
    local obs : list uniq obs
    local p : word count `obs'
    if `p' == 0 {
        display as error "dwlsmv: the model names no observed variable in the data"
        exit 198
    }

    local rng `c(rngstate)'
    tempname est
    _estimates hold `est', restore nullok
    preserve
    capture noisily {
        drop _all
        quietly set obs `= 200 * max(`G', 1)'
        set seed 20261003
        * under invariance(), the group variable carries the user's values
        * (sorted, O13), 200 dummy observations each, so sem's stripe names
        * the groups as the user's data will
        if `G' {
            quietly generate long `group' = 0
            local i 0
            foreach v of local gvals {
                quietly replace `group' = `v' if mod(_n - 1, `G') == `i++'
            }
        }
        * integer codes 0-3 from a one-factor latent: a group() variable the
        * scan picks up is then a valid sem group variable with four groups
        * of about 50, so Module 11's own multi-group refusal is what the user
        * sees, not sem's complaint about non-integer dummy values (spec D1)
        tempvar f z
        quietly generate double `f' = rnormal()
        foreach v of local obs {
            quietly generate double `z' = 0.6 * `f' + 0.8 * rnormal()
            quietly generate double `v' = (`z' > -0.6) + (`z' > 0) + (`z' > 0.6)
            drop `z'
        }
        drop `f'
        mata: dwlsmv_benign(tokens(st_local("obs")))
        local inj
        local soinj
        local sofac
        if `G' {
            * multiple groups, invariance() only (decided 2026-10-12): sem's
            * group grammar with ginvariant(none) (spec principle 4); every
            * factor without var(F@1) in group 1 gets it injected, in every
            * group; a loading fixed at 1 is refused on the final pass
            local gopt group(`group') ginvariant(none)
            quietly sem `paths', `semopts' `gopt' forcenoanchor noestimate
            mata: st_local("inj", dwlsmv_inject_g(dwlsmv_translate_g(0, ///
                "`group'", strtoreal(tokens("`gvals'")))))
            quietly sem `paths', `semopts' `gopt' `inj' noestimate
            mata: st_local("canon", invtokens(dwlsmv_translate_g(1, ///
                "`group'", strtoreal(tokens("`gvals'")))', ";"))
            mata: st_local("semnames", invtokens(dwlsmv_semnames()'))
            local ind `e(oyvars)'
        }
        else {
            if "`metric'" == "unitvariance" {
                quietly sem `paths', `semopts' forcenoanchor noestimate
                mata: st_local("inj", dwlsmv_inject(dwlsmv_translate(0)))
            }
            else {
                * release B (B3.2): a factor with no observed indicator is scaled
                * by its variance. The ordinary pass is read first; sem runs again
                * only when something is injected
                quietly sem `paths', `semopts' noestimate
                mata: st_local("inj", dwlsmv_inject_so(dwlsmv_translate(0)))
                local soinj `inj'
            }
            if "`metric'" == "unitvariance" | `"`inj'"' != "" {
                quietly sem `paths', `semopts' `inj' noestimate
            }
            mata: st_local("canon", invtokens(dwlsmv_translate(1)', ";"))
            mata: st_local("sofac", invtokens(dwlsmv_tr_sofac_t(st_local("canon"))))
            local cov `e(oxvars)'
            mata: st_local("semnames", invtokens(dwlsmv_semnames()'))
            local ind `e(oyvars)'
        }
    }
    local rc = _rc
    restore
    set rngstate `rng'
    if `rc' exit `rc'
    * B3.2's line, only after a successful translation that injected var(G@1)
    if `"`soinj'"' != "" {
        local sonames
        foreach t of local soinj {
            local t : subinstr local t "var(e." "", all
            local t : subinstr local t "var(" "", all
            local t : subinstr local t "@1)" "", all
            local sonames `sonames' `t'
        }
        if `: word count `sonames'' == 1 {
            display as text `"Note: `sonames' has no observed indicator and is scaled by its variance, `soinj'; sem's default fixes the path to its first child at 1 instead. Standardized estimates are the same."'
        }
        else {
            display as text `"Note: `sonames' have no observed indicator and are scaled by their variances, `soinj'; sem's default fixes each factor's path to its first child at 1 instead. Standardized estimates are the same."'
        }
    }

    return local canon    `"`canon'"'
    return local ind      `"`ind'"'
    return local semnames `"`semnames'"'
    return local injected `"`inj'"'
    return local cov      `"`cov'"'
    return local sofac    `"`sofac'"'
end
