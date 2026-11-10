*! version 0.3.0  difftest -- scaled-and-shifted difference test between nested dwlsmv models
*!
*!   difftest name1 [name2 ...]
*!
*! Written 2026-10-23 (spec Module 12, Architecture item 3 and the validation
*! contract; design docs/notes/level-chain-frontend/design.md, U1-U8; user
*! 2026-10-23: single-group comparisons only). lrtest-style: names are stored
*! estimation results (estimates store), "." the active results; one name is
*! compared with the active results; three or more form a chain, ordered by
*! check 9 and tested as adjacent steps (an extension of lrtest, which takes
*! two). Only dwlsmv results are accepted. Every refusal is rc 3498.
*!
*! This file runs checks 1-5 on the stored results; dwlsmv_dtfront() (Mata)
*! refits every model over one shared Stage 1 and runs checks 6-11, the
*! statistic and the output. The active results are restored afterwards.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

program difftest, rclass
    version 16
    syntax anything(name=names id="model names")
    local names : list retokenize names
    if `: word count `names'' == 1 local names `names' .
    local dups : list dups names
    if "`dups'" != "" {
        display as error "difftest: model `: word 1 of `dups'' is named twice."
        exit 3498
    }

    * the installed dwlsmv release: line 1 of dwlsmv.ado (check 1)
    quietly findfile dwlsmv.ado
    tempname fh
    file open `fh' using `"`r(fn)'"', read text
    file read `fh' line1
    file close `fh'
    local release : word 3 of `line1'

    local act
    if `: list posof "." in names' {
        if "`e(cmd)'" == "" {
            display as error "difftest: there are no active estimation results to use for {bf:.}."
            exit 3498
        }
        tempname act
    }
    * hold before storing "." so the restored results keep their own
    * e(_estimates_name), not the tempname dropped below
    tempname hold
    _estimates hold `hold', copy restore nullok
    if "`act'" != "" quietly estimates store `act'
    capture noisily _difftest_run, names(`names') act(`act') release(`release')
    local rc = _rc
    * taken before the clean-up, which resets r()
    if !`rc' return add
    if "`act'" != "" {
        capture estimates drop `act'
    }
    capture mata: dwlsmv_m11_clear()
    if `rc' exit `rc'
end

program _difftest_run, rclass
    syntax , names(string) release(string) [act(string)]

    mata: dwlsmv_dtfront_init()
    local n : word count `names'
    forvalues k = 1/`n' {
        local nm : word `k' of `names'
        local en = cond("`nm'" == ".", "`act'", "`nm'")
        capture quietly estimates restore `en'
        if _rc {
            display as error "difftest: no stored estimation results named `nm'."
            exit 3498
        }
        if "`e(cmd)'" != "dwlsmv" {
            display as error "difftest: `nm' was not fitted by dwlsmv; difftest compares dwlsmv models only."
            exit 3498
        }
        * EFA (0.3.0, spec section 5 Module 12, Surface): rc 3498,
        * difftest's own convention
        if "`e(efa)'" != "" {
            display as error "difftest: `nm' is an EFA fit (efa(`e(efa)'))."
            display as error "Difference tests with EFA models are not available in this version."
            exit 3498
        }
        * check 1: the installed release
        if "`e(version)'" != "`release'" {
            display as error "difftest check 1 (`nm'): fitted by dwlsmv `e(version)', but dwlsmv `release' is installed."
            display as error "Refit the model with the installed version and store it again."
            exit 3498
        }
        if `"`e(canon)'"' == "" | "`e(datasignature)'" == "" {
            display as error "difftest check 1 (`nm'): the stored results do not carry the specification"
            display as error "difftest refits. Refit the model with the installed dwlsmv and store it again."
            exit 3498
        }
        * check 2: the same data
        capture quietly checkestimationsample
        if _rc {
            display as error "difftest check 2 (`nm'): the data in memory are not the data the model was"
            display as error "fitted to. Load those data, or refit the model on these."
            exit 3498
        }
        tempvar s`k'
        quietly generate byte `s`k'' = e(sample)
        local ind`k' `e(indicators)'
        local kb`k'  `e(keepbound)'
        local met`k' `e(metric)'
        local grp`k' `e(groupvar)'
        mata: dwlsmv_dtfront_add("`nm'")
    }

    local nm1 : word 1 of `names'
    forvalues k = 2/`n' {
        local nmk : word `k' of `names'
        * check 3: the same estimation sample
        capture assert `s`k'' == `s1'
        if _rc {
            display as error "difftest check 3 (`nm1', `nmk'): the models were fitted to different estimation samples."
            exit 3498
        }
        * check 4: the same indicators, keepbound and metric
        if !`: list ind1 === ind`k'' {
            display as error "difftest check 4 (`nm1', `nmk'): the models have different indicators"
            display as error "(`ind1' against `ind`k'')."
            exit 3498
        }
        if "`kb1'" != "`kb`k''" {
            display as error "difftest check 4 (`nm1', `nmk'): keepbound was specified for one model only."
            exit 3498
        }
        if "`met1'" != "`met`k''" {
            display as error "difftest check 4 (`nm1', `nmk'): the models use different metrics (`met1', `met`k'')."
            exit 3498
        }
        if "`grp1'" != "`grp`k''" {
            if "`grp1'" == "" | "`grp`k''" == "" {
                display as error "difftest check 4 (`nm1', `nmk'): one model is a multiple-group fit and the other is not."
            }
            else {
                display as error "difftest check 4 (`nm1', `nmk'): the models use different group variables (`grp1', `grp`k'')."
            }
            exit 3498
        }
        * check 5: the same indicator order
        if "`ind1'" != "`ind`k''" {
            display as error "difftest check 5 (`nm1', `nmk'): the indicators are in a different order"
            display as error "(`ind1' against `ind`k''). Write the paths so that the indicators appear in the"
            display as error "same order in every model; difftest does not reorder them."
            exit 3498
        }
    }

    local gvals
    if "`grp1'" != "" quietly levelsof `grp1' if `s1', local(gvals)
    mata: st_rclear()
    mata: dwlsmv_dtfront("`s1'", "`ind1'", "`kb1'" != "", "`grp1'", "`gvals'")
    * Mata wrote r() directly; an rclass program returns only its return list
    return add
end
