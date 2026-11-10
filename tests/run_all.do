*! Run the dwlsmv test suite. Execute from the repository root:
*!
*!     Windows, Git Bash -- the doubled slash is required, because MSYS rewrites a single /e
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tests/run_all.do
*!     cat run_all.log
*!
*!     macOS or Linux
*!     stata -b do tests/run_all.do
*!
*! Verified under StataMP 18 on Windows. List-building, harness-first ordering,
*! and the pass and skip paths ran and counted correctly on 2026-09-17.
*!
*! The FAIL path (rc 9) was exercised on 2026-09-18 with a throwaway test that
*! did nothing but `exit 9`. It behaved: nfail incremented, the test was named
*! in the `failing:` list, and the runner exited 9. An earlier version of this
*! header said that path was not worth manufacturing a failure to prove. It
*! cost one three-line file and took a minute, so that judgement was wrong and
*! is not worth preserving.
*!
*! Still untested: the harness's own abort branch -- the two `exit 9`s that
*! fire when test_harness.do itself skips or fails. Reaching those means
*! breaking the comparator on purpose, which is a larger and less honest
*! forgery than a test that just exits 9, so it is left for the first real
*! occurrence.
*!
*! Tests are ordinary do-files that signal their outcome through the return
*! code, because it is the only channel the runner can read without parsing
*! the log:
*!
*!     0    pass
*!     8    skip -- the module under test does not exist yet
*!     9    fail
*!
*! A skip is not a pass. They are counted separately so an unimplemented
*! module can never be mistaken for a working one; conflating them is how a
*! suite reports green while testing nothing.
*!
*! Return-code 9 is also what `assert` raises, which is deliberate: a bare
*! failed assertion and a deliberate exit 9 should be indistinguishable. 8 was
*! chosen for skip because no natural Stata failure produces it.
*!
*! That claim is now partially checked (StataMP 18, 2026-09-18): `error 8`
*! emits no message text, where an assigned code emits one -- `error 4` prints
*! "no; dataset in memory has changed since last saved". So 8 is unassigned in
*! Stata's error table. That is consistent with nothing natural raising it but
*! does not prove it, since unassigned is not the same as unreachable and no
*! enumeration can settle that. If a test ever skips without an `exit 8`
*! written for it, suspect this before suspecting the runner, and move the
*! skip signal rather than leaving an ambiguous one.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

* ---------------------------------------------------------------------------
* Build the test list.
*
* `: dir` returns filenames ALREADY WRAPPED IN DOUBLE QUOTES. `foreach ... of
* local` strips those when it iterates, but the macro-list operators may not
* match a quoted element against an unquoted literal, so the harness is pulled
* to the front by explicit name comparison rather than by `: list -`. That is
* longer than the one-liner it replaces and does not depend on behaviour this
* file has never been run to confirm.
* ---------------------------------------------------------------------------

local found : dir "tests/unit" files "test_*.do"
local found : list sort found

local harness "test_harness.do"
local rest ""
local sawharness 0

foreach f of local found {
    if "`f'" == "`harness'" {
        local sawharness 1
    }
    else {
        local rest `"`rest' "`f'""'
    }
}

if !`sawharness' {
    display as error "tests/unit/`harness' not found -- refusing to run a suite"
    display as error "whose comparator has not been self-tested"
    exit 9
}

local tests `""`harness'"`rest'"'

local npass 0
local nfail 0
local nskip 0
local failed ""
local skipped ""

foreach t of local tests {
    display _n as text "{hline 70}"
    display as text "running: `t'"
    display as text "{hline 70}"

    capture noisily do "tests/unit/`t'"
    local rc = _rc

    if `rc' == 0 {
        local npass = `npass' + 1
    }
    else if `rc' == 8 {
        local nskip = `nskip' + 1
        local skipped `"`skipped' `t'"'

        * the harness has no dependencies, so it has nothing to skip for
        if "`t'" == "`harness'" {
            display as error _n "harness self-test reported SKIP -- it has no"
            display as error "dependencies and cannot legitimately skip"
            exit 9
        }
    }
    else if `rc' == 9 {
        local nfail = `nfail' + 1
        local failed `"`failed' `t'"'
        if "`t'" == "`harness'" {
            display as error _n "harness self-test failed -- aborting suite"
            exit 9
        }
    }
    else {
        display as error "unexpected return code `rc' from `t'"
        local nfail = `nfail' + 1
        local failed `"`failed' `t'"'
    }
}

display _n as text "{hline 70}"
display as text "passed: `npass'    skipped: `nskip'    failed: `nfail'"

if `nskip' {
    display as text "skipped:`skipped'"
}

if `nfail' {
    display as error "failing:`failed'"
    exit 9
}

if `nskip' {
    display as result "suite green, `nskip' module(s) not yet implemented"
}
else {
    display as result "suite green"
}
