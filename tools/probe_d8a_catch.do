*! dwlsmv -- Module 11 D8a probe: catching a Mata rc 3498 abort at ado level
*!
*! A probe, not a test (tools/). D8a (open) chooses the rc Module 11 exits with
*! after catching a Stage 1 abort, and the spec already has Module 11 catch rc
*! 3498 and print the label map. This probe records, on the machine it runs on,
*! what an ado-level program sees when a Mata routine calls errprintf() and then
*! _error() (rc 3498, the project's abort):
*!   bare     no capture: the caller's capture sees rc 3498
*!   cn       capture noisily: errprintf text and _error text visible? _rc?
*!   cq       capture (quiet): is the errprintf text hidden?
*!   rethrow  capture noisily, print a legend, exit _rc: caller sees 3498?
*!   recode   capture noisily, print a legend, exit 459: caller sees 459?
*!   quiet    quietly { capture noisily ... }: is the text suppressed?
*! Output lines prefixed D8A are the probe's record; message lines carry the
*! markers D8AMSG (errprintf) and D8AERR (_error). Whether they appear in the
*! log is the visibility result. Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_d8a_catch.do      (log: probe_d8a_catch.log)
*!     grep -E "^(====|D8A)|D8AMSG|D8AERR" probe_d8a_catch.log | diff - tools/probe_d8a_catch_stata18.txt
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
clear all
set more off
set linesize 200
display "D8AENV stata_version " c(stata_version) " processors " c(processors)

mata:
void d8a_abort()
{
    errprintf("D8AMSG detailed guidance from errprintf\n")
    _error("D8AERR short abort text")
}
end

program d8a_inner
    args mode
    if "`mode'" == "bare" {
        mata: d8a_abort()
    }
    else if "`mode'" == "cn" {
        capture noisily mata: d8a_abort()
        display "D8A cn inner _rc " _rc
    }
    else if "`mode'" == "cq" {
        capture mata: d8a_abort()
        display "D8A cq inner _rc " _rc
    }
    else if "`mode'" == "rethrow" {
        capture noisily mata: d8a_abort()
        local rc = _rc
        display as error "D8A legend: 1 = y1, 2 = y2"
        exit `rc'
    }
    else if "`mode'" == "recode" {
        capture noisily mata: d8a_abort()
        display as error "D8A legend: 1 = y1, 2 = y2"
        exit 459
    }
    else if "`mode'" == "quiet" {
        quietly {
            capture noisily mata: d8a_abort()
        }
        display "D8A quiet inner _rc " _rc
    }
end

foreach m in bare cn cq rethrow recode quiet {
    display "==== `m'"
    capture noisily d8a_inner `m'
    display "D8A `m' outer _rc " _rc
}
display "==== done"
