*! dwlsmv -- Module 11 D8d probe: does a matrix stripe accept varA~~varB names?
*!
*! A probe, not a test (tools/). D8d (2026-10-04, the user's) posts stage1.bound
*! as a 1 x P matrix e(bound) with variable-name stripes, one column per pair,
*! named varA~~varB, posted always. This probe records, on the machine it runs on,
*! whether matrix colnames accepts such names, whether they survive ereturn
*! matrix, and whether colnumb() finds a column by name:
*!   short  y1~~y2
*!   long   two 32-character variable names, a 66-character element
*!   mixed  underscores and digits, q7_anx~~q9_dep
*!   len32, len33  a varA~~varB element of exactly 32 and 33 characters
*!   slash  /varA~~varB with 32-character names (D5's 46-character stripe had
*!          a leading slash)
*!   eqform varA:varB, the two names as equation and column (no ~~)
*! The first run (StataMP 18) refused the 66-character element as "invalid
*! name", rc 198; the length cases were added after seeing that.
*! Output lines prefixed D8D are the record. Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_d8d_stripe.do      (log: probe_d8d_stripe.log)
*!     grep -E "^(====|D8D)" probe_d8d_stripe.log | diff - tools/probe_d8d_stripe_stata18.txt
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
clear all
set more off
set linesize 200
display "D8DENV stata_version " c(stata_version) " processors " c(processors)
set obs 10
generate byte touse = 1

local a32 "abcdefghijklmnopqrstuvwxyz012345"
local b32 "ABCDEFGHIJKLMNOPQRSTUVWXYZ678901"

program d8dpost, eclass
    args names
    tempname b V B
    * esample() consumes touse; recreate it on every call
    capture drop touse
    generate byte touse = 1
    matrix `b' = (0.5)
    matrix colnames `b' = theta1
    matrix `V' = (0.01)
    matrix rownames `V' = theta1
    matrix colnames `V' = theta1
    ereturn post `b' `V', esample(touse)
    local k : word count `names'
    matrix `B' = J(1, `k', 0)
    matrix `B'[1, 1] = 1
    matrix colnames `B' = `names'
    ereturn matrix bound = `B'
end

local a15 "abcdefghijklmno"
local b15 "ABCDEFGHIJKLMNO"
local b16 "ABCDEFGHIJKLMNOP"
foreach case in short long mixed len32 len33 slash eqform {
    display "==== `case'"
    if "`case'" == "short"  local names "y1~~y2 y1~~y3 y2~~y3"
    if "`case'" == "long"   local names "`a32'~~`b32' y1~~`a32'"
    if "`case'" == "mixed"  local names "q7_anx~~q9_dep q7_anx~~x2 q9_dep~~x2"
    * 15 + 2 + 15 = 32 characters, and 15 + 2 + 16 = 33
    if "`case'" == "len32"  local names "`a15'~~`b15' y1~~y2"
    if "`case'" == "len33"  local names "`a15'~~`b16' y1~~y2"
    * the forms D5's 46-character stripe used: a leading slash, an equation
    if "`case'" == "slash"  local names "/`a32'~~`b32' /y1~~y2"
    if "`case'" == "eqform" local names "`a32':`b32' y1:y2"
    capture noisily d8dpost "`names'"
    display "D8D `case' post rc " _rc
    if _rc == 0 {
        tempname M
        matrix `M' = e(bound)
        local got : colfullnames `M'
        display "D8D `case' names round-trip " ("`got'" == "`names'")
        display "D8D `case' got `got'"
        local first : word 1 of `names'
        display "D8D `case' colnumb first " colnumb(`M', "`first'")
        display "D8D `case' value first " el(`M', 1, colnumb(`M', "`first'"))
    }
}
display "==== done"
