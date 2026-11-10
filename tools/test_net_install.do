*! dwlsmv -- test a net install of the released package, from a URL or folder.
*!
*! Executed 2026-10-06 from the local folder dist/ at aca3ac8 on StataMP 18 and
*! StataMP 19: every rc 0, T3 chi2(8) = 5.1542, p 0.7410. NOT YET EXECUTED
*! against the tag URL. Written 2026-10-06 from the 2026-10-05 probes (net install
*! then a fresh session; PLUS set before any Mata call; `net get` for the example
*! data). Not a suite test: no pass/fail, and the runner ignores it.
*!
*! Run in a FRESH session, and set the source first. For a tag vX.Y.Z:
*!     global DWLSMV_FROM "https://raw.githubusercontent.com/bbachilles/dwlsmv/vX.Y.Z/dist/"
*!     do tools/test_net_install.do
*! or, from a local staging folder: global DWLSMV_FROM "C:\path\to\dist"
*!
*! It points PLUS and PERSONAL at a throwaway folder BEFORE touching Mata, because
*! Mata builds its library list on first use (a session that used Mata before
*! `sysdir set` gets rc 111 until `mata: mata mlib index`). It installs, checks
*! that the command and the library are found, runs `net get` for the example
*! data into a throwaway working directory, fits the help file's example model
*! and `standardized`, and prints NETINST| lines. Read them: every rc should be 0,
*! and the example fit should give T3 chi2(8) = 5.154, p = 0.7410 (the 2026-10-05
*! smoke runs on Stata 16.1 and 18).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
if `"$DWLSMV_FROM"' == "" {
    display as error "set the source first: global DWLSMV_FROM ""<url or folder>"""
    exit 198
}
local tmp "`c(tmpdir)'/dwlsmv_netinst"
capture mkdir "`tmp'"
capture mkdir "`tmp'/plus"
capture mkdir "`tmp'/personal"
capture mkdir "`tmp'/work"
sysdir set PLUS     "`tmp'/plus"
sysdir set PERSONAL "`tmp'/personal"
cd "`tmp'/work"
display "NETINST|env|stata_version|`c(stata_version)'|MP|`c(MP)'|os|`c(os)'"
display `"NETINST|from|$DWLSMV_FROM"'

capture noisily net install dwlsmv, from(`"$DWLSMV_FROM"') replace
display "NETINST|install_rc|" _rc
capture noisily which dwlsmv
display "NETINST|which_rc|" _rc
capture noisily which difftest
display "NETINST|which_difftest_rc|" _rc
capture noisily findfile difftest.sthlp
display "NETINST|findfile_difftest_help_rc|" _rc
capture noisily mata: mata which dwlsmv_benign()
display "NETINST|mata_which_rc|" _rc
capture noisily findfile ldwlsmv.mlib
display "NETINST|findfile_mlib_rc|" _rc
capture noisily net get dwlsmv, from(`"$DWLSMV_FROM"')
display "NETINST|net_get_rc|" _rc
capture noisily use dwlsmv_example.dta, clear
display "NETINST|use_example_rc|" _rc
capture noisily dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)
display "NETINST|fit_rc|" _rc
if !_rc {
    display "NETINST|N|" e(N)
    display "NETINST|T3|" %9.4f e(T3) " df " e(df_T3) " p " %6.4f e(p_T3)
}
capture noisily dwlsmv, standardized
display "NETINST|standardized_rc|" _rc
* difftest ships in the package from 0.3.0: a nested pair through the installed command
capture noisily estimates store h1
capture noisily dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1@0)
capture noisily estimates store h0
capture noisily difftest h0 h1
display "NETINST|difftest_rc|" _rc
