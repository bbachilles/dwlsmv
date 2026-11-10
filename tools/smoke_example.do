*! dwlsmv -- smoke test: build the library and the example data, then fit.
*!
*! Not a test (no pass/fail, and the suite runner ignores it): a minimal
*! end-to-end check that this Stata builds ldwlsmv.mlib, finds it, and fits the
*! help file's example model, plain and with `standardized`. It rebuilds
*! build/ldwlsmv.mlib, so do not run it where a library built by another Stata
*! version must be kept.
*!
*! Run from the repository root:
*!     StataMP-64.exe /e do tools/smoke_example.do
*!
*! The `mata: mata mlib index` after the adopath lines is required: Mata finds a
*! library only after indexing, and a session that adds build/ to the adopath
*! after startup must run it. Without it the fit stops at rc 3499
*! (dwlsmv_benign() not found), which is what the first version of this file did.
*!
*! RESULT, 2026-10-05, StataMP 16.1 (serial [milestone], 4 processors): the
*! library built (92 functions), `mata which dwlsmv_benign()` reports ldwlsmv,
*! the example fit ran (N = 5000, T3 chi2(8) = 5.154, p = 0.7410) and so did
*! `dwlsmv, standardized`.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
clear all
do tools/build_mlib.do
do tools/build_example.do
adopath ++ "`c(pwd)'/build"
adopath ++ "`c(pwd)'/src/ado"
mata: mata mlib index
mata: mata which dwlsmv_benign()
use build/dwlsmv_example.dta, clear
dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)
dwlsmv, standardized
