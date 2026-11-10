*! dwlsmv -- build the Mata library ldwlsmv.mlib from src/mata
*!
*! Packaging (decided 2026-10-04, the user's): the command is distributed, so
*! its Mata code ships as a compiled library that dwlsmv.ado finds on the
*! adopath with no per-call loading. This file is the only place the library
*! is built. It clears Mata (it is the caller; no definition file clears),
*! loads dwlsmv_types.mata first -- every struct must exist before a function
*! that uses it compiles -- then every other .mata file in src/mata, and writes
*! build/ldwlsmv.mlib, which .gitignore excludes.
*!
*! Run from the repository root:
*!     StataMP-64.exe /e do tools/build_mlib.do
*! The Module 11 tests run it themselves before anything else, so a stale
*! library cannot pass a test. A library built by Stata 18 does NOT load in
*! Stata 16 (measured 2026-10-05 by tools/probe_mlib.do on StataMP 16.1: it
*! is ignored as "compiled by Stata 18.0, is too new", rc 3499), so the
*! distributed library is built on the Stata 16 milestone machine. A Stata 16
*! build does load and run in Stata 18, with fits identical to Stata 18's own
*! build (tools/probe_mlib16.do, 2026-10-05); Stata 17 is not yet tested. Mata finds a
*! library only after "mata: mata mlib index", which a session that puts
*! build/ on the adopath after startup must run itself. See docs/decisions.md
*! and docs/milestone.md (2026-10-05).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek
version 16
mata: mata clear
local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
quietly do src/mata/dwlsmv_types.mata
foreach f of local srcs {
    quietly do "src/mata/`f'"
}
capture mkdir build
mata: mata mlib create ldwlsmv, dir("build") replace
mata: mata mlib add ldwlsmv *(), dir("build") complete
mata: mata mlib index
display "build_mlib: build/ldwlsmv.mlib from " (`: word count `srcs'' + 1) " source file(s)"
