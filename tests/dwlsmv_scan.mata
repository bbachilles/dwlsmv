*! dwlsmv -- source scanner for the version-floor allowlist check
*!
*! Extracts the names of functions CALLED in a Mata source file, so
*! tests/unit/test_version_floor.do can compare them against a reviewed
*! allowlist. See that file and E-072 for why this exists: the
*! Stata 16 floor cannot be enforced on the development machine, and the
*! milestone machine is only synced periodically, so a post-16 routine can sit
*! green here for a whole milestone. This closes that gap at write time.
*!
*! This file is deliberately NOT scanned by the check it supports. It is test
*! infrastructure rather than shipped code, and a scanner that must appear in
*! its own allowlist is a circularity with no benefit.
*!
*! The scanner is textual and therefore approximate. That is the right trade:
*! a false positive costs one allowlist line, while the thing it must not do
*! is miss a real call. It handles // and /* */ comments, Stata-style * lines,
*! and double-quoted strings. It does not parse Mata.
*!
*! Only pre-16 Mata is used here -- fopen/fget/fclose, strlen, substr, strpos,
*! strtrim, uniqrows, sort -- because this file runs on the milestone machine
*! too, and a floor check that breaks the floor would be a poor joke.
*!
*! This file does NOT issue `mata clear`, and must not. A definition file that
*! clears can only ever be loaded first: whatever was loaded before it is
*! silently wiped. Clearing belongs to the caller, and every test here opens
*! with `clear all`, which clears Mata too -- verified, not assumed. See
*! E-072, 2026-09-18.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Is c the start of an identifier?
real scalar dwlsmv_isidstart(string scalar c)
{
    if (c == "") return(0)
    return(strpos("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_", c) > 0)
}

// Is c a continuation character of an identifier?
real scalar dwlsmv_isidchar(string scalar c)
{
    if (c == "") return(0)
    return(strpos("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789", c) > 0)
}

// Read a file, dropping lines whose first non-blank character is * (a Stata
// comment, which is how these files carry their headers). Returns one string
// with the surviving lines joined by newlines.
string scalar dwlsmv_readsrc(string scalar fn)
{
    real scalar fh
    string scalar line, out, t

    out = ""
    fh  = fopen(fn, "r")
    while ((line = fget(fh)) != J(0, 0, "")) {
        t = strtrim(line)
        if (substr(t, 1, 1) != "*") {
            out = out + line + char(10)
        }
    }
    fclose(fh)
    return(out)
}

// Blank out // comments, /* */ comments and "..." string literals, so that a
// parenthesis inside any of them cannot be read as a call.
string scalar dwlsmv_decomment(string scalar src)
{
    real scalar i, n, st
    string scalar out, c, c2, q, nl

    q   = char(34)
    nl  = char(10)
    out = ""
    st  = 0
    n   = strlen(src)
    i   = 1

    while (i <= n) {
        c  = substr(src, i, 1)
        c2 = substr(src, i + 1, 1)

        if (st == 0) {
            if (c == q) {
                st = 1
            }
            else if (c == "/" & c2 == "/") {
                st = 2
                i  = i + 1
            }
            else if (c == "/" & c2 == "*") {
                st = 3
                i  = i + 1
            }
            else {
                out = out + c
            }
        }
        else if (st == 1) {
            if (c == q) st = 0
        }
        else if (st == 2) {
            if (c == nl) {
                st  = 0
                out = out + c
            }
        }
        else {
            if (c == "*" & c2 == "/") {
                st = 0
                i  = i + 1
            }
        }
        i = i + 1
    }
    return(out)
}

// Every identifier that is immediately followed by "(", sorted and deduped.
// Whitespace between the name and the paren is allowed, since Mata permits it.
string colvector dwlsmv_callnames(string scalar src)
{
    real scalar i, n, j
    string scalar c, name
    string colvector found

    found = J(0, 1, "")
    n     = strlen(src)
    i     = 1

    while (i <= n) {
        c = substr(src, i, 1)
        if (dwlsmv_isidstart(c)) {
            name = ""
            while (i <= n & dwlsmv_isidchar(substr(src, i, 1))) {
                name = name + substr(src, i, 1)
                i    = i + 1
            }
            j = i
            while (j <= n & (substr(src, j, 1) == " " | substr(src, j, 1) == char(9))) {
                j = j + 1
            }
            if (substr(src, j, 1) == "(") {
                found = found \ name
            }
        }
        else {
            i = i + 1
        }
    }

    if (rows(found) == 0) return(found)
    return(uniqrows(sort(found, 1)))
}

// Convenience: scan one file.
string colvector dwlsmv_scanfile(string scalar fn)
{
    return(dwlsmv_callnames(dwlsmv_decomment(dwlsmv_readsrc(fn))))
}

end
