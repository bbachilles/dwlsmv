*! Version-floor check: every external Mata routine called in the scanned
*! sources must appear in tests/mata_allowlist.txt.
*!
*! This does NOT prove Stata 16 compatibility -- only a run on the milestone
*! machine does that, and the development record holds the protocol. What it does
*! is shrink the window. Without it, a post-16 call compiles clean here and
*! stays green until the next sync, so the blast radius is a whole milestone's
*! work. With it, the blast radius is one commit.
*!
*! A failure here is not a bug. It means a routine is being used that nobody
*! has checked against the floor yet. Look it up once, record the answer in
*! the allowlist, move on.
*!
*! KNOWN BEHAVIOUR, verified rather than assumed: the scanner is textual and
*! cannot tell a definition from a call, so a function DEFINED in a scanned
*! file is reported as a call to itself. The project's own names are excluded
*! by prefix (dwls_, dwlsmv_, _dwlsmv_), so in practice this is invisible --
*! but define a function in src/mata/ under some other name and this test will
*! flag it. That is a feature rather than a wart: it makes the check enforce
*! the naming convention as a side effect. Do not "fix" it by adding the
*! function to the allowlist, which is for external routines; rename it to
*! carry the project prefix instead.
*!
*! Verified 2026-09-18 by adding a throwaway src/mata file defining
*! _tmp_probe() and calling sqrt(exp(x)). All three names were reported, the
*! test exited 9, and the runner counted a failure. The pass path and the fail
*! path have both now run.
*!
*! This test has no module dependency, so it never legitimately skips: it
*! passes (0) or fails (9).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

capture noisily do tests/dwlsmv_scan.mata
if _rc {
    display as error "test_version_floor: scanner failed to load"
    exit 9
}

* Scan every Mata source under src/, plus the comparator, which must also run
* on the milestone machine. tests/dwlsmv_scan.mata is deliberately excluded --
* see its header.
local targets ""
local found : dir "src/mata" files "*.mata"
foreach f of local found {
    local targets `"`targets' src/mata/`f'"'
}
local targets `"`targets' tests/dwlsmv_assert.mata"'

mata:
// Mata keywords and type words that the textual scanner cannot distinguish
// from calls, because `if (`, `for (` and `return(` look identical to one.
KW = ("if", "else", "for", "while", "do", "return", "break", "continue",
      "switch", "case", "default", "goto", "pragma", "version", "external",
      "static", "void", "real", "complex", "string", "numeric",
      "transmorphic", "pointer", "struct", "class", "scalar", "vector",
      "rowvector", "colvector", "matrix")'

// Read the allowlist: first field of each non-comment line.
ALLOW = J(0, 1, "")
fh = fopen("tests/mata_allowlist.txt", "r")
while ((line = fget(fh)) != J(0, 0, "")) {
    t = strtrim(line)
    if (t != "" & substr(t, 1, 1) != "#") {
        p = strpos(t, "|")
        if (p > 0) t = strtrim(substr(t, 1, p - 1))
        ALLOW = ALLOW \ t
    }
}
fclose(fh)

files = tokens(st_local("targets"))
bad   = J(0, 1, "")
nscan = 0

for (f = 1; f <= cols(files); f++) {
    nm = dwlsmv_scanfile(files[f])
    nscan = nscan + rows(nm)
    for (i = 1; i <= rows(nm); i++) {
        n = nm[i]
        if (sum(KW :== n) > 0)                continue
        if (substr(n, 1, 5)  == "dwls_")      continue
        if (substr(n, 1, 7)  == "dwlsmv_")    continue
        if (substr(n, 1, 8)  == "_dwlsmv_")   continue
        if (sum(ALLOW :== n) > 0)             continue
        bad = bad \ (n + "   [" + files[f] + "]")
    }
}

printf("scanned %f file(s), %f distinct call site name(s)\n", cols(files), nscan)
printf("allowlist holds %f entr(y/ies)\n", rows(ALLOW))

if (rows(bad) > 0) {
    printf("\nNOT ON THE ALLOWLIST:\n")
    for (i = 1; i <= rows(bad); i++) printf("    %s\n", bad[i])
}
st_local("nbad", strofreal(rows(bad)))
end

if `nbad' {
    display as error _n "test_version_floor: `nbad' routine(s) not on the allowlist."
    display as error "Check `help mata <name>' for a 'new in Stata NN' note, then add"
    display as error "a line to tests/mata_allowlist.txt recording what you found."
    exit 9
}

display as result "test_version_floor: all external Mata calls are on the allowlist"
exit 0
