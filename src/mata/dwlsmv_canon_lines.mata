*! dwlsmv -- split a canonical model specification into directive lines (Mod. 3)
*!
*! Executed under StataMP 18 on 2026-09-22. tests/unit/test_identify.do check 7 feeds it comments, blank lines and trailing comments in a shuffled specification.
*!
*! Spec Section 4.8: one directive per line, '#' begins a comment, blank lines
*! are ignored. This returns the surviving lines, trimmed, with comments
*! removed; everything else about the grammar is dwls_identify()'s.
*!
*! CARRIAGE RETURNS ARE STRIPPED EXPLICITLY. strtrim() removes blanks, not
*! char(13), so a specification that reached Mata from a Windows text file
*! would otherwise end its last token with an invisible character and every
*! keyword comparison on that line would fail -- "the directive 'free' is not
*! recognised" on a line that plainly reads 'free'. The specification is
*! normally built in Mata from string literals joined by char(10) and carries
*! no CR at all, which is exactly why this would be found the hard way.
*!
*! TABS ARE CONVERTED TO SPACES, for a sharper version of the same reason, and
*! this one is measured rather than anticipated: **Mata's tokens() does not
*! treat char(9) as a separator.** Verified at a Stata 18 prompt --
*! tokens("lambda" + char(9) + "F1" + char(9) + "y1" + char(9) + "free")
*! returns ONE token, where the space-separated form returns four. Spec
*! Section 4.8 invites canonical specifications written by hand in a do-file,
*! and its own worked example is column-aligned, so an author who aligns with
*! tabs instead of spaces writes a directive that is correct in every visible
*! respect. Without this line dwls_identify() sees the whole line as the
*! directive name and rejects it as unrecognized -- rc 3498, so it aborts
*! cleanly, but the message names a "directive" that is the entire line with
*! tabs in it and nothing in it suggests whitespace is the problem.
*!
*! strtrim() does not remove char(9) either, so trimming does not cover this.
*!
*! Newlines are char(10). '\n' is an escape inside a printf FORMAT only, never
*! inside a string being scanned, and a specification built with "\n" between
*! its lines arrives here as one line containing two literal characters.
*!
*! src/mata/dwlsmv_types.mata is not required by this file, which touches no
*! struct; it is loaded with the rest of Module 3 for convenience.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string colvector dwlsmv_canon_lines(string scalar modelspec)
{
    string colvector out
    string scalar rest, line
    real scalar nl, hash

    out  = J(0, 1, "")
    rest = subinstr(modelspec, char(13), "",  .)
    rest = subinstr(rest,      char(9),  " ", .)

    while (1) {
        nl = strpos(rest, char(10))
        if (nl > 0) {
            line = substr(rest, 1, nl - 1)
            rest = substr(rest, nl + 1, strlen(rest) - nl)
        }
        else {
            line = rest
            rest = ""
        }

        hash = strpos(line, "#")
        if (hash > 0) {
            line = substr(line, 1, hash - 1)
        }

        line = strtrim(line)
        if (line != "") {
            out = out \ line
        }

        if (nl == 0) {
            break
        }
    }

    return(out)
}

end
