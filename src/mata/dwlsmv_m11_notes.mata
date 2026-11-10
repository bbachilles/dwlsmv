*! dwlsmv -- Module 11: prints the notes block (D8e)
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_dwlsmv.do.
*! Line by line: one long printf() wraps mid-word in a log.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_notes()
{
    pointer() scalar p
    string scalar s
    real scalar i

    // line by line: one long printf() wraps mid-word in a log
    p = findexternal("dwlsmv_m11_txt")
    if (p != NULL) {
        s = *p
        while ((i = strpos(s, char(10))) > 0) {
            printf("%s\n", substr(s, 1, i - 1))
            s = substr(s, i + 1, .)
        }
        if (s != "") {
            printf("%s\n", s)
        }
    }
}

end
