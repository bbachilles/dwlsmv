*! dwlsmv -- Module 11: D2's var(F@1) / var(e.F@1) injection
*!
*! Executed under StataMP 18 on 2026-10-04: tests/unit/test_translate.do.
*! Called by _dwlsmv_translate.ado on the forcenoanchor pass's canonical text.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// D2: factors the user normalised -- a metric line, a fixed loading, or a
// fixed variance -- from a forcenoanchor pass, P; every other factor gets
// var(F@1), or var(e.F@1) when it has a beta line (endogenous).
string scalar dwlsmv_inject(string colvector P)
{
    string rowvector fac, endo, t
    string scalar inj
    real scalar i, k, norm

    fac  = tokens(substr(P[1], 8, .))
    endo = J(1, 0, "")
    for (i = 3; i <= rows(P); i++) {
        t = tokens(P[i])
        if (t[1] == "beta") {
            endo = endo, t[2]
        }
    }
    // sem's own reading too: a path fixed at 0 leaves no beta line, but sem
    // still treats the factor as endogenous and refuses var(F@1) on it
    endo = endo, tokens(st_global("e(lyvars)"))
    inj = ""
    for (k = 1; k <= cols(fac); k++) {
        norm = 0
        for (i = 3; i <= rows(P); i++) {
            // nested, not &: Mata's & does not short-circuit, and a metric
            // line has three tokens where lambda and psi lines have four
            t = tokens(P[i])
            if (t[1] == "metric") {
                if (t[3] == fac[k]) {
                    norm = 1
                }
            }
            else if (t[1] == "lambda" | t[1] == "psi") {
                if (t[2] == fac[k] & t[4] == "fixed") {
                    if (t[1] == "lambda" | t[3] == fac[k]) {
                        norm = 1
                    }
                }
            }
        }
        if (!norm) {
            inj = inj + (anyof(endo, fac[k]) ? " var(e." : " var(") + fac[k] + "@1)"
        }
    }
    return(strtrim(inj))
}

end
