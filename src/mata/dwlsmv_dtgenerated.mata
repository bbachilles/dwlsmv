*! dwlsmv -- check 8's exemption for generated invariance chains (Module 12)
*!
*!     real scalar dwlsmv_dtgenerated(string scalar cc0, string scalar can0,
*!                                    string scalar inv0, real matrix K0,
*!                                    string scalar cc1, string scalar can1,
*!                                    string scalar inv1, real matrix K1,
*!                                    string scalar name0, string scalar name1)
*!
*! Built 2026-10-13 (exemption gate step 2b; it was a refusing stub from
*! 2026-10-08). Spec 4.8 G8 (what each fit posts) and Module 12 check 8
*! (*Generated invariance chains are exempt*, decided 2026-10-08; non-adjacent
*! pairs (a), condition (4) option (ii), a failed generated claim (a), all
*! 2026-10-13). For each model: cc = e(canon_configural), can = e(canon),
*! inv = e(invariance), K = e(K), all from stored results. name0 and name1 only
*! label the note.
*!
*! Returns 1 (exempt: nesting holds by construction) only if both models post a
*! level and
*!   (1) both post the same e(canon_configural), byte for byte;
*!   (2) both post the same e(K), shape and every element;
*!   (3) dwlsmv_levelgen(cc, K, inv) reproduces each model's e(canon) exactly.
*! Returns 0 otherwise. Condition (4), check 6's refit, is NOT checked here: it
*! is check 6's, which the front end runs on every model before checks 8 and 9
*! (option (ii)). Any two levels from one base qualify, adjacent or not.
*!
*! A FAILED GENERATED CLAIM (option (a)): when both models post a level but a
*! condition fails, one line names the first failed condition, then 0 is
*! returned and check 8 applies as written (at G > 1, its tie refusal). No line
*! when only one model, or neither, posts a level. A configural text the
*! generator rejects aborts there, rc 3498, with the generator's own message.
*!
*! WHAT IT DOES NOT DO. No parameter table, no data, no fit; it does not
*! compute d (check 8 does) or order models (check 9).
*!
*! Load order: dwlsmv_canon_lines, dwlsmv_levelgen.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_dtg_note(string scalar name0, string scalar name1, string scalar why)
{
    errprintf("dwlsmv check 8 (%s, %s): the pair is not treated as a generated invariance\n",
              name0, name1)
    errprintf("chain: %s. Check 8 applies as written.\n", why)
}

// 1 iff the configural text fixes a loading or factor covariance at a nonzero
// value: a fixed row cannot absorb the rescaling that nests the levels (user, 2026-11-02)
real scalar dwlsmv_dtg_fixedrow(string scalar cc)
{
    string colvector L
    string rowvector w
    real scalar i

    L = dwlsmv_canon_lines(cc)
    for (i = 1; i <= rows(L); i++) {
        w = tokens(L[i])
        if (cols(w) >= 5) {
            if (w[4] == "fixed" & strtoreal(w[5]) != 0) {
                if (w[1] == "lambda") {
                    return(1)
                }
                if (w[1] == "psi") {
                    if (w[2] != w[3]) {
                        return(1)
                    }
                }
            }
        }
    }
    return(0)
}

real scalar dwlsmv_dtgenerated(string scalar cc0, string scalar can0,
                               string scalar inv0, real matrix K0,
                               string scalar cc1, string scalar can1,
                               string scalar inv1, real matrix K1,
                               string scalar name0, string scalar name1)
{
    if (inv0 == "" | inv1 == "") {
        return(0)
    }
    if (cc0 != cc1) {
        dwlsmv_dtg_note(name0, name1, "the two fits post different e(canon_configural)")
        return(0)
    }
    if (rows(K0) != rows(K1) | cols(K0) != cols(K1)) {
        dwlsmv_dtg_note(name0, name1, "the two fits post different e(K)")
        return(0)
    }
    if (K0 != K1) {
        dwlsmv_dtg_note(name0, name1, "the two fits post different e(K)")
        return(0)
    }
    if (dwlsmv_dtg_fixedrow(cc0)) {
        dwlsmv_dtg_note(name0, name1, "a fixed nonzero loading or factor covariance is present, and nesting is not established")
        return(0)
    }
    if (dwlsmv_levelgen(cc0, K0, inv0) != can0) {
        dwlsmv_dtg_note(name0, name1, "e(canon) of " + name0 +
                        " is not the generator's text for invariance(" + inv0 + ")")
        return(0)
    }
    if (dwlsmv_levelgen(cc1, K1, inv1) != can1) {
        dwlsmv_dtg_note(name0, name1, "e(canon) of " + name1 +
                        " is not the generator's text for invariance(" + inv1 + ")")
        return(0)
    }
    return(1)
}

end
