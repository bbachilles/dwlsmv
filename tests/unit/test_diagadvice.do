*! Module 10: the prose, dwlsmv_diagadvice(G, C, B, ind, slab, pname).
*!
*! Executed under StataMP 18 on 2026-09-26 and it passes. Before the routine
*! existed it failed with r(3499). Its first real run failed three checks on
*! a bug in its own helper (td_lines consumed its by-reference argument),
*! fixed there, not in the routine under test.
*!
*! Spec Section 5, Module 10; E-060, E-071 and E-075.
*! The flags/prose split of dwlsmv_admissadvice(): a program branches on the
*! matrices, a person reads this. Flag notes by family, then the category
*! counts, ALWAYS, as descriptive diagnostics (2026-10-07, the user's). The
*! singular-value ratios are no longer printed (Module 11 stores them).
*!
*!   1. CASE D, healthy: no flag note, and nothing about a ratio, Gamma, the
*!      bread, a cutoff or "v1": the advice is the category counts, carrying
*!      C's values.
*!   2. SMALL COUNTS ARE NEVER GAMMA CONDITIONING, NOR A WARNING: with
*!      countdiag rows of 1, no count line mentions Gamma, conditioning or
*!      singularity, and the advice has the same number of lines as check 1.
*!   3. EVERY FLAG HAS ITS NOTE, in the 2026-10-07 wording, naming its
*!      location through the labels: Gamma rows 2-5 and bread rows 1-3, one
*!      at a time; a flagged fit prints no ratio either.
*!   4. ASYMMETRY: an internal-check report ("Please report") that says the
*!      remaining Gamma checks were not run.
*!   5. HYGIENE: no literal backslash-n anywhere (a project coding rule), every line ends
*!      in char(10).
*!   6. CONTRACT, rc 3498: wrong shapes, a location outside its labels, and
*!      (after the 2026-09-26 review) a bad parameter location and invalid
*!      category positions. Check 3b, also from the review: the count values,
*!      the constant-column wording, and "evaluated, nothing to report" kept
*!      apart from "not evaluated".
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load tests/dwlsmv_assert.mata src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_diagadvice: `f' failed to load"
        exit 9
    }
}

import delimited using "fixtures/data/case_D.csv", clear varnames(1) delimiters(",")
mata: YD = st_data(., .)

mata:

string scalar td_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The lines of s, split on char(10). Splits a COPY: Mata passes by
// reference, and consuming s itself emptied the caller's string on the
// first run of this test (the advice then printed blank and checks 4-5
// failed on a string that was no longer there).
string colvector td_lines(string scalar s)
{
    string colvector L
    string scalar t
    real scalar k

    t = s
    L = J(0, 1, "")
    while ((k = strpos(t, char(10))) > 0) {
        L = L \ substr(t, 1, k - 1)
        t = substr(t, k + 1, .)
    }
    if (t != "") L = L \ t
    return(L)
}

// Does any line containing `a' also contain `b' (case-insensitive)?
real scalar td_cooccur(string colvector L, string scalar a, string scalar b)
{
    real scalar i

    for (i = 1; i <= rows(L); i++) {
        if (strpos(strlower(L[i]), strlower(a)) & strpos(strlower(L[i]), strlower(b))) {
            return(1)
        }
    }
    return(0)
}

// The index of the first line containing a, or 0.
real scalar td_find(string colvector L, string scalar a)
{
    real scalar i

    for (i = 1; i <= rows(L); i++) {
        if (strpos(L[i], a)) return(i)
    }
    return(0)
}

real scalar dwlsmv_t_diagadvice(real matrix YD)
{
    struct dwlsmv_stage1 scalar st
    struct dwlsmv_parmap scalar par
    struct dwlsmv_fit scalar fit
    string colvector FAIL, L, L0, ind, slab, pname
    string scalar out, out0, cav
    real matrix G, C, B, Gx, Bx, Cx
    real scalar i, k, np

    FAIL = J(0, 1, "")
    cav  = "A small ratio is not by itself evidence that the standard errors are wrong."

    st  = dwls_sample_stats(YD, (4 \ 4 \ 3 \ 3))
    par = dwls_identify(td_spec(("factor f", "ind y1 y2 y3 y4",
                                  "lambda f y2 free", "lambda f y3 free",
                                  "lambda f y4 free", "psi f f free",
                                  "metric marker f y1")), st.map)
    fit = dwls_fit(st, par)
    G   = dwlsmv_gammadiag(st, YD)
    C   = dwlsmv_countdiag(YD, st.map.K)
    B   = dwlsmv_breaddiag(fit.Delta, st.Gamma)
    ind  = ("y1" \ "y2" \ "y3" \ "y4")
    slab = st.map.label
    np   = par.npar
    pname = J(np, 1, "")
    for (i = 1; i <= par.nall; i++) {
        if (par.free[i]) pname[par.fidx[i]] = par.name[i]
    }

    // ---- 1. healthy: no flag note, no ratio, the counts --------------------------
    // The singular-value ratios are stored by Module 11, not printed, since
    // 2026-10-07 (the user's): nothing about them may appear on a clean fit.
    out0 = dwlsmv_diagadvice(G, C, B, ind, slab, pname)
    L0   = td_lines(out0)
    printf("{txt}  1. case D, the whole advice:\n%s", out0)
    for (i = 1; i <= rows(L0); i++) {
        if (strpos(L0[i], "flag") | strpos(L0[i], "not symmetric") |
            strpos(L0[i], "not positive") | strpos(L0[i], "Too few") |
            strpos(L0[i], "move together") | strpos(L0[i], "is constant") |
            strpos(L0[i], "estimating equation") | strpos(L0[i], "insufficient") |
            strpos(L0[i], "zero") | strpos(L0[i], "internal")) {
            FAIL = FAIL \ ("1: a flag note on healthy case D: " + L0[i])
        }
        if (strpos(L0[i], "singular value") | strpos(L0[i], "ratio") |
            strpos(L0[i], "Gamma") | strpos(L0[i], "Bread") | strpos(L0[i], "bread") |
            strpos(L0[i], "cutoff") | strpos(L0[i], "v1")) {
            FAIL = FAIL \ ("1: a diagnostic value or internal label printed on a clean fit: " + L0[i])
        }
    }
    if (td_find(L0, "Category counts") != 1) FAIL = FAIL \ "1: the category counts are not the whole of a clean fit's advice"
    if (!strpos(out0, "lowest or highest category: " + strofreal(C[1, 2]) + " (")) FAIL = FAIL \ "1: the boundary count does not carry C[1, 2]"
    if (!strpos(out0, "middle category:            " + strofreal(C[2, 2]) + " (")) FAIL = FAIL \ "1: the interior count does not carry C[2, 2]"

    // ---- 2. small counts are never Gamma conditioning, nor a warning ---------------
    Cx = (.n, 1, 1, 1 \ .n, 2, 2, 3)
    out = dwlsmv_diagadvice(G, Cx, B, ind, slab, pname)
    L   = td_lines(out)
    if (td_cooccur(L, "count", "gamma") | td_cooccur(L, "count", "condition") |
        td_cooccur(L, "count", "singular")) {
        FAIL = FAIL \ "2: a category-count line mentions Gamma, conditioning or singularity"
    }
    if (rows(L) != rows(L0)) FAIL = FAIL \ "2: small counts added or removed a line (counts are descriptive, never a warning)"
    if (!td_find(L, "y1, lowest)")) FAIL = FAIL \ "2: the boundary count's location is not named"

    // ---- 3. every flag has its note ------------------------------------------------
    Gx = G; Gx[2, .] = (1, -1e-10, 0, 0)
    out = dwlsmv_diagadvice(Gx, C, B, ind, slab, pname)
    if (!(strpos(out, "not" + char(10) + "    positive semidefinite") & strpos(out, "internal numerical check failed") &
          strpos(out, "Please report"))) FAIL = FAIL \ "3: the PSD note is not an internal-check report"
    Gx = G; Gx[3, .] = (1, 12, 16, 0)
    if (!strpos(dwlsmv_diagadvice(Gx, C, B, ind, slab, pname), "Too few distinct response patterns: 12 patterns for 16")) FAIL = FAIL \ "3: no pattern note"
    Gx = G; Gx[4, .] = (1, 1, 1, 4)
    out = dwlsmv_diagadvice(Gx, C, B, ind, slab, pname)
    if (!(strpos(out, slab[1] + " and " + slab[4] + " move together") & strpos(out, "exactly singular"))) FAIL = FAIL \ "3: aliasing note does not name y1|t1 and y2|t1"
    Gx = G; Gx[5, .] = (1, 0.58, 11, 0)
    out = dwlsmv_diagadvice(Gx, C, B, ind, slab, pname)
    if (!(strpos(out, "The correlation " + slab[11]) & strpos(out, "estimating equation"))) FAIL = FAIL \ "3: centring note does not name its statistic"
    Bx = B; Bx[1, .] = (1, 20, 16, 0)
    if (!strpos(dwlsmv_diagadvice(G, C, Bx, ind, slab, pname), "More free parameters")) FAIL = FAIL \ "3: no npar > D note"
    Bx = B; Bx[2, .] = (1, 1, np, 0)
    out = dwlsmv_diagadvice(G, C, Bx, ind, slab, pname)
    if (!(strpos(out, pname[np] + " has no effect") & strpos(out, "derivative matrix"))) FAIL = FAIL \ "3: zero-column note does not name its parameter"
    Bx = B; Bx[3, .] = (1, np - 1, np, 0)
    out = dwlsmv_diagadvice(G, C, Bx, ind, slab, pname)
    if (!(strpos(out, "insufficient to determine all parameters") &
          strpos(out, "rank " + strofreal(np - 1) + " of " + strofreal(np)))) FAIL = FAIL \ "3: no numerical-rank note"
    // a flagged fit still prints no ratio
    if (strpos(out, "singular value") | strpos(out, "ratio")) FAIL = FAIL \ "3: a ratio printed with a flag note"

    // ---- 3b. added after the 2026-09-26 review ----------------------------------------
    // count values (M33), with distinct numbers so the rows cannot swap unseen
    out = dwlsmv_diagadvice(G, (.n, 7, 1, 1 \ .n, 11, 2, 3), B, ind, slab, pname)
    if (!(strpos(out, "category: 7 (") & strpos(out, "category:            11 ("))) FAIL = FAIL \ "3b: the count lines do not carry C's values"
    // a constant column names ITS statistic (M31)
    Gx = G; Gx[4, .] = (1, 1, 5, 0)
    out = dwlsmv_diagadvice(Gx, C, B, ind, slab, pname)
    if (!strpos(out, slab[5] + " is constant")) FAIL = FAIL \ "3b: the constant-column note does not name its statistic"
    // "evaluated, nothing to report" is not "not evaluated" (D3)
    out = dwlsmv_diagadvice(G, (.n, 706, 1, 1 \ ., ., 0, 0), B, ind, slab, pname)
    if (!strpos(out, "middle category:            not evaluated")) FAIL = FAIL \ "3b: an unevaluated interior row reads as all binary"
    out = dwlsmv_diagadvice(G, (.n, 706, 1, 1 \ .n, ., 0, 0), B, ind, slab, pname)
    if (!strpos(out, "none (every item is binary)")) FAIL = FAIL \ "3b: an all-binary interior row is not said so"

    // ---- 4. asymmetry stops the Gamma checks, and is an internal check -------------
    Gx = (1, 3e-14, 1, 2 \ J(6, 2, .), J(6, 2, 0))
    out = dwlsmv_diagadvice(Gx, C, B, ind, slab, pname)
    if (!(strpos(out, "not" + char(10) + "    symmetric") & strpos(out, "not run") &
          strpos(out, "internal numerical check failed") & strpos(out, "Please report"))) {
        FAIL = FAIL \ "4: asymmetry note is not an internal-check report saying the rest was not run"
    }
    // ---- 5. hygiene ----------------------------------------------------------------
    if (strpos(out0, "\" + "n")) FAIL = FAIL \ "5: a literal backslash-n in the advice"
    if (substr(out0, strlen(out0), 1) != char(10)) FAIL = FAIL \ "5: the advice does not end in char(10)"

    printf("{txt}test_diagadvice: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }
    return(rows(FAIL))
}

// ---- 6's contract violations ----------------------------------------------------------
void dwlsmv_t_da_bad(real scalar which)
{
    real matrix G, C, B
    string colvector ind, slab, pname

    G = (0, 0, 0, 0 \ 0, 0.1, 0, 0 \ 0, 139, 16, 0 \ 0, 0, 0, 0 \ 0, 1e-15, 1, 0 \
         .n, 0.1, 16, 0 \ .n, 0.1, 16, 0)
    C = (.n, 706, 1, 1 \ .n, 1500, 2, 3)
    B = (0, 14, 16, 0 \ 0, 0, 0, 0 \ 0, 14, 14, 0 \ .n, 0.12, 14, 0)
    ind = ("y1" \ "y2" \ "y3" \ "y4"); slab = J(16, 1, "s"); pname = J(14, 1, "p")
    if (which == 1) G = G[1..6, .]
    if (which == 2) C = C[., 1..3]
    if (which == 3) B = B[1..3, .]
    if (which == 4) C[1, 3] = 5
    if (which == 5) G[4, .] = (1, 1, 1, 17)
    if (which == 6) B[2, .] = (1, 1, 20, 0)
    if (which == 7) C[1, 4] = 0
    if (which == 8) C[2, 4] = 1
    (void) dwlsmv_diagadvice(G, C, B, ind, slab, pname)
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_diagadvice(YD)))

local nerr 0
forvalues b = 1/8 {
    capture mata: dwlsmv_t_da_bad(`b')
    if _rc != 3498 {
        display as error "test_diagadvice: contract violation #`b' gave rc " _rc ", not 3498"
        local nerr = `nerr' + 1
    }
}
capture mata: dwlsmv_t_da_bad(0)
if _rc {
    display as error "test_diagadvice: the valid base of check 6 was rejected, rc " _rc "; its violations prove nothing"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_diagadvice: `nfail' numeric failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_diagadvice: flags become notes in plain words, no ratio is printed, counts are always shown and never Gamma"
exit 0
