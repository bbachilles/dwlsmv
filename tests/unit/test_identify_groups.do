*! Module 3 under multiple groups: dwls_identify() against spec 4.8 G0-G7 and O1-O9.
*!
*! Drafted 2026-10-10 by a read-only agent from 003819a (cloud task M), with no
*! Stata available. It cannot load until the building commit lands, because it is
*! written against code that does not exist yet:
*!   (1) O14: dwls_identify(string, struct dwlsmv_statmap vector), dwls_delta(theta,
*!       par, struct dwlsmv_statmap vector). A single group is a 1-vector.
*!   (2) a group key on dwlsmv_parmap. This file assumes  real colvector grp  (one
*!       entry per row, 1..G, tied members keeping their own group). Spec G6: "the
*!       struct field lands with the code". The name is used in ONE place, tig_grp().
*!   (3) Mata syntax never run here, all confined to tig_mapsK(): the struct-vector
*!       constructor dwlsmv_statmap(1, G), element assignment m[g] = ..., and passing
*!       maps[g] to a function.
*!
*! WHAT MUST BE RUN BY THE MAIN SESSION, in this order:
*!   a. tests/unit/test_identify.do after editing ONLY its dwls_identify() call
*!      sites to 1-vectors (the expectations stay unedited; this is G0's first
*!      requirement).
*!   b. this file. Read the log of every rejection: the judge below checks the
*!      _error() text and tokens of the errprintf() text, not their wording.
*!   c. the mutants in the list below, one at a time.
*!
*! ADAPTED 2026-10-10 by the main session to the decisions of that day (each change
*! marked in the code): calls go through tig_ident(), which passes synthetic Stage 1
*! results to a groups text; expected tau starts are those sample values (tied: group
*! 1's), so t2 also tests O5(a); n13, n17a/b, n22a/b assert one message (precedence);
*! n29a-c test the Stage 1 argument; eq1-eq3 compare every field but tau starts, and
*! assert those equal group 1's sample thresholds. M27 should now be killed by n17.
*!
*! HOW A REJECTION IS JUDGED. Each case runs inside a named log; the log text is read
*! back. A case passes when _rc == 3498 AND the log holds "dwls_identify: " + one of
*! the expected _error() texts AND, after that _error() text is removed from the log,
*! every detail token still appears (so a token that only the _error() text carries
*! cannot satisfy it). msg "?" means rc 3498 only; msg "" means the text must be
*! accepted (a control). msg "a|b" means either text passes (the spec leaves the
*! order of two checks open; see the choices list).
*!
*! CASE IDS (the mutant list cites them)
*!   g0a-g0f  a single-group text carrying group / tie / nu / omega / alpha / tau
*!   g0ok     the base text, accepted
*!   eq1-eq3  G0's second requirement: T vs a groups-2 text of T and T + psi F F fixed 1
*!   n1a-n1g  'groups' malformed        n2a-n2b  groups twice
*!   n3a-n3b  group outside 1..G        n4       a section twice
*!   n5a-n5b  numbering gap / missing   n6a-n6b  directive outside a section
*!   n7a-n7b  factors differ            n8a-n8b  indicators differ (order, name)
*!   n9a-n9b  tie with < 2 operands     n10a-b   unknown / incomplete operand
*!   n11a-d   O4 whitelist (psi, nu, omega, alpha)
*!   n12a-d   O3, one key position each (lambda F-same y-diff; F-diff y-same;
*!            tau y-same k-diff; y-diff k-same)
*!   n13      lambda joined to tau      n14a-b   tau with different K
*!   n15a-c   absent element            n16a-b   fixed element (1st / 2nd operand)
*!   n17a-b   two operands, one group   n18a-b   duplicate / overlapping ties
*!   n19a-b   tied lambda, differing starts (2 groups; 3 groups, third differs)
*!   n20*     O9: nu / omega / alpha in group 1 x {ref, else, free}
*!   n21a-d   omega must be positive    n22a-b   metric outside group 1
*!   n23a-b   marker metric (O2)        n24a-s   structural path (O2, lifted 2026-10-10, widened twice 2026-10-11)
*!   n25a-d   map count != G            n26a-b   later-group factor, no psi
*!   n27a-c   per-group existing rules  n28      standalone tau under groups
*!   k1-k8    controls (must be ACCEPTED)
*!   t1  design E configural: nall 34, npar 32, every row's (ptype, group, row, col,
*!       fidx, free, value) against a hand-written table, in table order   [O6]
*!   t2  design E intercepts: nall 40, npar 23, ties share fidx, no group-1 nu/omega/
*!       alpha rows and no nu rows [O1], block order tau, lambda, psi, omega, alpha
*!   t3  configural Delta: threshold rows x threshold columns is EXACTLY I(22), and
*!       threshold rows x structural columns is exactly 0               [O6, G6]
*!   t4  three groups, tie lambda 2 F y1 = lambda 3 F y1: the shared fidx is the
*!       group-2 member's position; nall 51, npar 47                    [G6]
*!   t5-t8  groups line last / between sections / sections reversed / directives
*!       shuffled inside sections: each gives t1's table field for field  [G1]
*!
*! MUTANTS (each must turn this file red; "killed by" = the first case that fails)
*!   M1  a missing groups line read as groups 1 (so 'group' opens a section)  g0a
*!   M2  nu, omega, alpha or tie accepted without groups                      g0b-g0e
*!   M3  standalone tau accepted, with or without groups                      g0f, n28
*!   M4  groups 1 accepted                                                    n1c
*!   M5  groups 0, 2.5, abc or -2 accepted (an integer-only test dropped)     n1d-n1g
*!   M6  second groups line ignored (last wins, or first wins)                n2a, n2b
*!   M7  group number outside 1..G accepted                                   n3a, n3b
*!   M8  duplicate section merged or last-wins                                n4
*!   M9  gap / missing section not detected                                   n5a, n5b
*!   M10 a directive after a tie line attributed to the last section          n6b
*!   M11 a directive before the first group accepted                          n6a
*!   M12 factors not compared with group 1 (name, count)                      n7a, n7b
*!   M13 indicators compared as a set, or only by count                       n8a, n8b
*!   M14 tie with 0 or 1 operands accepted                                    n9a, n9b
*!   M15 unknown tie operand type skipped as a comment                        n10a
*!   M16 incomplete trailing operand ignored                                  n10b
*!   M17 O4 whitelist dropped; or applied only to generator-written ties      n11a-n11d
*!   M18 O3 compares only the y / the F / the k / the indicator position      n12a-n12d
*!         (a mutant that compares only one position survives three of four)
*!   M19 operand types not compared                                           n13
*!   M20 K-mismatch check removed                                             n14a, n14b
*!   M21 K-mismatch compares only the first threshold's existence             n14b
*!   M22 absent-element check removed                                         n15a-n15c
*!   M23 absent check on the first operand only                               n15b
*!   M24 fixed check removed                                                  n16a, n16b
*!   M25 fixed check on the first operand only                                n16b
*!   M26 fixed check on the last operand only                                 n16a
*!   M27 distinct-group check removed                                         n17a, n17b (see choices item 5: may be equivalent)
*!   M28 duplicate check across ties removed                                  n18a, n18b
*!   M29 duplicate check compares identical lines only                        n18b
*!   M30 differing tied-lambda starts accepted                                n19a
*!   M31 lowest group's start wins (no refusal)                               n19a
*!   M32 start comparison on the first two operands only                      n19b
*!   M33 start comparison applied to tau as well (tau starts are all 0, so it  k6, t2
*!         must not fire; a data-dependent tau start would)
*!   M34 O9 accepts fixed at the reference value (nu 0, omega 1, alpha 0)     n20nuR, n20omR, n20alR
*!   M35 O9 refuses only 'free'                                               n20*E, n20*R
*!   M36 O9 for omega only (nu and alpha accepted)                            n20nu*, n20al*
*!   M37 O9 reports the omega-positivity text for omega fixed at 1 or elsewhere   n20omR/E (text assertion)
*!   M38 O9 message omits the indicator/factor, the status, or "reference"    n20* details
*!   M39 omega fixed 0 / start 0 accepted (positivity as > -eps, or only < 0) n21a, n21c
*!   M40 positivity checked for fixed only                                    n21c, n21d
*!   M41 positivity checked for free only                                     n21a, n21b
*!   M42 omega positivity refuses fixed 1 (the reference value in group >= 2) k5
*!   M43 metric accepted in group >= 2                                        n22a, n22b
*!   M44 marker refused only with a tie / nu / omega / alpha present          n23a
*!   M45 RETIRED 2026-10-11: it accepted a one-group beta placed before the lambda
*!         lines at G = 2, which is now the correct behaviour (n24c, n24d accepted).
*!   M46 the free-in-every-group rule (G > 2) checks group 1 only               n24p
*!   M47 beta refused only when the text has a tie (not configural)           n24p, n24q, n24r
*!   M48 a fixed nonzero beta's value not compared across groups (fixed in
*!         every group at any nonzero values accepted)                         n24m
*!   M67 fixed nonzero beta still refused (the 2026-10-11 widening not applied) n24i, n24o
*!   M68 a group where the beta is free counted as agreeing                   n24j
*!   (M48, M67, M68 executed 2026-10-11: each turned the test red.)
*!   M64 "free in every group" read as "free in some group" (G > 2)          n24p, n24r
*!   M65 a beta free in every group still refused (the lift not applied)      n24g, n24k, n24s
*!   M66 beta fixed at 0 treated as nonzero (absent path refused)             n24l, n24h
*!   M69 the G > 2 restriction dropped (one-group beta accepted at any G)     n24p, n24q, n24r
*!   M70 one-group beta still refused at G = 2 (the 2026-10-11 widening not   n24a-e, n24h
*!         applied)
*!   (M69, M70 executed 2026-10-11: each turned the test red. M46, M47 and M64 were
*!   re-targeted to the G = 3 cases and were NOT re-executed in that form.)
*!   M49 map count not compared with G                                        n25a-n25d
*!   M50 later-group factor without psi accepted (defaults to 1)              n26a, n26b
*!   M51 the per-group factor/metric/duplicate rules applied to group 1 only  n27a-n27c
*!   M52 theta ordered (row, col, group) in tau / lambda                      t1, t3, t4
*!   M53 theta fully group-major (all of g1's blocks, then g2's)              t1, t3
*!   M54 nu, omega, alpha placed before psi, or ordered (row, group)          t2
*!   M55 a tie's fidx taken from the highest-group member                     t4
*!   M56 tied members given separate fidx (tie ignored)                       t2, t4
*!   M57 a tie adds a row, or changes a member's free flag                    t2 (nall 40)
*!   M58 a reference-valued nu/omega/alpha written as a fixed row (O1)        t2 (nall 40)
*!   M59 group 1 loses its group key (all rows grp 1)                         t1, eq1-eq3
*!   M60 order of sections or of lines inside a section changes the table     t5-t8
*!   M61 a groups line anywhere but first refused                             t5, t6
*!   M62 the Delta rank check run per group, not once on the stacked table    needs a case where only the stacked
*!         table is deficient; NOT covered here (choices item 14)
*!   M63 single-group behaviour changed (a missing groups line altering the   eq1-eq3, g0ok, and test_identify.do
*!         single-group theta, starts or order)
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_identify_groups: `f' failed to load"
        exit 9
    }
}

mata:

// ---------------------------------------------------------------------------
// The one place the group key is named.
// ---------------------------------------------------------------------------
real colvector tig_grp(struct dwlsmv_parmap scalar p)
{
    return(p.grp)
}

// ---------------------------------------------------------------------------
// String plumbing
// ---------------------------------------------------------------------------
string scalar tig_join(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

string rowvector tig_none()
{
    return(J(1, 0, ""))
}

string rowvector tig_free(real scalar n)
{
    return(J(1, n, "free"))
}

string rowvector tig_IE()
{
    return(("y1", "y2", "y3", "y4", "y5"))
}

// n = 5 status vector with entry j replaced
string rowvector tig_with(real scalar j, string scalar s)
{
    string rowvector v

    v    = tig_free(5)
    v[j] = s
    return(v)
}

string rowvector tig_split(string scalar s, string scalar sep)
{
    string rowvector out
    string scalar t
    real scalar k

    t   = s
    out = J(1, 0, "")
    while ((k = strpos(t, sep)) > 0) {
        out = out, substr(t, 1, k - 1)
        t   = substr(t, k + strlen(sep), .)
    }
    return((out, t))
}

string scalar tig_slurp(string scalar path)
{
    real scalar fh
    string scalar line, out

    fh  = fopen(path, "r")
    out = ""
    while ((line = fget(fh)) != J(0, 0, "")) {
        out = out + line + char(10)
    }
    fclose(fh)
    return(out)
}

// ---------------------------------------------------------------------------
// Section builders. st[j] is the whole status text after the indicator name
// ("free", "free 0.6", "fixed 0.8") or "OMIT" for no lambda line at all.
// top = "" means the default: the metric line in group 1, psi F F fixed 1 after.
// ---------------------------------------------------------------------------
string rowvector tig_sec(real scalar g, string scalar fn, string rowvector inds,
                         string rowvector st, string scalar top, string rowvector extra)
{
    string rowvector out
    real scalar j

    out = ("group " + strofreal(g), "factor " + fn, "ind " + invtokens(inds))
    for (j = 1; j <= cols(inds); j++) {
        if (st[j] != "OMIT") {
            out = out, ("lambda " + fn + " " + inds[j] + " " + st[j])
        }
    }
    if (top != "") {
        out = out, top
    }
    else if (g == 1) {
        out = out, ("metric unitvariance " + fn)
    }
    else {
        out = out, ("psi " + fn + " " + fn + " fixed 1")
    }
    return((out, extra))
}

// design E section: one factor F, y1-y5
string rowvector tig_s(real scalar g, string rowvector extra)
{
    return(tig_sec(g, "F", tig_IE(), tig_free(5), "", extra))
}

string rowvector tig_sx(real scalar g, string rowvector st, string scalar top,
                        string rowvector extra)
{
    return(tig_sec(g, "F", tig_IE(), st, top, extra))
}

// two factors, y1-y6 (F1: y1-y3, F2: y4-y6), no structural path
string rowvector tig_sec2(real scalar g, string rowvector extra)
{
    string rowvector out

    out = ("group " + strofreal(g), "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
           "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
           "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free")
    if (g == 1) {
        out = out, ("metric unitvariance F1", "metric unitvariance F2")
    }
    else {
        out = out, ("psi F1 F1 fixed 1", "psi F2 F2 fixed 1")
    }
    return((out, extra))
}

// design E's tie lines, spec G8: tau ties first (y1-y3 k = 1..3, y4 1, y5 1), then lambda
string rowvector tig_ties_int()
{
    string rowvector out
    real rowvector C
    real scalar j, k

    C   = (3, 3, 3, 1, 1)
    out = J(1, 0, "")
    for (j = 1; j <= 5; j++) {
        for (k = 1; k <= C[j]; k++) {
            out = out, ("tie tau 1 y" + strofreal(j) + " " + strofreal(k) +
                        " tau 2 y" + strofreal(j) + " " + strofreal(k))
        }
    }
    for (j = 1; j <= 5; j++) {
        out = out, ("tie lambda 1 F y" + strofreal(j) + " lambda 2 F y" + strofreal(j))
    }
    return(out)
}

string rowvector tig_Econf_lines()
{
    return(("groups 2", tig_s(1, tig_none()), tig_s(2, tig_none())))
}

string rowvector tig_Eint_lines()
{
    string rowvector ex

    ex = ("omega y1 free", "omega y2 free", "omega y3 free", "omega y4 free",
          "omega y5 free", "alpha F free")
    return(("groups 2", tig_s(1, tig_none()),
            tig_sx(2, tig_free(5), "psi F F free", ex), tig_ties_int()))
}

// ---------------------------------------------------------------------------
// Index maps. Codes: 1 design E x 2 groups; 2 design E x 3; 3 design E and a
// group 2 whose y1 has K = 3; 4 design E x 1; 5 case D x 1; 6 six K = 3 x 2.
// ---------------------------------------------------------------------------
struct dwlsmv_statmap vector tig_mapsK(real matrix Kc)
{
    struct dwlsmv_statmap vector m
    real scalar g

    m = dwlsmv_statmap(1, cols(Kc))
    for (g = 1; g <= cols(Kc); g++) {
        m[g] = dwlsmv_statmap_build(Kc[., g])
    }
    return(m)
}

struct dwlsmv_statmap vector tig_maps(real scalar code)
{
    real colvector KE

    KE = (4 \ 4 \ 4 \ 2 \ 2)
    if (code == 1) {
        return(tig_mapsK((KE, KE)))
    }
    if (code == 2) {
        return(tig_mapsK((KE, KE, KE)))
    }
    if (code == 3) {
        return(tig_mapsK((KE, (3 \ 4 \ 4 \ 2 \ 2))))
    }
    if (code == 4) {
        return(tig_mapsK(KE))
    }
    if (code == 5) {
        return(tig_mapsK((4 \ 4 \ 3 \ 3)))
    }
    if (code == 7) {
        return(tig_mapsK((J(6, 1, 3), J(6, 1, 3), J(6, 1, 3))))
    }
    return(tig_mapsK((J(6, 1, 3), J(6, 1, 3))))
}

// ---------------------------------------------------------------------------
// Synthetic Stage 1 results (adapted 2026-10-10, main session). Under groups
// dwls_identify() takes the Stage 1 results and starts thresholds at the sample
// values (spec 4.8 O14 C3 note, O5). Thresholds here are distinct, nonzero and
// differ by group: qnorm(k/K_j) + 0.02 j + 0.05 (g - 1), so a test can tell group
// 1's value from group 2's. Correlations are 0.3; nothing else is read.
// ---------------------------------------------------------------------------
real scalar tig_tau(real scalar g, real scalar j, real scalar k, real scalar Kj)
{
    return(invnormal(k / Kj) + 0.02 * j + 0.05 * (g - 1))
}

struct dwlsmv_stage1 vector tig_st(struct dwlsmv_statmap vector m)
{
    struct dwlsmv_stage1 vector st
    real scalar g, i

    st = dwlsmv_stage1(1, cols(m))
    for (g = 1; g <= cols(m); g++) {
        st[g].map = m[g]
        st[g].s   = J(m[g].D, 1, 0.3)
        for (i = 1; i <= m[g].T; i++) {
            st[g].s[i] = tig_tau(g, m[g].var_i[i], m[g].thr_k[i], m[g].K[m[g].var_i[i]])
        }
    }
    return(st)
}

// a text with a 'groups' line gets the Stage 1 results; one without does not
struct dwlsmv_parmap scalar tig_ident(string scalar spec, struct dwlsmv_statmap vector m)
{
    if (strpos(char(10) + spec, char(10) + "groups") > 0) {
        return(dwls_identify(spec, m, tig_st(m)))
    }
    return(dwls_identify(spec, m))
}

// tau start values of design E's threshold pairs tp in group g
real colvector tig_tv(real scalar g, real matrix tp)
{
    real colvector KE, v
    real scalar i

    KE = (4 \ 4 \ 4 \ 2 \ 2)
    v  = J(rows(tp), 1, .)
    for (i = 1; i <= rows(tp); i++) {
        v[i] = tig_tau(g, tp[i, 1], tp[i, 2], KE[tp[i, 1]])
    }
    return(v)
}

// code < 100: tig_ident. 100 + c: a groups text with NO Stage 1. 200 + c: one
// Stage 1 result only. 300 + c: Stage 1 built on code 1's maps, passed with c's.
void tig_try(string scalar spec, real scalar code)
{
    struct dwlsmv_parmap scalar p
    struct dwlsmv_statmap vector m

    if (code >= 300) {
        m = tig_maps(code - 300)
        p = dwls_identify(spec, m, tig_st(tig_maps(1)))
        return
    }
    if (code >= 200) {
        m = tig_maps(code - 200)
        p = dwls_identify(spec, m, tig_st(m)[1])
        return
    }
    if (code >= 100) {
        m = tig_maps(code - 100)
        p = dwls_identify(spec, m)
        return
    }
    m = tig_maps(code)
    p = tig_ident(spec, m)
}

// ---------------------------------------------------------------------------
// The judge. row = (id, why, spec, code, msg, det).
// ---------------------------------------------------------------------------
string scalar tig_judge(string rowvector row, real scalar rc, string scalar txt)
{
    string rowvector alts, dets
    string scalar rest
    real scalar a, d, hit

    if (row[5] == "") {
        if (rc == 0) {
            return("")
        }
        return("a control was refused, rc " + strofreal(rc))
    }
    if (rc != 3498) {
        return("rc " + strofreal(rc) + ", wanted 3498")
    }
    if (row[5] == "?") {
        return("")
    }
    alts = tig_split(row[5], "|")
    hit  = 0
    for (a = 1; a <= cols(alts); a++) {
        if (hit == 0) {
            if (strpos(txt, "dwls_identify: " + alts[a]) > 0) {
                hit = a
            }
        }
    }
    if (hit == 0) {
        return("none of the expected _error() texts is in the log: " + row[5])
    }
    rest = subinstr(txt, "dwls_identify: " + alts[hit], "", .)
    if (row[6] != "") {
        dets = tig_split(row[6], ";;")
        for (d = 1; d <= cols(dets); d++) {
            if (strpos(rest, dets[d]) == 0) {
                return("the errprintf() text lacks: " + dets[d])
            }
        }
    }
    return("")
}

void tig_add(string matrix T, string scalar id, string scalar why, string scalar spec,
             real scalar code, string scalar msg, string scalar det)
{
    T = T \ (id, why, spec, strofreal(code), msg, det)
}

// ---------------------------------------------------------------------------
// The case table.
// ---------------------------------------------------------------------------
void tig_cases(string matrix T)
{
    string rowvector none, s1, s2, s3, w, tie1
    string scalar D, M

    none = tig_none()
    s1   = tig_s(1, none)
    s2   = tig_s(2, none)
    s3   = tig_s(3, none)
    tie1 = ("tie lambda 1 F y1 lambda 2 F y1")

    // ---- G0: a single-group text with any of the five tokens -----------------
    D = tig_join(("factor f", "ind y1 y2 y3 y4",
                  "lambda f y1 free 0.75", "lambda f y2 free 0.65",
                  "lambda f y3 free 0.55", "lambda f y4 free 0.45",
                  "metric unitvariance f"))
    M = "unrecognized directive"
    tig_add(T, "g0ok", "case D single factor, no groups line (control)", D, 5, "", "")
    tig_add(T, "g0a", "'group 1' opens a section in a text with no groups line",
            "group 1" + char(10) + D, 5, M, "")
    tig_add(T, "g0b", "'tie' without groups",
            D + "tie lambda 1 f y1 lambda 2 f y1" + char(10), 5, M, "")
    tig_add(T, "g0c", "'nu' without groups", D + "nu y1 free" + char(10), 5, M, "")
    tig_add(T, "g0d", "'omega' without groups", D + "omega y1 free" + char(10), 5, M, "")
    tig_add(T, "g0e", "'alpha' without groups", D + "alpha f free" + char(10), 5, M, "")
    tig_add(T, "g0f", "'tau' stays rejected without groups",
            D + "tau y1 1 free" + char(10), 5, M, "")

    // ---- G7: groups line --------------------------------------------------------
    w = ("'groups' needs one integer, 2 or more")
    tig_add(T, "n1a", "groups with no integer",
            tig_join(("groups", s1, s2)), 1, w, "")
    tig_add(T, "n1b", "groups with two integers",
            tig_join(("groups 2 3", s1, s2)), 1, w, "")
    tig_add(T, "n1c", "groups 1", tig_join(("groups 1", s1)), 4, w, "")
    tig_add(T, "n1d", "groups 0", tig_join(("groups 0", s1)), 4, w, "")
    tig_add(T, "n1e", "groups 2.5", tig_join(("groups 2.5", s1, s2)), 1, w, "")
    tig_add(T, "n1f", "groups abc", tig_join(("groups abc", s1, s2)), 1, w, "")
    tig_add(T, "n1g", "groups -2", tig_join(("groups -2", s1, s2)), 1, w, "")
    tig_add(T, "n2a", "groups twice, same value",
            tig_join(("groups 2", "groups 2", s1, s2)), 1,
            "more than one groups directive", "")
    tig_add(T, "n2b", "groups twice, different values",
            tig_join(("groups 2", s1, s2, "groups 3")), 1,
            "more than one groups directive", "")

    // ---- G7: sections -----------------------------------------------------------
    tig_add(T, "n3a", "group 3 in a groups-2 text (sections 1, 2, 3)",
            tig_join(("groups 2", s1, s2, s3)), 1, "group number outside 1..G", "")
    tig_add(T, "n3b", "group 0 in a groups-2 text (sections 0, 1, 2)",
            tig_join(("groups 2", tig_s(0, none), s1, s2)), 1,
            "group number outside 1..G", "")
    tig_add(T, "n4", "group 1 section twice",
            tig_join(("groups 2", s1, s1, s2)), 1, "a group section appears twice", "")
    tig_add(T, "n5a", "groups 3 with sections 1 and 3",
            tig_join(("groups 3", s1, s3)), 2, "group sections are not numbered 1..G", "")
    tig_add(T, "n5b", "groups 2 with section 1 only",
            tig_join(("groups 2", s1)), 1, "group sections are not numbered 1..G", "")
    tig_add(T, "n6a", "a factor line before the first group line",
            tig_join(("groups 2", "factor F", s1, s2)), 1,
            "a directive lies outside any group section", "")
    tig_add(T, "n6b", "an omega line after the tie line (a tie ends the section)",
            tig_join(("groups 2", s1, s2, tie1, "omega y1 free")), 1,
            "a directive lies outside any group section", "")
    tig_add(T, "n7a", "group 2 declares factor G, group 1 declares F",
            tig_join(("groups 2", s1, tig_sec(2, "G", tig_IE(), tig_free(5), "", none))), 1,
            "a group declares different factors from group 1", "")
    tig_add(T, "n7b", "group 1 declares F1 F2, group 2 only F1",
            tig_join(("groups 2", tig_sec2(1, none),
                      tig_sec(2, "F1", ("y1", "y2", "y3", "y4", "y5", "y6"),
                              tig_free(6), "", none))), 6,
            "a group declares different factors from group 1", "")
    w     = s2
    w[3]  = "ind y1 y2 y3 y5 y4"
    tig_add(T, "n8a", "group 2 lists the same indicators in another order",
            tig_join(("groups 2", s1, w)), 1,
            "a group declares different indicators from group 1", "")
    tig_add(T, "n8b", "group 2 declares z5 where group 1 declares y5",
            tig_join(("groups 2", s1,
                      tig_sec(2, "F", ("y1", "y2", "y3", "y4", "z5"), tig_free(5), "", none))), 1,
            "a group declares different indicators from group 1", "")

    // ---- G7: tie operands -------------------------------------------------------
    tig_add(T, "n9a", "tie with no operand",
            tig_join(("groups 2", s1, s2, "tie")), 1, "a tie needs at least two elements", "")
    tig_add(T, "n9b", "tie with one operand",
            tig_join(("groups 2", s1, s2, "tie lambda 1 F y1")), 1,
            "a tie needs at least two elements", "")
    tig_add(T, "n10a", "an unknown element type",
            tig_join(("groups 2", s1, s2, "tie foo 1 F y1 foo 2 F y1")), 1,
            "unrecognized element type in a tie", "")
    tig_add(T, "n10b", "the second operand is cut short",
            tig_join(("groups 2", s1, s2, "tie lambda 1 F y1 lambda 2 F")), 1,
            "an element in a tie is incomplete", "")

    w = ("a tie may join only tau or lambda elements")
    tig_add(T, "n11a", "O4: psi operands, both rows present and free",
            tig_join(("groups 2", tig_sec2(1, ("psi F2 F1 free")), tig_sec2(2, ("psi F2 F1 free")),
                      "tie psi 1 F2 F1 psi 2 F2 F1")), 6, w, "psi")
    tig_add(T, "n11b", "O4: nu operands (groups 2 and 3)",
            tig_join(("groups 3", s1, tig_s(2, ("nu y1 free")), tig_s(3, ("nu y1 free")),
                      "tie nu 2 y1 nu 3 y1")), 2, w, "nu")
    tig_add(T, "n11c", "O4: omega operands (groups 2 and 3)",
            tig_join(("groups 3", s1, tig_s(2, ("omega y1 free")), tig_s(3, ("omega y1 free")),
                      "tie omega 2 y1 omega 3 y1")), 2, w, "omega")
    tig_add(T, "n11d", "O4: alpha operands (groups 2 and 3)",
            tig_join(("groups 3", s1, tig_s(2, ("alpha F free")), tig_s(3, ("alpha F free")),
                      "tie alpha 2 F alpha 3 F")), 2, w, "alpha")

    w = ("a tie joins different elements across groups")
    tig_add(T, "n12a", "O3 lambda: same F, different y",
            tig_join(("groups 2", s1, s2, "tie lambda 1 F y1 lambda 2 F y2")), 1, w, "y1;;y2")
    tig_add(T, "n12b", "O3 lambda: different F, same y (y1 cross-loads in both groups)",
            tig_join(("groups 2", tig_sec2(1, ("lambda F2 y1 free")),
                      tig_sec2(2, ("lambda F2 y1 free")),
                      "tie lambda 1 F1 y1 lambda 2 F2 y1")), 6, w, "F1;;F2")
    tig_add(T, "n12c", "O3 tau: same y, different k",
            tig_join(("groups 2", s1, s2, "tie tau 1 y1 1 tau 2 y1 2")), 1, w, "y1")
    tig_add(T, "n12d", "O3 tau: different y, same k",
            tig_join(("groups 2", s1, s2, "tie tau 1 y1 1 tau 2 y2 1")), 1, w, "y1;;y2")
    tig_add(T, "n13", "a lambda joined to a tau (precedence: types before elements)",
            tig_join(("groups 2", s1, s2, "tie lambda 1 F y1 tau 2 y1 1")), 1,
            "a tie joins elements of different types", "")
    w = ("a tie joins thresholds of indicators with different K")
    tig_add(T, "n14a", "tau 1 of y1, K = 4 in group 1 and 3 in group 2",
            tig_join(("groups 2", s1, s2, "tie tau 1 y1 1 tau 2 y1 1")), 3, w, "y1")
    tig_add(T, "n14b", "tau 2 of y1 (exists in both groups), same K mismatch",
            tig_join(("groups 2", s1, s2, "tie tau 1 y1 2 tau 2 y1 2")), 3, w, "y1")
    w = ("a tie names an element that is not in the table")
    tig_add(T, "n15a", "tau k = 5 of y1 (y1 has 3 thresholds in both groups)",
            tig_join(("groups 2", s1, s2, "tie tau 1 y1 5 tau 2 y1 5")), 1, w, "")
    tig_add(T, "n15b", "lambda y5 present in group 1, absent in group 2",
            tig_join(("groups 2", s1, tig_sx(2, tig_with(5, "OMIT"), "", none),
                      "tie lambda 1 F y5 lambda 2 F y5")), 1, w, "")
    tig_add(T, "n15c", "lambda y5 absent in both groups",
            tig_join(("groups 2", tig_sx(1, tig_with(5, "OMIT"), "", none),
                      tig_sx(2, tig_with(5, "OMIT"), "", none),
                      "tie lambda 1 F y5 lambda 2 F y5")), 1, w, "")
    w = ("a tie names a fixed element")
    tig_add(T, "n16a", "the first operand (group 1 y1) is fixed",
            tig_join(("groups 2", tig_sx(1, tig_with(1, "fixed 0.8"), "", none), s2, tie1)),
            1, w, "")
    tig_add(T, "n16b", "the second operand (group 2 y1) is fixed",
            tig_join(("groups 2", s1, tig_sx(2, tig_with(1, "fixed 0.8"), "", none), tie1)),
            1, w, "")
    tig_add(T, "n17a", "one element named twice in a tie (precedence: two-in-one-group first)",
            tig_join(("groups 2", s1, s2, "tie lambda 1 F y1 lambda 1 F y1")), 1,
            "a tie names two elements of one group", "")
    tig_add(T, "n17b", "three groups, group 1 named twice (precedence: two-in-one-group first)",
            tig_join(("groups 3", s1, s2, s3,
                      "tie lambda 1 F y1 lambda 2 F y1 lambda 1 F y1")), 2,
            "a tie names two elements of one group", "")
    w = ("an element appears in a tie more than once")
    tig_add(T, "n18a", "the same tie line twice",
            tig_join(("groups 2", s1, s2, tie1, tie1)), 1, w, "")
    tig_add(T, "n18b", "two ties overlapping at group 2's y1",
            tig_join(("groups 3", s1, s2, s3, tie1, "tie lambda 2 F y1 lambda 3 F y1")), 2, w, "")
    w = ("tied elements carry different start values")
    tig_add(T, "n19a", "tied lambda y1 starts 0.6 in group 1, 0.8 in group 2",
            tig_join(("groups 2", tig_sx(1, tig_with(1, "free 0.6"), "", none),
                      tig_sx(2, tig_with(1, "free 0.8"), "", none), tie1)), 1, w, "")
    tig_add(T, "n19b", "a three-way tie, starts 0.6, 0.6, 0.8",
            tig_join(("groups 3", tig_sx(1, tig_with(1, "free 0.6"), "", none),
                      tig_sx(2, tig_with(1, "free 0.6"), "", none),
                      tig_sx(3, tig_with(1, "free 0.8"), "", none),
                      "tie lambda 1 F y1 lambda 2 F y1 lambda 3 F y1")), 2, w, "")

    // ---- G7 / O9: nu, omega, alpha in group 1, all nine -------------------------
    M = "nu, omega and alpha are not allowed in group 1"
    tig_add(T, "n20nuR", "O9 nu fixed 0 (the reference value) in group 1",
            tig_join(("groups 2", tig_s(1, ("nu y2 fixed 0")), s2)), 1, M,
            "group 1;;nu;;y2;;fixed;;reference")
    tig_add(T, "n20nuE", "O9 nu fixed 0.3 in group 1",
            tig_join(("groups 2", tig_s(1, ("nu y2 fixed 0.3")), s2)), 1, M,
            "group 1;;nu;;y2;;fixed;;reference")
    tig_add(T, "n20nuF", "O9 nu free in group 1",
            tig_join(("groups 2", tig_s(1, ("nu y2 free")), s2)), 1, M,
            "group 1;;nu;;y2;;free;;reference")
    tig_add(T, "n20omR", "O9 omega fixed 1 (the reference value) in group 1",
            tig_join(("groups 2", tig_s(1, ("omega y2 fixed 1")), s2)), 1, M,
            "group 1;;omega;;y2;;fixed;;reference")
    tig_add(T, "n20omE", "O9 omega fixed 1.5 in group 1",
            tig_join(("groups 2", tig_s(1, ("omega y2 fixed 1.5")), s2)), 1, M,
            "group 1;;omega;;y2;;fixed;;reference")
    tig_add(T, "n20omF", "O9 omega free in group 1",
            tig_join(("groups 2", tig_s(1, ("omega y2 free")), s2)), 1, M,
            "group 1;;omega;;y2;;free;;reference")
    tig_add(T, "n20alR", "O9 alpha fixed 0 (the reference value) in group 1",
            tig_join(("groups 2", tig_s(1, ("alpha F fixed 0")), s2)), 1, M,
            "group 1;;alpha;;F;;fixed;;reference")
    tig_add(T, "n20alE", "O9 alpha fixed 0.2 in group 1",
            tig_join(("groups 2", tig_s(1, ("alpha F fixed 0.2")), s2)), 1, M,
            "group 1;;alpha;;F;;fixed;;reference")
    tig_add(T, "n20alF", "O9 alpha free in group 1",
            tig_join(("groups 2", tig_s(1, ("alpha F free")), s2)), 1, M,
            "group 1;;alpha;;F;;free;;reference")

    // ---- G7: omega positivity (a group >= 2 line; Module 3's own refusal) -------
    w = ("omega must be positive")
    tig_add(T, "n21a", "omega fixed 0 in group 2",
            tig_join(("groups 2", s1, tig_s(2, ("omega y1 fixed 0")))), 1, w, "")
    tig_add(T, "n21b", "omega fixed -0.5 in group 2",
            tig_join(("groups 2", s1, tig_s(2, ("omega y1 fixed -0.5")))), 1, w, "")
    tig_add(T, "n21c", "omega free with start 0 in group 2",
            tig_join(("groups 2", s1, tig_s(2, ("omega y1 free 0")))), 1, w, "")
    tig_add(T, "n21d", "omega free with start -1 in group 2",
            tig_join(("groups 2", s1, tig_s(2, ("omega y1 free -1")))), 1, w, "")

    // ---- G5: metric outside group 1; O2 -----------------------------------------
    tig_add(T, "n22a", "metric unitvariance in group 2, no psi row (precedence: metric first)",
            tig_join(("groups 2", s1, tig_sx(2, tig_free(5), "metric unitvariance F", none))), 1,
            "metric appears in group 1 only", "")
    tig_add(T, "n22b", "metric unitvariance in group 2 beside psi F F fixed 1 (precedence: metric first)",
            tig_join(("groups 2", s1, tig_s(2, ("metric unitvariance F")))), 1,
            "metric appears in group 1 only", "")
    w = ("marker metric is not supported with groups")
    tig_add(T, "n23a", "O2: marker metric in group 1, configural, no ties",
            tig_join(("groups 2",
                      tig_sx(1, tig_with(1, "OMIT"), "metric marker F y1", ("psi F F free")), s2)),
            1, w, "")
    tig_add(T, "n23b", "O2: marker metric in group 1 with a tie and an omega row",
            tig_join(("groups 2",
                      tig_sx(1, tig_with(1, "OMIT"), "metric marker F y1", ("psi F F free")),
                      tig_s(2, ("omega y1 free")), "tie lambda 1 F y2 lambda 2 F y2")), 1, w, "")

    // O2: free in every group admitted (2026-10-10); with exactly two groups, free
    // in one group and absent or fixed at 0 in the other admitted too (2026-10-11,
    // after mgb-absent-prereg passed); at G > 2 free in some groups only refused.
    // n24a-e, h were refusals before 2026-10-11 and are now ACCEPTED controls.
    s2 = tig_sec2(2, none)
    tig_add(T, "n24a", "O2 (G = 2): beta F2 F1 free in group 2 only, after the lambda lines: ACCEPTED",
            tig_join(("groups 2", tig_sec2(1, none), tig_sec2(2, ("beta F2 F1 free")))), 6, "", "")
    tig_add(T, "n24b", "O2 (G = 2): the same text, sections in reverse textual order: ACCEPTED",
            tig_join(("groups 2", tig_sec2(2, ("beta F2 F1 free")), tig_sec2(1, none))), 6, "", "")
    tig_add(T, "n24c", "O2 (G = 2): beta free in group 2 only, directly after the group 2 line: ACCEPTED",
            tig_join(("groups 2", tig_sec2(1, none), s2[1], "beta F2 F1 free",
                      s2[2..cols(s2)])), 6, "", "")
    tig_add(T, "n24d", "O2 (G = 2): beta free in group 2 only, before the lambda lines: ACCEPTED",
            tig_join(("groups 2", tig_sec2(1, none), s2[1..3], "beta F2 F1 free",
                      s2[4..cols(s2)])), 6, "", "")
    tig_add(T, "n24e", "O2 (G = 2): beta free in group 1 only: ACCEPTED",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 free")), tig_sec2(2, none))), 6, "", "")
    tig_add(T, "n24h", "O2 (G = 2): beta free in group 1, fixed at 0 in group 2: ACCEPTED",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 free")),
                      tig_sec2(2, ("beta F2 F1 fixed 0")))), 6, "", "")
    w = ("with more than two groups, a beta must be free in every group")
    tig_add(T, "n24p", "O2 (G = 3): beta free in groups 1 and 2, absent in group 3",
            tig_join(("groups 3", tig_sec2(1, ("beta F2 F1 free")), tig_sec2(2, ("beta F2 F1 free")),
                      tig_sec2(3, none))), 7, w, "beta[F2,F1];;group 1;;group 3;;with 3 groups")
    tig_add(T, "n24q", "O2 (G = 3): beta free in group 2 only",
            tig_join(("groups 3", tig_sec2(1, none), tig_sec2(2, ("beta F2 F1 free")),
                      tig_sec2(3, none))), 7, w, "beta[F2,F1];;group 2;;group 1")
    tig_add(T, "n24r", "O2 (G = 3): beta free in groups 1-2, fixed at 0 in group 3",
            tig_join(("groups 3", tig_sec2(1, ("beta F2 F1 free")), tig_sec2(2, ("beta F2 F1 free")),
                      tig_sec2(3, ("beta F2 F1 fixed 0")))), 7, w, "beta[F2,F1];;group 3")
    tig_add(T, "n24s", "O2 (G = 3): beta free in every group: ACCEPTED (a control)",
            tig_join(("groups 3", tig_sec2(1, ("beta F2 F1 free")), tig_sec2(2, ("beta F2 F1 free")),
                      tig_sec2(3, ("beta F2 F1 free")))), 7, "", "")
    // Widened 2026-10-11 (mgb-fixed-prereg passed): fixed at the SAME nonzero value in
    // every group admitted; any other fixed mix refused.
    w = ("a fixed beta must have the same value in every group")
    tig_add(T, "n24f", "O2: beta fixed at a nonzero value in group 2 (group 1 none)",
            tig_join(("groups 2", tig_sec2(1, none), tig_sec2(2, ("beta F2 F1 fixed 0.3")))),
            6, w, "group 2 fixes beta[F2,F1] at .3, but group 1")
    tig_add(T, "n24j", "O2: beta free in group 1, fixed 0.3 in group 2 (the fixed refusal comes first)",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 free")),
                      tig_sec2(2, ("beta F2 F1 fixed 0.3")))), 6, w,
            "group 2 fixes beta[F2,F1] at .3, but group 1")
    tig_add(T, "n24m", "O2: beta fixed at 0.3 in group 1 and 0.4 in group 2",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 fixed 0.3")),
                      tig_sec2(2, ("beta F2 F1 fixed 0.4")))), 6, w,
            "group 1 fixes beta[F2,F1] at .3, but group 2")
    tig_add(T, "n24n", "O2: beta fixed at 0.3 in group 1 and at 0 in group 2",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 fixed 0.3")),
                      tig_sec2(2, ("beta F2 F1 fixed 0")))), 6, w,
            "group 1 fixes beta[F2,F1] at .3, but group 2")
    tig_add(T, "n24o", "O2: beta fixed at 0.3 in group 1, written differently (.30) in group 2: ACCEPTED (a control)",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 fixed 0.3")),
                      tig_sec2(2, ("beta F2 F1 fixed .30")))), 6, "", "")
    tig_add(T, "n24i", "O2 widened: beta fixed at 0.3 in both groups is ACCEPTED (a control)",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 fixed 0.3")),
                      tig_sec2(2, ("beta F2 F1 fixed 0.3")))), 6, "", "")
    tig_add(T, "n24g", "O2 lifted: beta free in both groups is ACCEPTED (a control)",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 free")),
                      tig_sec2(2, ("beta F2 F1 free")))), 6, "", "")
    tig_add(T, "n24k", "O2 lifted: beta free in both, written before the lambda lines in group 2 (a control)",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 free")), s2[1..3], "beta F2 F1 free",
                      s2[4..cols(s2)])), 6, "", "")
    tig_add(T, "n24l", "O2: beta fixed at 0 in both groups is an absent path, ACCEPTED (a control)",
            tig_join(("groups 2", tig_sec2(1, ("beta F2 F1 fixed 0")),
                      tig_sec2(2, ("beta F2 F1 fixed 0")))), 6, "", "")

    // ---- G7: map count; later-group psi; the existing per-group rules -----------
    w = ("index maps do not match the number of groups")
    tig_add(T, "n25a", "groups 2 with one map", tig_join(("groups 2", s1, tig_s(2, none))), 4, w, "")
    tig_add(T, "n25b", "groups 2 with three maps", tig_join(("groups 2", s1, tig_s(2, none))), 2, w, "")
    tig_add(T, "n25c", "groups 3 with two maps",
            tig_join(("groups 3", s1, tig_s(2, none), tig_s(3, none))), 1, w, "")
    s2 = tig_s(1, none)
    tig_add(T, "n25d", "no groups line (one group) with two maps (an inference; choices item 6)",
            tig_join(s2[2..cols(s2)]), 1, w, "")
    w = ("a group other than 1 has a factor with no psi diagonal")
    tig_add(T, "n26a", "group 2 writes no psi F F",
            tig_join(("groups 2", s1, tig_sx(2, tig_free(5), "# no psi row", none))), 1, w, "")
    tig_add(T, "n26b", "group 2 has psi F1 F1 but no psi F2 F2",
            tig_join(("groups 2", tig_sec2(1, none),
                      "group 2", "factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
                      "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                      "psi F1 F1 fixed 1")), 6, w, "")
    tig_add(T, "n27a", "group 2 factor with no loadings (any rc-3498 text; the wording may name the group)",
            tig_join(("groups 2", s1,
                      tig_sx(2, J(1, 5, "OMIT"), "", none))), 1, "?", "")
    tig_add(T, "n27b", "group 2 names lambda F y1 twice",
            tig_join(("groups 2", s1, tig_s(2, ("lambda F y1 free")))), 1, "?", "")
    tig_add(T, "n27c", "group 1 has no metric line",
            tig_join(("groups 2", tig_sx(1, tig_free(5), "# no metric", none), tig_s(2, none))), 1, "?", "")
    tig_add(T, "n28", "G3: a standalone tau directive inside a group section",
            tig_join(("groups 2", s1, tig_s(2, ("tau y1 1 free")))), 1, "?", "")

    // ---- Stage 1 argument under groups (user, 2026-10-10) --------------------------
    tig_add(T, "n29a", "a groups text with no Stage 1 results",
            tig_join(tig_Econf_lines()), 101, "groups need the Stage 1 results", "")
    tig_add(T, "n29b", "one Stage 1 result for two groups",
            tig_join(tig_Econf_lines()), 201, "one Stage 1 result is needed per group", "")
    tig_add(T, "n29c", "Stage 1 built on other maps (group 2's y1 K differs)",
            tig_join(tig_Econf_lines()), 303,
            "a Stage 1 result does not match its group's index map", "")

    // ---- controls: all ACCEPTED ----------------------------------------------------
    tig_add(T, "k1", "tied lambda y1, equal explicit starts 0.6 and 0.6",
            tig_join(("groups 2", tig_sx(1, tig_with(1, "free 0.6"), "", none),
                      tig_sx(2, tig_with(1, "free 0.6"), "", none), tie1)), 1, "", "")
    tig_add(T, "k2", "tied lambda y1, default start against an explicit 0.7",
            tig_join(("groups 2", s1, tig_sx(2, tig_with(1, "free 0.7"), "", none), tie1)), 1, "", "")
    tig_add(T, "k3", "omega fixed 1.5 in group 2",
            tig_join(("groups 2", s1, tig_s(2, ("omega y1 fixed 1.5")))), 1, "", "")
    tig_add(T, "k4", "omega free, start 1.2, in group 2, with y1's thresholds tied and nu y1 free (identified)",
            tig_join(("groups 2", s1, tig_s(2, ("nu y1 free", "omega y1 free 1.2")),
                      "tie tau 1 y1 1 tau 2 y1 1", "tie tau 1 y1 2 tau 2 y1 2",
                      "tie tau 1 y1 3 tau 2 y1 3")), 1, "", "")
    tig_add(T, "k5", "omega y1 fixed 1 in group 2 (the check-8 spelling of the reference value)",
            tig_join(("groups 2", s1, tig_s(2, ("omega y1 fixed 1")))), 1, "", "")
    tig_add(T, "k6", "tau tie on y2, same K in both groups of the K-mismatch maps",
            tig_join(("groups 2", s1, tig_s(2, none), "tie tau 1 y2 1 tau 2 y2 1")), 3, "", "")
    tig_add(T, "k7", "a three-operand lambda tie, equal starts",
            tig_join(("groups 3", s1, tig_s(2, none), s3,
                      "tie lambda 1 F y1 lambda 2 F y1 lambda 3 F y1")), 2, "", "")
    tig_add(T, "k8", "identical cross-loading tie lambda 1 F2 y1 = lambda 2 F2 y1",
            tig_join(("groups 2", tig_sec2(1, ("lambda F2 y1 free")),
                      tig_sec2(2, ("lambda F2 y1 free")),
                      "tie lambda 1 F2 y1 lambda 2 F2 y1")), 6, "", "")
}

// ---------------------------------------------------------------------------
// Hand-written expectations. Columns: ptype grp row col fidx free value.
// The tau pairs are the 4.2 index-map order of design E: y1-y3 k = 1..3, y4 1, y5 1.
// ---------------------------------------------------------------------------
real matrix tig_taupairs()
{
    return((1, 1 \ 1, 2 \ 1, 3 \ 2, 1 \ 2, 2 \ 2, 3 \ 3, 1 \ 3, 2 \ 3, 3 \ 4, 1 \ 5, 1))
}

real matrix tig_X(real scalar kind)
{
    real matrix tp, X
    real colvector fidx

    tp = tig_taupairs()
    if (kind == 1) {
        // configural: tau g1 1-11, tau g2 12-22, lambda g1 23-27, lambda g2 28-32, two fixed psi
        fidx = (1::22) \ (23::32) \ 0 \ 0
        X = (J(22, 1, 4) \ J(10, 1, 1) \ J(2, 1, 2)),
            (J(11, 1, 1) \ J(11, 1, 2) \ J(5, 1, 1) \ J(5, 1, 2) \ 1 \ 2),
            (tp[., 1] \ tp[., 1] \ (1::5) \ (1::5) \ 1 \ 1),
            (tp[., 2] \ tp[., 2] \ J(10, 1, 1) \ 1 \ 1),
            fidx,
            (J(32, 1, 1) \ 0 \ 0),
            (tig_tv(1, tp) \ tig_tv(2, tp) \ J(10, 1, 0.7) \ 1 \ 1)
        return(X)
    }
    if (kind == 2) {
        // intercepts: tau 1-11 tied, lambda 12-16 tied, psi g2 17, omega 18-22, alpha 23
        fidx = (1::11) \ (1::11) \ (12::16) \ (12::16) \ 0 \ 17 \ (18::22) \ 23
        X = (J(22, 1, 4) \ J(10, 1, 1) \ J(2, 1, 2) \ J(5, 1, 6) \ 7),
            (J(11, 1, 1) \ J(11, 1, 2) \ J(5, 1, 1) \ J(5, 1, 2) \ 1 \ 2 \ J(5, 1, 2) \ 2),
            (tp[., 1] \ tp[., 1] \ (1::5) \ (1::5) \ 1 \ 1 \ (1::5) \ 1),
            (tp[., 2] \ tp[., 2] \ J(10, 1, 1) \ 1 \ 1 \ J(5, 1, 1) \ 1),
            fidx,
            (J(32, 1, 1) \ 0 \ 1 \ J(5, 1, 1) \ 1),
            (tig_tv(1, tp) \ tig_tv(1, tp) \ J(10, 1, 0.7) \ 1 \ 1 \ J(5, 1, 1) \ 0)
        return(X)
    }
    // kind 3: three groups, tie lambda 2 F y1 = lambda 3 F y1
    fidx = (1::33) \ (34::38) \ (39::43) \ (39 \ 44 \ 45 \ 46 \ 47) \ 0 \ 0 \ 0
    X = (J(33, 1, 4) \ J(15, 1, 1) \ J(3, 1, 2)),
        (J(11, 1, 1) \ J(11, 1, 2) \ J(11, 1, 3) \ J(5, 1, 1) \ J(5, 1, 2) \ J(5, 1, 3) \ 1 \ 2 \ 3),
        (tp[., 1] \ tp[., 1] \ tp[., 1] \ (1::5) \ (1::5) \ (1::5) \ 1 \ 1 \ 1),
        (tp[., 2] \ tp[., 2] \ tp[., 2] \ J(15, 1, 1) \ 1 \ 1 \ 1),
        fidx,
        (J(48, 1, 1) \ 0 \ 0 \ 0),
        (tig_tv(1, tp) \ tig_tv(2, tp) \ tig_tv(3, tp) \ J(15, 1, 0.7) \ 1 \ 1 \ 1)
    return(X)
}

string scalar tig_cmp(struct dwlsmv_parmap scalar p, real matrix X)
{
    real colvector g
    real scalar i

    g = tig_grp(p)
    if (p.nall != rows(X)) {
        return("nall is " + strofreal(p.nall) + ", expected " + strofreal(rows(X)))
    }
    for (i = 1; i <= rows(X); i++) {
        if (p.ptype[i] != X[i, 1] | g[i] != X[i, 2] | p.row[i] != X[i, 3] |
            p.col[i] != X[i, 4] | p.fidx[i] != X[i, 5] | p.free[i] != X[i, 6] |
            p.value[i] != X[i, 7]) {
            return("table row " + strofreal(i) + " is (ptype " + strofreal(p.ptype[i]) +
                   ", group " + strofreal(g[i]) + ", row " + strofreal(p.row[i]) +
                   ", col " + strofreal(p.col[i]) + ", fidx " + strofreal(p.fidx[i]) +
                   ", free " + strofreal(p.free[i]) + ", value " + strofreal(p.value[i]) +
                   "), expected (" + strofreal(X[i, 1]) + ", " + strofreal(X[i, 2]) + ", " +
                   strofreal(X[i, 3]) + ", " + strofreal(X[i, 4]) + ", " + strofreal(X[i, 5]) +
                   ", " + strofreal(X[i, 6]) + ", " + strofreal(X[i, 7]) + ")")
        }
    }
    return("")
}

string scalar tig_same(struct dwlsmv_parmap scalar p, struct dwlsmv_parmap scalar q)
{
    if (p.nall != q.nall | p.npar != q.npar) {
        return("nall or npar differs")
    }
    if (p.name != q.name | p.ptype != q.ptype | tig_grp(p) != tig_grp(q) |
        p.row != q.row | p.col != q.col | p.free != q.free | p.value != q.value |
        p.fidx != q.fidx) {
        return("a table field differs")
    }
    return("")
}

real colvector tig_theta0(struct dwlsmv_parmap scalar p)
{
    real colvector th
    real scalar i

    th = J(p.npar, 1, .)
    for (i = 1; i <= p.nall; i++) {
        if (p.free[i]) {
            th[p.fidx[i]] = p.value[i]
        }
    }
    return(th)
}

// Delta restricted to the threshold rows (a selection of the group-major stacked s)
// and the threshold columns is exactly the identity; the threshold rows carry exactly
// zero in every structural column. Exactness follows test_delta.do check 1.
string scalar tig_delta(struct dwlsmv_parmap scalar p, struct dwlsmv_statmap vector m)
{
    real colvector g, trow, tcol, scol
    real matrix Dm, S
    real scalar i, k, off

    g    = tig_grp(p)
    Dm   = dwls_delta(tig_theta0(p), p, m)
    trow = J(0, 1, .)
    tcol = J(0, 1, .)
    scol = J(0, 1, .)
    off  = 0
    for (k = 1; k <= cols(m); k++) {
        for (i = 1; i <= p.nall; i++) {
            if (p.ptype[i] == 4 & g[i] == k) {
                trow = trow \ (off + dwlsmv_pos_tau(m[k], p.row[i], p.col[i]))
                tcol = tcol \ p.fidx[i]
            }
        }
        off = off + m[k].D
    }
    for (i = 1; i <= p.nall; i++) {
        if (p.ptype[i] != 4 & p.free[i]) {
            scol = scol \ p.fidx[i]
        }
    }
    if (cols(Dm) != p.npar) {
        return("Delta has " + strofreal(cols(Dm)) + " columns, npar is " + strofreal(p.npar))
    }
    S = Dm[trow, tcol]
    if (S != I(rows(S))) {
        return("the threshold submatrix of Delta is not exactly the identity")
    }
    if (Dm[trow, scol] != J(rows(trow), rows(scol), 0)) {
        return("a threshold row of Delta is not exactly zero in a structural column")
    }
    return("")
}

// reversed order inside a section, with a comment and a blank line (group line stays first)
string rowvector tig_shuf(string rowvector sec)
{
    string rowvector out
    real scalar j

    out = (sec[1], "# the order inside a section is irrelevant", "")
    for (j = cols(sec); j >= 2; j--) {
        out = out, sec[j]
    }
    return(out)
}

// ---------------------------------------------------------------------------
// G0's second requirement, per corpus text: T alone against a groups-2 text whose
// sections are T and T with an explicit psi F F fixed 1 (the metric lines dropped).
// Compared: ptype, row, col, free, value, and the ORDER of fidx (fidx itself differs
// by construction under O6). Corpus texts have no beta and no marker metric (O2) and
// are chosen to pass the rank check at theta0 (no second-point start; choices item 8).
// ---------------------------------------------------------------------------
string scalar tig_g0eq(real scalar k)
{
    struct dwlsmv_parmap scalar p, q
    string rowvector T, fac, g2
    real colvector Kv, g, ix, fq, fp, fqf
    real scalar j, n4

    if (k == 1) {
        T   = ("factor f", "ind y1 y2 y3 y4",
               "lambda f y1 free 0.75", "lambda f y2 free 0.65",
               "lambda f y3 free 0.55", "lambda f y4 free 0.45",
               "metric unitvariance f")
        fac = ("f")
        Kv  = (4 \ 4 \ 3 \ 3)
    }
    else if (k == 2) {
        T   = ("factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
               "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
               "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
               "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
               "psi F1 F3 free", "psi F3 F2 free", "psi F1 F2 free",
               "metric unitvariance F1", "metric unitvariance F2", "metric unitvariance F3")
        fac = ("F1", "F2", "F3")
        Kv  = J(9, 1, 3)
    }
    else {
        T   = ("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
               "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
               "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
               "psi F1 F2 free",
               "metric unitvariance F1", "metric unitvariance F2")
        fac = ("F1", "F2")
        Kv  = J(6, 1, 3)
    }
    g2 = J(1, 0, "")
    for (j = 1; j <= cols(T); j++) {
        if (substr(T[j], 1, 6) != "metric") {
            g2 = g2, T[j]
        }
    }
    for (j = 1; j <= cols(fac); j++) {
        g2 = g2, ("psi " + fac[j] + " " + fac[j] + " fixed 1")
    }
    p = tig_ident(tig_join(T), tig_mapsK(Kv))
    q = tig_ident(tig_join(("groups 2", "group 1", T, "group 2", g2)),
                      tig_mapsK((Kv, Kv)))
    g  = tig_grp(q)
    ix = select(1::q.nall, g :== 1)
    if (rows(ix) != p.nall) {
        return("group 1 has " + strofreal(rows(ix)) + " rows, the single-group table " +
               strofreal(p.nall))
    }
    // Threshold STARTS differ by decision (2026-10-10): 0 for a single group, the
    // sample thresholds under groups (O5, the C3 note). Every other field matches.
    if (q.ptype[ix] != p.ptype | q.row[ix] != p.row | q.col[ix] != p.col |
        q.free[ix] != p.free | (q.value[ix] :* (p.ptype :!= 4)) != (p.value :* (p.ptype :!= 4))) {
        return("group 1's ptype, row, col, free or value differ from the single-group table")
    }
    for (j = 1; j <= rows(ix); j++) {
        if (p.ptype[j] == 4) {
            if (q.value[ix[j]] != tig_tau(1, p.row[j], p.col[j], Kv[p.row[j]])) {
                return("group 1's threshold start is not its sample threshold")
            }
        }
    }
    fq  = q.fidx[ix]
    fp  = p.fidx
    if ((fq :== 0) != (fp :== 0)) {
        return("group 1 and the single-group table disagree on which rows are free")
    }
    n4 = sum(p.ptype :== 4)
    if (fq[1::n4] != fp[1::n4]) {
        return("group 1's threshold fidx is not the single-group threshold fidx")
    }
    fqf = select(fq, fq :> 0)
    if (rows(fqf) > 1) {
        if (any(fqf[2::rows(fqf)] :<= fqf[1::(rows(fqf) - 1)])) {
            return("group 1's free rows are not numbered in the single-group order")
        }
    }
    if (q.npar != 2 * p.npar) {
        return("npar is " + strofreal(q.npar) + ", expected twice " + strofreal(p.npar))
    }
    return("")
}

real scalar tig_tables()
{
    struct dwlsmv_parmap scalar p, q
    struct dwlsmv_statmap vector m2, m3
    string colvector FAIL
    string rowvector s1, s2
    string scalar r
    real scalar k

    FAIL = J(0, 1, "")
    m2   = tig_maps(1)
    m3   = tig_maps(2)
    s1   = tig_s(1, tig_none())
    s2   = tig_s(2, tig_none())

    // ---- t1: design E configural -------------------------------------------------
    p = tig_ident(tig_join(tig_Econf_lines()), m2)
    printf("{txt}  t1 configural: nall = %f, npar = %f (spec: 34, 32)\n", p.nall, p.npar)
    if (p.nall != 34 | p.npar != 32) {
        FAIL = FAIL \ "t1: design E configural is not nall 34, npar 32"
    }
    r = tig_cmp(p, tig_X(1))
    if (r != "") {
        FAIL = FAIL \ ("t1: " + r)
    }

    // ---- t3: its Delta -----------------------------------------------------------
    r = tig_delta(p, m2)
    printf("{txt}  t3 configural Delta, threshold rows x columns: %s\n", (r == "" ? "identity" : r))
    if (r != "") {
        FAIL = FAIL \ ("t3: " + r)
    }

    // ---- t5-t8: textual freedom --------------------------------------------------
    q = tig_ident(tig_join((s1, s2, "groups 2")), m2)
    r = tig_same(p, q)
    if (r != "") {
        FAIL = FAIL \ ("t5 groups line last: " + r)
    }
    q = tig_ident(tig_join((s1, "groups 2", s2)), m2)
    r = tig_same(p, q)
    if (r != "") {
        FAIL = FAIL \ ("t6 groups line between the sections: " + r)
    }
    q = tig_ident(tig_join(("groups 2", s2, s1)), m2)
    r = tig_same(p, q)
    if (r != "") {
        FAIL = FAIL \ ("t7 sections in reverse textual order: " + r)
    }
    q = tig_ident(tig_join(("groups 2", tig_shuf(s1), tig_shuf(s2))), m2)
    r = tig_same(p, q)
    if (r != "") {
        FAIL = FAIL \ ("t8 lines reversed inside each section: " + r)
    }

    // ---- t2: design E, intercepts level -------------------------------------------
    p = tig_ident(tig_join(tig_Eint_lines()), m2)
    printf("{txt}  t2 intercepts: nall = %f, npar = %f (spec: 40, 23)\n", p.nall, p.npar)
    if (p.nall != 40 | p.npar != 23) {
        FAIL = FAIL \ "t2: design E intercepts is not nall 40, npar 23"
    }
    r = tig_cmp(p, tig_X(2))
    if (r != "") {
        FAIL = FAIL \ ("t2: " + r)
    }
    if (rows(tig_ties_int()') != 16) {
        FAIL = FAIL \ "t2: the helper did not write the 11 + 5 tie lines of spec G8"
    }

    // ---- t4: a tie whose lowest member is not in group 1 ----------------------------
    p = tig_ident(tig_join(("groups 3", s1, s2, tig_s(3, tig_none()),
                                "tie lambda 2 F y1 lambda 3 F y1")), m3)
    printf("{txt}  t4 three groups, one lambda tie: nall = %f, npar = %f (expected 51, 47)\n",
           p.nall, p.npar)
    if (p.nall != 51 | p.npar != 47) {
        FAIL = FAIL \ "t4: nall or npar is not 51, 47"
    }
    r = tig_cmp(p, tig_X(3))
    if (r != "") {
        FAIL = FAIL \ ("t4: " + r)
    }

    // ---- eq1-eq3 -----------------------------------------------------------------------
    for (k = 1; k <= 3; k++) {
        r = tig_g0eq(k)
        printf("{txt}  eq%f G0 equivalence: %s\n", k, (r == "" ? "same rows" : r))
        if (r != "") {
            FAIL = FAIL \ ("eq" + strofreal(k) + ": " + r)
        }
    }

    printf("{txt}test_identify_groups tables: %f failure(s)\n", rows(FAIL))
    for (k = 1; k <= rows(FAIL); k++) {
        printf("{err}    %s\n", FAIL[k])
    }
    return(rows(FAIL))
}

end

* ===========================================================================
* The case table: every row of G7, G0's token check, O2-O9, and the controls.
* ===========================================================================
mata: TAB = J(0, 6, "")
mata: tig_cases(TAB)
mata: st_local("ncase", strofreal(rows(TAB)))

local nbadfail 0
tempfile lg
forvalues i = 1/`ncase' {
    mata: st_local("id", TAB[`i', 1])
    mata: st_local("why", TAB[`i', 2])
    quietly log using "`lg'", text replace name(tig)
    capture noisily mata: tig_try(TAB[`i', 3], strtoreal(TAB[`i', 4]))
    local rc = _rc
    quietly log close tig
    mata: st_local("verdict", tig_judge(TAB[`i', .], `rc', tig_slurp(st_local("lg"))))
    if `"`verdict'"' == "" {
        display as text "    `id' (rc `rc'): `why'"
    }
    else {
        display as error `"    `id' FAILED (rc `rc'): `why' -- `verdict'"'
        local nbadfail = `nbadfail' + 1
    }
}

mata: st_local("nfail", strofreal(tig_tables()))

display as text "    O5(a) is covered by t2: tied thresholds start at group 1's synthetic sample threshold, which differs from group 2's."
display as text "    M62, the stacked-only rank check: design under numerical verification (user, 2026-10-10); not yet a case."

local ntot = `nfail' + `nbadfail'
if `ntot' {
    display as error _n "test_identify_groups: `nfail' table failure(s), `nbadfail' case failure(s) of `ncase'"
    exit 9
}

display as result "test_identify_groups: every G7 rejection is rc 3498 with its text, a text without groups refuses all five tokens, design E configural and intercepts give nall 34/40 and npar 32/23 with the O6 positions and an exactly-identity threshold submatrix of Delta, ties share one fidx at the lowest member, group 1 refuses nu/omega/alpha in every status, a beta is admitted free in every group, fixed at the same value in every group, or at G = 2 free in one group and absent in the other (at G > 2 one-group, and at any G mixed and differing fixed betas, refused), and a marker metric is refused in every groups text"
exit 0
