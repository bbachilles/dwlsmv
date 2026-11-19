*! Drafted by a read-only agent 2026-11-18 from E-110 and the owner decisions of
*! 2026-11-18 (C1 (a), C2 (a), M1-M13); run and fixed by the main session the same day: C2 at every G, user-fixed
*! models dwlsmv accepts, MathAb <- MathAtt, equal numbers on two rows, a tempfile log. 26 mutants killed
*! (E-110).
*!
*! The compact grouped coefficient table (src/mata/dwlsmv_m11_ctab.mata, _dwlsmv_table and _dwlsmv_nctable in
*! dwlsmv.ado). End to end through build/ldwlsmv.mlib and src/ado, as test_efa_display.do. Display checks read the
*! command's own log text (tct_run keeps it and r(table) of the replay); a phrase test says what is in the output, a
*! number test builds the expected string here from e(b), e(V) and the normal quantile.
*!
*! Tolerances, fixed before the first run (2026-11-18, chosen before measuring):
*!   r(table) against e(b) and e(V): b, se, crit, eform exact (the same copy, the same sqrt); z, pvalue, ll, ul 1e-14
*!     relative to the value recomputed here (the DECIDING criterion for the numbers);
*!   r(table) against Stata's own ereturn display r(table) on the free columns: 1e-12 relative (a DIAGNOSTIC where
*!     ereturn display posts no r(table); a value mismatch is a failure);
*!   printed numbers: the exact string sprintf("%11.4f%10.4f%8.2f%8.3f%12.4f%12.4f") of b, se, z, p, ll, ul built here.
*!
*! Checks
*!   1  display only: e() (every macro, scalar and matrix) is byte-identical after replays with each display option, and
*!      lincom still runs
*!   2  numbers: r(table) at levels 90, 95, 99 on a marker CFA, a unit-variance CFA, a structural fit, a constrained
*!      fit, G = 2 at every invariance level, a covariate fit and the standardized table; the printed rows (label and
*!      numbers on one line) at level 90, 95 and 99; the independent comparator (ereturn display)
*!   3  row-set accounting: every e(b) column is shown exactly once or omitted with a documented reason; a row's number
*!      is the one at its own e(b) index (the reorder trap); only fixed rows are omitted for a reason other than
*!      thresholds
*!   4  widths at linesize 78, 80 and 120 with a 28-character name: 78 columns wide, the label column grows up to
*!      linesize - 62, a long label is cut with a tilde
*!   5  the layouts, cases (a) to (g): (a) marker CFA, (b) unit-variance CFA, (c) structural fit, marker and unit
*!      variance, and its standardized table, (d) user constraints, (e) G = 2 at four invariance levels, standardized
*!      at G = 2, (f) a covariate, (g) a fit that did not converge (real, and hand-set at G = 2)
*!   6  threshold defaults and overrides, G = 1 and G = 2
*!   7  error paths: rc 198 exactly for contradictory or foreign display options and for standardized on a fit that did
*!      not converge; rc 3498 exactly for malformed planner input, each with its own message
*!   8  replay = estimation: the table text (header to the thresholds note) is the same at estimation and on replay,
*!      raw, standardized, and at G = 2
*!
*! What this does not certify: that the table reads well. The numbers are certified, the layout is checked by phrase.
*!
*! Return codes: 0 pass, 9 fail.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 120

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

global TCT_NF 0
tempfile tctlg
global TCT_LG "`tctlg'"

* runs a command with its output in the log file $TCT_LG; returns the command's rc, and keeps r(table) of the
* command (copied while the log is open) in the matrix tct_RT; r(hasrt) says whether there was one
capture program drop tct_run
program tct_run, rclass
    capture matrix drop tct_RT
    capture quietly log close tctlog
    quietly log using "$TCT_LG", text replace name(tctlog)
    capture noisily `0'
    local rc = _rc
    capture matrix tct_RT = r(table)
    local hasrt = (_rc == 0)
    quietly log close tctlog
    return scalar rc = `rc'
    return scalar hasrt = `hasrt'
end

capture program drop tct_chk
program tct_chk
    args want phrase label
    mata: st_local("got", strofreal(tct_has(st_local("phrase"))))
    if `got' != `want' {
        display as error "test_ctab: `label' (phrase [`phrase'] found `got', expected `want')"
        global TCT_NF = $TCT_NF + 1
    }
end

* the number of times a phrase occurs in the log
capture program drop tct_n
program tct_n
    args want phrase label
    mata: st_local("got", strofreal(tct_count(st_local("phrase"))))
    if `got' != `want' {
        display as error "test_ctab: `label' (phrase [`phrase'] occurs `got' times, expected `want')"
        global TCT_NF = $TCT_NF + 1
    }
end

* the phrases occur in the log in this order (up to six)
capture program drop tct_ord
program tct_ord
    args label p1 p2 p3 p4 p5 p6
    mata: st_local("ok", strofreal(tct_ordered(st_local("p1"), st_local("p2"), st_local("p3"), st_local("p4"), st_local("p5"), st_local("p6"))))
    if !`ok' {
        display as error "test_ctab: `label' (the phrases are not in the order [`p1'] [`p2'] [`p3'] [`p4'] [`p5'] [`p6'])"
        global TCT_NF = $TCT_NF + 1
    }
end

capture program drop tct_ok
program tct_ok
    args cond label
    if !(`cond') {
        display as error "test_ctab: `label'"
        global TCT_NF = $TCT_NF + 1
    }
end

* a fit with e(converged) set to 0 (display only), as test_unconv_post.do
capture program drop tct_unconv
program tct_unconv, eclass
    ereturn scalar converged = 0
end

mata:
string colvector tct_lines()
{
    return(cat(st_global("TCT_LG")))
}

real scalar tct_has(string scalar phrase)
{
    return(strpos(invtokens(tct_lines()', " "), phrase) > 0)
}

real scalar tct_count(string scalar phrase)
{
    string scalar s

    s = invtokens(tct_lines()', char(10))
    return((strlen(s) - strlen(subinstr(s, phrase, "", .))) / strlen(phrase))
}

real scalar tct_ordered(string scalar a, string scalar b, string scalar c, string scalar d,
                        string scalar e, string scalar f)
{
    string rowvector p
    string scalar s
    real scalar i, at, prev

    p    = (a, b, c, d, e, f)
    s    = invtokens(tct_lines()', " ")
    prev = 0
    for (i = 1; i <= cols(p); i++) {
        if (p[i] != "") {
            at = strpos(s, p[i])
            if (at == 0 | at <= prev) {
                return(0)
            }
            prev = at
        }
    }
    return(1)
}

// everything e() holds, as one string: macros, scalars (hex), matrices (hex), so that equality is byte-exact
string scalar tct_snap()
{
    string colvector nm
    string scalar s
    real scalar i
    real matrix M

    s  = ""
    nm = st_dir("e()", "macro", "*")
    for (i = 1; i <= rows(nm); i++) {
        s = s + nm[i] + "=" + st_global("e(" + nm[i] + ")") + char(10)
    }
    nm = st_dir("e()", "numscalar", "*")
    for (i = 1; i <= rows(nm); i++) {
        s = s + nm[i] + "=" + strofreal(st_numscalar("e(" + nm[i] + ")"), "%21x") + char(10)
    }
    nm = st_dir("e()", "matrix", "*")
    for (i = 1; i <= rows(nm); i++) {
        M = st_matrix("e(" + nm[i] + ")")
        s = s + nm[i] + ":" + strofreal(rows(M)) + "x" + strofreal(cols(M)) + ":"
        if (rows(M) * cols(M) > 0) {
            s = s + invtokens(strofreal(vec(M), "%21x")')
        }
        s = s + char(10)
    }
    return(s)
}

real scalar tct_reldif(real scalar x, real scalar y)
{
    return(abs(x - y) / max((1, abs(y))))
}

// r(table), kept in tct_RT, against the matrices bn, vn and the level lv; returns failures
real scalar tct_numeq(real scalar lv, string scalar bn, string scalar vn, string scalar tag)
{
    real matrix T, V
    real rowvector b
    real scalar i, nb, se, z, crit, bad

    T   = st_matrix("tct_RT")
    b   = st_matrix(bn)
    V   = st_matrix(vn)
    nb  = cols(b)
    bad = 0
    if (rows(T) != 9 | cols(T) != nb) {
        printf("{err}  %s: r(table) is %f x %f, wanted 9 x %f\n", tag, rows(T), cols(T), nb)
        return(1)
    }
    if (st_matrixrowstripe("tct_RT")[., 2]' != ("b", "se", "z", "pvalue", "ll", "ul", "df", "crit", "eform")) {
        printf("{err}  %s: r(table) row names\n", tag)
        bad++
    }
    if (st_matrixcolstripe("tct_RT") != st_matrixcolstripe(bn)) {
        printf("{err}  %s: r(table) column names are not the stripe of b\n", tag)
        bad++
    }
    crit = invnormal(.5 + lv / 200)
    for (i = 1; i <= nb; i++) {
        se = sqrt(V[i, i])
        if (T[1, i] != b[i] | T[2, i] != se | T[8, i] != crit | T[9, i] != 0 | !missing(T[7, i])) {
            printf("{err}  %s: column %f b, se, df, crit or eform\n", tag, i)
            bad++
        }
        if (se > 0) {
            z = b[i] / se
            if (tct_reldif(T[3, i], z) > 1e-14 | tct_reldif(T[4, i], 2 * normal(-abs(z))) > 1e-14 |
                tct_reldif(T[5, i], b[i] - crit * se) > 1e-14 | tct_reldif(T[6, i], b[i] + crit * se) > 1e-14) {
                printf("{err}  %s: column %f z, pvalue or interval at level %f\n", tag, i, lv)
                bad++
            }
        }
        else {
            if (!missing(T[3, i]) | !missing(T[4, i]) | !missing(T[5, i]) | !missing(T[6, i])) {
                printf("{err}  %s: column %f has a zero se and a z, p or interval\n", tag, i)
                bad++
            }
        }
    }
    return(bad)
}

// the same r(table) against Stata's own (ereturn display's), kept in tct_ED, on the free columns; returns failures,
// or -1 when tct_ED does not exist (diagnostic)
real scalar tct_vsed(string scalar vn, string scalar tag)
{
    real matrix T, E, V
    string matrix sT, sE
    real scalar i, j, r, nb, bad, hit

    if (rows(st_matrix("tct_ED")) == 0) {
        return(-1)
    }
    T  = st_matrix("tct_RT")
    E  = st_matrix("tct_ED")
    V  = st_matrix(vn)
    sT = st_matrixcolstripe("tct_RT")
    sE = st_matrixcolstripe("tct_ED")
    nb = cols(T)
    bad = 0
    for (i = 1; i <= nb; i++) {
        if (V[i, i] > 0) {
            hit = 0
            for (j = 1; j <= cols(E); j++) {
                if (sE[j, 1] == sT[i, 1] & sE[j, 2] == sT[i, 2]) {
                    hit = j
                }
            }
            if (hit == 0) {
                printf("{err}  %s: ereturn display has no column %s:%s\n", tag, sT[i, 1], sT[i, 2])
                bad++
            }
            else {
                for (r = 1; r <= 6; r++) {
                    if (tct_reldif(T[r, i], E[r, hit]) > 1e-12) {
                        printf("{err}  %s: column %f row %f differs from ereturn display (%g, %g)\n", tag, i, r, T[r, i], E[r, hit])
                        bad++
                    }
                }
            }
        }
    }
    return(bad)
}

// the printed rows of the current raw fit (G = 1): every free column shown has a line with its numbers, whose label is the
// column's label; at a linesize where the label column is 16 wide
real scalar tct_rows(real scalar lv, real scalar showthr)
{
    real rowvector b
    real matrix V
    string matrix stp
    string rowvector fac, ind
    string colvector L
    real scalar i, k, nb, se, z, crit, bad, hit, isthr, p, ndup
    string scalar body, lab, eq, nm, got, x

    b   = st_matrix("e(b)")
    V   = st_matrix("e(V)")
    stp = st_matrixcolstripe("e(b)")
    ind = tokens(st_global("e(indicators)"))
    nb  = cols(b)
    fac = J(1, 0, "")
    for (i = 1; i <= nb; i++) {
        if (stp[i, 1] == "/" & substr(stp[i, 2], 1, 4) == "var(") {
            x = substr(stp[i, 2], 5, strlen(stp[i, 2]) - 5)
            if (substr(x, 1, 2) == "e.") {
                x = substr(x, 3, .)
            }
            if (!anyof(ind, x) & !anyof(fac, x)) {
                fac = fac, x
            }
        }
    }
    L    = tct_lines()
    crit = invnormal(.5 + lv / 200)
    bad  = 0
    for (i = 1; i <= nb; i++) {
        eq    = stp[i, 1]
        nm    = stp[i, 2]
        isthr = (substr(eq, 1, 1) == "/" & eq != "/")
        if (V[i, i] > 0 & (showthr | !isthr)) {
            se = sqrt(V[i, i])
            z  = b[i] / se
            p  = 2 * normal(-abs(z))
            if (eq == "/" | isthr) {
                lab = nm
            }
            else if (anyof(ind, eq) & anyof(fac, nm)) {
                lab = eq
            }
            else {
                lab = eq + " <- " + nm
            }
            body = sprintf("%11.4f%10.4f%8.2f%8.3f%12.4f%12.4f", b[i], se, z, p, b[i] - crit * se, b[i] + crit * se)
            hit = 0
            for (k = 1; k <= rows(L); k++) {
                if (strpos(L[k], body) > 0) {
                    hit++
                    got = strtrim(substr(L[k], 1, 15))
                    if (got != lab) {
                        printf("{err}  row %s:%s printed under the label [%s], wanted [%s]\n", eq, nm, got, lab)
                        bad++
                    }
                }
            }
            // two shown parameters with the same numbers (y5 and y6's cut2 on the example data) are on
            // as many lines as there are such parameters (main session, first run 2026-11-18)
            ndup = 0
            for (k = 1; k <= nb; k++) {
                if (V[k, k] > 0 & (showthr | !(substr(stp[k, 1], 1, 1) == "/" & stp[k, 1] != "/"))) {
                    x = sprintf("%11.4f%10.4f%8.2f%8.3f%12.4f%12.4f", b[k], sqrt(V[k, k]), b[k] / sqrt(V[k, k]),
                        2 * normal(-abs(b[k] / sqrt(V[k, k]))), b[k] - crit * sqrt(V[k, k]), b[k] + crit * sqrt(V[k, k]))
                    if (x == body) {
                        ndup++
                    }
                }
            }
            if (hit != ndup) {
                printf("{err}  row %s:%s: its numbers are on %f lines, wanted %f\n", eq, nm, hit, ndup)
                bad++
            }
        }
    }
    return(bad)
}

// plan from the matrices bn, vn ("" for none), cn ("" for none) and the current e() metadata
void tct_plan(string scalar bn, string scalar vn, string scalar cn, real scalar showthr, real scalar std)
{
    external string matrix TP_ptxt, TP_otxt
    external real matrix TP_pnum, TP_tab
    external real colvector TP_oidx
    real matrix V, C, K
    string matrix K_s
    string rowvector gv
    real scalar hasv, hasc

    hasv = (vn != "")
    hasc = (cn != "")
    V    = J(0, 0, .)
    C    = J(0, 0, .)
    if (hasv) {
        V = st_matrix(vn)
    }
    if (hasc) {
        C = st_matrix(cn)
    }
    gv = J(1, 0, "")
    if (st_global("e(invariance)") != "") {
        K_s = st_matrixrowstripe("e(K)")
        gv  = K_s[., 2]'
    }
    dwlsmv_ct_plan(st_matrixcolstripe(bn), st_matrix(bn), V, hasv, C, hasc, st_global("e(metric)"), st_global("e(canon)"),
                   tokens(st_global("e(indicators)")), gv, st_global("e(groupvar)"), showthr, std, 95,
                   TP_ptxt, TP_pnum, TP_otxt, TP_oidx, TP_tab)
}

// accounting on the plan last made by tct_plan. allowed: the reasons an omission may have; nmk, ncs: the wanted
// numbers of marker and constrained rows (-1: not checked); nomit: the wanted number of omissions (-1: not checked).
// Every b column once; a shown row's numbers are tab's at ITS OWN index (b and the square root of V's diagonal there).
real scalar tct_acct(string scalar bn, string scalar vn, string scalar cn, string scalar allowed, real scalar nmk,
                     real scalar ncs, real scalar nomit, string scalar tag)
{
    external string matrix TP_ptxt, TP_otxt
    external real matrix TP_pnum, TP_tab
    external real colvector TP_oidx
    real rowvector b, seen, fixed
    real matrix V, C
    real scalar i, j, k, nb, bad, bi, nz
    string rowvector ok

    b    = st_matrix(bn)
    nb   = cols(b)
    V    = st_matrix(vn)
    seen = J(1, nb, 0)
    bad  = 0
    fixed = J(1, nb, 0)
    if (cn != "") {
        C = st_matrix(cn)
        for (k = 1; k <= rows(C); k++) {
            for (j = 1; j <= nb; j++) {
                if (C[k, j] != 0) {
                    fixed[j] = 1
                }
            }
        }
    }
    for (i = 1; i <= rows(TP_pnum); i++) {
        bi = TP_pnum[i, 1]
        seen[bi] = seen[bi] + 1
        if (TP_ptxt[i, 3] == "" | TP_ptxt[i, 1] == "") {
            printf("{err}  %s: row %f has no label or section\n", tag, i)
            bad++
        }
        if (TP_pnum[i, 3] != b[bi] | TP_tab[1, bi] != b[bi]) {
            printf("{err}  %s: row %f (b index %f) is not that estimate\n", tag, i, bi)
            bad++
        }
        if (TP_pnum[i, 2] == 1) {
            if (fixed[bi]) {
                printf("{err}  %s: a fixed column (%f) is shown as free\n", tag, bi)
                bad++
            }
            if (TP_pnum[i, 4] != sqrt(V[bi, bi]) | TP_pnum[i, 4] != TP_tab[2, bi]) {
                printf("{err}  %s: row %f (b index %f) carries a standard error that is not its own\n", tag, i, bi)
                bad++
            }
        }
        else {
            if (!fixed[bi]) {
                printf("{err}  %s: a free column (%f) is shown as fixed\n", tag, bi)
                bad++
            }
        }
    }
    ok = tokens(allowed)
    for (i = 1; i <= rows(TP_oidx); i++) {
        bi = TP_oidx[i]
        seen[bi] = seen[bi] + 1
        if (!anyof(ok, TP_otxt[i, 2])) {
            printf("{err}  %s: omitted %s for the reason [%s], not in [%s]\n", tag, TP_otxt[i, 1], TP_otxt[i, 2], allowed)
            bad++
        }
        if (TP_otxt[i, 2] != "thresholds" & !fixed[bi]) {
            printf("{err}  %s: a free column (%s) is omitted for [%s]\n", tag, TP_otxt[i, 1], TP_otxt[i, 2])
            bad++
        }
    }
    if (any(seen :!= 1)) {
        printf("{err}  %s: %f columns are not accounted for exactly once\n", tag, sum(seen :!= 1))
        bad++
    }
    if (nmk >= 0) {
        if (sum(TP_pnum[., 2] :== 2) != nmk) {
            printf("{err}  %s: %f marker rows, wanted %f\n", tag, sum(TP_pnum[., 2] :== 2), nmk)
            bad++
        }
    }
    if (ncs >= 0) {
        if (sum(TP_pnum[., 2] :== 3) != ncs) {
            printf("{err}  %s: %f constrained rows, wanted %f\n", tag, sum(TP_pnum[., 2] :== 3), ncs)
            bad++
        }
    }
    if (nomit >= 0) {
        if (rows(TP_oidx) != nomit) {
            printf("{err}  %s: %f omitted rows, wanted %f\n", tag, rows(TP_oidx), nomit)
            bad++
        }
    }
    return(bad)
}

// the number of omissions with a reason
real scalar tct_nreason(string scalar reason)
{
    external string matrix TP_otxt

    return(sum(TP_otxt[., 2] :== reason))
}

// the number of threshold columns of e(b)
real scalar tct_nthr()
{
    string matrix stp

    stp = st_matrixcolstripe("e(b)")
    return(sum((substr(stp[., 1], 1, 1) :== "/") :& (stp[., 1] :!= "/")))
}

// the number of fixed columns of the matrix cn
real scalar tct_nfix(string scalar cn)
{
    return(rows(st_matrix(cn)))
}

// the table lines of the log: the header line to the end, without the notes under the table
string colvector tct_tablelines()
{
    string colvector L, out
    real scalar k, at

    L  = tct_lines()
    at = 0
    for (k = 1; k <= rows(L); k++) {
        if (at == 0 & strpos(L[k], "Estimate") > 0) {
            at = k
        }
    }
    out = J(0, 1, "")
    if (at == 0) {
        return(out)
    }
    for (k = at; k <= rows(L); k++) {
        if (substr(L[k], 1, 10) != "Thresholds" & substr(L[k], 1, 8) != "Fixed by" & substr(L[k], 1, 5) != "Note:") {
            out = out \ L[k]
        }
    }
    return(out)
}

// the widest table line and whether every one is within the linesize
real rowvector tct_width(real scalar ls)
{
    string colvector T
    real scalar k, wmax

    T    = tct_tablelines()
    wmax = 0
    for (k = 1; k <= rows(T); k++) {
        wmax = max((wmax, ustrlen(T[k])))
    }
    return((wmax, rows(T) > 0 & wmax <= ls))
}

// the text of the table, from its column titles to the thresholds note, of the log; "" if there is no such span
string scalar tct_segment()
{
    string colvector L
    real scalar k, a, z
    string scalar s

    L = tct_lines()
    a = 0
    z = 0
    for (k = 1; k <= rows(L); k++) {
        if (a == 0 & strpos(L[k], "Estimate") > 0) {
            a = k
        }
        if (a > 0 & z == 0 & substr(L[k], 1, 24) == "Thresholds are not shown") {
            z = k
        }
    }
    if (a == 0 | z == 0) {
        return("")
    }
    s = ""
    for (k = a; k <= z; k++) {
        s = s + L[k] + char(10)
    }
    return(s)
}

// the i-th malformed planner input (G = 1 fit in e()); the planner must abort at 3498 with its own message
void tct_bad(real scalar c)
{
    real rowvector b
    real matrix V, C, pnum, tab
    string matrix stp, ptxt, otxt
    string rowvector ind, gv
    real colvector oidx
    real scalar lv, il, nb, iv

    b   = st_matrix("e(b)")
    V   = st_matrix("e(V)")
    C   = st_matrix("e(Cns)")
    stp = st_matrixcolstripe("e(b)")
    ind = tokens(st_global("e(indicators)"))
    gv  = J(1, 0, "")
    lv  = 95
    nb  = cols(b)
    il  = selectindex(stp[., 1] :== "y2")[1]
    iv  = selectindex(stp[., 2] :== "var(F1)")[1]
    if (c == 1) {
        stp[il, 1] = "zzz"
    }
    else if (c == 2) {
        stp[il, .] = stp[il - 1, .]
    }
    else if (c == 3) {
        V = V[1..nb - 1, 1..nb - 1]
    }
    else if (c == 4) {
        b[il] = .
    }
    else if (c == 5) {
        C[1, il] = 1
    }
    else if (c == 6) {
        C[1, .] = J(1, cols(C), 0)
    }
    else if (c == 7) {
        C = C \ C[1, .]
    }
    else if (c == 8) {
        C = C[., 1..cols(C) - 1]
    }
    else if (c == 9) {
        lv = 100
    }
    else if (c == 10) {
        stp[iv, 2] = "foo(F1)"
    }
    else if (c == 11) {
        stp = stp[1..nb - 1, .]
    }
    else if (c == 12) {
        V[il, il] = .
    }
    dwlsmv_ct_plan(stp, b, V, 1, C, 1, st_global("e(metric)"), st_global("e(canon)"), ind, gv, "", 0, 0, lv,
                   ptxt, pnum, otxt, oidx, tab)
}

// the i-th malformed planner input on the G = 2 fit in e()
void tct_badg(real scalar c)
{
    real rowvector b
    real matrix V, C, pnum, tab, K
    string matrix stp, ptxt, otxt, K_s
    string rowvector ind, gv
    string scalar gvar
    real colvector oidx
    real scalar il

    b    = st_matrix("e(b)")
    V    = st_matrix("e(V)")
    C    = st_matrix("e(Cns)")
    stp  = st_matrixcolstripe("e(b)")
    ind  = tokens(st_global("e(indicators)"))
    K_s  = st_matrixrowstripe("e(K)")
    gv   = K_s[., 2]'
    gvar = st_global("e(groupvar)")
    il   = selectindex(stp[., 1] :== "y2")[1]
    if (c == 1) {
        stp[il, 2] = "c.F"
    }
    else if (c == 2) {
        gvar = ""
    }
    else if (c == 3) {
        gv[2] = ""
    }
    dwlsmv_ct_plan(stp, b, V, 1, C, 1, st_global("e(metric)"), st_global("e(canon)"), ind, gv, gvar, 1, 0, 95,
                   ptxt, pnum, otxt, oidx, tab)
}

// the loadings of group 1 and group 2 equal, label by label (tied at metric and scalar), from the plan last made
real scalar tct_tied()
{
    external string matrix TP_ptxt
    external real matrix TP_pnum
    real scalar i, j, bad

    bad = 0
    for (i = 1; i <= rows(TP_pnum); i++) {
        if (TP_pnum[i, 11] == 1 & TP_pnum[i, 10] == 1) {
            for (j = 1; j <= rows(TP_pnum); j++) {
                if (TP_pnum[j, 11] == 1 & TP_pnum[j, 10] == 2 & TP_ptxt[j, 3] == TP_ptxt[i, 3]) {
                    if (abs(TP_pnum[i, 3] - TP_pnum[j, 3]) > 1e-12 | abs(TP_pnum[i, 4] - TP_pnum[j, 4]) > 1e-12) {
                        bad++
                    }
                }
            }
        }
    }
    return(bad)
}
end

* a matrix of the free-parameter display numbers for the fits below is read from r(table); the checks of 2 and 3
* run on whatever e() holds, through this program. tag names the fit; ng1: 1 at one group (the printed rows are
* checked at one group only: the label column of the rows check is the 16-wide one)
capture program drop tct_numbers
program tct_numbers
    args tag ng1
    tempname hold
    foreach lv in 90 95 99 {
        set linesize 80
        tct_run dwlsmv, level(`lv') thresholds
        tct_ok "r(rc) == 0 & r(hasrt) == 1" "2 `tag': replay at level `lv' rc and r(table)"
        mata: st_local("f", strofreal(tct_numeq(`lv', "e(b)", "e(V)", "2 `tag' `lv'")))
        tct_ok "`f' == 0" "2 `tag': r(table) at level `lv' against e(b) and e(V), `f' failures"
        tct_chk 1 "[`lv'% conf. interval]" "2 `tag': the interval header names level `lv'"
        if "`ng1'" == "1" {
            mata: st_local("f", strofreal(tct_rows(`lv', 1)))
            tct_ok "`f' == 0" "2 `tag': printed rows at level `lv', `f' failures"
        }
        * the comparator that shares no code with ours: Stata's own table of the same e()
        capture matrix drop tct_ED
        capture quietly ereturn display, level(`lv')
        capture matrix tct_ED = r(table)
        tct_run dwlsmv, level(`lv') thresholds
        mata: st_local("f", strofreal(tct_vsed("e(V)", "2 `tag' vs ereturn display, level `lv'")))
        if `f' < 0 {
            display as text "  2 `tag' level `lv': ereturn display posted no r(table); comparator skipped (diagnostic)"
        }
        else tct_ok "`f' == 0" "2 `tag': r(table) against ereturn display at level `lv', `f' failures"
    }
    set linesize 120
end

import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
local mA "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6)"
local mE "(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)"

* ============================================================================================
* 1  display only
* ============================================================================================
quietly dwlsmv `mA'
mata: TCT_S0 = tct_snap()
tempname b0 V0
matrix `b0' = e(b)
matrix `V0' = e(V)
foreach opts in "" "level(90)" "standardized" "thresholds" "nothresholds" "standardized level(99) thresholds" {
    tct_run dwlsmv, `opts'
    tct_ok "r(rc) == 0" "1: replay [`opts'] rc"
    mata: st_local("same", strofreal(tct_snap() == TCT_S0))
    tct_ok "`same'" "1: e() is not the same after replay [`opts']"
}
tct_ok "mreldif(e(b), `b0') == 0 & mreldif(e(V), `V0') == 0" "1: e(b) and e(V) after the replays"
capture lincom _b[y2:F1]
tct_ok "_rc == 0" "1: lincom still runs after the display"
capture test _b[y2:F1] = _b[y3:F1]
tct_ok "_rc == 0" "1: test still runs after the display"

* ============================================================================================
* 2  numbers: levels 90, 95, 99, on every fit shape
* ============================================================================================
quietly dwlsmv `mA'
tct_numbers "marker CFA" 1
quietly dwlsmv `mA', metric(unitvariance)
tct_numbers "unit-variance CFA" 1
quietly dwlsmv `mE'
tct_numbers "structural, marker" 1
quietly dwlsmv `mE', metric(unitvariance)
tct_numbers "structural, unit variance" 1
quietly dwlsmv (F1 -> y1@1 y2@.8 y3) (F2 -> y4 y5 y6), cov(F1*F2@.2)
tct_numbers "user constraints" 1

* the standardized table: r(table) against dwlsmv_m11_std's own matrices
quietly dwlsmv `mE'
tempname bs Vs Cs
mata: dwlsmv_m11_std("`bs'", "`Vs'", "`Cs'")
foreach lv in 90 95 99 {
    tct_run dwlsmv, standardized level(`lv') thresholds
    tct_ok "r(rc) == 0 & r(hasrt) == 1" "2 standardized: level `lv' rc and r(table)"
    mata: st_local("f", strofreal(tct_numeq(`lv', "`bs'", "`Vs'", "2 standardized `lv'")))
    tct_ok "`f' == 0" "2 standardized: r(table) at level `lv', `f' failures"
}

* G = 2, every invariance level
preserve
import delimited using "fixtures/data/mg_E.csv", clear
local mG "(F -> y1 y2 y3 y4 y5)"
foreach lev in configural thresholds metric scalar {
    quietly dwlsmv `mG', group(g) invariance(`lev') metric(unitvariance)
    tct_numbers "G = 2 `lev'" 0
}
restore

* a covariate fit
preserve
set seed 20261114
clear
quietly set obs 1500
quietly generate byte x = runiform() < .5
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
drop eta z1 z2 z3 z4
quietly dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- x)
tct_numbers "covariate" 1
restore

* ============================================================================================
* 3  row-set accounting
* ============================================================================================
quietly dwlsmv `mA'
tempname b3 V3 C3
matrix `b3' = e(b)
matrix `V3' = e(V)
matrix `C3' = e(Cns)
mata: tct_plan("`b3'", "`V3'", "`C3'", 0, 0)
mata: st_local("nt", strofreal(tct_nthr()))
mata: st_local("f", strofreal(tct_acct("`b3'", "`V3'", "`C3'", "thresholds", 2, 0, `nt', "3 marker CFA, thresholds hidden")))
tct_ok "`f' == 0" "3: marker CFA, thresholds hidden: `f' failures"
mata: tct_plan("`b3'", "`V3'", "`C3'", 1, 0)
mata: st_local("f", strofreal(tct_acct("`b3'", "`V3'", "`C3'", "thresholds", 2, 0, 0, "3 marker CFA, thresholds shown")))
tct_ok "`f' == 0" "3: marker CFA, thresholds shown, nothing omitted: `f' failures"

quietly dwlsmv `mA', metric(unitvariance)
matrix `b3' = e(b)
matrix `V3' = e(V)
matrix `C3' = e(Cns)
mata: tct_plan("`b3'", "`V3'", "`C3'", 0, 0)
mata: st_local("nt", strofreal(tct_nthr()))
mata: st_local("f", strofreal(tct_acct("`b3'", "`V3'", "`C3'", "thresholds metric", 0, 0, `nt' + 2, "3 unit variance")))
tct_ok "`f' == 0" "3: unit-variance CFA: thresholds and the two metric fixings omitted: `f' failures"
mata: st_local("nm", strofreal(tct_nreason("metric")))
tct_ok "`nm' == 2" "3: unit-variance CFA: two rows omitted as metric fixings (`nm')"

quietly dwlsmv `mE', metric(unitvariance)
matrix `b3' = e(b)
matrix `V3' = e(V)
matrix `C3' = e(Cns)
mata: tct_plan("`b3'", "`V3'", "`C3'", 1, 0)
mata: st_local("nm", strofreal(tct_nreason("metric")))
tct_ok "`nm' == 2" "3: structural, unit variance: var(F1) and var(e.F2) omitted as metric fixings (`nm')"
mata: st_local("f", strofreal(tct_acct("`b3'", "`V3'", "`C3'", "metric", 0, 0, 2, "3 structural unit variance")))
tct_ok "`f' == 0" "3: structural, unit variance: `f' failures"

quietly dwlsmv (F1 -> y1@1 y2@.8 y3) (F2 -> y4 y5 y6), cov(F1*F2@.2)
matrix `b3' = e(b)
matrix `V3' = e(V)
matrix `C3' = e(Cns)
mata: tct_plan("`b3'", "`V3'", "`C3'", 1, 0)
mata: st_local("f", strofreal(tct_acct("`b3'", "`V3'", "`C3'", "", 2, 2, 0, "3 user constraints")))
tct_ok "`f' == 0" "3: user constraints: two markers, two constrained, none omitted: `f' failures"

* standardized, G = 1: the constant exogenous variances are omitted (C2 (a), owner 2026-11-18)
quietly dwlsmv `mA'
mata: dwlsmv_m11_std("`bs'", "`Vs'", "`Cs'")
mata: tct_plan("`bs'", "`Vs'", "`Cs'", 1, 1)
mata: st_local("f", strofreal(tct_acct("`bs'", "`Vs'", "`Cs'", "standardized", 0, 0, 2, "3 standardized G = 1")))
tct_ok "`f' == 0" "3: standardized G = 1: var(F1) and var(F2) are omitted, nothing constrained: `f' failures"

* G = 2: the omissions are fixed rows, the level conventions
preserve
import delimited using "fixtures/data/mg_E.csv", clear
foreach lev in configural thresholds metric scalar {
    quietly dwlsmv `mG', group(g) invariance(`lev') metric(unitvariance)
    matrix `b3' = e(b)
    matrix `V3' = e(V)
    matrix `C3' = e(Cns)
    mata: tct_plan("`b3'", "`V3'", "`C3'", 1, 0)
    mata: st_local("f", strofreal(tct_acct("`b3'", "`V3'", "`C3'", "metric level", 0, 0, `=rowsof(`C3')', "3 G = 2 `lev'")))
    tct_ok "`f' == 0" "3: G = 2 `lev': every fixed row is omitted as a convention, nothing constrained: `f' failures"
    mata: tct_plan("`b3'", "`V3'", "`C3'", 0, 0)
    mata: st_local("nt", strofreal(tct_nthr()))
    mata: st_local("f", strofreal(tct_acct("`b3'", "`V3'", "`C3'", "metric level thresholds", 0, 0, `=rowsof(`C3')' + `nt', "3 G = 2 `lev' hidden")))
    tct_ok "`f' == 0" "3: G = 2 `lev', thresholds hidden: `f' failures"
    mata: dwlsmv_m11_std_g("`bs'", "`Vs'", "`Cs'")
    mata: tct_plan("`bs'", "`Vs'", "`Cs'", 1, 1)
    mata: st_local("f", strofreal(tct_acct("`bs'", "`Vs'", "`Cs'", "standardized", 0, 0, `=rowsof(`Cs')', "3 G = 2 `lev' standardized")))
    tct_ok "`f' == 0" "3: G = 2 `lev' standardized: every constant row omitted, nothing constrained: `f' failures"
    if inlist("`lev'", "metric", "scalar") {
        mata: tct_plan("`b3'", "`V3'", "`C3'", 1, 0)
        mata: st_local("f", strofreal(tct_tied()))
        tct_ok "`f' == 0" "3: G = 2 `lev': a tied loading is repeated in both groups with the same estimate and se: `f' differ"
    }
}
restore

* ============================================================================================
* 4  widths
* ============================================================================================
preserve
rename y1 abcdefghijabcdefghijabcdefgh
foreach ls in 78 80 120 {
    set linesize `ls'
    quietly dwlsmv (F1 -> abcdefghijabcdefghijabcdefgh y2 y3) (F2 -> y4 y5 y6)
    tct_run dwlsmv
    set linesize 120
    tct_ok "r(rc) == 0" "4: replay at linesize `ls'"
    * the label column: max(16, longest + 2) up to linesize - 62; the longest is the indented 28-character name, 30
    local w = min(max(16, 30 + 2), max(`ls' - 62, 16))
    mata: st_local("wd", strofreal(tct_width(`ls')[1]))
    mata: st_local("fits", strofreal(tct_width(`ls')[2]))
    tct_ok "`wd' == `w' + 62" "4: linesize `ls': the table is `wd' wide, wanted `=`w' + 62'"
    tct_ok "`fits'" "4: linesize `ls': no table line is wider than the linesize"
    if `ls' < 94 {
        * the label column w-1 wide holds the two-space indent, the first w-4 characters of the name and the tilde
        local cutp = substr("abcdefghijabcdefghijabcdefgh", 1, `w' - 4) + "~"
        tct_chk 1 "`cutp'" "4: linesize `ls': the long label is cut with a tilde, [`cutp']"
        tct_chk 0 "abcdefghijabcdefghijabcdefgh" "4: linesize `ls': and not printed in full"
    }
    else {
        tct_chk 1 "abcdefghijabcdefghijabcdefgh" "4: linesize `ls': the long label is printed in full"
        tct_chk 0 "~" "4: linesize `ls': and not cut"
    }
}
set linesize 120
restore

* ============================================================================================
* 5  the layouts
* ============================================================================================
* ---- (a) marker CFA ----
quietly dwlsmv `mA'
tct_run dwlsmv
foreach p in "Loadings on F1" "Loadings on F2" "Variances" "var(F1)" "var(F2)" "Covariances" "cov(F1,F2)" "Estimate" ///
    "Std. err." "P>|z|" "[95% conf. interval]" "Thresholds are not shown" "Metric: marker" {
    tct_chk 1 "`p'" "5a: marker CFA prints [`p']"
}
foreach p in "Structural" "Covariates" "cut1" "Thresholds: y1" "Coefficient" "Fixed" "(constrained)" "(omitted)" "Group" {
    tct_chk 0 "`p'" "5a: marker CFA does not print [`p']"
}
tct_n 2 "1  (marker)" "5a: two markers"
tct_ord "5a: section order" "Loadings on F1" "Loadings on F2" "Variances" "Covariances" "Thresholds are not shown"

* ---- (b) unit-variance CFA ----
quietly dwlsmv `mA', metric(unitvariance)
tct_run dwlsmv
foreach p in "Loadings on F1" "Loadings on F2" "Covariances" "cov(F1,F2)" "Metric: factor variances fixed at 1" {
    tct_chk 1 "`p'" "5b: unit-variance CFA prints [`p']"
}
foreach p in "Variances" "var(F1)" "var(F2)" "(marker)" "(constrained)" {
    tct_chk 0 "`p'" "5b: unit-variance CFA does not print [`p']"
}
tct_chk 0 "Fixed by the level" "5b: no convention line at one group (the Metric line says it)"

* ---- (c) structural, marker and unit variance, and the standardized table ----
quietly dwlsmv `mE'
tct_run dwlsmv
foreach p in "Loadings on F1" "Loadings on F2" "Structural" "F2 <- F1" "Variances" "var(F1)" "var(e.F2)" {
    tct_chk 1 "`p'" "5c: structural marker fit prints [`p']"
}
tct_chk 0 "Covariances" "5c: structural marker fit has no Covariances section"
tct_n 2 "1  (marker)" "5c: two markers"
tct_ord "5c: section order" "Loadings on F1" "Loadings on F2" "Structural" "Variances" "Thresholds are not shown"
quietly dwlsmv `mE', metric(unitvariance)
tct_run dwlsmv
tct_chk 1 "Structural" "5c: unit-variance structural fit prints Structural"
foreach p in "Variances" "var(F1)" "var(e.F2)" "(marker)" "(constrained)" {
    tct_chk 0 "`p'" "5c: unit-variance structural fit does not print [`p']"
}
quietly dwlsmv `mE'
tct_run dwlsmv, standardized
foreach p in "Standardized (std.all); /var(e.y) rows are the derived residual variances" "var(e.y1)" "var(e.y6)" "var(e.F2)" ///
    "Thresholds are not shown" {
    tct_chk 1 "`p'" "5c: the standardized table prints [`p']"
}
tct_chk 1 "Constant rows (fixed by the level or by the standardization) are not shown." "5c: the omission clause at one group (C2 (a))"
tct_n 0 "(constrained)" "5c: standardized: var(F1) of the one exogenous factor is omitted, nothing constrained"
tct_chk 0 "var(F1)" "5c: standardized: no var(F1) row"
tct_chk 0 "(marker)" "5c: standardized: the marker has a standard error and is no marker row"

* ---- (d) user constraints ----
quietly dwlsmv (F1 -> y1@1 y2@.8 y3) (F2 -> y4 y5 y6), cov(F1*F2@.2)
tct_run dwlsmv
tct_chk 1 "0.8  (constrained)" "5d: a fixed loading shows its value and (constrained)"
tct_chk 1 "0.2  (constrained)" "5d: a fixed covariance shows its value and (constrained), never (omitted)"
tct_n 2 "(constrained)" "5d: exactly two constrained rows"
tct_n 2 "1  (marker)" "5d: and the two markers"
tct_chk 0 "(omitted)" "5d: no (omitted)"
* a user-fixed factor variance is not a metric fixing under the marker metric: shown with its value
quietly dwlsmv (F1 -> y1@1 y2 y3) (F2 -> y4 y5 y6), var(F1@.5)
tct_run dwlsmv
tct_chk 1 "0.5  (constrained)" "5d: a user-fixed variance shows its value and (constrained)"
tct_chk 1 "var(F1)" "5d: in the Variances section"
tct_n 2 "1  (marker)" "5d: the markers stay markers"

* ---- (e) G = 2 at four invariance levels, and standardized ----
preserve
import delimited using "fixtures/data/mg_E.csv", clear
foreach lev in configural thresholds metric scalar {
    quietly dwlsmv `mG', group(g) invariance(`lev') metric(unitvariance)
    tct_run dwlsmv
    tct_ok "r(rc) == 0" "5e `lev': replay rc"
    foreach p in "Group 1 (g = 1)" "Group 2 (g = 2)" "Loadings on F" "Estimate" "Std. err." {
        tct_chk 1 "`p'" "5e `lev': prints [`p']"
    }
    tct_n 2 "Loadings on F" "5e `lev': one loadings block per group"
    tct_n 1 "Fixed by the level, not shown: variances 1, intercepts 0, scales 1, means 0." "5e `lev': the convention line, once"
    foreach p in "(constrained)" "(omitted)" "(marker)" "Group 2 (g = 1)" {
        tct_chk 0 "`p'" "5e `lev': does not print [`p']"
    }
    local wantvar = inlist("`lev'", "metric", "scalar")
    tct_chk `wantvar' "Variances" "5e `lev': the Variances section is there only where group 2's variance is free"
    if "`lev'" == "configural" {
        tct_chk 0 "cut1" "5e `lev': thresholds hidden by default"
        tct_chk 1 "Thresholds are not shown" "5e `lev': and the note says so"
    }
    else {
        foreach p in "Thresholds: y1" "Thresholds: y5" "cut1 (Group 1)" "cut1 (Group 2)" {
            tct_chk 1 "`p'" "5e `lev': thresholds shown by default, [`p']"
        }
        tct_chk 0 "Thresholds are not shown" "5e `lev': no hidden note when shown"
        tct_ord "5e `lev': thresholds come after the group tables" "Group 1 (g = 1)" "Group 2 (g = 2)" "Thresholds: y1" "Thresholds: y5"
    }
    if "`lev'" == "scalar" {
        tct_chk 1 "Factor means" "5e `lev': the free factor mean of group 2 is a row"
    }
    * the table of G = 2 is the same on replay as the notes-block group form is different (M9): no colon after the title
    tct_chk 0 "Group 2 (g = 2):" "5e `lev': the group title has no colon (the notes-block form)"
}
quietly dwlsmv `mG', group(g) invariance(metric) metric(unitvariance)
tct_run dwlsmv, standardized
foreach p in "Standardized (std.all); /var(e.y) rows are the derived residual variances" ///
    "Constant rows (fixed by the level or by the standardization) are not shown." "Note: equal raw parameters can have different standardized values across groups." ///
    "Group 1 (g = 1)" "var(e.y1)" {
    tct_chk 1 "`p'" "5e standardized G = 2: prints [`p']"
}
tct_n 0 "(constrained)" "5e standardized G = 2: no constant row is shown"
tct_n 0 "Fixed by the level" "5e standardized G = 2: the title clause says it, not the convention line"
restore

* ---- (f) a covariate ----
preserve
set seed 20261114
clear
quietly set obs 1500
quietly generate byte x = runiform() < .5
quietly generate double eta = .4 * x + rnormal()
forvalues j = 1/4 {
    quietly generate double z`j' = .7 * eta + sqrt(.51) * rnormal() + (`j' == 4) * .3 * x
    quietly generate byte y`j' = (z`j' > -.5) + (z`j' > .5)
}
drop eta z1 z2 z3 z4
quietly dwlsmv (F -> y1 y2 y3 y4) (F <- x) (y4 <- x)
tct_run dwlsmv
foreach p in "Loadings on F" "Covariates" "F <- x" "y4 <- x" "Variances" "var(e.F)" "Thresholds are not shown" {
    tct_chk 1 "`p'" "5f: covariate fit prints [`p']"
}
tct_n 1 "1  (marker)" "5f: one marker"
tct_ord "5f: section order" "Loadings on F" "Covariates" "Variances" "Thresholds are not shown"
quietly dwlsmv (F -> y1 y2 y3 y4) (F <- x@.4) (y4 <- x@-.3)
tct_run dwlsmv
tct_chk 1 "0.4  (constrained)" "5f: a fixed gamma shows its value and (constrained)"
tct_chk 1 "-0.3  (constrained)" "5f: a fixed kappa shows its value and (constrained)"
restore

* ---- (g) a fit that did not converge: real (the marker trap), and hand-set at G = 2 ----
preserve
import delimited using "fixtures/data/marker_trap.csv", clear varnames(1) delimiters(",") stringcols(_all)
quietly destring q* att*, replace
quietly dwlsmv (MathAb -> q1@1 q2@1 q3-q8) (MathAtt -> att1-att5) (MathAtt -> MathAb)
tct_ok "e(converged) == 0" "5g: PREMISE: the marker-trap fit does not converge (the data no longer provoke the case)"
tct_run dwlsmv
foreach p in "Last iterate, not a solution; no standard errors" "Loadings on MathAb" "Loadings on MathAtt" "Structural" "MathAb <- MathAtt" "Estimate" {
    tct_chk 1 "`p'" "5g: last iterate prints [`p']"
}
foreach p in "Std. err." "P>|z|" "conf. interval" "Fixed" "Prob > chi2" {
    tct_chk 0 "`p'" "5g: last iterate does not print [`p']"
}
* q1 is the marker; q2 is a user loading fixed at 1 and is NOT the marker (M5); att1 is the second factor's marker
tct_n 2 "1  (marker)" "5g: the markers named by e(canon): q1 and att1"
tct_n 1 "1  (constrained)" "5g: q2@1 is constrained, not a marker"
tct_run dwlsmv, thresholds
tct_chk 1 "Thresholds: q1" "5g: thresholds shows the item blocks of the last iterate"
tct_chk 0 "Thresholds are not shown" "5g: and no hidden note"
mata: st_local("hasrt", strofreal(rows(st_matrix("tct_RT"))))
tct_ok "`hasrt' == 9" "5g: r(table) of a last-iterate table has its nine rows"
restore
preserve
import delimited using "fixtures/data/case_E.csv", clear varnames(1) delimiters(",")
generate g = mod(_n, 2)
quietly dwlsmv (F1 -> y1-y3) (F2 -> y4-y6), group(g) invariance(metric) metric(unitvariance)
tct_unconv
tct_run dwlsmv
foreach p in "Last iterate, not a solution; no standard errors" "Group 1 (g = 0)" "Group 2 (g = 1)" "Loadings on F1" "Estimate" ///
    "Fixed by the level, not shown" "Thresholds: y1" {
    tct_chk 1 "`p'" "5g: hand-set G = 2 prints [`p']"
}
foreach p in "Std. err." "P>|z|" "(constrained)" {
    tct_chk 0 "`p'" "5g: hand-set G = 2 does not print [`p']"
}
restore

* ============================================================================================
* 6  threshold defaults
* ============================================================================================
quietly dwlsmv `mA'
tct_run dwlsmv
tct_chk 0 "cut1" "6: G = 1 hides the thresholds by default"
tct_chk 1 "Thresholds are not shown; the option thresholds shows them (they are in e(b))." "6: and the note says so, in its words"
tct_run dwlsmv, thresholds
foreach p in "Thresholds: y1" "Thresholds: y6" "cut1" "cut2" {
    tct_chk 1 "`p'" "6: thresholds shows [`p']"
}
tct_chk 0 "Thresholds are not shown" "6: and says nothing"
tct_ord "6: thresholds last" "Loadings on F2" "Variances" "Covariances" "Thresholds: y1" "Thresholds: y6"
tct_run dwlsmv, nothresholds
tct_chk 0 "cut1" "6: nothresholds hides them"
tct_run dwlsmv `mA', thresholds
tct_chk 1 "Thresholds: y1" "6: thresholds at estimation"
tct_run dwlsmv `mA', nothresholds
tct_chk 1 "Thresholds are not shown" "6: nothresholds at estimation"
preserve
import delimited using "fixtures/data/mg_E.csv", clear
foreach lev in configural thresholds metric scalar {
    quietly dwlsmv `mG', group(g) invariance(`lev') metric(unitvariance)
    local want = ("`lev'" != "configural")
    tct_run dwlsmv
    tct_chk `want' "Thresholds: y1" "6: G = 2 `lev' default"
    tct_chk `=1 - `want'' "Thresholds are not shown" "6: G = 2 `lev' default note"
    tct_run dwlsmv, nothresholds
    tct_chk 0 "Thresholds: y1" "6: G = 2 `lev' nothresholds"
    tct_chk 1 "Thresholds are not shown" "6: G = 2 `lev' nothresholds note"
    tct_run dwlsmv, thresholds
    tct_chk 1 "Thresholds: y1" "6: G = 2 `lev' thresholds"
    tct_chk 0 "Thresholds are not shown" "6: G = 2 `lev' thresholds note"
}
restore

* ============================================================================================
* 7  error paths
* ============================================================================================
quietly dwlsmv `mA'
tct_run dwlsmv, thresholds nothresholds
tct_ok "r(rc) == 198" "7: thresholds with nothresholds is rc 198 exactly"
tct_run dwlsmv, se
tct_ok "r(rc) == 198" "7: se without efa() is rc 198 exactly"
tct_run dwlsmv, ci
tct_ok "r(rc) == 198" "7: ci without efa() is rc 198 exactly"
tct_run dwlsmv `mA', thresholds nothresholds
tct_ok "r(rc) == 198" "7: thresholds with nothresholds at estimation is rc 198 exactly"
tct_unconv
tct_run dwlsmv, standardized
tct_ok "r(rc) == 198" "7: standardized on a fit that did not converge is rc 198 exactly"
tct_chk 1 "standardized is not available for a fit that did not converge: the estimates are" "7: and says why"
tct_chk 0 "Estimate" "7: and prints no table"
quietly dwlsmv `mA'
local msg1 "not in e(b)'s grammar"
local msg2 "appears twice"
local msg3 "e(V) does not fit e(b)"
local msg4 "e(b) has a missing value"
local msg5 "does not select exactly one parameter"
local msg6 "does not select exactly one parameter"
local msg7 "is fixed twice"
local msg8 "fixed-value rows do not fit e(b)"
local msg9 "confidence level is not between 0 and 100"
local msg10 "an unreadable row"
local msg11 "column stripe does not fit e(b)"
local msg12 "e(V) has a missing diagonal entry"
forvalues c = 1/12 {
    tct_run mata: tct_bad(`c')
    tct_ok "r(rc) == 3498" "7.`c': malformed planner input is rc 3498 exactly"
    tct_chk 1 "dwlsmv: internal error, coefficient table input is malformed" "7.`c': the short message"
    tct_chk 1 "`msg`c''" "7.`c': and errprintf names what (`msg`c'')"
}
preserve
import delimited using "fixtures/data/mg_E.csv", clear
quietly dwlsmv `mG', group(g) invariance(metric) metric(unitvariance)
local gmsg1 "a loading without its group tag or factor"
local gmsg2 "group values without the group variable"
local gmsg3 "an empty group value"
forvalues c = 1/3 {
    tct_run mata: tct_badg(`c')
    tct_ok "r(rc) == 3498" "7g.`c': malformed G = 2 planner input is rc 3498 exactly"
    tct_chk 1 "`gmsg`c''" "7g.`c': and says what (`gmsg`c'')"
}
restore

* ============================================================================================
* 8  replay = estimation
* ============================================================================================
quietly dwlsmv `mA'
tct_run dwlsmv `mA'
mata: TCT_T1 = tct_segment()
tct_run dwlsmv
mata: st_local("same", strofreal(TCT_T1 != "" & TCT_T1 == tct_segment()))
tct_ok "`same'" "8: raw table at estimation = replay"
tct_run dwlsmv `mA', standardized
mata: TCT_T1 = tct_segment()
tct_run dwlsmv, standardized
mata: st_local("same", strofreal(TCT_T1 != "" & TCT_T1 == tct_segment()))
tct_ok "`same'" "8: standardized table at estimation = replay"
preserve
import delimited using "fixtures/data/mg_E.csv", clear
tct_run dwlsmv `mG', group(g) invariance(configural) metric(unitvariance)
mata: TCT_T1 = tct_segment()
tct_run dwlsmv
mata: st_local("same", strofreal(TCT_T1 != "" & TCT_T1 == tct_segment()))
tct_ok "`same'" "8: G = 2 configural table at estimation = replay"
tct_run dwlsmv `mG', group(g) invariance(metric) metric(unitvariance) nothresholds
mata: TCT_T1 = tct_segment()
tct_run dwlsmv, nothresholds
mata: st_local("same", strofreal(TCT_T1 != "" & TCT_T1 == tct_segment()))
tct_ok "`same'" "8: G = 2 metric, nothresholds, table at estimation = replay"
restore
* r(table) at estimation: a diagnostic only (the notes block runs after the table; whether r() survives it is not
* what this file asserts)
quietly dwlsmv `mA'
tct_run dwlsmv `mA'
display as text "  8 diagnostic: r(table) after an estimation call: r(hasrt) = " r(hasrt)

adopath - "`here'/src/ado"
adopath - "`here'/build"
if $TCT_NF {
    display as error _n "test_ctab: $TCT_NF failure(s)"
    exit 9
}
display as result "test_ctab: the compact table's numbers, row accounting, widths, layouts, threshold defaults, error paths and replay hold"
exit 0
