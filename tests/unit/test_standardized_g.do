*! Module 11: the standardized solution at G > 1 (spec 4.8 O17, derivation
*! @sec-standardize-mg, tier 2 fixtures/expected/lavaan_mg_std*.csv).
*!
*! Executed under StataMP 18 on 2026-11-01 and it passes. Mutants: E-074, 2026-11-01.
*!
*! Two comparisons with lavaan 0.7.2, kept apart because their errors differ:
*!   1  the map alone: lavaan's own theta-hat and vcov go through
*!      dwlsmv_m11_std_gc, so no optimizer difference enters. Points 1e-12
*!      absolute, SEs 1e-7 relative, both carried from std-mg-frozen; a row
*!      lavaan shows with SE 0 must be below 1e-10 here (std-mg-frozen's rule).
*!   2  end to end: the command's own fit against lavaan's solution, after
*!      per-group sign equivalence. Points at the theta tier (1e-5 abs + 1e-4
*!      rel), SEs at the SE tier (1e-8 abs + 1e-5 rel): the tiers of
*!      test_invariance_post.do check 5, not re-chosen. Any revision needs its
*!      own measurement, per quantity, frozen before it is evaluated.
*! Independent of lavaan:
*!   3  Vs against J V J', J from central differences of the map, h = 1e-6:
*!      1e-6 relative on nonzero diagonals (chosen before the first run);
*!   4  lambda*^2 + Theta* = 1 on one-factor designs, 1e-12; group 1's tau and
*!      lambda rows equal the raw estimate exactly (omega = 1, psi = 1);
*!   5  display: constant rows kept (Q-D (a)); Stata prints a zero-valued one as
*!      "0 (omitted)", the rest as (constrained), as the raw table; the tie note
*!      (Q-C (a)) past configural only, e() restored, replay the same;
*!   6  groups 0 and 7 with shuffled rows, and G = 3;
*!   7  fail closed: each malformed input gives rc 3498 and its own message.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 255

quietly do tools/build_mlib.do
quietly do tests/dwlsmv_assert.mata
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

foreach f in std std_theta std_vcov {
    import delimited using "fixtures/expected/lavaan_mg_`f'.csv", ///
        clear varnames(1) stringcols(_all) delimiters(",")
    mata: F_`f' = st_sdata(., .)
}

mata:
// the fixture's key of an output row, from spec 4.8's name rule
string scalar tsg_key(string scalar eq, string scalar nm, string rowvector ind,
                      string rowvector fac, string rowvector gv)
{
    real scalar p, a
    string scalar h, t, x, y, gs

    p = strpos(nm, "#")
    h = (p ? substr(nm, 1, p - 1) : "")
    t = (p ? substr(nm, p + 1, .) : nm)
    if (substr(h, 2, .) != "" & strpos(h, ".") & !strpos(h, "(")) {
        x = h
        h = t
        t = x
    }
    gs = strofreal(selectindex(gv :== subinstr(substr(t, 1, strpos(t, ".") - 1), "bn", "")))
    if (substr(eq, 1, 1) == "/" & eq != "/") {
        return("4:" + gs + ":" + strofreal(selectindex(ind :== substr(eq, 2, .))) + ":" + substr(h, 6, .))
    }
    if (eq != "/" & h == "") {
        return("5:" + gs + ":" + strofreal(selectindex(ind :== eq)) + ":1")
    }
    if (eq != "/") {
        return("1:" + gs + ":" + strofreal(selectindex(ind :== eq)) + ":" +
               strofreal(selectindex(fac :== substr(h, 3, .))))
    }
    x = substr(h, strpos(h, "(") + 1, strlen(h) - strpos(h, "(") - 1)
    if (substr(h, 1, 6) == "scale(") {
        return("6:" + gs + ":" + strofreal(selectindex(ind :== x)) + ":1")
    }
    if (substr(h, 1, 5) == "mean(") {
        return("7:" + gs + ":" + strofreal(selectindex(fac :== x)) + ":1")
    }
    if (substr(h, 1, 4) == "cov(") {
        a = strpos(x, ",")
        y = substr(x, a + 1, .)
        x = substr(x, 1, a - 1)
        return("2:" + gs + ":" + strofreal(max((selectindex(fac :== x), selectindex(fac :== y)))) +
               ":" + strofreal(min((selectindex(fac :== x), selectindex(fac :== y)))))
    }
    if (substr(x, 1, 2) == "e.") {
        return("8:" + gs + ":" + strofreal(selectindex(ind :== substr(x, 3, .))) + ":1")
    }
    return("2:" + gs + ":" + strofreal(selectindex(fac :== x)) + ":" + strofreal(selectindex(fac :== x)))
}

string rowvector tsg_keys(string matrix st, string rowvector ind, string rowvector fac)
{
    string rowvector gv, k
    string matrix K_s
    real scalar i

    K_s = st_matrixrowstripe("e(K)")
    gv  = K_s[., 2]'
    k   = J(1, rows(st), "")
    for (i = 1; i <= rows(st); i++) {
        k[i] = tsg_key(st[i, 1], st[i, 2], ind, fac, gv)
    }
    return(k)
}

// the map's value with row j of the estimate moved by h
real rowvector tsg_g(real rowvector bp, real scalar j, real scalar h, string matrix cs,
                     string rowvector ind, real matrix K, string rowvector gv, string scalar gvar)
{
    real rowvector b, g
    real matrix Vs
    string matrix out

    b    = bp
    b[j] = bp[j] + h
    // no ties, a free identity V: only the map is under test
    dwlsmv_m11_std_gc(b, I(cols(b)), cs, cs, ind, K, gv, gvar, "", cols(b), g, Vs, out)
    return(g)
}

// each tie class exact from its first member, as e(V) shares theta
void tsg_share(real rowvector b, real matrix Vl, real matrix V)
{
    real scalar i, k

    for (i = 2; i <= cols(b); i++) {
        for (k = 1; k < i; k++) {
            if (V[i, i] > 0) {
                if (V[i, .] == V[k, .]) {
                    if (abs(b[i] - b[k]) > 1e-6) {
                        _error("test: tied lavaan values differ")
                    }
                    b[i]     = b[k]
                    Vl[i, .] = Vl[k, .]
                    Vl[., i] = Vl[., k]
                    break
                }
            }
        }
    }
}

// checks 1, 3, 4 on the current e(); returns failures
real scalar tsg_map(string scalar dn, string scalar lv, string rowvector fac)
{
    external string matrix F_std, F_std_theta, F_std_vcov
    real rowvector b, g, bp
    real matrix V, Vs, K, Vl, Jn, Vfd
    string matrix cs, out, K_s
    string rowvector ind, gv, kb, ko, kth
    string colvector keys
    real colvector d
    real scalar i, j, m, q, h, se, ptmax, semax, zmax, hit, rank, gr
    string scalar canon, gvar, st

    b     = st_matrix("e(b)")
    V     = st_matrix("e(V)")
    cs    = st_matrixcolstripe("e(b)")
    K     = st_matrix("e(K)")
    K_s   = st_matrixrowstripe("e(K)")
    gv    = K_s[., 2]'
    ind   = tokens(st_global("e(indicators)"))
    canon = st_global("e(canon)")
    gvar  = st_global("e(groupvar)")
    rank  = st_numscalar("e(rank)")
    q     = cols(b)
    kb    = tsg_keys(cs, ind, fac)
    keys  = kb'

    d   = selectindex(F_std_theta[., 1] :== dn :& F_std_theta[., 2] :== lv)
    kth = F_std_theta[d, 3]'
    bp  = b
    for (i = 1; i <= cols(kth); i++) {
        j = selectindex(kb :== kth[i])
        if (cols(j) * rows(j) != 1) {
            printf("{err}  %s %s: lavaan key %s not in e(b)\n", dn, lv, kth[i])
            return(1)
        }
        bp[j] = strtoreal(F_std_theta[d[i], 4])
    }
    d  = selectindex(F_std_vcov[., 1] :== dn :& F_std_vcov[., 2] :== lv)
    Vl = J(q, q, 0)
    for (i = 1; i <= rows(d); i++) {
        Vl[selectindex(keys :== F_std_vcov[d[i], 3]), selectindex(keys :== F_std_vcov[d[i], 4])] =
            strtoreal(F_std_vcov[d[i], 5])
    }
    Vl = (Vl + Vl') / 2
    tsg_share(bp, Vl, V)

    dwlsmv_m11_std_gc(bp, Vl, cs, cs, ind, K, gv, gvar, canon, rank, g, Vs, out)
    ko = tsg_keys(out, ind, fac)

    // check 1
    d = selectindex(F_std[., 1] :== dn :& F_std[., 2] :== lv)
    ptmax = 0
    semax = 0
    zmax  = 0
    hit   = 0
    for (i = 1; i <= rows(d); i++) {
        j = selectindex(ko :== F_std[d[i], 3])
        if (cols(j) * rows(j) != 1) {
            printf("{err}  1 %s %s: fixture key %s not in the output\n", dn, lv, F_std[d[i], 3])
            return(1)
        }
        hit++
        ptmax = max((ptmax, abs(g[j] - strtoreal(F_std[d[i], 4]))))
        se = strtoreal(F_std[d[i], 5])
        if (se > 0) {
            semax = max((semax, abs(sqrt(Vs[j, j]) - se) / se))
        }
        else {
            zmax = max((zmax, sqrt(Vs[j, j])))
        }
    }
    printf("{txt}  1 %s %s: %f rows, points %9.2e, SEs %9.2e rel, SE-0 rows %9.2e\n",
           dn, lv, hit, ptmax, semax, zmax)
    if (hit != cols(ko) | ptmax > 1e-12 | semax > 1e-7 | zmax > 1e-10) {
        printf("{err}  1 %s %s FAILED\n", dn, lv)
        return(1)
    }

    // check 3
    h  = 1e-6
    m  = cols(g)
    Jn = J(m, q, 0)
    for (j = 1; j <= q; j++) {
        if (Vl[j, j] > 0) {
            Jn[., j] = (tsg_g(bp, j, h, cs, ind, K, gv, gvar) -
                        tsg_g(bp, j, -h, cs, ind, K, gv, gvar))' / (2 * h)
        }
    }
    Vfd = Jn * Vl * Jn'
    semax = 0
    for (i = 1; i <= m; i++) {
        if (Vs[i, i] > 1e-20) {
            semax = max((semax, abs(Vfd[i, i] - Vs[i, i]) / Vs[i, i]))
        }
    }
    printf("{txt}  3 %s %s: J V J' from differences, %9.2e rel\n", dn, lv, semax)
    if (semax > 1e-6) {
        printf("{err}  3 %s %s FAILED\n", dn, lv)
        return(1)
    }

    // check 4
    if (cols(fac) == 1) {
        semax = 0
        for (gr = 1; gr <= cols(gv); gr++) {
            for (j = 1; j <= cols(ind); j++) {
                st = strofreal(gr) + ":" + strofreal(j) + ":1"
                semax = max((semax, abs(g[selectindex(ko :== "1:" + st)]^2 +
                                        g[selectindex(ko :== "8:" + st)] - 1)))
            }
        }
        if (semax > 1e-12) {
            printf("{err}  4 %s %s FAILED, identity %g\n", dn, lv, semax)
            return(1)
        }
    }
    for (j = 1; j <= q; j++) {
        if (substr(ko[j], 3, 2) == "1:") {
            if (substr(ko[j], 1, 1) == "4" | substr(ko[j], 1, 1) == "1") {
                if (g[j] != bp[j]) {
                    printf("{err}  4 %s %s: group 1 row %s moved\n", dn, lv, ko[j])
                    return(1)
                }
            }
        }
    }
    return(0)
}

// check 2: the command's own solution against lavaan, by key
real scalar tsg_e2e(string scalar dn, string scalar lv, string rowvector fac,
                    string scalar bs, string scalar vs)
{
    external string matrix F_std
    real rowvector g, sgn
    real matrix Vs
    string rowvector ind, ko, w
    real colvector d
    real scalar i, j, se, est, ptmax, semax, nf, gr, kd

    ind = tokens(st_global("e(indicators)"))
    nf  = cols(fac)
    g   = st_matrix(bs)
    Vs  = st_matrix(vs)
    ko  = tsg_keys(st_matrixcolstripe(bs), ind, fac)
    d   = selectindex(F_std[., 1] :== dn :& F_std[., 2] :== lv)
    sgn = J(1, 2 * nf, 0)
    for (i = 1; i <= rows(d); i++) {
        w = tokens(subinstr(F_std[d[i], 3], ":", " "))
        if (w[1] == "1") {
            gr = strtoreal(w[2])
            j  = (gr - 1) * nf + strtoreal(w[4])
            if (sgn[j] == 0) {
                sgn[j] = sign(g[selectindex(ko :== F_std[d[i], 3])] * strtoreal(F_std[d[i], 4]))
            }
        }
    }
    ptmax = 0
    semax = 0
    for (i = 1; i <= rows(d); i++) {
        w   = tokens(subinstr(F_std[d[i], 3], ":", " "))
        kd  = strtoreal(w[1])
        gr  = strtoreal(w[2])
        j   = selectindex(ko :== F_std[d[i], 3])
        est = strtoreal(F_std[d[i], 4])
        if (kd == 1) {
            est = est * sgn[(gr - 1) * nf + strtoreal(w[4])]
        }
        else if (kd == 7) {
            est = est * sgn[(gr - 1) * nf + strtoreal(w[3])]
        }
        else if (kd == 2) {
            est = est * sgn[(gr - 1) * nf + strtoreal(w[3])] * sgn[(gr - 1) * nf + strtoreal(w[4])]
        }
        ptmax = max((ptmax, abs(g[j] - est) / (1e-5 + 1e-4 * abs(est))))
        se    = strtoreal(F_std[d[i], 5])
        semax = max((semax, abs(sqrt(Vs[j, j]) - se) / (1e-8 + 1e-5 * se)))
    }
    printf("{txt}  2 %s %s: points at %6.3f of the tier, SEs at %6.3f of the tier\n", dn, lv, ptmax, semax)
    if (ptmax > 1 | semax > 1) {
        printf("{err}  2 %s %s FAILED\n", dn, lv)
        return(1)
    }
    return(0)
}

// check 5: every omega row is a constraint at 1, and every constraint has SE 0
real scalar tsg_const(string scalar bs, string scalar vs, string scalar cs_n, real scalar nfix)
{
    real rowvector g
    real matrix Vs, C
    string matrix out
    real scalar i, nom, nz

    g   = st_matrix(bs)
    Vs  = st_matrix(vs)
    C   = st_matrix(cs_n)
    out = st_matrixcolstripe(bs)
    nom = 0
    for (i = 1; i <= cols(g); i++) {
        if (substr(out[i, 2], 1, 6) == "scale(") {
            nom++
            if (g[i] != 1 | Vs[i, i] != 0) {
                return(1)
            }
        }
    }
    nz = sum(diagonal(Vs) :== 0)
    return(nom == 0 | nz != nfix | rows(C) != nfix)
}

// check 7: the i-th malformed input
void tsg_bad(real scalar c)
{
    real rowvector b, g
    real matrix V, Vs, K
    string matrix cs, rs, K_s, out
    string rowvector ind, gv
    string scalar canon, gvar
    real scalar rank, q, at, i

    b     = st_matrix("e(b)")
    V     = st_matrix("e(V)")
    cs    = st_matrixcolstripe("e(b)")
    rs    = cs
    K     = st_matrix("e(K)")
    K_s   = st_matrixrowstripe("e(K)")
    gv    = K_s[., 2]'
    ind   = tokens(st_global("e(indicators)"))
    canon = st_global("e(canon)")
    gvar  = st_global("e(groupvar)")
    rank  = st_numscalar("e(rank)")
    q     = cols(b)
    // a tied threshold row of group 2 and a lambda row of group 1
    at = 0
    for (i = 1; i <= q; i++) {
        if (cs[i, 1] == "/y1" & cs[i, 2] == "2.g#c.cut1") {
            at = i
        }
    }
    if (at == 0) {
        _error("test: no tied threshold row")
    }
    i = 0
    for (q = 1; q <= cols(b); q++) {
        if (cs[q, 1] == "y1" & cs[q, 2] == "1bn.g#c.F") {
            i = q
        }
    }
    q = cols(b)
    if (c == 1) {
        cs[i, 1] = "F"
        rs[i, 1] = "F"
    }
    else if (c == 2) {
        cs[i, 2] = "c.F"
        rs[i, 2] = "c.F"
    }
    else if (c == 3) {
        cs[i, 2] = "9.g#c.F"
        rs[i, 2] = "9.g#c.F"
    }
    else if (c == 4) {
        cs[at, .] = cs[at - 1, .]
        rs[at, .] = rs[at - 1, .]
    }
    else if (c == 5) {
        V = V[1..q - 1, 1..q - 1]
    }
    else if (c == 6) {
        rs[i, 2] = "2.g#c.F"
    }
    else if (c == 7) {
        b[at] = b[at] + 1e-9
    }
    else if (c == 8) {
        V[at, 1] = V[at, 1] + 1e-9
        V[1, at] = V[at, 1]
    }
    else if (c == 9) {
        rank = rank + 1
    }
    else if (c == 10) {
        b[i] = .
    }
    else if (c == 11) {
        gvar = ""
    }
    else if (c == 12) {
        K = K[1, .]
    }
    else if (c == 13) {
        b  = b[1, 1..at - 1], b[1, at + 1..q]
        V  = V[(1..at - 1, at + 1..q), (1..at - 1, at + 1..q)]
        cs = cs[(1..at - 1, at + 1..q), .]
        rs = cs
    }
    else if (c == 14) {
        canon = canon + char(10) + "tie tau 1 y1 9 tau 2 y1 9"
    }
    else if (c == 15) {
        ind[2] = ind[1]
    }
    else if (c == 16) {
        gv[2] = gv[1]
    }
    else if (c == 17) {
        canon = canon + char(10) + "tie tau 1"
    }
    else if (c == 18) {
        // a psi_ff at or below 0 in group 2
        for (i = 1; i <= q; i++) {
            if (cs[i, 2] == "var(F)#2.g") {
                b[i] = -.5
            }
        }
    }
    dwlsmv_m11_std_gc(b, V, cs, rs, ind, K, gv, gvar, canon, rank, g, Vs, out)
}
// the standardized table in a captured log: its title through its closing note
string scalar tsg_table(string scalar s)
{
    real scalar p, q
    string scalar r

    p = strpos(s, "Standardized (std.all)")
    if (!p) return("")
    r = substr(s, p, .)
    q = strpos(r, "equal raw parameters")
    if (!q) return("")
    return(substr(r, 1, q))
}
end

local nfail 0
local m1 (F -> y1 y2 y3 y4 y5)
local m2 (F -> y1 y2 y3) (H -> y4 y5)

capture program drop tsg_log
program tsg_log
    args cmd
    tempfile lg
    capture log close tsg
    quietly log using "`lg'", text replace name(tsg)
    capture noisily `cmd'
    local rc = _rc
    quietly log close tsg
    mata: TSGLOG = invtokens(cat(st_local("lg"))')
    c_local tsg_rc `rc'
end

* ---- checks 1 to 4, designs x levels --------------------------------------
local d1 E
local d2 B
local d3 D
local d4 B2F
forvalues di = 1/4 {
    local dn `d`di''
    local file = cond("`dn'" == "B2F", "B", "`dn'")
    import delimited using "fixtures/data/mg_`file'.csv", clear
    local model = cond("`dn'" == "B2F", "`m2'", "`m1'")
    local fac   = cond("`dn'" == "B2F", "F H", "F")
    foreach lv in configural thresholds metric scalar {
        local lav = cond("`lv'" == "metric", "loadings", cond("`lv'" == "scalar", "intercepts", "`lv'"))
        capture noisily quietly dwlsmv `model', group(g) invariance(`lv') metric(unitvariance)
        if _rc {
            display as error "`dn' `lv': the fit failed, rc " _rc
            local ++nfail
            continue
        }
        mata: st_local("f", strofreal(tsg_map("`dn'", "`lav'", tokens("`fac'"))))
        local nfail = `nfail' + `f'
        tempname bs Vs Cs
        mata: dwlsmv_m11_std_g("`bs'", "`Vs'", "`Cs'")
        mata: st_local("f", strofreal(tsg_e2e("`dn'", "`lav'", tokens("`fac'"), "`bs'", "`Vs'")))
        local nfail = `nfail' + `f'
    }
}

* ---- check 5, display ------------------------------------------------------
import delimited using "fixtures/data/mg_E.csv", clear
foreach lv in configural metric {
    quietly dwlsmv `m1', group(g) invariance(`lv') metric(unitvariance)
    tempname b0 V0 bs Vs Cs
    matrix `b0' = e(b)
    matrix `V0' = e(V)
    mata: dwlsmv_m11_std_g("`bs'", "`Vs'", "`Cs'")
    mata: st_local("f", strofreal(tsg_const("`bs'", "`Vs'", "`Cs'", `nfixs')))
    if `f' {
        display as error "  5 `lv': constant rows wrong"
        local ++nfail
    }
    tsg_log "dwlsmv, standardized"
    mata: st_local("nc", strofreal((strlen(TSGLOG) - strlen(subinstr(TSGLOG, "(constrained)", ""))) / 13))
    mata: st_local("nt", strofreal(strpos(TSGLOG, "Equal raw") + strpos(TSGLOG, "equal raw parameters") > 0))
    mata: st_local("nb", strofreal(strpos(TSGLOG, "var(e.y5)") > 0))
    mata: st_local("no", strofreal((strlen(TSGLOG) - strlen(subinstr(TSGLOG, "(omitted)", ""))) / 9))
    local want = cond("`lv'" == "configural", 0, 1)
    if `tsg_rc' | `nc' + `no' != `nfixs' | `nc' == 0 | `nt' != `want' | !`nb' | "`e(cmd)'" != "dwlsmv" | ///
        mreldif(`b0', e(b)) != 0 | mreldif(`V0', e(V)) != 0 {
        display as error "  5 `lv' FAILED (rc `tsg_rc', constrained `nc' plus omitted `no' of `nfixs', note `nt', Theta `nb')"
        local ++nfail
    }
    else display as text "  5 `lv': ok, `nc' constrained and `no' zero rows, note `nt'"
}
tsg_log "dwlsmv `m1', group(g) invariance(scalar) metric(unitvariance) standardized"
mata: TSGFIT = TSGLOG
tsg_log "dwlsmv, standardized"
* the table printed at estimation, from its title to its note, equals replay's
mata: st_local("same", strofreal(tsg_table(TSGFIT) != "" & tsg_table(TSGFIT) == tsg_table(TSGLOG)))
if `tsg_rc' | !`same' {
    display as error "  5 standardized at estimation FAILED (rc `tsg_rc', same table as replay `same')"
    local ++nfail
}
else display as text "  5 standardized at estimation: ok, same table as replay"

* ---- check 6, group values and G = 3 ---------------------------------------
quietly dwlsmv `m1', group(g) invariance(metric) metric(unitvariance)
tempname b12 V12 Cs
mata: dwlsmv_m11_std_g("`b12'", "`V12'", "`Cs'")
preserve
replace g = cond(g == 1, 0, 7)
set seed 20261101
generate double u = runiform()
sort u
drop u
quietly dwlsmv `m1', group(g) invariance(metric) metric(unitvariance)
tempname b07 V07
mata: dwlsmv_m11_std_g("`b07'", "`V07'", "`Cs'")
local cn : colnames `b07'
local hit0 : list posof "0bn.g#c.F" in cn
local hit7 : list posof "7.g#c.F" in cn
if abs(`b07'[1, 1] - `b12'[1, 1]) > 1e-7 | mreldif(`b07', `b12') > 1e-6 | ///
    mreldif(`V07', `V12') > 1e-4 | !`hit0' | !`hit7' {
    display as error "  6a group values 0, 7 FAILED (" mreldif(`b07', `b12') ", " mreldif(`V07', `V12') ")"
    local ++nfail
}
else display as text "  6a group values 0, 7: ok, reldif " mreldif(`b07', `b12')
restore
preserve
replace g = 3 if g == 2 & mod(_n, 2)
quietly dwlsmv `m1', group(g) invariance(metric) metric(unitvariance)
tempname b3 V3
mata: dwlsmv_m11_std_g("`b3'", "`V3'", "`Cs'")
local nI 5
local rows3 = colsof(e(b)) + 3 * `nI'
mata: st_local("f", strofreal(any(diagonal(st_matrix("`V3'")) :< 0) | missing(st_matrix("`b3'"))))
if colsof(`b3') != `rows3' | `f' {
    display as error "  6b G = 3 FAILED"
    local ++nfail
}
else display as text "  6b G = 3: ok, " colsof(`b3') " rows"
restore

* ---- check 7, fail closed (design E, metric: has ties) ---------------------
quietly dwlsmv `m1', group(g) invariance(metric) metric(unitvariance)
local bad1  "a structural path reached"
local bad2  "e(b) row has no group"
local bad3  "e(b) row has no group"
local bad4  "e(b) has a row twice"
local bad5  "e(V) does not fit e(b)"
local bad6  "e(V) stripe differs"
local bad7  "tied rows differ"
local bad8  "tied rows differ"
local bad9  "e(rank) does not fit"
local bad10 "missing value in e(b)"
local bad11 "group metadata missing"
local bad12 "e(K) does not fit"
local bad13 "e(b) lacks a threshold"
local bad14 "tie names no e(b) row"
local bad15 "indicator named twice"
local bad16 "group values not unique"
local bad17 "tie line not recognized"
local bad18 "at or below 0"
forvalues c = 1/18 {
    tsg_log "mata: tsg_bad(`c')"
    mata: st_local("hit", strofreal(strpos(TSGLOG, st_local("bad`c'")) > 0))
    if `tsg_rc' != 3498 | !`hit' {
        display as error "  7.`c' FAILED (rc `tsg_rc', message `hit': `bad`c'')"
        local ++nfail
    }
}
display as text "  7 done"

adopath - "`here'/src/ado"
adopath - "`here'/build"
if `nfail' {
    display as error "test_standardized_g: `nfail' failure(s)"
    exit 9
}
display as text "test_standardized_g: all checks passed"
