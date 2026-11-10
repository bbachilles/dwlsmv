*! Module 0-1: the sample-statistic index map (Gate 1, spec Section 4.2).
*!
*! Four things are checked, and they are not redundant:
*!
*!   1. ORACLE. The exact label sequence lavInspect(fit, "gamma") produces for
*!      the four-indicator design of Section 4.2, pinned as a literal. This is
*!      the only check here with an external source; the rest is internal
*!      consistency, which cannot detect the whole map being shifted or
*!      ordered wrongly in a self-consistent way.
*!   2. ORDER. Positions from the closed-form arithmetic agree with labels
*!      built by sequential construction. The formula and the counter are
*!      genuinely independent routes to the same answer.
*!   3. BIJECTIVITY. Every position 1..D filled exactly once.
*!   4. INVERSE ROUND-TRIP. pos_inv(s) recovers what s indexes, and feeding
*!      that back through pos_tau/pos_rho returns s.
*!
*! Designs are the derivation file's: the four fixture shapes plus four they do
*! not cover, including p == 1 (no correlation block at all), a 12-variable
*! binary design, and a 2-variable 10-category design.
*!
*! Two Mata mechanics worth knowing, both learned by hitting them here:
*!
*!   The body is a function rather than interactive statements. At Mata's
*!   interactive level every variable is transmorphic, member access on a
*!   transmorphic fails with "struct expected", and the obvious fix does not
*!   work either: a bare declaration at the prompt is parsed as the start of a
*!   function definition and raises "nothing found where ( expected". Inside a
*!   function, declarations behave normally.
*!
*!   This test does NOT load tests/dwlsmv_assert.mata, for one reason only:
*!   positions are exact integers, so a tolerance comparator is the wrong
*!   instrument. It is no longer a loading constraint -- no definition file
*!   issues `mata clear` any more, and tests/unit/test_loading.do holds the
*!   comparator and the structs in one session to keep it that way.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

capture noisily do "src/mata/dwlsmv_types.mata"
if _rc {
    display as error "test_statmap: dwlsmv_types.mata failed to load"
    exit 9
}
foreach f in dwlsmv_statmap_build dwlsmv_pos_tau dwlsmv_pos_rho dwlsmv_pos_inv {
    capture noisily do "src/mata/`f'.mata"
    if _rc {
        display as error "test_statmap: `f'.mata failed to load"
        exit 9
    }
}

mata:

// Reference labels by sequential construction -- a counter, not the closed
// form. This is the independent route the order check needs.
string colvector dwlsmv_t_reflabels(real colvector K)
{
    real scalar p, i, j, c
    real colvector C
    string colvector L

    p = rows(K)
    C = K :- 1
    L = J(0, 1, "")
    for (j = 1; j <= p; j++) {
        for (c = 1; c <= C[j]; c++) {
            L = L \ ("y" + strofreal(j) + "|t" + strofreal(c))
        }
    }
    for (i = 1; i <= p - 1; i++) {
        for (j = i + 1; j <= p; j++) {
            L = L \ ("y" + strofreal(i) + "~~y" + strofreal(j))
        }
    }
    return(L)
}

real scalar dwlsmv_t_statmap()
{
    struct dwlsmv_statmap scalar mD, m
    string colvector FAIL, ORACLE, L
    pointer(real colvector) rowvector DES
    string rowvector NAME
    real colvector K
    real rowvector z
    real scalar d, i, j, c, s, pt, pr, back
    string scalar nm, want

    FAIL = J(0, 1, "")

    // ---- 1. ORACLE --------------------------------------------------------
    // Section 4.2: two 4-category and two 3-category indicators, T = 10, P = 6.
    ORACLE = ("y1|t1" \ "y1|t2" \ "y1|t3" \ "y2|t1" \ "y2|t2" \ "y2|t3" \
              "y3|t1" \ "y3|t2" \ "y4|t1" \ "y4|t2" \
              "y1~~y2" \ "y1~~y3" \ "y1~~y4" \ "y2~~y3" \ "y2~~y4" \ "y3~~y4")

    mD = dwlsmv_statmap_build((4 \ 4 \ 3 \ 3))

    if (mD.D != 16) FAIL = FAIL \ ("oracle: D is " + strofreal(mD.D) + ", expected 16")
    if (mD.T != 10) FAIL = FAIL \ ("oracle: T is " + strofreal(mD.T) + ", expected 10")
    if (mD.P != 6)  FAIL = FAIL \ ("oracle: P is " + strofreal(mD.P) + ", expected 6")

    if (rows(mD.label) == rows(ORACLE)) {
        for (i = 1; i <= rows(ORACLE); i++) {
            if (mD.label[i] != ORACLE[i]) {
                FAIL = FAIL \ ("oracle: position " + strofreal(i) + " is " +
                               mD.label[i] + ", lavaan says " + ORACLE[i])
            }
        }
    }
    else {
        FAIL = FAIL \ ("oracle: label vector has " + strofreal(rows(mD.label)) +
                       " rows, expected " + strofreal(rows(ORACLE)))
    }

    // ---- 2-4. the design sweep -------------------------------------------
    DES  = (&((2 \ 2)), &((2 \ 2 \ 2)), &((3 \ 3)), &((4 \ 4 \ 3 \ 3)),
            &(J(1, 1, 2)), &((6 \ 2 \ 4 \ 3 \ 5 \ 2)), &(J(12, 1, 2)),
            &((10 \ 10)))
    NAME = ("A", "B", "C", "D", "E", "F", "G", "H")

    for (d = 1; d <= cols(DES); d++) {
        K  = *DES[d]
        nm = NAME[d]
        m  = dwlsmv_statmap_build(K)
        L  = dwlsmv_t_reflabels(K)

        if (m.D != rows(L)) {
            FAIL = FAIL \ (nm + ": D is " + strofreal(m.D) +
                           ", construction gives " + strofreal(rows(L)))
            continue
        }
        if (m.D != m.T + m.P) {
            FAIL = FAIL \ (nm + ": D does not equal T + P")
        }

        // 2. ORDER
        for (i = 1; i <= m.D; i++) {
            if (m.label[i] != L[i]) {
                FAIL = FAIL \ (nm + ": position " + strofreal(i) + " is " +
                               m.label[i] + ", construction gives " + L[i])
            }
        }

        // 3. BIJECTIVITY -- every slot filled, no label used twice.
        for (i = 1; i <= m.D; i++) {
            if (m.label[i] == "") {
                FAIL = FAIL \ (nm + ": position " + strofreal(i) + " unfilled")
            }
        }
        if (rows(uniqrows(m.label)) != m.D) {
            FAIL = FAIL \ (nm + ": labels are not unique")
        }

        // Positions computed directly must agree with where the builder put them.
        for (j = 1; j <= m.p; j++) {
            for (c = 1; c <= m.C[j]; c++) {
                pt = dwlsmv_pos_tau(m, j, c)
                if (m.label[pt] != "y" + strofreal(j) + "|t" + strofreal(c)) {
                    FAIL = FAIL \ (nm + ": pos_tau(" + strofreal(j) + "," +
                                   strofreal(c) + ") points at " + m.label[pt])
                }
            }
        }
        for (i = 1; i <= m.p - 1; i++) {
            for (j = i + 1; j <= m.p; j++) {
                pr = dwlsmv_pos_rho(m, i, j)
                if (m.label[pr] != "y" + strofreal(i) + "~~y" + strofreal(j)) {
                    FAIL = FAIL \ (nm + ": pos_rho(" + strofreal(i) + "," +
                                   strofreal(j) + ") points at " + m.label[pr])
                }
            }
        }

        // 4. INVERSE ROUND-TRIP
        for (s = 1; s <= m.D; s++) {
            z = dwlsmv_pos_inv(m, s)
            if (z[1] == 1) {
                back = dwlsmv_pos_tau(m, z[2], z[4])
                want = "y" + strofreal(z[2]) + "|t" + strofreal(z[4])
            }
            else {
                back = dwlsmv_pos_rho(m, z[2], z[3])
                want = "y" + strofreal(z[2]) + "~~y" + strofreal(z[3])
            }
            if (back != s) {
                FAIL = FAIL \ (nm + ": pos_inv(" + strofreal(s) +
                               ") does not round-trip, got " + strofreal(back))
            }
            if (want != m.label[s]) {
                FAIL = FAIL \ (nm + ": pos_inv(" + strofreal(s) + ") says " +
                               want + ", label says " + m.label[s])
            }
        }
    }

    printf("test_statmap: %f design(s) swept, %f failure(s)\n",
           cols(DES), rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) printf("    %s\n", FAIL[i])

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(dwlsmv_t_statmap()))

* Error paths. Each is a caller bug and must abort rather than return a
* plausible number. The map is rebuilt inside each call because a struct
* variable cannot be declared at Mata's interactive level.
local nerr 0
capture mata: dwlsmv_pos_rho(dwlsmv_statmap_build((4 \ 4 \ 3 \ 3)), 3, 2)
if _rc != 3498 {
    display as error "test_statmap: pos_rho on i > j gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_pos_rho(dwlsmv_statmap_build((4 \ 4 \ 3 \ 3)), 2, 2)
if _rc != 3498 {
    display as error "test_statmap: pos_rho on i == j gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_pos_tau(dwlsmv_statmap_build((4 \ 4 \ 3 \ 3)), 1, 99)
if _rc != 3498 {
    display as error "test_statmap: pos_tau on an out-of-range threshold gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_pos_inv(dwlsmv_statmap_build((4 \ 4 \ 3 \ 3)), 0)
if _rc != 3498 {
    display as error "test_statmap: pos_inv on position 0 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}
capture mata: dwlsmv_statmap_build((2 \ 1))
if _rc != 3498 {
    display as error "test_statmap: statmap_build on K < 2 gave rc " _rc ", not 3498"
    local nerr = `nerr' + 1
}

if `nfail' | `nerr' {
    display as error _n "test_statmap: `nfail' map failure(s), `nerr' error-path failure(s)"
    exit 9
}

display as result "test_statmap: index map verified -- oracle, order, bijectivity, inverse"
exit 0
