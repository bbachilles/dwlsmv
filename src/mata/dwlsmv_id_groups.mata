*! dwlsmv -- Module 3: a 'groups' canonical specification (spec 4.8 G1-G7)
*!
*! Executed under StataMP 18 on 2026-10-10 (O2 lift for beta free in every group:
*! 2026-10-10, mutants M48, M64-M66 each turned the test red; widened to a beta
*! fixed at the same nonzero value in every group 2026-10-11, mutants M48, M67,
*! M68 each red; one-group beta at G = 2 admitted 2026-10-11, mutants M69, M70
*! each red). tests/unit/test_identify_groups.do.
*!
*! Called by dwls_identify() when the text has a 'groups' line. Builds ONE
*! stacked parameter table and runs ONE rank check on it (contract C3, spec
*! 4.8 O14). In order:
*!
*!   1. 'groups G' (once, an integer >= 2), then the sections: 'group g' opens
*!      section g; 'groups' and 'tie' lines belong to none (G1).
*!   2. one map per group, and the Stage 1 results: required at G > 1, one per
*!      map, K matching (user, 2026-10-10).
*!   3. each section through dwlsmv_id_parse(), the single-group parser, in
*!      group mode; every group must declare group 1's factors and
*!      indicators in group 1's order.
*!   4. O2 on the assembled rows, whatever group or line order put them there:
*!      a beta free in every group is admitted (user, 2026-10-10, option (a)),
*!      and a beta fixed at the same nonzero value in every group (2026-10-11,
*!      after mgb-fixed-prereg passed); with exactly two groups, a beta free
*!      in one group and absent (or fixed at 0) in the other, either
*!      orientation, any path (2026-10-11, after mgb-absent-prereg passed).
*!      Refused: fixed values differing, fixed in some groups only, free beside
*!      fixed nonzero, and at G > 2 a beta free in some groups only.
*!      (The marker refusal is the metric line's, in the parser.)
*!   5. 'tie' lines, judged in the precedence decided 2026-10-10: type
*!      whitelist (O4), different types, group in range, two operands in one
*!      group, different elements (O3), thresholds with different K (G3),
*!      absent, fixed, an element in two ties, differing lambda starts (O5b).
*!   6. the table, group-major within each block (user, 2026-10-10, option
*!      (a)): blocks tau, lambda, beta, psi, nu, omega, alpha; inside a block
*!      (group, row, col), fixed rows in place, a tied member in its own
*!      group's position. fidx numbers the free rows in table order (= theta
*!      order, O6), a tie taking the number of its first, lowest-group member.
*!   7. starts: an untied tau at its own group's sample threshold, a tied tau
*!      at its lowest-numbered member group's (O5a); every other row as written.
*!   8. dwlsmv_id_rank() on the stacked table.
*!
*! Row names carry the group, "g2:lambda[y1,F]"; the e() naming (D5, D6,
*! D8d) is Module 11's and lands with it.
*!
*! Load order: dwlsmv_types, dwlsmv_pos_tau, dwlsmv_canon_find,
*! dwlsmv_canon_status, dwlsmv_id_parse, dwlsmv_id_tie, dwlsmv_id_rank.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

struct dwlsmv_parmap scalar dwlsmv_id_groups(string colvector lines,
                                             struct dwlsmv_statmap vector maps,
                                             struct dwlsmv_stage1 vector st,
                                             real scalar hasst)
{
    struct dwlsmv_parmap scalar par
    string colvector fac1, ind1, fac, ind, sl
    string rowvector w
    string scalar d, tn
    real colvector secof, seen, ept, erw, ecl, efr, evl
    real colvector Eg, Ept, Er, Ec, Ef, Ev, blk, tf
    real matrix ops, Tl, S
    real scalar G, nl, g, i, j, k, t, q, nt, ngl, cur, pt, b, r, c
    real scalar nall, pos, npar, hit, tid, v1, nfac, p
    string rowvector tyn

    tyn = ("lambda", "psi", "beta", "tau", "nu", "omega", "alpha")

    // ---- 1. groups and sections (G1) ---------------------------------------
    ngl = 0
    for (nl = 1; nl <= rows(lines); nl++) {
        w = tokens(lines[nl])
        if (w[1] == "groups") {
            ngl++
            if (ngl > 1) {
                errprintf("dwls_identify: a second 'groups' line (line %f). Write 'groups G' once.\n", nl)
                _error("dwls_identify: more than one groups directive")
            }
            G = (cols(w) == 2 ? strtoreal(w[2]) : .)
            if (missing(G) | G != floor(G) | G < 2) {
                errprintf("dwls_identify: 'groups' takes one integer, 2 or more; a single group is written\n")
                errprintf("without a 'groups' line (spec 4.8 G1).\n")
                _error("dwls_identify: 'groups' needs one integer, 2 or more")
            }
        }
    }

    secof = J(rows(lines), 1, 0)
    seen  = J(G, 1, 0)
    cur   = 0
    for (nl = 1; nl <= rows(lines); nl++) {
        w = tokens(lines[nl])
        d = w[1]
        if (d == "groups" | d == "tie") {
            cur = 0
            continue
        }
        if (d == "group") {
            g = (cols(w) == 2 ? strtoreal(w[2]) : .)
            if (missing(g) | g != floor(g) | g < 1 | g > G) {
                errprintf("dwls_identify: line %f: 'group' needs one integer in 1..%f.\n", nl, G)
                _error("dwls_identify: group number outside 1..G")
            }
            if (seen[g]) {
                errprintf("dwls_identify: 'group %f' opens a second section (line %f).\n", g, nl)
                _error("dwls_identify: a group section appears twice")
            }
            seen[g] = 1
            cur     = g
            continue
        }
        if (cur == 0) {
            errprintf("dwls_identify: line %f ('%s') lies outside any group section. Every directive but\n",
                      nl, d)
            errprintf("'groups' and 'tie' belongs to the section a 'group g' line opens.\n")
            _error("dwls_identify: a directive lies outside any group section")
        }
        secof[nl] = cur
    }
    if (sum(seen) != G) {
        errprintf("dwls_identify: 'groups %f' but sections were opened for %f group(s); write\n", G, sum(seen))
        errprintf("one 'group g' section for each g = 1..%f.\n", G)
        _error("dwls_identify: group sections are not numbered 1..G")
    }

    // ---- 2. maps and Stage 1 (O14; user 2026-10-10) ------------------------
    if (length(maps) != G) {
        errprintf("dwls_identify: 'groups %f' but %f index map(s) were passed; one per group.\n",
                  G, length(maps))
        _error("dwls_identify: index maps do not match the number of groups")
    }
    if (!hasst) {
        errprintf("dwls_identify: a 'groups' text needs the Stage 1 results, one per group: threshold\n")
        errprintf("starts come from the data (O5), and a data-free start can make an identified\n")
        errprintf("model look rank-deficient (spec 4.8 O14, C3 note).\n")
        _error("dwls_identify: groups need the Stage 1 results")
    }
    if (length(st) != G) {
        errprintf("dwls_identify: %f Stage 1 result(s) for %f groups.\n", length(st), G)
        _error("dwls_identify: one Stage 1 result is needed per group")
    }
    for (g = 1; g <= G; g++) {
        hit = (rows(st[g].map.K) != rows(maps[g].K) | rows(st[g].s) != maps[g].D)
        if (!hit) {
            hit = (st[g].map.K != maps[g].K)    // nested: Mata's | does not short-circuit
        }
        if (hit) {
            errprintf("dwls_identify: group %f's Stage 1 result was not built on its index map.\n", g)
            _error("dwls_identify: a Stage 1 result does not match its group's index map")
        }
    }

    // ---- 3. each section, by the single-group parser --------------------------
    Eg = J(0, 1, .)
    Ept = Eg
    Er  = Eg
    Ec  = Eg
    Ef  = Eg
    Ev  = Eg
    for (g = 1; g <= G; g++) {
        sl = J(0, 1, "")
        for (nl = 1; nl <= rows(lines); nl++) {
            if (secof[nl] == g) {
                sl = sl \ lines[nl]
            }
        }
        dwlsmv_id_parse(sl, maps[g], g, 1, fac, ind, ept, erw, ecl, efr, evl)
        if (g == 1) {
            fac1 = fac
            ind1 = ind
        }
        else {
            hit = (rows(fac) != rows(fac1))
            if (!hit) {
                hit = (fac != fac1)
            }
            if (hit) {
                errprintf("dwls_identify: group %f declares factors that differ from group 1's (names or\n", g)
                errprintf("order). Every group declares group 1's factors in group 1's order (spec 4.8 G1).\n")
                _error("dwls_identify: a group declares different factors from group 1")
            }
            hit = (rows(ind) != rows(ind1))
            if (!hit) {
                hit = (ind != ind1)
            }
            if (hit) {
                errprintf("dwls_identify: group %f declares indicators that differ from group 1's (names or\n", g)
                errprintf("order). Every group declares group 1's indicators in group 1's order (spec 4.8 G1).\n")
                _error("dwls_identify: a group declares different indicators from group 1")
            }
        }
        Eg  = Eg  \ J(rows(ept), 1, g)
        Ept = Ept \ ept
        Er  = Er  \ erw
        Ec  = Ec  \ ecl
        Ef  = Ef  \ efr
        Ev  = Ev  \ evl
    }
    nfac = rows(fac1)
    p    = rows(ind1)

    // ---- 4. O2: structural paths, decided on the assembled rows ---------------
    // Admitted: a beta free in EVERY group (user, 2026-10-10, option (a)), and a
    // beta fixed at the SAME nonzero value in every group (widened 2026-10-11 after
    // mgb-fixed-prereg passed). Any other mix (a value differing across groups,
    // fixed in some groups and absent or free in others) is refused until validated;
    // fixed at 0 everywhere is an absent path. Widened 2026-10-11 after
    // mgb-absent-prereg passed (S1, S2, SM): with exactly TWO groups a free beta
    // may be absent (or fixed at 0) in the other group, in either orientation and
    // on any path; with more than two groups it must still be free in every group.
    // A free beta beside a fixed nonzero one is refused by the first loop.
    for (i = 1; i <= rows(Ept); i++) {
        if (Ept[i] != 3 | Ef[i] | Ev[i] == 0) {
            continue
        }
        for (g = 1; g <= G; g++) {
            hit = 0
            for (k = 1; k <= rows(Ept); k++) {
                if (Eg[k] == g & Ept[k] == 3 & Er[k] == Er[i] & Ec[k] == Ec[i]) {
                    if (!Ef[k] & Ev[k] == Ev[i]) {
                        hit = 1
                    }
                }
            }
            if (!hit) {
                errprintf("dwls_identify: group %f fixes beta[%s,%s] at %g, but group %f does not fix it at\n",
                          Eg[i], fac1[Er[i]], fac1[Ec[i]], Ev[i], g)
                errprintf("the same value (absent, free, or another value). With groups, a fixed structural\n")
                errprintf("path must be fixed at the same value in every group (spec 4.8 O2).\n")
                _error("dwls_identify: a fixed beta must have the same value in every group")
            }
        }
    }
    for (i = 1; i <= rows(Ept); i++) {
        if (Ept[i] != 3 | !Ef[i] | G == 2) {
            continue
        }
        for (g = 1; g <= G; g++) {
            hit = 0
            for (k = 1; k <= rows(Ept); k++) {
                if (Eg[k] == g & Ept[k] == 3 & Er[k] == Er[i] & Ec[k] == Ec[i]) {
                    hit = Ef[k]
                }
            }
            if (!hit) {
                errprintf("dwls_identify: beta[%s,%s] is free in group %f but not in group %f (absent or fixed\n",
                          fac1[Er[i]], fac1[Ec[i]], Eg[i], g)
                errprintf("at 0). A structural path absent from some groups is validated for two groups only;\n")
                errprintf("with %f groups it must be free in every group (spec 4.8 O2).\n", G)
                _error("dwls_identify: with more than two groups, a beta must be free in every group")
            }
        }
    }

    // ---- 5. ties (G2, O3, O4, O5b, G3, G7; precedence 2026-10-10) -----------------
    Tl = J(0, 5, .)          // tie id, ptype, group, row, col
    nt = 0
    for (nl = 1; nl <= rows(lines); nl++) {
        w = tokens(lines[nl])
        if (w[1] != "tie") {
            continue
        }
        ops = dwlsmv_id_tie(w, fac1, ind1)
        nt++
        tn = "tie line " + strofreal(nl)

        for (q = 1; q <= rows(ops); q++) {
            if (ops[q, 1] != 1 & ops[q, 1] != 4) {
                errprintf("dwls_identify: %s joins a '%s' element. Only tau and lambda may be tied in\n",
                          tn, tyn[ops[q, 1]])
                errprintf("0.3.0, the types the invariance levels use (spec 4.8 O4).\n")
                _error("dwls_identify: a tie may join only tau or lambda elements")
            }
        }
        for (q = 2; q <= rows(ops); q++) {
            if (ops[q, 1] != ops[1, 1]) {
                errprintf("dwls_identify: %s joins a tau and a lambda element.\n", tn)
                _error("dwls_identify: a tie joins elements of different types")
            }
        }
        for (q = 1; q <= rows(ops); q++) {
            if (ops[q, 2] < 1 | ops[q, 2] > G) {
                errprintf("dwls_identify: %s names group %f; the groups are 1..%f.\n", tn, ops[q, 2], G)
                _error("dwls_identify: a tie names an element that is not in the table")
            }
        }
        for (q = 1; q <= rows(ops); q++) {
            for (k = q + 1; k <= rows(ops); k++) {
                if (ops[q, 2] == ops[k, 2]) {
                    errprintf("dwls_identify: %s names two elements of group %f. A tie joins one element\n",
                              tn, ops[q, 2])
                    errprintf("per group (spec 4.8 G2).\n")
                    _error("dwls_identify: a tie names two elements of one group")
                }
            }
        }
        for (q = 2; q <= rows(ops); q++) {
            if (ops[q, 3] != ops[1, 3] | ops[q, 4] != ops[1, 4]) {
                errprintf("dwls_identify: %s joins different elements: '%s' and '%s'.\n", tn,
                          dwlsmv_id_opname(ops[1, .], fac1, ind1), dwlsmv_id_opname(ops[q, .], fac1, ind1))
                errprintf("Operands name the same element in each group, differing only in the group\n")
                errprintf("(spec 4.8 O3).\n")
                _error("dwls_identify: a tie joins different elements across groups")
            }
        }
        if (ops[1, 1] == 4 & ops[1, 3] > 0) {
            r = ops[1, 3]
            for (q = 2; q <= rows(ops); q++) {
                if (maps[ops[q, 2]].K[r] != maps[ops[1, 2]].K[r]) {
                    errprintf("dwls_identify: %s ties thresholds of '%s', which has %f categories in group %f\n",
                              tn, ind1[r], maps[ops[1, 2]].K[r], ops[1, 2])
                    errprintf("and %f in group %f. That is an incompatible restriction, not a partial tie, and\n",
                              maps[ops[q, 2]].K[r], ops[q, 2])
                    errprintf("partial threshold invariance is not in 0.3.0 (spec 4.8 G3).\n")
                    _error("dwls_identify: a tie joins thresholds of indicators with different K")
                }
            }
        }
        for (q = 1; q <= rows(ops); q++) {
            g = ops[q, 2]
            r = ops[q, 3]
            c = ops[q, 4]
            hit = 0
            if (ops[q, 1] == 4) {
                hit = (r >= 1 & r <= p)
                if (hit) {
                    hit = (c >= 1 & c <= maps[g].C[r])
                }
            }
            else {
                for (i = 1; i <= rows(Ept); i++) {
                    if (Eg[i] == g & Ept[i] == 1 & Er[i] == r & Ec[i] == c) {
                        hit = i
                    }
                }
            }
            if (!hit) {
                errprintf("dwls_identify: %s, operand %f, names an element that group %f does not have.\n",
                          tn, q, g)
                _error("dwls_identify: a tie names an element that is not in the table")
            }
            if (ops[q, 1] == 1) {
                if (!Ef[hit]) {
                    errprintf("dwls_identify: %s, operand %f, is fixed in group %f; only free elements are\n",
                              tn, q, g)
                    errprintf("tied (spec 4.8 G2).\n")
                    _error("dwls_identify: a tie names a fixed element")
                }
            }
        }
        for (q = 1; q <= rows(ops); q++) {
            for (i = 1; i <= rows(Tl); i++) {
                if (Tl[i, 2] == ops[q, 1] & Tl[i, 3] == ops[q, 2] & Tl[i, 4] == ops[q, 3] &
                    Tl[i, 5] == ops[q, 4]) {
                    errprintf("dwls_identify: %s names an element (operand %f) already in another tie; an\n",
                              tn, q)
                    errprintf("equality class is one line (spec 4.8 G2).\n")
                    _error("dwls_identify: an element appears in a tie more than once")
                }
            }
        }
        // O5(b): tied lambda operands must carry one start. tau never: its start
        // is its lowest member group's sample threshold.
        if (ops[1, 1] == 1) {
            v1 = .
            for (q = 1; q <= rows(ops); q++) {
                for (i = 1; i <= rows(Ept); i++) {
                    if (Eg[i] == ops[q, 2] & Ept[i] == 1 & Er[i] == ops[q, 3] & Ec[i] == ops[q, 4]) {
                        if (q == 1) {
                            v1 = Ev[i]
                        }
                        else if (Ev[i] != v1) {
                            errprintf("dwls_identify: %s ties lambda[%s,%s] with start %g in group %f to start %g in\n",
                                      tn, ind1[ops[q, 3]], fac1[ops[q, 4]], v1, ops[1, 2], Ev[i])
                            errprintf("group %f. Tied elements share one start; nothing picks a winner (spec 4.8 O5).\n",
                                      ops[q, 2])
                            _error("dwls_identify: tied elements carry different start values")
                        }
                    }
                }
            }
        }
        Tl = Tl \ (J(rows(ops), 1, nt), ops)
    }

    // ---- 6. the table: blocks tau, lambda, beta, psi, nu, omega, alpha ------------
    nall = rows(Ept)
    for (g = 1; g <= G; g++) {
        nall = nall + maps[g].T
    }
    par.nall  = nall
    par.name  = J(nall, 1, "")
    par.ptype = J(nall, 1, 0)
    par.row   = J(nall, 1, 0)
    par.col   = J(nall, 1, 0)
    par.free  = J(nall, 1, 0)
    par.value = J(nall, 1, 0)
    par.fidx  = J(nall, 1, 0)
    par.grp   = J(nall, 1, 0)

    blk = (4 \ 1 \ 3 \ 2 \ 5 \ 6 \ 7)
    pos = 0
    for (b = 1; b <= rows(blk); b++) {
        pt = blk[b]
        for (g = 1; g <= G; g++) {
            if (pt == 4) {
                for (i = 1; i <= maps[g].T; i++) {
                    pos++
                    par.name[pos]  = "g" + strofreal(g) + ":" + maps[g].label[i]
                    par.ptype[pos] = 4
                    par.row[pos]   = maps[g].var_i[i]
                    par.col[pos]   = maps[g].thr_k[i]
                    par.free[pos]  = 1
                    par.grp[pos]   = g
                }
                continue
            }
            S = J(0, 3, .)
            for (i = 1; i <= rows(Ept); i++) {
                if (Eg[i] == g & Ept[i] == pt) {
                    S = S \ (Er[i], Ec[i], i)
                }
            }
            if (rows(S) == 0) {
                continue
            }
            S = sort(S, (1, 2))
            for (k = 1; k <= rows(S); k++) {
                i = S[k, 3]
                r = Er[i]
                c = Ec[i]
                pos++
                if (pt == 1)      par.name[pos] = "lambda[" + ind1[r] + "," + fac1[c] + "]"
                else if (pt == 3) par.name[pos] = "beta[" + fac1[r] + "," + fac1[c] + "]"
                else if (pt == 2) par.name[pos] = "psi[" + fac1[r] + "," + fac1[c] + "]"
                else if (pt == 5) par.name[pos] = "nu[" + ind1[r] + "]"
                else if (pt == 6) par.name[pos] = "omega[" + ind1[r] + "]"
                else              par.name[pos] = "alpha[" + fac1[r] + "]"
                par.name[pos]  = "g" + strofreal(g) + ":" + par.name[pos]
                par.ptype[pos] = pt
                par.row[pos]   = r
                par.col[pos]   = c
                par.free[pos]  = Ef[i]
                par.value[pos] = Ev[i]
                par.grp[pos]   = g
            }
        }
    }
    if (pos != nall) {
        _error("dwls_identify: the parameter table lost a row while being ordered")
    }

    // fidx in table order; a tie takes its first (lowest-group) member's number
    tf   = J(nt, 1, 0)
    npar = 0
    for (i = 1; i <= nall; i++) {
        if (!par.free[i]) {
            continue
        }
        tid = 0
        for (k = 1; k <= rows(Tl); k++) {
            if (Tl[k, 2] == par.ptype[i] & Tl[k, 3] == par.grp[i] & Tl[k, 4] == par.row[i] &
                Tl[k, 5] == par.col[i]) {
                tid = Tl[k, 1]
            }
        }
        if (tid > 0) {
            if (tf[tid] == 0) {
                npar++
                tf[tid] = npar
            }
            par.fidx[i] = tf[tid]
        }
        else {
            npar++
            par.fidx[i] = npar
        }
    }
    par.npar = npar

    // ---- 7. threshold starts from Stage 1 (O5a) -------------------------------------
    for (i = 1; i <= nall; i++) {
        if (par.ptype[i] != 4) {
            continue
        }
        g   = par.grp[i]
        tid = 0
        for (k = 1; k <= rows(Tl); k++) {
            if (Tl[k, 2] == 4 & Tl[k, 3] == g & Tl[k, 4] == par.row[i] & Tl[k, 5] == par.col[i]) {
                tid = Tl[k, 1]
            }
        }
        if (tid) {
            g = min(select(Tl[., 3], Tl[., 1] :== tid))
            par.value[i] = st[g].s[dwlsmv_pos_tau(maps[g], par.row[i], par.col[i])]
        }
        else {
            par.value[i] = st[g].s[dwlsmv_pos_tau(maps[g], par.row[i], par.col[i])]
        }
    }

    // ---- 8. one rank check, on the stacked table (C3) -------------------------------
    dwlsmv_id_rank(par, maps)
    return(par)
}

end
