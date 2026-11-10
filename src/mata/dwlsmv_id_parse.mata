*! dwlsmv -- Module 3: one section of a canonical specification (passes 1-3)
*!
*! Executed under StataMP 18 on 2026-10-10. tests/unit/test_identify.do (unedited), tests/unit/test_identify_groups.do.
*!
*! Moved out of dwls_identify() on 2026-10-10 so that a single-group text and
*! each group section of a 'groups' text are parsed by the same code (spec 4.8
*! G1: "each section is a complete single-group text under the existing
*! rules"). With grpmode = 0 this is the old body of dwls_identify() passes
*! 1-3 and its whole-table checks, line for line, so a text without 'groups'
*! parses exactly as before (G0). grpmode = 1 adds only what spec 4.8 G4, G5
*! and G7 require of a section:
*!
*!   - 'nu y', 'omega y', 'alpha F' (ptypes 5, 6, 7; element (row, 1); the
*!     status keyword is the THIRD token). In group 1 every one is refused,
*!     whatever its status (O9). A fixed value or start for omega must be
*!     positive. Without 'groups' the three stay unrecognized directives.
*!   - 'metric' in a group g >= 2 is refused first (G5; precedence decided
*!     2026-10-10: before the missing-psi and collision messages), and
*!     'metric marker' in any group (O2).
*!   - a group g >= 2 needs no metric line; a factor there without a psi
*!     diagonal gets G7's later-group message.
*!
*! Outputs: the section's factor and indicator lists and its element rows
*! (ptype, row, col, free, value), in directive order, metric expanded.
*! Nothing here builds the table, numbers theta or checks rank.
*!
*! Load order: dwlsmv_types, dwlsmv_canon_find, dwlsmv_canon_status.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_id_parse(string colvector lines,
                     struct dwlsmv_statmap scalar map,
                     real scalar g,
                     real scalar grpmode,
                     string colvector fac,
                     string colvector ind,
                     real colvector eptype,
                     real colvector erow,
                     real colvector ecol,
                     real colvector efree,
                     real colvector eval)
{
    string rowvector w
    string scalar d, st
    real colvector mkind
    real matrix Rch
    string scalar cyc
    real scalar nl, k, i, j, r, c, a, b
    real scalar nfac, hit, ptc, defval, isfree, val

    // =====================================================================
    // Pass 1 -- declarations. Position in these lists IS the matrix index.
    // =====================================================================
    fac = J(0, 1, "")
    ind = J(0, 1, "")

    for (nl = 1; nl <= rows(lines); nl++) {
        w = tokens(lines[nl])
        d = w[1]

        if (d != "factor" & d != "ind") {
            continue
        }
        if (cols(w) < 2) {
            errprintf("dwls_identify: line %f: '%s' declares no names.\n", nl, d)
            _error("dwls_identify: a declaration names nothing")
        }

        // A SECOND declaration directive of the same kind is refused, and the
        // reason is Section 4.8's "order within the file is irrelevant".
        // Declaration ORDER is what fixes B's and Psi's row and column order
        // and the observed-variable order behind the index map, so splitting
        // a declaration across two lines would make the order of those two
        // lines matter -- the one thing the section promises does not. One
        // line per kind keeps the promise literally true rather than nearly
        // true, and the grammar is written that way: 'factor <F> [<F> ...]'
        // takes the whole list.
        if (d == "factor" & rows(fac) > 0) {
            errprintf("dwls_identify: a second 'factor' directive on line %f.\n", nl)
            errprintf("Declare every factor on one line: its order fixes B's and Psi's indices,\n")
            errprintf("so two lines would make the order of the lines matter.\n")
            _error("dwls_identify: more than one factor directive")
        }
        if (d == "ind" & rows(ind) > 0) {
            errprintf("dwls_identify: a second 'ind' directive on line %f.\n", nl)
            errprintf("Declare every indicator on one line: its order is the observed-variable\n")
            errprintf("order behind the Section 4.2 index map.\n")
            _error("dwls_identify: more than one ind directive")
        }

        for (k = 2; k <= cols(w); k++) {
            if (d == "factor") {
                if (dwlsmv_canon_find(fac, w[k]) > 0) {
                    errprintf("dwls_identify: factor '%s' is declared twice.\n", w[k])
                    _error("dwls_identify: a factor is declared twice")
                }
                fac = fac \ w[k]
            }
            else {
                if (dwlsmv_canon_find(ind, w[k]) > 0) {
                    errprintf("dwls_identify: indicator '%s' is declared twice.\n", w[k])
                    _error("dwls_identify: an indicator is declared twice")
                }
                ind = ind \ w[k]
            }
        }
    }

    if (rows(fac) == 0) {
        errprintf("dwls_identify: no 'factor' directive; the specification declares no factors.\n")
        _error("dwls_identify: no factor is declared")
    }
    if (rows(ind) == 0) {
        errprintf("dwls_identify: no 'ind' directive; the specification declares no indicators.\n")
        _error("dwls_identify: no indicator is declared")
    }
    if (rows(ind) != map.p) {
        errprintf("dwls_identify: %f indicator(s) declared, but the index map was built for %f.\n",
                  rows(ind), map.p)
        errprintf("The 'ind' order must be the observed-variable order behind the Section 4.2 map.\n")
        _error("dwls_identify: indicator count does not match the index map")
    }

    nfac   = rows(fac)
    eptype = J(0, 1, .)
    erow   = J(0, 1, .)
    ecol   = J(0, 1, .)
    efree  = J(0, 1, .)
    eval   = J(0, 1, .)

    // =====================================================================
    // Pass 2 -- the element directives. One matrix element each.
    // =====================================================================
    for (nl = 1; nl <= rows(lines); nl++) {
        w = tokens(lines[nl])
        d = w[1]

        if (d == "factor" | d == "ind" | d == "metric") {
            continue
        }

        // nu, omega, alpha: groups only (spec 4.8 G4). Status is the THIRD token.
        if (grpmode & (d == "nu" | d == "omega" | d == "alpha")) {
            st = (cols(w) >= 3 ? w[3] : "(none)")
            if (g == 1) {
                errprintf("dwls_identify: group 1, '%s %s %s': %s is not allowed in group 1, whatever its status\n",
                          d, (cols(w) >= 2 ? w[2] : ""), st, d)
                errprintf("(here '%s'). Group 1 is the reference group (nu 0, omega 1, alpha 0); the line is\n", st)
                errprintf("refused because of where it is, not because of its value.\n")
                _error("dwls_identify: nu, omega and alpha are not allowed in group 1")
            }
            if (cols(w) < 3) {
                errprintf("dwls_identify: group %f, line %f: '%s' needs a name then 'free' or 'fixed'.\n",
                          g, nl, d)
                _error("dwls_identify: a parameter directive is incomplete")
            }
            if (cols(w) > 4) {
                errprintf("dwls_identify: group %f: '%s %s' carries %f token(s); at most 4 are used.\n",
                          g, d, w[2], cols(w))
                _error("dwls_identify: a directive has extra tokens")
            }
            if (d == "alpha") {
                a = dwlsmv_canon_find(fac, w[2])
                if (a == 0) {
                    errprintf("dwls_identify: group %f: 'alpha %s' names an undeclared factor.\n", g, w[2])
                    _error("dwls_identify: undeclared factor in an alpha directive")
                }
                ptc    = 7
                defval = 0
            }
            else {
                a = dwlsmv_canon_find(ind, w[2])
                if (a == 0) {
                    errprintf("dwls_identify: group %f: '%s %s' names an undeclared indicator.\n",
                              g, d, w[2])
                    _error("dwls_identify: undeclared indicator in a nu or omega directive")
                }
                ptc    = (d == "nu" ? 5 : 6)
                defval = (d == "nu" ? 0 : 1)
            }
            r = a
            c = 1
            if (w[3] == "free") {
                isfree = 1
                val    = defval
                if (cols(w) == 4) {
                    val = strtoreal(w[4])
                    if (missing(val)) {
                        errprintf("dwls_identify: group %f: '%s' is not a number, as a start for '%s %s'.\n",
                                  g, w[4], d, w[2])
                        _error("dwls_identify: a start value is not a number")
                    }
                }
            }
            else if (w[3] == "fixed") {
                isfree = 0
                if (cols(w) < 4) {
                    errprintf("dwls_identify: group %f: '%s %s fixed' does not say what to fix it to.\n",
                              g, d, w[2])
                    _error("dwls_identify: 'fixed' needs a value")
                }
                val = strtoreal(w[4])
                if (missing(val)) {
                    errprintf("dwls_identify: group %f: '%s' is not a number, as a fixed value for '%s %s'.\n",
                              g, w[4], d, w[2])
                    _error("dwls_identify: a fixed value is not a number")
                }
            }
            else {
                errprintf("dwls_identify: group %f: '%s %s %s' -- expected 'free' or 'fixed' as the third token.\n",
                          g, d, w[2], w[3])
                _error("dwls_identify: expected 'free' or 'fixed'")
            }
            if (ptc == 6 & val <= 0) {
                errprintf("dwls_identify: group %f: 'omega %s' %s at %g; omega is a scale and must be positive.\n",
                          g, w[2], (isfree ? "starts" : "is fixed"), val)
                _error("dwls_identify: omega must be positive")
            }
            for (i = 1; i <= rows(eptype); i++) {
                if (eptype[i] == ptc & erow[i] == r & ecol[i] == c) {
                    errprintf("dwls_identify: group %f: '%s %s' names an element that is already specified.\n",
                              g, d, w[2])
                    _error("dwls_identify: a parameter is named twice")
                }
            }
            eptype = eptype \ ptc
            erow   = erow   \ r
            ecol   = ecol   \ c
            efree  = efree  \ isfree
            eval   = eval   \ val
            continue
        }

        if (d != "lambda" & d != "beta" & d != "psi") {
            if (d == "tau") {
                errprintf("dwls_identify: 'tau' is not a directive. Threshold rows are generated\n")
                errprintf("from the index map, which knows K. Constrained thresholds are measurement\n")
                errprintf("invariance, which is Module 12 and not v1.\n")
            }
            else if (d == "phi") {
                errprintf("dwls_identify: 'phi' is not a directive and cannot become one. Phi_eta is\n")
                errprintf("DERIVED as (I-B)^-1 Psi (I-B)^-T for every factor. Write 'psi', which is the\n")
                errprintf("disturbance (co)variance and the thing that is actually free.\n")
            }
            else if (d == "theta") {
                errprintf("dwls_identify: 'theta' is not a directive. Under the delta parameterization\n")
                errprintf("theta_jj = 1 - (Lambda Phi Lambda')_jj is derived, and a free residual\n")
                errprintf("variance is unidentified from correlations alone.\n")
            }
            else {
                errprintf("dwls_identify: line %f: '%s' is not a directive of the canonical form.\n",
                          nl, d)
                errprintf("Accepted: factor, ind, lambda, beta, psi, metric.\n")
            }
            _error("dwls_identify: unrecognized directive")
        }

        if (cols(w) < 4) {
            errprintf("dwls_identify: line %f: '%s' needs two names then 'free' or 'fixed'.\n",
                      nl, d)
            _error("dwls_identify: a parameter directive is incomplete")
        }

        if (d == "lambda") {
            a = dwlsmv_canon_find(fac, w[2])
            b = dwlsmv_canon_find(ind, w[3])
            if (a == 0) {
                errprintf("dwls_identify: 'lambda %s %s' names an undeclared factor '%s'.\n",
                          w[2], w[3], w[2])
                _error("dwls_identify: undeclared factor in a lambda directive")
            }
            if (b == 0) {
                errprintf("dwls_identify: 'lambda %s %s' names an undeclared indicator '%s'.\n",
                          w[2], w[3], w[3])
                _error("dwls_identify: undeclared indicator in a lambda directive")
            }
            ptc    = 1
            r      = b                      // Lambda is p x nfac
            c      = a
            defval = 0.7
        }
        else if (d == "beta") {
            a = dwlsmv_canon_find(fac, w[2])
            b = dwlsmv_canon_find(fac, w[3])
            if (a == 0 | b == 0) {
                errprintf("dwls_identify: 'beta %s %s' names an undeclared factor.\n", w[2], w[3])
                _error("dwls_identify: undeclared factor in a beta directive")
            }
            if (a == b) {
                errprintf("dwls_identify: 'beta %s %s' regresses a factor on itself.\n", w[2], w[3])
                _error("dwls_identify: a beta directive names one factor twice")
            }
            ptc    = 3
            r      = a                      // row = TO, col = FROM
            c      = b
            defval = 0
        }
        else {
            a = dwlsmv_canon_find(fac, w[2])
            b = dwlsmv_canon_find(fac, w[3])
            if (a == 0 | b == 0) {
                errprintf("dwls_identify: 'psi %s %s' names an undeclared factor.\n", w[2], w[3])
                _error("dwls_identify: undeclared factor in a psi directive")
            }
            ptc    = 2
            r      = max((a, b))            // Psi is symmetric; keep the
            c      = min((a, b))            // lower triangle
            defval = (r == c ? 1 : 0)
        }

        for (i = 1; i <= rows(eptype); i++) {
            if (eptype[i] == ptc & erow[i] == r & ecol[i] == c) {
                errprintf("dwls_identify: '%s %s %s' names an element that is already specified.\n",
                          d, w[2], w[3])
                _error("dwls_identify: a parameter is named twice")
            }
        }

        dwlsmv_canon_status(w, defval, isfree, val)

        eptype = eptype \ ptc
        erow   = erow   \ r
        ecol   = ecol   \ c
        efree  = efree  \ isfree
        eval   = eval   \ val
    }

    // =====================================================================
    // Pass 3 -- metric. Expanded here, and only here.
    // =====================================================================
    mkind = J(nfac, 1, 0)

    for (nl = 1; nl <= rows(lines); nl++) {
        w = tokens(lines[nl])
        if (w[1] != "metric") {
            continue
        }
        // Spec 4.8 G5 / O2. Group g >= 2 first: the precedence decided 2026-10-10.
        if (grpmode & g > 1) {
            errprintf("dwls_identify: group %f carries a 'metric' line. Under groups the metric is group 1's;\n", g)
            errprintf("a later group's factor variance is its own 'psi F F' row, fixed at 1 or free.\n")
            _error("dwls_identify: metric appears in group 1 only")
        }
        if (grpmode & cols(w) >= 2) {
            if (w[2] == "marker") {
                errprintf("dwls_identify: 'metric marker' with groups. The 0.3.0 invariance levels use the\n")
                errprintf("unit-variance metric only (spec 4.8 G5, O2).\n")
                _error("dwls_identify: marker metric is not supported with groups")
            }
        }
        if (cols(w) < 3) {
            errprintf("dwls_identify: line %f: 'metric' needs 'marker <F> <y>' or 'unitvariance <F>'.\n",
                      nl)
            _error("dwls_identify: a metric directive is incomplete")
        }

        if (w[2] == "marker") {
            if (cols(w) != 4) {
                errprintf("dwls_identify: 'metric marker' takes a factor and its marker indicator.\n")
                _error("dwls_identify: 'metric marker' needs a factor and an indicator")
            }
            a = dwlsmv_canon_find(fac, w[3])
            b = dwlsmv_canon_find(ind, w[4])
            if (a == 0) {
                errprintf("dwls_identify: 'metric marker %s %s' names an undeclared factor.\n",
                          w[3], w[4])
                _error("dwls_identify: undeclared factor in a metric directive")
            }
            if (b == 0) {
                errprintf("dwls_identify: 'metric marker %s %s' names an undeclared indicator.\n",
                          w[3], w[4])
                _error("dwls_identify: undeclared indicator in a metric directive")
            }
            ptc = 1
            r   = b
            c   = a
        }
        else if (w[2] == "unitvariance") {
            if (cols(w) != 3) {
                errprintf("dwls_identify: 'metric unitvariance' takes one factor.\n")
                _error("dwls_identify: 'metric unitvariance' takes one factor")
            }
            a = dwlsmv_canon_find(fac, w[3])
            if (a == 0) {
                errprintf("dwls_identify: 'metric unitvariance %s' names an undeclared factor.\n",
                          w[3])
                _error("dwls_identify: undeclared factor in a metric directive")
            }
            ptc = 2
            r   = a
            c   = a
        }
        else {
            errprintf("dwls_identify: 'metric %s' -- expected 'marker' or 'unitvariance'.\n", w[2])
            _error("dwls_identify: unrecognized metric method")
        }

        if (mkind[a] != 0) {
            errprintf("dwls_identify: factor '%s' carries two metric directives.\n", fac[a])
            _error("dwls_identify: a factor carries two metric constraints")
        }
        mkind[a] = (w[2] == "marker" ? 1 : 2)

        for (i = 1; i <= rows(eptype); i++) {
            if (eptype[i] == ptc & erow[i] == r & ecol[i] == c) {
                errprintf("dwls_identify: the metric on '%s' and an explicit directive both constrain\n",
                          fac[a])
                errprintf("the same element. Two constraints on one parameter is the overparameterization\n")
                errprintf("this check exists to catch; drop one of them.\n")
                _error("dwls_identify: a metric directive collides with an explicit one")
            }
        }

        eptype = eptype \ ptc
        erow   = erow   \ r
        ecol   = ecol   \ c
        efree  = efree  \ 0
        eval   = eval   \ 1
    }

    // =====================================================================
    // Identification checks that need the whole section.
    // =====================================================================
    for (a = 1; a <= nfac; a++) {
        // A later group's factor needs no metric line (spec 4.8 G1, G5).
        if (mkind[a] == 0 & !(grpmode & g > 1)) {
            errprintf("dwls_identify: factor '%s' has no metric directive, so its scale is not fixed.\n",
                      fac[a])
            errprintf("Write 'metric marker %s <y>' or 'metric unitvariance %s'.\n", fac[a], fac[a])
            _error("dwls_identify: a factor has no metric constraint")
        }

        hit = 0
        for (i = 1; i <= rows(eptype); i++) {
            if (eptype[i] == 2 & erow[i] == a & ecol[i] == a) {
                hit = 1
            }
        }
        if (hit == 0) {
            if (grpmode & g > 1) {
                errprintf("dwls_identify: group %f: factor '%s' has no psi diagonal. A later group's factor\n",
                          g, fac[a])
                errprintf("variance is its own 'psi %s %s' row, fixed at 1 or free (spec 4.8 G1).\n",
                          fac[a], fac[a])
                _error("dwls_identify: a group other than 1 has a factor with no psi diagonal")
            }
            errprintf("dwls_identify: factor '%s' has no psi diagonal, so its disturbance variance\n",
                      fac[a])
            errprintf("would be zero. Under a marker metric write 'psi %s %s free'.\n",
                      fac[a], fac[a])
            _error("dwls_identify: a factor has no disturbance variance")
        }
        if (grpmode & g > 1) {
            for (i = 1; i <= rows(eptype); i++) {
                if (eptype[i] == 2 & erow[i] == a & ecol[i] == a & !efree[i]) {
                    if (eval[i] != 1) {
                        errprintf("dwls_identify: group %f: psi %s %s is fixed at %g, but a later group's factor\n",
                                  g, fac[a], fac[a], eval[i])
                        errprintf("variance is fixed at 1 or free (spec 4.8 G5).\n")
                        _error("dwls_identify: a later group's factor variance is fixed at a value other than 1")
                    }
                }
            }
        }

        hit = 0
        for (i = 1; i <= rows(eptype); i++) {
            if (eptype[i] == 1 & ecol[i] == a) {
                hit = 1
            }
        }
        if (hit == 0) {
            errprintf("dwls_identify: factor '%s' has no loadings, so nothing observed measures it.\n",
                      fac[a])
            _error("dwls_identify: a factor has no indicators")
        }
    }

    // =====================================================================
    // v1 is RECURSIVE: no beta path may lead from a factor back to itself.
    // Decided 2026-09-23 (spec Section 4.8). B's directed graph has an edge
    // FROM col TO row for every beta that is free or fixed at a nonzero value
    // -- a fixed zero is no path. Its transitive closure (Warshall) has a 1 on
    // the diagonal exactly for the factors that lie on a cycle. Self-loops
    // were refused above, so any cycle here involves two or more factors.
    // =====================================================================
    Rch = J(rows(fac), rows(fac), 0)
    for (i = 1; i <= rows(eptype); i++) {
        if (eptype[i] == 3 & (efree[i] | eval[i] != 0)) {
            Rch[erow[i], ecol[i]] = 1
        }
    }
    for (k = 1; k <= rows(fac); k++) {
        for (i = 1; i <= rows(fac); i++) {
            for (j = 1; j <= rows(fac); j++) {
                if (Rch[i, k] & Rch[k, j]) {
                    Rch[i, j] = 1
                }
            }
        }
    }
    cyc = ""
    for (a = 1; a <= rows(fac); a++) {
        if (Rch[a, a]) {
            cyc = cyc + " " + fac[a]
        }
    }
    if (cyc != "") {
        errprintf("dwls_identify: the beta directives form a cycle through these factors:%s.\n", cyc)
        errprintf("Reciprocal and feedback paths (non-recursive models) are not supported in v1.\n")
        errprintf("At the default starts every reciprocal pair has identical Delta columns, so the\n")
        errprintf("rank check would misreport it as overparameterized; and a feedback loop also\n")
        errprintf("needs a stability condition v1 does not check. See spec Section 4.8.\n")
        _error("dwls_identify: cyclic beta paths are not supported in v1")
    }
}

end
