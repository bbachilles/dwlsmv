*! dwlsmv -- Module 3: the conditional-covariate directives (spec, Conditional covariates)
*!
*! Executed 2026-11-13 by the main session (Stata 18 MP, suite green). Drafted by a read-only agent 2026-11-13 (revision 2, to the owner's
*! decisions of 2026-11-13); the main session loads and runs it
*! (tests/unit/test_condx_identify.do).
*!
*! Two functions, called by dwls_identify() only:
*!
*!   dwlsmv_id_cxsplit()  reads the RAW model text, line by line, and separates the 'covariate',
*!                        'gamma' and 'kappa' lines from the rest, so dwlsmv_id_parse() (and
*!                        dwlsmv_id_groups()) never see them and are unchanged in what they
*!                        accept. Every line it returns carries its ORIGINAL line number in the
*!                        text (blank and comment lines count), so a message from here or from
*!                        dwlsmv_id_parse() (which receives the numbers as its optional last
*!                        argument) points at the line the user wrote. Each line is cleaned by
*!                        dwlsmv_canon_lines() itself, one raw line at a time, so there is no
*!                        second copy of the comment/blank rules to drift. Refuses, at rc 3498:
*!                        more than one covariate (one directive line naming several, or a second
*!                        'covariate' line), any covariate under a 'groups' text, a gamma or kappa
*!                        line with no 'covariate' line. One binary covariate, one group, is all
*!                        that is specified (spec, Conditional covariates).
*!   dwlsmv_id_cx()       parses the gamma and kappa lines against the declared factors,
*!                        indicators and covariate; refuses a covariate named like a factor or an
*!                        indicator; and applies the by-name anchor screen below.
*!
*! gmat and kmat are n x 4: row, col, free, value. gamma: row = factor, col = covariate; kappa:
*! row = indicator, col = covariate (ptypes 8 and 9). Starts are 0, as for beta.
*!
*! PRESENT (owner's decisions, 2026-11-13; @sec-condx-anchor, "What present means"). A loading
*! or a beta is present when it is free, or fixed at a value other than 0; fixed at 0, or with no
*! directive, it is absent. A KAPPA is different: only a FREE kappa absorbs the covariate's effect,
*! so an indicator whose kappa is absent or fixed at ANY value (zero or not) is kappa-free.
*!
*! ANCHOR SCREEN N (owner's decision, 2026-11-13; @sec-condx-anchor, "The by-name check"). For
*! every factor f whose gamma is FREE, let R(f) be f together with every factor reachable from f
*! along present betas (`beta F2 F1` means F1 predicts F2, so F2 is reachable from F1). Some
*! kappa-free indicator must have a present loading on a factor of R(f); otherwise the model is
*! refused at rc 3498, naming f. Factors with no free gamma (absent, or fixed at any value) are not
*! checked, even when every one of their indicators carries a kappa.
*!
*! This is a NECESSARY screen, not the full condition: a failure means the factor's column of
*! (Lambda (I-B)^-1) restricted to the kappa-free rows is identically zero at every parameter
*! value. A model that passes can still be unidentified (one anchor shared by two free-gamma
*! factors, a shared descendant, a free beta or loading sitting at 0). The rank check of Delta
*! that follows in dwls_identify() is the exact test and stays authoritative.
*!
*! Requires dwlsmv_types, dwlsmv_canon_lines, dwlsmv_canon_find, dwlsmv_canon_status.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_id_cxsplit(string scalar modelspec,
                       string colvector rest,
                       real colvector lno,
                       string colvector lcx,
                       real colvector lcxno,
                       string colvector cvn,
                       real scalar covno)
{
    string rowvector w
    string colvector c
    string scalar txt, line, d
    real scalar n, k, hasgroups, grpno

    rest      = J(0, 1, "")
    lno       = J(0, 1, .)
    lcx       = J(0, 1, "")
    lcxno     = J(0, 1, .)
    cvn       = J(0, 1, "")
    covno     = .
    hasgroups = 0
    grpno     = .

    // The same carriage-return rule as dwlsmv_canon_lines(); the lines are then split on
    // char(10) exactly as it does, and each is cleaned by it.
    txt = subinstr(modelspec, char(13), "", .)
    n   = 0

    while (1) {
        n = n + 1
        k = strpos(txt, char(10))
        if (k > 0) {
            line = substr(txt, 1, k - 1)
            txt  = substr(txt, k + 1, strlen(txt) - k)
        }
        else {
            line = txt
            txt  = ""
        }

        c = dwlsmv_canon_lines(line)
        if (rows(c) == 1) {
            w = tokens(c[1])
            d = w[1]
            if (d == "groups") {
                if (!hasgroups) {
                    grpno = n
                }
                hasgroups = 1
            }
            if (d == "covariate") {
                if (rows(cvn) > 0) {
                    errprintf("dwls_identify: line %f: a second 'covariate' directive (the first is on line %f).\n",
                              n, covno)
                    errprintf("One binary covariate is specified so far; several covariates are not\n")
                    errprintf("specified yet (spec, Conditional covariates).\n")
                    _error("dwls_identify: more than one covariate is not specified yet")
                }
                if (cols(w) < 2) {
                    errprintf("dwls_identify: line %f: 'covariate' names no covariate.\n", n)
                    _error("dwls_identify: a covariate directive names nothing")
                }
                if (cols(w) > 2) {
                    errprintf("dwls_identify: line %f: 'covariate' names %f covariates. One binary covariate is\n",
                              n, cols(w) - 1)
                    errprintf("specified so far; several covariates are not specified yet (spec,\n")
                    errprintf("Conditional covariates).\n")
                    _error("dwls_identify: more than one covariate is not specified yet")
                }
                cvn   = w[2]
                covno = n
            }
            else if (d == "gamma" | d == "kappa") {
                lcx   = lcx \ c[1]
                lcxno = lcxno \ n
            }
            else {
                rest = rest \ c[1]
                lno  = lno \ n
            }
        }

        if (k == 0) {
            break
        }
    }

    if (rows(lcx) > 0 & rows(cvn) == 0) {
        errprintf("dwls_identify: line %f: a 'gamma' or 'kappa' directive, but no 'covariate' directive\n",
                  lcxno[1])
        errprintf("declares the covariate it names.\n")
        _error("dwls_identify: gamma or kappa without a covariate directive")
    }
    if (rows(cvn) > 0 & hasgroups) {
        errprintf("dwls_identify: line %f: a covariate together with a 'groups' text (line %f). Covariates\n",
                  covno, grpno)
        errprintf("under several groups are not specified yet (spec, Conditional covariates).\n")
        _error("dwls_identify: covariates with several groups are not specified yet")
    }
}

void dwlsmv_id_cx(string colvector lcx,
                  real colvector lcxno,
                  string colvector cvn,
                  real scalar covno,
                  string colvector fac,
                  string colvector ind,
                  real colvector eptype,
                  real colvector erow,
                  real colvector ecol,
                  real colvector efree,
                  real colvector eval,
                  real matrix gmat,
                  real matrix kmat)
{
    string rowvector w
    string scalar d
    real scalar nl, ln, i, k, a, b, h, nfac, isfree, val, hasfreeg, anchored, freekappa
    real matrix reach

    gmat = J(0, 4, .)
    kmat = J(0, 4, .)

    if (rows(cvn) == 0) {
        return
    }
    if (dwlsmv_canon_find(fac, cvn[1]) > 0 | dwlsmv_canon_find(ind, cvn[1]) > 0) {
        errprintf("dwls_identify: line %f: the covariate '%s' has the name of a factor or an indicator.\n",
                  covno, cvn[1])
        errprintf("A covariate needs a name of its own (spec, Conditional covariates).\n")
        _error("dwls_identify: a covariate shares a name with a factor or indicator")
    }

    for (nl = 1; nl <= rows(lcx); nl++) {
        w  = tokens(lcx[nl])
        d  = w[1]
        ln = lcxno[nl]

        // The shape of the line, checked here (not in dwlsmv_canon_status(), whose messages
        // carry no line number) so that every refusal names the line the user wrote. The
        // later call to dwlsmv_canon_status() therefore never fires for these lines.
        if (cols(w) < 4) {
            errprintf("dwls_identify: line %f: '%s' needs two names then 'free' or 'fixed'.\n", ln, d)
            _error("dwls_identify: a parameter directive is incomplete")
        }
        if (cols(w) > 5) {
            errprintf("dwls_identify: line %f: '%s %s %s' carries %f token(s); at most 5 are used.\n",
                      ln, d, w[2], w[3], cols(w))
            _error("dwls_identify: a directive has extra tokens")
        }
        if (w[4] != "free" & w[4] != "fixed") {
            errprintf("dwls_identify: line %f: '%s %s %s' -- expected 'free' or 'fixed' as the fourth token.\n",
                      ln, d, w[2], w[3])
            _error("dwls_identify: expected 'free' or 'fixed'")
        }
        if (w[4] == "fixed" & cols(w) < 5) {
            errprintf("dwls_identify: line %f: '%s %s %s fixed' does not say what to fix it to.\n",
                      ln, d, w[2], w[3])
            _error("dwls_identify: 'fixed' needs a value")
        }
        if (cols(w) == 5) {
            if (missing(strtoreal(w[5]))) {
                errprintf("dwls_identify: line %f: '%s' is not a number, as a start or fixed value for '%s %s %s'.\n",
                          ln, w[5], d, w[2], w[3])
                _error("dwls_identify: a value is not a number")
            }
        }

        b = dwlsmv_canon_find(cvn, w[3])
        if (b == 0) {
            errprintf("dwls_identify: line %f: '%s %s %s' names an undeclared covariate '%s'.\n",
                      ln, d, w[2], w[3], w[3])
            _error("dwls_identify: undeclared covariate in a gamma or kappa directive")
        }
        if (d == "gamma") {
            a = dwlsmv_canon_find(fac, w[2])
            if (a == 0) {
                errprintf("dwls_identify: line %f: 'gamma %s %s' names an undeclared factor '%s'.\n",
                          ln, w[2], w[3], w[2])
                _error("dwls_identify: undeclared factor in a gamma directive")
            }
        }
        else {
            a = dwlsmv_canon_find(ind, w[2])
            if (a == 0) {
                errprintf("dwls_identify: line %f: 'kappa %s %s' names an undeclared indicator '%s'.\n",
                          ln, w[2], w[3], w[2])
                _error("dwls_identify: undeclared indicator in a kappa directive")
            }
        }

        // Duplicate: the same (row, col) in the same block.
        if (d == "gamma") {
            for (i = 1; i <= rows(gmat); i++) {
                if (gmat[i, 1] == a & gmat[i, 2] == b) {
                    errprintf("dwls_identify: line %f: 'gamma %s %s' names an element that is already specified.\n",
                              ln, w[2], w[3])
                    _error("dwls_identify: a parameter is named twice")
                }
            }
        }
        else {
            for (i = 1; i <= rows(kmat); i++) {
                if (kmat[i, 1] == a & kmat[i, 2] == b) {
                    errprintf("dwls_identify: line %f: 'kappa %s %s' names an element that is already specified.\n",
                              ln, w[2], w[3])
                    _error("dwls_identify: a parameter is named twice")
                }
            }
        }

        dwlsmv_canon_status(w, 0, isfree, val)
        if (d == "gamma") {
            gmat = gmat \ (a, b, isfree, val)
        }
        else {
            kmat = kmat \ (a, b, isfree, val)
        }
    }

    // Identification, screen N (@sec-condx-anchor; spec, Conditional covariates). Necessary only:
    // the rank check of Delta that follows in dwls_identify() is the exact test.
    // reach[f, h] = 1 when h is f or is reachable from f along present betas. A beta entry
    // (eptype 3) has erow = to and ecol = from, so it gives reach[ecol, erow] = 1.
    nfac  = rows(fac)
    reach = I(nfac)
    for (i = 1; i <= rows(eptype); i++) {
        if (eptype[i] == 3 & (efree[i] | eval[i] != 0)) {
            reach[ecol[i], erow[i]] = 1
        }
    }
    // Transitive closure (Warshall, k outermost).
    for (k = 1; k <= nfac; k++) {
        for (a = 1; a <= nfac; a++) {
            for (b = 1; b <= nfac; b++) {
                if (reach[a, k] == 1 & reach[k, b] == 1) {
                    reach[a, b] = 1
                }
            }
        }
    }

    for (a = 1; a <= nfac; a++) {
        // Only a factor with a FREE gamma is checked.
        hasfreeg = 0
        for (k = 1; k <= rows(gmat); k++) {
            if (gmat[k, 1] == a & gmat[k, 3] == 1) {
                hasfreeg = 1
            }
        }
        if (hasfreeg == 0) {
            continue
        }

        // An anchor: an indicator with a present loading on some factor of R(a) and no FREE kappa.
        anchored = 0
        for (i = 1; i <= rows(eptype); i++) {
            if (eptype[i] == 1 & (efree[i] | eval[i] != 0)) {
                h = ecol[i]
                if (reach[a, h] == 1) {
                    freekappa = 0
                    for (k = 1; k <= rows(kmat); k++) {
                        if (kmat[k, 1] == erow[i] & kmat[k, 3] == 1) {
                            freekappa = 1
                        }
                    }
                    if (freekappa == 0) {
                        anchored = 1
                    }
                }
            }
        }
        if (anchored == 0) {
            errprintf("dwls_identify: factor %s has a free gamma on %s, but no indicator without a free kappa\n",
                      fac[a], cvn[1])
            errprintf("loads on it or on a factor it predicts (by a free beta or a beta fixed at a nonzero\n")
            errprintf("value). A kappa fixed at any value does not count as free. Free one kappa less on\n")
            errprintf("such an item, or fix the gamma, or add a beta from this factor to one that has\n")
            errprintf("such an item (spec, Conditional covariates, Identification).\n")
            _error("dwls_identify: a free-gamma factor has no indicator without a free kappa")
        }
    }
}

end
