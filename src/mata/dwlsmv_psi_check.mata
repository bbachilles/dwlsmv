*! dwlsmv -- the disturbance-variance check (release B, gate B3b; spec Module 12 Release B, B3.6, B3.7, B3.11)
*!
*! Drafted by an agent 2026-11-20, revised by the main session (wording to B3.6), run 2026-11-20.
*! B3.8 change (dwlsmv_psi_so3note(), dwlsmv_psi_hasso(), dwlsmv_psi_unitnote(), the call in
*! dwlsmv_psi_note()) drafted by an agent and executed under StataMP 18 on 2026-11-20:
*! tests/unit/test_notes_so.do (all checks pass; mutants N1-N5 killed).
*!
*! B3.8 (spec Module 12 Release B, "Notes (B3.8)"; SO2.3, SO3.3). Two notes, in plain language.
*!   (a) dwlsmv_psi_so3note(): after the psi lines, once per second-order factor G that has no lambda row,
*!       exactly three free outgoing beta rows (distinct parameters), and no other row but psi_GG (no
*!       covariance, no path into G, no fixed or tied path), when the check classes (class 2) the free
*!       diagonal psi of at least one of its three children as negative. It is added to the psi lines
*!       before they are merged into the admissibility note, so it sits inside that note after the lines.
*!       Class 3 (negative, band not computed) is not counted: the spec says the classes of the check, and
*!       B3.8's sentence rests on a negative that the band has confirmed. Several groups: no note (a
*!       second-order factor is refused there).
*!   (b) dwlsmv_psi_unitnote(fit, metric): for a fit that did not converge, under metric "unitvariance",
*!       when the table has a second-order factor (a factor with no lambda row and a free outgoing beta).
*!       The caller (dwlsmv_m11_fit) appends it after dwlsmv_psi_note().
*!
*! A routine below the locked boundary: dwlsmv_admissible()'s return, every existing row and every e()
*! result are unchanged, and this check does not feed D11's trigger. It reads existing fit fields only
*! (theta, par, stage1, Delta, converged); no struct field is added.
*!
*! THE RULE (owner, gate B3b, 2026-11-20, option (2)). For every FREE DIAGONAL psi (ptype 2, row == col)
*! with theta entry k, at theta-hat, with r = s - sigma(theta), W the stacked weights as the fit uses
*! them (weighting (a) under groups), g = -2 Delta'W r and H = 2 Delta'W Delta (Gauss-Newton):
*!
*!     band_k = 2 * sqrt( [H^-1]_kk * g H^-1 g' )
*!
*! THE FACTOR 2 WAS CHOSEN AFTER MEASURING (gate B3b, 2026-11-20; E-114,
*! the psi measurement, six fits, StataMP 18). It covers the difference-quotient term (<= 0.41 of the
*! first-order bound), Gauss-Newton against the exact H (<= 0.05) and the processor-count spread (<= 1.0).
*! The bound itself is SO6.5(1). The band's evidence is that measurement only; it is independent of the
*! start rule for an item-anchored higher-order factor (dwlsmv_hof_start).
*!
*! CLASSES (dwlsmv_psi_check(), column 4):
*!     0  no line: psi-hat > band (proper), or band missing and psi-hat > 0
*!     1  |psi-hat| <= band: zero within numerical error; NOT called inadmissible
*!     2  psi-hat < -band: negative
*!     3  band missing (H singular or [H^-1]_kk <= 0) and psi-hat < 0: sign only, band not computed
*!     4  band missing and psi-hat == 0: sign only, band not computed
*! A tied psi is ONE theta entry, so it has one band and one line that names every factor/group it covers.
*! A fit that did not converge gets no line at all (dwlsmv_psi_note() returns the note it was given).
*!
*! WORDING (B3.6, adopted 2026-11-19 as shown for veto). Endogenous: "F1: disturbance variance -0.3200 is
*! negative; its standardized loading on G is 1.2600, above 1." Exogenous: "F: variance -0.0500 is
*! negative." In the band: "F1: disturbance variance is zero within numerical error; its standardized
*! loading on G is 1." The standardized-loading clause is printed only where it is defined by one path:
*! one group, exactly one structural parent (a beta row, no gamma), and both factor variances in Phi
*! positive. Elsewhere (several parents, a covariate, several groups) the line stops after the
*! classification; B3.6 shows no wording for those cases (main session, 2026-11-20, mechanical).
*!
*! Load order: after dwlsmv_admissadvice (dwlsmv_admiss_num) and dwls_fit (dwlsmv_fit_stack); Mata
*! resolves both at run time.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// The band, npar x 1: missing at every entry that is not a free diagonal psi, and missing everywhere when
// g and H cannot be formed or H is singular. Delta is the D x npar matrix at theta (a fit passes fit.Delta).
// Nothing passed in is modified.
real colvector dwlsmv_psi_band(real colvector theta,
                               struct dwlsmv_parmap  scalar par,
                               struct dwlsmv_stage1 vector sts,
                               real matrix Delta)
{
    struct dwlsmv_statmap vector maps
    real colvector s, w, r, band
    real rowvector g
    real matrix H, Hi
    real scalar N, i, k, q

    band = J(par.npar, 1, .)
    dwlsmv_fit_stack(sts, s, w, N, maps, 1)
    if (rows(Delta) != rows(s) | cols(Delta) != par.npar | rows(theta) != par.npar) {
        return(band)
    }
    if (missing(Delta) | missing(theta)) {
        return(band)
    }
    r = s - dwls_model(theta, par, maps)
    if (missing(r)) {
        return(band)
    }
    g  = -2 * quadcross(w :* r, Delta)
    H  = 2 * quadcross(Delta, w, Delta)
    Hi = luinv(H)
    if (missing(Hi)) {
        return(band)
    }
    q = g * Hi * g'
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.row[i] == par.col[i] & par.free[i]) {
            k = par.fidx[i]
            // [H^-1]_kk <= 0 means H is not positive definite here: no band, not a guess
            if (Hi[k, k] > 0) {
                band[k] = 2 * sqrt(Hi[k, k] * abs(q))
            }
        }
    }
    return(band)
}

// One value's class (see the header). abs(p) <= b, not <: a value ON the band is inside it.
real scalar dwlsmv_psi_class(real scalar p, real scalar b)
{
    if (missing(b)) {
        if (p < 0) {
            return(3)
        }
        if (p == 0) {
            return(4)
        }
        return(0)
    }
    if (abs(p) <= b) {
        return(1)
    }
    if (p < -b) {
        return(2)
    }
    return(0)
}

// One row per free diagonal-psi theta entry (a tied entry once): (k, psi-hat, band, class).
real matrix dwlsmv_psi_check(real colvector theta,
                             struct dwlsmv_parmap scalar par,
                             real colvector band)
{
    real matrix C
    real colvector seen
    real scalar i, k

    if (rows(band) != par.npar | rows(theta) != par.npar) {
        _error("dwlsmv_psi_check: theta or band is not npar x 1")
    }
    C    = J(0, 4, .)
    seen = J(par.npar, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.row[i] == par.col[i] & par.free[i]) {
            k = par.fidx[i]
            if (seen[k] == 0) {
                seen[k] = 1
                C = C \ (k, theta[k], band[k], dwlsmv_psi_class(theta[k], band[k]))
            }
        }
    }
    return(C)
}

// The standardized loading clause for factor f (B3.6), or "" where it is not defined by one path.
// cl is the class (1 or 2). p is map.p, the item count; theta is the point the classes were read at.
string scalar dwlsmv_psi_stdclause(real colvector theta, struct dwlsmv_parmap scalar par,
                                   real scalar p, real scalar f, string colvector fac,
                                   real scalar cl)
{
    real colvector vals
    real matrix Lam, Psi, B, Ai, Phi
    real scalar i, npar_in, ip, P, nf, bs

    if (max(par.grp) != 1) {
        return("")
    }
    vals    = dwlsmv_parvals(theta, par)
    npar_in = 0
    ip      = 0
    for (i = 1; i <= par.nall; i++) {
        if ((par.ptype[i] == 3 | par.ptype[i] == 8) & par.row[i] == f) {
            if (par.free[i] | vals[i] != 0) {
                npar_in++
                ip = i
            }
        }
    }
    if (npar_in != 1) {
        return("")
    }
    if (par.ptype[ip] != 3) {
        return("")
    }
    P  = par.col[ip]
    nf = dwlsmv_nfac(par)
    dwlsmv_matrices(vals, par, p, nf, Lam, Psi, B)
    Ai = luinv(I(nf) - B)
    if (missing(Ai)) {
        return("")
    }
    Phi = Ai * Psi * Ai'
    if (!(Phi[f, f] > 0) | !(Phi[P, P] > 0)) {
        return("")
    }
    if (cl == 1) {
        return("; its standardized loading on " + fac[P] + " is 1")
    }
    bs = vals[ip] * sqrt(Phi[P, P] / Phi[f, f])
    if (bs > 1) {
        return("; its standardized loading on " + fac[P] + " is " + dwlsmv_admiss_num(bs) + ", above 1")
    }
    if (bs < -1) {
        return("; its standardized loading on " + fac[P] + " is " + dwlsmv_admiss_num(bs) + ", below -1")
    }
    return("")
}

// The text, one line per row of C with class > 0, "" if none. fac: factor names in declaration order;
// p: map.p. hasneg is set to 1 when some line says "negative" (class 2 or 3), else 0.
string scalar dwlsmv_psi_lines(real matrix C, real colvector theta,
                               struct dwlsmv_parmap scalar par,
                               real scalar p, string colvector fac,
                               real scalar hasneg)
{
    string scalar out, nl, who, noun, why, std
    real scalar r, i, k, G, endo, cl, v, first

    nl     = char(10)
    out    = ""
    hasneg = 0
    G      = max(par.grp)
    why    = " (the numerical band could not be computed:" + nl +
             "    the information matrix is singular or not positive definite)"
    for (r = 1; r <= rows(C); r++) {
        k  = C[r, 1]
        v  = C[r, 2]
        cl = C[r, 4]
        if (cl > 0) {
            who   = ""
            first = 0
            for (i = 1; i <= par.nall; i++) {
                if (par.ptype[i] == 2 & par.row[i] == par.col[i] & par.free[i]) {
                    if (par.fidx[i] == k) {
                        if (first == 0) {
                            first = i
                        }
                        if (who != "") {
                            who = who + ", "
                        }
                        who = who + fac[par.row[i]]
                        if (G > 1) {
                            who = who + " (group " + strofreal(par.grp[i]) + ")"
                        }
                    }
                }
            }
            // an endogenous factor (a beta or gamma row into it) has a DISTURBANCE variance
            endo = any(((par.ptype :== 3) :| (par.ptype :== 8)) :& (par.row :== par.row[first]))
            noun = "variance"
            std  = ""
            if (endo) {
                noun = "disturbance variance"
                if (cl == 1 | cl == 2) {
                    std = dwlsmv_psi_stdclause(theta, par, p, par.row[first], fac, cl)
                }
            }
            if (cl == 1) {
                out = out + "    " + who + ": " + noun + " is zero within numerical error" + std + "." + nl
            }
            else if (cl == 2) {
                out = out + "    " + who + ": " + noun + " " + dwlsmv_admiss_num(v) +
                      " is negative" + std + "." + nl
            }
            else if (cl == 3) {
                out = out + "    " + who + ": " + noun + " " + dwlsmv_admiss_num(v) +
                      " is negative" + why + "." + nl
            }
            else {
                out = out + "    " + who + ": " + noun + " is 0" + why + "." + nl
            }
            if (cl == 2 | cl == 3) {
                hasneg = 1
            }
        }
    }
    return(out)
}

// Word-wrap s at 76 columns and indent each line by two spaces (the admissibility note's paragraph style).
string scalar dwlsmv_psi_para(string scalar s)
{
    string scalar w, out
    real scalar i

    w   = dwlsmv_m11_recwrap(s, 76)
    out = ""
    while ((i = strpos(w, char(10))) > 0) {
        out = out + "  " + substr(w, 1, i - 1) + char(10)
        w   = substr(w, i + 1, .)
    }
    return(out)
}

// 1 if the table has a second-order factor: a factor with no lambda row (any row, free or fixed) and at
// least one FREE outgoing beta row (a beta row has the parent in col). One group only is read as one.
real scalar dwlsmv_psi_hasso(struct dwlsmv_parmap scalar par)
{
    real scalar i, j, f, hasl

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 3 & par.free[i]) {
            f    = par.col[i]
            hasl = 0
            for (j = 1; j <= par.nall; j++) {
                if (par.ptype[j] == 1 & par.col[j] == f) {
                    hasl = 1
                }
            }
            if (hasl == 0) {
                return(1)
            }
        }
    }
    return(0)
}

// B3.8 note (a): "" or one paragraph per factor G that has exactly three free, distinct outgoing paths and
// nothing else but its own variance, when one of the three children's free diagonal psi is class 2 in C.
string scalar dwlsmv_psi_so3note(real matrix C, struct dwlsmv_parmap scalar par, string colvector fac)
{
    string scalar out, kids, txt
    real scalar G, nf, i, j, hasl, nb, other, neg, k, r, tied

    out = ""
    if (max(par.grp) != 1) {
        return(out)
    }
    nf = dwlsmv_nfac(par)
    for (G = 1; G <= nf; G++) {
        hasl  = 0
        nb    = 0
        other = 0
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] == 1 & par.col[i] == G) {
                hasl = 1
            }
            else if (par.ptype[i] == 3 & par.col[i] == G) {
                if (par.free[i]) {
                    nb++
                }
                else {
                    other = 1
                }
            }
            else if ((par.ptype[i] == 3 | par.ptype[i] == 8) & par.row[i] == G) {
                other = 1
            }
            else if (par.ptype[i] == 2 & par.row[i] != par.col[i] & (par.row[i] == G | par.col[i] == G)) {
                other = 1
            }
        }
        // free outgoing betas: three, each its own parameter (a tie is a restriction on G)
        tied = 0
        if (hasl == 0 & other == 0 & nb == 3) {
            for (i = 1; i <= par.nall; i++) {
                for (j = i + 1; j <= par.nall; j++) {
                    if (par.ptype[i] == 3 & par.ptype[j] == 3 & par.col[i] == G & par.col[j] == G) {
                        if (par.free[i] & par.free[j] & par.fidx[i] == par.fidx[j]) {
                            tied = 1
                        }
                    }
                }
            }
        }
        if (hasl == 0 & other == 0 & nb == 3 & tied == 0) {
            kids = ""
            neg  = 0
            for (i = 1; i <= par.nall; i++) {
                if (par.ptype[i] == 3 & par.col[i] == G & par.free[i]) {
                    if (kids != "") {
                        kids = kids + ", "
                    }
                    kids = kids + fac[par.row[i]]
                    for (j = 1; j <= par.nall; j++) {
                        if (par.ptype[j] == 2 & par.row[j] == par.row[i] & par.col[j] == par.row[i] & par.free[j]) {
                            k = par.fidx[j]
                            for (r = 1; r <= rows(C); r++) {
                                if (C[r, 1] == k & C[r, 4] == 2) {
                                    neg = 1
                                }
                            }
                        }
                    }
                }
            }
            if (neg) {
                txt = fac[G] + " has exactly three first-order factors (" + kids + "), so this model fits " +
                      "exactly as well as one in which those three factors are simply correlated. The " +
                      "negative disturbance variance means their correlations cannot be represented by a " +
                      "proper second-order factor."
                out = out + dwlsmv_psi_para(txt)
            }
        }
    }
    return(out)
}

// B3.8 note (b): "" unless the fit did not converge, the metric is unitvariance and the table has a
// second-order factor. converged is tested as != 1 (0 and any other value are not a solution).
string scalar dwlsmv_psi_unitnote(struct dwlsmv_fit scalar fit, string scalar metric)
{
    if (fit.converged == 1) {
        return("")
    }
    if (metric != "unitvariance") {
        return("")
    }
    if (dwlsmv_psi_hasso(fit.par) == 0) {
        return("")
    }
    return(dwlsmv_psi_para("This model has a second-order factor. With metric(unitvariance) its " +
                           "first-order factors cannot have a negative disturbance variance; with " +
                           "metric(marker) they can, so a fit that needs one may converge there."))
}

// Put the lines inside the admissibility note adv (dwlsmv_admissadvice()'s or dwlsmv_cx_admissnote()'s).
// adv non-empty: before its closing "  Point estimates are reported." paragraph (appended if that text
// is absent). adv empty: a negative line gets the note's own header and closing paragraph; a zero-only
// line gets the short header "Disturbance variances:" (it is not a statement of inadmissibility).
string scalar dwlsmv_psi_merge(string scalar adv, string scalar lines, real scalar hasneg)
{
    string scalar nl, trailer
    real scalar pos

    if (lines == "") {
        return(adv)
    }
    nl = char(10)
    if (adv != "") {
        pos = strpos(adv, "  Point estimates are reported.")
        if (pos > 0) {
            return(substr(adv, 1, pos - 1) + lines + substr(adv, pos, .))
        }
        return(adv + lines)
    }
    if (hasneg) {
        trailer = "  Point estimates are reported. Standard errors and the WLSMV test assume a" + nl +
                  "  solution inside the allowed range and should not be relied on." + nl
        return("The fitted solution is not admissible (some estimates are outside the range" + nl +
               "a valid model allows):" + nl + lines + trailer)
    }
    return("Disturbance variances:" + nl + lines)
}

// The entry point the Module 11 callers use: the note adv, with this check's lines inside it.
// A fit that did not converge: adv unchanged. fac: factor names in declaration order (a colvector).
string scalar dwlsmv_psi_note(struct dwlsmv_fit scalar fit,
                              string colvector fac,
                              string scalar adv)
{
    struct dwlsmv_stage1 vector sv
    real colvector band
    real matrix C
    string scalar lines
    real scalar hasneg

    if (fit.converged != 1) {
        return(adv)
    }
    sv    = fit.stage1
    band  = dwlsmv_psi_band(fit.theta, fit.par, sv, fit.Delta)
    C     = dwlsmv_psi_check(fit.theta, fit.par, band)
    lines = dwlsmv_psi_lines(C, fit.theta, fit.par, sv[1].map.p, fac, hasneg)
    if (lines != "") {
        lines = lines + dwlsmv_psi_so3note(C, fit.par, fac)
    }
    return(dwlsmv_psi_merge(adv, lines, hasneg))
}

end
