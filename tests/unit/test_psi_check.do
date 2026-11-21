*! Gate B3b, the disturbance-variance check (release B, unit B5): dwlsmv_psi_band / _check / _lines / _note.
*!
*! Drafted by an agent 2026-11-20, revised and run by the main session the same day. The tolerance in
*! check 6 (1e-3 relative) was set before any run (measured: 1.0e-5 on Case E, 4e-7 on Case I).
*!
*! Spec Module 12 Release B, B3.6, B3.7 (settled at gate B3b, 2026-11-20, option (2)): band_k =
*! 2 sqrt([H^-1]_kk g H^-1 g'), |psi| <= band "zero within numerical error", psi < -band "negative".
*! Checks (E-114, the B3b package, item B3b.7):
*!   1 a hand-built improper point (psi_F2 < 0 on Case E) reports negative, with its standardized loading
*!   2 a proper fit reports nothing and returns the note it was given
*!   3 psi_k = 0 exactly reports zero-within-error (class 1; class 4 when no band)
*!   4 values around the band: +-3 band, +-0.5 band, +-band exactly and just outside, with the band
*!     held at Case E's fitted value for psi_F2 and the values taken about 0. At psi_F2 = 0.5 band the
*!     command's own wording is asserted (zero, standardized loading on F1 is 1). Then the documented
*!     edge (owner, 2026-11-20): at a residual that is exactly zero g = 0 and the band is 0, so a
*!     psi-hat of +-1e-17 is classed as nothing to report / negative. Asserted with the band set to 0,
*!     not from a fit. A first version fitted H4's known-truth design (s = sigma(theta*) exactly): at
*!     six processors its band was 1.6e-16 and at four exactly 0, so the test depended on the
*!     processor count.
*!   5 the existing rows unchanged: admissible() is 6 x 4, and its note survives byte for byte around
*!     the inserted lines (Case I, the improper fixture)
*!   6 the band's value against an independent recomputation (the formula of tools/b3b_psi_error.do:
*!     cross() and a fresh dwls_delta(), not quadcross() and fit.Delta)
*!   7 singular H: the band is missing and the line says the sign only
*!   8 a non-converged fit gets no line
*! Mutants this is meant to kill: the check removed (1, 3, 4, 5, 8); the factor 2 dropped (6: the band
*! would be half the reference; check 4 cannot see it, its points are relative to the production band);
*! '<=' written '<' (4: +-band exactly); the sign test reversed (1, 4, 5); the convergence guard removed
*! (8); the missing-band branch dropped (7).
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
    quietly do "`f'"
}

foreach c in E I {
    import delimited using "fixtures/data/case_`c'.csv", clear varnames(1) ///
        delimiters(",") stringcols(_all)
    destring _all, replace
    mata: Y`c' = st_data(., .)
}

mata:

string scalar tp_join(string colvector L)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= rows(L); i++) {
        s = s + L[i] + char(10)
    }
    return(s)
}

// Case E, marker metric: F2 is endogenous (its psi is a DISTURBANCE variance), F1 exogenous.
string scalar tp_canon_e()
{
    return(tp_join(("factor F1 F2" \ "ind y1 y2 y3 y4 y5 y6" \
        "metric marker F1 y1" \ "lambda F1 y2 free" \ "lambda F1 y3 free" \
        "metric marker F2 y4" \ "lambda F2 y5 free" \ "lambda F2 y6 free" \
        "psi F1 F1 free" \ "psi F2 F2 free" \ "beta F2 F1 free")))
}

// Metric V (the tool's bp_v): nk first-order factors of three items under G, 'unitvariance G'. Case I uses
// nk = 3 (the tool's I_V), Case H's first twelve items nk = 4 (H4).
string scalar tp_canon_v(real scalar nk)
{
    string colvector L
    string scalar fl, il
    real scalar k, j, y

    fl = ""
    for (k = 1; k <= nk; k++) {
        fl = fl + " F" + strofreal(k)
    }
    il = ""
    for (y = 1; y <= 3 * nk; y++) {
        il = il + " y" + strofreal(y)
    }
    L = ("factor" + fl + " G") \ ("ind" + il)
    y = 0
    for (k = 1; k <= nk; k++) {
        for (j = 1; j <= 3; j++) {
            y++
            if (j == 1) {
                L = L \ ("metric marker F" + strofreal(k) + " y" + strofreal(y))
            }
            else {
                L = L \ ("lambda F" + strofreal(k) + " y" + strofreal(y) + " free")
            }
        }
        L = L \ ("psi F" + strofreal(k) + " F" + strofreal(k) + " free")
        L = L \ ("beta F" + strofreal(k) + " G free")
    }
    L = L \ "metric unitvariance G"
    return(tp_join(L))
}

// theta entry of the free diagonal psi of factor f (row f); 0 if none
real scalar tp_psik(struct dwlsmv_parmap scalar par, real scalar f)
{
    real scalar i

    for (i = 1; i <= par.nall; i++) {
        if (par.ptype[i] == 2 & par.row[i] == f & par.col[i] == f & par.free[i]) {
            return(par.fidx[i])
        }
    }
    return(0)
}

real scalar tp_cls(real matrix C, real scalar k)
{
    real scalar r

    for (r = 1; r <= rows(C); r++) {
        if (C[r, 1] == k) {
            return(C[r, 4])
        }
    }
    return(.)
}

real scalar tp_chk(string scalar what, real scalar ok)
{
    if (ok) {
        printf("{txt}  ok    %s\n", what)
        return(0)
    }
    printf("{err}  FAIL  %s\n", what)
    return(1)
}

// DWLS fit of Y.
struct dwlsmv_fit scalar tp_fit(real matrix Y, real colvector K, string scalar canon)
{
    struct dwlsmv_stage1 vector sts
    struct dwlsmv_parmap scalar par

    sts     = dwlsmv_stage1(1, 1)
    sts[1]  = dwls_sample_stats(Y, K)
    par     = dwls_identify(canon, sts[1].map)
    return(dwls_fit(sts, par))
}

// The band again, by the formula of tools/b3b_psi_error.do: a fresh dwls_delta(), cross() not quadcross(),
// explicit g and H. Same inputs, independent code path.
real colvector tp_bandref(struct dwlsmv_fit scalar f)
{
    struct dwlsmv_stage1 vector sv
    struct dwlsmv_statmap vector maps
    real colvector w, r, g, th, out
    real matrix D, H, Hi
    real scalar i, k, q

    sv      = f.stage1
    maps    = dwlsmv_statmap(1, 1)
    maps[1] = sv[1].map
    th      = f.theta
    w       = 1 :/ diagonal(sv[1].Gamma)
    r       = sv[1].s - dwls_model(th, f.par, maps)
    D       = dwls_delta(th, f.par, maps)
    g       = -2 * D' * (w :* r)
    H       = 2 * cross(D, w, D)
    Hi      = luinv(H)
    q       = g' * Hi * g
    out     = J(rows(th), 1, .)
    for (i = 1; i <= f.par.nall; i++) {
        if (f.par.ptype[i] == 2 & f.par.row[i] == f.par.col[i] & f.par.free[i]) {
            k = f.par.fidx[i]
            out[k] = 2 * sqrt(abs(Hi[k, k] * q))
        }
    }
    return(out)
}

real scalar tp_main(real matrix YE, real matrix YI)
{
    struct dwlsmv_fit scalar fE, fI, fn
    struct dwlsmv_parmap scalar parE, parI, parZ
    struct dwlsmv_stage1 vector svE, svI
    struct dwlsmv_statmap scalar mapE, mapI
    real colvector th, bandE, bandI, bandZ, bandS, bandM, bandRef, vals
    real matrix C, A, Dz
    real scalar nf, kF1, kF2, kZ, kI1, hasneg, b, i, pos, ok, rel, v
    string colvector facE, facI, indI
    string scalar lines, adv, note, trailer
    string colvector tags

    nf   = 0
    facE = ("F1" \ "F2")
    facI = ("F1" \ "F2" \ "F3" \ "G")
    indI = ("y1" \ "y2" \ "y3" \ "y4" \ "y5" \ "y6" \ "y7" \ "y8" \ "y9")

    // ---- the proper fit: Case E ----------------------------------------------
    fE    = tp_fit(YE, J(6, 1, 3), tp_canon_e())
    parE  = fE.par
    svE   = fE.stage1
    mapE  = svE[1].map
    kF1   = tp_psik(parE, 1)
    kF2   = tp_psik(parE, 2)
    nf    = nf + tp_chk("E: converged, both psi are free diagonals", fE.converged == 1 & kF1 > 0 & kF2 > 0)
    bandE = dwlsmv_psi_band(fE.theta, parE, fE.stage1, fE.Delta)
    nf    = nf + tp_chk("E: band is defined at exactly the two free diagonal psi, positive",
                        rows(bandE) == parE.npar & missing(bandE) == parE.npar - 2 &
                        !missing(bandE[kF1]) & !missing(bandE[kF2]) & bandE[kF1] > 0 & bandE[kF2] > 0)
    printf("{txt}  Case E: psi-hat (%g, %g), band (%g, %g)\n", fE.theta[kF1], fE.theta[kF2], bandE[kF1], bandE[kF2])

    // ---- 2: a proper point reports nothing ------------------------------------
    C = dwlsmv_psi_check(fE.theta, parE, bandE)
    nf = nf + tp_chk("2: proper fit, every class 0", rows(C) == 2 & all(C[., 4] :== 0))
    nf = nf + tp_chk("2: note is exactly the note it was given (empty)", dwlsmv_psi_note(fE, facE, "") == "")
    nf = nf + tp_chk("2: note is exactly the note it was given (non-empty)",
                     dwlsmv_psi_note(fE, facE, "sentinel" + char(10)) == "sentinel" + char(10))

    // ---- 1: hand-built improper point -----------------------------------------
    // Row 3 fires at any psi < 0 (SO6.3, Sylvester); the existing rows are not touched by this check.
    th      = fE.theta
    th[kF2] = -0.02
    A       = dwlsmv_admissible(th, parE, mapE)
    nf      = nf + tp_chk("1: the existing row 3 fires at the hand-built point (SO6.3)", rows(A) == 6 & cols(A) == 4 & A[3, 1] == 1)
    C       = dwlsmv_psi_check(th, parE, bandE)
    nf      = nf + tp_chk("1: psi_F2 = -0.02 is class 2 (negative)", tp_cls(C, kF2) == 2)
    nf      = nf + tp_chk("1: psi_F1 unchanged is class 0", tp_cls(C, kF1) == 0)
    lines   = dwlsmv_psi_lines(C, th, parE, mapE.p, facE, hasneg)
    printf("{txt}%s", lines)
    nf      = nf + tp_chk("1: the line names F2's negative disturbance variance and its standardized loading on F1 above 1",
                          hasneg == 1 & strpos(lines, "F2: disturbance variance -0.0200 is negative; its standardized loading on F1 is ") > 0 &
                          strpos(lines, ", above 1.") > 0 & strpos(lines, "zero within") == 0)
    note = dwlsmv_psi_merge("", lines, hasneg)
    nf   = nf + tp_chk("1: merged into an empty note it carries the not-admissible header and closing paragraph",
                       strpos(note, "The fitted solution is not admissible") == 1 &
                       strpos(note, "  Point estimates are reported.") > strpos(note, "F2: disturbance variance"))

    // ---- 3: psi_k = 0 exactly --------------------------------------------------
    th      = fE.theta
    th[kF1] = 0
    C       = dwlsmv_psi_check(th, parE, bandE)
    nf      = nf + tp_chk("3: psi_F1 = 0 exactly is class 1", tp_cls(C, kF1) == 1)
    lines   = dwlsmv_psi_lines(C, th, parE, mapE.p, facE, hasneg)
    printf("{txt}%s", lines)
    nf      = nf + tp_chk("3: wording: an exogenous variance, zero within numerical error, no loading clause",
                          hasneg == 0 & strpos(lines, "F1: variance is zero within numerical error.") > 0 &
                          strpos(lines, "negative") == 0 & strpos(lines, "standardized") == 0)
    bandM   = J(parE.npar, 1, .)
    C       = dwlsmv_psi_check(th, parE, bandM)
    nf      = nf + tp_chk("3: with no band, psi = 0 is class 4", tp_cls(C, kF1) == 4)

    // ---- 4: around the band (Case E's band for psi_F2, values about 0) ------------
    parZ  = parE
    kZ    = kF2
    bandZ = bandE
    b     = bandZ[kZ]
    th      = fE.theta
    th[kZ]  = 0.5 * b
    C       = dwlsmv_psi_check(th, parZ, bandZ)
    lines   = dwlsmv_psi_lines(C, th, parZ, mapE.p, facE, hasneg)
    printf("{txt}%s", lines)
    nf    = nf + tp_chk("4: at 0.5 band the line says zero within numerical error, standardized loading on F1 is 1",
                        hasneg == 0 & strpos(lines, "F2: disturbance variance is zero within numerical error; its standardized loading on F1 is 1.") > 0)
    vals = (3 * b \ -3 * b \ 0.5 * b \ -0.5 * b \ b \ -b \ b * (1 + 1e-9) \ -b * (1 + 1e-9))
    tags = ("+3 band" \ "-3 band" \ "+0.5 band" \ "-0.5 band" \ "+band exactly" \ "-band exactly" \ "just above band" \ "just below -band")
    for (i = 1; i <= rows(vals); i++) {
        th     = fE.theta
        th[kZ] = vals[i]
        C      = dwlsmv_psi_check(th, parZ, bandZ)
        v      = tp_cls(C, kZ)
        if (i == 1 | i == 7) {
            ok = (v == 0)
        }
        else if (i == 2 | i == 8) {
            ok = (v == 2)
        }
        else {
            ok = (v == 1)
        }
        nf = nf + tp_chk("4: psi at " + tags[i] + " gets its class", ok)
    }
    // The documented edge (owner, 2026-11-20; spec Module 12 B3.7): an exactly zero residual gives g = 0
    // and a band of 0, so +-1e-17 is classed by its sign.
    bandM     = bandE
    bandM[kZ] = 0
    th        = fE.theta
    th[kZ]    = 1e-17
    C         = dwlsmv_psi_check(th, parZ, bandM)
    nf        = nf + tp_chk("4: band 0 (exactly zero residual): psi +1e-17 has no line", tp_cls(C, kZ) == 0)
    th[kZ]    = -1e-17
    C         = dwlsmv_psi_check(th, parZ, bandM)
    nf        = nf + tp_chk("4: band 0 (exactly zero residual): psi -1e-17 is negative", tp_cls(C, kZ) == 2)
    th[kZ]    = 0
    C         = dwlsmv_psi_check(th, parZ, bandM)
    nf        = nf + tp_chk("4: band 0 (exactly zero residual): psi 0 exactly is zero within numerical error", tp_cls(C, kZ) == 1)

    // ---- 5: existing rows unchanged, on the improper fixture -------------------
    fI    = tp_fit(YI, J(9, 1, 4), tp_canon_v(3))
    parI  = fI.par
    svI   = fI.stage1
    mapI  = svI[1].map
    kI1   = tp_psik(parI, 1)
    nf    = nf + tp_chk("5: Case I converged", fI.converged == 1 & kI1 > 0)
    bandI = dwlsmv_psi_band(fI.theta, parI, fI.stage1, fI.Delta)
    C     = dwlsmv_psi_check(fI.theta, parI, bandI)
    nf    = nf + tp_chk("5: Case I: three free diagonal psi (G is fixed, not checked), psi_F1 negative (class 2)",
                        rows(C) == 3 & tp_cls(C, kI1) == 2)
    A     = dwlsmv_admissible(fI.theta, parI, mapI)
    nf    = nf + tp_chk("5: dwlsmv_admissible is still 6 x 4", rows(A) == 6 & cols(A) == 4)
    adv   = dwlsmv_admissadvice(A, indI, facI)
    note  = dwlsmv_psi_note(fI, facI, adv)
    printf("{txt}%s", note)
    nf    = nf + tp_chk("5: the note names F1's disturbance variance as negative, standardized loading on G 1.3259 in modulus (G's sign is free until B3.10)",
                        strpos(note, "F1: disturbance variance -0.3742 is negative; its standardized loading on G is ") > 0 &
                        (strpos(note, "is 1.3259, above 1.") > 0 | strpos(note, "is -1.3259, below -1.") > 0))
    trailer = "  Point estimates are reported."
    pos     = strpos(adv, trailer)
    if (adv == "") {
        ok = (strpos(note, "The fitted solution is not admissible") == 1)
    }
    else {
        ok = (pos > 0 & substr(note, 1, pos - 1) == substr(adv, 1, pos - 1) &
              substr(note, strlen(note) - strlen(substr(adv, pos, .)) + 1, .) == substr(adv, pos, .) &
              strlen(note) > strlen(adv))
    }
    nf = nf + tp_chk("5: admissibility note's own lines survive byte for byte around the inserted lines", ok)

    // ---- 6: the band's value against the independent recomputation ------------
    // 1e-3 relative, set BEFORE any run: g is a sum of order-1 terms cancelling to ~1e-11, so cross()
    // against quadcross() can differ in g at about 1e-5 relative. Not to be widened to make a failure pass.
    bandRef = tp_bandref(fE)
    rel = max((abs(bandE[kF1] - bandRef[kF1]) / bandRef[kF1] \ abs(bandE[kF2] - bandRef[kF2]) / bandRef[kF2]))
    printf("{txt}  check 6, Case E: relative difference %g\n", rel)
    nf = nf + tp_chk("6: Case E band (factor 2 applied) agrees with the tool's formula to 1e-3 relative", rel < 1e-3)
    bandRef = tp_bandref(fI)
    rel = abs(bandI[kI1] - bandRef[kI1]) / bandRef[kI1]
    printf("{txt}  check 6, Case I: relative difference %g\n", rel)
    nf = nf + tp_chk("6: Case I band for psi_F1 agrees to 1e-3 relative", rel < 1e-3)

    // ---- 7: singular H -------------------------------------------------------------
    Dz         = fE.Delta
    Dz[., kF2] = J(rows(Dz), 1, 0)
    bandS      = dwlsmv_psi_band(fE.theta, parE, fE.stage1, Dz)
    nf         = nf + tp_chk("7: a zero column in Delta makes every band missing", missing(bandS) == rows(bandS))
    th         = fE.theta
    th[kF2]    = -0.02
    C          = dwlsmv_psi_check(th, parE, bandS)
    lines      = dwlsmv_psi_lines(C, th, parE, mapE.p, facE, hasneg)
    nf         = nf + tp_chk("7: negative psi, no band: class 3, sign stated, band said not computed",
                             tp_cls(C, kF2) == 3 & hasneg == 1 &
                             strpos(lines, "F2: disturbance variance -0.0200 is negative") > 0 &
                             strpos(lines, "could not be computed") > 0)
    nf         = nf + tp_chk("7: positive psi, no band: no line", tp_cls(C, kF1) == 0)

    // ---- 8: a fit that did not converge gets no line --------------------------------
    fn           = fI
    fn.converged = 0
    nf           = nf + tp_chk("8: non-converged improper fit: note is the note it was given",
                               dwlsmv_psi_note(fn, facI, adv) == adv & dwlsmv_psi_note(fn, facI, "") == "")
    nf           = nf + tp_chk("8: the same fit, converged, does produce a line", dwlsmv_psi_note(fI, facI, "") != "")

    // ---- 9: the header with no admissibility note (spec Module 12, B3.8 as built) ---
    // B6 review (2026-11-20; added after its mutant P18 survived): with adv empty, a
    // zero-only line gets "Disturbance variances:", not the inadmissibility header.
    nf = nf + tp_chk("9: zero-only line, no adv: short header, no inadmissibility text",
                     dwlsmv_psi_merge("", "  L" + char(10), 0) == "Disturbance variances:" + char(10) + "  L" + char(10))
    nf = nf + tp_chk("9: negative line, no adv: the inadmissibility header and closing paragraph",
                     strpos(dwlsmv_psi_merge("", "  L" + char(10), 1), "is not admissible") > 0 &
                     strpos(dwlsmv_psi_merge("", "  L" + char(10), 1), "  Point estimates are reported.") > 0)

    return(nf)
}

end

mata: st_numscalar("tp_nfail", tp_main(YE, YI))
if scalar(tp_nfail) > 0 {
    display as error "test_psi_check: " scalar(tp_nfail) " check(s) failed"
    exit 9
}
display as text "test_psi_check: all checks passed"
