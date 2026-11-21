*! dwlsmv -- canonicalize factor signs after convergence (Module 7)
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_fit.do check 13. The first run aborted at rc 3301: a '|' in the apply loop subscripted flip[] with a threshold row's variable index, so the loop now branches on ptype first.
*! Pass 2 (B3.10) executed under StataMP 18 on 2026-11-20: tests/unit/test_signcanon_so.do (checks a-k, eight mutants killed) and test_translate_second_order.do check 8.
*!
*!     real colvector dwlsmv_signcanon(real colvector theta,
*!                                     struct dwlsmv_parmap scalar par)
*!
*! Returns theta with every REFLECTABLE factor turned so that its anchor
*! loading is non-negative. Decided 2026-09-23 (spec Section 5, Module 7):
*! optimize the unrestricted parameterization, with no sign constraint and
*! above all no bound -- that would put an artificial boundary into a smooth
*! problem -- and canonicalize once, after convergence, before anything is
*! derived from theta-hat.
*!
*! THE REFLECTION. Reversing factor f is Lambda -> Lambda D, Psi -> D Psi D,
*! B -> D B D with D = diag(+-1). On the parameter table that negates:
*!
*!     lambda   in column f
*!     beta     with exactly one of (to, from) equal to f
*!     psi      off-diagonal, with exactly one index equal to f
*!
*! A beta or psi joining two reflected factors is negated twice and so left
*! alone, which is D B D and D Psi D exactly. Diagonal psi -- the unit
*! variances -- never move, and Phi -> D Phi D, so Lambda Phi Lambda', R_y*,
*! sigma(theta), F and T are all unchanged. Thresholds are not touched.
*!
*! WHICH FACTORS ARE REFLECTABLE. Only those for which the reflection is a
*! symmetry of THIS table: no FIXED nonzero parameter in any of the positions
*! above. A marker factor fails this by construction -- its fixed lambda = 1
*! already identifies the sign -- and so does a factor with any other fixed
*! nonzero loading, path or covariance, since negating a fixed value would be
*! fitting a different model. Such factors are left exactly as they are, even
*! with a negative first free loading; that sign is data, not an arbitrary
*! choice. The criterion is deliberately this and not "unit-variance": it is
*! the reason the sign is arbitrary, stated directly.
*!
*! THE ANCHOR is the factor's first free loading in indicator order. Exactly
*! zero leaves the factor alone. For v1 CFA the anchor is one of the factor's
*! own indicators; ESEM (Module 12) frees cross-loadings, where the first may
*! sit near zero and a sum-of-loadings anchor would be the sturdier rule. Not
*! built: nothing in v1 reaches it.
*!
*! GROUPS (spec 4.8 O16, decided 2026-10-10, option (a)). A factor is reflected
*! in EVERY group at once: factor indices are shared across groups, so the loops
*! above already cover every group's rows, and a tie, which joins the same
*! element across groups (O3), is never broken. The anchor is the factor's
*! first free loading in theta order (the smallest fidx: group 1's when tied,
*! theta order fixed by O6); at G = 1 that is the first in indicator order, as
*! before. An alpha (ptype 7) of the factor is negated too: the factor mean
*! changes sign with the factor. tau, nu and omega never move. In a configural
*! fit a later group may keep a negative first loading.
*!
*! Negation is exact in floating point, so a reflected theta canonicalizes
*! back to the original bit for bit. test_fit.do check 13 asserts that.
*!
*! CONDITIONAL MODE (2026-11-13). The slopes are b = Lambda (I-B)^-1 Gamma_x + K, so reversing
*! factor f is Lambda -> Lambda D, B -> D B D AND Gamma_x -> D Gamma_x: the gamma rows (ptype 8,
*! row = factor) of a reversed factor change sign, and a FIXED nonzero gamma on the factor makes
*! the reflection not a symmetry of the table, as a fixed nonzero loading does. kappa (ptype 9)
*! never moves: K is not scaled by a factor. A table without gamma rows takes no new branch.
*!
*! SECOND-ORDER FACTORS (B3.10, 2026-11-20). A higher-order factor G with no free loading of
*! its own (no anchor: pass 1 leaves it alone) is reflectable, so its sign is arbitrary too.
*! Spec Module 12, Release B: after the first-order factors are anchored (pass 1, above,
*! unchanged), G is reflected if sum_k b*_k < 0 over its factor children k, the STANDARDIZED
*! paths b*_k = b_k sqrt(Phi[G,G] / Phi[k,k]); a tie at exactly 0 leaves G as fitted. Each G
*! is decided separately. A "G" here is a factor that is reflectable, has no free loading, and
*! has at least one FREE outgoing beta (ptype 3, col = G); its children are the rows of those
*! betas. b_k is the path as fitted and s_k = 1 - 2 flip[k] turns it into the path the child's
*! own reflection leaves: the score uses the first-order factors AS ANCHORED, which is why pass
*! 2 follows pass 1. Phi = (I-B)^-1 Psi (I-B)^-T is built here from the table (value =
*! theta where free, the fixed value otherwise; ptype 2 and 3 only); its diagonal is invariant
*! under any reflection (Phi -> D Phi D), so the theta passed in serves. THE PHI GUARD: a
*! child with Phi[k,k] not > 0 (an improper fit, or a missing value) contributes nothing, and
*! a G with Phi[G,G] not > 0 is left as fitted; tested as !(x > 0), true for missing. Without
*! it sqrt of a negative is missing, the score is missing, and "missing < 0" is false: G
*! would stay unflipped, silently, because a missing value is not an error. The score is
*! accumulated in an explicit loop so no missing value is dropped by sum(). ORDER: a G whose
*! child is itself a G (a third level) waits for that child; G's are decided in sweeps, each
*! deciding every unresolved G whose G-children are all decided, until none is unresolved.
*! B is acyclic (Module 3), so this ends; a sweep that decides none leaves the rest unflipped.
*! ONE GROUP ONLY: pass 2 runs when max(par.grp) == 1. Several groups cannot reach it, since
*! B3.5 refuses group() with a second-order factor (a G with no loading has no later-group
*! metric), and a standardization that mixed groups would be a decision nobody has taken.
*! MD (a G with a free item loading, anchor > 0) is covered by pass 1 and skipped by pass 2.
*!
*! Load order: dwlsmv_types, dwlsmv_nfac.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_signcanon(real colvector theta,
                                struct dwlsmv_parmap scalar par)
{
    real colvector th, flip, neg, canv, anchv, nout, pending
    real matrix B, Psi, Ai, Phi
    real scalar nf, f, i, r, c, t, can, anchor, afid, v
    real scalar progress, ready, score, k

    if (rows(theta) != par.npar) {
        _error("dwlsmv_signcanon: theta is not npar rows long")
    }

    th    = theta
    nf    = dwlsmv_nfac(par)
    flip  = J(nf, 1, 0)
    canv  = J(nf, 1, 0)
    anchv = J(nf, 1, 0)
    nout  = J(nf, 1, 0)

    // Pass 1: every reflectable factor with a free loading anchor.
    for (f = 1; f <= nf; f++) {
        can    = 1
        anchor = 0
        afid   = .
        for (i = 1; i <= par.nall; i++) {
            r = par.row[i]
            c = par.col[i]
            t = par.ptype[i]
            if (!((t == 1 & c == f) |
                  (t == 3 & ((r == f) != (c == f))) |
                  (t == 2 & r != c & ((r == f) != (c == f))) |
                  (t == 7 & r == f) |
                  (t == 8 & r == f))) {
                continue
            }
            if (!par.free[i] & par.value[i] != 0) {
                can = 0
            }
            if (t == 1 & par.free[i]) {
                if (par.fidx[i] < afid) {
                    afid   = par.fidx[i]
                    anchor = par.fidx[i]
                }
            }
        }
        canv[f]  = can
        anchv[f] = anchor
        if (can & anchor > 0) {
            flip[f] = (th[anchor] < 0)
        }
    }

    // Pass 2 (B3.10): a reflectable factor with no free loading and a free
    // outgoing path. One group only; see the header.
    if (max(par.grp) == 1) {
        pending = J(nf, 1, 0)
        for (i = 1; i <= par.nall; i++) {
            if (par.ptype[i] == 3 & par.free[i]) {
                nout[par.col[i]] = nout[par.col[i]] + 1
            }
        }
        for (f = 1; f <= nf; f++) {
            if (canv[f] & anchv[f] == 0 & nout[f] > 0) {
                pending[f] = 1
            }
        }

        if (sum(pending) > 0) {
            B   = J(nf, nf, 0)
            Psi = J(nf, nf, 0)
            for (i = 1; i <= par.nall; i++) {
                t = par.ptype[i]
                if (t != 2 & t != 3) {
                    continue
                }
                if (par.free[i]) {
                    v = th[par.fidx[i]]
                }
                else {
                    v = par.value[i]
                }
                if (t == 3) {
                    B[par.row[i], par.col[i]] = v
                }
                else {
                    Psi[par.row[i], par.col[i]] = v
                    Psi[par.col[i], par.row[i]] = v
                }
            }
            Ai = luinv(I(nf) - B)

            if (!missing(Ai)) {
                Phi      = Ai * Psi * Ai'
                progress = 1
                while (progress & sum(pending) > 0) {
                    progress = 0
                    for (f = 1; f <= nf; f++) {
                        if (!pending[f]) {
                            continue
                        }
                        // Wait while a child is itself a G not yet decided.
                        ready = 1
                        for (i = 1; i <= par.nall; i++) {
                            if (par.ptype[i] == 3 & par.free[i] & par.col[i] == f) {
                                if (pending[par.row[i]]) {
                                    ready = 0
                                }
                            }
                        }
                        if (!ready) {
                            continue
                        }
                        if (Phi[f, f] > 0) {
                            score = 0
                            for (i = 1; i <= par.nall; i++) {
                                if (par.ptype[i] == 3 & par.free[i] & par.col[i] == f) {
                                    k = par.row[i]
                                    if (Phi[k, k] > 0) {
                                        score = score + (1 - 2 * flip[k]) * th[par.fidx[i]] * sqrt(Phi[f, f] / Phi[k, k])
                                    }
                                }
                            }
                            flip[f] = (score < 0)
                        }
                        pending[f] = 0
                        progress   = 1
                    }
                }
            }
        }
    }

    if (sum(flip) == 0) {
        return(th)
    }

    // Branch on the type BEFORE subscripting flip: a threshold row's 'row' is
    // a variable index, which can exceed the factor count. Under groups a tied
    // element is one theta entry shared by several rows (spec 4.8 G2), so the
    // rows MARK the entries and each is negated once below: negating per row
    // would flip a tied entry twice, i.e. not at all.
    neg = J(par.npar, 1, 0)
    for (i = 1; i <= par.nall; i++) {
        if (!par.free[i]) continue
        r = par.row[i]
        c = par.col[i]
        t = par.ptype[i]
        if (t == 1) {
            if (flip[c]) {
                neg[par.fidx[i]] = 1
            }
        }
        else if (t == 3 | (t == 2 & r != c)) {
            if (flip[r] != flip[c]) {
                neg[par.fidx[i]] = 1
            }
        }
        else if (t == 7) {
            if (flip[r]) {
                neg[par.fidx[i]] = 1
            }
        }
        else if (t == 8) {
            if (flip[r]) {
                neg[par.fidx[i]] = 1
            }
        }
    }
    for (i = 1; i <= par.npar; i++) {
        if (neg[i]) {
            th[i] = -th[i]
        }
    }

    return(th)
}

end
