*! dwlsmv -- canonicalize factor signs after convergence (Module 7)
*!
*! Executed under StataMP 18 on 2026-09-23. tests/unit/test_fit.do check 13. The first run aborted at rc 3301: a '|' in the apply loop subscripted flip[] with a threshold row's variable index, so the loop now branches on ptype first.
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
*! Load order: dwlsmv_types, dwlsmv_nfac.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_signcanon(real colvector theta,
                                struct dwlsmv_parmap scalar par)
{
    real colvector th, flip, neg
    real scalar nf, f, i, r, c, t, can, anchor, afid

    if (rows(theta) != par.npar) {
        _error("dwlsmv_signcanon: theta is not npar rows long")
    }

    th   = theta
    nf   = dwlsmv_nfac(par)
    flip = J(nf, 1, 0)

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
                  (t == 7 & r == f))) {
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
        if (can & anchor > 0) {
            flip[f] = (th[anchor] < 0)
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
    }
    for (i = 1; i <= par.npar; i++) {
        if (neg[i]) {
            th[i] = -th[i]
        }
    }

    return(th)
}

end
