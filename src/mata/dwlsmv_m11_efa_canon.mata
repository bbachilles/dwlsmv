*! dwlsmv -- Module 11: the EFA block's echelon canonical text (0.3.0, Unit 6)
*!
*! Executed under StataMP 18 on 2026-10-29: tests/unit/test_efa_surface.do.
*!
*!     string scalar dwlsmv_m11_efa_canon(string scalar canon,
*!                                        string scalar efalist)
*!
*! canon is _dwlsmv_translate's text for the user's paths under
*! metric(unitvariance), lines separated by ';'. efalist is efa()'s namelist.
*! Spec section 5 Module 12, Surface (2026-10-29): in 0.3.0 the model must be
*! exactly one complete EFA block. Accepted lines, and nothing else:
*!
*!     factor ... / ind ...                  the block's factors and items
*!     lambda F y free                       one for every factor x item
*!     metric unitvariance F                 one per factor (injected)
*!     psi F G free, F != G                  sem's default factor covariance
*!
*! Returns the echelon text dwls_identify() fits (as test_efa_unrotated.do):
*! factor k (in the factor line's order) fixed at 0 on the first k-1 items (in
*! the ind line's order), every psi covariance fixed at 0, so Phi = I. Rows are
*! fixed, never dropped, so every name in sem's stripe is still built (D6).
*!
*! On a refusal it returns "" and sets the local efamsg; dwlsmv.ado prints it
*! and exits 198 before anything is fitted. Refused: a name in efa() that is
*! not a latent variable; a latent variable not in efa(); a missing loading
*! (a partial block); any other line (paths, constraints, a fixed covariance
*! or loading); more factors than leave df >= 0, df = p(p-1)/2 - (pm -
*! m(m-1)/2) (thresholds are saturated). df = 0 and m = 1 are accepted.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

string scalar dwlsmv_m11_efa_canon(string scalar canon, string scalar efalist)
{
    string rowvector lines, w, fac, ind, efa
    string scalar out, msg
    real matrix seen
    real scalar i, j, k, f, y, p, m, df, nl, npsi

    lines = select(tokens(canon, ";"), tokens(canon, ";") :!= ";")
    nl    = cols(lines)
    fac   = tokens(substr(lines[1], 8, .))
    ind   = tokens(substr(lines[2], 5, .))
    efa   = tokens(efalist)
    p     = cols(ind)
    m     = cols(fac)

    for (k = 1; k <= cols(efa); k++) {
        if (!anyof(fac, efa[k])) {
            st_local("efamsg", "efa(): " + efa[k] + " is not a latent variable of the model.")
            return("")
        }
        if (sum(efa :== efa[k]) > 1) {
            st_local("efamsg", "efa(): " + efa[k] + " is listed twice.")
            return("")
        }
    }
    for (k = 1; k <= m; k++) {
        if (!anyof(efa, fac[k])) {
            st_local("efamsg", "efa(): every latent variable must be in the EFA block in this version; " +
                     fac[k] + " is not.")
            return("")
        }
    }

    seen = J(m, p, 0)
    npsi = 0
    for (i = 3; i <= nl; i++) {
        w = tokens(lines[i])
        if (w[1] == "lambda" & cols(w) == 4) {
            if (w[4] == "free") {
                f = selectindex(fac :== w[2])
                y = selectindex(ind :== w[3])
                seen[f, y] = 1
                continue
            }
        }
        else if (w[1] == "metric" & cols(w) == 3) {
            if (w[2] == "unitvariance") {
                continue
            }
        }
        else if (w[1] == "psi" & cols(w) == 4) {
            if (w[2] != w[3] & w[4] == "free") {
                npsi++
                continue
            }
            if (w[2] != w[3]) {
                st_local("efamsg", "efa(): constraints on the factor covariances are not available with efa();" +
                         " the rotation sets the factor correlations.")
                return("")
            }
        }
        st_local("efamsg", "efa(): only the EFA block itself is available with efa() in this version;" +
                 " the model also has '" + lines[i] + "'.")
        return("")
    }
    for (k = 1; k <= m; k++) {
        for (j = 1; j <= p; j++) {
            if (!seen[k, j]) {
                st_local("efamsg", "efa(): every factor in the block must load on every indicator;" +
                         " " + fac[k] + " -> " + ind[j] + " is missing.")
                return("")
            }
        }
    }
    // sem frees every covariance among exogenous factors by default, so a
    // missing psi row is a user constraint (cov(F1*F2@0) leaves no row at all)
    if (npsi != m * (m - 1) / 2) {
        st_local("efamsg", "efa(): constraints on the factor covariances are not available with efa();" +
                 " the rotation sets the factor correlations.")
        return("")
    }
    df = p * (p - 1) / 2 - (p * m - m * (m - 1) / 2)
    if (df < 0) {
        st_local("efamsg", "efa(): " + strofreal(m) + " factors on " + strofreal(p) +
                 " indicators leave df = " + strofreal(df) + "; the model is not identified.")
        return("")
    }

    out = lines[1] + ";" + lines[2]
    for (i = 3; i <= nl; i++) {
        w = tokens(lines[i])
        if (w[1] == "lambda") {
            k = selectindex(fac :== w[2])
            j = selectindex(ind :== w[3])
            out = out + ";lambda " + w[2] + " " + w[3] + (j < k ? " fixed 0" : " free")
        }
        else if (w[1] == "psi") {
            out = out + ";psi " + w[2] + " " + w[3] + " fixed 0"
        }
        else {
            out = out + ";" + lines[i]
        }
    }
    msg = ""
    st_local("efamsg", msg)
    return(out)
}

end
