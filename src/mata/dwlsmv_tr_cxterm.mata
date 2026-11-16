*! dwlsmv -- Module 11: which of sem's parameters concern the covariate (dwlsmv_tr_cxterm)
*!
*! Executed by the main session (Stata 18 MP) and mutant-tested, 2026-11-14 at the latest. Drafted by a read-only agent 2026-11-14; target path
*! src/mata/dwlsmv_tr_cxterm.mata. Remove this banner only after the file has been run.
*!
*!     real scalar dwlsmv_tr_cxterm(string scalar co, string scalar cvn)
*!
*! co is a name in the second column of sem's e(b) stripe; cvn is the covariate's name (""
*! when the model has none). Returns
*!     0  the name does not concern the covariate
*!     1  var(x)        the covariate's variance
*!     2  mean(x)       the covariate's mean
*!     3  cov(x,.) or cov(.,x)   a covariance of the covariate with another variable
*! The paths x -> F and x -> y are not names of this kind: they are stripe rows whose second
*! column is the covariate, and dwlsmv_tr_core reads them as gamma and kappa.
*!
*! One definition, two callers: dwlsmv_tr_core (translate) and dwlsmv_semnames (D6's stripe),
*! so the rows the translator skips and the rows the D6 check leaves out cannot differ.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_tr_cxterm(string scalar co, string scalar cvn)
{
    string scalar t
    real scalar ln, k

    if (cvn == "") {
        return(0)
    }
    if (co == "var(" + cvn + ")") {
        return(1)
    }
    if (co == "mean(" + cvn + ")") {
        return(2)
    }
    if (substr(co, 1, 4) == "cov(") {
        t  = substr(co, 5, strlen(co) - 5)
        ln = strlen(t)
        k  = strlen(cvn) + 1
        if (substr(t, 1, k) == cvn + ",") {
            return(3)
        }
        // the suffix ",x": its first character is at ln - k + 1
        if (ln >= k) {
            if (substr(t, ln - k + 1, k) == "," + cvn) {
                return(3)
            }
        }
    }
    return(0)
}

end
