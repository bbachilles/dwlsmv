*! dwlsmv -- read 'free [<start>]' / 'fixed <value>' off a directive (Module 3)
*!
*! Executed under StataMP 18 on 2026-09-22. tests/unit/test_identify.do rejects a bare 'fixed', a non-keyword in the fourth slot and a non-numeric value, each at rc 3498.
*!
*! Spec Section 4.8's three element directives share one tail:
*!
*!     lambda  <F> <y>      free [<start>] | fixed <value>
*!     beta    <to> <from>  free [<start>] | fixed <value>
*!     psi     <F> <F>      free [<start>] | fixed <value>
*!
*! so the status keyword is always token 4 and the value, when present, token
*! 5. 'isfree' and 'val' are filled by reference; Mata passes by reference, and
*! this is the same shape the struct-filling routines use.
*!
*! 'defval' is the caller's default start value for a FREE parameter whose
*! start was omitted, and the caller chooses it by parameter kind:
*!
*!     lambda          0.7   a plausible standardized loading
*!     beta            0     no effect
*!     psi, diagonal   1     unit disturbance variance
*!     psi, off-diag   0     uncorrelated disturbances
*!
*! THESE DEFAULTS ARE NOT COSMETIC, because Module 3's rank check runs at
*! exactly these values. A default of zero for every loading would make Delta's
*! correlation block vanish and report every model as rank deficient; that is
*! why lambda's default is not 0. Conversely the overparameterization the check
*! exists to catch is detected AT psi_12 = 0 -- see dwls_identify()'s header,
*! which works through why the null direction is exactly (beta, psi_12) there
*! and rotates into the loadings away from it. Changing these numbers changes
*! what the rank check sees.
*!
*! A 'fixed' with no value is an error rather than a fix-to-default: the two
*! readings of 'psi F1 F1 fixed' -- fix at 1, or fix at the default 0 -- differ
*! by an entire identification convention, and guessing between them is how a
*! silently misidentified model gets fitted.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_canon_status(string rowvector w,
                         real scalar defval,
                         real scalar isfree,
                         real scalar val)
{
    real scalar v

    if (cols(w) > 5) {
        errprintf("dwls_identify: '%s %s %s' carries %f token(s); at most 5 are used.\n",
                  w[1], w[2], w[3], cols(w))
        _error("dwls_identify: a directive has extra tokens")
    }

    if (w[4] == "free") {
        isfree = 1
        if (cols(w) == 5) {
            v = strtoreal(w[5])
            if (missing(v)) {
                errprintf("dwls_identify: '%s' is not a number, as a start value for '%s %s %s'.\n",
                          w[5], w[1], w[2], w[3])
                _error("dwls_identify: a start value is not a number")
            }
            val = v
        }
        else {
            val = defval
        }
        return
    }

    if (w[4] == "fixed") {
        isfree = 0
        if (cols(w) < 5) {
            errprintf("dwls_identify: '%s %s %s fixed' does not say what to fix it to.\n",
                      w[1], w[2], w[3])
            _error("dwls_identify: 'fixed' needs a value")
        }
        v = strtoreal(w[5])
        if (missing(v)) {
            errprintf("dwls_identify: '%s' is not a number, as a fixed value for '%s %s %s'.\n",
                      w[5], w[1], w[2], w[3])
            _error("dwls_identify: a fixed value is not a number")
        }
        val = v
        return
    }

    errprintf("dwls_identify: '%s %s %s %s' -- expected 'free' or 'fixed' as the fourth token.\n",
              w[1], w[2], w[3], w[4])
    _error("dwls_identify: expected 'free' or 'fixed'")
}

end
