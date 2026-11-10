*! dwlsmv -- position of a declared name in the factor or indicator list (Mod. 3)
*!
*! Executed under StataMP 18 on 2026-09-22. tests/unit/test_identify.do, whose undeclared-factor and undeclared-indicator rejections are this routine returning 0.
*!
*! Returns the 1-based position, or 0 when the name was never declared. Zero is
*! the signal dwls_identify() turns into "unknown factor" / "unknown indicator"
*! diagnostics, so this routine deliberately does not raise an error itself:
*! the caller knows which of the two lists it searched and can say so, and a
*! message naming the wrong kind of object is worse than no message.
*!
*! Position IS the index. Spec Section 4.8 makes the `factor` declaration order
*! the row and column order of B and Psi, and the `ind` declaration order the
*! order of the observed variables behind the Section 4.2 index map. So a
*! lookup here is how a user-facing name becomes a matrix subscript, and it is
*! the only place that translation happens.
*!
*! Linear search, deliberately. These lists have a handful of entries and are
*! searched a few dozen times while parsing; an index would be a second copy of
*! the same fact for no measurable gain.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_canon_find(string colvector pool, string scalar want)
{
    real scalar i

    for (i = 1; i <= rows(pool); i++) {
        if (pool[i] == want) {
            return(i)
        }
    }

    return(0)
}

end
