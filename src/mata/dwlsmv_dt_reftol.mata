*! dwlsmv -- check 6's refit-reproduction tolerance (Module 12)
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9 as a placeholder (missing);
*! set 2026-10-09 by the main session from the measurement below.
*!
*! Spec Module 12 check 6: the refit's e(b), e(T3) and e(bound) agree with the
*! stored values "within a relative tolerance set from the measured
*! processor-count residue ... labelled as chosen after measuring."
*!
*! 1e-10, CHOSEN AFTER MEASURING, by a rule frozen before the measurement
*! (E-026, 2026-10-09): 100 x the observed maximum, rounded up to a
*! power of ten. Observed maximum of max|x - y|/(1 + |y|) over e(b) and e(T3):
*! 8.05e-13 (the D11-recovered marker-trap fit, e(b) element 35, 1 against 2
*! processors), over case D, case E's four fits, the trap fit and synthetic
*! p = 20 and p = 40 designs at processor counts 1, 2, 4, 6 on StataMP 18; each
*! count repeated bitwise identically. Covers a different processor count on the
*! same Stata version only; not a different Stata version, and not G > 1
*! (re-checked when the group code exists). Rerun on polished code 2026-11-02:
*! maximum 9.93e-13, value unchanged (E-088).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real scalar dwlsmv_dt_reftol()
{
    return(1e-10)
}

end
