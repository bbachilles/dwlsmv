*! Drafted 2026-09-29 by a read-only agent. Executed by the main session on
*! StataMP 18, 2026-09-29. Fixed: the cell thresholds' selection mask used
*! missing(), which returns a count, not an elementwise mask (rc 3203).
*! Added: argument B (default 0.999). delta_Phi's BUDGET is the maximum over
*! |rho| <= B and thresholds in [-2.5, 2.5], chosen after measuring (user,
*! 2026-09-29); the whole-grid ENVELOPE is printed alongside. Thresholds
*! beyond [-2.5, 2.5] are unmeasured; whether to extend is open.
*!
*! Usage: do tools/measure_binormal.do [csvdir [B]]
*!
*! Measurement of delta_Phi for Stata's binormal(). NOT a test: the suite
*! runner globs tests/unit/test_*.do only. It asserts its INPUT (row counts,
*! parse, reference range) and exits 9 if that is broken; the measurement
*! itself has no pass/fail -- it prints DPHI| lines for a human to judge.
*!
*! Why: the draft @sec-rootfind builds tau_S, tau_rho and the usability cut
*!   delta_P = 4 delta_Phi + 3 eps   (eps = 2^-52)
*! from delta_Phi = 2^-52, the per-corner absolute error of Phi2, which was
*! measured only in R. Production uses binormal(), so it is measured here.
*!
*! Run (repo root, Git Bash):
*!     "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/measure_binormal.do
*!     grep "^DPHI|" measure_binormal.log
*! Optional argument: directory holding the two reference CSVs
*! (default tools/data). Diff the DPHI| lines between Stata 18 and the Stata
*! 16 milestone machine; tools/probe_machine.do found binormal() bitwise equal
*! on its 35 values, so any difference here is a new finding.
*!
*! Reference: binormal_ref_corners.csv and binormal_ref_cells.csv, made by
*! make_binormal_ref.R (160-bit MPFR tanh-sinh quadrature, two independent
*! integration ranges; agreement reported by the REF| lines). Each reference is a double-double
*! ref_hi + ref_lo; the error is formed as (computed - ref_hi) - ref_lo so the
*! reference's own rounding does not enter the measurement.
*!
*! CHOSEN AFTER MEASURING / LOOKING: the (h, k, rho) grid and the cell
*! threshold sets were chosen after looking at the lavaan threshold estimates
*! in fixtures/replication/{mc_c_998, mc_c_999, mc_B_995, mc_B_9975}. The
*! regime cut |rho| >= 0.97 was chosen after the rho~-0.98 underflow was seen.
*! delta_Phi = 2^-52 (the value being tested) was chosen after measuring R.
*! Nothing in this file is tuned to Stata's output.

version 16
clear all
set more off
capture set processors 1

args csvdir B
if `"`csvdir'"' == "" local csvdir "tools/data"
if `"`B'"' == "" local B 0.999
scalar dphi_B = `B'

capture noisily do "src/mata/dwlsmv_cellprob.mata"
if _rc {
    display as error "measure_binormal: src/mata/dwlsmv_cellprob.mata failed to load"
    exit 9
}

display "DPHI|env|stata=" c(stata_version) "|born=" c(born_date) "|os=" c(os) ///
        "|machine=" c(machine_type) "|procs=" c(processors)

* ---- corners ---------------------------------------------------------------
import delimited using `"`csvdir'/binormal_ref_corners.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) clear

mata:
void dphi_line(string scalar scope, string scalar reg, real colvector err,
               real colvector rel, real colvector h, real colvector k,
               real colvector rho)
{
    real scalar i, eps, mr
    real colvector a, rr

    eps = epsilon(1)
    a = abs(err)
    i = dphi_argmax(a)
    rr = select(rel, !missing(rel))
    mr = (rows(rr) ? max(rr) : .)
    printf("DPHI|%s|%s|n=%g|max_abs=%21.13e|units_2^-52=%9.4f|mean_abs=%21.13e|max_rel=%21.13e|at_h=%g|at_k=%g|at_rho=%g\n",
           scope, reg, rows(err), a[i], a[i] / eps, sum(a) / rows(a), mr,
           h[i], k[i], rho[i])
}

real scalar dphi_argmax(real colvector a)
{
    real scalar i, best

    best = 1
    for (i = 2; i <= rows(a); i++) {
        if (a[i] > a[best]) {
            best = i
        }
    }
    return(best)
}

void dphi_corners()
{
    string matrix S
    string colvector kind, reg, ureg
    real colvector h, k, rho, hi, lo, cmp, err, rel, sel, ismarg
    real scalar n, i, j, eps, dphi, dP, m, B, denv

    eps = epsilon(1)
    S   = st_sdata(., ("kind", "regime", "h", "k", "rho", "ref_hi", "ref_lo"))
    n   = rows(S)
    kind = S[., 1]
    reg  = S[., 2]
    h   = strtoreal(S[., 3])
    k   = strtoreal(S[., 4])
    rho = strtoreal(S[., 5])
    hi  = strtoreal(S[., 6])
    lo  = strtoreal(S[., 7])
    ismarg = (kind :== "margin")

    // Row counts as written by make_binormal_ref.R (REF|corners=12501; 27 margins).
    if (n != 12501 | sum(ismarg) != 27) {
        errprintf("measure_binormal: corners CSV has %g rows, %g margins; expected 12501, 27\n",
                  n, sum(ismarg))
        exit(9)
    }

    // Input assertions: every row parses, kinds are known, refs lie in [0, 1].
    if (sum(ismarg :| (kind :== "corner")) != n) {
        errprintf("measure_binormal: unknown kind in corners CSV" + char(10))
        exit(9)
    }
    if (missing(h) | missing(hi) | missing(lo)) {
        errprintf("measure_binormal: unparsed h / ref in corners CSV" + char(10))
        exit(9)
    }
    if (missing(select(k, !ismarg)) | missing(select(rho, !ismarg))) {
        errprintf("measure_binormal: unparsed k / rho on a corner row" + char(10))
        exit(9)
    }
    if (sum((hi :< 0) :| (hi :> 1)) > 0) {
        errprintf("measure_binormal: reference outside [0, 1]" + char(10))
        exit(9)
    }
    if (sum(ismarg) < 1 | sum(!ismarg) < 1) {
        errprintf("measure_binormal: corners CSV lacks corner or margin rows" + char(10))
        exit(9)
    }

    cmp = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        if (ismarg[i]) {
            cmp[i] = normal(h[i])
        }
        else {
            cmp[i] = binormal(h[i], k[i], rho[i])
        }
    }
    if (missing(cmp)) {
        errprintf("measure_binormal: binormal()/normal() returned missing" + char(10))
        exit(9)
    }
    err = (cmp - hi) - lo
    // Relative error only where the reference exceeds 1e-18: below that the
    // reference's own absolute accuracy (<= 5e-22) makes it meaningless. x/0 is
    // missing, so dividing by the 0/1 mask blanks those rows.
    rel = (abs(err) :/ (hi + lo)) :/ ((hi + lo) :> 1e-18)

    ureg = uniqrows(reg)
    for (j = 1; j <= rows(ureg); j++) {
        sel = (reg :== ureg[j])
        dphi_line("corner", ureg[j], select(err, sel), select(rel, sel),
                  select(h, sel), select(k, sel), select(rho, sel))
    }
    sel = !ismarg
    dphi_line("corner", "ALL_binormal", select(err, sel), select(rel, sel),
              select(h, sel), select(k, sel), select(rho, sel))
    sel = ismarg
    dphi_line("corner", "ALL_normal", select(err, sel), select(rel, sel),
              select(h, sel), select(k, sel), select(rho, sel))

    // Budget: max over the searched range |rho| <= B (margins always in).
    // Envelope: max over the whole grid. Both chosen after measuring.
    B = st_numscalar("dphi_B")
    if (sum(abs(rho) :== B) < 2 | sum(rho :== B) < 1 | sum(rho :== -B) < 1) {
        errprintf("measure_binormal: B = %g is not a grid point on both sides\n", B)
        exit(9)
    }
    sel  = ismarg :| (abs(rho) :<= B)
    dphi = max(abs(select(err, sel)))
    denv = max(abs(err))
    dP   = 4 * dphi + 3 * eps
    printf("DPHI|delta_Phi|budget|B=%g|thresholds [-2.5, 2.5]|max_abs=%21.13e|units_2^-52=%9.4f\n",
           B, dphi, dphi / eps)
    printf("DPHI|delta_Phi|envelope|whole grid|max_abs=%21.13e|units_2^-52=%9.4f\n",
           denv, denv / eps)
    printf("DPHI|delta_P|implied=4*measured+3eps=%21.13e|assumed=%21.13e|ratio=%9.4f\n",
           dP, 7 * eps, dP / (7 * eps))
    if (dphi <= eps) {
        printf("DPHI|verdict|delta_Phi=2^-52 COVERS binormal on this grid|headroom=%9.4f\n",
               eps / dphi)
    }
    else {
        printf("DPHI|verdict|delta_Phi=2^-52 EXCEEDED on this grid|excess=%9.4f\n",
               dphi / eps)
    }
    st_numscalar("dphi_measured", dphi)
    m = sum(hi :== 0 :& lo :== 0)
    printf("DPHI|corner|zero_ref_rows=%g\n", m)
}

end

mata: dphi_corners()

* ---- cells -----------------------------------------------------------------
import delimited using `"`csvdir'/binormal_ref_cells.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) clear

mata:
void dphi_cells()
{
    string matrix S
    string colvector setn, reg, ureg
    real matrix T, P
    real colvector rho, r, s, hi, lo, cmp, err, rel, a, b, sel, pos, eh, ek
    real scalar n, i, i0, eps, dphi, dP, dP0, nusable, nlost, ngain, nneg

    eps  = epsilon(1)
    dphi = st_numscalar("dphi_measured")
    dP   = 4 * dphi + 3 * eps          // implied by this machine's binormal()
    dP0  = 7 * eps                     // the draft's value (delta_Phi = 2^-52)

    S = st_sdata(., ("set", "regime", "rho", "a1", "a2", "a3", "b1", "b2", "b3",
                     "r", "s", "ref_hi", "ref_lo"))
    n    = rows(S)
    if (n != 3333) {
        errprintf("measure_binormal: cells CSV has %g rows; expected 3333\n", n)
        exit(9)
    }
    setn = S[., 1]
    reg  = S[., 2]
    rho  = strtoreal(S[., 3])
    T    = strtoreal(S[., 4..9])
    r    = strtoreal(S[., 10])
    s    = strtoreal(S[., 11])
    hi   = strtoreal(S[., 12])
    lo   = strtoreal(S[., 13])
    if (missing(rho) | missing(r) | missing(s) | missing(hi) | missing(lo)) {
        errprintf("measure_binormal: unparsed field in cells CSV" + char(10))
        exit(9)
    }
    if (missing(T[., 1]) | missing(T[., 4])) {
        errprintf("measure_binormal: a cell row has no thresholds" + char(10))
        exit(9)
    }

    cmp = J(n, 1, .)
    i0 = 1
    while (i0 <= n) {
        // A group is a run of rows sharing (set, rho): one call, many cells.
        i = i0
        while (i < n) {
            if (setn[i + 1] != setn[i0] | rho[i + 1] != rho[i0]) {
                break
            }
            i++
        }
        a = select(T[i0, 1..3]', T[i0, 1..3]' :< .)
        b = select(T[i0, 4..6]', T[i0, 4..6]' :< .)
        P = dwlsmv_cellprob_raw(a, b, rho[i0])
        if (i - i0 + 1 != rows(P) * cols(P)) {
            errprintf("measure_binormal: cell group size mismatch" + char(10))
            exit(9)
        }
        for (; i0 <= i; i0++) {
            cmp[i0] = P[r[i0], s[i0]]
        }
    }
    if (missing(cmp)) {
        errprintf("measure_binormal: a cell was never filled" + char(10))
        exit(9)
    }

    err = (cmp - hi) - lo
    rel = (abs(err) :/ (hi + lo)) :/ ((hi + lo) :> 1e-18)
    eh = J(n, 1, .)
    ek = J(n, 1, .)

    ureg = uniqrows(reg)
    for (i = 1; i <= rows(ureg); i++) {
        sel = (reg :== ureg[i])
        dphi_line("cell", ureg[i], select(err, sel), select(rel, sel),
                  select(eh, sel), select(ek, sel), select(rho, sel))
    }
    dphi_line("cell", "ALL", err, rel, eh, ek, rho)
    // Relative error only where the cell is resolvable at all.
    pos = ((hi + lo) :> dP0)
    dphi_line("cell", "ALL_ref_gt_7eps", select(err, pos), select(rel, pos),
              select(eh, pos), select(ek, pos), select(rho, pos))

    // The draft's error model says |cell error| <= delta_P. Check it both ways.
    printf("DPHI|cell_budget|max_abs_over_implied_delta_P=%9.4f|max_abs_over_7eps=%9.4f\n",
           max(abs(err)) / dP, max(abs(err)) / dP0)

    // Usability rule (computed pi > delta_P) against the truth.
    nusable = sum(cmp :> dP)
    nlost   = sum((cmp :<= dP) :& ((hi + lo) :> 2 * dP))
    ngain   = sum((cmp :> dP) :& ((hi + lo) :<= 0))
    nneg    = sum(cmp :<= 0)
    printf("DPHI|usable|delta_P=%21.13e|n=%g|usable=%g|computed_le_0=%g|true_gt_2dP_but_unusable=%g|true_le_0_but_usable=%g\n",
           dP, n, nusable, nneg, nlost, ngain)
    nusable = sum(cmp :> dP0)
    printf("DPHI|usable|delta_P=%21.13e|n=%g|usable=%g|(draft 7eps)\n", dP0, n, nusable)
}
end

mata: dphi_cells()

display "DPHI|done"
