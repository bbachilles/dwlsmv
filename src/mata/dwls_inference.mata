*! dwlsmv -- Modules 8-9: inference on a fit
*!
*! Executed under StataMP 18 on 2026-09-24; passed on the first run. tests/unit/test_inference.do: T3 end to end within 9.8e-3 of spec 7.2's test-statistic tier on all five fits, and b within 1.4e-2 of its own (check 3); every error path aborts at rc 3498 and leaves the fit unfilled, one of them after acov was computed (check 8). Adversarially reviewed 2026-09-24, 23 mutations run: see test_inference.do's header.
*!
*! The locked Section 5.1 entry point:
*!
*!     void dwls_inference(struct dwlsmv_fit scalar f)
*!
*! Fills acov, shift_a, shift_b, T3 and T3_df, and nothing else:
*!
*!     acov      dwlsmv_acov(Delta, Gamma) -- ACOV of sqrt(N)(theta-hat -
*!               theta_0), never e(V); dwlsmv_ev() is Module 11's call
*!     shift_a   a = sqrt(D / Tr(M^2))           (dwlsmv_scaleshift)
*!     shift_b   b = D - a Tr(M)
*!     T3        a T_raw + b, with T_raw = (N-1) F-hat as Module 7 stored it
*!     T3_df     D_stat. T3 keeps df = D; stored explicitly because T2 would
*!               not (spec 4.9)
*!
*! with the traces from dwlsmv_traces(). Tr(M) and Tr(M^2) are NOT stored:
*! they are (D - b)/a and D/a^2, and a stored copy could disagree.
*!
*! ALL OR NOTHING (E-038, 2026-09-23). Everything is computed into locals
*! and the fields are assigned only once every step has succeeded. f is passed
*! by reference, so assigning acov before the traces were computed would leave
*! a half-filled fit visible to any caller that ran this under capture and an
*! error was raised in between. Unfilled is rows(acov) == 0 with T3, T3_df,
*! shift_a and shift_b missing, which is what dwls_fit() returns.
*!
*! D_STAT = 0 (decided 2026-09-24, spec 4.9): acov is filled, T3 = T_raw,
*! T3_df = 0, and shift_a and shift_b are missing -- lavaan 0.7.2's behaviour.
*! That is a complete state, not a partial one: rows(acov) > 0 tells it apart
*! from an unfilled fit.
*!
*! CONVERGENCE IS NOT CHECKED (decided 2026-09-24): a fit with converged = 0 is
*! processed like any other, as Module 7 reports non-convergence and never
*! enforces it. The results then describe the last iterate, which is not an
*! optimum; whether to show them is the caller's decision, and
*! f.converged remains the record.
*!
*! GROUPS (built 2026-10-10). The fit carries one Stage 1 per group and Delta
*! stacked. dwlsmv_omega() forms Omega = blockdiag(Gamma_g / pi_g); since
*! diag(Omega)^-1 is the fit's weighting (a), the single-group primitives fed
*! (Delta, Omega) are the pooled ones: ACOV = E^-1 B E^-1 on the sqrt(N) scale
*! with N the total (MG6; e(V) = ACOV / (N - G) is Module 11's), and the traces
*! of MG4. G = 1 passes Gamma itself, so nothing changes there.
*!
*! Load order: dwlsmv_types, dwlsmv_acov, dwlsmv_traces, dwlsmv_scaleshift.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Omega = blockdiag(Gamma_g / pi_g), pi_g = N_g / N (@sec-multigroup-scaling
// MG1). diag(Omega)^-1 is weighting (a), so dwlsmv_acov() and dwlsmv_traces()
// fed (Delta, Omega) give MG6's E^-1 B E^-1 and MG4's traces. G = 1 returns
// Gamma itself.
real matrix dwlsmv_omega(struct dwlsmv_stage1 vector sts)
{
    real matrix Om
    real scalar G, g, N, o, d

    G = length(sts)
    if (G == 1) {
        return(sts[1].Gamma)
    }
    N = 0
    d = 0
    for (g = 1; g <= G; g++) {
        N = N + sts[g].N
        d = d + rows(sts[g].Gamma)
    }
    Om = J(d, d, 0)
    o  = 0
    for (g = 1; g <= G; g++) {
        d = rows(sts[g].Gamma)
        Om[(o + 1)::(o + d), (o + 1)::(o + d)] = sts[g].Gamma :/ (sts[g].N / N)
        o = o + d
    }
    return(Om)
}

void dwls_inference(struct dwlsmv_fit scalar f)
{
    real matrix A, Om
    real rowvector tr, ab
    real scalar T3, D, g

    if (length(f.stage1) < 1) {
        _error("dwls_inference: the fit carries no Stage 1 result")
    }
    D = 0
    for (g = 1; g <= length(f.stage1); g++) {
        D = D + f.stage1[g].map.D
        if (rows(f.stage1[g].Gamma) != f.stage1[g].map.D) {
            errprintf("dwls_inference: group %f's Gamma has %f row(s); its map has D = %f\n",
                      g, rows(f.stage1[g].Gamma), f.stage1[g].map.D)
            _error("dwls_inference: a group's Gamma does not match its map")
        }
    }

    if (rows(f.Delta) != D | cols(f.Delta) != f.par.npar) {
        errprintf("dwls_inference: Delta is %f x %f; the map(s) have D = %f and the table npar = %f\n",
                  rows(f.Delta), cols(f.Delta), D, f.par.npar)
        _error("dwls_inference: Delta does not match the fit's map and parameter table")
    }
    if (f.D_stat != D - f.par.npar) {
        errprintf("dwls_inference: D_stat = %f but D - npar = %f\n",
                  f.D_stat, D - f.par.npar)
        _error("dwls_inference: D_stat is not D - npar")
    }
    if (missing(f.T_raw) | f.T_raw < 0) {
        _error("dwls_inference: T_raw is missing or negative")
    }

    Om = dwlsmv_omega(f.stage1)
    A  = dwlsmv_acov(f.Delta, Om)
    tr = dwlsmv_traces(f.Delta, Om)
    ab = dwlsmv_scaleshift(tr[1], tr[2], f.D_stat)
    T3 = (f.D_stat == 0 ? f.T_raw : ab[1] * f.T_raw + ab[2])

    f.acov    = A
    f.shift_a = ab[1]
    f.shift_b = ab[2]
    f.T3      = T3
    f.T3_df   = f.D_stat
}

end
