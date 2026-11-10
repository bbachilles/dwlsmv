*! test_mg_admissible.do -- admissibility at G > 1 (spec O15, decided
*! 2026-10-31, user): dwlsmv_mg_admissible(), from the hand-built tables of
*! tests/unit/mg_model_ref.mata (case S: two groups, every row type 1-7,
*! omega != 1 in group 2; design E configural).
*!
*! Checks (tolerances fixed BEFORE the first run):
*!   a1  E configural: each group's result is bitwise dwlsmv_admissible() on
*!       that group's own rows (omega = 1, nu = alpha = 0: the transformation
*!       is the identity)
*!   a2  S, both groups, against Module 5 (dwls_model) and mg_parts, a code
*!       path the helper does not use for these values:
*!         row 1 value = min_j 1 - omega_j^2 (Lambda Phi Lambda')_jj, <= 1e-14
*!         row 4 value = min eigenvalue of R*_g built from dwls_model's rho
*!                       block with 1 on the diagonal, <= 1e-12
*!         row 6 value = min within-indicator gap of dwls_model's thresholds,
*!                       <= 1e-14
*!   a3  a Heywood case through omega alone: group 2's omega_1 set so that
*!       1 - V_1 > 0 but 1 - omega_1^2 V_1 < 0; row 1 must flag indicator 1
*!       in group 2 and not in group 1
*!   a4  a negative omega_1 in group 2: row 6 flags indicator 1 (order of
*!       sigma, not tau); rows 1 and 4 are unchanged (omega enters squared
*!       or as D_omega R D_omega); O12's row 7 flags omega_1 there and not
*!       at the positive omegas (a1 also asserts row 7 at omega = 1)
*!   a5  dwlsmv_admissadvice() on a4's result names y1 with its scale factor
*!       and its threshold order
*!
*! Mutants, run 2026-10-31:
*!   M1 loadings not scaled by omega                     killed  a2
*!   M2 thresholds not transformed                       killed  a2, a4
*!   M3 nu enters the threshold shift with +             SURVIVES, EQUIVALENT:
*!      the shift nu_j + (Lambda kappa)_j is constant within an indicator,
*!      so it cancels in every threshold gap, row 6's only use of it; kept
*!      so that the table holds sigma_g literally (MS1)
*!   M4 absent omega read as 0                           killed  a1
*!   M5 ptypes 5-7 kept in the table                     SURVIVES, EQUIVALENT:
*!      dwlsmv_matrices() and dwlsmv_admissible() read ptypes 1-4 only
*!   M6 row 7 never flags                                killed  a4, a5
*!   M7 row 7 takes the largest omega                    killed  a4, a5
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
local load `load' tests/unit/mg_model_ref.mata
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_mg_admissible: `f' failed to load"
        exit 9
    }
}

mata:

struct dwlsmv_statmap vector ta_maps(real matrix K)
{
    struct dwlsmv_statmap vector m
    real scalar g
    m = dwlsmv_statmap(cols(K))
    for (g = 1; g <= cols(K); g++) {
        m[g] = dwlsmv_statmap_build(K[., g])
    }
    return(m)
}

// a1
real scalar ta_a1()
{
    struct dwlsmv_parmap scalar par, sp
    struct dwlsmv_statmap vector maps
    real colvector th, sig, idx, fr
    real matrix K, Dl, A, A1
    real scalar g, i, n, nf

    nf = 0
    mgref_E_configural(par, th, K, sig, Dl)
    maps = ta_maps(K)
    for (g = 1; g <= 2; g++) {
        sp = dwlsmv_mg_sub(par, g, idx)
        fr = J(0, 1, .)
        n  = 0
        for (i = 1; i <= sp.nall; i++) {
            if (sp.free[i]) {
                n++
                fr = fr \ sp.fidx[i]
                sp.fidx[i] = n
            }
        }
        sp.npar = n
        A  = dwlsmv_mg_admissible(th, par, g, maps[g])
        A1 = dwlsmv_admissible(th[fr], sp, maps[g])
        if (A[1::6, .] != A1 | A[7, .] != (0, 1, 1, 0)) {
            printf("{err}a1: group %f is not bitwise the single-group result\n", g)
            A, A1
            nf++
        }
    }
    return(nf)
}

// the three reference values of a2 for group g
void ta_ref(real colvector th, struct dwlsmv_parmap scalar par,
            struct dwlsmv_statmap vector maps, real scalar g,
            real scalar r1, real scalar r4, real scalar r6, real colvector V)
{
    struct dwlsmv_parmap scalar sp
    struct dwlsmv_statmap scalar m
    real matrix Lam, Psi, B, Ai, R
    real colvector idx, v, nu, om, al, kap, s, sg
    real scalar off, k, j, c, d

    m   = maps[g]
    sp  = dwlsmv_mg_sub(par, g, idx)
    v   = dwlsmv_parvals(th, par)[idx]
    dwlsmv_mg_parts(v, sp, m.p, dwlsmv_nfac(sp), Lam, Psi, B, Ai, nu, om, al, kap)
    V   = diagonal(Lam * Ai * Psi * Ai' * Lam')
    r1  = min(1 :- (om :^ 2) :* V)

    s   = dwls_model(th, par, maps)
    off = (g == 1 ? 0 : maps[1].D)
    sg  = s[(off + 1)::(off + m.D)]
    R   = I(m.p)
    for (k = m.T + 1; k <= m.D; k++) {
        R[m.var_i[k], m.var_j[k]] = sg[k]
        R[m.var_j[k], m.var_i[k]] = sg[k]
    }
    r4 = min(symeigenvalues(R))

    r6 = .
    for (j = 1; j <= m.p; j++) {
        for (c = 2; c <= m.C[j]; c++) {
            d = sg[m.Coff[j] + c] - sg[m.Coff[j] + c - 1]
            if (d < r6) {
                r6 = d
            }
        }
    }
}

// a2
real scalar ta_a2()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    real colvector th, sig, V
    real matrix K, Dl, A
    real scalar g, r1, r4, r6, nf

    nf = 0
    mgref_S(par, th, K, sig, Dl)
    maps = ta_maps(K)
    for (g = 1; g <= 2; g++) {
        A = dwlsmv_mg_admissible(th, par, g, maps[g])
        ta_ref(th, par, maps, g, r1, r4, r6, V)
        printf("{txt}  S group %f: row1 %9.2e  row4 %9.2e  row6 %9.2e\n", g,
               abs(A[1, 2] - r1), abs(A[4, 2] - r4), abs(A[6, 2] - r6))
        if (!(abs(A[1, 2] - r1) <= 1e-14)) {
            printf("{err}a2: group %f row 1 %g, reference %g\n", g, A[1, 2], r1)
            nf++
        }
        if (!(abs(A[4, 2] - r4) <= 1e-12)) {
            printf("{err}a2: group %f row 4 %g, reference %g\n", g, A[4, 2], r4)
            nf++
        }
        if (!(abs(A[6, 2] - r6) <= 1e-14)) {
            printf("{err}a2: group %f row 6 %g, reference %g\n", g, A[6, 2], r6)
            nf++
        }
        if (any(A[., 1] :== 1)) {
            printf("{err}a2: group %f of case S flagged inadmissible\n", g)
            nf++
        }
    }
    return(nf)
}

// a3, a4. Group 2's omega_1 is theta[47] (mgref_S: row "omega g2 1 1").
real scalar ta_a34()
{
    struct dwlsmv_parmap scalar par
    struct dwlsmv_statmap vector maps
    real colvector th, th2, sig, V
    real matrix K, Dl, A0, A, A1
    real scalar r1, r4, r6, nf
    string scalar msg

    nf = 0
    mgref_S(par, th, K, sig, Dl)
    maps = ta_maps(K)
    if (par.name[47] != "omega g2 1 1" | par.fidx[47] != 47) {
        printf("{err}a3: theta[47] is not group 2's omega_1; the fixture changed\n")
        return(1)
    }

    // a3
    ta_ref(th, par, maps, 2, r1, r4, r6, V)
    if (!(V[1] > 0 & V[1] < 1)) {
        printf("{err}a3: V_1 = %g, so the case cannot separate Theta* from 1 - V\n", V[1])
        return(1)
    }
    th2     = th
    th2[47] = 1.05 / sqrt(V[1])
    A  = dwlsmv_mg_admissible(th2, par, 2, maps[2])
    A1 = dwlsmv_mg_admissible(th2, par, 1, maps[1])
    printf("{txt}  a3: 1 - V_1 = %g, row 1 = %g at indicator %f\n", 1 - V[1], A[1, 2], A[1, 3])
    if (A[1, 1] != 1 | A[1, 3] != 1) {
        printf("{err}a3: the Heywood case through omega is not flagged on indicator 1\n")
        nf++
    }
    if (!(abs(A[1, 2] - (1 - th2[47]^2 * V[1])) <= 1e-14)) {
        printf("{err}a3: row 1 value %g, expected %g\n", A[1, 2], 1 - th2[47]^2 * V[1])
        nf++
    }
    if (A1[1, 1] != 0) {
        printf("{err}a3: group 1 flagged by group 2's omega\n")
        nf++
    }

    // a4
    A0      = dwlsmv_mg_admissible(th, par, 2, maps[2])
    th2     = th
    th2[47] = -th[47]
    A  = dwlsmv_mg_admissible(th2, par, 2, maps[2])
    if (A[6, 1] != 1 | A[6, 3] != 1) {
        printf("{err}a4: a negative omega_1 does not flag indicator 1's threshold order\n")
        nf++
    }
    if (!(abs(A[1, 2] - A0[1, 2]) <= 1e-15 & abs(A[4, 2] - A0[4, 2]) <= 1e-12)) {
        printf("{err}a4: rows 1 or 4 moved with the sign of omega_1\n")
        nf++
    }
    // O12's row 7: flagged at the negative omega, not at the positive ones
    if (A0[7, 1] != 0 | A[7, .] != (1, th2[47], 1, 0)) {
        printf("{err}a4: row 7 (omega <= 0) wrong: before %g, after %g at %g\n",
               A0[7, 1], A[7, 2], A[7, 3])
        nf++
    }
    // a5: the note names the indicator and the value
    msg = dwlsmv_admissadvice(A, ("y1", "y2", "y3", "y4", "y5", "y6")', ("F1", "F2")')
    if (!strpos(msg, "y1: scale factor -0.8900, not above 0") |
        !strpos(msg, "y1: thresholds are not in increasing order")) {
        printf("{err}a5: the note does not report omega_1 and the threshold order\n")
        printf("{txt}%s\n", msg)
        nf++
    }
    return(nf)
}

real scalar ta_all()
{
    real scalar nf
    nf = ta_a1()
    nf = nf + ta_a2()
    nf = nf + ta_a34()
    return(nf)
}

end

mata: st_local("nf", strofreal(ta_all()))
if `nf' {
    display as error "test_mg_admissible: `nf' failure(s)"
    exit 9
}
display as text "test_mg_admissible: all checks passed"
