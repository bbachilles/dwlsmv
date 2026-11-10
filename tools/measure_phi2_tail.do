*! measure_phi2_tail.do -- the Phi2 tail-validation study: Arms A, B, C and the
*! univariate add-on (docs/notes/phi2-tail/prereg.md).
*!
*! Executed on the frozen grid 2026-10-29 (docs/notes/phi2-tail/result.md).
*!
*! Usage: do tools/measure_phi2_tail.do csvdir [smoke]
*!   csvdir  directory holding the phi2tail_*.csv files of make_phi2_tail_ref.R
*!   smoke   1 = smoke mode: checks the pipeline and prints NO error figures
*!
*! A measurement tool, not a test: the suite runner globs tests/unit/test_*.do
*! only. It asserts its INPUT (files, row counts, parse) and exits 9 if that is
*! broken. The study's outcome is printed as PT| lines and judged by the
*! preregistered rule.
*!
*! What is measured:
*!   C0  production: binormal() and normal() corners, differenced as in
*!       dwlsmv_cellprob_raw(); the density as in dwlsmv_dcellprob_rho.
*!   C1  candidate: the same, after reflecting each axis whose cell lies mostly
*!       above 0 (x -> -x, rho -> -rho): the corners differenced are small.
*!   The ratio terms u, du/drho, du/dtheta are a TEST-ONLY Mata transcription
*!   of cxh_pair_d (@sec-condx-hess). Production code for conditional mode does
*!   not exist; when it does, it is tested against this transcription.
*!
*! Constants (prereg §3; user decisions 2026-10-25): L = 4, B = 0.999,
*! delta_Phi budget 4.23 * 2^-52, Arm B domain P_ref >= 1e3 * delta_P, Arm C
*! target 1e-8. Arm B tests |err| <= tau_U2 (@sec-condx-hess-2nd), with no
*! allowance (user, 2026-10-28).

version 16
clear all
set more off
capture set processors 1

args csvdir smoke
if `"`csvdir'"' == "" {
    display as error "measure_phi2_tail: give the directory of the reference CSVs"
    exit 9
}
if `"`smoke'"' == "" local smoke 0

display "PT|env|stata=" c(stata_version) "|born=" c(born_date) "|os=" c(os) ///
        "|machine=" c(machine_type) "|procs=" c(processors) "|smoke=`smoke'"

mata:
// Corner with flags: f = -1 for a -inf bound, +1 for +inf, 0 finite.
real scalar pt_phi2(real scalar h, real scalar hf, real scalar k,
                    real scalar kf, real scalar rho)
{
    if (hf == -1 | kf == -1) {
        return(0)
    }
    if (hf == 1 & kf == 1) {
        return(1)
    }
    if (hf == 1) {
        return(normal(k))
    }
    if (kf == 1) {
        return(normal(h))
    }
    return(binormal(h, k, rho))
}

// Phi((y - rho x) / s), x finite
real scalar pt_cond(real scalar x, real scalar y, real scalar yf,
                    real scalar rho, real scalar s)
{
    if (yf == 1) {
        return(1)
    }
    if (yf == -1) {
        return(0)
    }
    return(normal((y - rho * x) / s))
}

// Gradient of a bound in theta = (tau_1..tau_{K-1}, b), as cxh_grad.
real rowvector pt_grad(real scalar cat, real scalar K, real scalar x)
{
    real rowvector g

    g = J(1, K, 0)
    if (cat >= 1 & cat <= K - 1) {
        g[cat] = 1
    }
    g[K] = -x
    return(g)
}

// Bounds of category cat of an item with theta th, at covariate x.
void pt_bounds(real rowvector th, real scalar cat, real scalar x,
               real rowvector a, real rowvector af, real matrix A)
{
    real scalar K

    K = cols(th)
    a = (0, 0)
    af = (0, 0)
    if (cat == 0) {
        af[1] = -1
    }
    else {
        a[1] = th[cat] - th[K] * x
    }
    if (cat == K - 1) {
        af[2] = 1
    }
    else {
        a[2] = th[cat + 1] - th[K] * x
    }
    A = pt_grad(cat, K, x) \ pt_grad(cat + 1, K, x)
}

// The closed forms of @sec-condx-hess for one cell: (P, u, v, w).
real rowvector pt_cell(real rowvector a, real rowvector af, real rowvector b,
                       real rowvector bf, real scalar rho, real matrix A,
                       real matrix B)
{
    real scalar s2, s, P, Pr, Prr, k, l, sg, x, y, f, Q, u
    real rowvector dPi, dPri, dPj, dPrj

    s2 = 1 - rho * rho
    s = sqrt(s2)
    P = 0
    Pr = 0
    Prr = 0
    dPi = J(1, cols(A), 0)
    dPri = dPi
    dPj = J(1, cols(B), 0)
    dPrj = dPj
    for (k = 1; k <= 2; k++) {
        for (l = 1; l <= 2; l++) {
            sg = (-1)^(k + l)
            P = P + sg * pt_phi2(a[k], af[k], b[l], bf[l], rho)
            if (af[k] == 0 & bf[l] == 0) {
                x = a[k]
                y = b[l]
                f = exp(-(x * x - 2 * rho * x * y + y * y) /
                        (2 * (1 - rho * rho))) / (2 * pi() * s)
                Q = (x * x - 2 * rho * x * y + y * y) / s2
                Pr = Pr + sg * f
                Prr = Prr + sg * f * (rho * (1 - Q) + x * y) / s2
                dPri = dPri - sg * f * (x - rho * y) / s2 * A[k, .]
                dPrj = dPrj - sg * f * (y - rho * x) / s2 * B[l, .]
            }
            if (af[k] == 0) {
                dPi = dPi + sg * normalden(a[k]) *
                      pt_cond(a[k], b[l], bf[l], rho, s) * A[k, .]
            }
            if (bf[l] == 0) {
                dPj = dPj + sg * normalden(b[l]) *
                      pt_cond(b[l], a[k], af[k], rho, s) * B[l, .]
            }
        }
    }
    u = Pr / P
    return((P, u, Prr / P - u^2, (dPri, dPrj) / P - u * (dPi, dPj) / P))
}

// C1: reflect an axis whose cell lies mostly above 0, i.e. a1 + a2 > 0
// (a +inf upper bound counts as +inf, a -inf lower bound as -inf).
real scalar pt_reflect(real rowvector a, real rowvector af, real matrix A)
{
    real scalar up

    if (af[1] == -1) {
        up = 0
    }
    else if (af[2] == 1) {
        up = 1
    }
    else {
        up = (a[1] + a[2] > 0)
    }
    if (!up) {
        return(1)
    }
    a = (-a[2], -a[1])
    af = (-af[2], -af[1])
    A = (-A[2, .]) \ (-A[1, .])
    return(-1)
}
real rowvector pt_cell_c1(real rowvector a0, real rowvector af0,
                          real rowvector b0, real rowvector bf0,
                          real scalar rho, real matrix A0, real matrix B0)
{
    real rowvector a, af, b, bf, r
    real matrix A, B
    real scalar sg

    a = a0
    af = af0
    b = b0
    bf = bf0
    A = A0
    B = B0
    sg = pt_reflect(a, af, A) * pt_reflect(b, bf, B)
    r = pt_cell(a, af, b, bf, sg * rho, A, B)
    r[2] = sg * r[2]
    r[|4 \ cols(r)|] = sg * r[|4 \ cols(r)|]
    return(r)
}

real matrix pt_num(string rowvector v)
{
    return(strtoreal(st_sdata(., v)))
}

void pt_require(real scalar ok, string scalar what)
{
    if (!ok) {
        errprintf("measure_phi2_tail: input check failed: %s\n", what)
        exit(9)
    }
}

// tau_U2 (@sec-condx-hess-2nd): first-order tau_U plus the exact remainder,
// formed from the reference's columns as tools/make_phi2_tail_ref.R (sec_cols)
// documents. Missing where eta-bar >= 1 (the bound is undefined there).
real rowvector pt_tau2(real rowvector TT, real rowvector SS, real scalar du)
{
    real rowvector tau, pw, wz
    real scalar eb, pu, pv, d

    tau = (TT[1] + du * TT[2], TT[3] + du * TT[4],
           TT[|5 \ 12|] + du * TT[|13 \ 20|])
    eb = SS[1] + du * SS[2]
    if (eb >= 1 | eb < 0) {
        return(J(1, 10, .))
    }
    d = (1 - eb)^2
    pu = SS[3] + du * SS[4]
    pv = SS[5] + du * SS[6]
    pw = SS[|11 \ 18|] + du * SS[|19 \ 26|]
    wz = SS[|43 \ 50|] + du * SS[|51 \ 58|]
    tau[1] = tau[1] + pu * eb / (1 - eb)
    tau[2] = tau[2] + (pv + SS[7] + eb * SS[8] + SS[9] + eb^2 * SS[10]) / d - pv
    tau[|3 \ 10|] = tau[|3 \ 10|] +
        (pw + SS[|27 \ 34|] + eb * SS[|35 \ 42|] + wz + eb^2 * SS[|59 \ 66|]) :/ d - pw
    return(tau)
}

// Arm B: per cell and representation, max over u, v, w_k of |err| / tau,
// tau = tau_U2 (missing where the evaluation is not finite; -1 where tau_U2
// is undefined, eta-bar >= 1, which never occurs inside Arm B).
real matrix pt_armB(real matrix C, real matrix R, real matrix T0,
                    real matrix T1, real matrix S0, real matrix S1,
                    string colvector tmpl, real scalar du)
{
    real matrix ratio, A, B, TT, SS
    real rowvector thi, thj, a, af, b, bf, r, tau, e, tq, rr
    real scalar i, c, q, nq, h

    ratio = J(rows(C), 2, 0)
    for (i = 1; i <= rows(C); i++) {
        thi = (strtoreal(tokens(tmpl[C[i, 1]])), C[i, 4])
        thj = (strtoreal(tokens(tmpl[C[i, 2]])), C[i, 5])
        pt_bounds(thi, C[i, 7], C[i, 3], a, af, A)
        pt_bounds(thj, C[i, 8], C[i, 3], b, bf, B)
        nq = 2 + cols(thi) + cols(thj)
        for (c = 1; c <= 2; c++) {
            if (c == 1) {
                r = pt_cell(a, af, b, bf, C[i, 6], A, B)
                TT = T0[i, .]
                SS = S0[i, .]
            }
            else {
                r = pt_cell_c1(a, af, b, bf, C[i, 6], A, B)
                TT = T1[i, .]
                SS = S1[i, .]
            }
            tau = pt_tau2(TT, SS, du)
            if (hasmissing(r)) {
                ratio[i, c] = .
            }
            else if (hasmissing(tau)) {
                ratio[i, c] = -1
            }
            else {
                e = J(1, nq, 0)
                for (q = 1; q <= nq; q++) {
                    e[q] = abs((r[q + 1] - R[i, 2 * q + 1]) - R[i, 2 * q + 2])
                }
                // into the cell's units of 2^k, as tau: exact power-of-2 steps
                // of at most 2^1000, since 2^-k (and its half) can overflow
                // (k <= -2181 on the grid; changed after the frozen run,
                // 2026-10-29, prereg §6)
                h = -SS[67]
                while (abs(h) > 1000) {
                    e = e :* 2^(sign(h) * 1000)
                    h = h - sign(h) * 1000
                }
                e = e :* 2^h
                tq = tau[|1 \ nq|]
                rr = (e :/ (tq + (tq :== 0))) :* (tq :> 0) + (tq :== 0) :* (e :> 0) * 1e300
                ratio[i, c] = max(rr)
            }
        }
    }
    return(ratio)
}

// Arm D: per cell with D_dom 1 or 2, the relative error of P under C0 and C1
// against D_P (hi + lo); missing where the evaluation is not finite, and for
// cells outside those classes.
real matrix pt_armD(real matrix C, real matrix DD, string colvector tmpl)
{
    real matrix rel, A, B
    real rowvector thi, thj, a, af, b, bf, r
    real scalar i, c

    rel = J(rows(C), 2, .)
    for (i = 1; i <= rows(C); i++) {
        if (DD[i, 3] == 1 | DD[i, 3] == 2) {
            thi = (strtoreal(tokens(tmpl[C[i, 1]])), C[i, 4])
            thj = (strtoreal(tokens(tmpl[C[i, 2]])), C[i, 5])
            pt_bounds(thi, C[i, 7], C[i, 3], a, af, A)
            pt_bounds(thj, C[i, 8], C[i, 3], b, bf, B)
            for (c = 1; c <= 2; c++) {
                if (c == 1) {
                    r = pt_cell(a, af, b, bf, C[i, 6], A, B)
                }
                else {
                    r = pt_cell_c1(a, af, b, bf, C[i, 6], A, B)
                }
                if (!missing(r[1])) {
                    rel[i, c] = abs((r[1] - DD[i, 1]) - DD[i, 2]) / DD[i, 1]
                }
            }
        }
    }
    return(rel)
}

// Arm C: relative error of U_rho and normwise of U_theta, per block, C0 and C1.
real matrix pt_armC(real matrix D, real matrix TH, real matrix K)
{
    real matrix relC, A, B
    real colvector sel
    real rowvector thi, thj, a, af, b, bf, r, Ut, reft, et
    real scalar m, c, z, i, q, Ur, refr

    relC = J(rows(K), 4, .)
    for (m = 1; m <= rows(K); m++) {
        sel = selectindex(D[., 1] :== K[m, 1] :& D[., 2] :== K[m, 2])
        pt_require(rows(sel) > 0, "block without cells")
        for (c = 1; c <= 2; c++) {
            Ur = 0
            Ut = J(1, 8, 0)
            for (z = 1; z <= rows(sel); z++) {
                i = sel[z]
                thi = TH[i, (1..D[i, 7])]
                thj = TH[i, (5..(4 + D[i, 8]))]
                pt_bounds(thi, D[i, 3], D[i, 5], a, af, A)
                pt_bounds(thj, D[i, 4], D[i, 5], b, bf, B)
                if (c == 1) {
                    r = pt_cell(a, af, b, bf, D[i, 9], A, B)
                }
                else {
                    r = pt_cell_c1(a, af, b, bf, D[i, 9], A, B)
                }
                Ur = Ur + D[i, 6] * r[3]
                Ut[|1 \ cols(r) - 3|] = Ut[|1 \ cols(r) - 3|] + D[i, 6] * r[|4 \ cols(r)|]
            }
            refr = K[m, 3] + K[m, 4]
            reft = J(1, 8, 0)
            et = J(1, 8, 0)
            for (q = 1; q <= 8; q++) {
                reft[q] = K[m, 3 + 2 * q] + K[m, 4 + 2 * q]
                et[q] = abs((Ut[q] - K[m, 3 + 2 * q]) - K[m, 4 + 2 * q])
            }
            relC[m, 2 * c - 1] = abs((Ur - K[m, 3]) - K[m, 4]) / abs(refr)
            relC[m, 2 * c] = max(et) / max(abs(reft))
        }
    }
    return(relC)
}
end

* ---- Arm A: corners --------------------------------------------------------
import delimited using `"`csvdir'/phi2tail_corners.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) case(preserve) clear
mata:
X = pt_num(("h", "k", "rho", "cert", "ref_hi", "ref_lo"))
pt_require(rows(X) > 0 & !hasmissing(X), "corners")
err = J(rows(X), 1, .)
for (i = 1; i <= rows(X); i++) {
    err[i] = (binormal(X[i, 1], X[i, 2], X[i, 3]) - X[i, 5]) - X[i, 6]
}
pt_require(!hasmissing(err), "corner evaluation")
eps52 = epsilon(1)
cert = X[., 4]
A_cert = max(abs(select(err, cert)))
A_env  = (sum(!cert) ? max(abs(select(err, !cert))) : 0)
end

* ---- univariate add-on (also normal()'s margins for Arm A) -----------------
import delimited using `"`csvdir'/phi2tail_uni.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) case(preserve) clear
mata:
ls = st_sdata(., "l")
us = st_sdata(., "u")
U = pt_num(("cert", "P_hi", "P_lo", "rl_hi", "rl_lo", "ru_hi", "ru_lo"))
pt_require(rows(U) > 0 & !hasmissing(U), "univariate rows")
lf = -(ls :== "-inf")
uf = (us :== "inf")
lv = strtoreal(ls)
uv = strtoreal(us)
pt_require(!hasmissing(select(lv, !lf)) & !hasmissing(select(uv, !uf)), "univariate bounds")
uerr = J(rows(U), 3, 0)
merr = J(0, 1, .)
for (i = 1; i <= rows(U); i++) {
    if (lf[i] == 0 & lv[i] > 0) {
        P = normal(-lv[i]) - (uf[i] ? 0 : normal(-uv[i]))
    }
    else {
        P = (uf[i] ? 1 : normal(uv[i])) - (lf[i] ? 0 : normal(lv[i]))
    }
    rl = (lf[i] ? 0 : normalden(lv[i]) / P)
    ru = (uf[i] ? 0 : normalden(uv[i]) / P)
    uerr[i, 1] = abs(((P - U[i, 2]) - U[i, 3]) / U[i, 2])
    uerr[i, 2] = (lf[i] ? 0 : abs(((rl - U[i, 4]) - U[i, 5]) / U[i, 4]))
    uerr[i, 3] = (uf[i] ? 0 : abs(((ru - U[i, 6]) - U[i, 7]) / U[i, 6]))
    if (lf[i] & !uf[i]) {
        merr = merr \ ((normal(uv[i]) - U[i, 2]) - U[i, 3])
    }
}
pt_require(!hasmissing(uerr), "univariate evaluation")
A_norm = max(abs(merr))
ucert = U[., 1]
Uni_cert = max(select(uerr, ucert))
Uni_env = (sum(!ucert) ? max(select(uerr, !ucert)) : 0)

// delta_Phi for Arms B and C: the budget if Arm A passes, else the measured
// certification maximum (prereg §5), in units of 2^-52.
A_pass = (max((A_cert, A_norm)) <= 4.23 * eps52)
dphi_units = (A_pass ? 4.23 : max((A_cert, A_norm)) / eps52)
deltaP = 4 * dphi_units * eps52 + 3 * eps52
end

* ---- Arm B: the cell design ------------------------------------------------
import delimited using `"`csvdir'/phi2tail_templates.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) case(preserve) clear
mata:
tmpl = st_sdata(., "tau")
end
import delimited using `"`csvdir'/phi2tail_cells.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) case(preserve) clear
mata:
C = pt_num(("ti", "tj", "x", "bi", "bj", "rho", "ci", "cj", "cert"))
R = pt_num(("P_hi", "P_lo", "u_hi", "u_lo", "v_hi", "v_lo",
            "w1_hi", "w1_lo", "w2_hi", "w2_lo", "w3_hi", "w3_lo", "w4_hi", "w4_lo",
            "w5_hi", "w5_lo", "w6_hi", "w6_lo", "w7_hi", "w7_lo", "w8_hi", "w8_lo"))
T0 = pt_num(("c0_u_t0", "c0_u_t1", "c0_v_t0", "c0_v_t1",
             "c0_w1_t0", "c0_w2_t0", "c0_w3_t0", "c0_w4_t0", "c0_w5_t0", "c0_w6_t0", "c0_w7_t0", "c0_w8_t0",
             "c0_w1_t1", "c0_w2_t1", "c0_w3_t1", "c0_w4_t1", "c0_w5_t1", "c0_w6_t1", "c0_w7_t1", "c0_w8_t1"))
T1 = pt_num(("c1_u_t0", "c1_u_t1", "c1_v_t0", "c1_v_t1",
             "c1_w1_t0", "c1_w2_t0", "c1_w3_t0", "c1_w4_t0", "c1_w5_t0", "c1_w6_t0", "c1_w7_t0", "c1_w8_t0",
             "c1_w1_t1", "c1_w2_t1", "c1_w3_t1", "c1_w4_t1", "c1_w5_t1", "c1_w6_t1", "c1_w7_t1", "c1_w8_t1"))
// A cell whose reference P underflows double (P_hi == 0) has no defined u, v or
// w (0/0), and the reference leaves them blank. Such a cell is outside Arm B
// (P >= 1e3 delta_P) and must not be Arm D certified (checked below). Every
// other cell must be complete. (Changed after the frozen run stopped here on
// input, 2026-10-29, before any figure was read; prereg §6.)
underP = (R[., 1] :== 0)
pt_require(rows(C) > 0 & !hasmissing(C) & !hasmissing(T0) & !hasmissing(T1), "cells")
pt_require(!hasmissing(select(R, !underP)), "cells")
printf("PT|cells|P_underflows_double=%g|of=%g\n", sum(underP), rows(C))
// tau_U2's pieces, in make_phi2_tail_ref.R's sec_cols order (66 per representation),
// then k: tau_U and these pieces are in the cell's units of 2^k
sn = ("eb_t0", "eb_t1", "pu_t0", "pu_t1", "pv_t0", "pv_t1", "vx", "vy", "vz", "va")
for (q = 1; q <= 8; q++) sn = sn, "pw" + strofreal(q) + "_t0"
for (q = 1; q <= 8; q++) sn = sn, "pw" + strofreal(q) + "_t1"
for (q = 1; q <= 8; q++) sn = sn, "wx" + strofreal(q)
for (q = 1; q <= 8; q++) sn = sn, "wy" + strofreal(q)
for (q = 1; q <= 8; q++) sn = sn, "wz" + strofreal(q) + "_t0"
for (q = 1; q <= 8; q++) sn = sn, "wz" + strofreal(q) + "_t1"
for (q = 1; q <= 8; q++) sn = sn, "wa" + strofreal(q)
sn = sn, "k"
S0 = pt_num("c0_" :+ sn)
S1 = pt_num("c1_" :+ sn)
pt_require(cols(S0) == 67 & !hasmissing(S0) & !hasmissing(S1), "tau_U2 columns")
ratio = pt_armB(C, R, T0, T1, S0, S1, tmpl, dphi_units)
nonfin = sum(missing(ratio))
inB = C[., 9] :& (R[., 1] :>= 1e3 * deltaP)
// tau = -1 marks "not defined" (double P underflows); never inside Arm B
pt_require(!sum(inB :& (rowmin(T0) :< 0 :| rowmin(T1) :< 0)), "tau undefined inside Arm B")
pt_require(!sum(inB :& (rowmin(S0[|1, 1 \ rows(S0), 66|]) :< 0 :| rowmin(S1[|1, 1 \ rows(S1), 66|]) :< 0)), "tau_U2 pieces undefined inside Arm B")
pt_require(!sum(inB :& (ratio[., 1] :< 0 :| ratio[., 2] :< 0)), "tau_U2 undefined inside Arm B")
// Arm D: C1's relative error of P against the M4-based reference D_P, on the
// tail cells; D_dom 1 certified, 2 diagnostic, 3 floor, 0 not in Arm D
DD = pt_num(("D_P_hi", "D_P_lo", "D_dom"))
pt_require(!hasmissing(DD), "Arm D columns")
pt_require(!sum(underP :& (DD[., 3] :== 1)), "underflowed P inside Arm D certified")
relD = pt_armD(C, DD, tmpl)
end

* ---- Arm C: the design datasets --------------------------------------------
import delimited using `"`csvdir'/phi2tail_design.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) case(preserve) clear
mata:
D = pt_num(("dataset", "pair", "ci", "cj", "x", "w", "Ki", "Kj", "rho"))
TH = strtoreal(st_sdata(., ("thi1", "thi2", "thi3", "thi4", "thj1", "thj2", "thj3", "thj4")))
pt_require(rows(D) > 0 & !hasmissing(D), "design cells")
end
import delimited using `"`csvdir'/phi2tail_blocks.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) case(preserve) clear
mata:
K = pt_num(("dataset", "pair", "Ur_hi", "Ur_lo",
            "Ut1_hi", "Ut1_lo", "Ut2_hi", "Ut2_lo", "Ut3_hi", "Ut3_lo", "Ut4_hi", "Ut4_lo",
            "Ut5_hi", "Ut5_lo", "Ut6_hi", "Ut6_lo", "Ut7_hi", "Ut7_lo", "Ut8_hi", "Ut8_lo"))
pt_require(rows(K) > 0 & !hasmissing(K), "blocks")
relC = pt_armC(D, TH, K)
end

* ---- report ----------------------------------------------------------------
mata:
if (`smoke') {
    printf("PT|smoke|corners=%g|uni=%g|cells=%g|cells_in_B=%g|nonfinite=%g|nonfinite_in_B=%g|design=%g|blocks=%g|missing_relC=%g\n",
           rows(X), rows(U), rows(C), sum(inB), nonfin,
           sum(inB :& rowmissing(ratio)), rows(D), rows(K), sum(missing(relC)))
    printf("PT|smoke|D|certified=%g|diagnostic=%g|floor=%g|nonfinite_C1_in_D=%g\n",
           sum(DD[., 3] :== 1), sum(DD[., 3] :== 2), sum(DD[., 3] :== 3),
           sum(((DD[., 3] :== 1) :| (DD[., 3] :== 2)) :& rowmissing(relD[., 2])))
    printf("PT|smoke|no error figures printed in smoke mode\n")
}
else {
    printf("PT|A|binormal_cert_units=%9.4f|normal_units=%9.4f|env_units=%9.4f|budget=4.23|pass=%g\n",
           A_cert / eps52, A_norm / eps52, A_env / eps52, A_pass)
    printf("PT|A|dphi_units_used=%9.4f|deltaP=%21.13e\n", dphi_units, deltaP)
    for (c = 1; c <= 2; c++) {
        rc = select(ratio[., c], inB)
        printf("PT|B|C%g|cells=%g|nonfinite=%g|max_ratio=%12.6f|violations=%g|bound=tau_U2|pass=%g\n",
               c - 1, rows(rc), sum(missing(rc)), max(rc), sum(rc :> 1),
               (sum(rc :> 1) == 0 & sum(missing(rc)) == 0))
        re = select(ratio[., c], !inB)
        if (rows(re)) {
            printf("PT|B|C%g|envelope_cells=%g|tau2_undefined=%g|max_ratio=%12.6e\n",
                   c - 1, rows(re), sum(re :< 0), max(re))
        }
    }
    // Arm D: C1 decides on the certified cells; C0, and the diagnostic cells,
    // are reported only. max() skips missing, so non-finite cells are counted.
    dc = (DD[., 3] :== 1)
    dg = (DD[., 3] :== 2)
    rc = select(relD[., 2], dc)
    printf("PT|D|C1|certified_cells=%g|nonfinite=%g|max_rel=%12.4e|violations=%g|target=1e-8|pass=%g\n",
           rows(rc), sum(missing(rc)), max(rc), sum(rc :> 1e-8),
           (rows(rc) > 0 & sum(rc :> 1e-8) == 0 & sum(missing(rc)) == 0))
    re = select(relD[., 2], dg)
    printf("PT|D|C1|diagnostic_cells=%g|nonfinite=%g|max_rel=%12.4e|over_target=%g|not certified\n",
           rows(re), sum(missing(re)), max(re), sum(re :> 1e-8))
    re = select(relD[., 1], dc :| dg)
    printf("PT|D|C0|cells=%g|nonfinite=%g|max_rel=%12.4e|over_target=%g|diagnostic\n",
           rows(re), sum(missing(re)), max(re), sum(re :> 1e-8))
    printf("PT|D|floor_cells=%g|envelope\n", sum(DD[., 3] :== 3))
    for (m = 1; m <= rows(K); m++) {
        printf("PT|C|dataset=%g|pair=%g|C0_rel_r=%12.4e|C0_rel_t=%12.4e|C1_rel_r=%12.4e|C1_rel_t=%12.4e\n",
               K[m, 1], K[m, 2], relC[m, 1], relC[m, 2], relC[m, 3], relC[m, 4])
    }
    printf("PT|C|C0_max=%12.4e|C1_max=%12.4e|target=1e-8|C0_pass=%g|C1_pass=%g\n",
           max(relC[., 1..2]), max(relC[., 3..4]),
           max(relC[., 1..2]) <= 1e-8, max(relC[., 3..4]) <= 1e-8)
    printf("PT|uni|cert_max_rel=%12.4e|env_max_rel=%12.4e|target=1e-8|pass=%g\n",
           Uni_cert, Uni_env, Uni_cert <= 1e-8)
}
end
