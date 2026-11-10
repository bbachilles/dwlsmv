*! Module 12, DIFFTEST validation check 6, the comparison: dwlsmv_dtrefit() accepts
*! a refit that reproduces the stored e(b), e(T3), e(bound) and e(df_T3) and refuses
*! (rc 3498) anything else, including an unset tolerance.
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9; the T3 and df cases
*! rewritten by the main session (user, 2026-10-08). Executed 2026-10-08: 22/22
*! green, each refusal read for its intended rule; mutants below each turned it red.
*!
*! Rules (spec Module 12 check 6; user, 2026-10-08): measure |x - y| / (1 + |y|),
*! y stored, max over elements; e(bound) identical; e(df_T3) equal; a missing
*! e(T3) in either result refused (dwlsmv posts T3 = T_raw for a saturated fit).
*!
*! THE PRODUCTION TOLERANCE, dwlsmv_dt_reftol(), is 1e-10 since 2026-10-09, chosen
*! after measuring (observed maximum 8.05e-13 across processor counts; frozen rule
*! 100 x, rounded up to a power of ten; E-026). TOL = 1e-9 in c1-c22 is
*! TEST-LOCAL: it only puts constructed perturbations on each side of it, each at
*! least a factor 10 away. Cases t1-t3 test the production value.
*!
*! Constructed values: RB, 20 elements including an exact 0; RT = 1234.5678;
*! RK, a 5 x 1 0/1 vector; DF = 7.
*!   c1  identical                                          pass
*!   c2  one e(b) element within TOL (1e-12 under the measure)  pass
*!   c3  e(T3) within TOL at large magnitude (relative 1e-11,
*!       absolute 1.2e-8)                                   pass; kills an ABSOLUTE measure
*!   c4  saturated (df 0), T3 = 1e-20 stored, 3e-20 refit   pass; kills a pure ratio
*!   c5  e(b) off by 1e-6 at an element that is exactly 0   refuse; kills a pure ratio
*!   c6  ONE of 20 e(b) elements off by 1e-8                refuse; kills a MEAN
*!   c7  e(T3) off by relative 1e-6                         refuse; kills b-only
*!   c8  one e(bound) flag flipped, all else identical      refuse; kills ignoring bound
*!   c9  e(T3) missing in the refit only                    refuse
*!   c10 e(T3) missing in both                              refuse
*!   c11 tol missing (the former placeholder's value)      refuse
*!   c12 tol = 0                                            refuse
*!   c13 tol < 0                                            refuse
*!   c14 e(b) lengths differ                                refuse
*!   c15 e(bound) lengths differ                            refuse
*!   c16 e(b) missing in the refit                          refuse
*!   c17 e(b) missing in the stored result                  refuse
*!   c18 e(bound) missing                                   refuse
*!   c19 empty e(b)                                         refuse
*!   c20 e(df_T3) 7 stored, 8 refit, all else identical     refuse
*!   c21 e(df_T3) missing in both                           refuse
*!   c22 saturated, T3 0 stored and 1e-6 refit              refuse
*!   t1  dwlsmv_dt_reftol() is exactly 1e-10 (the measured, labelled value)
*!   t2  at the production tolerance, e(T3) off by the observed maximum
*!       (relative 8.1e-13)                                 pass
*!   t3  at the production tolerance, e(T3) off by relative 1e-9   refuse
*!
*! MUTANTS, each of which must turn this file red (killed by):
*!   R1  only e(b) compared                          c7, c22
*!   R2  only e(T3) compared                         c5, c6
*!   R3  e(bound) ignored                            c8
*!   R4  absolute difference instead of the measure  c3
*!   R5  mean over elements (sum / n) instead of max c6
*!   R6  pure ratio |x - y| / |y|                    c5
*!   R7  missing tol accepted                        c11
*!   R8  missing T3 read as agreement                c9, c10
*!   R8c both-missing T3 as agreement (draft C rule) c10
*!   (Removing the missing-T3 guard alone is rc-equivalent: a missing T3 then
*!   fails the measure and is still refused, with a less specific message.)
*!   R9  df not compared                             c20
*!   R10 dimension checks removed                    c14, c15
*!   R11 missing e(b) guard removed                  c16, c17
*!   R12 missing df guard removed                    c21
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
foreach f of local load {
    capture noisily do "`f'"
    if _rc {
        display as error "test_dtrefit: `f' failed to load"
        exit 9
    }
}

mata:

real scalar tdr_rc(real matrix b0, real scalar T0, real matrix k0, real scalar d0,
                   real matrix b1, real scalar T1, real matrix k1, real scalar d1,
                   real scalar tol)
{
    dwlsmv_dtrefit(b0, T0, k0, d0, b1, T1, k1, d1, tol, "M")
    return(0)
}

end

mata: RB = (0.5, -1.25, 0.8, 0.1, -0.3, 1.2, 0.75, 0, 2, 0.6, 0.65, 0.7, 0.4, 0.45, 0.5, 0.55, 0.6, 0.62, 0.64, 0.66)
mata: RK = (0 \ 0 \ 1 \ 0 \ 0)
mata: RT = 1234.5678
mata: DF = 7
mata: TOL = 1e-9
mata: B_OK = RB
mata: B_OK[7] = RB[7] + 1e-12 * (1 + abs(RB[7]))
mata: B_BIG = RB
mata: B_BIG[8] = RB[8] + 1e-6
mata: B_ONE = RB
mata: B_ONE[3] = RB[3] + 1e-8 * (1 + abs(RB[3]))
mata: B_MISS = RB
mata: B_MISS[2] = .
mata: B_SHORT = RB[(1..19)]
mata: B_EMPTY = J(1, 0, .)
mata: T_OK = RT * (1 + 1e-11)
mata: T_BAD = RT * (1 + 1e-6)
mata: K_FLIP = RK
mata: K_FLIP[2] = 1
mata: K_SHORT = RK[(1::4)]
mata: K_MISS = RK
mata: K_MISS[5] = .

* indexed locals, one case each: word-of extraction splits a quoted case on its spaces
local c1  tdr_rc(RB, RT, RK, DF, RB, RT, RK, DF, TOL)
local l1  identical
local w1  0
local c2  tdr_rc(RB, RT, RK, DF, B_OK, RT, RK, DF, TOL)
local l2  one b element within tol
local w2  0
local c3  tdr_rc(RB, RT, RK, DF, RB, T_OK, RK, DF, TOL)
local l3  T3 within tol at large magnitude
local w3  0
local c4  tdr_rc(RB, 1e-20, RK, 0, RB, 3e-20, RK, 0, TOL)
local l4  saturated, T3 near 0 in both
local w4  0
local c5  tdr_rc(RB, RT, RK, DF, B_BIG, RT, RK, DF, TOL)
local l5  b off by 1e-6 at an exact zero
local w5  3498
local c6  tdr_rc(RB, RT, RK, DF, B_ONE, RT, RK, DF, TOL)
local l6  one of 20 b elements off by 1e-8
local w6  3498
local c7  tdr_rc(RB, RT, RK, DF, RB, T_BAD, RK, DF, TOL)
local l7  T3 off by 1e-6 relative
local w7  3498
local c8  tdr_rc(RB, RT, RK, DF, RB, RT, K_FLIP, DF, TOL)
local l8  bound flag flipped
local w8  3498
local c9  tdr_rc(RB, RT, RK, DF, RB, ., RK, DF, TOL)
local l9  T3 missing in the refit only
local w9  3498
local c10 tdr_rc(RB, ., RK, DF, RB, ., RK, DF, TOL)
local l10 T3 missing in both
local w10 3498
local c11 tdr_rc(RB, RT, RK, DF, RB, RT, RK, DF, .)
local l11 tol missing
local w11 3498
local c12 tdr_rc(RB, RT, RK, DF, RB, RT, RK, DF, 0)
local l12 tol = 0
local w12 3498
local c13 tdr_rc(RB, RT, RK, DF, RB, RT, RK, DF, -1e-9)
local l13 tol < 0
local w13 3498
local c14 tdr_rc(RB, RT, RK, DF, B_SHORT, RT, RK, DF, TOL)
local l14 b lengths differ
local w14 3498
local c15 tdr_rc(RB, RT, RK, DF, RB, RT, K_SHORT, DF, TOL)
local l15 bound lengths differ
local w15 3498
local c16 tdr_rc(RB, RT, RK, DF, B_MISS, RT, RK, DF, TOL)
local l16 b missing in the refit
local w16 3498
local c17 tdr_rc(B_MISS, RT, RK, DF, RB, RT, RK, DF, TOL)
local l17 b missing in the stored result
local w17 3498
local c18 tdr_rc(RB, RT, K_MISS, DF, RB, RT, RK, DF, TOL)
local l18 bound missing
local w18 3498
local c19 tdr_rc(B_EMPTY, RT, RK, DF, B_EMPTY, RT, RK, DF, TOL)
local l19 empty b
local w19 3498
local c20 tdr_rc(RB, RT, RK, DF, RB, RT, RK, DF + 1, TOL)
local l20 df differ
local w20 3498
local c21 tdr_rc(RB, RT, RK, ., RB, RT, RK, ., TOL)
local l21 df missing in both
local w21 3498
local c22 tdr_rc(RB, 0, RK, 0, RB, 1e-6, RK, 0, TOL)
local l22 saturated, T3 0 against 1e-6
local w22 3498

local nfail 0
local ncase 22

forvalues i = 1/`ncase' {
    local c `c`i''
    local l `l`i''
    local w `w`i''
    capture mata: (void) `c'
    local rc = _rc
    display as text "  c`i' `l': rc `rc'"
    if `rc' != `w' {
        capture noisily mata: (void) `c'
        display as error "test_dtrefit: c`i' `l' gave rc `rc', wanted `w'"
        local nfail = `nfail' + 1
    }
}

* t1: the production value, set 2026-10-09 from the processor-count measurement.
mata: st_local("tolok", strofreal(dwlsmv_dt_reftol() == 1e-10))
display as text "  t1 dwlsmv_dt_reftol() == 1e-10: `tolok'"
if "`tolok'" != "1" {
    display as error "test_dtrefit: dwlsmv_dt_reftol() is not the measured 1e-10"
    local nfail = `nfail' + 1
}
* t2, t3: the production value against the observed residue and beyond it.
capture mata: (void) tdr_rc(RB, RT, RK, DF, RB, RT * (1 + 8.1e-13), RK, DF, dwlsmv_dt_reftol())
local rc2 = _rc
display as text "  t2 observed maximum at the production tolerance: rc `rc2'"
if `rc2' != 0 {
    display as error "test_dtrefit: t2 refused the observed processor-count residue"
    local nfail = `nfail' + 1
}
capture mata: (void) tdr_rc(RB, RT, RK, DF, RB, RT * (1 + 1e-9), RK, DF, dwlsmv_dt_reftol())
local rc3 = _rc
display as text "  t3 relative 1e-9 at the production tolerance: rc `rc3'"
if `rc3' != 3498 {
    display as error "test_dtrefit: t3 gave rc `rc3', wanted 3498"
    local nfail = `nfail' + 1
}

if `nfail' {
    display as error _n "test_dtrefit: `nfail' failure(s)"
    exit 9
}

display as result "test_dtrefit: dwlsmv_dtrefit() accepts a refit within the (test-local) tolerance and refuses with rc 3498 any difference beyond it, a flipped bound flag, a df or T3 mismatch, an unset tolerance and malformed input"
exit 0
