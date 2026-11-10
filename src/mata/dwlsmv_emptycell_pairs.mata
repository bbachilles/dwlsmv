*! dwlsmv -- the one-empty-cell pair: membership from (Y, map) (O6, spec Module 2)
*!
*! Executed under StataMP 18 on 2026-10-06: tests/unit/test_emptycell.do,
*! test_sample_stats.do and test_dwlsmv.do pass; suite green.
*!
*!     real scalar   dwlsmv_emptycell1(real matrix n)          (own file)
*!     real matrix   dwlsmv_emptycell_repair(real matrix n)    (own file)
*!     real colvector dwlsmv_emptycell_pairs(real matrix Y,
*!                                           struct dwlsmv_statmap scalar map)
*!
*! The one place the repair rule lives (decided 2026-10-06, spec Module 2,
*! *One-empty-cell repair*). Stage 1 repairs the pairs it names; Module 10
*! excludes their rho entries from the centring row's bug-report z; Module 11
*! lists them. All three call dwlsmv_emptycell_pairs() on the same Y, so they
*! cannot disagree, and no struct field stores what is derivable from (Y, K):
*! if a repaired pair fails, Stage 1 aborts and no struct exists, so "has
*! exactly one empty cell" and "was repaired" are the same set whenever a
*! stage1 struct is in hand.
*!
*! Membership is read from the TABLE, never from a status. A one-empty-cell
*! pair is a binary 2x2 table with exactly one empty cell. NOT members:
*!   - two empty cells on a diagonal (the perfect-association pair): no repair,
*!     as in lavaan; Stage 1 aborts it;
*!   - two empty cells in one row or column: an empty category, which
*!     dwlsmv_thresholds() rejects before any pair is fitted;
*!   - any larger table with an empty cell: lavaan's rule there sets empty
*!     cells to zero, which changes nothing.
*!
*! The repair is lavaan 0.7.2's keep-margins rule with c = 0.5 (fixed, not
*! user-settable): +c on the empty cell e and its diagonal partner, -c on the
*! other two, i.e. n+ = n + c sigma_e sigma with sigma = (1, -1 \ -1, 1). Both
*! margins are unchanged, so Module 1's thresholds are unchanged
*! (derivation @sec-o6).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

real colvector dwlsmv_emptycell_pairs(real matrix Y, struct dwlsmv_statmap scalar map)
{
    real colvector r
    real scalar k, pos, vi, vj

    r = J(map.P, 1, 0)
    for (k = 1; k <= map.P; k++) {
        pos = map.T + k
        vi  = map.var_i[pos]
        vj  = map.var_j[pos]
        if (map.K[vi] == 2) {
            if (map.K[vj] == 2) {
                r[k] = dwlsmv_emptycell1(dwlsmv_crosstab(Y[., vi], Y[., vj], 2, 2))
            }
        }
    }
    return(r)
}

end
