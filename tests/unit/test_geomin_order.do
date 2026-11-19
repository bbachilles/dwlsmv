*! test_geomin_order.do -- the complete column order of a rotated EFA solution (spec Module 12,
*! owner 2026-11-18, A3.7 revisited): dwlsmv_geomin_colorder() and dwlsmv_geomin_tiekeys().
*!
*! Constructed solutions, each one deciding at a known key:
*!   1  key 1 (index statistic) decides
*!   2  key 1 ties exactly, key 2 (sum of squared loadings) decides
*!   3  keys 1 and 2 tie exactly, key 3 (loadings item by item) decides
*!   4  identical loading columns at m = 3, Phi differs: key 4 decides
*!   5  identical loading columns at m = 2: nothing decides, and the posted solution is the same
*! For every case, every column arrangement of the input (Phi permuted with it) must give a
*! bitwise identical posted (Lambda, Phi): the order depends on the solution, not on the order
*! the starts happened to produce. The posted solution must be a fixed point (idempotence), and
*! tiekeys() must report the deciding key.
*! Executed under StataMP 18 on 2026-11-18 (main session).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later

version 16
clear all
set more off
set linesize 120

quietly do src/mata/dwlsmv_geomin.mata

mata:

// all column arrangements of (L, P) must post the same solution; returns failures
real scalar tgo_case(string scalar nm, real matrix L, real matrix P, real rowvector wantkeys)
{
    real scalar nf, m, first
    real rowvector ord, k
    real colvector pc
    real matrix L0, P0, Lc, Pc
    transmorphic info

    nf = 0
    m = cols(L)
    first = 1
    info = cvpermutesetup((1::m))
    while ((pc = cvpermute(info)) != J(0, 1, .)) {
        Lc = L[., pc']
        Pc = P[pc', pc']
        ord = dwlsmv_geomin_colorder(Lc, Pc)
        Lc = Lc[., ord]
        Pc = Pc[ord, ord]
        if (first) {
            L0 = Lc
            P0 = Pc
            first = 0
        }
        else if ((L0 != Lc) | (P0 != Pc)) {
            printf("{err}  %s: arrangement %s posts a different solution\n", nm, invtokens(strofreal(pc')))
            nf++
        }
    }
    if (dwlsmv_geomin_colorder(L0, P0) != (1..m)) {
        printf("{err}  %s: the posted solution is not a fixed point\n", nm)
        nf++
    }
    k = dwlsmv_geomin_tiekeys(L0, P0)
    if (k != wantkeys) {
        printf("{err}  %s: tie keys %s, expected %s\n", nm, invtokens(strofreal(k)), invtokens(strofreal(wantkeys)))
        nf++
    }
    printf("{txt}  %s: tie keys %s, %s\n", nm, invtokens(strofreal(k)), (nf ? "FAIL" : "ok"))
    return(nf)
}

real scalar tgo_all()
{
    real scalar nf
    real matrix L, P

    nf = 0
    // 1: key 1 decides (column 1 salient on rows 1-2, column 2 on rows 3-4)
    L = (0.8, 0.1 \ 0.7, 0 \ 0, 0.8 \ 0.1, 0.7)
    P = (1, 0.3 \ 0.3, 1)
    nf = nf + tgo_case("key 1", L, P, (1))
    // 2: index 2.5 for both (rows {1,4} and {2,3}); sums of squares 1.13 and 0.51
    L = (0.8, 0 \ 0, 0.5 \ 0, 0.5 \ 0.7, 0.1)
    P = (1, 0.3 \ 0.3, 1)
    nf = nf + tgo_case("key 2", L, P, (2))
    // 3: index 2.5 for both, sums of squares exactly 0.5 for both; item 1 decides
    L = (0.5, 0 \ 0, 0.5 \ 0, 0.5 \ 0.5, 0)
    P = (1, 0.3 \ 0.3, 1)
    nf = nf + tgo_case("key 3", L, P, (3))
    // 4: m = 3, columns 1 and 2 identical, phi13 != phi23 (column 3 distinct, index 4.5)
    L = (0.7, 0.7, 0 \ 0.6, 0.6, 0 \ 0, 0, 0.1 \ 0.1, 0.1, 0.8 \ 0, 0, 0.7 \ 0, 0, 0)
    P = (1, 0.5, 0.2 \ 0.5, 1, 0.4 \ 0.2, 0.4, 1)
    nf = nf + tgo_case("key 4", L, P, (4, 1))
    // 5: m = 2, identical columns: either order posts the same Lambda and Phi
    L = (0.7, 0.7 \ 0.6, 0.6 \ 0.1, 0.1)
    P = (1, 0.5 \ 0.5, 1)
    nf = nf + tgo_case("key 0", L, P, (0))
    return(nf)
}

end

local nfail = 0
mata: st_local("n1", strofreal(tgo_all()))
local nfail = `nfail' + `n1'

if `nfail' {
    display as error _n "test_geomin_order: `nfail' failure(s)"
    exit 9
}
display as result "test_geomin_order: the column order is complete and start-independent; each key decides where it should; tie keys reported"
exit 0
