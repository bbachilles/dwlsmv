*! dwlsmv -- DIFFTEST validation check 9: order the models, test the chain
*! (Module 12)
*!
*! Drafted 2026-10-08 by a read-only agent from ddcd1d9. Executed 2026-10-08:
*! test_dtnest.do green, mutants C1-C6 each turned it red. Generated chains
*! added 2026-10-13 (exemption gate step 2b): tests/unit/test_dtgenerated.do.
*!
*!     real colvector dwlsmv_dtchain(struct dwlsmv_parmap vector pars,
*!                                   real colvector df, string colvector names
*!                                   | string colvector cc, string colvector can,
*!                                   string colvector inv,
*!                                   pointer(real matrix) colvector Kp,
*!                                   real matrix info)
*!
*! Spec Module 12, Validation contract, check 9 (decided 2026-10-08): the
*! models are ordered by their free-parameter count, so the user may name them
*! in any order; more than two must form a chain under check 8, and adjacent
*! pairs are tested. pars holds the n >= 2 refit parameter tables, df their
*! e(df_T3), names their labels, all in the user's order. The optional cc,
*! can, inv and Kp are each model's e(canon_configural), e(canon),
*! e(invariance) and a pointer to its e(K), passed to check 8's exemption.
*!
*! Returns the permutation (n x 1) that puts the models in order of increasing
*! par.npar: ord[1] is the most restricted model (the first H0), ord[n] the
*! least restricted. Ties in npar (decided 2026-10-13, user, option (a)): when
*! every model of the tie posts a level, by level rank, most restricted first
*! (scalar, metric, thresholds, configural), so a d = 0 generated step leaves
*! the order, and which pairs are tested, independent of the order the models
*! are named in; otherwise by input position, and check 8 then refuses the
*! equal-sized pair (d >= 1). Nesting is transitive, so the adjacent pairs
*! suffice; an endpoint-only test would miss a middle model that is nested in
*! neither neighbour.
*!
*! info (optional output, (n - 1) x 2): row k is the pair (ord[k], ord[k + 1]),
*! column 1 its d, column 2 1 if it passed through the generated-chain
*! exemption. The caller reads d = 0 exempt pairs as "no test" (spec check 8,
*! zero-df pairs).
*!
*! WHAT IT DOES NOT DO. Everything dwlsmv_dtnest does not do. It does not run
*! check 11 on each H0 or check 10 on each pair: those are the caller's, per
*! pair, on the returned order (dwlsmv_dtstep). It orders by par.npar, never by
*! par.nall (fixed rows are not parameters).
*!
*! Requires dwlsmv_dtnest. Load order: dwlsmv_types.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

// Level rank, most restricted first; 0 for a model that posts no level.
real scalar dwlsmv_dtc_rank(string scalar lev)
{
    if (lev == "scalar")     return(1)
    if (lev == "metric")     return(2)
    if (lev == "thresholds") return(3)
    if (lev == "configural") return(4)
    return(0)
}

real colvector dwlsmv_dtchain(struct dwlsmv_parmap vector pars,
                              real colvector df, string colvector names,
                              | string colvector cc, string colvector can,
                              string colvector inv,
                              pointer(real matrix) colvector Kp,
                              real matrix info)
{
    struct dwlsmv_parmap scalar pa, pb
    real colvector npar, sec, ord, tie
    real matrix S
    real scalar n, k, j, d, gen, ex, a, b

    n = length(pars)
    if (n < 2) {
        errprintf("dwlsmv check 9: a comparison needs at least two models, got %g.\n", n)
        _error("dwlsmv check 9: need at least two models")
    }
    if (rows(df) != n | cols(df) != 1 | rows(names) != n | cols(names) != 1) {
        errprintf("dwlsmv check 9: %g models but %g df and %g names.\n",
                  n, rows(df) * cols(df), rows(names) * cols(names))
        _error("dwlsmv check 9: models, df and names differ in number")
    }
    gen = (args() >= 7)
    if (args() != 3 & args() < 7) {
        errprintf("dwlsmv check 9: the generated-chain inputs are incomplete.\n")
        _error("dwlsmv check 9: incomplete generated-chain inputs")
    }
    if (gen) {
        if (rows(cc) != n | rows(can) != n | rows(inv) != n | rows(Kp) != n) {
            errprintf("dwlsmv check 9: %g models but the generated-chain inputs differ in number.\n", n)
            _error("dwlsmv check 9: generated-chain inputs differ in number")
        }
    }

    npar = J(n, 1, .)
    for (k = 1; k <= n; k++) {
        npar[k] = pars[k].npar
    }

    // Secondary key: level rank when every model sharing this npar posts a
    // level, else 0 (input order).
    sec = J(n, 1, 0)
    if (gen) {
        for (k = 1; k <= n; k++) {
            tie = selectindex(npar :== npar[k])
            a = 1
            for (j = 1; j <= rows(tie); j++) {
                if (dwlsmv_dtc_rank(inv[tie[j]]) == 0) {
                    a = 0
                }
            }
            if (a) {
                sec[k] = dwlsmv_dtc_rank(inv[k])
            }
        }
    }

    S   = sort((npar, sec, (1::n)), (1, 2, 3))
    ord = S[., 3]

    if (args() == 8) {
        info = J(n - 1, 2, .)
    }
    for (k = 1; k < n; k++) {
        a  = ord[k]
        b  = ord[k + 1]
        pa = pars[a]
        pb = pars[b]
        ex = 0
        if (gen) {
            d = dwlsmv_dtnest(pa, pb, df[a], df[b], names[a], names[b],
                              cc[a], can[a], inv[a], *Kp[a],
                              cc[b], can[b], inv[b], *Kp[b], ex)
        }
        else {
            d = dwlsmv_dtnest(pa, pb, df[a], df[b], names[a], names[b])
        }
        if (args() == 8) {
            info[k, .] = (d, ex)
        }
    }

    return(ord)
}

end
