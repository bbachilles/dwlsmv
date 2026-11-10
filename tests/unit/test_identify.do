*! Module 3: dwls_identify() -- the canonical form, the parameter table, and
*! the Section 4.8 rank check.
*!
*! Executed under StataMP 18 on 2026-09-22 and it passes.
*!
*! WHAT THIS FILE IS FOR. Module 3 turns the CANONICAL model specification of
*! spec Section 4.8 into a free/fixed parameter table and a theta ordering. It
*! never sees the user's text -- Module 11 owns the sem-style surface syntax
*! and normalizes into this form -- which is exactly what makes this file
*! writable today: every specification below is hand-written canonical form,
*! and none of it waits on a grammar that is still open.
*!
*! The checks, in the order they run:
*!
*!   1-4   the Section 4.8 worked example (case E) parses, and the table has
*!         the shape theta's ordering requires: thresholds first in index-map
*!         order, then lambda, beta, psi, each ascending by (row, col).
*!   5-6   'metric unitvariance F' fixes psi_ff -- ptype 2, the DISTURBANCE
*!         variance -- for an ENDOGENOUS factor as well as an exogenous one.
*!         This is the 2026-09-21 decision asserted against our own code
*!         rather than against lavaan's output; test_structural_metric.do
*!         carries the lavaan measurement behind it.
*!   7     order within the file is irrelevant, asserted by shuffling.
*!   8     the marker metric is the other half, and needs an explicit psi
*!         diagonal, because the canonical form has no defaults.
*!   9     THE RANK CHECK. The overparameterized case-E model -- psi_12 free
*!         alongside beta -- is rejected. This is the check Section 4.8 is
*!         built around and the reason Modules 3, 5 and 6 are one unit: it
*!         cannot run without dwls_delta().
*!   10    the error paths abort with rc 3498 and not merely "nonzero". Both a
*!         bare _error() and a failed assert abort; only the code says whether
*!         the user was actually told anything.
*!
*!         COVERAGE, stated precisely because an earlier version of this line
*!         said "every documented error path" and that was not true.
*!         dwls_identify() has 33 _error() sites and dwlsmv_canon_status()
*!         five -- the 33rd is the cycle check of 2026-09-23 (rejects 16-18:
*!         a reciprocal pair at NONZERO starts, where the rank check alone
*!         would not fire; a three-factor loop; a loop closed by a fixed
*!         nonzero beta). The list below covers the reachable ones. The two it omits
*!         are marked in the source as backstops no canonical specification
*!         can reach: npar == 0 (every threshold becomes a free tau row, so
*!         npar >= T >= 1) and sv[1] <= 0 (Delta's leading block is an
*!         identity of order T). Do not "cover" those by weakening the code
*!         around them.
*!
*!         WHY THE COUNT MATTERS, rather than being box-ticking: a
*!         format/argument mismatch in an errprintf on an unvisited path
*!         aborts at rc 3001 AT THE PRINTF, not 3498, and nothing short of a
*!         runtime visit distinguishes the two. Not hypothetical -- see
*!         the backtick-in-Mata-string trap, which took out ten of these paths
*!         simultaneously on their first run.
*!   11    v1 IS RECURSIVE, and the cycle check does not over-reject: a reverse
*!         beta fixed at ZERO is no path, and a three-factor chain with a
*!         shortcut (F1 -> F2 -> F3, F1 -> F3) is connected but acyclic. Both
*!         are accepted. Added 2026-09-23.
*!   12    PSI ORDER (spec 4.8, settled 2026-09-28): the lower triangle,
*!         row-major in factor declaration order, whichever spelling is
*!         written; both spellings give one identical table, and writing both
*!         is rejected (reject 38), with the one-spelling control accepted. A
*!         regression pin -- the code already did this. Proved able to fire: a
*!         mutant reversing the column loop is killed by check 12 alone; the
*!         two mutants that break the lower-triangle normalization die earlier,
*!         at dwls_identify()'s own "lost a row" backstop (one-off scratch
*!         runs, 2026-09-28).
*!   13    THE SECOND POINT (spec 4.8, derivation Finding 5, 2026-10-04): two
*!         models deficient only at the symmetric starts are accepted -- a
*!         random intercept with psi_21 free, and two factors of two
*!         indicators with psi_21 free -- while check 9's model, deficient
*!         everywhere, is still rejected. A dwls_identify() without the
*!         second point aborts check 13 (a); one that shifts loadings only
*!         aborts (b).
*!
*!         Check 11 (b) uses three indicators per factor because, before the
*!         second point, two made the rank check reject it at B = 0.
*!
*! Reference fixtures are not read here. Module 3 produces a table, not
*! numbers, and the numbers it would be compared on belong to Modules 5 and 6.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off

foreach f in tests/dwlsmv_assert.mata               ///
             src/mata/dwlsmv_types.mata             ///
             src/mata/dwlsmv_mg_check.mata          ///
             src/mata/dwlsmv_mg_sub.mata            ///
             src/mata/dwlsmv_mg_parts.mata          ///
             src/mata/dwlsmv_mg_sigma.mata          ///
             src/mata/dwlsmv_mg_dcols.mata          ///
             src/mata/dwlsmv_id_parse.mata          ///
             src/mata/dwlsmv_id_rank.mata           ///
             src/mata/dwlsmv_id_tie.mata            ///
             src/mata/dwlsmv_id_groups.mata         ///
             src/mata/dwlsmv_id_opname.mata         ///
             src/mata/dwlsmv_statmap_build.mata     ///
             src/mata/dwlsmv_pos_tau.mata           ///
             src/mata/dwlsmv_pos_rho.mata           ///
             src/mata/dwlsmv_nfac.mata              ///
             src/mata/dwlsmv_parvals.mata           ///
             src/mata/dwlsmv_matrices.mata          ///
             src/mata/dwlsmv_rystar.mata            ///
             src/mata/dwlsmv_corrvec.mata           ///
             src/mata/dwls_model.mata               ///
             src/mata/dwlsmv_delta_cols.mata        ///
             src/mata/dwls_delta.mata               ///
             src/mata/dwlsmv_canon_lines.mata       ///
             src/mata/dwlsmv_canon_find.mata        ///
             src/mata/dwlsmv_canon_status.mata      ///
             src/mata/dwlsmv_numrank.mata           ///
             src/mata/dwls_identify.mata {
    capture noisily do "`f'"
    if _rc {
        display as error "test_identify: `f' failed to load"
        exit 9
    }
}

mata:

// Join directive lines with char(10). '\n' is an escape inside a printf
// FORMAT only -- written into a string literal it stays two characters and
// the whole specification arrives as one unparseable line.
string scalar tid_spec(string rowvector lines)
{
    string scalar s
    real scalar i

    s = ""
    for (i = 1; i <= cols(lines); i++) {
        s = s + lines[i] + char(10)
    }
    return(s)
}

// The Section 4.8 worked example, verbatim from the spec.
string scalar tid_caseE()
{
    return(tid_spec((
        "factor  F1 F2",
        "ind     y1 y2 y3 y4 y5 y6",
        "lambda  F1 y1 free",
        "lambda  F1 y2 free",
        "lambda  F1 y3 free",
        "lambda  F2 y4 free",
        "lambda  F2 y5 free",
        "lambda  F2 y6 free",
        "beta    F2 F1 free",
        "metric  unitvariance F1",
        "metric  unitvariance F2")))
}

// The same model with the directives in a different order, and with comments
// and blank lines thrown in.
string scalar tid_caseE_shuffled()
{
    return(tid_spec((
        "metric  unitvariance F2",
        "# the order of directives is irrelevant (spec 4.8)",
        "lambda  F2 y6 free",
        "",
        "beta    F2 F1 free",
        "lambda  F1 y3 free",
        "ind     y1 y2 y3 y4 y5 y6",
        "lambda  F2 y4 free    # trailing comment",
        "lambda  F1 y1 free",
        "metric  unitvariance F1",
        "lambda  F2 y5 free",
        "factor  F1 F2",
        "lambda  F1 y2 free")))
}

// Case E under the marker metric. The psi diagonals must be written out:
// the canonical form is a literal serialization and carries no defaults.
string scalar tid_caseE_marker()
{
    return(tid_spec((
        "factor  F1 F2",
        "ind     y1 y2 y3 y4 y5 y6",
        "lambda  F1 y2 free",
        "lambda  F1 y3 free",
        "lambda  F2 y5 free",
        "lambda  F2 y6 free",
        "beta    F2 F1 free",
        "psi     F1 F1 free",
        "psi     F2 F2 free",
        "metric  marker F1 y1",
        "metric  marker F2 y4")))
}

// The error Section 4.8 is built around: psi_12 free alongside beta.
string scalar tid_caseE_over()
{
    return(tid_caseE() + "psi     F1 F2 free" + char(10))
}

end

* ===========================================================================
* The specifications that must be REJECTED, each with the reason it exists.
* Kept at Stata level because the return code is the assertion: both a bare
* _error() and a failed assert abort, and only rc 3498 says the user was
* actually handed a message.
* ===========================================================================

mata:

void tid_add(string colvector why, string colvector spec, real colvector isE,
             string scalar w, string scalar s, real scalar e)
{
    why  = why  \ w
    spec = spec \ s
    isE  = isE  \ e
}

void tid_badlist(string colvector why, string colvector spec, real colvector isE)
{
    why  = J(0, 1, "")
    spec = J(0, 1, "")
    isE  = J(0, 1, 0)

    tid_add(why, spec, isE,
            "psi_12 free alongside beta -- the 4.8 rank check",
            tid_caseE_over(), 1)

    tid_add(why, spec, isE,
            "a factor with no metric directive",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y1 free", "lambda F1 y2 free",
                      "lambda F1 y3 free", "lambda F1 y4 free",
                      "psi F1 F1 free")), 0)

    tid_add(why, spec, isE,
            "two metric directives on one factor",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1", "metric unitvariance F1")), 0)

    tid_add(why, spec, isE,
            "a metric colliding with an explicit directive on the same element",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y1 fixed 1", "lambda F1 y2 free",
                      "lambda F1 y3 free", "lambda F1 y4 free",
                      "psi F1 F1 free", "metric marker F1 y1")), 0)

    tid_add(why, spec, isE,
            "a marker metric with no psi diagonal written out",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F1 y4 free", "metric marker F1 y1")), 0)

    tid_add(why, spec, isE,
            "a 'tau' directive -- thresholds come from the index map",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1", "tau y1 1 free")), 0)

    tid_add(why, spec, isE,
            "a 'phi' directive -- Phi_eta is derived and can never be a row",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1", "phi F1 F1 free")), 0)

    tid_add(why, spec, isE,
            "a 'theta' directive -- residual variances are derived under delta",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1", "theta y1 y1 free")), 0)

    tid_add(why, spec, isE,
            "an undeclared factor name",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F9 y2 free", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1")), 0)

    tid_add(why, spec, isE,
            "an undeclared indicator name",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y7 free", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1")), 0)

    tid_add(why, spec, isE,
            "an indicator count that does not match the index map",
            tid_spec(("factor F1", "ind y1 y2 y3",
                      "lambda F1 y2 free", "lambda F1 y3 free",
                      "psi F1 F1 free", "metric marker F1 y1")), 0)

    tid_add(why, spec, isE,
            "'fixed' with no value",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 fixed", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1")), 0)

    tid_add(why, spec, isE,
            "neither 'free' nor 'fixed'",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 maybe", "lambda F1 y3 free",
                      "lambda F1 y4 free", "psi F1 F1 free",
                      "metric marker F1 y1")), 0)

    tid_add(why, spec, isE,
            "the same element specified twice",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "lambda F1 y2 free 0.5",
                      "lambda F1 y3 free", "lambda F1 y4 free",
                      "psi F1 F1 free", "metric marker F1 y1")), 0)


    tid_add(why, spec, isE,
            "a factor regressed on itself",
            tid_spec(("factor F1 F2", "ind y1 y2 y3 y4",
                      "lambda F1 y1 free", "lambda F1 y2 free",
                      "lambda F2 y3 free", "lambda F2 y4 free",
                      "beta F1 F1 free",
                      "metric unitvariance F1", "metric unitvariance F2")), 0)

    // ---- cycles in B: v1 is recursive-only (spec 4.8, 2026-09-23) ---------
    // The reciprocal pair is given NONZERO starts on purpose: at the default
    // starts (both 0) the rank check would also reject it, for the wrong
    // reason, and this case must be refused by the cycle check alone. The
    // three-factor loop has no reciprocal pair at all, and the last closes
    // its loop with a FIXED nonzero beta, which is still a path.
    tid_add(why, spec, isE,
            "a reciprocal pair of beta paths, at nonzero starts",
            tid_spec(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
                      "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                      "beta F2 F1 free 0.3", "beta F1 F2 free 0.4",
                      "metric unitvariance F1", "metric unitvariance F2")), 1)

    tid_add(why, spec, isE,
            "a three-factor feedback loop F1 -> F2 -> F3 -> F1",
            tid_spec(("factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6",
                      "lambda F1 y1 free", "lambda F1 y2 free",
                      "lambda F2 y3 free", "lambda F2 y4 free",
                      "lambda F3 y5 free", "lambda F3 y6 free",
                      "beta F2 F1 free 0.3", "beta F3 F2 free 0.3", "beta F1 F3 free 0.3",
                      "metric unitvariance F1", "metric unitvariance F2",
                      "metric unitvariance F3")), 1)

    tid_add(why, spec, isE,
            "a loop closed by a FIXED nonzero beta",
            tid_spec(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
                      "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                      "beta F2 F1 free", "beta F1 F2 fixed 0.2",
                      "metric unitvariance F1", "metric unitvariance F2")), 1)

    // ---- the declaration pass, which nothing above reached ----------------
    // Added 2026-09-22 after a review observed that the cases above exercise
    // 15 of the 37 _error() sites and all of them sit downstream of the
    // declarations. The gap matters for one specific reason: a format/argument
    // mismatch on an untested path aborts at rc 3001 at the printf, NOT 3498,
    // and only a runtime visit distinguishes the two. A project coding rule.

    tid_add(why, spec, isE, "an empty specification", "", 0)

    tid_add(why, spec, isE, "a specification that is only comments",
            tid_spec(("# nothing but a comment", "   ", "# and another")), 0)

    tid_add(why, spec, isE, "no 'factor' directive",
            tid_spec(("ind y1 y2 y3 y4", "lambda F1 y2 free")), 0)

    tid_add(why, spec, isE, "no 'ind' directive",
            tid_spec(("factor F1", "lambda F1 y2 free")), 0)

    tid_add(why, spec, isE, "'factor' declaring no names",
            tid_spec(("factor", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "metric unitvariance F1")), 0)

    tid_add(why, spec, isE, "the same factor declared twice on one line",
            tid_spec(("factor F1 F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "metric unitvariance F1")), 0)

    tid_add(why, spec, isE, "the same indicator declared twice",
            tid_spec(("factor F1", "ind y1 y2 y3 y3",
                      "lambda F1 y2 free", "metric unitvariance F1")), 0)

    // Split declarations are refused so that 4.8's "order within the file is
    // irrelevant" stays literally true: declaration ORDER fixes Psi's and B's
    // indices, so two factor lines would make the order of those lines matter.
    tid_add(why, spec, isE, "a second 'factor' directive",
            tid_spec(("factor F1", "factor F2", "ind y1 y2 y3 y4",
                      "lambda F1 y1 free", "lambda F2 y3 free",
                      "metric unitvariance F1", "metric unitvariance F2")), 0)

    tid_add(why, spec, isE, "a second 'ind' directive",
            tid_spec(("factor F1", "ind y1 y2", "ind y3 y4",
                      "lambda F1 y2 free", "metric unitvariance F1")), 0)

    tid_add(why, spec, isE, "'metric' with nothing after it",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "psi F1 F1 free", "metric")), 0)

    tid_add(why, spec, isE, "'metric marker' missing its indicator",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "psi F1 F1 free",
                      "metric marker F1")), 0)

    tid_add(why, spec, isE, "'metric unitvariance' given an extra token",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "metric unitvariance F1 y1")), 0)

    tid_add(why, spec, isE, "an unrecognized metric method",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "psi F1 F1 free",
                      "metric effectscoding F1")), 0)

    tid_add(why, spec, isE, "'metric marker' naming an undeclared factor",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "psi F1 F1 free",
                      "metric marker F9 y1")), 0)

    tid_add(why, spec, isE, "'metric marker' naming an undeclared indicator",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free", "psi F1 F1 free",
                      "metric marker F1 y9")), 0)

    tid_add(why, spec, isE, "a directive with too many tokens",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free 0.7 extra",
                      "psi F1 F1 free", "metric marker F1 y1")), 0)

    tid_add(why, spec, isE, "a non-numeric start value",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2 free abc",
                      "psi F1 F1 free", "metric marker F1 y1")), 0)

    tid_add(why, spec, isE, "a lambda directive with only three tokens",
            tid_spec(("factor F1", "ind y1 y2 y3 y4",
                      "lambda F1 y2", "psi F1 F1 free",
                      "metric marker F1 y1")), 0)

    // More free parameters than sample statistics, i.e. negative df. Caught
    // before the rank check, which would otherwise report a deficiency and
    // blame the wrong thing: rank(Delta) <= D < npar necessarily, so every
    // such model looks overparameterized in the 4.8 sense whether or not it
    // is. Eight free cross-loadings on four indicators does it -- npar = 18
    // against D = 16.
    //
    // Note there is no companion case for npar == 0: Module 3 generates one
    // FREE tau row per threshold and the index map guarantees T >= p >= 1, so
    // npar >= 1 always. That branch of dwls_identify() is a backstop and says
    // so; it is not reachable from any canonical specification.
    tid_add(why, spec, isE, "more free parameters than sample statistics",
            tid_spec(("factor F1 F2", "ind y1 y2 y3 y4",
                      "lambda F1 y1 free", "lambda F1 y2 free",
                      "lambda F1 y3 free", "lambda F1 y4 free",
                      "lambda F2 y1 free", "lambda F2 y2 free",
                      "lambda F2 y3 free", "lambda F2 y4 free",
                      "metric unitvariance F1", "metric unitvariance F2")), 0)

    // Spec 4.8 (2026-09-28): both spellings of an off-diagonal psi name ONE
    // element. No beta here, so the rank check has nothing to catch; check 12's
    // control accepts the same model with one spelling alone. Appended last so
    // the numbers of the rejects above, which docs cite, do not move.
    tid_add(why, spec, isE,
            "an off-diagonal psi written in both spellings",
            tid_spec(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
                      "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                      "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                      "psi F2 F1 free", "psi F1 F2 free",
                      "metric unitvariance F1", "metric unitvariance F2")), 1)
}

// Each bad specification is tried against the map whose p it claims -- except
// the one whose whole point is that the two disagree. The map is rebuilt here
// rather than passed in so that no struct has to survive at the interactive
// level, where a bare declaration is transmorphic and member access fails.
void tid_try(string scalar spec, real scalar isE)
{
    struct dwlsmv_parmap  scalar p
    struct dwlsmv_statmap scalar map

    if (isE) {
        map = dwlsmv_statmap_build((3 \ 3 \ 3 \ 3 \ 3 \ 3))
    }
    else {
        map = dwlsmv_statmap_build((4 \ 4 \ 3 \ 3))
    }

    p = dwls_identify(spec, map)
}

end

mata: BADWHY = J(0, 1, "")
mata: BADSPEC = J(0, 1, "")
mata: BADE = J(0, 1, 0)
mata: tid_badlist(BADWHY, BADSPEC, BADE)
mata: st_local("nbad", strofreal(rows(BADSPEC)))

local nbadfail 0
forvalues i = 1/`nbad' {
    mata: st_local("why", BADWHY[`i'])
    capture noisily mata: tid_try(BADSPEC[`i'], BADE[`i'])
    local rc = _rc
    if `rc' == 3498 {
        display as text "    reject `i' (rc 3498): `why'"
    }
    else {
        display as error "    reject `i' FAILED with rc `rc': `why'"
        local nbadfail = `nbadfail' + 1
    }
}

* ===========================================================================
* The specifications that must be ACCEPTED, and what their tables must say.
* ===========================================================================

mata:

real scalar tid_checks()
{
    struct dwlsmv_parmap  scalar p, q, m
    struct dwlsmv_statmap scalar mapE, mapD
    string colvector FAIL
    real scalar i, nT, ok, run, seenl, seenb, seenp
    string colvector pnames

    mapE = dwlsmv_statmap_build((3 \ 3 \ 3 \ 3 \ 3 \ 3))
    mapD = dwlsmv_statmap_build((4 \ 4 \ 3 \ 3))

    FAIL = J(0, 1, "")
    nT   = mapE.T

    p = dwls_identify(tid_caseE(), mapE)

    // ---- 1. the table's size -------------------------------------------
    // 12 thresholds + 6 loadings + 1 beta + 2 metric-fixed psi = 21 rows,
    // of which the two fixed psi are not free.
    printf("{txt}  1. case E table              : nall = %f, npar = %f\n", p.nall, p.npar)
    if (p.nall != nT + 9) {
        FAIL = FAIL \ ("case E nall is " + strofreal(p.nall) + ", expected " + strofreal(nT + 9))
    }
    if (p.npar != nT + 7) {
        FAIL = FAIL \ ("case E npar is " + strofreal(p.npar) + ", expected " + strofreal(nT + 7))
    }

    // ---- 2. thresholds come FIRST, in index-map order --------------------
    // This is what makes Delta's leading block a literal identity. A
    // permutation appearing here is the tell that something reordered theta.
    ok = 1
    for (i = 1; i <= nT; i++) {
        if (p.ptype[i] != 4)            ok = 0
        if (p.row[i]   != mapE.var_i[i]) ok = 0
        if (p.col[i]   != mapE.thr_k[i]) ok = 0
        if (p.name[i]  != mapE.label[i]) ok = 0
        if (p.free[i]  != 1)            ok = 0
        if (p.fidx[i]  != i)            ok = 0
    }
    printf("{txt}  2. thresholds first          : %f free tau rows in map order, fidx 1..%f\n",
           nT, nT)
    if (!ok) {
        FAIL = FAIL \ "the leading T rows are not the thresholds in index-map order"
    }

    // ---- 3. then lambda, then beta, then psi -----------------------------
    // Section 4.8 fixes the GROUP order as lambda, beta, psi. ptype's codes
    // are 1 lambda, 2 psi, 3 beta, so the group order is NOT ascending ptype
    // and a check written on ptype alone would pass on the wrong ordering.
    seenl = 0
    seenb = 0
    seenp = 0
    ok    = 1
    for (i = nT + 1; i <= p.nall; i++) {
        if (p.ptype[i] == 1) {
            seenl = i
            if (seenb | seenp) ok = 0
        }
        else if (p.ptype[i] == 3) {
            seenb = i
            if (seenp) ok = 0
        }
        else if (p.ptype[i] == 2) {
            seenp = i
        }
        else {
            ok = 0
        }
    }
    printf("{txt}  3. group order               : last lambda %f, last beta %f, last psi %f\n",
           seenl, seenb, seenp)
    if (!ok) {
        FAIL = FAIL \ "the structural rows are not grouped lambda, beta, psi in that order"
    }

    // ---- 4. within a group, ascending by (row, col) ----------------------
    ok  = 1
    run = 0
    for (i = nT + 1; i <= p.nall; i++) {
        if (p.ptype[i] == 1) {
            if (p.row[i] * 1000 + p.col[i] <= run) ok = 0
            run = p.row[i] * 1000 + p.col[i]
        }
    }
    printf("{txt}  4. lambda ascending (row,col): %s\n", (ok ? "yes" : "NO"))
    if (!ok) {
        FAIL = FAIL \ "lambda rows are not ascending by (row, col)"
    }

    // ---- 5. 'metric unitvariance' fixes PSI, for both factors ------------
    // ptype 2 is psi. Phi_eta has no rows at all and cannot acquire one.
    for (i = 1; i <= p.nall; i++) {
        if (p.ptype[i] == 2) {
            printf("{txt}  5. metric row                : %-14s ptype %f row %f col %f free %f value %g\n",
                   p.name[i], p.ptype[i], p.row[i], p.col[i], p.free[i], p.value[i])
            if (p.row[i] != p.col[i]) {
                FAIL = FAIL \ "unitvariance produced an off-diagonal psi row"
            }
            if (p.free[i] != 0 | p.value[i] != 1) {
                FAIL = FAIL \ ("unitvariance did not fix " + p.name[i] + " at 1")
            }
        }
    }

    // ---- 6. the endogenous factor is NOT treated differently -------------
    // F2 is predicted by F1, and what unitvariance fixes for it is still
    // psi_22, the DISTURBANCE variance.
    //
    // WHERE THE DISCRIMINATION ACTUALLY LIVES, because it is not here.
    // The project's case-E trap warns about checks that "pass identically under
    // both candidate answers", and this is one: dwls_identify() can only ever
    // write ptype 1-4, and dwlsmv_parmap has no phi code at all, so a wrong
    // implementation has nothing to emit that this would reject. What it does
    // establish is structural -- exactly two psi rows, both diagonal, both
    // fixed at 1 -- which is worth asserting but does not separate the
    // disturbance reading from the total-variance one.
    //
    // The separating checks are in test_model.do: check 6 compares DERIVED
    // Theta against lavaan's fitted residual variances at lavaan's own
    // loadings, which diverge under the total-variance reading, and check 7
    // rejects an implied Var(F2) of 1. Both are numerical and both would fail.
    ok = 0
    k  = 0
    for (i = 1; i <= p.nall; i++) {
        if (p.ptype[i] == 2) {
            k = k + 1
            if (p.row[i] == 2 & p.col[i] == 2) ok = 1
        }
    }
    printf("{txt}  6. endogenous F2             : %f psi row(s), psi[F2,F2] present as a fixed row: %s\n",
           k, (ok ? "yes" : "NO"))
    if (!ok) {
        FAIL = FAIL \ "no psi[F2,F2] row for the endogenous factor"
    }
    if (k != 2) {
        FAIL = FAIL \ ("the unitvariance table carries " + strofreal(k) +
                       " psi row(s), expected exactly 2")
    }

    // ---- 7. order within the file is irrelevant --------------------------
    q  = dwls_identify(tid_caseE_shuffled(), mapE)
    ok = (q.nall == p.nall & q.npar == p.npar)
    if (ok) {
        for (i = 1; i <= p.nall; i++) {
            if (q.name[i]  != p.name[i])  ok = 0
            if (q.ptype[i] != p.ptype[i]) ok = 0
            if (q.row[i]   != p.row[i])   ok = 0
            if (q.col[i]   != p.col[i])   ok = 0
            if (q.free[i]  != p.free[i])  ok = 0
            if (q.value[i] != p.value[i]) ok = 0
            if (q.fidx[i]  != p.fidx[i])  ok = 0
        }
    }
    printf("{txt}  7. order-independence        : shuffled spec gives an identical table: %s\n",
           (ok ? "yes" : "NO"))
    if (!ok) {
        FAIL = FAIL \ "reordering the directives changed the parameter table"
    }

    // ---- 8. the marker metric --------------------------------------------
    // Same model, other metric: two loadings fixed at 1 and two psi diagonals
    // free, so npar is unchanged. That it is unchanged is the point -- the
    // two metrics spend the same budget.
    m  = dwls_identify(tid_caseE_marker(), mapE)
    ok = 0
    for (i = 1; i <= m.nall; i++) {
        if (m.ptype[i] == 1 & m.free[i] == 0 & m.value[i] == 1) ok = ok + 1
    }
    printf("{txt}  8. marker metric             : nall = %f, npar = %f, %f loading(s) fixed at 1\n",
           m.nall, m.npar, ok)
    if (m.npar != p.npar) {
        FAIL = FAIL \ ("marker npar " + strofreal(m.npar) + " differs from unitvariance npar " +
                       strofreal(p.npar))
    }
    if (ok != 2) {
        FAIL = FAIL \ "the marker metric did not fix exactly two loadings at 1"
    }

    // ---- 9. a single-factor case D model, as the simplest accepted form --
    q = dwls_identify(tid_spec(("factor f",
                                "ind y1 y2 y3 y4",
                                "lambda f y1 free 0.75",
                                "lambda f y2 free 0.65",
                                "lambda f y3 free 0.55",
                                "lambda f y4 free 0.45",
                                "metric unitvariance f")), mapD)
    printf("{txt}  9. case D one factor         : nall = %f, npar = %f (T = %f)\n",
           q.nall, q.npar, mapD.T)
    if (q.nall != mapD.T + 5 | q.npar != mapD.T + 4) {
        FAIL = FAIL \ "case D single-factor table is the wrong size"
    }

    // ---- 11. the cycle check does not over-reject --------------------------
    // Two recursive models it must ACCEPT. (a) A reverse beta FIXED AT ZERO is
    // no path, so F2 <- F1 with 'beta F1 F2 fixed 0' is recursive; the table
    // still carries the fixed row. (b) F1 -> F2 -> F3 with the shortcut
    // F1 -> F3: every factor reaches F3 two ways, but none reaches itself.
    // A check that treated "connected" as "cyclic" would refuse (b). (b) has
    // three indicators per factor, on its own nine-indicator map: with two,
    // B = 0 at the starts leaves each factor an isolated two-indicator model
    // whose loadings are identified only as a product, and the RANK check
    // rejects it -- which was this check's first version.
    q = dwls_identify(tid_spec(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
                                "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                                "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                                "beta F2 F1 free", "beta F1 F2 fixed 0",
                                "metric unitvariance F1", "metric unitvariance F2")), mapE)
    m = dwls_identify(tid_spec(("factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
                                "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                                "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                                "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
                                "beta F2 F1 free", "beta F3 F2 free", "beta F3 F1 free",
                                "metric unitvariance F1", "metric unitvariance F2",
                                "metric unitvariance F3")),
                      dwlsmv_statmap_build(J(9, 1, 3)))
    printf("{txt} 11. recursive models accepted : fixed-zero reverse path, %f beta row(s); three-factor chain with shortcut, %f free beta(s)\n",
           sum(q.ptype :== 3), sum(m.ptype :== 3 :& m.free))
    if (sum(q.ptype :== 3) != 2 | sum(q.ptype :== 3 :& q.free) != 1) {
        FAIL = FAIL \ "a fixed-zero reverse beta was not carried as one free and one fixed beta row"
    }
    if (sum(m.ptype :== 3 :& m.free) != 3) {
        FAIL = FAIL \ "the recursive three-factor chain did not keep its three free betas"
    }

    // ---- 12. psi order: lower triangle, row-major (spec 4.8, 2026-09-28) ----
    // Three correlated factors, covariances written in mixed spellings and out
    // of order. The psi rows must come out psi[F1,F1], [F2,F1], [F2,F2],
    // [F3,F1], [F3,F2], [F3,F3] -- the lower triangle, row by row in factor
    // declaration order -- and the all-lower spelling must give the identical
    // table. A regression pin: dwls_identify() already did this when the rule
    // was written down; no test held it.
    q = dwls_identify(tid_spec(("factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
                                "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                                "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                                "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
                                "psi F1 F3 free", "psi F3 F2 free", "psi F1 F2 free",
                                "metric unitvariance F1", "metric unitvariance F2",
                                "metric unitvariance F3")),
                      dwlsmv_statmap_build(J(9, 1, 3)))
    m = dwls_identify(tid_spec(("factor F1 F2 F3", "ind y1 y2 y3 y4 y5 y6 y7 y8 y9",
                                "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                                "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                                "lambda F3 y7 free", "lambda F3 y8 free", "lambda F3 y9 free",
                                "psi F2 F1 free", "psi F3 F1 free", "psi F3 F2 free",
                                "metric unitvariance F1", "metric unitvariance F2",
                                "metric unitvariance F3")),
                      dwlsmv_statmap_build(J(9, 1, 3)))
    pnames = select(q.name, q.ptype :== 2)
    printf("{txt} 12. psi order, mixed spellings: %s\n", invtokens(pnames'))
    if (pnames != ("psi[F1,F1]" \ "psi[F2,F1]" \ "psi[F2,F2]" \
                   "psi[F3,F1]" \ "psi[F3,F2]" \ "psi[F3,F3]")) {
        FAIL = FAIL \ "psi rows are not the lower triangle, row-major in factor order"
    }
    if (sum((q.ptype :== 2) :& (q.row :< q.col))) {
        FAIL = FAIL \ "a psi row is stored in the upper triangle"
    }
    if (q.name != m.name | q.row != m.row | q.col != m.col | q.free != m.free |
        q.nall != m.nall | q.npar != m.npar) {
        FAIL = FAIL \ "the upper and lower spellings of psi gave different tables"
    }
    // The control for the both-spellings rejection above: one spelling alone,
    // on the same two-factor model, is accepted.
    q = dwls_identify(tid_spec(("factor F1 F2", "ind y1 y2 y3 y4 y5 y6",
                                "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                                "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                                "psi F1 F2 free",
                                "metric unitvariance F1", "metric unitvariance F2")), mapE)
    printf("{txt} 12. one spelling alone accepted: %f free psi covariance row(s), named %s\n",
           sum((q.ptype :== 2) :& q.free), invtokens(select(q.name, (q.ptype :== 2) :& q.free)'))
    if (select(q.name, (q.ptype :== 2) :& q.free) != "psi[F2,F1]") {
        FAIL = FAIL \ "the upper spelling alone did not give the one row psi[F2,F1]"
    }

    // ---- 13. the second point (spec 4.8, derivation Finding 5, 2026-10-04) ----
    // Two models that lose rank at the symmetric starts alone must be ACCEPTED;
    // check 9's model, deficient everywhere, is still rejected there.
    // (a) A random intercept G (all loadings fixed 1, psi_GG free, orthogonal)
    // with psi_21 free: rank 19 of 20 at equal loadings. (b) Two factors of
    // two indicators with psi_21 free: rank 11 of 13 at psi_21 = 0, whatever
    // the loadings, so a second point that shifted loadings only would still
    // refuse it.
    q = dwls_identify(tid_spec(("factor F1 F2 G", "ind y1 y2 y3 y4 y5 y6",
                                "lambda F1 y1 free", "lambda F1 y2 free", "lambda F1 y3 free",
                                "lambda F2 y4 free", "lambda F2 y5 free", "lambda F2 y6 free",
                                "metric marker G y1",
                                "lambda G y2 fixed 1", "lambda G y3 fixed 1", "lambda G y4 fixed 1",
                                "lambda G y5 fixed 1", "lambda G y6 fixed 1",
                                "metric unitvariance F1", "metric unitvariance F2",
                                "psi F2 F1 free", "psi G G free")), mapE)
    m = dwls_identify(tid_spec(("factor F1 F2", "ind y1 y2 y3 y4",
                                "lambda F1 y1 free", "lambda F1 y2 free",
                                "lambda F2 y3 free", "lambda F2 y4 free",
                                "metric unitvariance F1", "metric unitvariance F2",
                                "psi F2 F1 free")),
                      dwlsmv_statmap_build(J(4, 1, 3)))
    printf("{txt} 13. deficient at the starts only, accepted: random intercept %f free; two-by-two %f free\n",
           q.npar, m.npar)
    if (q.npar != mapE.T + 8) {
        FAIL = FAIL \ "the random-intercept model with psi_21 free does not have T + 8 free parameters"
    }
    if (m.npar != 8 + 5) {
        FAIL = FAIL \ "the two-by-two model with psi_21 free does not have 13 free parameters"
    }

    printf("{txt}test_identify: %f failure(s)\n", rows(FAIL))
    for (i = 1; i <= rows(FAIL); i++) {
        printf("{err}    %s\n", FAIL[i])
    }

    return(rows(FAIL))
}

end

mata: st_local("nfail", strofreal(tid_checks()))

local ntot = `nfail' + `nbadfail'
if `ntot' {
    display as error _n "test_identify: `nfail' table failure(s), `nbadfail' error-path failure(s)"
    exit 9
}

display as result "test_identify: the canonical form parses, theta orders thresholds first then lambda/beta/psi, unitvariance fixes psi for an endogenous factor, directive order is irrelevant, and every reachable error path aborts with rc 3498 including the 4.8 rank check"
exit 0
