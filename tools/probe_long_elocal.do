*! dwlsmv -- probe: do long e() locals survive estimates store/restore and
*! save/use intact? (spec 4.8 G8, "Not established: that Stata 16's estimates
*! store/restore keep long e() locals intact (a milestone probe)")
*!
*! A probe, not a test (tools/). Check 8's generated-chain exemption compares
*! each fit's e(canon) byte for byte with dwlsmv_levelgen's regeneration from
*! e(canon_configural) and e(K) (condition (3)); the comparison command will
*! read them from stored estimates. If a Stata version truncated or altered
*! either text on restore, every generated pair would fail the exemption there
*! (fail-closed: refused by check 8's tie rule, with the failed-claim note).
*!
*! Records, for e(canon) and e(canon_configural) of an invariance(scalar) fit
*! on fixtures/data/mg_E.csv: the length in bytes, the number of newlines, an
*! exact checksum sum_i i * ascii(s_i) (integers, exact in double at these
*! lengths), and whether e(canon) equals the regeneration; at four points:
*!   fit       right after the fit
*!   restore   after estimates store, a second (configural) fit, estimates restore
*!   use       after estimates save to a file, estimates clear, estimates use
*!   same      1 iff each text is identical to its copy taken right after the fit
*! Output lines prefixed LONGE are the record. Run from the repository root:
*!     StataMP-64.exe /e do tools/probe_long_elocal.do      (log: probe_long_elocal.log)
*!     grep -E "^(====|LONGE)" probe_long_elocal.log | diff - tools/probe_long_elocal_stata18.txt
*! On Stata 16 every LONGE line except LONGEENV must match the Stata 18 record.
*!
*! Written and run on StataMP 18 by the main session, 2026-10-13: every text
*! intact at every point (record tools/probe_long_elocal_stata18.txt). The
*! synthetic case was added after the first run showed design E's text is
*! only 815 bytes; its first run stopped at estimates store (r(301), the
*! minimal program set no e(cmd)), fixed in the probe. Run on Stata 16.1 by
*! the user on the milestone machine, 2026-10-13, at f4818f5: identical to the
*! Stata 18 record except LONGEENV (tools/probe_long_elocal_stata16.txt).
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16
clear all
set more off
set linesize 200

quietly do tools/build_mlib.do
local here `c(pwd)'
adopath ++ "`here'/build"
adopath ++ "`here'/src/ado"
mata: mata mlib index

display "LONGEENV stata_version " c(stata_version) " processors " c(processors)

mata:
void le_record(string scalar point)
{
    string scalar c, cc
    real scalar k

    c  = st_global("e(canon)")
    cc = st_global("e(canon_configural)")
    printf("LONGE %s canon length %g newlines %g checksum %21.0g\n", point,
           strlen(c), le_nl(c), le_sum(c))
    printf("LONGE %s canon_configural length %g newlines %g checksum %21.0g\n", point,
           strlen(cc), le_nl(cc), le_sum(cc))
    printf("LONGE %s regenerates %g\n", point,
           dwlsmv_levelgen(cc, st_matrix("e(K)"), st_global("e(invariance)")) == c)
    if (point == "fit") {
        *crexternal("LE_c")  = c
        *crexternal("LE_cc") = cc
    }
    else {
        printf("LONGE %s same canon %g canon_configural %g\n", point,
               c == *findexternal("LE_c"), cc == *findexternal("LE_cc"))
    }
}

real scalar le_nl(string scalar s)
{
    return(strlen(s) - strlen(subinstr(s, char(10), "")))
}

real scalar le_sum(string scalar s)
{
    real rowvector a

    if (strlen(s) == 0) {
        return(0)
    }
    a = ascii(s)
    return(sum((1..cols(a)) :* a))
}
end

import delimited using "fixtures/data/mg_E.csv", clear varnames(1) delimiters(",")
local mdl (F -> y1 y2 y3 y4 y5), group(g) metric(unitvariance)

display "==== fit"
quietly dwlsmv `mdl' invariance(scalar)
display "LONGE fit invariance `e(invariance)'"
mata: le_record("fit")

display "==== restore"
estimates store A
quietly dwlsmv `mdl' invariance(configural)
display "LONGE second fit invariance `e(invariance)'"
estimates restore A
display "LONGE restore invariance `e(invariance)'"
mata: le_record("restore")

display "==== use"
tempfile est
quietly estimates save "`est'"
estimates clear
ereturn clear
quietly estimates use "`est'"
display "LONGE use invariance `e(invariance)'"
mata: le_record("use")

* A synthetic long text (design E's e(canon) is only 815 bytes; a 30-item,
* 3-group model runs to tens of kilobytes): 200,000 bytes, 10,000 lines of
* 19 characters plus a newline, posted as e(canon) by a minimal eclass program.
display "==== synthetic"
mata:
void le_synth()
{
    string scalar s, line
    real scalar k

    s = ""
    for (k = 1; k <= 10000; k++) {
        line = "lambda F" + strofreal(mod(k, 97), "%02.0f") + " y" +
               strofreal(k, "%05.0f") + " free    "
        s = s + substr(line, 1, 19) + char(10)
    }
    *crexternal("LE_s") = s
    st_global("e(canon)", s)
}

void le_synth_record(string scalar point)
{
    string scalar c

    c = st_global("e(canon)")
    printf("LONGE synthetic %s length %g newlines %g checksum %21.0g same %g\n", point,
           strlen(c), le_nl(c), le_sum(c), c == *findexternal("LE_s"))
}
end
capture program drop le_post
program le_post, eclass
    tempname b
    matrix `b' = (1)
    matrix colnames `b' = x
    ereturn post `b'
    ereturn local cmd "le_post"
    mata: le_synth()
end
le_post
mata: le_synth_record("posted")
estimates store S
quietly dwlsmv `mdl' invariance(configural)
estimates restore S
mata: le_synth_record("restore")
tempfile est2
quietly estimates save "`est2'"
estimates clear
ereturn clear
quietly estimates use "`est2'"
mata: le_synth_record("use")
