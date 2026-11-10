*! measure_normal_r2.do -- the revised normal() study (A-L revision 2),
*! docs/notes/lnnormal-cert-r2/prereg.md. An extension of measure_lnnormal.do, whose
*! frozen file and hash are untouched; its helpers and normal() section are copied here.
*!
*! EXECUTED on the frozen design: StataMP 18 and 16.1 gates (docs/notes/lnnormal-cert-r2/result.md).
*! Banner changed 2026-10-30, after execution; frozen hash in prereg §6.
*!
*! Usage: do tools/measure_normal_r2.do csvdir [smoke]
*!   csvdir  directory holding normal_ref_r2.csv
*!   smoke   1 = smoke mode: checks the pipeline and prints NO error figures
*!
*! Rule: normal(): max |((yhat - hi) / hi) - lorel| / u <= 16, u = 2^-53, on
*! [beta_h, 0], beta_h = -8.3743889230674569; per Stata version. Sensitivity at 8 and 32.
version 16
clear all
set more off
capture set processors 1

args csvdir smoke
if `"`csvdir'"' == "" {
    display as error "measure_normal_r2: give the directory of the reference CSVs"
    exit 9
}
if `"`smoke'"' == "" local smoke 0

display "LN|env|stata=" c(stata_version) "|born=" c(born_date) "|os=" c(os) ///
        "|machine=" c(machine_type) "|procs=" c(processors) "|smoke=`smoke'"

mata:
void ln_require(real scalar ok, string scalar what)
{
    if (!ok) {
        errprintf("measure_normal_r2: input check failed: %s\n", what)
        exit(9)
    }
}

// beta, hi, lo from the hex columns; reports hex/decimal agreement of beta
real matrix ln_read(string scalar hicol, string scalar locol)
{
    real colvector b, bd, hi, lo

    b = strtoreal(st_sdata(., "beta_x"))
    bd = strtoreal(st_sdata(., "beta_d"))
    hi = strtoreal(st_sdata(., hicol))
    lo = strtoreal(st_sdata(., locol))
    ln_require(rows(b) > 0, "rows")
    ln_require(!missing(b) & !missing(hi) & !missing(lo), "hex parse")
    printf("LN|input|%s|rows=%g|hex_vs_decimal_beta_disagree=%g\n", hicol, rows(b),
           sum(b :!= bd))
    return((b, hi, lo))
}

// per design component: n, missing, max normalized error, and its beta
void ln_report(string scalar what, real colvector z, string colvector comp,
               real colvector b, real scalar smoke)
{
    real scalar k, i, j
    string colvector cs
    real colvector sel, zz

    cs = uniqrows(comp)
    for (k = 1; k <= rows(cs); k++) {
        sel = (comp :== cs[k])
        zz = select(z, sel)
        if (smoke) {
            printf("LN|%s|comp=%s|n=%g|missing=%g\n", what, cs[k], rows(zz), missing(zz))
        }
        else {
            printf("LN|%s|comp=%s|n=%g|missing=%g|max_norm=%9.4f\n", what, cs[k],
                   rows(zz), missing(zz), max(zz))
        }
    }
    if (smoke) {
        printf("LN|%s|all|n=%g|missing=%g|no error figures in smoke mode\n", what,
               rows(z), missing(z))
        return
    }
    // missing counts as an exceed (prereg §3); max() skips missing, so it is counted
    // Mata's | does not short-circuit: z[i] must not be read while i is missing
    i = .
    for (j = 1; j <= rows(z); j++) {
        if (z[j] < .) {
            if (i >= .) {
                i = j
            }
            else if (z[j] > z[i]) {
                i = j
            }
        }
    }
    ln_require(i < ., "no finite error")
    printf("LN|%s|all|n=%g|missing=%g|max_norm=%9.4f|at_beta=%21x (%s)\n", what, rows(z),
           missing(z), z[i], b[i], comp[i])
    printf("LN|%s|gate|c=16|exceed=%g|pass=%g|sensitivity: c=8 exceed=%g, c=32 exceed=%g\n",
           what, sum(z :> 16) , (sum(z :> 16) == 0), sum(z :> 8), sum(z :> 32))
}
end

* ---- normal() ----------------------------------------------------------------
import delimited using `"`csvdir'/normal_ref_r2.csv"', ///
    delimiters(",") stringcols(_all) varnames(1) case(preserve) clear
mata:
X = ln_read("P_hi_x", "P_lorel_x")
u = epsilon(1) / 2
ln_require(!sum(X[., 2] :<= 0), "normal reference positive")
comp = st_sdata(., "comp")
yP = normal(X[., 1])
zP = abs(((yP - X[., 2]) :/ X[., 2]) - X[., 3]) :/ u
ln_report("normal", zP, comp, X[., 1], `smoke')
end
