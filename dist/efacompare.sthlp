{smcl}
{* every example run 2026-11-21 on the 0.6.0 release tree (cdaaa78 plus the 0.6.0 help text), StataMP 18, linesize 80, on build/dwlsmv_example.dta: rc 0 throughout}
{* *! version 0.6.0  21nov2026}{...}
{vieweralsosee "dwlsmv" "help dwlsmv"}{...}
{vieweralsosee "[MV] screeplot" "help screeplot"}{...}
{viewerjumpto "Syntax" "efacompare##syntax"}{...}
{viewerjumpto "Description" "efacompare##description"}{...}
{viewerjumpto "Options" "efacompare##options"}{...}
{viewerjumpto "Remarks" "efacompare##remarks"}{...}
{viewerjumpto "Examples" "efacompare##examples"}{...}
{viewerjumpto "Stored results" "efacompare##results"}{...}
{viewerjumpto "Also see" "efacompare##alsosee"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:efacompare} {hline 2}}Fit and over-extraction diagnostics by number of factors{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:efacompare} {varlist} {ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt max:factors(#)}}largest number of factors to fit; the default is the largest
with at least one degree of freedom{p_end}
{synopt:{opt eps:ilon(#)}}geomin epsilon, the same for every number of factors{p_end}
{synopt:{opt rstarts(#)}}number of geomin starts{p_end}
{synopt:{opt keepbound}}as for {helpb dwlsmv}{p_end}
{synoptline}
{p2colreset}{...}

{p 4 4 2}
{it:varlist} is the list of the indicators: at least four numeric variables, which must not be
named {cmd:F1}, {cmd:F2}, and so on (those are the names {cmd:efacompare} gives its factors).


{marker description}{...}
{title:Description}

{p 4 4 2}
{cmd:efacompare} fits the exploratory model of {cmd:dwlsmv}'s {opt efa()} (oblique geomin) with
one factor, two factors, and so on, and shows one row per number of factors, m. Each row is
the result of the standalone command

{phang2}{cmd:. dwlsmv (F1 F2 -> }{it:varlist}{cmd:), efa(F1 F2)}{p_end}

{p 4 4 2}
for that m, on the same sample, and reads its stored results. The command reports the fit,
and diagnostics from the rotated solution that bear on whether m is too large. It gives no
verdict and applies no cutoff: it does not name a number of factors, and it marks no row as best.

{p 4 4 2}
{bf:Columns, first block.}

{p2colset 8 20 22 2}{...}
{p2col:{cmd:m}}the number of factors{p_end}
{p2col:{cmd:df}}degrees of freedom, {it:p}({it:p} - 1)/2 - ({it:p}m - m(m - 1)/2) for {it:p} items;
the thresholds are saturated{p_end}
{p2col:{cmd:status}}{cmd:fitted}: the fit converged. {cmd:not conv}: the unrotated fit did not
converge. {cmd:rot fail}: the geomin rotation did not converge from any start{p_end}
{p2col:{cmd:T3}, {cmd:p}}the scaled-and-shifted test statistic and its p-value, {cmd:e(T3)} and
{cmd:e(p_T3)}{p_end}
{p2col:{cmd:RMSEA}, {cmd:CFI}, {cmd:TLI}}the scaled fit indices, {cmd:e(rmsea_s)}, {cmd:e(cfi_s)},
{cmd:e(tli_s)}{p_end}
{p2colreset}{...}

{p 4 4 2}
{bf:Columns, second block.}

{p2colset 8 20 22 2}{...}
{p2col:{cmd:m}}the number of factors{p_end}
{p2col:{cmd:minima}}the number of distinct minima of the geomin criterion that the starts
reached, {cmd:e(rot_minima)}{p_end}
{p2col:{cmd:gap}}the difference in the criterion between the best and the second-best distinct
minimum, {cmd:e(rot_gap)}; missing when only one minimum was found, and at m = 1{p_end}
{p2col:{cmd:lam3}}for each factor, the third-largest absolute value among its loadings; the
smallest of these over the factors. A factor defined by fewer than three loadings gives
a small value{p_end}
{p2col:{cmd:phimax}}the largest absolute off-diagonal element of the factor correlation matrix
{cmd:e(Phi)}; missing at m = 1{p_end}
{p2col:{cmd:comm1}}1 if the unrotated fit is inadmissible because an indicator's communality
exceeds 1 ({cmd:e(rot_nose)} is {cmd:inadmissible}), 0 otherwise{p_end}
{p2col:{cmd:ev_m}}the m-th eigenvalue of the sample polychoric (tetrachoric) correlation matrix,
from {cmd:e(Ev)}; the same eigenvalues {helpb screeplot} plots{p_end}
{p2colreset}{...}

{p 4 4 2}
A value that does not exist for a row is shown as {cmd:.}.

{p 4 4 2}
{cmd:efacompare} gives no difference test between m and m + 1 (a design decision; {helpb difftest}
refuses EFA results). Parallel analysis is not provided. The eigenvalue column is descriptive,
as is {helpb screeplot}.


{marker options}{...}
{title:Options}

{phang}
{opt maxfactors(#)} sets the largest number of factors fitted; the table has rows for
m = 1, ..., {it:#}. By default the table ends at the largest m that leaves at least one degree of
freedom. {opt maxfactors()} may reach the one m with df = 0, which is shown with a note. A
larger value is refused (rc 198), and the message names the largest m allowed.
Fewer than four items are refused: with three, one factor already has no degrees of freedom.

{phang}
{opt epsilon(#)} and {opt rstarts(#)} are passed to every fit and mean what they mean in
{helpb dwlsmv}. Without them each row uses what the standalone fit would use, which includes
a geomin epsilon that depends on m; see {help efacompare##remarks:Remarks}.

{phang}
{opt keepbound} is passed to every fit; see {helpb dwlsmv}.

{p 4 4 2}
{cmd:group()}, {cmd:invariance()}, {cmd:metric()}, {cmd:standardized}, and the options for
standard errors, intervals and thresholds are not options of {cmd:efacompare} (they belong to one
fit, or are not available with {cmd:efa()}). Using one gives Stata's error 198.


{marker remarks}{...}
{title:Remarks}

{p 4 4 2}
{bf:The geomin epsilon changes with m.} Unless you give {opt epsilon()}, the rotation uses 1e-4
at m = 2, 1e-3 at m = 3 and 1e-2 at m of 4 or more (m = 1 has no rotation). The rotated columns
({cmd:minima}, {cmd:gap}, {cmd:lam3}, {cmd:phimax}) are then not computed on one footing across
rows. The table says so. With {opt epsilon()} one value is used for every m.

{p 4 4 2}
{bf:The rotated diagnostics are conditional on the selected minimum.} Geomin can have several
local minima. Each row's {cmd:lam3} and {cmd:phimax} describe the minimum that the fit
selected, the one with the smallest criterion among the starts that converged. {cmd:minima} and
{cmd:gap} say whether other minima were found.

{p 4 4 2}
{bf:Over-extraction signs are descriptive.} A small {cmd:lam3}, a large {cmd:phimax}, several
minima with a small {cmd:gap}, {cmd:comm1} equal to 1, or a failed rotation can accompany an m
that is too large, and none is a test of it. No column was shown to identify the number
of factors, and none comes with a cutoff.

{p 4 4 2}
{bf:Conventional cutoffs were not calibrated for WLSMV EFA.} The usual RMSEA and CFI thresholds
have not been established for this estimator and these models, and {cmd:dwlsmv} has not checked
them. Do not read them as validated for these columns.

{p 4 4 2}
{bf:Rows that did not converge.} If the unrotated fit at an m does not converge, the row shows
m, df and {cmd:status} {cmd:not conv} and nothing else. {cmd:dwlsmv} posts the last iterate of such
a fit, but a last iterate is not a solution, so no fit statistic or rotated diagnostic is read from
it. If the geomin rotation does not converge from any start, {cmd:dwlsmv} posts nothing, and the row
shows m, df and {cmd:status} {cmd:rot fail}. The command goes on to the next m after either. Any
other failure of {cmd:dwlsmv} at an m stops {cmd:efacompare} with that return code, and the
message gives the {cmd:dwlsmv} command to run to see the failure itself.

{p 4 4 2}
{bf:Time.} A row for an m that is too large for the data can take far longer than the others,
because the fit runs to its iteration limit before it is declared not converged. The command prints
one dot per m as it goes.

{p 4 4 2}
{bf:Results left behind.} {cmd:efacompare} leaves your estimation results as they were:
nothing is left in {cmd:e()}.


{marker examples}{...}
{title:Examples}

{p 4 4 2}
The example data have six three-category indicators, {cmd:y1}-{cmd:y6}, and 5,000 observations
({cmd:dwlsmv_example.dta}, simulated from two correlated factors):

{phang2}{cmd:. use dwlsmv_example, clear}{p_end}
{phang2}{cmd:. efacompare y1-y6}{p_end}

{p 4 4 2}
One epsilon for every m, and a subsample:

{phang2}{cmd:. efacompare y1-y6 if y1 > 0, epsilon(1e-3)}{p_end}

{p 4 4 2}
With six items, m = 3 leaves df = 0; {opt maxfactors()} reaches it. On these data the geomin
rotation at m = 3 converges from none of its starts, as {cmd:dwlsmv} itself reports for that model,
so the m = 3 row is blank and footnoted:

{phang2}{cmd:. efacompare y1-y6, maxfactors(3)}{p_end}

{p 4 4 2}
The table as a matrix, and the eigenvalues:

{phang2}{cmd:. return list}{p_end}
{phang2}{cmd:. matrix list r(table)}{p_end}


{marker results}{...}
{title:Stored results}

{p 4 4 2}
{cmd:efacompare} stores the following in {cmd:r()}:

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}
{synopt:{cmd:r(p)}}number of items{p_end}
{synopt:{cmd:r(maxfactors)}}largest number of factors shown{p_end}

{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}one row per m (row names {cmd:m1}, {cmd:m2}, ...), columns
{cmd:m df status conv T3 p RMSEA CFI TLI minima gap lam3 phimax comm1 ev_m}; {cmd:status} is 1 (fitted),
2 (not converged) or 3 (rotation failed); {cmd:conv} is 1 or 0, and missing when the rotation
failed{p_end}
{synopt:{cmd:r(Ev)}}1 x p eigenvalues of the sample correlation matrix, as {cmd:e(Ev)} of any one
of the fits{p_end}


{marker alsosee}{...}
{title:Also see}

{p 4 13 2}
Help:  {helpb dwlsmv}, {helpb screeplot}, {helpb difftest}
{p_end}
