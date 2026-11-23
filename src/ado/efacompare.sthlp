{smcl}
{* every example run 2026-11-22 on the 0.7.0 release tree (fb6386e plus the version bump), residual diagnostics included, StataMP 18, linesize 80, on build/dwlsmv_example.dta: rc 0 throughout;}
{* the screeplot example included (rc 0); the Scree plot remark checked 2026-11-21 on StataMP 18 (r(301) and the stale e(Ev) as stated).}{...}
{* *! version 0.7.0  22nov2026}{...}
{vieweralsosee "dwlsmv" "help dwlsmv"}{...}
{vieweralsosee "[MV] screeplot" "help screeplot"}{...}
{viewerjumpto "Syntax" "efacompare##syntax"}{...}
{viewerjumpto "Description" "efacompare##description"}{...}
{viewerjumpto "Options" "efacompare##options"}{...}
{viewerjumpto "Remarks" "efacompare##remarks"}{...}
{viewerjumpto "Examples" "efacompare##examples"}{...}
{viewerjumpto "Stored results" "efacompare##results"}{...}
{viewerjumpto "References" "efacompare##references"}{...}
{viewerjumpto "Author" "efacompare##author"}{...}
{viewerjumpto "Acknowledgments" "efacompare##acknowledgments"}{...}
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
diagnostics from the rotated solution that bear on whether m is too large, and the largest
residuals of each fit. It gives no
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
{p2col:{cmd:SRMR}}the standardized root mean square residual, {cmd:e(srmr)}; see
{help efacompare##residuals:Residuals}{p_end}
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
{bf:Third block: largest residuals.} For each fitted m, the three pairs of items with the
largest absolute normalized residual, largest first; ties keep the order of the items.

{p2colset 8 20 22 2}{...}
{p2col:{cmd:m}, {cmd:rank}}the number of factors, and the pair's rank within that m{p_end}
{p2col:{cmd:pair}}the two items, each name abbreviated to 10 characters{p_end}
{p2col:{cmd:residual}}the sample polychoric (tetrachoric) correlation minus the correlation the
model implies, signed: positive when the sample correlation is larger than the implied one; from
{cmd:e(Resid)}{p_end}
{p2col:{cmd:norm.res.}}the residual divided by the estimated standard error of the sample
correlation; from {cmd:e(Nresid)}{p_end}
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

{marker residuals}{...}
{p 4 4 2}
{bf:Residuals: SRMR and the largest residuals answer different questions.} A residual is the
sample polychoric (tetrachoric) correlation of a pair minus the correlation the model implies.
SRMR summarizes all of them in one number, in correlation units: the square root of (the sum of
the squared residuals over the {it:p}({it:p} - 1)/2 pairs, divided by {it:p}({it:p} + 1)/2),
as lavaan computes {cmd:srmr} for categorical indicators. The thresholds do not enter. The
largest-residuals block shows where the most conspicuous discrepancies are, each relative to its
own sampling precision. The two can disagree: a small SRMR can come with one pair the model
reproduces badly, and a pair whose correlation is estimated precisely can have a large normalized
residual and a small residual. SRMR need not fall as m rises, because the fit minimizes a weighted
discrepancy, not the sum of squares SRMR uses.

{p 4 4 2}
{bf:The normalized residual is descriptive, not a test.} It is the residual divided by
sqrt(gamma_kk/N), where gamma_kk is the estimated asymptotic variance of the pair's sample
correlation: the asymptotic standard error of the sample correlation, not of the residual. So
even in large samples it is not standard normal when the model is correct. It serves to rank the
eligible pairs within one row. It is not a formal z-test: no critical value applies to it, 1.96
included, and the table marks none. In sparse data and small samples that standard error can
itself be poorly estimated, for ordinary pairs as well as the ones footnoted.

{p 4 4 2}
{bf:Pairs not ranked.} A pair whose two-way table had one empty cell repaired, or whose
correlation is held at the bound under {opt keepbound}, has no normalized residual (it is
missing in {cmd:e(Nresid)}). It is left out of the ranking and named in a footnote. Its residual
still counts in SRMR and is in {cmd:e(Resid)}. Which pairs these are depends on the data only,
so they are the same in every row.

{p 4 4 2}
{bf:Rows left out, and rows with df = 0.} A row that did not converge, or whose rotation failed,
has no residuals and is left out of the block, with a footnote. A fitted row with df = 0 is kept.
There the model has as many parameters as there are correlations, and T3 and the fit indices are
not available. Residuals near zero mean the model reproduces the correlations exactly; larger
ones mean the model does not reproduce them, as happens when no exact solution exists, which T3
and the fit indices cannot show. No cutoff separates the two, and a footnote says so.

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

{p 4 4 2}
{bf:Scree plot.} Typed without an argument, {helpb screeplot} reads {cmd:e(Ev)}, so after
{cmd:efacompare} it either stops with r(301) or, if an earlier {cmd:dwlsmv} EFA fit is still in
{cmd:e()}, plots that fit's eigenvalues, which can come from other items or another sample. Plot
{cmd:r(Ev)} instead, as in the last example.


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

{p 4 4 2}
A scree plot of the eigenvalues, from {cmd:r(Ev)}:

{phang2}{cmd:. matrix ev = r(Ev)}{p_end}
{phang2}{cmd:. screeplot ev}{p_end}


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
{cmd:m df status conv T3 p RMSEA CFI TLI minima gap lam3 phimax comm1 ev_m srmr}; {cmd:status} is 1 (fitted),
2 (not converged) or 3 (rotation failed); {cmd:conv} is 1 or 0, and missing when the rotation
failed. The largest residuals are shown, not stored; they are in the standalone fit's
{cmd:e(Resid)} and {cmd:e(Nresid)}{p_end}
{synopt:{cmd:r(Ev)}}1 x p eigenvalues of the sample correlation matrix, as {cmd:e(Ev)} of any one
of the fits{p_end}


{marker references}{...}
{title:References}

{p 4 8 2}
Asparouhov, T., and B. Muthén. 2010. Simple second order chi-square correction.
Mplus technical note, May 3, 2010.
{browse "https://www.statmodel.com/download/WLSMV_new_chi21.pdf"}.
{p_end}

{p 4 8 2}
Browne, M. W. 2001. An overview of analytic rotation in exploratory factor analysis.
{it:Multivariate Behavioral Research} 36: 111{c -}150.
{p_end}

{p 4 8 2}
Muthén, B., S. H. C. du Toit, and D. Spisic. 1997. Robust inference using
weighted least squares and quadratic estimating equations in latent variable
modeling with categorical and continuous outcomes. Unpublished manuscript.
{p_end}

{p 4 8 2}
Yates, A. 1987. {it:Multivariate Exploratory Data Analysis: A Perspective on Exploratory}
{it:Factor Analysis}. Albany, NY: State University of New York Press.
{p_end}


{marker author}{...}
{title:Author}

{p 4 4 2}
Erik Ruzek{break}
ORCID: {browse "https://orcid.org/0000-0002-0931-1951"}{break}
Source and issues: {browse "https://github.com/bbachilles/dwlsmv"}
{p_end}


{marker acknowledgments}{...}
{title:Acknowledgments}

{p 4 4 2}
Development of {cmd:dwlsmv} was assisted by Claude Code (Anthropic), predominantly
using Claude Opus models. The author reviewed and tested all code and is responsible
for any errors.
{p_end}


{marker alsosee}{...}
{title:Also see}

{p 4 13 2}
Help:  {helpb dwlsmv}, {helpb screeplot}, {helpb difftest}
{p_end}
