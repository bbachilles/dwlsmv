{smcl}
{* Executed under StataMP 18 on 2026-11-22 on the 0.6.1 release tree (2d22697 plus the version bump): every example in this file and in}{...}
{* difftest.sthlp and efacompare.sthlp run on the source build, net get excepted, including esttab, the two covariate examples and the}{...}
{* second-order example on dwlsmv_example_so.dta (its Note printed, the fit proper); the random-intercept abort at rc 3498 as documented,}{...}
{* every other command rc 0, every fit posting e(version) 0.6.1, no fit triggering recovery; tests/unit/test_dwlsmv_gaps.do G11-G15.}{...}
{* The notes quoted under "Second-order factors" are copied from the source strings; the examples do not print all of them.}{...}
{* *! version 0.6.1  22nov2026}{...}
{vieweralsosee "[SEM] sem" "help sem"}{...}
{vieweralsosee "[R] nlcom" "help nlcom"}{...}
{vieweralsosee "[SEM] gsem" "help gsem"}{...}
{vieweralsosee "difftest" "help difftest"}{...}
{vieweralsosee "efacompare" "help efacompare"}{...}
{viewerjumpto "Syntax" "dwlsmv##syntax"}{...}
{viewerjumpto "Description" "dwlsmv##description"}{...}
{viewerjumpto "Options" "dwlsmv##options"}{...}
{viewerjumpto "What sem reads and what dwlsmv refuses" "dwlsmv##refusals"}{...}
{viewerjumpto "The metric" "dwlsmv##metric"}{...}
{viewerjumpto "Second-order factors" "dwlsmv##secondorder"}{...}
{viewerjumpto "Categories and indicator order" "dwlsmv##categories"}{...}
{viewerjumpto "Measurement invariance across groups" "dwlsmv##invariance"}{...}
{viewerjumpto "Output" "dwlsmv##output"}{...}
{viewerjumpto "Standardized estimates" "dwlsmv##standardized"}{...}
{viewerjumpto "Comparison with lavaan and Mplus" "dwlsmv##comparison"}{...}
{viewerjumpto "Examples" "dwlsmv##examples"}{...}
{viewerjumpto "Stored results" "dwlsmv##results"}{...}
{viewerjumpto "References" "dwlsmv##references"}{...}
{viewerjumpto "Author" "dwlsmv##author"}{...}
{viewerjumpto "Acknowledgments" "dwlsmv##acknowledgments"}{...}
{viewerjumpto "Also see" "dwlsmv##alsosee"}{...}
{p2colset 1 15 17 2}{...}
{p2col:{bf:dwlsmv} {hline 2}}SEM/CFA for ordinal indicators by DWLS with the WLSMV (T3) test{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:dwlsmv} {it:paths} {ifin}
[{cmd:,} {it:options}]

{p 4 4 2}
{it:paths} are {helpb sem} path specifications, for example
{cmd:(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)}.

{p 4 4 2}
A factor may also be measured only by other factors, a {it:second-order factor},
as in {cmd:(G -> F1 F2 F3)} with each {cmd:F} measured by indicators; see
{help dwlsmv##secondorder:Second-order factors}.

{p 4 4 2}
One observed variable may appear as a {it:covariate}: a binary variable coded 0 and 1, written as a
{cmd:sem} path to a factor, {cmd:(}{it:F}{cmd: <- }{it:x}{cmd:)}, or to an indicator,
{cmd:(}{it:y}{cmd: <- }{it:x}{cmd:)}. The model is then fitted conditional on {it:x}: the
thresholds are those at {it:x} = 0, {cmd:(}{it:F}{cmd: <- }{it:x}{cmd:)} is the effect of {it:x}
on the factor, and {cmd:(}{it:y}{cmd: <- }{it:x}{cmd:)} is its direct effect on the indicator's
latent response (minus a threshold shift). A path written with {cmd:@}, as in {cmd:(}{it:F}{cmd: <- }{it:x}{cmd:@.4)},
fixes that coefficient at the value, which is a number in the model's own units, not a
standardized coefficient: {cmd:(}{it:F}{cmd: <- }{it:x}{cmd:@}{it:c)} fixes the effect of {it:x} = 1 against
{it:x} = 0 on the factor at {it:c} factor units (with {cmd:metric(unitvariance)}, {it:c} disturbance standard
deviations; with the marker metric, units of the marker indicator's latent response), and
{cmd:(}{it:y}{cmd: <- }{it:x}{cmd:@}{it:c)} fixes the direct effect on {it:y}'s latent response, whose
conditional variance is 1, at {it:c}, so the threshold shift is {it:-c}. A fixed path is shown in the table, under
Covariates, as its value followed by {cmd:(constrained)}, with no standard error and a row in {cmd:e(Cns)}; it is not an estimated parameter, so it
adds one degree of freedom. {cmd:(}{it:F}{cmd: <- }{it:x}{cmd:@0)} is dropped by {cmd:sem} and adds no row: it is the
same model as no effect of {it:x} on {it:F}. A path fixed at a value other than 0 also fixes the factor's
orientation under {cmd:metric(unitvariance)}: the loadings are not reflected to make the first one positive. A limitation follows: if the sign of the fixed value
is opposite to the orientation the starting values imply, the fit can report convergence at a local optimum whose
test statistic is much larger than that of the same model with the opposite sign, whose mirror-image solution
(every loading reversed) would fit equally well. To avoid it, fix the value with the other sign and reverse the
interpretation of the factor, or use the marker metric, where the orientation is set by the marker loading, fixed
at 1. This was observed on one test dataset; how often it happens has not been studied.
{cmd:sem} lets an observed exogenous {it:x} covary
with every exogenous factor; {cmd:dwlsmv} does not take that covariance as a parameter, so for each
exogenous factor write the path from {it:x} or fix the covariance at zero with the {cmd:sem} option
{cmd:cov(}{it:x}{cmd:*}{it:F}{cmd:@0)}.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt metric(marker|unitvariance)}}how factors the user did not
normalise are scaled; default {cmd:metric(marker)}{p_end}
{synopt:{opth group(varname)}}fit the model in each group defined by
{it:varname}; requires {opt invariance()}{p_end}
{synopt:{opt invariance(level)}}equality across groups: {cmd:configural},
{cmd:thresholds}, {cmd:metric} or {cmd:scalar}; requires {opt group()} and
{cmd:metric(unitvariance)}; see
{help dwlsmv##invariance:Measurement invariance}{p_end}
{synopt:{opt efa(namelist)}}fit the model as an exploratory factor analysis
of these factors, with oblique geomin rotation{p_end}
{synopt:{opt eps:ilon(#)}}geomin epsilon; default depends on the number of
factors{p_end}
{synopt:{opt rstarts(#)}}number of rotation starts; default is
{cmd:rstarts(30)}{p_end}
{synopt:{opt se}}with {opt efa()}, display the standard errors of the rotated
loadings and factor correlations{p_end}
{synopt:{opt ci}}with {opt efa()}, display Wald confidence intervals for them{p_end}
{synopt:{opt thresholds}}display the thresholds, which are hidden by default
unless {opt invariance()} constrains them{p_end}
{synopt:{opt nothresholds}}hide the thresholds{p_end}
{synopt:{opt keepbound}}keep polychoric pairs at the bound instead of
aborting; see {help dwlsmv##keepbound:keepbound}{p_end}
{synopt:{opt standardized}}display the standardized solution, with the
residual variances; see {help dwlsmv##standardized:Standardized estimates}{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{it:sem_options}}options that change how {cmd:sem} reads the model;
see {help dwlsmv##semopts:sem options}{p_end}
{synoptline}

{p 8 15 2}
{cmd:dwlsmv} [{cmd:,} {opt standardized} {opt l:evel(#)} {opt se} {opt ci} {opt thresholds} {opt nothresholds}]

{p 4 6 2}
Typing {cmd:dwlsmv} without a model redisplays the header and coefficient
table of the last {cmd:dwlsmv} fit, standardized if you ask, with the display
options above (the notes block is not redisplayed).{p_end}


{marker description}{...}
{title:Description}

{p 4 4 2}
{cmd:dwlsmv} fits structural equation and confirmatory factor models to ordinal
and binary indicators by diagonally weighted least squares on thresholds and
polychoric correlations, with the scaled-and-shifted (T3) mean- and
variance-adjusted test statistic. The asymptotic covariance matrix of the
sample statistics is estimated empirically.

{p 4 4 2}
The model is written in {helpb sem}'s path syntax, and {cmd:sem} itself reads
it: {cmd:dwlsmv} runs {cmd:sem} {it:paths}{cmd:,} {it:sem_options}
{cmd:noestimate} on a dummy dataset that holds only the model's observed
variables (N = 200, generated from one factor), and translates what {cmd:sem}
reports into the model it fits. {cmd:sem} therefore never sees your data, and
its defaults apply exactly as in {cmd:sem}: the first-listed indicator anchors
a factor, exogenous factor covariances are free, varlist shorthand and either
arrow direction work. Your data, sort order, random-number state and active
estimation results are restored afterwards.

{p 4 4 2}
Estimation uses the observations selected by {ifin} with no missing value on
any indicator (listwise deletion).


{marker options}{...}
{title:Options}

{phang}
{opt metric(marker|unitvariance)} chooses how a factor is scaled when you did
not scale it yourself. {cmd:marker}, the default, is {cmd:sem}'s own
normalisation: the first-listed indicator's loading is fixed at 1.
{cmd:unitvariance} fixes the factor's variance at 1 for an exogenous factor
and its {it:disturbance} variance at 1 for an endogenous factor. It never
overrides a constraint you wrote. See {help dwlsmv##metric:The metric}.

{phang}
{opth group(varname)} names the variable whose values define the groups. Its
values must be nonnegative integers; for a string variable, create one with
{cmd:egen newgrp = group(}{it:varname}{cmd:)} or {helpb encode}. Observations
with a missing value of {it:varname} are excluded. Groups are ordered by value,
and the lowest value is the first, or reference, group. At least two groups
must remain in the estimation sample. In this version {opt group()} is
available only together with {opt invariance()}.

{phang}
{opt invariance(level)} fits the model at one level of measurement invariance
across the groups: {cmd:configural} (the same model in every group, nothing
constrained equal), {cmd:thresholds} (thresholds equal), {cmd:metric}
(thresholds and loadings equal) or {cmd:scalar} (thresholds and loadings
equal, factor means estimated in the groups after the first). It requires
{cmd:metric(unitvariance)}. Each command fits one level. These levels are not
{cmd:sem}'s {cmd:ginvariant()} classes, and {cmd:ginvariant()} is not
available. See {help dwlsmv##invariance:Measurement invariance}.

{phang}
{opt efa(namelist)} fits an exploratory factor analysis (EFA). The paths must
be one block in which every factor in {it:namelist} loads on every indicator,
as in {cmd:dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)}, and {it:namelist} must name
every latent variable. dwlsmv fits the block, then rotates the loadings by
oblique geomin. The rotated loadings and factor correlations are shown without
standard errors unless you ask for them with {opt se} or {opt ci}; the
covariance is stored in {cmd:e(V)} either way. See
{help dwlsmv##efa:Exploratory factor analysis}.

{phang}
{opt epsilon(#)} sets the geomin epsilon, a small positive number. The default
follows Mplus: .0001 with two factors, .001 with three and .01 with four or
more. It requires {opt efa()}.

{phang}
{opt rstarts(#)} sets how many starting rotations are tried; the default is
30. It requires {opt efa()}.

{phang}
{opt se} displays the standard errors of the rotated loadings and factor
correlations, as matrices after the estimates. It requires {opt efa()}. The
caveat in {help dwlsmv##efase:Remarks on rotated standard errors} applies to
every number it shows.

{phang}
{opt ci} displays Wald confidence intervals, the estimate plus or minus z times
its standard error with z from {opt level()}, in brackets beside each estimate.
It requires {opt efa()}. With {opt efa()} no intervals are shown unless you give
{opt ci}. The intervals are not truncated, so the interval for a factor
correlation can extend beyond -1 or 1, and the output says so when it does. The
same caveat applies. {opt se} and {opt ci} may be combined.

{phang}
{opt thresholds} and {opt nothresholds} show or hide the thresholds in the
output of any fit. By default they are hidden, and shown when {opt invariance()}
constrains them across groups ({cmd:thresholds}, {cmd:metric} or
{cmd:scalar}). Hiding changes only what is displayed: {cmd:e(b)} and
{cmd:e(V)} are the same.

{marker keepbound}{...}
{phang}
{opt keepbound} changes how the first stage treats a polychoric correlation
whose estimate stops at the bound of the search. By default such a pair makes
{cmd:dwlsmv} abort with rc 3498, naming the pair. Under {opt keepbound} the
remaining status-2 pairs (see the first-stage message) are kept in the fit at
the bound, as lavaan does, and listed by variable name in the notes block. It
does not admit a 2x2 table with two empty cells on a diagonal or any other
failure status, and it has nothing to do with the one-empty-cell repair, which
is applied whether or not {opt keepbound} is set (see
{help dwlsmv##comparison:Comparison with lavaan and Mplus}). Its effect on standard
errors and on T3 is not established. If a fit aborts without {opt keepbound}
and the same data would pass the first stage with it, {cmd:dwlsmv} prints
{cmd:keep_bound = 1 above is the dwlsmv option keepbound.}

{marker semopts}{...}
{phang}
{it:sem_options} are the {cmd:sem} options that change how the model is read,
such as {cmd:covariance()}, {cmd:variance()}, {cmd:latent()} and
{cmd:nocapslatent}. They are passed unchanged to {cmd:sem} when it reads the
model on the dummy data (with {cmd:noestimate}, and under
{cmd:metric(unitvariance)} also with {cmd:forcenoanchor} in a first pass).
{ifin} is not passed to {cmd:sem}.

{phang}
{bf:sem's estimation options are refused} with rc 198: {cmd:vce()},
{cmd:method()}, {cmd:iterate()}, {cmd:technique()}, {cmd:tolerance()},
{cmd:ltolerance()}, {cmd:nrtolerance()}, {cmd:nonrtolerance}, {cmd:from()},
{cmd:startvalues()}, {cmd:startgrouping()}, {cmd:noivstart} and
{cmd:difficult}. {cmd:dwlsmv} estimates by DWLS, and its standard errors are
always the robust (sandwich) DWLS ones. In particular {cmd:vce(cluster} {it:id}{cmd:)}
is not available: clustering and sampling weights are planned for a later version.
{cmd:nolog} is accepted and has no effect.

{phang}
{opt level(#)} sets the confidence level of the intervals in the table, on
the fit and on replay, with or without {opt standardized}; see
{helpb estimation options##level():[R] Estimation options}.

{phang}
{opt standardized} replaces the coefficient table by the standardized
solution, with delta-method standard errors and one residual-variance row
{cmd:/var(e.}{it:y}{cmd:)} per indicator. It works on the fit and on replay.
Nothing is stored: {cmd:e()} is the same with or without it. See
{help dwlsmv##standardized:Standardized estimates}.


{marker refusals}{...}
{title:What sem reads and what dwlsmv refuses}

{p 4 4 2}
{bf:sem's own refusals} reach you as {cmd:sem}'s message and return code
(rc 498). Two of them enforce rules of this command: {cmd:var(F2@1)} on an
endogenous factor (a total variance) and any covariance with an error term,
such as {cmd:cov(F1*e.F2)}. The latter is {cmd:sem}'s general rule that no
error term may covary with an exogenous variable; the reasons are {cmd:sem}'s.
(The covariance of {cmd:F1} with {cmd:F2}'s disturbance alongside the path
{cmd:F2 <- F1} would also be rank-deficient, but that is not why {cmd:sem}
refuses it.) A disturbance covariance between two endogenous factors,
{cmd:cov(e.F2*e.F3)}, is accepted.

{p 4 4 2}
{bf:dwlsmv's own refusals} exit with rc 198 and a message. {cmd:dwlsmv}
refuses, among models that {cmd:sem} accepts:

{p2colset 8 12 12 2}{...}
{p2col:o}{opt group()} without {opt invariance()};{p_end}
{p2col:o}observed exogenous variables other than one binary 0/1 covariate
({cmd:(}{it:F}{cmd: <- }{it:x}{cmd:)} on a factor, {cmd:(}{it:y}{cmd: <- }{it:x}{cmd:)} on an
indicator), which fits the model conditional on {it:x}: there the statistic,
degrees of freedom and p-value are reported (not when the fit did not converge, see
{help dwlsmv##output:Output}), but no baseline model and no RMSEA, CFI or TLI (these would need their own derivation and validation),
{opt standardized}, {opt keepbound}, {opt efa()} or {helpb difftest} (difference tests with a covariate are refused until their finite-sample size has been studied);{p_end}
{p2col:o}a model with no latent factor;{p_end}
{p2col:o}a constraint that ties several parameters together (for example a
shared {cmd:@}{it:name}); only constraints fixing one parameter at a value are
accepted;{p_end}
{p2col:o}a constrained intercept ({cmd:_cons}); thresholds replace intercepts
for ordinal indicators;{p_end}
{p2col:o}a path between observed variables (an indicator on an indicator), or
a path from an observed variable that is not the covariate to a factor;{p_end}
{p2col:o}a constrained indicator residual variance {cmd:var(e.}{it:y}{cmd:)}:
residual variances are not parameters here, being
1 - (Lambda Phi Lambda')_jj, derived from the model;{p_end}
{p2col:o}covariances among indicator residuals (planned for a later
version);{p_end}
{p2col:o}a variance or covariance of an observed variable; of the covariate,
a constrained variance or mean, a covariance with a factor that is free (see
the covariate in the syntax) or fixed at a value other than 0;{p_end}
{p2col:o}more than one observed exogenous variable, a covariate together with
{opt group()}, and a covariate that is not numeric, not coded 0 and 1 in the
estimation sample, or constant;{p_end}
{p2col:o}a factor with no observed indicator and no path to another factor;
a path out of such a factor fixed at a nonzero value (a second-order factor is
scaled by its variance, {cmd:var(G@1)}); and a second-order factor together
with {opt group()}, {opt efa()} or a covariate;{p_end}
{p2col:o}a factor scaled by a loading or variance fixed at a value other than 1
(for example {cmd:y2@.5} or {cmd:var(F@2)}): a factor's scale is set by a
loading fixed at 1 or a variance fixed at 1;{p_end}
{p2col:o}a parameter {cmd:sem} reports that {cmd:dwlsmv} cannot read;{p_end}
{p2col:o}a model naming no observed variable in the data, and an invalid
{opt metric()};{p_end}
{p2col:o}a non-numeric (string) indicator; the message names it;{p_end}
{p2col:o}{cmd:sem}'s estimation options, {cmd:vce()} among them (see
{help dwlsmv##semopts:sem options}).{p_end}
{p2colreset}{...}

{p 4 4 2}
An estimation sample with no observations gives rc 2000. A loading fixed at 1
together with {cmd:var(F@1)} is accepted.

{p 4 4 2}
{bf:Aborts during estimation} keep their own messages. A first-stage
(threshold and polychoric) abort is rethrown as rc 3498, followed by a legend
mapping the positional labels {cmd:y1}, {cmd:y2}, ... in the message to your
variable names, a line pointing to
{help dwlsmv##comparison:Comparison with lavaan and Mplus}, and, when it
applies, the {opt keepbound} line above. When the model-identification step
aborts because a factor has no scale, {cmd:dwlsmv} adds, for each such factor,

{p 8 8 2}
{cmd:Factor F has no fixed scale: fix one of its loadings at 1,}{break}
{cmd:e.g. (F -> y@1 ...), or specify metric(unitvariance).}

{marker unanchored}{...}
{p 4 4 2}
{bf:A factor that sem leaves unanchored.} In

{p 8 8 2}
{cmd:(F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1)}

{p 4 4 2}
(a random intercept {cmd:G}), {cmd:sem} anchors neither {cmd:F1} nor
{cmd:F2}: its constraints are only {cmd:G}'s six. {cmd:dwlsmv} keeps that
reading, so the fit aborts with the no-scale message for {cmd:F1} and
{cmd:F2}. Either anchor {cmd:F1} and {cmd:F2} yourself, or specify
{cmd:metric(unitvariance)}, which leaves {cmd:G} alone (you normalised it) and
fixes the variances of {cmd:F1} and {cmd:F2}. {cmd:G} orthogonal to {cmd:F1}
and {cmd:F2} must be written ({cmd:cov(G*F1@0) cov(G*F2@0)}); {cmd:sem}'s
default frees those covariances. With unit variance and the {cmd:F1}-{cmd:F2}
covariance free, this model loses rank at the default starting values only;
{cmd:dwlsmv} accepts it and starts the fit at a second point.


{marker metric}{...}
{title:The metric: marker versus unit variance}

{p 4 4 2}
{cmd:metric(marker)} is {cmd:sem}'s normalisation, and lavaan's and Mplus's
default too: each factor's first-listed indicator has its loading fixed at 1.

{p 4 4 2}
{cmd:metric(unitvariance)} applies only to factors you did not normalise (by a
fixed loading or a fixed variance; a fixed covariance is not a normalisation).
For each of them {cmd:dwlsmv} adds {cmd:var(F@1)}, or {cmd:var(e.F@1)} when
{cmd:F} is endogenous, and has {cmd:sem} read the model again.

{p 4 4 2}
{bf:Unit variance fixes the disturbance variance psi of an endogenous factor,}
{bf:not its total variance.} The factor covariance matrix is derived,
Phi = (I-B)^-1 Psi (I-B)^-T. For an exogenous factor psi and Phi are the same
number. For an endogenous factor they differ by the explained variance, and
its total variance is not 1. On the example data ({cmd:dwlsmv_example.dta}, two factors,
{cmd:F2 <- F1}), lavaan 0.7.2 with unit-variance identification gives
Var(F2) = beta^2 + 1 = 1.2800 at beta = 0.5291, while Var(F1) = 1. This is
also {cmd:sem}'s convention: it refuses {cmd:var(F2@1)} for an endogenous
factor and requires {cmd:var(e.F2@1)}. Consequently {cmd:metric(unitvariance)}
is not a standardized solution: it changes which parameters are estimated,
not how a fitted solution is displayed, and an endogenous factor's loadings
still change on standardization.

{p 4 4 2}
A factor with no observed indicator (a second-order factor) is scaled by its
variance under either metric; see
{help dwlsmv##secondorder:Second-order factors}.

{p 4 4 2}
The header states the metric: "marker (first-listed indicator)" or "factor
variances fixed at 1 (for a factor predicted by other factors, its residual
variance)".

{marker recovery}{...}
{p 4 4 2}
{bf:Recovery of a marker-metric fit.} Under {cmd:metric(marker)} the free
loadings start from the same positive value. When an indicator loads
negatively (a reverse-keyed item), the fit can stop at a degenerate point that
has not converged, although the same model under {cmd:metric(unitvariance)}
reaches the proper solution.

{p 4 4 2}
{cmd:dwlsmv} therefore checks the first marker-metric fit. If it did not
converge, if the admissibility check flags a condition, or if the information
in the fit was insufficient to determine all parameters reliably (the bread of
the sandwich is numerically rank deficient), {cmd:dwlsmv} fits the same model
under unit variance on the same first-stage results, transforms that solution
to the marker metric, and fits the marker model again from that point. With
c_F the marker indicator's loading under unit variance, loadings are divided
by c_F, a regression of F on G is multiplied by c_F/c_G, and a (co)variance
psi_FG by c_F*c_G. The unit-variance fit only supplies starting values: the
reported fit is the final marker fit, and nothing is chosen by comparing fits.

{p 4 4 2}
Recovery applies only to {cmd:metric(marker)}. It is attempted only when the
transformation is exact for the model: every factor has exactly one fixed
nonzero loading, equal to 1; no factor (disturbance) variance is fixed; and
every other fixed loading, regression or covariance is 0. Otherwise the model
is ineligible, and the note names the reasons.

{p 4 4 2}
{bf:With a covariate} recovery is the same, with the effect of the covariate on
a factor (the {cmd:gamma} of {cmd:(F <- x)}) multiplied by c_F, and the effect
of the covariate on an indicator (the {cmd:kappa} of {cmd:(y <- x)}) and the
thresholds unchanged. A model that fixes such an effect on a factor at a
nonzero value is ineligible, because the value means something different in
the two metrics; one fixed at 0, and a fixed effect on an indicator of any
value, stay eligible. For a model with a covariate the first marker fit is
recovered when it did not converge, the admissibility check flags a
condition, or the weighted derivative matrix is numerically rank deficient
at the estimate: the same three causes as without a covariate.

{p 4 4 2}
Recovery succeeds when the final marker fit converges; its admissibility
check and diagnostic notes are reported as for any fit, so a genuine Heywood
case remains flagged. It fails when the unit-variance fit does not converge,
when some c_F is 0 or not finite, when the final fit does not converge, or
when any step of it stops with an error. The first marker fit is then
reported unchanged, with its own warnings.

{p 4 4 2}
{cmd:e(recovery)} records the outcome: {cmd:none} (not needed, or
{cmd:metric(unitvariance)}), {cmd:recovered}, {cmd:ineligible} or
{cmd:failed}. {cmd:e(converged)} and {cmd:e(iterations)} are those of the
reported fit. There is no option to turn recovery off.

{p 4 4 2}
When recovery is {cmd:ineligible} or {cmd:failed}, refit with
{cmd:metric(unitvariance)}, or reverse-code the negatively loading items, or
choose as marker an indicator that loads positively (list it first).


{marker secondorder}{...}
{title:Second-order factors}

{p 4 4 2}
{bf:What one is.} A factor with no observed indicator and a path to other
factors is a second-order factor: {cmd:(G -> F1 F2 F3)}, where each {cmd:F} is
measured by indicators, as in {cmd:(F1 -> y1 y2 y3)}. Write it as any {cmd:sem}
path. Second-order factors are available for a single group only.

{p 4 4 2}
{bf:Scale.} A second-order factor is scaled by its variance, not by a
loading. Under {cmd:metric(marker)}, the first-order factors are scaled by their
marker indicators as usual, and {cmd:dwlsmv} adds {cmd:var(G@1)} for each
second-order factor you did not normalise ({cmd:var(e.G@1)} if {cmd:G} is itself
predicted by another factor), has {cmd:sem} read the model again, and prints one line:

{p 8 8 2}
{cmd:Note: G has no observed indicator and is scaled by its variance, var(G@1); sem's default fixes the path to its first child at 1 instead. Standardized estimates are the same.}

{p 4 4 2}
With several such factors the line names them together, in the plural. If you
write {cmd:var(G@1)} yourself, the model is the same and there is no line.
{cmd:sem}'s default, the path to the first child fixed at 1, is not used. Under
{cmd:metric(unitvariance)} the factor's variance is fixed at 1 as for any factor
you did not normalise, and there is no such line. A path out of {cmd:G} fixed at a nonzero
value is refused (below), because it would be a second scale for {cmd:G}.

{p 4 4 2}
{bf:What is accepted.} {cmd:G} over three or more first-order factors; two
correlated second-order factors; a first-order factor outside the hierarchy,
which may covary with {cmd:G} or be regressed on it; a first-order factor with
two indicators; and an indicator that also loads on {cmd:G}, as in 0.5.0 (a
factor with an indicator is scaled as any such factor is, by its marker or
under {cmd:metric(unitvariance)}).

{p 4 4 2}
{bf:What is refused, with rc 198.} A second-order factor together with
{opt group()}, {opt efa()} or a covariate; a path out of a factor with no
indicator fixed at a nonzero value, for example {cmd:(G -> F1@1 F2 F3)}, with a message that names {cmd:var(G@1)};
and a factor with no indicator and no path to another factor. The same message
for a nonzero fixed path is given under both metrics, and whether or not you
wrote {cmd:var(G@1)}.

{p 4 4 2}
{bf:What is refused as not identified.} The model-identification check refuses
some models that {cmd:dwlsmv} otherwise accepts: {cmd:G} over exactly two
first-order factors with no further constraint, a first-order factor with a
single indicator, and two uncorrelated second-order factors. The refusal exits
with rc 3498; it says that the model appears not to be identified, how many free
parameters it has and how many the data seem to determine, and lists the
parameters that seem to trade off against one another. These directions are a diagnostic
at trial values, not a proof that the model is unidentified.

{p 4 4 2}
{bf:G over exactly three.} Such a model fits exactly as well as the model in
which the three first-order factors are simply correlated, since three
correlations are carried by three paths. When the fit is improper, that is,
some disturbance variance of a first-order factor is negative, {cmd:dwlsmv}
adds this to the admissibility note, after the lines about the variances:

{p 8 8 2}
{cmd:G has exactly three first-order factors (F1, F2, F3), so this model fits exactly as well as one in which those three factors are simply}
{cmd:correlated. The negative disturbance variance means their correlations cannot be represented by a proper second-order factor.}

{p 4 4 2}
Under {cmd:metric(unitvariance)}, a fit with a second-order factor that does not
converge is followed by this note:

{p 8 8 2}
{cmd:This model has a second-order factor. With metric(unitvariance) its first-order factors cannot have a negative disturbance variance; with}
{cmd:metric(marker) they can, so a fit that needs one may converge there.}

{p 4 4 2}
{bf:The disturbance-variance check.} For every free diagonal element of Psi, in
every model and not only second-order ones, {cmd:dwlsmv} classifies the
estimate against a numerical band computed at the solution: a value within the
band, of either sign, is reported as zero within numerical error and is not
called inadmissible; below the band it is reported as negative. A fixed variance and a fit that did not
converge are not checked. The lines are printed inside the admissibility note.
Where the standardized loading is defined by a single path (one group, one
structural parent of the factor), the line gives it:

{p 8 8 2}
{cmd:F1: disturbance variance -0.3200 is negative; its standardized loading on G is 1.2600, above 1.}{break}
{cmd:F1: disturbance variance is zero within numerical error; its standardized loading on G is 1.}

{p 4 4 2}
An exogenous factor's line says "variance", not "disturbance variance". If the band cannot be computed (a
singular information matrix), the line gives the sign only and says so. With no other
admissibility finding, a negative line brings the note's own header, and a line that is only
zero within numerical error gets the header {cmd:Disturbance variances:}. The check does not trigger
{help dwlsmv##recovery:recovery}.

{p 4 4 2}
{bf:Sign.} After the first-order factors are oriented, a second-order factor is
reflected when the sum of its standardized paths to its first-order factors
is negative, so that the paths are reported with a positive sum under both
metrics; a sum of exactly 0 leaves it as fitted. Each second-order factor is
treated separately. A second-order factor with an indicator keeps the rule for any factor with indicators.

{p 4 4 2}
{bf:Recovery.} Under {cmd:metric(marker)} a second-order factor takes part in
{help dwlsmv##recovery:recovery} like any other, with c_G = 1: its variance is 1
under both metrics. When recovery fails and the converged first fit has a
disturbance variance classed negative, the note says why:

{p 8 8 2}
{cmd:The marker-metric fit produced an invalid solution: F1 has a negative disturbance variance. A restart from the unit-variance solution,}
{cmd:transformed to the marker metric, was attempted but failed. The first fit is retained. Its negative disturbance variance cannot be represented}
{cmd:by the unit-variance parameterization, although the first fit itself converged.}

{p 4 4 2}
It is followed by a line naming the step that failed. "Variance" replaces
"disturbance variance" for a factor with no path into it, and several factors are listed.

{p 4 4 2}
{bf:Standardized estimates} are available for a second-order model, with the
paths from {cmd:G} to its first-order factors standardized as any structural
path (see {help dwlsmv##standardized:Standardized estimates}); they agree with
lavaan's {cmd:std.all}.

{p 4 4 2}
{bf:Comparing with {cmd:difftest}.} {helpb difftest} between a model with a
second-order factor and one without, for example its correlated-factor parent,
is refused, with rc 3498:

{p 8 8 2}
{cmd:dwlsmv check 8 (}{it:name0} {cmd:in} {it:name1}{cmd:): one model has a second-order factor and the other does not. DIFFTEST across}
{cmd:these structures (for example a second-order model against its correlated-factor parent) is not available in this version; it is planned}
{cmd:once its reference distribution is established (spec Module 12, B3.9).}

{p 4 4 2}
This is a limit of this version. Each such comparison needs its nesting, its
identification and its reference distribution established before it is offered.
A pair of second-order models is judged as any pair is.

{p 4 4 2}
{bf:Not yet available.} {opt group()} and {opt invariance()}, {opt efa()}, a
covariate, and equality constraints on the paths from a second-order factor.


{marker categories}{...}
{title:Categories and indicator order}

{p 4 4 2}
{bf:Category recoding.} For each indicator, the distinct values observed in
the estimation sample are sorted and recoded to 0, 1, ..., K-1, where K is the
number of distinct observed values. Codes 1, 2, 4 become 0, 1, 2. Value labels
are not read, and a gap in the codes is not treated as an empty category.

{p 4 4 2}
{bf:Indicator order} is {cmd:sem}'s {cmd:e(oyvars)}: observed indicators in
order of first appearance in the paths as written, not dataset order. A
cross-loaded indicator keeps its first position, and a factor whose paths are
split across pieces is not regrouped. This order fixes the positional labels
{cmd:y}{it:j} in first-stage messages, the threshold block of {cmd:e(b)}, and
{cmd:e(indicators)}.


{marker invariance}{...}
{title:Measurement invariance across groups}

{p 4 4 2}
With {opth group(varname)} and {opt invariance(level)}, {cmd:dwlsmv} fits the
same model in every group, with the equalities across groups that the level
sets. Groups are ordered by the value of {it:varname}; the first (lowest
value) is the reference group. In the reference group every factor has
variance 1 and mean 0, every indicator's latent response has intercept 0 and
scale 1, and nothing else is fixed.

{p 4 4 2}
What each level makes equal, and what the groups after the first then
estimate:

{p2colset 8 24 26 2}{...}
{p2col:{cmd:configural}}nothing is equal. Factor variances 1 and means 0 in
every group.{p_end}
{p2col:{cmd:thresholds}}thresholds equal. Each indicator's latent response
may shift (every item) and change scale (items with three or more
categories). Factor variances 1 and means 0.{p_end}
{p2col:{cmd:metric}}thresholds and loadings equal. Latent-response shifts and
scales for every item; a factor's variance is estimated if it has an
indicator with three or more categories, else fixed at 1. Means 0.{p_end}
{p2col:{cmd:scalar}}thresholds and loadings equal; no latent-response shifts.
Latent-response scales, factor variances and factor means are
estimated.{p_end}
{p2colreset}{...}

{p 4 4 2}
These are the levels of Wu and Estabrook's identification for ordinal
indicators, as {cmd:dwlsmv} generates them; they are not {cmd:sem}'s
{cmd:ginvariant()} classes. Each command fits one level; fit the levels you
need as separate commands and store each with {helpb estimates store}.

{p 4 4 2}
{bf:Comparing levels.} Each level reports its own T3 test of fit. The
difference between two T3 statistics is not a valid test of the difference
between two levels. To test one level against another, use {helpb difftest},
which also states the conditions under which its check of the refit has been
validated for multiple-group models.

{p 4 4 2}
{bf:Requirements.} {cmd:metric(unitvariance)} must be specified, and no
loading may be fixed at 1. The model must be the same in every group: a path
written for some groups only, or fixed in one group and free in another, is
refused. Structural paths between factors, {cmd:partial()} and
{cmd:ginvariant()} are not available with {cmd:invariance()}.

{p 4 4 2}
{bf:Categories.} At every level past {cmd:configural} each indicator's
thresholds are equal across groups, so each indicator must take the same
values in every group. If an indicator has a different number of categories,
or the same number but different values (for example 1, 2, 3 in one group
and 1, 2, 4 in another), {cmd:dwlsmv} lists each such indicator and its
values by group and fits nothing. {cmd:invariance(configural)} is not
affected.

{p 4 4 2}
{bf:Scalar invariance with only binary indicators.} {cmd:invariance(scalar)}
is refused when any factor has only binary indicators; this is checked from
the model and the category counts before anything is estimated. For a factor
with only binary indicators, the scalar model is not identified when the
items are equally difficult (thresholds proportional to loadings). Near this
point, studies on generated data found substantial standard-error inflation,
even in fits that converged. A rule to tell safe data from harmful data
before fitting was studied, but no single value could yet be chosen by a
pre-specified procedure, so this version refuses rather than apply an
arbitrary one; this may change in a later version. {cmd:configural},
{cmd:thresholds} and {cmd:metric} are not affected, nor is {cmd:scalar} when
every factor has an indicator with three or more categories.

{p 4 4 2}
{bf:Levels with the same degrees of freedom.} When a level has as many
degrees of freedom as the level before it, the notes say so. When every
factor has at most one indicator with three or more categories, the metric
level is equivalent to the thresholds level and adds no testable restriction,
and the note says that. In other such cases the note says only that the
degrees of freedom are equal: equal degrees of freedom alone do not make two
models equivalent.

{p 4 4 2}
{bf:Notes with groups.} Notes are printed group by group, on each group's own
estimates and scale. The other diagnostic notes and the category counts are
not printed with groups in this version; see {help dwlsmv##output:Output}.


{marker efa}{...}
{title:Exploratory factor analysis}

{pstd}
{cmd:dwlsmv (F1 F2 -> y1-y6), efa(F1 F2)} fits a two-factor exploratory factor
analysis of six ordinal indicators. The factor names are up to the user
({cmd:F1} and {cmd:F2} are used here as examples); the block
{cmd:(F1 F2 -> y1-y6)} says that both factors load on all six indicators, and
{cmd:efa(F1 F2)} marks the block as exploratory. It is the counterpart of
Mplus's {cmd:F1-F2 BY y1-y6 (*1)}.

{pstd}
{it:What is fitted.} The factors' variances are fixed at 1, and enough loadings
are fixed at 0 to identify the model; that solution is then rotated. The
rotation does not change the model-implied correlations, so the test
statistic, its degrees of freedom and the fit indices are those of the
unrotated fit, and do not depend on the rotation.

{pstd}
{it:The rotation.} Oblique geomin, with no row normalization. The criterion
uses a small constant epsilon, by default .0001 with two factors, .001 with
three and .01 with four or more, as Mplus does; lavaan uses .001 throughout,
so its loadings differ slightly from dwlsmv's with two or with four or more
factors unless you set {cmd:epsilon(.001)}. The rotation is run from 30
starting points (see {opt rstarts()}), and the best converged result is
reported.

{pstd}
The starting points are a fixed set, not random draws: the same data give the
same rotated solution every time, and {cmd:set seed} has no effect on it.
dwlsmv never uses or changes the random-number state. A larger
{opt rstarts()} keeps the first 30 starts and adds more, so it can only find
an equal or better solution.

{marker efaorder}{...}
{pstd}
{it:Sign and order.} A rotation determines the factors only up to their signs
and order. dwlsmv reverses any factor whose loadings sum to a negative number
(a sum of exactly zero is left as it is), then orders the factors by where their
largest loadings fall in the indicator list, as lavaan does: for each factor it
takes the mean position of the indicators whose loading is at least 0.8 times the
factor's largest absolute loading, and puts the factor with the smaller mean first.
List the indicators in the order you expect the factors, and the factors come out in
that order. The names in {opt efa()} are attached to the factors in that final order.

{pstd}
{it:Ties in the order.} The order is fully determined. If two factors have exactly
the same mean position, dwlsmv decides by further rules, each used only when the
previous one does not separate them: the factor with the larger sum of squared
loadings first; then the factor with the larger loading at the first indicator
where their loadings differ; then, if their loadings are identical, the order that
gives the larger factor correlations, compared along the lower triangle row by row.
If even that does not separate them, the factors have the same loadings and
correlations, so the order changes no estimate. Whenever the first rule ties, the
output names the two factors and the rule that decided. The order never depends on
which starting point of the rotation was reported. Reordering the factors reorders the
rows and columns of the loadings, the factor correlations and {cmd:e(V)} accordingly
and changes nothing else.

{pstd}
{it:Scree plot.} After {cmd:efa()}, {cmd:e(Ev)} holds the eigenvalues of the
polychoric (or tetrachoric) correlation matrix of the indicators, largest first,
and {cmd:e(properties)} includes {cmd:eigen}, so Stata's {helpb screeplot} plots
them: type {cmd:screeplot} after the fit. They are not printed. They come from the
correlation matrix alone, so they are the same for any number of factors and do not
depend on the rotation. The matrix is estimated pair by pair and need not be
positive definite, so an eigenvalue can be negative; none is clipped. The plot is
a picture to look at when choosing the number of factors; dwlsmv applies no rule
to it. To see the fit and the rotated diagnostics for several numbers of factors side by side,
see {helpb efacompare}.

{pstd}
{it:What is reported.} The rotated loadings and the factor correlations, the test
statistic and the fit indices. The thresholds are stored but not displayed unless
you give {opt thresholds}. Standard errors and intervals are displayed with
{opt se} and {opt ci}. The covariance of the rotated estimates is stored in
{cmd:e(V)}, so {cmd:test} and {cmd:lincom} can be used after an EFA; see
{help dwlsmv##efase:Remarks on rotated standard errors} for what those results
mean. When the rotation finds more than one distinct solution among its starts,
the output says how many, and how many starts reached the reported one.

{pstd}
If the covariance cannot be computed, because the unrotated fit did not converge
or is inadmissible, because the rotated model's weighted Jacobian is rank
deficient (smallest to largest singular value below 1e-8), or because the geomin
criterion has a degenerate critical point at the solution, the rotated estimates
are shown without standard errors, {cmd:e(V)} is not stored, the output says why
in error colour, and {cmd:e(rot_nose)} names the reason. The test statistic and the
fit indices are still reported, except when the unrotated fit did not converge:
then, as for any fit that does not converge, no inference is reported (see
{help dwlsmv##unconverged:A fit that does not converge}).

{pstd}
{it:When the rotation does not converge.} If none of the starting points
converges, dwlsmv stops with an error and posts nothing. This usually means
more factors than the data support, often with an unrotated fit that is
inadmissible (an indicator's communality above 1): try fewer factors.

{pstd}
{it:Degrees of freedom.} With p indicators and m factors the model has
p(p-1)/2 - (pm - m(m-1)/2) degrees of freedom. A model with negative degrees
of freedom is refused. With zero, it fits exactly: the test is not informative,
and the fit indices are shown as missing.

{marker efase}{...}
{pstd}
{it:Remarks on rotated standard errors.} What is and is not established.

{pstd}
{it:What is computed.} The covariance of the rotated loadings and factor
correlations is obtained from the covariance of the unrotated estimates by the
delta method, evaluated at the rotated solution that dwlsmv reports. It is a
local first-order approximation, and it is conditional on that solution: the
rotation is not re-run for each sample, and the covariance does not describe the
variation of the rotation itself. Rotation does not change the thresholds, T3 or
the fit indices; their standard errors and test results are those of the
unrotated fit, and the cautions here do not apply to them.

{pstd}
{it:How good the approximation is.} Its accuracy depends on how far a loading is
from zero, measured against the square root of the geomin epsilon (.01 at the
default with two factors, about .032 with three), and on the sample size. Loadings
near zero are where the approximation is least accurate. Beyond this, no rule is
established for deciding which standard errors can be trusted: dwlsmv does not
flag any loading and applies no cutoff.

{pstd}
{it:What was tested.} The standard errors were checked in a simulation study of
six cases: three populations, each at N = 500 and N = 2000. The populations were a
two-factor model with ten ordinal indicators of three categories, a three-factor
model with twelve ordinal indicators of three categories, and a two-factor model
with ten binary indicators, all at the default epsilon. Coverage below is over
about 4,000 samples per estimate. Nothing is claimed beyond these cases: not for
larger N, for more than three factors, for other epsilon values, or for designs
where every indicator has sizeable cross-loadings.

{pstd}
{it:What the simulation showed.}

{phang2}
{c -} Of 84 clearly nonzero loadings and factor correlations, 95% intervals from
the standard errors covered the true value in 92.6% to 97.5% of the samples. Seven
of them, which covered 92.6% to 93.8%, could not be classified as calibrated or
not by the study's rule; the others were classified as calibrated. None was
classified as miscalibrated.

{phang2}
{c -} Near zero, coverage was lower and uneven. Of 36 loadings that were exactly
zero in the population, 34 were covered in 93% or more of the samples and two
(three-factor model, N = 2000) in 88%. Of 36 loadings set at 0.5 to 6 times the
square root of epsilon, coverage ranged from 79% to 98%, and 18 were below 92.5%.
In the two-factor populations, loadings at 2, 4 and 6 times the square root of
epsilon (.02, .04 and .06) were covered in 79% to 91% of the samples, at both
N = 500 and N = 2000.

{phang2}
{c -} Four candidate measures of closeness to zero were tried as a way of saying
which standard errors to trust. None was shown to work: three were rejected and
the fourth could not be decided, so dwlsmv makes no flag, cutoff or narrowing of
the claim.

{pstd}
So: for loadings and factor correlations that are clearly nonzero, intervals from
these standard errors performed close to nominal in the cases tested. For loadings
near zero, treat standard errors and intervals as approximate, with coverage
sometimes well below nominal. The size of the shortfall in the tested cases does
not carry over to other models.

{pstd}
{it:Several minima.} Geomin rotation can have more than one local minimum. When
the starts reach more than one, the output says how many, and
{cmd:e(rot_minima)}, {cmd:e(rot_atbest)} and {cmd:e(rot_gap)} hold the same
facts. The standard errors are conditional on the reported solution: another
minimum, when it exists, has its own standard errors, which can differ from those
reported. Whether having several minima goes with poorer coverage has not been
studied.

{pstd}
{it:Not available with efa() in this version:} other paths, a second or
partial block, constraints on the loadings or factor covariances,
{opt group()}, {opt invariance()}, {opt metric()} (the block fixes the factor
variances itself), {opt standardized} (the rotated loadings are already
standardized), and {helpb difftest}.

{marker output}{...}
{title:Output}

{p 4 4 2}
The header gives the number of observations, the metric line, T3 with its
degrees of freedom and p-value, the scaled RMSEA, CFI and TLI (see
{help dwlsmv##fitindices:Fit indices}). When the fit did not converge it says so
first and reports no T3, p-value or fit index (see
{help dwlsmv##unconverged:A fit that does not converge}). The coefficient table
follows, in sections: {cmd:Loadings on} {it:F} (one per factor), {cmd:Structural}
({it:F2} {cmd:<-} {it:F1}), {cmd:Covariates} ({it:F} {cmd:<-} {it:x},
{it:y4} {cmd:<-} {it:x}), {cmd:Variances}, {cmd:Covariances}, and, with
{opt invariance()}, the intercepts, scales and factor means of the groups, then
the thresholds, one block per item; an empty section is not printed. The columns
are Estimate, Std. err., z, P>|z| and the Wald interval at {opt level()}
(normal theory). Each row is labelled by the parameter's name without the
leading slash, so it can be found in {cmd:_b[]}. The marker loading shows as
{cmd:1  (marker)}; any other parameter you fixed (for example
{cmd:y2@.8} beside the marker's {cmd:y1@1}, or {cmd:cov(F1*F2@.2)}) shows as its value followed by
{cmd:(constrained)}. Parameters that the metric fixes are not listed: under
{cmd:metric(unitvariance)} the variance of each factor fixed at 1 (the
disturbance variance of a factor predicted by other factors), which the metric
line of the header states. The thresholds are not shown unless
{opt invariance()} constrains them or you give {opt thresholds}; a one-line note
says when they are hidden, and they are always in {cmd:e(b)}. The table is
display only: {cmd:e(b)}, {cmd:e(V)} and {cmd:e(Cns)} are unchanged, and
{cmd:r(table)} (the layout of {cmd:ereturn display}: rows
{cmd:b se z pvalue ll ul df crit eform}) holds every parameter of {cmd:e(b)}, hidden ones included,
after a replay. The table is 78 columns wide, and the label column grows with
the longest label up to {cmd:linesize} - 62; a longer label is cut with a
{cmd:~}.

{p 4 4 2}
After the table comes a single notes block, in this order:

{p2colset 8 12 12 2}{...}
{p2col:1.}if the first marker-metric fit triggered
{help dwlsmv##recovery:recovery}, what triggered it and whether it was
restarted from the unit-variance solution; if recovery was ineligible or
failed, why, and that the first fit is retained;{p_end}
{p2col:2.}under {opt keepbound}, the pairs included at the bound, by variable
name ({it:varA}{cmd:~~}{it:varB});{p_end}
{p2col:}the pairs given the one-empty-cell repair, by variable name with their
repaired correlation, the rule, and a note that its effect on standard errors
and T3 at small N is not established;{p_end}
{p2col:3.}polychoric correlations at or beyond 0.99 in absolute value (not
included at the bound), by variable name with their value, and a note that
such a value has several possible readings (near-duplicate items, strong skew,
too few observations) that the data cannot tell apart;{p_end}
{p2col:4.}the admissibility check of the solution and the diagnostic
notes, printed only when a condition is found;{p_end}
{p2col:5.}on every fit, the category counts: the fewest observations in any
lowest or highest category, and in any middle category, with their items.
These are descriptive only: no threshold is applied, and they are not a
warning or an admissibility check.{p_end}
{p2colreset}{...}

{p 4 4 2}
{bf:With {cmd:invariance()}.} The header names the level and the number of
groups and says that the factor variances are fixed at 1 in the first group
(the reference group). The table has one sub-table per group, titled
{cmd:Group} {it:k} {cmd:(}{it:varname} {cmd:=} {it:value}{cmd:)}, each with
its loadings, variances and covariances and, where the level estimates them,
the latent-response intercepts, scales and factor means; the thresholds
(hidden under {cmd:configural} unless you give {opt thresholds}) follow as one
table, one block per item, each threshold labelled with its group. A parameter
constrained equal across groups appears in every group with the same estimate
and standard error. The parameters a level fixes, which are conventions and not
results, are not listed, and a line under the table says so: the factor
variances fixed at 1, the intercepts fixed at 0, the scales fixed at 1 and the
means fixed at 0. They are in {cmd:e(b)} and {cmd:e(Cns)}. Names follow
{cmd:sem}'s {cmd:group()} form:
{cmd:y1:2.g#c.F} is y1's loading on F in the group with g = 2,
{cmd:/y1:2.g#c.cut1} its first threshold, {cmd:y1:2.g} its latent-response
intercept, {cmd:/scale(y1)#2.g} its latent-response scale and
{cmd:/mean(F)#2.g} the mean of F; the first group's tag reads {cmd:1bn.g}.

{p 4 4 2}
The notes block then lists, under {cmd:Group} {it:g} {cmd:(}{it:varname}
{cmd:=} {it:value}{cmd:):} and only for groups where they occur, the pairs
included at the bound, the one-empty-cell repairs, the correlations at or
beyond 0.99 and the admissibility check, followed by the note on equal degrees
of freedom, if any. The other diagnostic notes and the category counts are not
printed with groups in this version.

{p 4 4 2}
Further numerical diagnostics, ratios of singular values with no established
cutoff, are not printed; they are stored in {cmd:e()} (see
{help dwlsmv##results:Stored results}).

{p 4 4 2}
A one-line-per-indicator legend ({cmd:y3 = q7_anx}) follows any first-stage
message, since only the first stage uses positional labels: after an abort, and
before a fit that continues when a pair was kept at the bound or repaired. The
residual variances (Theta) appear only in the standardized table.

{p 4 4 2}
If the fit does not converge from its ULS starting values, it is run once
more from the model's initial starting values, and a note says so. The criteria
do not change. If that also fails, the first attempt is reported. A
{cmd:metric(marker)} fit takes no such retry; it has
{help dwlsmv##recovery:recovery} instead.

{marker unconverged}{...}
{p 4 4 2}
A fit that does not converge is posted for inspection, with no inference:
inference is reported only for a converged solution, whether the model is a CFA,
an {cmd:efa()} fit, has several groups or a covariate. Under
{cmd:metric(marker)} such a fit is first put through
{help dwlsmv##recovery:recovery}; a fit whose recovery succeeds converges and is
reported as any fit.

{p 4 4 2}
{bf:A fit that does not converge} is reported as follows. A
line in error colour directly under the title says that the optimizer did not
converge, that the estimates are the last iterate and not a solution, and, under
{cmd:metric(marker)}, to refit with {cmd:metric(unitvariance)}, adding that the
recovery refit was not attempted (the model is ineligible) or failed, as the notes
say. Otherwise {cmd:dwlsmv} has no option for other
starting values. No standard errors, T3, p-value or fit index is reported. Only {cmd:e(b)}
(and {cmd:e(Cns_b)}, which marks the fixed parameters) is posted, with
{cmd:e(df_T3)}, {cmd:e(converged)}, {cmd:e(iterations)} and what identifies the model
and the sample; {cmd:e(V)}, {cmd:e(T3)}, {cmd:e(p_T3)}, the fit indices and the
baseline's T3 are not, so {cmd:test}, {cmd:lincom} and {cmd:nlcom} have nothing
to work from. The table lists the last iterate, in the sections of the coefficient table but
with an Estimate column only: the marker loading shows as {cmd:(marker)} and any other fixed parameter as
{cmd:(constrained)}, under the title {cmd:Last iterate, not a solution; no standard errors}. The replay says the same;
{cmd:standardized} is not available for such a fit.

{p 4 4 2}
{bf:Admissibility with a covariate.} The admissibility check is run on every fit
conditional on a covariate, converged or not, at the reported estimates, and its
notes appear in the notes block as for any fit, after the note on
{help dwlsmv##recovery:recovery} when it was triggered.

{p 4 4 2}
{bf:Diagnostics with a covariate.} The notes on the sampling covariance of the
thresholds, slopes and correlations (Gamma), on the identification of the parameters
at the estimate (the weighted derivative matrix) and the category counts are
produced for a model with a covariate as for any fit, and the four stored values
{cmd:e(gamma_ratio)}, {cmd:e(gamma_ratio_raw)}, {cmd:e(gamma_rank)} and
{cmd:e(bread_ratio)} are posted. A response pattern is then the set of responses
together with the value of the covariate. The category counts are shown pooled and,
in addition, within each value of the covariate, with the two group sizes: the
standard error of a slope depends on the counts within each group, which the pooled
counts do not show.
{bf:Limitation.} The check for thresholds that move together on every observation
is not carried out with a covariate. An item that is an exact function of another
item and the covariate (for example equal to it at covariate = 0 and to its reverse
at covariate = 1) is accepted, makes the sampling covariance lose rank exactly, and
is not named in a note; it shows only as a value near 1e-16 of
{cmd:e(gamma_ratio)}. Gamma is never inverted, so standard errors and the test
statistic are still produced.
Whether the standard errors and the test statistic are accurate when Gamma has lost
rank in this way has not been established.


{marker fitindices}{...}
{title:Fit indices}

{p 4 4 2}
RMSEA, CFI and TLI are computed from T3 and its degrees of freedom, as
lavaan's {cmd:rmsea.scaled}, {cmd:cfi.scaled} and {cmd:tli.scaled}. The
baseline model for CFI and TLI has every threshold free in every group and
every polychoric correlation fixed at 0; it is computed from the sample
statistics, with no additional fit, and is stored in {cmd:e(T3_b)} and
{cmd:e(df_T3_b)}. With groups, RMSEA includes a factor of sqrt(G), the
number of groups, as in lavaan, and its sample-size divisor is N - G.

{p 4 4 2}
When T3 is at or below its degrees of freedom, RMSEA is 0 and CFI is 1; TLI
is not truncated and can exceed 1. When the model has 0 degrees of freedom,
all three are missing. When the baseline model's T3 is at or below its own
degrees of freedom, CFI and TLI are missing, since there is no excess misfit
to compare against; RMSEA does not depend on the baseline and is still
reported.

{p 4 4 2}
The baseline is the same for every invariance level fitted to the same
data, so CFI and TLI at different levels are relative to one common
reference. It is not, however, a restricted version of every level: a level
that equates thresholds across groups can fit the thresholds less well than
the baseline does.

{p 4 4 2}
Differences in CFI or RMSEA between invariance levels (ΔCFI, ΔRMSEA) are
descriptive, not a test. Each level's T3 carries its own scaling, so the
difference between two levels' T3 values is not a chi-squared difference
test, and ΔCFI can be positive when restrictions are added. To test one level
against another, use {helpb difftest}. No published ΔCFI or ΔRMSEA cutoff has
been validated for DWLS estimation with the T3 statistic; cutoffs developed for maximum
likelihood should not be assumed to carry over.


{marker standardized}{...}
{title:Standardized estimates}

{p 4 4 2}
Every latent response has unit variance in this model, so standardizing
rescales by the factor standard deviations only, sqrt(Phi_FF), where
Phi = (I-B)^-1 Psi (I-B)^-T is the factor covariance matrix. A loading
becomes lambda*sqrt(Phi_FF), a path beta*sqrt(Phi_FF/Phi_TT), and
thresholds are unchanged. The table shown by {opt standardized} has the same
sections, rows and labels as the coefficient table (every loading, structural
path, variance and covariance, with its standard error), plus:

{p2colset 8 12 12 2}{...}
{p2col:o}{cmd:/var(e.}{it:y}{cmd:)}, the residual variance of each
indicator's latent response, 1 - (Lambda Phi Lambda')_jj, already
standardized;{p_end}
{p2col:o}for an endogenous factor, {cmd:/var(e.}{it:F}{cmd:)} is 1 - R^2,
the share of its variance not explained by its predictors;{p_end}
{p2col:o}for two exogenous factors, {cmd:/cov(}{it:Fa}{cmd:,}{it:Fb}{cmd:)}
is their correlation. The correlation of a pair involving an endogenous
factor is not shown; the path is.{p_end}
{p2colreset}{...}

{p 4 4 2}
An exogenous factor's {cmd:/var(}{it:F}{cmd:)} is 1, and a covariance fixed
at 0 stays 0. Rows the standardization makes constant are not shown; a line above
the table says so. Standard errors are the
delta method's, with the analytic Jacobian of the map. This is lavaan's
{cmd:std.all}; on the example data the standardized values and standard errors match
lavaan 0.7.2's, both metrics, structural and correlated models.

{p 4 4 2}
{bf:With groups}, each group is standardized on its own scale, and its scaling
factor omega enters: a loading becomes omega*lambda*sqrt(Phi_FF), a threshold
omega*tau, an intercept omega*nu, a factor mean alpha/sqrt(Phi_FF), and
{cmd:/var(e.}{it:y}{cmd:)} is 1 - omega^2*(Lambda Phi Lambda')_jj. The
{cmd:/scale()} row is the constant 1. So equal raw parameters can have different
standardized values across groups, and the table says so wherever parameters are
tied. The constant rows are not shown (an exogenous factor's {cmd:/var(}{it:F}{cmd:)},
the {cmd:/scale()} rows, and the intercepts and means the level fixes at 0; the
title line of the table says so); each has a standard error of exactly 0. These match lavaan 0.7.2's {cmd:std.all}; Mplus's
STDYX with groups has not been compared. If a factor variance is at or below 0 in
a group, the solution is not defined and {opt standardized} stops with an error.

{p 4 4 2}
{helpb nlcom} on the posted {cmd:e(b)} and {cmd:e(V)} gives the same values
for a single quantity you write yourself. For a marker-identified exogenous
factor:

{p 8 8 2}
{cmd:nlcom _b[y2:F1]*sqrt(_b[/var(F1)])}

{p 4 4 2}
For an endogenous factor the total variance must be written out. For
{cmd:F2 <- F1}, Phi_F2F2 = beta^2 psi_F1F1 + psi_F2F2:

{p 8 8 2}
{cmd:nlcom _b[y5:F2]*sqrt(_b[F2:F1]^2*_b[/var(F1)] + _b[/var(e.F2)])}

{p 4 4 2}
Rescaling by psi, or skipping the rescaling under {cmd:metric(unitvariance)},
is wrong for every endogenous factor (see {help dwlsmv##metric:The metric}).
These two expressions were checked against lavaan's standardized solution on
the example data; any other expression is yours to write.


{marker comparison}{...}
{title:Comparison with lavaan and Mplus}

{p 4 4 2}
{bf:Pairs at the bound.} lavaan 0.7.2 and Mplus 8.10 fit data that contain
polychoric pairs at the bound, with a warning, where {cmd:dwlsmv} aborts by
default. lavaan 0.7.2 caps at 0.999 and includes the pair with a
"(nearly) 1.0" warning; Mplus 8.10 includes it with a
"zero cells ... can be used to create one new variable" warning (whether Mplus
caps is unknown, as it prints three decimals). {cmd:dwlsmv} aborts because
this version establishes neither an interior estimator for such pairs nor that
truncation at the bound preserves the DWLS inference.

{p 4 4 2}
{bf:Empty-cell 2x2 tables.} For a binary pair with one empty cell the
likelihood has no interior maximum. {cmd:dwlsmv} repairs the table as lavaan
0.7.2 does by default ({cmd:zero.add = c(0.5, 0)},
{cmd:zero.keep.margins = TRUE}): 0.5 is added to the empty cell and its
diagonal partner and taken from the other two, so the margins and thresholds
are unchanged. The correlation is estimated from the repaired table; the
weight matrix is built from scores on the observed data, as in lavaan and in
Mplus 8.10, whose results agree with the same rule (on one such table all
three give 0.9656). The pair enters the fit and is listed in the notes block.
A pre-registered simulation could not show that this choice of weight matrix
is better or worse than the alternative built from the repaired table; it was
adopted because both reference programs use it. Its effect on standard errors
and T3 at small N is not established. A binary pair with two empty cells on a
diagonal (perfect association) is not repaired, and {cmd:dwlsmv} aborts; so does a repaired pair whose correlation reaches the
bound (0.999). On the one dataset compared directly (N = 40, three binary
variables), every standard error from Mplus was lavaan's times sqrt((N-1)/N),
repaired pair or not: a divisor convention, not a consequence of repair.

{p 4 4 2}
{bf:keepbound} admits only the remaining status-2 pairs, at the bound, as
lavaan does. It does not admit 2x2 tables with two empty cells on a
diagonal or other failure statuses, and its effect on standard errors and T3 is not established.

{p 4 4 2}
{bf:Metric.} {cmd:metric(marker)}, the first-listed indicator, is the default
of all three programs.

{p 4 4 2}
{bf:Standardized solutions.} Under this model's parameterization lavaan's
{cmd:std.all} equals its {cmd:std.lv} (lavaan 0.7.2, the example data, both metrics:
estimates identically, standard errors to 2.5e-12).


{marker examples}{...}
{title:Examples}

{p 4 4 2}
The example data, {cmd:dwlsmv_example.dta}, {cmd:dwlsmv_example_mg.dta} and
{cmd:dwlsmv_example_so.dta}, are ancillary files of the package: {cmd:net install}
does not copy them. {cmd:net get dwlsmv}, with the same {cmd:from()} location used
to install, copies them into the current working directory (see {helpb net});
change to the directory where you want the files before running it.

{pstd}Setup{p_end}
{phang2}{cmd:. net get dwlsmv, from(}{it:location}{cmd:)}{p_end}
{phang2}{cmd:. use dwlsmv_example.dta, clear}{p_end}

{pstd}Two factors, F2 predicted by F1, marker metric{p_end}
{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1)}{p_end}

{pstd}Standardized loadings by nlcom{p_end}
{phang2}{cmd:. nlcom _b[y2:F1]*sqrt(_b[/var(F1)])}{p_end}
{phang2}{cmd:. nlcom _b[y5:F2]*sqrt(_b[F2:F1]^2*_b[/var(F1)] + _b[/var(e.F2)])}{p_end}

{pstd}Unit variance: var(F1) = 1 and var(e.F2) = 1{p_end}
{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F2 <- F1), metric(unitvariance)}{p_end}

{pstd}Stored results and the bound marker{p_end}
{phang2}{cmd:. ereturn list}{p_end}
{phang2}{cmd:. matrix list e(bound)}{p_end}

{pstd}Random intercept: aborts (no scale for F1, F2), then the remedy{p_end}
{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1), cov(G*F1@0) cov(G*F2@0)}{p_end}
{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (G -> y1@1 y2@1 y3@1 y4@1 y5@1 y6@1), cov(G*F1@0) cov(G*F2@0) metric(unitvariance)}{p_end}

{pstd}Standardized solution of the last fit, with residual variances{p_end}
{phang2}{cmd:. dwlsmv, standardized}{p_end}

{pstd}Correlated factors, standardized at the fit{p_end}
{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), standardized}{p_end}

{pstd}Measurement invariance: the second example dataset, also copied by
{cmd:net get dwlsmv}{p_end}
{phang2}{cmd:. use dwlsmv_example_mg.dta, clear}{p_end}

{pstd}Configural and scalar invariance across the two groups of {cmd:grp}{p_end}
{phang2}{cmd:. dwlsmv (F -> y1-y5), group(grp) invariance(configural) metric(unitvariance)}{p_end}
{phang2}{cmd:. estimates store configural}{p_end}
{phang2}{cmd:. dwlsmv (F -> y1-y5), group(grp) invariance(scalar) metric(unitvariance)}{p_end}
{phang2}{cmd:. estimates store scalar}{p_end}

{pstd}A table of both fits with {cmd:esttab}, if installed
({cmd:ssc install estout}): the rows fixed in both fits dropped, the fit
statistics added. {cmd:esttab} writes the first group's tag as {cmd:1.grp},
not {cmd:1bn.grp}, and {cmd:drop()} needs that form; {cmd:/:} names the
parameters outside the equations.{p_end}
{phang2}{cmd:. esttab configural scalar, se nostar varwidth(24) mtitles(configural scalar) drop(y*:1.grp y*:2.grp /:var(F)#1.grp /:scale(*)#1.grp /:mean(F)#1.grp) stats(N T3 df_T3 p_T3)}{p_end}

{pstd}A binary covariate, on the same data: the factor regressed on an indicator
of the second group, coded 0 and 1. The thresholds are those of the group coded
0, and the slope is the difference in the mean of the factor between the groups. A covariate
cannot be combined with {opt group()}.{p_end}
{phang2}{cmd:. generate byte grp2 = (grp == 2)}{p_end}
{phang2}{cmd:. dwlsmv (F -> y1-y5) (F <- grp2)}{p_end}

{pstd}A direct effect of the covariate on one item, {cmd:y5}: a shift of its
thresholds between the groups beyond what the factor explains{p_end}
{phang2}{cmd:. dwlsmv (F -> y1-y5) (F <- grp2) (y5 <- grp2)}{p_end}

{pstd}A second-order factor: the third example dataset, also copied by
{cmd:net get dwlsmv}, holds twelve ordinal indicators with four categories each
({cmd:y1}-{cmd:y12}: four factors of three indicators) and a thirteenth
({cmd:y13}), 1,000 observations. One second-order factor {cmd:G} stands over
the four. {cmd:G} is scaled by its variance, and a note says so{p_end}
{phang2}{cmd:. use dwlsmv_example_so.dta, clear}{p_end}
{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6) (F3 -> y7 y8 y9) (F4 -> y10 y11 y12) (G -> F1 F2 F3 F4)}{p_end}

{pstd}The standardized solution of that fit: the paths from {cmd:G}
to the four factors are the standardized second-order loadings{p_end}
{phang2}{cmd:. dwlsmv, standardized}{p_end}


{marker results}{...}
{title:Stored results}

{p 4 4 2}
{cmd:dwlsmv} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(rank)}}number of free parameters{p_end}
{synopt:{cmd:e(T3)}}scaled-and-shifted test statistic; not posted, nor {cmd:e(p_T3)}, the fit
indices, {cmd:e(T3_b)}, {cmd:e(df_T3_b)} or {cmd:e(V)}, for a fit that did not converge{p_end}
{synopt:{cmd:e(df_T3)}}degrees of freedom of T3{p_end}
{synopt:{cmd:e(p_T3)}}p-value of T3 (missing if df is 0){p_end}
{synopt:{cmd:e(rmsea_s)}}RMSEA from T3 (missing if df is 0){p_end}
{synopt:{cmd:e(cfi_s)}}CFI from T3 and the baseline (missing if df is 0 or
the baseline's T3 does not exceed its df){p_end}
{synopt:{cmd:e(tli_s)}}TLI, likewise; not truncated at 1{p_end}
{synopt:{cmd:e(T3_b)}}T3 of the baseline model (thresholds free, correlations
fixed at 0){p_end}
{synopt:{cmd:e(df_T3_b)}}degrees of freedom of the baseline model{p_end}
{synopt:{cmd:e(converged)}}1 if the reported fit converged, 0 otherwise{p_end}
{synopt:{cmd:e(iterations)}}number of iterations of the reported fit{p_end}
{synopt:{cmd:e(gamma_ratio)}}smallest/largest singular value of the sampling
covariance of the thresholds and correlations, scale-free, within its rank{p_end}
{synopt:{cmd:e(gamma_ratio_raw)}}the same ratio, raw{p_end}
{synopt:{cmd:e(gamma_rank)}}the rank it is computed within{p_end}
{synopt:{cmd:e(bread_ratio)}}smallest/largest singular value of the weighted
derivative matrix{p_end}
{p2col 5 20 24 2:}diagnostic values with no established cutoff; missing when
not evaluated{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:dwlsmv}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(version)}}version of {cmd:dwlsmv} that produced the results{p_end}
{synopt:{cmd:e(metric)}}{cmd:marker} or {cmd:unitvariance}{p_end}
{synopt:{cmd:e(recovery)}}{cmd:none}, {cmd:recovered}, {cmd:ineligible} or
{cmd:failed}; see {help dwlsmv##recovery:recovery}{p_end}
{synopt:{cmd:e(indicators)}}indicators, in {cmd:sem}'s {cmd:e(oyvars)} order{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}; {cmd:b} for a fit that did not converge;
with {cmd:efa()}, {cmd:b V eigen} or {cmd:b eigen}{p_end}
{synopt:{cmd:e(canon)}}the model that was fitted, written out; {helpb difftest} refits from it{p_end}
{synopt:{cmd:e(keepbound)}}{cmd:keepbound} if specified{p_end}
{synopt:{cmd:e(datasignature)}}the signature of the indicators in the estimation sample, which {helpb difftest} checks{p_end}
{synopt:{cmd:e(datasignaturevars)}}the variables the signature covers{p_end}
{synopt:{cmd:e(covariates)}}the binary covariate, for a model conditional on it; such a fit
posts no {cmd:e(rmsea_s)}, {cmd:e(cfi_s)}, {cmd:e(tli_s)}, {cmd:e(T3_b)}, {cmd:e(df_T3_b)}
or {cmd:e(bound)}; it does post {cmd:e(gamma_ratio)}, {cmd:e(gamma_ratio_raw)},
{cmd:e(gamma_rank)} and {cmd:e(bread_ratio)}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}every parameter, free and fixed, in the order thresholds,
loadings, regressions, factor (co)variances; a fixed parameter carries its
value{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimates, exactly zero
in fixed rows and columns; not posted for a fit that did not converge{p_end}
{synopt:{cmd:e(Cns)}}one row per fixed parameter, selecting that column at its
value (posted only when some parameter is fixed){p_end}
{synopt:{cmd:e(Cns_b)}}the same rows, posted instead of {cmd:e(Cns)} for a fit
that did not converge (Stata keeps {cmd:e(Cns)} for results with {cmd:e(V)}){p_end}
{synopt:{cmd:e(bound)}}1 x P, one column per polychoric pair; 1 for a pair
included at the bound under {opt keepbound}, else 0{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{p 4 4 2}
{bf:Names in e(b).} A loading is {it:y}{cmd::}{it:F} ({cmd:_b[y2:F1]}); a
regression is {it:To}{cmd::}{it:From} ({cmd:_b[F2:F1]}); a variance is
{cmd:/var(}{it:F}{cmd:)} for an exogenous factor and
{cmd:/var(e.}{it:F}{cmd:)} for an endogenous one, fixed or free; a covariance
is {cmd:/cov(}{it:Fa}{cmd:,}{it:Fb}{cmd:)} with {cmd:e.} on each endogenous
factor; thresholds are {cmd:/}{it:y}{cmd::cut1}, ..., {cmd:/}{it:y}{cmd::cut}{it:(K-1)}
({cmd:_b[/y1:cut1]}). Structural names are those {cmd:sem} gives the same
model. Residual variances of the indicators are derived and not stored, and
neither is the standardized solution.

{p 4 4 2}
{bf:e(bound)'s stripe} is {it:varA}{cmd::}{it:varB}, so
{cmd:colnumb(e(bound), "q7_anx:q9_dep")} finds a pair. (A stripe element is
one name of at most 32 characters, so the {cmd:~~} form shown in the output
cannot be used there.)

{p 4 4 2}
{bf:With {cmd:efa()}.} {cmd:e(b)} holds the rotated loadings, the factor
correlations and the thresholds. {cmd:e(V)} is the covariance of these rotated
estimates, over every column of {cmd:e(b)}; it is posted except when
{cmd:e(rot_nose)} is set, and then {cmd:e(properties)} is {cmd:b eigen} rather
than {cmd:b V eigen}. {cmd:e(rank)} is the number of free parameters of the
unrotated model that was fitted, which is also the rank of {cmd:e(V)}. In addition:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(rot_starts)}}rotation starts tried{p_end}
{synopt:{cmd:e(rot_converged)}}rotation starts that converged{p_end}
{synopt:{cmd:e(rot_minima)}}distinct minima among the converged starts{p_end}
{synopt:{cmd:e(rot_atbest)}}converged starts that reached the reported solution{p_end}
{synopt:{cmd:e(rot_Q)}}geomin criterion at the reported solution{p_end}
{synopt:{cmd:e(rot_gpnorm)}}its gradient norm (converged below 1e-5){p_end}
{synopt:{cmd:e(rot_eps)}}geomin epsilon used{p_end}
{synopt:{cmd:e(rot_gap)}}criterion gap to the next distinct minimum; missing if
there is only one{p_end}
{synopt:{cmd:e(rot_cond)}}condition number of the matrix used to compute the rotated
covariance; stored for inspection, with no cutoff; it does not indicate how well
the intervals cover{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(Lambda)}}rotated loadings, indicators by factors{p_end}
{synopt:{cmd:e(Phi)}}factor correlations{p_end}
{synopt:{cmd:e(Ev)}}eigenvalues of the indicators' polychoric correlation matrix,
largest first, 1 x p; used by {helpb screeplot}{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(efa)}}the factors named in {opt efa()}{p_end}
{synopt:{cmd:e(rotation)}}{cmd:geomin oblique}{p_end}
{synopt:{cmd:e(rot_nose)}}set only when no covariance is stored: {cmd:unconverged},
{cmd:inadmissible}, {cmd:rank}, {cmd:critical} or {cmd:acov}, naming why{p_end}
{p2colreset}{...}

{p 4 4 2}
{bf:e(V) after efa().} The rotated covariance is singular by construction: it has
m(m-1) more rows than its rank, {cmd:e(rank)}, and it is not regularized.
{cmd:test} and {cmd:lincom} work as usual for ordinary combinations of the rotated
loadings and factor correlations. A combination that moves the rotated solution
along its own family of equivalent rotations (it mixes loadings and factor
correlations in a particular way) has variance zero; tests and intervals for such a
combination have no meaning, and Stata's handling of them depends on the sign of a
rounding residue, so it differs between computers: the constraint is dropped on
some and reported as decisive on others. Do not interpret them. The standard errors
in {cmd:e(V)} carry the caveat in
{help dwlsmv##efase:Remarks on rotated standard errors}.

{p 4 4 2}
{bf:With {cmd:invariance()}.} {cmd:e(b)} holds every group's parameters at
every level, in the order and with the names described in
{help dwlsmv##output:Output}, and adds each indicator's latent-response
intercept ({it:y}{cmd::}{it:tag}) and scale ({cmd:/scale(}{it:y}{cmd:)#}{it:tag})
and each factor's mean ({cmd:/mean(}{it:F}{cmd:)#}{it:tag}). Parameters the
level does not estimate, the first group's included, are posted at their
fixed values; {cmd:e(rank)} counts free parameters only. A parameter
constrained equal across groups appears once per group. {cmd:e(V)} divides by
N - G, G the number of groups, as lavaan does. {cmd:e(bread_ratio)} is
computed on all groups together, each weighted by its share of the
observations. {cmd:e(bound)}, {cmd:e(gamma_ratio)}, {cmd:e(gamma_ratio_raw)}
and {cmd:e(gamma_rank)} are not stored with groups in this version, and
{cmd:e(recovery)} is {cmd:none}. Also stored:

{synoptset 24 tabbed}{...}
{synopt:{cmd:e(invariance)}}the level: {cmd:configural}, {cmd:thresholds},
{cmd:metric} or {cmd:scalar}{p_end}
{synopt:{cmd:e(K)}}G x p matrix of category counts, rows named by group
value, columns by indicator{p_end}
{synopt:{cmd:e(canon)}}the model that was fitted, written out with the
level's equality constraints; {helpb difftest} refits from it{p_end}
{synopt:{cmd:e(canon_configural)}}the configural model the level was
generated from; equal to {cmd:e(canon)} for
{cmd:invariance(configural)}{p_end}
{synopt:{cmd:e(groupvar)}}the group variable{p_end}
{synopt:{cmd:e(keepbound)}}{cmd:keepbound} if specified{p_end}
{synopt:{cmd:e(datasignature)}}the signature of the indicators and the group
variable in the estimation sample, which {helpb difftest} checks{p_end}
{synopt:{cmd:e(datasignaturevars)}}the indicators and the group variable{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{p 4 8 2}
Asparouhov, T., and B. Muthén. 2006. Robust chi square difference testing with
mean and variance adjusted test statistics. Mplus Web Notes: No. 10.
{p_end}

{p 4 8 2}
Asparouhov, T., and B. Muthén. 2010. Simple second order chi-square correction.
Mplus technical note, May 3, 2010.
{browse "https://www.statmodel.com/download/WLSMV_new_chi21.pdf"}.
{p_end}

{p 4 8 2}
Muthén, B. 1984. A general structural equation model with dichotomous, ordered
categorical, and continuous latent variable indicators. {it:Psychometrika} 49: 115{c -}132.
{p_end}

{p 4 8 2}
Muthén, B., S. H. C. du Toit, and D. Spisic. 1997. Robust inference using
weighted least squares and quadratic estimating equations in latent variable
modeling with categorical and continuous outcomes. Unpublished manuscript.
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
Help:  {manhelp sem SEM}, {manhelp gsem SEM}, {manhelp nlcom R}, {helpb difftest}, {helpb efacompare}
{p_end}
