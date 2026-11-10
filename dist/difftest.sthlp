{smcl}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col:{cmd:difftest} {hline 2}}Difference test between nested {cmd:dwlsmv} models{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:difftest} {it:modelspec1} [{it:modelspec2} ...]

{p 4 4 2}
where {it:modelspec} is the name of results saved by {cmd:estimates store}, or
{cmd:.} for the active estimation results. With one name, that model is
compared with the active results.


{title:Description}

{p 4 4 2}
{cmd:difftest} tests a model against a less restricted model that contains it,
both fitted by {cmd:dwlsmv}, with the scaled-and-shifted difference test of
Asparouhov and Muthén (2006). With maximum likelihood, nested models are
often compared by a likelihood-ratio test, which is the difference between
their chi-squared statistics. That does not work for T3: each model's T3
carries its own scaling, so the difference between two T3 values does not
follow a chi-squared distribution.

{p 4 4 2}
Like {helpb lrtest}, it takes the names of stored results; the order in which
they are named does not matter. Unlike {cmd:lrtest}, which compares exactly
two models, {cmd:difftest} accepts three or more: they must form a chain, each
nested in the next, and each adjacent pair is tested.

{p 4 4 2}
{cmd:difftest} refits every model on the same data and checks that each refit
reproduces the stored results before computing anything. It compares
single-group models, and multiple-group models fitted with {cmd:group()} and
{cmd:invariance()} to the same data and group variable, for example the
configural, thresholds, metric and scalar levels as a chain.

{p 4 4 2}
Fit the models to be compared, and run {cmd:difftest}, in the same version
and edition of Stata (for example, Stata/MP 18). Results can differ slightly
across Stata versions and editions. For single-group models, differences in
the number of processors are allowed for. For multiple-group models, also keep
the number of processors the same ({help set processors}): that is the setting
under which {cmd:difftest}'s check of the refit has been validated for them.


{title:What is checked}

{p 4 4 2}
{cmd:difftest} refuses the comparison, with a message naming the check and
the models, in any of these cases:

{phang2}o  a model is an EFA fit ({cmd:efa()});{p_end}
{phang2}o  a model was fitted by another version of {cmd:dwlsmv};{p_end}
{phang2}o  the models were fitted to different data or different estimation samples;{p_end}
{phang2}o  the models have different indicators, or the same indicators in a different order;{p_end}
{phang2}o  {cmd:keepbound} was specified for only one of them;{p_end}
{phang2}o  the models use different metrics;{p_end}
{phang2}o  one model is a multiple-group fit and the other is not, or they use different group variables;{p_end}
{phang2}o  a refit does not reproduce the stored fit, or did not converge;{p_end}
{phang2}o  the models are not nested;{p_end}
{phang2}o  the fixed values in a model force a boundary, for example a factor correlation fixed at 1;{p_end}
{phang2}o  the more restricted model's unscaled statistic, T = (N - G)F with G groups, is below the less restricted model's by more than the precision of the fits, (N - G) times the convergence tolerance; a nested model cannot fit better than the model that contains it.{p_end}

{p 4 4 2}
Models nested only through a reparameterization (a different marker, for
example) are refused.

{p 4 4 2}
For multiple-group models, the comparison is available for levels of
{cmd:invariance()} fitted from the same configural specification. It is not
available for a model with a structural path, for other models with equality
constraints, or when a loading or factor covariance is fixed at a nonzero
value, as in {cmd:cov(F1*F2@.4)}. Nesting has not been established for that
last case, so a valid comparison may also be refused. A covariance fixed at 0
is allowed.


{title:Output}

{p 4 4 2}
A table of the models, most restricted first, with their T3, degrees of
freedom and scaled RMSEA, CFI and TLI; then a table of the tests, one row per
adjacent pair, with the difference in degrees of freedom d, the statistic,
its p-value, and the differences in CFI and RMSEA. Each difference is the
more restricted model's value minus the less restricted model's (H0 minus
H1), so a restriction that worsens fit shows a negative dCFI and a positive
dRMSEA; the sign shown is the actual one, and either can occur. A step
skipped because d = 0 shows no difference, and the next step is compared
with the last model that was fitted. A restriction that does not change the
fit to the precision the fits are computed to is shown as "within eps; not
tested". Differences in CFI or RMSEA are descriptive, not a
test; see {help dwlsmv##fitindices:Fit indices}.


{title:Stored results}

{p 4 4 2}
{cmd:difftest} stores the following in {cmd:r()}:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}

{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(models)}}one row per model, most restricted first: npar, T3,
df, rmsea, cfi, tli{p_end}
{synopt:{cmd:r(steps)}}one row per tested pair (H0:H1): d, T_d3, p, dCFI,
dRMSEA (each H0 minus H1), status (1 tested, 2 within eps, 0 skipped){p_end}

{p 4 4 2}
The models are ordered by their number of free parameters, not in the order
they were named, so use the row names to tell the rows apart: each row of
{cmd:r(models)} is named by the model's stored-results name, and each row of
{cmd:r(steps)} by its pair, as {it:H0}:{it:H1}. The row name used for the
active results ({cmd:.}) is not fixed and may change.


{title:Example}

{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), cov(F1*F2@0) metric(unitvariance)}{p_end}
{phang2}{cmd:. estimates store indep}{p_end}
{phang2}{cmd:. dwlsmv (F1 -> y1 y2 y3) (F2 -> y4 y5 y6), metric(unitvariance)}{p_end}
{phang2}{cmd:. difftest indep}{p_end}


{title:Reference}

{p 4 8 2}
Asparouhov, T., and B. Muthén. 2006. Robust chi square difference testing with
mean and variance adjusted test statistics. Mplus Web Notes: No. 10.
{p_end}
