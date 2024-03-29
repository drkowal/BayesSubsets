% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_subsel.R
\name{accept_family_randint}
\alias{accept_family_randint}
\title{Compute the acceptable family of linear subsets for the
random intercept model}
\usage{
accept_family_randint(
  post_y_pred,
  post_lpd,
  post_sigma_e,
  post_sigma_u,
  XX,
  YY,
  indicators,
  post_y_pred_sum = NULL,
  eps_level = 0.05,
  eta_level = 0,
  K = 10,
  sir_frac = 0.5,
  plot = TRUE
)
}
\arguments{
\item{post_y_pred}{\code{S x m x n} matrix of posterior predictive draws
at the given \code{XX} covariate values for \code{m} replicates per subject}

\item{post_lpd}{\code{S} evaluations of the log-likelihood computed
at each posterior draw of the parameters}

\item{post_sigma_e}{(\code{nsave}) draws from the posterior distribution
of the observation error SD}

\item{post_sigma_u}{(\code{nsave}) draws from the posterior distribution
of the random intercept SD}

\item{XX}{\code{n x p} matrix of covariates at which to evaluate}

\item{YY}{\code{m x n} matrix of response variables (optional)}

\item{indicators}{\code{L x p} matrix of inclusion indicators (booleans)
where each row denotes a candidate subset}

\item{post_y_pred_sum}{(\code{nsave x n}) matrix of the posterior predictive
draws summed over the replicates within each subject (optional)}

\item{eps_level}{probability required to match the predictive
performance of the "best" model (up to \code{eta_level})}

\item{eta_level}{allowable margin (%) between each acceptable model
and the "best" model}

\item{K}{number of cross-validation folds (optional)}

\item{sir_frac}{fraction of the posterior samples to use for SIR (optional)}

\item{plot}{logical; if TRUE, include a plot to summarize the predictive
performance across candidate subsets}
}
\value{
a list containing the following elements:
\itemize{
\item \code{all_accept}: indices (i.e., rows of \code{indicators})
that correspond to the acceptable subsets
\item \code{beta_hat_small} linear coefficients for the
smallest acceptable model
\item \code{beta_hat_min} linear coefficients for the
"best" acceptable model
\item \code{ell_small}: index (i.e., row of \code{indicators}) of the
smallest acceptable model
\item \code{ell_min}: index (i.e., row of \code{indicators}) of the
"best" acceptable model
}
}
\description{
Given output from a Bayesian random intercept model and a candidate of
subsets, compute the *acceptable family* of subsets that
match or nearly match the predictive accuracy of the "best" subset.
The acceptable family may be computed for any set of covariate values
\code{XX}; if \code{XX = X} are the in-sample points, then
cross-validation is used to assess out-of-sample predictive performance.
}
\details{
When \code{XX = X} is the observed covariate values,
then \code{post_lpd} and \code{yy} must be provided. These
are used to compute the cross-validated predictive and empirical
squared errors; the predictive version relies on a sampling importance-resampling
procedure.

When \code{XX} corresponds to a new set of covariate values, then set \code{post_lpd = NULL}
and \code{yy = NULL} (these are the default values).

Additional details on the predictive and empirical comparisons are
in \code{pp_loss_randint}.
}
