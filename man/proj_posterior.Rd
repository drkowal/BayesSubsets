% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_subsel.R
\name{proj_posterior}
\alias{proj_posterior}
\title{Projected predictive distribution for regression coefficients}
\usage{
proj_posterior(post_y_pred, XX, sub_x = 1:ncol(XX), use_ols = TRUE)
}
\arguments{
\item{post_y_pred}{\code{S x n} matrix of posterior predictive draws
at the given \code{XX} covariate values}

\item{XX}{\code{n x p} matrix of covariates}

\item{sub_x}{vector of inclusion indicators for the \code{p} covariates;
the remaining coefficients will be fixed at zero}

\item{use_ols}{logical; if TRUE, use ordinary least squares regression (default);
otherwise use logistic regression}
}
\value{
\code{post_beta}: the \code{S x p} matrix of
draws from the projected predictive distribution of  the regression coefficients.
}
\description{
Given draws from the predictive distribution, project these
draws onto (a  subset of) the covariates. This produces
many predictive draws for the regression coefficients,
which provides uncertainty quantification.
}
