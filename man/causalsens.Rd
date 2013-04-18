\name{causalsens}
\alias{causalsens}
\title{Calculate sensitivity of causal estimates to unmeasured confounding.}
\usage{
  causalsens(model.y, model.t, cov.form,
    confound = one.sided, data, alpha)
}
\arguments{
  \item{model.y}{outcome model object. Currently only
  handles \code{lm} objects.}

  \item{model.t}{propensity score model. Currently assumes
  a \code{glm} object.}

  \item{cov.form}{one-sided formula to describe any
  covariates to be included in the parital R^2
  calculations.}

  \item{confound}{function that calculates the confounding
  function. This function must take arguments \code{alpha},
  \code{pscores}, and \code{treat}. Defaults to
  \code{\link{one.sided}}. Other functions included with
  the package are \code{\link{one.sided.att}},
  \code{\link{alignment}}, and
  \code{\link{alignment.att}}.}

  \item{data}{data frame to find the covariates from
  \code{cov.form}.}

  \item{alpha}{vector of confounding values to pass the
  confounding function. Defaults to 11 points from -0.5 to
  0.5 for binary outcome variable, and 11 points covering
  the inter-quartile range for non-binary outcome
  variables.}
}
\value{
  Returns an object of class \code{causalsens}. \itemize{
  \item \code{sens} data frame containing alpha values,
  partial R^2s, estimates, and 95% confidence intervals
  \item \code{partial.r2} vector of partial R^2 values for
  the covariates to compare to sensitivity analysis
  results. }
}
\description{
  This function performs a sensitivity analysis of causal
  effects different assumptions about unmeasured
  confounding, as described by Blackwell (2013).
}

