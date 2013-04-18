\name{one.sided}
\alias{alignment}
\alias{alignment.att}
\alias{one.sided}
\alias{one.sided.att}
\title{Confounding functions}
\usage{
  one.sided(alpha, pscores, treat)

  alignment(alpha, pscores, treat)

  one.sided.att(alpha, pscores, treat)

  alignment.att(alpha, pscores, treat)
}
\arguments{
  \item{alpha}{vector of confounding values to use in the
  sensitivity analysis.}

  \item{pscores}{vector of propensity scores for each
  unit.}

  \item{treat}{vector of treatment values for each unit.}
}
\description{
  Various confounding functions for use with
  \code{\link{causalsens}}.
}

