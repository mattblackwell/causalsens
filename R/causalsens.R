#' Experimental data from the job training program first studied by LaLonde (1986)
#'
#' A dataset of units in an experimental evaluation of a job training
#' program. Subset to those units with two years of pre-treatment income data.
#'
#'
#' \itemize{
#'   \item \code{age} - age in years.
#'   \item \code{education} - number of years of schooling.
#'   \item \code{black} - 1 if black, 0 otherwise.
#'   \item \code{hispanic} - 1 if Hispanic, 0 otherwise.
#'   \item \code{married} - 1 if married, 0 otherwise.
#'   \item \code{nodegree} - 1 if no high school degree, 0 otherwise.
#'   \item \code{re74} - earnings ($) in 1974.
#'   \item \code{re75} - earnings ($) in 1975.
#'   \item \code{re78} - earnings ($) in 1978.
#'   \item \code{u74} - 1 if unemployed in 1974, 0 otherwise.
#'   \item \code{u75} - 1 if unemployed in 1975, 0 otherwise.
#'   \item \code{treat} - 1 if treated, 0 otherwise.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name lalonde.exp
#' @usage data(lalonde.exp)
#' @format A data frame with 445 rows and 12 variables
#' @references LaLonde, Robert J. (1986). Evaluating the Econometric
#' Evaluations of Training Programs with Experimental Data. The
#' American Economic Review, 76(4), 604--620.
NULL

#' Non-experimental data from Lalonde (1986)
#'
#' A dataset of experimental treated units and non-experimental
#' control units from the Panel Study of Income Dynamics (PSID).
#'
#' \itemize{
#'   \item \code{age} - age in years.
#'   \item \code{education} - number of years of schooling.
#'   \item \code{black} - 1 if black, 0 otherwise.
#'   \item \code{hispanic} - 1 if Hispanic, 0 otherwise.
#'   \item \code{married} - 1 if married, 0 otherwise.
#'   \item \code{nodegree} - 1 if no high school degree, 0 otherwise.
#'   \item \code{re74} - earnings ($) in 1974.
#'   \item \code{re75} - earnings ($) in 1975.
#'   \item \code{re78} - earnings ($) in 1978.
#'   \item \code{u74} - 1 if unemployed in 1974, 0 otherwise.
#'   \item \code{u75} - 1 if unemployed in 1975, 0 otherwise.
#'   \item \code{treat} - 1 if treated, 0 otherwise.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name lalonde.psid
#' @usage data(lalonde.psid)
#' @format A data frame with 2675 rows and 12 variables
#' @references LaLonde, Robert J. (1986). Evaluating the Econometric
#' Evaluations of Training Programs with Experimental Data. The
#' American Economic Review, 76(4), 604--620.
NULL


#' Calculate sensitivity of causal estimates to unmeasured confounding.
#'
#' This function performs a sensitivity analysis of causal effects
#' different assumptions about unmeasured confounding, as described
#' by Blackwell (2013).
#'
#'
#' @param model.y outcome model object. Currently only handles
#' \code{lm} objects.
#' @param model.t propensity score model. Currently assumes a
#' \code{glm} object.
#' @param cov.form one-sided formula to describe any covariates to be
#' included in the parital R^2 calculations.
#' @param confound function that calculates the confounding
#' function. This function must take arguments \code{alpha},
#' \code{pscores}, and \code{treat}. Defaults to
#' \code{\link{one.sided}}. Other functions included with the package are
#' \code{\link{one.sided.att}}, \code{\link{alignment}}, and \code{\link{alignment.att}}.
#' @param data data frame to find the covariates from \code{cov.form}.
#' @param alpha vector of confounding values to pass the confounding
#' function. Defaults to 11 points from -0.5 to 0.5 for binary
#' outcome variable, and 11 points covering the
#' a interval with width equal to the inter-quartile range and centered at 0 for
#' non-binary outcome variables.
#' @param level level of the confidence interval returned.
#' @return Returns an object of class \code{causalsens}.
#' \itemize{
#'   \item \code{sens} data frame containing alpha values, partial
#' R^2s, estimates, and 95% confidence intervals
#' \item \code{partial.r2} vector of partial R^2 values for the
#' covariates to compare to sensitivity analysis results.
#' }
#'
#' @examples
#' data(lalonde.exp)
#'
#' ymodel <- lm(re78 ~ treat+age + education + black + hispanic +
#' married + nodegree + re74 + re75 + u74 + u75, data = lalonde.exp)
#' 
#' pmodel <- glm(treat ~ age + education + black + hispanic + married
#' + nodegree + re74 + re75 + u74 + u75, data = lalonde.exp,
#' family = binomial())
#'
#' alpha <- seq(-4500, 4500, by = 250)
#'
#' ll.sens <- causalsens(ymodel, pmodel, ~ age + education, data =
#' lalonde.exp, alpha = alpha, confound = one.sided.att)
#' 
#' @export
#' @importFrom stats coef fitted lm.fit lm.wfit model.frame
#' model.matrix model.response model.weights qt quantile terms var
causalsens <- function(model.y, model.t, cov.form, confound = one.sided,
                       data, alpha, level = 0.95) {

  if (inherits(model.y, "glm")) {
    stop("Only works for linear outcome models right now. Check back soon.")
  }
  y.dat <- model.frame(model.y, data = data)
  t.dat <- model.frame(model.t, data = data)
  t.name <- colnames(t.dat)[1]
  treat <- model.response(t.dat)
  y <- model.response(y.dat)
  w <- model.weights(y.dat)
  X <- model.matrix(model.y, y.dat)
  pscores <- fitted(model.t)

  mt_c <- terms(cov.form, data = data)
  attr(mt_c, "intercept") <- 0
  overlap <- attr(mt_c, "term.labels") %in% colnames(X)
  Xc <- model.matrix(mt_c, data)[, !overlap]
  rn.y <- row.names(X)
  rn.t <- names(pscores)
  K <- ncol(X)
  N <- nrow(X)

  if (!identical(rn.y, rn.t)) {
    bothrows <- intersect(rn.y, rn.t)
    X <- X[bothrows, ]
    Xc <- Xc[bothrows, ]
    y <- y[bothrows]
    treat <- treat[bothrows]
    pscores <- pscores[bothrows]
    if (!is.null(w)) w <- w[bothrows]
  }
  Xall <- cbind(X, Xc)
  X0 <- X[treat == 0, -which(colnames(X) == t.name)]

  if (missing(alpha)) {
    if (length(unique(y)) == 2) {
      alpha <- seq(-0.5, 0.5, length = 11)
    } else {
      iqr <- quantile(y, 0.75) - quantile(y, 0.25)
      alpha <- seq(-iqr / 2, iqr / 2, length = 11)
    }
  }

  sens <- matrix(NA, nrow = length(alpha), ncol = 6)
  colnames(sens) <- c("rsqs", "alpha", "estimate", "lower", "upper", "se")
  sens[, "alpha"] <- alpha
  for (j in 1:length(alpha)) {
    adj <- confound(alpha = alpha[j], pscores = pscores, treat = treat)
    y.adj <- y - adj
    if (is.null(w)) {
      s.out <- lm.fit(x = X, y = y.adj, offset = model.y$offset)
      r.out <- lm.fit(x = X0, y = y.adj[treat == 0],
                      offset = model.y$offset[treat == 0])
      s.rss <- sum(s.out$residuals ^ 2)
      r.rss <- sum(r.out$residuals ^ 2)
      vt <- var(treat)
    } else {
      s.out <- lm.wfit(x = X, y = y.adj, w = w, offset = model.y$offset)
      r.out <- lm.wfit(x = X0, y = y.adj[treat == 0], w = w[treat == 0],
                       offset = model.y$offset[treat == 0])
      s.rss <- sum(w[treat == 0] * s.out$residuals ^ 2)
      r.rss <- sum(w[treat == 0] * r.out$residuals ^ 2)
      vt <- sum(w * (treat - mean(treat)) ^ 2) / (N - 1)
    }
    s.sigma2 <- s.rss / s.out$df.residual
    r.sigma2 <- r.rss / r.out$df.residual
    R <- s.sigma2  * chol2inv(s.out$qr$qr[1L:K, 1L:K, drop = FALSE])
    q.alpha <- abs(qt( (1 - level) / 2, df = s.out$df.residual))
    t.pos <- which(names(coef(s.out)) == t.name)
    sens.se <- sqrt(R[t.pos, t.pos])
    sens[j, 3] <- coef(s.out)[t.pos]
    sens[j, 4] <- coef(s.out)[t.pos] - q.alpha * sens.se
    sens[j, 5] <- coef(s.out)[t.pos] + q.alpha * sens.se
    sens[j, 6] <- sens.se
    sens[j, 1] <- alpha[j] ^ 2 * vt / r.sigma2
  }

  ## Calculate all partial R2 for covariates in untreated group.
  Kc <- ncol(Xall)
  prsqs <- rep(NA, Kc - 1)
  if (is.null(w)) {
    p.all <- lm.fit(Xall[treat == 0, ], y[treat == 0],
                    offset = model.y$offset[treat == 0])
    pa.rss <- sum(p.all$residuals ^ 2)
  } else {
    p.all <- lm.fit(Xall[treat == 0, ], y[treat == 0], w = w[treat == 0],
                    offset = model.y$offset[treat == 0])
    pa.rss <- sum(w[treat == 0] * p.all$residuals ^ 2)
  }
  for (k in 2:Kc) {
    jj <- setdiff(1L:Kc, k)
    if (is.null(w)) {
      p.out <- lm.fit(Xall[treat == 0, jj, drop = FALSE], y[treat == 0],
                      offset = model.y$offset[treat == 0])
      p.rss <- sum(p.out$residuals ^ 2)
    } else {
      p.out <- lm.fit(Xall[treat == 0, jj, drop = FALSE], y[treat == 0],
                      w = w[treat == 0], offset = model.y$offset[treat == 0])
      p.rss <- sum(w[treat == 0] * p.out$residuals ^ 2)
    }
    prsqs[k - 1] <- (p.rss - pa.rss) / p.rss
    names(prsqs)[k - 1] <- colnames(X)[k]
  }
  out <- list(sens = data.frame(sens), partial.r2 = prsqs)
  class(out) <- "causalsens"
  return(out)

}

print.causalsens <- function(x, ...) {
  cat("Sensitivity analysis output: \n")
  print(x$sens)
  cat("\nCovariate partial r-squared: \n")
  print(x$partial.r2)
  invisible()
}


summary.causalsens <- function(object, ...) {
  cat("Maximum estimate: ", max(object$sens$estimate), "\n")
  cat("Minimum estimate: ", min(object$sens$estimate), "\n")
  cat("Parial R^2 range: [", min(object$sens$rsqs), ", ",
      max(object$sens$rsqs), "]\n", sep = "")
  cat("Covariate R^2 range: [", min(object$partial.r2), ", ",
      max(object$partial.r2), "]\n", sep = "")
}

#' Plot a causal sensitivity analysis.
#'
#' Plot the results of a sensitivity analysis against unmeasured
#' confounding as perfomed by \code{\link{causalsens}}
#' @export
#' @param x \code{causalsens} object.
#' @param type a string taking either the value \code{"r.squared"}
#' (default), which plots the estimated effects as a function of the
#' partial R-squared values, or \code{"raw"}, which plots them as a
#' function of the raw confounding values, \code{alpha}.
#' @param ... other parameters to pass to the plot.
#' @importFrom grDevices rgb
#' @importFrom graphics abline points polygon
plot.causalsens <- function(x, type = "r.squared", ...) {
  m <- match.call(expand.dots = TRUE)
  m[[1L]] <- quote(graphics::plot)
  if (type == "r.squared") {
    m$x <- sign(x$sens$alpha) * x$sens$rsqs
    if (is.null(m$xlab)) {
      m$xlab <- "Variance explained by confounding"
    }
  } else if (type == "raw") {
    m$x <- x$sens$alpha
    if (is.null(m$xlab)) {
      m$xlab <- "Amount of confounding"
    }
  } else {
    stop("type must be 'r.squared' or 'raw'")
  }
  if (is.null(m$ylim)) {
    m$ylim <- c(min(x$sens$lower), max(x$sens$upper))
  }
  if (is.null(m$ylab)) {
    m$ylab <- "Estimated effect"
  }
  m$y <- x$sens$estimate
  m$type <- "l"
  eval(m, parent.frame())
  abline(h = 0, col = "grey")
  abline(v = 0, col = "grey")
  polygon(x = c(m$x, rev(m$x)),
          y = c(x$sens$lower, rev(x$sens$upper)),
          col = rgb(0.5, 0.5, 0.5, alpha = 0.5), border = NA)
  if (type == "r.squared") {
    points(x = x$partial.r2, y = rep(0, length(x$partial.r2)), pch = 4)
    points(x = -x$partial.r2, y = rep(0, length(x$partial.r2)), pch = 4)
  }
  invisible()
}


#' Confounding functions
#'
#' Various confounding functions for use with
#' \code{\link{causalsens}}.
#' @aliases one.sided alignment one.sided.att alignment.att
#' @param alpha vector of confounding values to use in the
#' sensitivity analysis.
#' @param pscores vector of propensity scores for each unit.
#' @param treat vector of treatment values for each unit.
#' @export
one.sided <- function(alpha, pscores, treat) {
  adj <- alpha * (1 - pscores) * treat - alpha * pscores * (1 - treat)
  return(adj)
}
#' @rdname one.sided
#' @export
alignment <- function(alpha, pscores, treat) {
  adj <- alpha * (1 - pscores) * treat + alpha * pscores * (1 - treat)
  return(adj)
}
#' @rdname one.sided
#' @export
one.sided.att <- function(alpha, pscores, treat) {
  adj <- -alpha  * (1 - treat)
  return(adj)
}
#' @rdname one.sided
#' @export
alignment.att <- function(alpha, pscores, treat) {
  adj <- alpha  * (1 - treat)
  return(adj)
}
