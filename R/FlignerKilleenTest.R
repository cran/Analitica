#' Fligner-Killeen Test for Homogeneity of Variances (Manual Implementation)
#'
#' Performs a non-parametric Fligner-Killeen test for equality of variances
#' across two or more groups, using raw vectors via a formula interface.
#'
#' This test is particularly useful when the assumption of normality is violated,
#' as it is robust to outliers and distributional deviations. It serves as a reliable
#' alternative to Bartlett’s test when data do not follow a normal distribution.
#'
#' Advantages:
#' - Non-parametric: No assumption of normality.
#' - Robust to outliers.
#' - Suitable for heterogeneous sample sizes.
#'
#' Disadvantages:
#' - Less powerful than parametric tests under normality.
#' - May be computationally intensive with large datasets.
#'
#' @param formula A formula of the form \code{y ~ group}, where \code{y} is numeric
#'        and \code{group} is a grouping variable (factor or coercible to factor).
#' @param data A data frame containing the variables in the formula.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"homocedasticidad"}, containing:
#' \describe{
#'   \item{Statistic}{The Fligner-Killeen chi-squared statistic.}
#'   \item{df}{Degrees of freedom.}
#'   \item{p_value}{The p-value for the test.}
#'   \item{Decision}{\code{"Homoscedastic"} or \code{"Heteroscedastic"} depending on the test result.}
#'   \item{Method}{A string indicating the method used ("Fligner-Killeen").}
#' }
#'
#' @references Fligner, M. A., & Killeen, T. J. (1976). "Distribution-free two-sample tests for scale." \emph{Journal of the American Statistical Association}, 71(353), 210–213. <https://doi.org/10.1080/01621459.1976.10480351>
#'
#' @examples
#' data(d_e, package = "Analitica")
#' res <- FKTest(Sueldo_actual ~ labor, data = d_e)
#' summary(res)
#'
#' @export
#' @importFrom stats pchisq model.frame
FKTest <- function(formula, data, alpha = 0.05) {

  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided.")
  }

  mf <- model.frame(formula, data)
  if (ncol(mf) != 2) {
    stop("Formula must be of the form 'response ~ group'.")
  }

  y <- mf[[1]]
  group <- as.factor(mf[[2]])
  k <- nlevels(group)

  if (!is.numeric(y)) {
    stop("The response variable must be numeric.")
  }

  if (k < 2) {
    stop("The grouping variable must have at least two levels.")
  }

  # Median per group
  medians <- tapply(y, group, median)

  # Compute absolute deviations from group medians
  z <- abs(y - medians[group])

  # Rank the deviations
  ranks <- rank(z)

  # Apply log transformation
  log_ranks <- log(ranks + 0.5)

  # Group means and overall mean
  overall_mean <- mean(log_ranks)
  group_means <- tapply(log_ranks, group, mean)
  n <- tapply(log_ranks, group, length)

  # Fligner-Killeen test statistic
  X2 <- sum(n * (group_means - overall_mean)^2)

  df <- k - 1
  p_val <- 1 - pchisq(X2, df)

  decision <- ifelse(p_val < alpha, "Heteroscedastic", "Homoscedastic")
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))

    out <- list(
    Statistic = round(X2, 4),
    df = df,
    p_value = round(p_val, 4),
    Significance = sig,
    Decision = decision,
    Method = "Fligner-Killeen"
  )
  class(out) <- "homocedasticidad"
  return(out)
}
