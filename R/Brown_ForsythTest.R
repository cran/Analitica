#' Brown-Forsythe Test for Homogeneity of Variances (Manual Implementation)
#'
#' Performs the Brown-Forsythe test using absolute deviations from the median
#' of each group, followed by a one-way ANOVA on those deviations.
#'
#' This test is a robust alternative to Bartlett's test, especially useful when
#' the assumption of normality is violated or when outliers are present.
#'
#' Advantages:
#' - More robust than Bartlett's test under non-normal distributions.
#' - Less sensitive to outliers due to the use of the median.
#'
#' Disadvantages:
#' - Lower power than Bartlett's test when normality strictly holds.
#' - Assumes that absolute deviations follow similar distributions across groups.
#'
#' @param formula A formula of the form \code{y ~ group}, where \code{y} is numeric
#'        and \code{group} is a factor.
#' @param data A data frame containing the variables.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"homocedasticidad"}, with:
#' \itemize{
#'   \item \code{Statistic}: F-statistic.
#'   \item \code{df1}: Numerator degrees of freedom.
#'   \item \code{df2}: Denominator degrees of freedom.
#'   \item \code{p_value}: P-value.
#'   \item \code{Decision}: "Heterocedastic" or "Homocedastic".
#'   \item \code{Method}: "Brown-Forsythe".
#' }
#'
#'@references Brown, M. B., & Forsythe, A. B. (1974). "Robust Tests for the Equality of Variances".
#' \emph{Journal of the American Statistical Association}, 69(346), 364–367.
#'
#'
#' @examples
#' data(d_e, package = "Analitica")
#' res <- BrownForsytheTest(Sueldo_actual ~ labor, data = d_e)
#' summary(res)
#'
#'
#' @export
#' @importFrom stats model.frame median lm anova
BrownForsytheTest <- function(formula, data, alpha = 0.05) {

  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided.")
  }

  mf <- model.frame(formula, data)
  if (ncol(mf) != 2) {
    stop("The formula must be of the form 'response ~ group'.")
  }

  y <- mf[[1]]
  group <- as.factor(mf[[2]])

  if (!is.numeric(y)) {
    stop("The response variable must be numeric.")
  }

  if (nlevels(group) < 2) {
    stop("The grouping variable must have at least two levels.")
  }

  # Absolute deviations from group medians
  medians <- tapply(y, group, median, na.rm = TRUE)
  deviations <- abs(y - medians[group])

  # ANOVA on deviations
  fit <- lm(deviations ~ group)
  aov_res <- anova(fit)

  F_stat <- aov_res$`F value`[1]
  df1 <- aov_res$Df[1]
  df2 <- aov_res$Df[2]
  p_val <- aov_res$`Pr(>F)`[1]
  decision <- ifelse(p_val < alpha, "Heterocedastic", "Homocedastic")

  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))

  out <- list(
    Statistic = round(F_stat, 4),
    df1 = df1,
    df2 = df2,
    p_value = round(p_val, 4),
    Significance = sig,
    Decision = decision,
    Method = "Brown-Forsythe"
  )
  class(out) <- "homocedasticidad"
  return(out)
}
