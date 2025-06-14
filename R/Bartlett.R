#' Bartlett's Test for Homogeneity of Variances (Manual Implementation)
#'
#' Conducts Bartlett's test to evaluate whether multiple groups have equal variances,
#' based on a formula interface and raw data vectors, without requiring a fitted model.
#' This implementation provides flexibility for exploratory variance testing in custom workflows.
#'
#' Bartlett’s test is appropriate when group distributions are approximately normal.
#' It tests the null hypothesis that all groups have equal variances (homoscedasticity).
#'
#' Advantages:
#' - Straightforward to compute.
#' - High sensitivity to variance differences under normality.
#'
#' Disadvantages:
#' - Highly sensitive to non-normal distributions.
#' - Less robust than alternatives like Levene’s test for skewed or heavy-tailed data.
#'
#' @param formula A formula of the form \code{y ~ group}, where \code{y} is a numeric response
#'        and \code{group} is a factor indicating group membership.
#' @param data A data frame containing the variables specified in the formula.
#' @param alpha Significance level for the test (default is 0.05).
#'
#' @return An object of class \code{"homocedasticidad"}, containing:
#' \itemize{
#'   \item \code{Statistic}: Bartlett's chi-squared test statistic.
#'   \item \code{df}: Degrees of freedom associated with the test.
#'   \item \code{p_value}: The p-value for the test statistic.
#'   \item \code{Decision}: A character string indicating the conclusion ("Heterocedastic" or "Homocedastic").
#'   \item \code{Method}: A character string indicating the method used ("Bartlett").
#' }
#'
#' @references Bartlett, M. S. (1937). "Properties of sufficiency and statistical tests." \emph{Proceedings of the Royal Society of London}, Series A, 160(901), 268–282.
#'
#' @examples
#' data(d_e, package = "Analitica")
#' res <- BartlettTest(Sueldo_actual ~ labor, data = d_e)
#' summary(res)
#'
#' summary(BartlettTest(Sueldo_actual ~ as.factor(labor), data = d_e))
#'
#'
#' @export
#' @importFrom stats pchisq var model.frame
BartlettTest <- function(formula, data, alpha = 0.05) {

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

  k <- nlevels(group)
  if (k < 2) {
    stop("The grouping variable must have at least two levels.")
  }

  n <- tapply(y, group, length)
  s2 <- tapply(y, group, var)
  N <- sum(n)

  # Bartlett's test calculation
  sum_num <- sum((n - 1) * log(s2))
  s2_pooled <- sum((n - 1) * s2) / (N - k)

  A <- (N - k) * log(s2_pooled) - sum_num
  C <- 1 + (1 / (3 * (k - 1))) * (sum(1 / (n - 1)) - 1 / (N - k))
  X2 <- A / C

  df <- k - 1
  p_val <- 1 - pchisq(X2, df)
  decision <- ifelse(p_val < alpha, "Heterocedastic", "Homocedastic")

  # Significance stars
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))

  out <- list(
    Statistic = round(X2, 4),
    df = df,
    p_value = round(p_val, 4),
    Significance = sig,
    Decision = decision,
    Method = "Bartlett"
  )
  class(out) <- "homocedasticidad"
  return(out)
}
