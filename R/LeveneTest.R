#' Levene's Test for Homogeneity of Variances (Manual Implementation)
#'
#' Performs Levene's test for equality of variances across groups using a formula interface.
#' This test evaluates the null hypothesis that the variances are equal across groups,
#' and is commonly used as a preliminary test before ANOVA or other parametric analyses.
#'
#' Levene’s test is based on an analysis of variance (ANOVA) applied to the absolute deviations
#' from each group’s center (either the mean or, more robustly, the median).
#' It is less sensitive to departures from normality than Bartlett’s test.
#'
#' Advantages:
#' - Robust to non-normality, especially when using the median.
#' - Suitable for equal or unequal sample sizes across groups.
#' - Widely used in practice for checking homoscedasticity.
#'
#' Disadvantages:
#' - Less powerful than parametric alternatives under strict normality.
#'
#' @param formula A formula of the form \code{y ~ group}, where \code{y} is numeric and \code{group} is a factor.
#' @param data A data frame containing the variables in the formula.
#' @param alpha Significance level (default is 0.05).
#' @param center Character string: use \code{"median"} (default) or \code{"mean"} as the center for deviations.
#'
#' @return An object of class \code{"homocedasticidad"}, containing:
#' \describe{
#'   \item{Statistic}{F statistic of the Levene test.}
#'   \item{df}{Degrees of freedom (between and within groups).}
#'   \item{p_value}{The p-value for the test.}
#'   \item{Decision}{\code{"Homoscedastic"} or \code{"Heteroscedastic"} depending on the test result.}
#'   \item{Method}{A string indicating the method used ("Levene").}
#' }
#'
#' @references Levene, H. (1960). "Robust Tests for Equality of Variances." In Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling (pp. 278–292). Stanford University Press.
#'
#' @examples
#' data(d_e, package = "Analitica")
#' res <- Levene.Test(Sueldo_actual ~ labor, data = d_e)
#' summary(res)
#'
#' @export
#' @importFrom stats model.frame pf
Levene.Test <- function(formula, data, alpha = 0.05, center = "median") {

  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided.")
  }

  if (!center %in% c("median", "mean")) {
    stop("The 'center' argument must be either 'median' or 'mean'.")
  }

  mf <- model.frame(formula, data)
  if (ncol(mf) != 2) {
    stop("Formula must be of the form 'response ~ group'.")
  }

  y <- mf[[1]]
  group <- as.factor(mf[[2]])
  group_levels <- levels(group)

  if (!is.numeric(y)) {
    stop("The response variable must be numeric.")
  }

  if (length(group_levels) < 2) {
    stop("The grouping variable must have at least two levels.")
  }

  # Compute absolute deviations from group center
  z <- numeric(length(y))
  for (g in group_levels) {
    idx <- group == g
    center_value <- if (center == "mean") mean(y[idx]) else median(y[idx])
    z[idx] <- abs(y[idx] - center_value)
  }

  # Degrees of freedom
  df_total <- length(z) - 1
  df_between <- length(group_levels) - 1
  df_within <- df_total - df_between

  # ANOVA on absolute deviations
  mean_z <- tapply(z, group, mean)
  n <- tapply(z, group, length)
  overall_mean <- mean(z)

  SS_between <- sum(n * (mean_z - overall_mean)^2)
  group_vector <- rep(mean_z, times = n)
  SS_within <- sum((z - group_vector)^2)

  MS_between <- SS_between / df_between
  MS_within <- SS_within / df_within

  F_stat <- MS_between / MS_within
  p_val <- 1 - pf(F_stat, df_between, df_within)

  decision <- ifelse(p_val < alpha, "Heteroscedastic", "Homoscedastic")
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))


  out <- list(
    Statistic = round(F_stat, 4),
    df = c(df_between = df_between, df_within = df_within),
    p_value = round(p_val, 4),
    Significance = sig,
    Decision = decision,
    Method = paste("Levene (", center, ")", sep = "")
  )
  class(out) <- "homocedasticidad"
  return(out)
}

