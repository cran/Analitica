#' Jarque-Bera Test with Glinskiy Corrections
#'
#' Performs the Jarque-Bera test for normality with optional corrections proposed by Glinskiy et al. (2024),
#' depending on whether the mean, variance, or both are known a priori.
#'
#' @param y A numeric vector to test for normality.
#' @param mu Optional known mean value. Default is \code{NULL}.
#' @param sigma2 Optional known variance value. Default is \code{NULL}.
#' @param alpha Significance level for the test (default is 0.05).
#'
#' @return An object of class \code{"normalidad"}, containing:
#' \itemize{
#'   \item \code{statistic}: Test statistic value.
#'   \item \code{df}: Degrees of freedom (always 2).
#'   \item \code{p_value}: P-value of the test.
#'   \item \code{decision}: Conclusion about normality.
#'   \item \code{variant}: Type of JB test applied.
#'   \item \code{method}: "Jarque-Bera (Glinskiy)"
#' }
#'
#' @references Glinskiy, Vladimir & Ismayilova, Yulia & Khrushchev, Sergey & Logachov, Artem & Logachova, Olga & Serga, Lyudmila & Yambartsev, Anatoly & Zaykov, Kirill. (2024). Modifications to the Jarque–Bera Test. Mathematics. 12. 2523. 10.3390/math12162523.
#'
#' @export
#' @importFrom stats pchisq quantile
#' @importFrom moments skewness kurtosis
#'
#' @examples
#'
#' data(d_e, package = "Analitica")
#' JBGTest(d_e$Sueldo_actual)
#' #output different of result
#' summary(JBGTest(d_e$Sueldo_actual))
#'
#'
JBGTest <- function(y, mu = NULL, sigma2 = NULL, alpha = 0.05) {
  if (missing(y)) stop("Argument 'y' must be provided.")

  if (!is.numeric(y)) stop("Input must be a numeric vector.")
  y <- y[!is.na(y)]

  n <- length(y)
  if (n < 8) stop("Sample size must be at least 8.")

  variant <- if (!is.null(mu) && !is.null(sigma2)) {
    "JB(mu, sigma2)"
  } else if (!is.null(mu)) {
    "JB(mu)"
  } else if (!is.null(sigma2)) {
    "JB(sigma2)"
  } else {
    "JB(Classic)"
  }

  if (variant == "JB(Classic)") {
    S <- moments::skewness(y)
    K <- moments::kurtosis(y)
  } else if (variant == "JB(mu)") {
    m2 <- mean((y - mu)^2)
    S <- mean((y - mu)^3) / m2^(3/2)
    K <- mean((y - mu)^4) / m2^2
  } else if (variant == "JB(sigma2)") {
    m <- mean(y)
    S <- mean((y - m)^3) / sigma2^(3/2)
    K <- mean((y - m)^4) / sigma2^2
  } else if (variant == "JB(mu, sigma2)") {
    S <- mean((y - mu)^3) / sigma2^(3/2)
    K <- mean((y - mu)^4) / sigma2^2
  }

  JB <- n/6 * (S^2 + (K - 3)^2 / 4)
  p_val <- 1 - pchisq(JB, df = 2)
  decision <- ifelse(p_val < alpha, "Not Normal", "Normal")

  result <- list(
    statistic = round(JB, 4),
    df = 2,
    p_value = round(p_val, 4),
    decision = decision,
    variant = variant,
    method = "Jarque-Bera (Glinskiy)"
  )
  class(result) <- "normalidad"
  return(result)
}

#' @export
print.normalidad <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("\nVariant:", x$variant)
  cat("\n\nStatistic:", x$statistic, ", df =", x$df, ", p-value =", x$p_value)
  invisible(x)
}

#' @export
summary.normalidad <- function(object, ...) {
  cat("\n\033[1mNormality Test Summary\033[0m\n")
  cat(rep("-", 30), "\n", sep = "")
  cat("Method:   ", object$method, "\n")
  cat("Variant:  ", object$variant, "\n")
  cat("Statistic:", object$statistic, "\n")
  cat("df:       ", object$df, "\n")
  cat("p-value:  ", object$p_value, "\n")
  cat("Decision: ", object$decision, "\n")
  invisible(object)
}





