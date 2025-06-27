#' Jonckheere-Terpstra Test for Ordered Alternatives (with Tie Correction)
#'
#' Performs the Jonckheere-Terpstra test to evaluate the presence of a monotonic trend
#' (increasing or decreasing) across three or more independent ordered groups.
#' This test is non-parametric and is particularly useful when the independent variable
#' is ordinal and the response is continuous or ordinal.
#'
#' The Jonckheere-Terpstra test compares all pairwise combinations of groups and counts
#' the number of times values in higher-ordered groups exceed those in lower-ordered groups.
#' This implementation includes a full correction for ties in the data, which ensures more accurate inference.
#'
#' Advantages:
#' - Non-parametric: does not assume normality or equal variances.
#' - More powerful than Kruskal-Wallis when there is an a priori ordering of groups.
#' - Tie correction included, improving robustness in real-world data.
#'
#' Disadvantages:
#' - Requires that the group variable be ordered (ordinal).
#' - Detects overall trend but not specific group differences.
#' - Sensitive to large numbers of ties or very unbalanced group sizes.
#'
#' @param formula A formula of the type y ~ group, where `group` is an ordered factor.
#' @param data A data.frame containing the variables in the formula.
#'
#' @return An object of class "jonckheere" with:
#' \itemize{
#'   \item \code{J}: Total Jonckheere-Terpstra statistic.
#'   \item \code{J_pares}: Pairwise J statistics between group combinations.
#'   \item \code{mu_J}: Expected value of J under the null hypothesis.
#'   \item \code{var_J}: Variance of J (with complete tie correction).
#'   \item \code{Z}: Standardized test statistic.
#'   \item \code{p_value}: Two-sided p-value.
#'   \item \code{Trend}: Detected trend ("increasing", "decreasing", or "none").
#'   \item \code{Method}: Description of the method.
#' }
#'
#' @references Hollander, M., Wolfe, D. A., & Chicken, E. (2014). Nonparametric statistical methods. p. 202 (3rd ed.). Wiley.
#'
#' @export
#' @importFrom stats pnorm
#'
#' @examples
#' df <- data.frame(
#'   group = factor(rep(1:3, each = 6), ordered = TRUE),
#'   y = c(40,35,38,43,44,41,38,40,47,44,40,42,48,40,45,43,46,44)
#' )
#' res <- JT_Test(y ~ group, data = df)
#'
JT_Test <- function(formula, data) {
  mf <- model.frame(formula, data)
  y <- mf[[1]]
  group <- mf[[2]]

  if (!is.ordered(group)) {
    warning("The group variable should be an ordered factor. Alphabetical order will be used.")
    group <- factor(group, ordered = TRUE)
  }

  levels_group <- levels(group)
  k <- length(levels_group)
  groups <- split(y, group)
  n_i <- sapply(groups, length)
  N <- length(y)

  # Compute total and pairwise J
  J <- 0
  J_pairs <- c()
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      g1 <- groups[[i]]
      g2 <- groups[[j]]
      Jij <- 0
      for (x in g1) {
        for (y_ in g2) {
          if (x < y_) {
            Jij <- Jij + 1
          } else if (x == y_) {
            Jij <- Jij + 0.5
          }
        }
      }
      J <- J + Jij
      pair_name <- paste(levels_group[i], "vs", levels_group[j])
      J_pairs[pair_name] <- Jij
    }
  }

  # Expected value
  mu_J <- 0.5 * sum(outer(n_i, n_i, "*")[lower.tri(matrix(0, k, k))])

  # Frequencies for ties
  ties <- table(y)

  # Full tie correction: A, B, C
  A <- (N * (N - 1) * (2 * N + 5) -
          sum(n_i * (n_i - 1) * (2 * n_i + 5)) -
          sum(ties * (ties - 1) * (2 * ties + 5))) / 72

  B <- (sum(n_i * (n_i - 1) * (n_i - 2)) *
          sum(ties * (ties - 1) * (ties - 2))) /
    (36 * N * (N - 1) * (N - 2))

  C <- (sum(n_i * (n_i - 1)) *
          sum(ties * (ties - 1))) /
    (8 * N * (N - 1))

  var_J <- A + C - B

  # Z statistic and p-value
  Z <- (J - mu_J) / sqrt(var_J)
  pval <- 2 * (1 - pnorm(abs(Z)))

  trend <- if (Z > 0 & pval < 0.05) {
    "increasing"
  } else if (Z < 0 & pval < 0.05) {
    "decreasing"
  } else {
    "none"
  }

  result <- list(
    J = J,
    J_pares = J_pairs,
    mu_J = mu_J,
    var_J = var_J,
    Z = Z,
    p_value = pval,
    Trend = trend,
    Method = "Jonckheere-Terpstra (full tie correction)"
  )
  class(result) <- "jonckheere"

  # Display summary
  cat("Jonckheere-Terpstra Test\n")
  cat("------------------------\n")
  cat("Method: ", result$Method, "\n")
  cat("Total J statistic:", round(result$J, 4), "\n")
  cat("Pairwise J decomposition:\n")
  print(round(result$J_pares, 4))
  cat("E(J):", round(result$mu_J, 4), "\n")
  cat("V(J):", round(result$var_J, 4), "\n")
  cat("Z:", round(result$Z, 4), "\n")
  cat("p-value:", format.pval(result$p_value, digits = 4), "\n")
  cat("Trend:", result$Trend, "\n\n")

  return(result)
}
