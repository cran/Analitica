#' Conover-Iman Test for Multiple Comparisons (Non-Parametric)
#'
#' Performs non-parametric pairwise comparisons based on rank-transformed data using
#' the Conover-Iman procedure. This method is typically applied as a post hoc test
#' following a significant Kruskal-Wallis test to identify specific group differences.
#'
#' The Conover-Iman test uses rank-based t-statistics, offering improved statistical
#' power over Dunn's test while maintaining flexibility in sample size.
#'
#' Advantages:
#' - More powerful than Dunn’s test, especially with moderate group differences.
#' - Robust to non-normal data and suitable for ordinal or skewed distributions.
#' - Allows for unequal sample sizes across groups.
#'
#' Disadvantages:
#' - Sensitive to heteroscedasticity (non-constant variances).
#' - Requires appropriate p-value adjustment to control the family-wise error rate.
#'
#' @param formula A formula of the form \code{y ~ group}, where \code{y} is a numeric variable
#'        and \code{group} is a factor indicating group membership.
#' @param data A data frame containing the variables specified in the formula.
#' @param alpha Significance level for hypothesis testing (default is 0.05).
#' @param method.p Method used to adjust p-values for multiple comparisons (default is \code{"holm"}).
#'
#' @return An object of class \code{"conover"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, t-statistics, raw and adjusted p-values, and significance markers.
#'   \item \code{Promedios}: A named numeric vector with mean ranks for each group.
#'   \item \code{Orden_Medias}: A character vector with group names sorted from highest to lowest rank.
#'   \item \code{Metodo}: A string describing the method used ("Conover (no parametrico)").
#' }
#'
#' @references Conover, W. J. & Iman, R. L. (1979). "Multiple comparisons using rank sums." \emph{Technometrics}, 21(4), 489–495.
#'
#' @examples
#'data(d_e, package = "Analitica")
#'ConoverTest(Sueldo_actual ~ labor, data = d_e)
#'
#'
#' @export
#' @importFrom stats pt p.adjust
#' @importFrom utils combn

ConoverTest <- function(formula, data, alpha = 0.05, method.p = "holm") {
  mf <- model.frame(formula, data)
  respuesta <- mf[[1]]
  grupo <- as.factor(mf[[2]])

  niveles <- levels(grupo)
  k <- length(niveles)
  N <- length(respuesta)
  ranks <- rank(respuesta)
  n <- table(grupo)
  Rj <- tapply(ranks, grupo, mean)

  comparaciones <- combn(niveles, 2, simplify = FALSE)
  resultados <- data.frame(
    Comparacion = character(),
    t_value = numeric(),
    p_value = numeric(),
    p_ajustada = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  S2 <- (N * (N + 1)) / 12

  t_vals <- numeric(length(comparaciones))
  p_vals <- numeric(length(comparaciones))

  for (i in seq_along(comparaciones)) {
    par <- comparaciones[[i]]
    g1 <- par[1]; g2 <- par[2]
    dif <- abs(Rj[g1] - Rj[g2])
    se <- sqrt(S2 * (1 / n[g1] + 1 / n[g2]))
    t <- dif / se
    df <- N - k  # Conservative approximation
    p <- 2 * pt(-abs(t), df)

    t_vals[i] <- t
    p_vals[i] <- p
  }

  p_ajustada <- p.adjust(p_vals, method = method.p)
  sig <- ifelse(p_ajustada < 0.001, "***",
                ifelse(p_ajustada < 0.01, "**",
                       ifelse(p_ajustada < 0.05, "*", "ns")))

  for (i in seq_along(comparaciones)) {
    comp <- paste(sort(comparaciones[[i]]), collapse = " - ")
    resultados[i, ] <- list(
      Comparacion = comp,
      t_value = round(t_vals[i], 4),
      p_value = round(p_vals[i], 4),
      p_ajustada = round(p_ajustada[i], 4),
      Significancia = sig[i]
    )
  }

  out <- list(
    Resultados = resultados,
    Promedios = Rj,
    Orden_Medias = names(sort(Rj, decreasing = TRUE)),
    Metodo = "Conover (no parametrico)"
  )
  class(out) <- c("comparaciones", "conover")
  return(out)
}
