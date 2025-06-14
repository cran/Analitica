#' Dunn's Test for Multiple Comparisons (Non-Parametric)
#'
#' Performs Dunn's test for pairwise comparisons following a Kruskal-Wallis test.
#' Suitable for non-parametric data (ordinal or non-normal), using rank sums.
#' Includes Holm correction by default for multiple comparisons.
#'
#' Advantages:
#' - Simple and widely used non-parametric alternative to Tukey's test.
#' - Handles unequal sample sizes.
#' - Compatible with various p-value corrections (e.g., Holm, Bonferroni).
#'
#' Disadvantages:
#' - Less powerful than DSCF or Conover when sample sizes vary widely.
#' - Requires ranking all data and can be conservative depending on adjustment.
#'
#' @param formula A formula of the form \code{y ~ group}.
#' @param data A data frame containing the variables.
#' @param alpha Significance level (default is 0.05).
#' @param method.p Method for p-value adjustment (default is "holm").
#'
#' @return An object of class \code{"dunn"} and \code{"comparaciones"}, including:
#' \itemize{
#'   \item \code{Resultados}: Data frame with group comparisons, z-values, raw and adjusted p-values, and significance.
#'   \item \code{Promedios}: Mean ranks of each group.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest rank.
#'   \item \code{Metodo}: "Dunn (no paramétrico)".
#' }
#'
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics*, 6(3), 241–252. \doi{10.1080/00401706.1964.10490181}
#'
#' @seealso \code{\link{kruskal.test}}, \code{\link[dunn.test]{dunn.test}}
#'
#'@examples
#'data(d_e, package = "Analitica")
#'DunnTest(Sueldo_actual ~ labor, data = d_e)
#'
#'
#' @export
#' @importFrom stats p.adjust pnorm kruskal.test
#' @importFrom utils combn

DunnTest <- function(formula, data, alpha = 0.05, method.p = "holm") {
  mf <- model.frame(formula, data)
  response <- mf[[1]]
  grupo <- as.factor(mf[[2]])

  # Kruskal-Wallis check (opcional, pero puede ayudarte a verificar)
  kruskal <- kruskal.test(formula, data = data)

  # Ranks y estructura
  ranks <- rank(response)
  niveles <- levels(grupo)
  k <- length(niveles)
  n <- table(grupo)
  N <- length(response)
  Rj <- tapply(ranks, grupo, mean)
  S2 <- (N * (N + 1)) / 12

  comparaciones <- combn(niveles, 2, simplify = FALSE)
  resultados <- data.frame(
    Comparacion = character(),
    z = numeric(),
    p_value = numeric(),
    p_ajustada = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  z_vals <- numeric(length(comparaciones))
  p_vals <- numeric(length(comparaciones))

  for (i in seq_along(comparaciones)) {
    par <- comparaciones[[i]]
    g1 <- par[1]; g2 <- par[2]

    dif <- abs(Rj[g1] - Rj[g2])
    se <- sqrt(S2 * (1 / n[g1] + 1 / n[g2]))
    z <- dif / se
    p <- 2 * (1 - pnorm(abs(z)))

    z_vals[i] <- z
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
      z = round(z_vals[i], 4),
      p_value = round(p_vals[i], 4),
      p_ajustada = round(p_ajustada[i], 4),
      Significancia = sig[i]
    )
  }

  ordenado <- names(sort(Rj, decreasing = TRUE))

  out <- list(
    Resultados = resultados,
    Promedios = Rj,
    Orden_Medias = ordenado,
    Metodo = "Dunn (no parametrico)"
  )
  class(out) <- c("comparaciones", "dunn")
  return(out)
}
