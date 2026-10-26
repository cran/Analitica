#' Bonferroni-Corrected Mann-Whitney Tests (Non-Parametric)
#'
#' Performs all pairwise comparisons using the Wilcoxon rank-sum test (Mann-Whitney)
#' with Bonferroni correction for multiple testing.
#'
#' Suitable for non-parametric data where ANOVA assumptions are violated.
#' Groups are ordered from highest to lowest mean rank, and each comparison is
#' reported as "higher - lower", so the W statistic always refers to the first
#' group of the pair.
#'
#' Advantages:
#' - Simple and intuitive non-parametric alternative to ANOVA post hoc tests.
#' - Strong control of Type I error via Bonferroni correction.
#' - Works with unequal group sizes.
#'
#' Disadvantages:
#' - Conservative with many groups.
#' - Only valid for pairwise comparisons; does not support complex contrasts.
#' - It is only useful in completely random or single-factor designs.
#'
#' @references
#' Wilcoxon, F. (1945). Individual Comparisons by Ranking Methods. \emph{Biometrics Bulletin}, 1(6), 80–83. \doi{10.2307/3001968}
#'
#' Dunn, O. J. (1964). Multiple Comparisons Using Rank Sums. \emph{Technometrics}, 6(3), 241–252. \doi{10.1080/00401706.1964.10490181}
#'
#' Shaffer, J. P. (1995). Multiple Hypothesis Testing. \emph{Annual Review of Psychology}, 46(1), 561–584. \doi{10.1146/annurev.ps.46.020195.003021}
#'
#' @param formula A formula of the form \code{y ~ group}.
#' @param data A data frame containing the variables.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"bonferroni_np"}, \code{"comparaciones_np"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: Data frame with comparisons (ordered from highest to lowest
#'         mean rank), W-statistics, raw and adjusted p-values, and significance levels.
#'   \item \code{Promedios}: Mean ranks of each group, sorted from highest to lowest.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest mean rank.
#'   \item \code{Metodo}: Name of the method used ("Bonferroni (non-parametric)").
#'   \item \code{Datos}: Data frame with the response and the group, used by \code{plot()}.
#' }
#'
#' @examples
#' data(iris)
#' BonferroniNPTest(Sepal.Length ~ Species, data = iris)
#'
#' @export
#' @importFrom stats wilcox.test p.adjust model.frame
#' @importFrom utils combn

BonferroniNPTest <- function(formula, data, alpha = 0.05) {
  mf <- model.frame(formula, data)
  respuesta <- mf[[1]]
  grupo <- droplevels(as.factor(mf[[2]]))
  if (nlevels(grupo) < 2) stop("At least 2 groups are required.")

  # Rangos globales y rango promedio por grupo
  ranks <- rank(respuesta)
  Rj <- tapply(ranks, grupo, mean)

  # Orden de mayor a menor rango promedio
  Rj_ord <- sort(Rj, decreasing = TRUE)
  niveles_ord <- names(Rj_ord)

  # Pares en ese orden: g1 siempre tiene mayor rango promedio que g2
  comparaciones <- combn(niveles_ord, 2, simplify = FALSE)
  k <- length(comparaciones)

  w_vals <- numeric(k)
  p_vals <- numeric(k)
  etiquetas <- character(k)

  for (i in seq_len(k)) {
    g1 <- comparaciones[[i]][1]
    g2 <- comparaciones[[i]][2]

    test <- wilcox.test(respuesta[grupo == g1], respuesta[grupo == g2],
                        exact = FALSE, correct = FALSE)
    w_vals[i] <- unname(test$statistic)
    p_vals[i] <- test$p.value
    etiquetas[i] <- paste(g1, "-", g2)   # sin sort(): coincide con la dirección de W
  }

  p_ajustada <- p.adjust(p_vals, method = "bonferroni")
  sig <- .sig_stars(p_ajustada, alpha)

  resultados <- data.frame(
    Comparacion   = etiquetas,
    W             = round(w_vals, 4),
    p_value       = round(p_vals, 4),
    p_ajustada    = round(p_ajustada, 4),
    Significancia = sig,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados   = resultados,
    Promedios    = Rj_ord,
    Orden_Medias = niveles_ord,
    Metodo       = "Bonferroni (non-parametric)",
    Datos        = data.frame(y = respuesta, grupo = grupo)
  )
  out <- .agregar_info_pares(out, vapply(comparaciones, `[`, character(1), 1),
                             vapply(comparaciones, `[`, character(1), 2),
                             p_vals, p_ajustada)
  class(out) <- c("bonferroni_np", "comparaciones_np", "comparaciones")
  out
}

