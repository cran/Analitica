#' Bonferroni-Corrected Mann-Whitney Tests (Non-Parametric)
#'
#' Performs all pairwise comparisons using the Wilcoxon rank-sum test (Mann-Whitney)
#' with Bonferroni correction for multiple testing.
#'
#' Suitable for non-parametric data where ANOVA assumptions are violated.
#'
#' Advantages:
#' - Simple and intuitive non-parametric alternative to ANOVA post hoc tests.
#' - Strong control of Type I error via Bonferroni correction.
#' - Works with unequal group sizes.
#'
#' Disadvantages:
#' - Conservative with many groups.
#' - Only valid for pairwise comparisons; does not support complex contrasts.
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
#' @return An object of class \code{"bonferroni_np"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: Data frame with comparisons, W-statistics, raw and adjusted p-values, and significance levels.
#'   \item \code{Promedios}: Mean ranks of each group.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest rank.
#'   \item \code{Metodo}: Name of the method used ("Bonferroni (non-parametric)").
#' }
#'
#' @examples
#' data(iris)
#' BonferroniNPTest(Sepal.Length ~ Species, data = iris)
#'
#' @export
#' @importFrom stats wilcox.test p.adjust
#' @importFrom utils combn

BonferroniNPTest <- function(formula, data, alpha = 0.05) {
  mf <- model.frame(formula, data)
  respuesta <- mf[[1]]
  grupo <- as.factor(mf[[2]])

  niveles <- levels(grupo)
  ranks <- rank(respuesta)
  Rj <- tapply(ranks, grupo, mean)

  comparaciones <- combn(niveles, 2, simplify = FALSE)
  resultados <- data.frame(
    Comparacion = character(),
    W = numeric(),
    p_value = numeric(),
    p_ajustada = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  w_vals <- numeric(length(comparaciones))
  p_vals <- numeric(length(comparaciones))

  for (i in seq_along(comparaciones)) {
    par <- comparaciones[[i]]
    g1 <- par[1]; g2 <- par[2]
    datos_g1 <- respuesta[grupo == g1]
    datos_g2 <- respuesta[grupo == g2]

    test <- wilcox.test(datos_g1, datos_g2, exact = FALSE, correct = FALSE)
    w_vals[i] <- test$statistic
    p_vals[i] <- test$p.value
  }

  p_ajustada <- p.adjust(p_vals, method = "bonferroni")
  sig <- ifelse(p_ajustada < 0.001, "***",
                ifelse(p_ajustada < 0.01, "**",
                       ifelse(p_ajustada < 0.05, "*", "ns")))

  for (i in seq_along(comparaciones)) {
    comp <- paste(sort(comparaciones[[i]]), collapse = " - ")
    resultados[i, ] <- list(
      Comparacion = comp,
      W = round(w_vals[i], 4),
      p_value = round(p_vals[i], 4),
      p_ajustada = round(p_ajustada[i], 4),
      Significancia = sig[i]
    )
  }

  out <- list(
    Resultados = resultados,
    Promedios = Rj,
    Orden_Medias = names(sort(Rj, decreasing = TRUE)),
    Metodo = "Bonferroni (non-parametric)"
  )
  class(out) <- c("comparaciones", "bonferroni_np")
  return(out)
}
