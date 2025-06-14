#' Bonferroni-Corrected Pairwise t-Tests
#'
#' Performs pairwise t-tests with Bonferroni adjustment for multiple comparisons.
#' This method controls the family-wise error rate by dividing the alpha level
#' by the number of comparisons.
#'
#' Advantages:
#' - Very simple and easy to implement.
#' - Strong control of Type I error.
#' - Applicable to any set of independent comparisons.
#'
#' Disadvantages:
#' - Highly conservative, especially with many groups.
#' - Can lead to low statistical power (increased Type II error).
#' - Does not adjust test statistics, only p-values.
#'
#' @references
#' Dunn, O. J. (1964). Multiple Comparisons Using Rank Sums. \emph{Technometrics}, 6(3), 241–252. \doi{10.1080/00401706.1964.10490181}
#'
#' Wilcoxon, F. (1945). Individual Comparisons by Ranking Methods. \emph{Biometrics Bulletin}, 1(6), 80–83. \doi{10.2307/3001968}
#'
#' @param modelo An object of class \code{aov} or \code{lm}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"bonferroni"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: Data frame with comparisons, mean differences, t-values, unadjusted and adjusted p-values, and significance.
#'   \item \code{Promedios}: Named numeric vector of group means.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: Name of the method used ("Bonferroni-adjusted t-test").
#' }
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- BonferroniTest(mod)
#' summary(resultado)
#'
#' @export
#' @importFrom stats pt p.adjust deviance
#' @importFrom utils combn

BonferroniTest <- function(modelo, alpha = 0.05) {
  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  s2 <- deviance(modelo) / modelo$df.residual
  nombres_grupos <- names(medias)
  df_error <- modelo$df.residual

  comparaciones <- combn(nombres_grupos, 2, simplify = FALSE)
  raw_pvals <- numeric(length(comparaciones))
  t_vals <- numeric(length(comparaciones))
  difs <- numeric(length(comparaciones))
  SEs <- numeric(length(comparaciones))
  nombres_comp <- character(length(comparaciones))

  for (i in seq_along(comparaciones)) {
    par <- comparaciones[[i]]
    g1 <- par[1]; g2 <- par[2]
    dif <- abs(medias[g1] - medias[g2])
    SE <- sqrt(s2 * (1 / n[g1] + 1 / n[g2]))
    t_val <- dif / SE
    p_val <- 2 * pt(-abs(t_val), df_error)

    t_vals[i] <- t_val
    raw_pvals[i] <- p_val
    difs[i] <- dif
    SEs[i] <- SE
    nombres_comp[i] <- paste(sort(c(g1, g2)), collapse = " - ")
  }

  pvals_ajustados <- p.adjust(raw_pvals, method = "bonferroni")

  significancia <- ifelse(pvals_ajustados < 0.001, "***",
                          ifelse(pvals_ajustados < 0.01, "**",
                                 ifelse(pvals_ajustados < 0.05, "*", "ns")))

  resultados <- data.frame(
    Comparacion = nombres_comp,
    Diferencia = round(difs, 4),
    t_value = round(t_vals, 4),
    p_value = round(raw_pvals, 4),
    p_ajustada = round(pvals_ajustados, 4),
    Significancia = significancia,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = names(sort(medias, decreasing = TRUE)),
    Metodo = "Bonferroni-adjusted t-test"
  )
  class(out) <- c("comparaciones", "bonferroni")

  return(out)
}
