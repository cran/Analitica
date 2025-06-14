#' Student-Newman-Keuls (SNK) Test for Multiple Comparisons
#'
#' Performs the Student-Newman-Keuls (SNK) post hoc test for pairwise comparisons
#' after fitting an ANOVA model. The test uses a stepwise approach where the
#' critical value depends on the number of means spanned between groups (range r).
#'
#' SNK is more powerful but less conservative than Tukey’s HSD, increasing the chance of
#' detecting real differences while slightly raising the Type I error rate.
#'
#' Assumptions: normality, homogeneity of variances, and independence of observations.
#'
#' Advantages:
#' - More powerful than Tukey when differences are large.
#' - Intermediate control of Type I error.
#'
#' Disadvantages:
#' - Error control is not family-wise.
#' - Type I error increases with more comparisons.
#'
#'
#' @param modelo An object of class \code{aov} or \code{lm} representing an ANOVA model.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"snk"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, including mean differences,
#'   critical values, p-values, and significance codes.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector with group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string indicating the test used ("SNK").
#' }
#'
#' @references Student, Newman, and Keuls (1952). "Student-Newman-Keuls Procedure".
#' See also: <https://doi.org/10.1002/bimj.200310019>
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- SNKTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' @export
#' @importFrom stats ptukey qtukey deviance qf
#' @importFrom utils combn
#'
SNKTest <- function(modelo, alpha = 0.05) {
  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  nombres_grupos <- names(medias)

  df_error <- modelo$df.residual
  MSerror <- deviance(modelo) / df_error

  # Función para el error estándar
  SE <- function(g1, g2) sqrt(MSerror * (1 / n[g1] + 1 / n[g2]))

  comparaciones <- combn(nombres_grupos, 2, simplify = FALSE)

  resultados <- data.frame(
    Comparacion = character(),
    Diferencia = numeric(),
    Valor_Critico = numeric(),
    p_value = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  # Ordenar medias de menor a mayor
  ordenado <- sort(medias)
  grupo_nombres_ordenados <- names(ordenado)

  for (par in comparaciones) {
    g1 <- par[1]
    g2 <- par[2]

    # r = número de medias entre los grupos (inclusive)
    pos1 <- which(grupo_nombres_ordenados == g1)
    pos2 <- which(grupo_nombres_ordenados == g2)
    r <- abs(pos2 - pos1) + 1

    dif <- abs(medias[g1] - medias[g2])
    se_val <- SE(g1, g2)
    q_crit <- qtukey(1 - alpha, r, df_error)
    valor_critico <- q_crit * se_val / sqrt(2)

    q_obs <- dif * sqrt(2) / se_val
    p_val <- 1 - ptukey(q_obs, r, df_error)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "ns")))

    comparacion <- paste(sort(c(g1, g2)), collapse = " - ")

    resultados <- rbind(resultados, data.frame(
      Comparacion = comparacion,
      Diferencia = round(dif, 4),
      Valor_Critico = round(valor_critico, 4),
      p_value = round(p_val, 4),
      Significancia = sig
    ))
  }

  ordenado_nombres <- names(sort(medias, decreasing = TRUE))

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = ordenado_nombres,
    Metodo = "SNK"
  )
  class(out) <- c("comparaciones", "snk")

  return(out)
}
