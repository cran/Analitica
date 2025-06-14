#' Least Significant Difference (LSD) Test
#'
#' Performs unadjusted pairwise t-tests following a significant ANOVA.
#'
#' Advantages:
#' - Very powerful when assumptions are met.
#' - Simple and easy to interpret.
#'
#' Disadvantages:
#' - High risk of Type I error without correction.
#' - Not recommended if many comparisons are made.
#'
#' @param modelo An object of class \code{aov} or \code{lm}.
#' @param alpha Significance level (default is 0.05).
#'
#' @references
#' Fisher, R. A. (1935). The Design of Experiments. Oliver & Boyd.
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- LSDTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' @return An object of class \code{"comparaciones"} with LSD results.
#' @export
LSDTest <- function(modelo, alpha = 0.05) {
  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  nombres_grupos <- names(medias)

  df_error <- modelo$df.residual
  MSerror <- deviance(modelo) / df_error

  comparaciones <- combn(nombres_grupos, 2, simplify = FALSE)
  resultados <- data.frame(
    Comparacion = character(),
    Diferencia = numeric(),
    Valor_LSD = numeric(),
    p_value = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  for (par in comparaciones) {
    g1 <- par[1]; g2 <- par[2]
    dif <- abs(medias[g1] - medias[g2])
    SE <- sqrt(MSerror * (1 / n[g1] + 1 / n[g2]))
    t_crit <- qt(1 - alpha / 2, df_error)
    LSD_val <- t_crit * SE
    t_val <- dif / SE
    p_val <- 2 * pt(-abs(t_val), df_error)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "ns")))

    resultados <- rbind(resultados, data.frame(
      Comparacion = paste(sort(c(g1, g2)), collapse = " - "),
      Diferencia = round(dif, 4),
      Valor_LSD = round(LSD_val, 4),
      p_value = round(p_val, 4),
      Significancia = sig
    ))
  }

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = names(sort(medias, decreasing = TRUE)),
    Metodo = "LSD"
  )
  class(out) <- c("comparaciones", "lsd")
  return(out)
}
