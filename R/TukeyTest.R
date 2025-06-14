#' Tukey HSD Test for Multiple Comparisons
#'
#' Performs Tukey's Honest Significant Difference (HSD) test for all pairwise
#' comparisons after fitting an ANOVA model. This post hoc method uses the
#' studentized range distribution and is appropriate when variances are equal
#' across groups and observations are independent.
#'
#' Tukey's test controls the family-wise error rate and is widely used when group
#' comparisons have not been planned in advance.
#'
#' Advantages:
#' - Strong control of Type I error rate.
#' - Ideal for balanced designs with equal variances.
#'
#' Disadvantages:
#' - Assumes equal variances and sample sizes.
#' - Less powerful with heteroscedasticity.
#'
#' @param modelo An object of class \code{aov} or \code{lm} representing an ANOVA model.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"tukey"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame of pairwise comparisons with mean differences, critical value, p-value, and significance level.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector with group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string indicating the method used ("Tukey").
#' }
#'
#' @references Tukey, J. W. (1949). "Comparing individual means in the analysis of variance." \emph{Biometrics}, 5(2), 99–114. <https://doi.org/10.2307/3001913>
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- TukeyTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' @export
#' @importFrom stats qtukey ptukey deviance qf
#' @importFrom utils combn
TukeyTest <- function(modelo, alpha = 0.05) {

  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  nombres_grupos <- names(medias)

  orden_medias <- order(medias, decreasing = TRUE)
  etiquetas_ordenadas <- nombres_grupos[orden_medias]

  df_error <- modelo$df.residual
  MSerror <- deviance(modelo) / df_error
  ng <- length(medias)

  comparaciones <- combn(nombres_grupos, 2, simplify = FALSE)

  resultados <- data.frame(
    Comparacion = character(),
    Diferencia = numeric(),
    Valor_Critico = numeric(),
    p_value = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  for (par in comparaciones) {
    g1 <- par[1]
    g2 <- par[2]

    dif <- abs(medias[g1] - medias[g2])
    SE <- sqrt(MSerror * (1 / n[g1] + 1 / n[g2]))

    q_crit <- qtukey(1 - alpha, ng, df_error)
    valor_critico <- q_crit * SE / sqrt(2)

    q_obs <- dif * sqrt(2) / SE
    p_val <- 1 - ptukey(q_obs, ng, df_error)

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

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = etiquetas_ordenadas,
    Metodo = "Tukey"
  )
  class(out) <- c("comparaciones", "tukey")

  return(out)
}

