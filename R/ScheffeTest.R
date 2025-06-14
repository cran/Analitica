#' Scheffé Test for Multiple Comparisons
#'
#' Performs Scheffé's post hoc test after fitting an ANOVA model. This test compares all possible
#' pairs of group means, using a critical value based on the F-distribution.
#'
#' The Scheffé test is a conservative method, making it harder to detect significant differences,
#' but reducing the likelihood of Type I errors (false positives). It is especially appropriate
#' when the comparisons were not pre-planned and the number of contrasts is large.
#'
#' Assumptions: normally distributed residuals and homogeneity of variances.
#'
#' Advantages:
#' - Very robust to violations of assumptions.
#' - Suitable for complex comparisons, not just pairwise.
#'
#' Disadvantages:
#' - Very conservative; reduced power.
#' - Not ideal for detecting small differences.
#'
#' @param modelo An object of class \code{aov} or \code{lm} representing an ANOVA model.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"scheffe"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame of pairwise comparisons with difference, critical value, p-value, and significance code.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector with group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string indicating the test name ("Scheffe").
#' }
#'
#' @references Scheffé, H. (1953). "A method for judging all contrasts in the analysis of variance." \emph{Biometrika}, 40(1/2), 87–104. <https://doi.org/10.1093/biomet/40.1-2.87>
#'
#' @importFrom  stats qf
#' @importFrom utils combn
#' @export
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- ScheffeTest(mod)
#' summary(resultado)
#' plot(resultado)
ScheffeTest <- function(modelo, alpha = 0.05) {

  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  nombres_grupos <- names(medias)

  orden_medias <- order(medias, decreasing = TRUE)
  etiquetas_ordenadas <- nombres_grupos[orden_medias]

  v1 <- modelo$rank - 1
  v2 <- modelo$df.residual
  MSerror <- deviance(modelo) / v2
  Fcrit <- qf(1 - alpha, v1, v2)

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
    SE <- MSerror * (1 / n[g1] + 1 / n[g2])
    Fobs <- (dif^2) / SE
    valor_critico <- sqrt(v1 * Fcrit * SE)
    p_val <- 1 - pf(Fobs, v1, v2)

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
    Metodo = "Scheffe"
  )
  class(out) <- c("comparaciones", "scheffe")

  return(out)
}

