#' Gabriel’s Post Hoc Test for Multiple Comparisons
#'
#' A modification of Tukey's test for use with moderately unequal sample sizes.
#'
#' Advantages:
#' - More powerful than Tukey for unequal group sizes.
#' - Controls error rates effectively with moderate imbalance.
#'
#' Disadvantages:
#' - Can be anti-conservative with large differences in group sizes.
#' - Less common in standard statistical software.
#'
#' @param modelo An object of class \code{aov} or \code{lm}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"gabriel"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: Data frame with comparisons, mean differences, adjusted critical value, p-value, and significance level.
#'   \item \code{Promedios}: Named numeric vector of group means.
#'   \item \code{Orden_Medias}: Vector of group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: Name of the method used ("Gabriel").
#' }
#'
#' @references Hochberg, Y., & Tamhane, A. C. (1987). Multiple Comparison Procedures.
#'
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- GabrielTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#'
#' @export
#' @importFrom stats qtukey ptukey deviance
#' @importFrom utils combn

GabrielTest <- function(modelo, alpha = 0.05) {
  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  nombres_grupos <- names(medias)

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

    # Gabriel SE adjustment
    ni <- n[g1]
    nj <- n[g2]
    N <- ni + nj
    MSE <- MSerror
    se_ij <- sqrt(MSE * (2 / N))

    # Critical value from studentized range
    q_crit <- qtukey(1 - alpha, ng, df_error)
    valor_critico <- q_crit * se_ij / sqrt(2)

    q_obs <- dif * sqrt(2) / se_ij
    p_val <- 1 - ptukey(q_obs, ng, df_error)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "ns")))

    resultados <- rbind(resultados, data.frame(
      Comparacion = paste(sort(c(g1, g2)), collapse = " - "),
      Diferencia = round(dif, 4),
      Valor_Critico = round(valor_critico, 4),
      p_value = round(p_val, 4),
      Significancia = sig
    ))
  }

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = names(sort(medias, decreasing = TRUE)),
    Metodo = "Gabriel"
  )
  class(out) <- c("comparaciones", "gabriel")

  return(out)
}
