#' Duncan Multiple Range Test (DMRT)
#'
#' Performs the Duncan test for pairwise comparisons after an ANOVA.
#' This method is more liberal than Tukey's HSD, using a stepwise approach
#' with critical values from the studentized range distribution.
#'
#' Advantages:
#' - High power for detecting differences.
#' - Simple to interpret and implement.
#'
#' Disadvantages:
#' - Inflates Type I error rate.
#' - Not recommended for confirmatory research.
#'
#' @param modelo An object of class \code{aov} or \code{lm}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"duncan"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, mean differences,
#'   critical values, p-values, and significance indicators.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector with group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string indicating the comparison method ("Duncan").
#' }
#'
#' @references Duncan, D. B. (1955). "Multiple range and multiple F tests." Biometrics, 11(1), 1-42.
#'
#' @export
#' @importFrom stats qtukey ptukey qf
#' @importFrom utils combn
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- DuncanTest(mod)
#' summary(resultado)
#' plot(resultado)
DuncanTest <- function(modelo, alpha = 0.05) {
  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  nombres_grupos <- names(medias)

  df_error <- modelo$df.residual
  MSerror <- deviance(modelo) / df_error

  # Función para error estándar
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

  # Ordenar las medias de menor a mayor
  ordenado <- sort(medias)
  grupo_nombres_ordenados <- names(ordenado)

  for (par in comparaciones) {
    g1 <- par[1]
    g2 <- par[2]

    # Número de medias intermedias + 2
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

  # Ordenar nombres de grupos por media (de mayor a menor)
  ordenado_nombres <- names(sort(medias, decreasing = TRUE))

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = ordenado_nombres,
    Metodo = "Duncan"
  )
  class(out) <- c("comparaciones", "duncan")

  return(out)
}

