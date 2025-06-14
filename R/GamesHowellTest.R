#' Games-Howell Post Hoc Test
#'
#' Performs the Games-Howell test for pairwise comparisons after ANOVA,
#' without assuming equal variances or sample sizes. It is suitable when
#' Levene or Bartlett test indicates heterogeneity of variances.
#'
#' Advantages:
#' - Excellent for heteroscedastic data.
#' - Controls Type I error across unequal group sizes.
#'
#' Disadvantages:
#' - Slightly conservative in small samples.
#' - More complex to compute than Tukey.

#'
#' @param modelo An object from \code{aov} or \code{lm}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"gameshowell"} and \code{"comparaciones"},
#' which contains:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, including mean differences, t-values, degrees of freedom, p-values, and significance labels.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector with group names ordered by their means.
#'   \item \code{Metodo}: A character string indicating the method used ("Games-Howell").
#' }
#'
#' @references Games, P. A., & Howell, J. F. (1976). "Pairwise Multiple Comparison Procedures with Unequal N's and/or Variances: A Monte Carlo Study". Journal of Educational Statistics, 1(2), 113–125. <https://doi.org/10.1002/j.2162-6057.1976.tb00211.x>
#'
#' @importFrom stats qf
#' @importFrom utils combn
#'
#' @export
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- GHTest(mod)
#' summary(resultado)
#' plot(resultado)
GHTest <- function(modelo, alpha = 0.05) {
  factor_name <- names(modelo$xlevels)[1]
  grupos <- modelo$model[[factor_name]]
  respuesta <- modelo$model[[1]]

  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  s2 <- tapply(respuesta, grupos, var)
  nombres_grupos <- names(medias)

  comparaciones <- combn(nombres_grupos, 2, simplify = FALSE)

  resultados <- data.frame(
    Comparacion = character(),
    Diferencia = numeric(),
    t_value = numeric(),
    gl = numeric(),
    p_value = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  for (par in comparaciones) {
    g1 <- par[1]
    g2 <- par[2]

    dif <- abs(medias[g1] - medias[g2])
    se_ij <- sqrt((s2[g1] / n[g1]) + (s2[g2] / n[g2]))

    # Welch-Satterthwaite degrees of freedom
    df_num <- (s2[g1] / n[g1] + s2[g2] / n[g2])^2
    df_den <- ((s2[g1]^2) / (n[g1]^2 * (n[g1] - 1))) + ((s2[g2]^2) / (n[g2]^2 * (n[g2] - 1)))
    gl <- df_num / df_den

    t_val <- dif / se_ij
    p_val <- 2 * pt(-abs(t_val), df = gl)  # two-tailed

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "ns")))


    # Ordenar nombres para coincidir con plot.comparacion()
    nombres_ordenados <- sort(c(g1, g2))
    comparacion <- paste(nombres_ordenados, collapse = " - ")

    resultados <- rbind(resultados, data.frame(
      Comparacion = comparacion,
      Diferencia = round(dif, 4),
      t_value = round(t_val, 4),
      gl = round(gl, 2),
      p_value = round(p_val, 4),
      Significancia = sig
    ))
  }

  # Ordenar los nombres de los grupos por media (de mayor a menor)
  ordenado <- names(sort(medias, decreasing = TRUE))

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = ordenado,
    Metodo = "Games-Howell"
  )
  class(out) <- c("comparaciones", "gameshowell")

  return(out)
}
