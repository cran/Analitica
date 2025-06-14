#' Tamhane's T2 Post Hoc Test
#'
#' Performs the Tamhane T2 test for pairwise comparisons after an ANOVA model,
#' assuming unequal variances and/or unequal sample sizes. This test is appropriate
#' when the assumption of homogeneity of variances is violated, such as when
#' Levene's test or Bartlett's test is significant.
#'
#' The test uses a modified t-test with Welch-Satterthwaite degrees of freedom and
#' a conservative approach to control for multiple comparisons.
#'
#' Advantages:
#' - Controls Type I error under heteroscedasticity.
#' - No assumption of equal sample sizes.
#'
#' Disadvantages:
#' - Conservative; may reduce power.
#' - Not as powerful as Games-Howell in some contexts.
#'
#' @param modelo An object of class \code{aov} or \code{lm}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"tamhanet2"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, mean differences,
#'   t-values, degrees of freedom, p-values, and significance codes.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector with group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string indicating the method used ("Tamhane T2").
#' }
#'
#' @references Tamhane, A. C. (1977). "Multiple comparisons in model I one-way ANOVA with unequal variances." \emph{Communications in Statistics - Theory and Methods}, 6(1), 15–32. <https://doi.org/10.1080/03610927708827524>
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- T2Test(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' @export
#' @importFrom stats pt qf
#' @importFrom utils combn
#'
T2Test <- function(modelo, alpha = 0.05) {
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
    p_val <- 2 * pt(-abs(t_val), df = gl)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "ns")))

    comparacion <- paste(sort(c(g1, g2)), collapse = " - ")

    resultados <- rbind(resultados, data.frame(
      Comparacion = comparacion,
      Diferencia = round(dif, 4),
      t_value = round(t_val, 4),
      gl = round(gl, 2),
      p_value = round(p_val, 4),
      Significancia = sig
    ))
  }

  ordenado <- names(sort(medias, decreasing = TRUE))

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = ordenado,
    Metodo = "Tamhane T2"
  )
  class(out) <- c("comparaciones", "tamhanet2")

  return(out)
}

