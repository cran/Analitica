#' Dunnett's T3 Post Hoc Test
#'
#' Performs Dunnett's T3 test for pairwise comparisons after an ANOVA model.
#' This test is recommended when group variances are unequal and sample sizes differ.
#' It is based on the studentized range distribution and provides conservative control
#' over Type I error without assuming homoscedasticity.
#'
#'
#' Advantages:
#' - More powerful than T2 when group sizes are small.
#' - Adjusted for unequal variances.
#'
#' Disadvantages:
#' - Complex critical value estimation.
#' - Less frequently used and harder to find in software.

#'
#' @param modelo An object of class \code{aov} or \code{lm}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"dunnettt3"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, mean differences,
#'   q-values, degrees of freedom, p-values, and significance indicators.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector of group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string with the test name ("Dunnett T3").
#' }
#'
#' @references Dunnett, C. W. (1980). "Pairwise multiple comparisons in the unequal variance case." \emph{Journal of the American Statistical Association}, 75(372), 796–800. <https://doi.org/10.1080/01621459.1980.10477558>
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- T3Test(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' @export
#' @importFrom stats ptukey qf
#' @importFrom utils combn
T3Test <- function(modelo, alpha = 0.05) {
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
    q_value = numeric(),
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

    # q observed
    q_val <- dif / se_ij * sqrt(2)

    # p-value from studentized range distribution
    p_val <- 1 - ptukey(q_val, nmeans = length(nombres_grupos), df = gl)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "ns")))

    comparacion <- paste(sort(c(g1, g2)), collapse = " - ")

    resultados <- rbind(resultados, data.frame(
      Comparacion = comparacion,
      Diferencia = round(dif, 4),
      q_value = round(q_val, 4),
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
    Metodo = "Dunnett T3"
  )
  class(out) <- c("comparaciones", "dunnettt3")

  return(out)
}
