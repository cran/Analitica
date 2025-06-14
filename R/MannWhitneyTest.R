#' Mann-Whitney U Test (Wilcoxon Rank-Sum, Manual Implementation)
#'
#' Performs the Mann-Whitney U test (Wilcoxon rank-sum) for two independent groups,
#' using a manual implementation. Suitable when the assumptions of parametric tests
#' (normality, homogeneity of variances) are not met.
#'
#' Advantages:
#' - Does not assume normality.
#' - More powerful than t-test for skewed distributions.
#'
#' Disadvantages:
#' - Only compares two groups at a time.
#' - Sensitive to unequal variances or shapes.
#'
#' This implementation allows one- or two-sided alternatives and optionally applies a continuity correction.
#'
#' @param grupo1 Numeric vector for the first group.
#' @param grupo2 Numeric vector for the second group.
#' @param alpha Significance level (default = 0.05).
#' @param alternative Character string specifying the alternative hypothesis.
#'        Options are \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#' @param continuity Logical indicating whether to apply continuity correction (default = TRUE).
#'
#' @return An object of class \code{"comparacion"} and \code{"mannwhitney"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with the comparison name, difference in means, p-value, and significance.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: A character vector of group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A string describing the test and hypothesis direction.
#' }
#'
#' @references Mann, H. B., & Whitney, D. R. (1947). "On a Test of Whether One of Two Random Variables is Stochastically Larger than the Other." \emph{Annals of Mathematical Statistics}, 18(1), 50–60.
#'
#' @examples
#' data(d_e, package = "Analitica")
#' g1 <- d_e$Sueldo_actual[d_e$labor == 1]
#' g2 <- d_e$Sueldo_actual[d_e$labor == 2]
#' resultado <- MWTest(g1, g2, alternative = "greater")
#' summary(resultado)
#'
#' @export
#' @importFrom stats pnorm
#'
MWTest <- function(grupo1, grupo2, alpha = 0.05,
                   alternative = c("two.sided", "less", "greater"),
                   continuity = TRUE) {

  alternative <- match.arg(alternative)

  if (!is.numeric(grupo1) || !is.numeric(grupo2)) {
    stop("Both groups must be numeric vectors.")
  }

  m <- length(grupo1)
  n <- length(grupo2)

  if (m < 2 || n < 2) {
    stop("Each group must contain at least two observations.")
  }

  datos <- c(grupo1, grupo2)
  rangos <- rank(datos)

  R1 <- sum(rangos[1:m])
  U1 <- R1 - m * (m + 1) / 2
  U <- U1
  mu_U <- m * n / 2
  sigma_U <- sqrt(m * n * (m + n + 1) / 12)

  z0 <- switch(alternative,
               "two.sided" = {
                 z <- (U - mu_U - if (continuity) 0.5 * sign(U - mu_U) else 0) / sigma_U
                 p <- 2 * (1 - pnorm(abs(z)))
                 p
               },
               "greater" = {
                 z <- (U - mu_U - if (continuity) 0.5 else 0) / sigma_U
                 1 - pnorm(z)
               },
               "less" = {
                 z <- (U - mu_U + if (continuity) 0.5 else 0) / sigma_U
                 pnorm(z)
               }
  )

  p_val <- z0

  nombres <- c("Grupo1", "Grupo2")
  medias <- c(mean(grupo1), mean(grupo2))
  names(medias) <- nombres
  orden <- order(medias, decreasing = TRUE)
  etiquetas_ordenadas <- names(medias)[orden]
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))

  resultados <- data.frame(
    Comparacion = paste(nombres[1], nombres[2], sep = " - "),
    Diferencia = round(abs(diff(medias)), 4),
    Valor_Critico = NA,
    p_value = round(p_val, 4),
    Significancia = sig,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = etiquetas_ordenadas,
    Metodo = paste("Mann-Whitney U (", alternative, ", manual)", sep = "")
  )
  class(out) <- c("comparacion", "mannwhitney")
  return(out)
}
