#' Dwass-Steel-Critchlow-Fligner (DSCF) Test (Non-Parametric)
#'
#' Robust non-parametric method for multiple comparisons after Kruskal-Wallis.
#' Uses rank-based pairwise tests with a pooled variance estimate.
#'
#' Advantages:
#' - Strong control of Type I error with unequal sample sizes.
#' - More powerful than Dunn in many conditions.
#'
#' Disadvantages:
#' - Computationally more complex.
#' - Less commonly available in standard software.
#'
#' @param formula A formula of the form \code{y ~ group}.
#' @param data A data frame containing the variables.
#' @param alpha Significance level (default is 0.05).
#' @param method.p Method for p-value adjustment (default is "holm").
#'
#' @return An object of class \code{"dscf"} and \code{"comparaciones"}, including:
#' \itemize{
#'   \item \code{Resultados}: Data frame with comparisons, z-statistics, p-values, adjusted p-values, and significance levels.
#'   \item \code{Promedios}: Mean ranks of each group.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest mean rank.
#'   \item \code{Metodo}: "DSCF (no paramétrico)".
#' }
#'
#' @references Dwass, M. (1960). Some k-sample rank-order tests. In I. Olkin et al. (Eds.), Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling (pp. 198–202). Stanford University Press.
#'
#'@examples
#'data(d_e, package = "Analitica")
#'DSCFTest(Sueldo_actual ~ labor, data = d_e)
#'
#'
#' @export
#' @importFrom stats pnorm p.adjust
#' @importFrom utils combn

DSCFTest <- function(formula, data, alpha = 0.05, method.p = "holm") {
  mf <- model.frame(formula, data)
  respuesta <- mf[[1]]
  grupo <- as.factor(mf[[2]])

  niveles <- levels(grupo)
  k <- length(niveles)
  N <- length(respuesta)
  ranks <- rank(respuesta)
  n <- table(grupo)
  Rj <- tapply(ranks, grupo, mean)

  comparaciones <- combn(niveles, 2, simplify = FALSE)
  resultados <- data.frame(
    Comparacion = character(),
    z = numeric(),
    p_value = numeric(),
    p_ajustada = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  S2 <- (sum(ranks^2) - N * (N + 1)^2 / 4) / (N - 1)

  z_vals <- numeric(length(comparaciones))
  p_vals <- numeric(length(comparaciones))

  for (i in seq_along(comparaciones)) {
    par <- comparaciones[[i]]
    g1 <- par[1]; g2 <- par[2]
    ni <- n[g1]; nj <- n[g2]
    Ri <- Rj[g1]; Rj_ <- Rj[g2]

    dif <- abs(Ri - Rj_)
    se <- sqrt(S2 * (1 / ni + 1 / nj))
    z <- dif / se
    p <- 2 * (1 - pnorm(z))

    z_vals[i] <- z
    p_vals[i] <- p
  }

  p_ajustada <- p.adjust(p_vals, method = method.p)
  sig <- ifelse(p_ajustada < 0.001, "***",
                ifelse(p_ajustada < 0.01, "**",
                       ifelse(p_ajustada < 0.05, "*", "ns")))

  for (i in seq_along(comparaciones)) {
    comp <- paste(sort(comparaciones[[i]]), collapse = " - ")
    resultados[i, ] <- list(
      Comparacion = comp,
      z = round(z_vals[i], 4),
      p_value = round(p_vals[i], 4),
      p_ajustada = round(p_ajustada[i], 4),
      Significancia = sig[i]
    )
  }

  out <- list(
    Resultados = resultados,
    Promedios = Rj,
    Orden_Medias = names(sort(Rj, decreasing = TRUE)),
    Metodo = "DSCF (no parametrico)"
  )
  class(out) <- c("comparaciones", "dscf")
  return(out)
}
