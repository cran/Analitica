#' Dunn's Test for Multiple Comparisons (Non-Parametric)
#'
#' Performs Dunn's test for pairwise comparisons following a Kruskal-Wallis test.
#' Suitable for non-parametric data (ordinal or non-normal), using rank sums.
#' Includes Holm correction by default for multiple comparisons.
#'
#' For each pair, \eqn{z_{ij} = |\bar R_i - \bar R_j| / \sqrt{\frac{N(N+1)}{12} C \left(\frac{1}{n_i} + \frac{1}{n_j}\right)}},
#' where \eqn{C = 1 - \sum (t^3 - t)/(N^3 - N)} corrects for ties.
#'
#' Advantages:
#' - Simple and widely used non-parametric alternative to Tukey's test.
#' - Handles unequal sample sizes.
#' - Compatible with various p-value corrections (e.g., Holm, Bonferroni).
#'
#' Disadvantages:
#' - Less powerful than DSCF or Conover when sample sizes vary widely.
#' - Requires ranking all data and can be conservative depending on adjustment.
#'
#' @param formula A formula of the form \code{y ~ group}.
#' @param data A data frame containing the variables.
#' @param alpha Significance level (default is 0.05).
#' @param method.p Method for p-value adjustment (default is "holm").
#'
#' @return An object of class \code{"dunn"}, \code{"comparaciones_np"} and \code{"comparaciones"}, including:
#' \itemize{
#'   \item \code{Resultados}: Data frame with group comparisons (higher - lower mean rank), z-values, raw and adjusted p-values, and significance.
#'   \item \code{Promedios}: Mean ranks of each group, sorted from highest to lowest.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest rank.
#'   \item \code{Metodo}: "Dunn (non-parametric)".
#'   \item \code{Datos}: Data frame with the response and the group, used by \code{plot()}.
#'   \item \code{KruskalWallis}: Result of \code{kruskal.test()}.
#' }
#'
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics*, 6(3), 241–252. \doi{10.1080/00401706.1964.10490181}
#'
#'@examples
#'data(d_e, package = "Analitica")
#'DunnTest(Sueldo_actual ~ labor, data = d_e)
#'
#'
#' @export
#' @importFrom stats p.adjust pnorm kruskal.test model.frame
#' @importFrom utils combn

DunnTest <- function(formula, data, alpha = 0.05, method.p = "holm") {
  mf <- model.frame(formula, data)
  response <- mf[[1]]
  grupo <- droplevels(as.factor(mf[[2]]))

  k <- nlevels(grupo)
  if (k < 2) stop("At least 2 groups are required.")
  N <- length(response)
  n <- table(grupo)

  ranks <- rank(response)
  Rj <- tapply(ranks, grupo, mean)
  S2 <- N * (N + 1) / 12 * .factor_empates(response)

  pares <- .pares_ordenados(Rj)
  g1 <- vapply(pares, `[`, character(1), 1)
  g2 <- vapply(pares, `[`, character(1), 2)

  dif <- as.numeric(Rj[g1] - Rj[g2])
  se  <- sqrt(S2 * (1 / as.numeric(n[g1]) + 1 / as.numeric(n[g2])))
  z_vals <- dif / se
  p_vals <- 2 * stats::pnorm(-abs(z_vals))
  p_ajustada <- stats::p.adjust(p_vals, method = method.p)

  resultados <- data.frame(
    Comparacion   = paste(g1, "-", g2),
    Diferencia    = round(dif, 4),
    z             = round(z_vals, 4),
    p_value       = round(p_vals, 4),
    p_ajustada    = round(p_ajustada, 4),
    Significancia = .sig_stars(p_ajustada, alpha),
    stringsAsFactors = FALSE
  )

  orden <- names(sort(Rj, decreasing = TRUE))
  out <- list(
    Resultados    = resultados,
    Promedios     = Rj[orden],
    Orden_Medias  = orden,
    Metodo        = "Dunn (non-parametric)",
    Ajuste        = method.p,
    Datos         = data.frame(y = response, grupo = grupo),
    KruskalWallis = stats::kruskal.test(response, grupo)
  )
  out <- .agregar_info_pares(out, g1, g2, p_vals, p_ajustada)
  class(out) <- c("dunn", "comparaciones_np", "comparaciones")
  out
}
