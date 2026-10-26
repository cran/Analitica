#' Conover-Iman Test for Multiple Comparisons (Non-Parametric)
#'
#' Performs non-parametric pairwise comparisons based on rank-transformed data using
#' the Conover-Iman procedure. This method is typically applied as a post hoc test
#' following a significant Kruskal-Wallis test to identify specific group differences.
#'
#' The Conover-Iman test uses rank-based t-statistics, offering improved statistical
#' power over Dunn's test while maintaining flexibility in sample size.
#'
#' Ranks are computed on the pooled sample (mid-ranks for ties). For each pair,
#' \deqn{t_{ij} = \frac{|\bar R_i - \bar R_j|}{\sqrt{S^2 \, \frac{N - 1 - H}{N - k} \left(\frac{1}{n_i} + \frac{1}{n_j}\right)}}}
#' with \eqn{S^2 = \frac{1}{N-1}\left(\sum R^2 - \frac{N(N+1)^2}{4}\right)} and \eqn{H} the
#' (tie-corrected) Kruskal-Wallis statistic. The statistic is compared with a
#' t distribution with \eqn{N - k} degrees of freedom.
#'
#' Advantages:
#' - More powerful than Dunn’s test, especially with moderate group differences.
#' - Robust to non-normal data and suitable for ordinal or skewed distributions.
#' - Allows for unequal sample sizes across groups.
#'
#' Disadvantages:
#' - Sensitive to heteroscedasticity (non-constant variances).
#' - Requires appropriate p-value adjustment to control the family-wise error rate.
#' - It is only useful in completely random or single-factor designs.
#'
#' @param formula A formula of the form \code{y ~ group}, where \code{y} is a numeric variable
#'        and \code{group} is a factor indicating group membership.
#' @param data A data frame containing the variables specified in the formula.
#' @param alpha Significance level for hypothesis testing (default is 0.05).
#' @param method.p Method used to adjust p-values for multiple comparisons (default is \code{"holm"}).
#'
#' @return An object of class \code{"conover"}, \code{"comparaciones_np"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons (higher - lower mean rank), t-statistics, raw and adjusted p-values, and significance markers.
#'   \item \code{Promedios}: A named numeric vector with mean ranks for each group, sorted from highest to lowest.
#'   \item \code{Orden_Medias}: A character vector with group names sorted from highest to lowest rank.
#'   \item \code{Metodo}: A string describing the method used ("Conover-Iman (non-parametric)").
#'   \item \code{Datos}: Data frame with the response and the group, used by \code{plot()}.
#'   \item \code{KruskalWallis}: Result of \code{kruskal.test()}.
#' }
#'
#' @references Conover, W. J. & Iman, R. L. (1979). "Multiple comparisons using rank sums." \emph{Technometrics}, 21(4), 489–495.
#'
#' @examples
#'data(d_e, package = "Analitica")
#'ConoverTest(Sueldo_actual ~ labor, data = d_e)
#'
#'
#' @export
#' @importFrom stats pt p.adjust kruskal.test model.frame
#' @importFrom utils combn

ConoverTest <- function(formula, data, alpha = 0.05, method.p = "holm") {
  mf <- model.frame(formula, data)
  respuesta <- mf[[1]]
  grupo <- droplevels(as.factor(mf[[2]]))

  k <- nlevels(grupo)
  if (k < 2) stop("At least 2 groups are required.")
  N <- length(respuesta)
  if (N - k < 1) stop("Not enough observations (N - k must be at least 1).")
  n <- table(grupo)

  ranks <- rank(respuesta)
  Rj <- tapply(ranks, grupo, mean)

  # Varianza de los rangos (incluye empates) y estadistico H de Kruskal-Wallis
  S2 <- (sum(ranks^2) - N * (N + 1)^2 / 4) / (N - 1)
  H  <- (sum(as.numeric(n) * Rj^2) - N * (N + 1)^2 / 4) / S2
  gl <- N - k

  pares <- .pares_ordenados(Rj)
  g1 <- vapply(pares, `[`, character(1), 1)
  g2 <- vapply(pares, `[`, character(1), 2)

  dif <- as.numeric(Rj[g1] - Rj[g2])
  se  <- sqrt(S2 * (N - 1 - H) / gl * (1 / as.numeric(n[g1]) + 1 / as.numeric(n[g2])))
  t_vals <- dif / se
  p_vals <- 2 * stats::pt(-abs(t_vals), gl)
  p_ajustada <- stats::p.adjust(p_vals, method = method.p)

  resultados <- data.frame(
    Comparacion   = paste(g1, "-", g2),
    Diferencia    = round(dif, 4),
    t_value       = round(t_vals, 4),
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
    Metodo        = "Conover-Iman (non-parametric)",
    Ajuste        = method.p,
    gl            = gl,
    Datos         = data.frame(y = respuesta, grupo = grupo),
    KruskalWallis = stats::kruskal.test(respuesta, grupo)
  )
  out <- .agregar_info_pares(out, g1, g2, p_vals, p_ajustada)
  class(out) <- c("conover", "comparaciones_np", "comparaciones")
  out
}
