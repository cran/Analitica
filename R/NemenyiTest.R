#' Nemenyi Test for Multiple Comparisons (Non-Parametric)
#'
#' Performs the Nemenyi (Nemenyi-Damico-Wolfe) all-pairs test after a significant
#' Kruskal-Wallis test, using mean ranks computed on the pooled sample and the
#' studentized range distribution.
#'
#' For each pair of groups the statistic is
#' \deqn{q_{ij} = \frac{|\bar R_i - \bar R_j|}{\sqrt{\frac{N(N+1)}{12} \, C \, \frac{1}{2}\left(\frac{1}{n_i} + \frac{1}{n_j}\right)}}}
#' where \eqn{C = 1 - \sum (t^3 - t)/(N^3 - N)} is the tie correction. The p-value is
#' obtained from the studentized range distribution with \eqn{k} means and infinite
#' degrees of freedom, so it already controls the family-wise error rate.
#' With unequal group sizes the Tukey-Kramer extension is used.
#'
#' Advantages:
#' - Controls the family-wise error rate without an additional p-value adjustment.
#' - Direct non-parametric analogue of Tukey's HSD.
#'
#' Disadvantages:
#' - Exact only for equal group sizes; with unequal sizes it is an approximation (Tukey-Kramer).
#' - Usually less powerful than Conover-Iman.
#' - It is only useful in completely random or single-factor designs
#'   (for blocked designs use a Friedman-type procedure).
#'
#' @param formula A formula of the form \code{y ~ group}.
#' @param data A data frame containing the variables.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"nemenyi"}, \code{"comparaciones_np"} and \code{"comparaciones"}, including:
#' \itemize{
#'   \item \code{Resultados}: Data frame with comparisons (higher - lower mean rank), mean-rank
#'         difference, q statistic, critical difference, p-values and significance codes.
#'   \item \code{Promedios}: Mean ranks of each group, sorted from highest to lowest.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest mean rank.
#'   \item \code{Metodo}: Name of the method ("Nemenyi (non-parametric)").
#'   \item \code{Datos}: Data frame with the response and the group, used by \code{plot()}.
#'   \item \code{KruskalWallis}: Result of \code{kruskal.test()}.
#' }
#'
#' @references
#' Nemenyi, P. (1963). \emph{Distribution-free Multiple Comparisons}. PhD thesis, Princeton University.
#'
#' Hollander, M., Wolfe, D. A. & Chicken, E. (2014). \emph{Nonparametric Statistical Methods} (3rd ed.). Wiley.
#'
#' @examples
#' set.seed(123)
#' datos <- data.frame(
#'   grupo = rep(c("A", "B", "C", "D"), each = 10),
#'   medida = c(rnorm(10, 10), rnorm(10, 12), rnorm(10, 15), rnorm(10, 11))
#' )
#' resultado <- NemenyiTest(medida ~ grupo, data = datos)
#' summary(resultado)
#' plot(resultado)
#'
#' @export
#' @importFrom stats ptukey qtukey kruskal.test model.frame
#' @importFrom utils combn
NemenyiTest <- function(formula, data, alpha = 0.05) {
  mf <- model.frame(formula, data)
  respuesta <- mf[[1]]
  grupo <- droplevels(as.factor(mf[[2]]))

  k <- nlevels(grupo)
  if (k < 2) stop("At least 2 groups are required.")
  N <- length(respuesta)
  n <- table(grupo)

  ranks <- rank(respuesta)
  Rj <- tapply(ranks, grupo, mean)
  C <- .factor_empates(respuesta)
  S2 <- N * (N + 1) / 12 * C

  pares <- .pares_ordenados(Rj)
  m <- length(pares)
  g1 <- vapply(pares, `[`, character(1), 1)
  g2 <- vapply(pares, `[`, character(1), 2)

  q_crit <- stats::qtukey(1 - alpha, k, Inf)
  dif <- as.numeric(Rj[g1] - Rj[g2])
  se <- sqrt(S2 / 2 * (1 / as.numeric(n[g1]) + 1 / as.numeric(n[g2])))
  q_obs <- dif / se
  p_val <- stats::ptukey(q_obs, k, Inf, lower.tail = FALSE)

  resultados <- data.frame(
    Comparacion   = paste(g1, "-", g2),
    Diferencia    = round(dif, 4),
    q_value       = round(q_obs, 4),
    Valor_Critico = round(q_crit * se, 4),
    p_value       = round(p_val, 4),
    p_ajustada    = round(p_val, 4),
    Significancia = .sig_stars(p_val, alpha),
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados    = resultados,
    Promedios     = Rj[names(sort(Rj, decreasing = TRUE))],
    Orden_Medias  = names(sort(Rj, decreasing = TRUE)),
    Metodo        = "Nemenyi (non-parametric)",
    Datos         = data.frame(y = respuesta, grupo = grupo),
    KruskalWallis = stats::kruskal.test(respuesta, grupo)
  )
  out <- .agregar_info_pares(out, g1, g2, p_val, p_val)
  class(out) <- c("nemenyi", "comparaciones_np", "comparaciones")
  out
}
