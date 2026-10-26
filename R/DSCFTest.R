#' Dwass-Steel-Critchlow-Fligner (DSCF) Test (Non-Parametric)
#'
#' All-pairs non-parametric multiple comparison procedure, usually applied after
#' a Kruskal-Wallis test. Unlike Dunn or Conover-Iman, each pair of groups is
#' ranked separately (only the observations of those two groups), and the
#' standardized Wilcoxon statistic is referred to the studentized range
#' distribution with \eqn{k} groups and infinite degrees of freedom.
#'
#' For groups \eqn{i} and \eqn{j}, with \eqn{W_i} the sum of the ranks of group
#' \eqn{i} in the joint ranking of both groups,
#' \deqn{W^*_{ij} = \sqrt{2}\,\frac{W_i - n_i(n_i + n_j + 1)/2}{\sqrt{V_{ij}}}}
#' \deqn{V_{ij} = \frac{n_i n_j}{12}\left(n_i + n_j + 1 - \frac{\sum (t^3 - t)}{(n_i + n_j)(n_i + n_j - 1)}\right)}
#' where \eqn{t} are the sizes of the tie groups within the pair. The p-value is
#' \eqn{P(Q_{k,\infty} \ge |W^*_{ij}|)}, which already controls the family-wise error rate.
#'
#' Advantages:
#' - Strong control of the family-wise error rate with unequal sample sizes.
#' - Each comparison only uses the data of the two groups involved.
#'
#' Disadvantages:
#' - Computationally more complex.
#' - Less commonly available in standard software.
#' - It is only useful in completely random or single-factor designs.
#'
#' @references
#' Dwass, M. (1960). Some k-sample rank-order tests. In I. Olkin et al. (Eds.), \emph{Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling} (pp. 198-202). Stanford University Press.
#'
#' Steel, R. G. D. (1960). A rank sum test for comparing all pairs of treatments. \emph{Technometrics}, 2(2), 197-207. \doi{10.1080/00401706.1960.10489894}
#'
#' Critchlow, D. E. & Fligner, M. A. (1991). On distribution-free multiple comparisons in the one-way analysis of variance. \emph{Communications in Statistics - Theory and Methods}, 20(1), 127-139. \doi{10.1080/03610929108830487}
#'
#' @examples
#' data(d_e, package = "Analitica")
#' DSCFTest(Sueldo_actual ~ labor, data = d_e)
#'
#' @param formula  A formula of the form \code{y ~ group}.
#' @param data     A data frame containing the variables.
#' @param alpha    Significance level (default 0.05).
#' @param method.p Ignored. Kept for backward compatibility: the studentized range
#'   distribution already controls the family-wise error rate.
#' @param na.rm    Remove NA (TRUE by default).
#' @param include_kw If TRUE, add the Kruskal-Wallis test to the output.
#'
#' @return An object of class \code{"dscf"}, \code{"comparaciones_np"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: Data frame with comparisons (higher - lower mean rank), W statistic,
#'         p-values and significance codes.
#'   \item \code{Promedios}: Mean ranks of each group (pooled ranking), sorted from highest to lowest.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest mean rank.
#'   \item \code{N_por_grupo}: Sample size of each group.
#'   \item \code{Metodo}: "DSCF (non-parametric)".
#'   \item \code{Datos}: Data frame with the response and the group, used by \code{plot()}.
#'   \item \code{KruskalWallis}: Result of \code{kruskal.test()} (if \code{include_kw = TRUE}).
#' }
#'
#' @export
#' @importFrom stats ptukey kruskal.test model.frame
#' @importFrom utils combn
DSCFTest <- function(formula, data, alpha = 0.05, method.p = NULL,
                     na.rm = TRUE, include_kw = TRUE) {
  if (!is.null(method.p)) {
    warning("'method.p' is ignored: DSCF p-values come from the studentized range ",
            "distribution and already control the family-wise error rate.")
  }

  mf <- model.frame(formula, data, na.action = if (na.rm) stats::na.omit else stats::na.pass)
  y <- mf[[1]]
  g <- droplevels(as.factor(mf[[2]]))

  if (!is.numeric(y)) stop("The response must be numeric.")
  k <- nlevels(g)
  if (k < 2L) stop("At least 2 groups are required.")

  n_j <- table(g)
  # Rangos sobre la muestra completa: solo para ordenar los grupos y graficar
  Rbar_j <- tapply(rank(y), g, mean)

  pares <- .pares_ordenados(Rbar_j)
  g1 <- vapply(pares, `[`, character(1), 1)
  g2 <- vapply(pares, `[`, character(1), 2)

  W_vals <- numeric(length(pares))
  for (i in seq_along(pares)) {
    yi <- y[g == g1[i]]; yj <- y[g == g2[i]]
    ni <- length(yi); nj <- length(yj); nn <- ni + nj
    r  <- rank(c(yi, yj))                 # rangos solo dentro del par
    Wi <- sum(r[seq_len(ni)])
    tt <- table(c(yi, yj))
    V  <- ni * nj / 12 * (nn + 1 - sum(tt^3 - tt) / (nn * (nn - 1)))
    W_vals[i] <- sqrt(2) * (Wi - ni * (nn + 1) / 2) / sqrt(V)
  }
  p_vals <- stats::ptukey(abs(W_vals), k, Inf, lower.tail = FALSE)

  resultados <- data.frame(
    Comparacion   = paste(g1, "-", g2),
    W             = round(W_vals, 4),
    p_value       = round(p_vals, 4),
    p_ajustada    = round(p_vals, 4),
    Significancia = .sig_stars(p_vals, alpha),
    stringsAsFactors = FALSE
  )

  orden <- names(sort(Rbar_j, decreasing = TRUE))
  out <- list(
    Resultados    = resultados,
    Promedios     = Rbar_j[orden],
    Orden_Medias  = orden,
    N_por_grupo   = n_j,
    Metodo        = "DSCF (non-parametric)",
    Datos         = data.frame(y = y, grupo = g)
  )
  if (isTRUE(include_kw)) out$KruskalWallis <- stats::kruskal.test(y, g)

  out <- .agregar_info_pares(out, g1, g2, p_vals, p_vals)
  class(out) <- c("dscf", "comparaciones_np", "comparaciones")
  out
}
