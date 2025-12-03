#' Dwass-Steel-Critchlow-Fligner (DSCF) Test (Non-Parametric) v2.0.3
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
#' - It is only useful in completely random or single-factor designs.
#'
#' @references Dwass, M. (1960). Some k-sample rank-order tests. In I. Olkin et al. (Eds.), Contribution1s to Probability and Statistics: Essays in Honor of Harold Hotelling (pp. 198 - 202). Stanford University Press.
#'
#'@examples
#'data(d_e, package = "Analitica")
#'DSCFTest(Sueldo_actual ~ labor, data = d_e)
#'
#' @param formula  y ~ group
#' @param data     data.frame con las variables
#' @param alpha    nivel (0.05 por defecto)  just for the little star
#' @param method.p adjustment method (default "holm")
#' @param na.rm    remove NA (TRUE by default)
#' @param include_kw if TRUE, add summary of Kruskal-Wallis test
#'
#' @return objeto con clases c("comparaciones","dscf")
#'
#' @export
#' @importFrom stats pnorm p.adjust kruskal.test model.frame
#' @importFrom utils combn
DSCFTest <- function(formula, data, alpha = 0.05, method.p = "holm",
                     na.rm = TRUE, include_kw = TRUE) {
  mf <- model.frame(formula, data)
  y <- mf[[1]]
  g <- as.factor(mf[[2]])

  if (na.rm) {
    keep <- !(is.na(y) | is.na(g))
    y <- y[keep]; g <- g[keep]
  }
  if (!is.numeric(y)) stop("The answer must be numerical.")
  if (nlevels(g) < 2L) stop("at least 2 groups are required.")

  niveles <- levels(g)
  k <- length(niveles)
  N <- length(y)

  ## Rangos globales (promedio en empates)
  r <- rank(y, ties.method = "average")

  ## Medias de rangos y tamanos por grupo
  split_r <- split(r, g)
  n_j <- vapply(split_r, length, integer(1))
  Rbar_j <- vapply(split_r, mean, numeric(1))

  ## Varianza pooled de rangos dentro de grupos:
  ## S^2 = sum_j sum_i (R_ij - Rbar_j)^2 / (N - k)
  ss_within <- sum(vapply(names(split_r), function(nm) {
    rx <- split_r[[nm]]
    sum( (rx - mean(rx))^2 )
  }, numeric(1)))
  S2 <- ss_within / (N - k)

  ## Correccion por empates (Critchlow–Fligner):
  ## T = 1 - sum(t^3 - t)/(N^3 - N)    donde t = tamanos de cada empate en y
  tab <- table(y)
  tie_term <- if (any(tab > 1L)) sum(tab^3 - tab) else 0
  T_corr <- 1 - tie_term / (N^3 - N)
  if (T_corr <= 0) warning("Degenerate tie correction; check the data.")
  S2_t <- S2 / T_corr

  ## Todas las comparaciones por pares
  pares <- utils::combn(niveles, 2, simplify = FALSE)

  z_vals <- numeric(length(pares))
  p_vals <- numeric(length(pares))
  comp_lab <- character(length(pares))

  for (i in seq_along(pares)) {
    g1 <- pares[[i]][1]; g2 <- pares[[i]][2]
    dif <- abs(Rbar_j[g1] - Rbar_j[g2])
    se  <- sqrt(S2_t * (1 / n_j[g1] + 1 / n_j[g2]))
    z   <- dif / se
    z_vals[i] <- z
    p_vals[i] <- 2 * (1 - pnorm(z))
    comp_lab[i] <- paste(g1, g2, sep = " - ")
  }

  p_adj <- p.adjust(p_vals, method = method.p)
  sig <- ifelse(p_adj < 0.001, "***",
                ifelse(p_adj < 0.01,  "**",
                       ifelse(p_adj < 0.05,  "*",  "ns")))

  resultados <- data.frame(
    Comparacion   = comp_lab,
    z             = round(z_vals, 4),
    p_value       = round(p_vals, 6),
    p_ajustada    = round(p_adj, 6),
    Significancia = sig,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados    = resultados[order(p_adj, z_vals, decreasing = FALSE), ],
    Promedios     = Rbar_j,
    Orden_Medias  = names(sort(Rbar_j, decreasing = TRUE)),
    N_por_grupo   = n_j,
    Parametros    = list(N = N, k = k, S2 = S2, T_corr = T_corr, S2_t = S2_t),
    Metodo        = "DSCF (in the parametric)",
    Ajuste        = method.p
  )

  if (isTRUE(include_kw)) {
    out$KruskalWallis <- stats::kruskal.test(y ~ g)
  }

  class(out) <- c("comparaciones", "dscf")
  out
}
