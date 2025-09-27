#' Holm-Adjusted Pairwise Comparisons v2.0
#'
#' Performs pairwise t-tests with p-values adjusted using Holm’s sequential method.
#'
#' Advantages:
#' - Controls family-wise error rate more efficiently than Bonferroni.
#' - Easy to apply over any set of p-values.
#'
#' Disadvantages:
#' - Does not adjust test statistics, only p-values.
#' - Slightly more conservative than false discovery rate (FDR) methods.
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{"holm"} and \code{"comparaciones"} containing:
#' \itemize{
#'   \item \code{Resultados}: a data.frame with columns \code{Comparacion}, \code{Diferencia}, \code{SE}, \code{t_value},
#'         \code{p_value} (unadjusted), \code{p_ajustada} (Holm), \code{Valor_Critico} (critical difference), and \code{Significancia}.
#'   \item \code{Promedios}: a named vector of group means as defined by \code{comparar}.
#'   \item \code{Orden_Medias}: group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: "Holm t-test".
#'   \item \code{Termino}: the term being compared (e.g., "A", "B", or "A:B").
#'   \item \code{MSerror}, \code{df_error}, \code{N}: useful for plots with error bars.
#' }
#'
#' @references
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.Scandinavian Journal of Statistics, 6(2), 65–70.
#'
#' @export
#' @importFrom stats pt p.adjust deviance
#' @importFrom utils combn
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- HolmTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' # RCBD
#' mod <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' res <- HolmTest(mod, comparar = "as.factor(labor)")
#' summary(res); plot(res)                # usa p_ajustada automaticamente
#'
#' # Factorial
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) * Sexo, data = d_e)
#' resAB <- HolmTest(mod2, comparar = c("as.factor(labor)","Sexo"))
#' summary(resAB, n = Inf); plot(resAB, horizontal = TRUE)
#'
HolmTest <- function(modelo, comparar = NULL, alpha = 0.05) {
  if (is.null(modelo$model)) {
    stop("The 'model' object must contain the data (use aov/lm with embedded data).")
  }

  # Factores disponibles
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("The model has no factors in 'xlevels'. It converts categoricals to 'factor'.")
  }

  # Si no se especifica, usa el primer factor
  if (is.null(comparar)) {
    comparar <- names(xlv)[1]
  }
  comparar <- as.character(comparar)

  mf <- modelo$model
  respuesta <- mf[[1]]

  # Verifica/convierte a factor los terminos a comparar
  for (nm in comparar) {
    if (!nm %in% names(mf)) stop(sprintf("The term '%s' is not in the model data.", nm))
    if (!is.factor(mf[[nm]])) mf[[nm]] <- factor(mf[[nm]])
  }

  # Factor de grupos (principal o interaccion)
  if (length(comparar) == 1) {
    grupos <- mf[[comparar]]
    term_label <- comparar
  } else {
    grupos <- interaction(mf[, comparar, drop = FALSE], drop = TRUE)
    term_label <- paste(comparar, collapse = ":")
  }

  # Medias y tamaños
  medias <- tapply(respuesta, grupos, mean)
  n      <- tapply(respuesta, grupos, length)
  if (length(medias) < 2) stop("At least two levels are needed in the term to be compared.")
  nombres_grupos <- names(medias)

  # MS de error y g.l. del modelo completo
  df_error <- modelo$df.residual
  MSerror  <- deviance(modelo) / df_error

  # Comparaciones pareadas
  pares <- combn(nombres_grupos, 2, simplify = FALSE)
  m <- length(pares)

  # Prealocar
  Comparacion <- character(m)
  Diferencia  <- numeric(m)
  SE          <- numeric(m)
  t_value     <- numeric(m)
  p_value     <- numeric(m)

  for (i in seq_along(pares)) {
    g1 <- pares[[i]][1]; g2 <- pares[[i]][2]
    dif <- abs(medias[g1] - medias[g2])
    se  <- sqrt(MSerror * (1 / n[g1] + 1 / n[g2]))
    t   <- dif / se
    p   <- 2 * stats::pt(-abs(t), df_error)

    Comparacion[i] <- paste(sort(c(g1, g2)), collapse = " - ")
    Diferencia[i]  <- dif
    SE[i]          <- se
    t_value[i]     <- t
    p_value[i]     <- p
  }

  p_ajustada <- p.adjust(p_value, method = "holm")
  Sig <- ifelse(p_ajustada < 0.001, "***",
                ifelse(p_ajustada < 0.01,  "**",
                       ifelse(p_ajustada < 0.05,   "*", "ns")))

  resultados <- data.frame(
    Comparacion = Comparacion,
    Diferencia  = round(Diferencia, 4),
    SE          = round(SE, 4),
    t_value     = round(t_value, 4),
    p_value     = round(p_value, 4),
    p_ajustada  = round(p_ajustada, 4),
    Significancia = Sig,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados   = resultados,
    Promedios    = medias,
    Orden_Medias = names(sort(medias, decreasing = TRUE)),
    Metodo       = "Holm-adjusted t-test",
    Termino      = term_label,
    MSerror      = MSerror,
    df_error     = df_error,
    N            = n
  )
  class(out) <- c("comparaciones", "holm")
  return(out)
}
