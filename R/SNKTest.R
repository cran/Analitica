#' Student-Newman-Keuls (SNK) Test for Multiple Comparisons v2.0
#'
#' Performs the Student-Newman-Keuls (SNK) post hoc test for pairwise comparisons
#' after fitting an ANOVA model. The test uses a stepwise approach where the
#' critical value depends on the number of means spanned between groups (range r).
#'
#' SNK is more powerful but less conservative than Tukey’s HSD, increasing the chance of
#' detecting real differences while slightly raising the Type I error rate.
#'
#' Assumptions: normality, homogeneity of variances, and independence of observations.
#'
#' Advantages:
#' - More powerful than Tukey when differences are large.
#' - Intermediate control of Type I error.
#'
#' Disadvantages:
#' - Error control is not family-wise.
#' - Type I error increases with more comparisons.
#'
#' @references
#' Student, Newman, and Keuls (1952). "Student-Newman-Keuls Procedure". See also: <https://doi.org/10.1002/bimj.200310019>
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{"SNK"} and \code{"comparaciones"} containing:
#' \itemize{
#'   \item \code{Resultados}: a data.frame with columns \code{Comparacion}, \code{Diferencia}, \code{SE}, \code{t_value},
#'         \code{p_value} (unadjusted), \code{p_ajustada} (SNK), \code{Valor_Critico} (critical difference), and \code{Significancia}.
#'   \item \code{Promedios}: a named vector of group means as defined by \code{comparar}.
#'   \item \code{Orden_Medias}: group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: "SNK t-test".
#'   \item \code{Termino}: the term being compared (e.g., "A", "B", or "A:B").
#'   \item \code{MSerror}, \code{df_error}, \code{N}: useful for plots with error bars.
#' }
#'
#' @export
#' @importFrom stats ptukey qtukey deviance qf
#' @importFrom utils combn
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- SNKTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' # RCBD
#' mod <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' res <- SNKTest(mod, comparar = "as.factor(labor)")
#' summary(res); plot(res)                      # plot usara p_value
#'
#' # Factorial
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) * Sexo, data = d_e)
#' resAB <- SNKTest(mod2, comparar = c("as.factor(labor)","Sexo"))
#' summary(resAB, n = Inf); plot(resAB, horizontal = TRUE)
#'
SNKTest <- function(modelo, comparar = NULL, alpha = 0.05) {
  if (is.null(modelo$model)) {
    stop("The 'model' object must contain the data (use aov/lm with embedded data).")
  }

  # Factores disponibles
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("The model has no factors in 'xlevels'. Convert categoricals to 'factor'.")
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

  # Grupos (efecto principal o interaccion)
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
  if (length(medias) < 2) stop("At least two levels are required in the term to be compared.")
  nombres_grupos <- names(medias)

  # MS de error y g.l. residuales del modelo completo
  df_error <- modelo$df.residual
  MSerror  <- deviance(modelo) / df_error

  # Error estandar por par (Tukey-Kramer para tamaños desiguales)
  SE_ij <- function(g1, g2) sqrt(MSerror * (1 / n[g1] + 1 / n[g2]))

  # Todas las comparaciones
  pares <- combn(nombres_grupos, 2, simplify = FALSE)

  # Orden de medias ascendente para calcular r (número de medias abarcadas)
  ordenado_asc <- names(sort(medias, decreasing = FALSE))

  # Prealocar
  m <- length(pares)
  Comparacion   <- character(m)
  Diferencia    <- numeric(m)
  SE            <- numeric(m)
  Valor_Critico <- numeric(m)
  p_value       <- numeric(m)
  Significancia <- character(m)

  for (i in seq_along(pares)) {
    g1 <- pares[[i]][1]
    g2 <- pares[[i]][2]

    # r = |pos2 - pos1| + 1 en el orden de medias ascendente
    pos1 <- match(g1, ordenado_asc)
    pos2 <- match(g2, ordenado_asc)
    r    <- abs(pos2 - pos1) + 1

    dif    <- abs(medias[g1] - medias[g2])
    se_val <- SE_ij(g1, g2)

    # Critico y p con rango studentizado y 'r' grados de rango
    q_crit       <- stats::qtukey(1 - alpha, r, df_error)
    valor_crit   <- q_crit * se_val / sqrt(2)

    q_obs <- dif * sqrt(2) / se_val
    p_val <- 1 - stats::ptukey(q_obs, r, df_error)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01,  "**",
                         ifelse(p_val < 0.05,   "*", "ns")))

    Comparacion[i]   <- paste(sort(c(g1, g2)), collapse = " - ")
    Diferencia[i]    <- round(dif, 4)
    SE[i]            <- round(se_val, 4)
    Valor_Critico[i] <- round(valor_crit, 4)
    p_value[i]       <- round(p_val, 4)
    Significancia[i] <- sig
  }

  resultados <- data.frame(
    Comparacion   = Comparacion,
    Diferencia    = Diferencia,
    SE            = SE,
    Valor_Critico = Valor_Critico,
    p_value       = p_value,
    Significancia = Significancia,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados   = resultados,
    Promedios    = medias,
    Orden_Medias = names(sort(medias, decreasing = TRUE)),
    Metodo       = "SNK",
    Termino      = term_label,
    MSerror      = MSerror,
    df_error     = df_error,
    N            = n
  )
  class(out) <- c("comparaciones", "snk")
  return(out)
}
