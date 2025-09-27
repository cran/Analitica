#' Gabriel’s Post Hoc Test for Multiple Comparisons v 2.0
#'
#' A modification of Tukey's test for use with moderately unequal sample sizes.
#'
#' Advantages:
#' - More powerful than Tukey for unequal group sizes.
#' - Controls error rates effectively with moderate imbalance.
#'
#' Disadvantages:
#' - Can be anti-conservative with large differences in group sizes.
#' - Less common in standard statistical software.
#'
#' @references
#' Hochberg, Y., & Tamhane, A. C. (1987). Multiple Comparison Procedures.
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{"gabriel"} and \code{"comparaciones"} containing:
#' \itemize{
#'   \item \code{Resultados}: a data.frame with columns \code{Comparacion}, \code{Diferencia}, \code{SE}, \code{t_value},
#'         \code{p_value} (unadjusted), \code{p_ajustada} (gabriel), \code{Valor_Critico} (critical difference), and \code{Significancia}.
#'   \item \code{Promedios}: a named vector of group means as defined by \code{comparar}.
#'   \item \code{Orden_Medias}: group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: "Gabriel t-test".
#'   \item \code{Termino}: the term being compared (e.g., "A", "B", or "A:B").
#'   \item \code{MSerror}, \code{df_error}, \code{N}: useful for plots with error bars.
#' }
#'
#' @export
#' @importFrom stats qtukey ptukey deviance
#' @importFrom utils combn
#'
#' @examples
#'
#' # DCA
#' data(d_e, package = "Analitica")
#' mod1 <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- GabrielTest(mod1)
#' summary(resultado)
#' plot(resultado)
#'
#' # RCBD
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' res <- GabrielTest(mod2, comparar = "as.factor(labor)")
#' summary(res); plot(res)
#'
#' # Factorial
#' mod3 <- aov(Sueldo_actual ~ as.factor(labor) * Sexo, data = d_e)
#' resAB <- GabrielTest(mod3, comparar = c("as.factor(labor)","Sexo"))  # celdas A:B
#' summary(resAB, n = Inf); plot(resAB, horizontal = TRUE)
#'
GabrielTest <- function(modelo, comparar = NULL, alpha = 0.05) {
  if (is.null(modelo$model)) {
    stop("The 'model' object must contain the data (use aov/lm with embedded data)..")
  }

  # Factores disponibles
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("The model has no factors in 'xlevels'. Convert categoricals to 'factor'.")
  }

  # Si no se especifica, toma el primer factor
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
  if (length(medias) < 2) stop("At least two levels are required in the term to be compared.")
  nombres_grupos <- names(medias)

  # MS de error y g.l. del modelo completo
  df_error <- modelo$df.residual
  MSerror  <- deviance(modelo) / df_error
  ng <- length(medias)

  # Todas las comparaciones
  pares <- combn(nombres_grupos, 2, simplify = FALSE)

  # Prealocar resultados
  m <- length(pares)
  Comparacion   <- character(m)
  Diferencia    <- numeric(m)
  SE            <- numeric(m)
  Valor_Critico <- numeric(m)
  p_value       <- numeric(m)
  Significancia <- character(m)

  for (i in seq_along(pares)) {
    g1 <- pares[[i]][1]; g2 <- pares[[i]][2]

    dif <- abs(medias[g1] - medias[g2])

    # SE Gabriel: sqrt(MSE * 2 / (n_i + n_j))
    ni <- n[g1]; nj <- n[g2]
    se_ij <- sqrt(MSerror * (2 / (ni + nj)))

    # Critico y p usando rango studentizado con ng y df_error
    q_crit       <- stats::qtukey(1 - alpha, ng, df_error)
    valor_crit   <- q_crit * se_ij / sqrt(2)

    q_obs <- dif * sqrt(2) / se_ij
    p_val <- 1 - stats::ptukey(q_obs, ng, df_error)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01,  "**",
                         ifelse(p_val < 0.05,   "*", "ns")))

    Comparacion[i]   <- paste(sort(c(g1, g2)), collapse = " - ")
    Diferencia[i]    <- round(dif, 4)
    SE[i]            <- round(se_ij, 4)
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
    Metodo       = "Gabriel",
    Termino      = term_label,
    MSerror      = MSerror,
    df_error     = df_error,
    N            = n
  )
  class(out) <- c("comparaciones", "gabriel")
  return(out)
}
