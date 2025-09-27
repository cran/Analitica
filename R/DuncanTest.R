#' Duncan Multiple Range Test (DMRT) v2.0
#'
#' Performs the Duncan test for pairwise comparisons after an ANOVA.
#' This method is more liberal than Tukey's HSD, using a stepwise approach
#' with critical values from the studentized range distribution.
#'
#' Advantages:
#' - High power for detecting differences.
#' - Simple to interpret and implement.
#'
#' Disadvantages:
#' - Inflates Type I error rate.
#' - Not recommended for confirmatory research.
#'
#' @references
#' Duncan, D. B. (1955). "Multiple range and multiple F tests." Biometrics, 11(1), 1-42.
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{"duncan"} and \code{"comparaciones"} containing:
#' \itemize{
#'   \item \code{Resultados}: a data.frame with columns \code{Comparacion}, \code{Diferencia}, \code{SE}, \code{t_value},
#'         \code{p_value} (unadjusted), \code{p_ajustada} (duncan), \code{Valor_Critico} (critical difference), and \code{Significancia}.
#'   \item \code{Promedios}: a named vector of group means as defined by \code{comparar}.
#'   \item \code{Orden_Medias}: group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: "Duncan t-test".
#'   \item \code{Termino}: the term being compared (e.g., "A", "B", or "A:B").
#'   \item \code{MSerror}, \code{df_error}, \code{N}: useful for plots with error bars.
#' }
#'
#' @export
#' @importFrom stats qtukey ptukey qf deviance
#' @importFrom utils combn
#' @examples
#'
#' # DCA
#' data(d_e, package = "Analitica")
#' mod1 <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- DuncanTest(mod1)
#' summary(resultado)
#' plot(resultado)
#'
#' # DBA
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' res <- DuncanTest(mod2, comparar = "as.factor(labor)")
#' summary(res); plot(res)
#'
#' # DFactorial
#' mod3 <- aov(Sueldo_actual ~as.factor(labor) * Sexo, data = d_e)
#' resAB <- DuncanTest(mod3, comparar = c("as.factor(labor)","Sexo"))  # celdas A:B
#' summary(resAB, n = Inf); plot(resAB, horizontal = TRUE)
#'
DuncanTest <- function(modelo, comparar = NULL, alpha = 0.05) {
  if (is.null(modelo$model)) {
    stop("The 'model' object must contain the data (use aov/lm with embedded data).")
  }

  # Factores disponibles en el modelo
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("The model has no factors in 'xlevels'. Be sure to convert categoricals to 'factor'.")
  }

  # Si no se especifica, toma el primer factor del modelo
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

  # Factor de grupos a comparar (principal o interaccion)
  if (length(comparar) == 1) {
    grupos <- mf[[comparar]]
    term_label <- comparar
  } else {
    grupos <- interaction(mf[, comparar, drop = FALSE], drop = TRUE)
    term_label <- paste(comparar, collapse = ":")
  }

  # Medias, tamaños y orden
  medias <- tapply(respuesta, grupos, mean)
  n      <- tapply(respuesta, grupos, length)
  if (length(medias) < 2) stop("At least two levels are required in the term to be compared.")
  nombres_grupos <- names(medias)

  # MS de error y g.l. residuales del modelo completo
  df_error <- modelo$df.residual
  MSerror  <- deviance(modelo) / df_error

  # Funcion de SE por par (permite tamaños desiguales)
  SE_ij <- function(g1, g2) sqrt(MSerror * (1 / n[g1] + 1 / n[g2]))

  # Todas las comparaciones
  pares <- combn(nombres_grupos, 2, simplify = FALSE)

  # Para DMRT se ordenan las medias (ascendente) para calcular r = |pos2 - pos1| + 1
  ordenado_asc <- names(sort(medias, decreasing = FALSE))

  resultados <- data.frame(
    Comparacion   = character(length(pares)),
    Diferencia    = numeric(length(pares)),
    Valor_Critico = numeric(length(pares)),
    p_value       = numeric(length(pares)),
    Significancia = character(length(pares)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(pares)) {
    g1 <- pares[[i]][1]
    g2 <- pares[[i]][2]

    pos1 <- which(ordenado_asc == g1)
    pos2 <- which(ordenado_asc == g2)
    r    <- abs(pos2 - pos1) + 1

    dif    <- abs(medias[g1] - medias[g2])
    se_val <- SE_ij(g1, g2)

    # Critico DMRT para rango r (studentized range)
    q_crit       <- stats::qtukey(1 - alpha, r, df_error)
    valor_crit   <- q_crit * se_val / sqrt(2)

    # Estadistico observado y p-valor (para rango r)
    q_obs <- dif * sqrt(2) / se_val
    p_val <- 1 - stats::ptukey(q_obs, r, df_error)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01,  "**",
                         ifelse(p_val < 0.05,   "*", "ns")))

    resultados$Comparacion[i]   <- paste(sort(c(g1, g2)), collapse = " - ")
    resultados$Diferencia[i]    <- round(dif, 4)
    resultados$Valor_Critico[i] <- round(valor_critico <- valor_crit, 4)
    resultados$p_value[i]       <- round(p_val, 4)
    resultados$Significancia[i] <- sig
  }

  out <- list(
    Resultados   = resultados,
    Promedios    = medias,
    Orden_Medias = names(sort(medias, decreasing = TRUE)),
    Metodo       = "Duncan",
    Termino      = term_label,
    MSerror      = MSerror,
    df_error     = df_error,
    N            = n
  )
  class(out) <- c("comparaciones", "duncan")
  return(out)
}
