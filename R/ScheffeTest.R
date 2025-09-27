#' Scheffe Test for Multiple Comparisons v2.0
#'
#' Performs Scheffe's post hoc test after fitting an ANOVA model. This test compares all possible
#' pairs of group means, using a critical value based on the F-distribution.
#'
#' The Scheffe test is a conservative method, making it harder to detect significant differences,
#' but reducing the likelihood of Type I errors (false positives). It is especially appropriate
#' when the comparisons were not pre-planned and the number of contrasts is large.
#'
#' Assumptions: normally distributed residuals and homogeneity of variances.
#'
#' Advantages:
#' - Very robust to violations of assumptions.
#' - Suitable for complex comparisons, not just pairwise.
#'
#' Disadvantages:
#' - Very conservative; reduced power.
#' - Not ideal for detecting small differences.
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default 0.05).
#'
#' @return Objeto de clase \code{"scheffe"} and \code{"comparaciones"} with:
#' \itemize{
#'   \item \code{Resultados}: data.frame with \code{Comparacion}, \code{Diferencia},
#'         \code{SE2} (= MSerror*(1/n_i+1/n_j)), \code{F_obs}, \code{Valor_Critico}, \code{p_value}, \code{Significancia}.
#'   \item \code{Promedios}, \code{Orden_Medias}, \code{Metodo}="Scheffe", \code{Termino}.
#'   \item \code{MSerror}, \code{df_error}, \code{N} (utiles para \code{plot.comparaciones()}).
#' }
#'
#' @references
#' Scheffe, H. (1953). "A method for judging all contrasts in the analysis of variance." \emph{Biometrika}, 40(1/2), 87–104. <https://doi.org/10.1093/biomet/40.1-2.87>
#'
#' @importFrom  stats qf pf deviance
#' @importFrom utils combn
#' @export
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- ScheffeTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' # RCBD
#' mod <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' res <- ScheffeTest(mod, comparar = "as.factor(labor)")
#' summary(res); plot(res)                      # plot usara p_value
#'
#' # Factorial
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) * Sexo, data = d_e)
#' resAB <- ScheffeTest(mod2, comparar = c("as.factor(labor)","Sexo"))
#' summary(resAB, n = Inf); plot(resAB, horizontal = TRUE)
#'
#'
ScheffeTest <- function(modelo, comparar = NULL, alpha = 0.05) {
  if (is.null(modelo$model)) {
    stop("The 'model' object must contain the data (use aov/lm with embedded data)..")
  }

  # Factores disponibles
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("The model has no factors in 'xlevels'. Convert the categoricals to 'factor'.")
  }

  # Si no se especifica, usa el primer factor
  if (is.null(comparar)) comparar <- names(xlv)[1]
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
  if (length(medias) < 2) stop("At least two levels are needed in the term to be compared.")
  nombres_grupos <- names(medias)

  # MS de error y g.l.
  df_error <- modelo$df.residual
  MSerror  <- deviance(modelo) / df_error
  g        <- length(medias)           # numero de grupos del termino
  df1      <- g - 1
  Fcrit    <- stats::qf(1 - alpha, df1, df_error)

  # Todas las comparaciones
  pares <- utils::combn(nombres_grupos, 2, simplify = FALSE)

  # Prealocar
  m <- length(pares)
  Comparacion   <- character(m)
  Diferencia    <- numeric(m)
  SE2           <- numeric(m)  # MSerror*(1/ni+1/nj)
  F_obs         <- numeric(m)
  Valor_Critico <- numeric(m)  # umbral de diferencia |yi-yj|
  p_value       <- numeric(m)
  Significancia <- character(m)

  for (i in seq_along(pares)) {
    g1 <- pares[[i]][1]; g2 <- pares[[i]][2]
    dif <- abs(medias[g1] - medias[g2])
    se2 <- MSerror * (1 / n[g1] + 1 / n[g2])

    # Fobs de Scheffe para pares: ((dif)^2) / ((g-1)*MSerror*(1/ni+1/nj))
    Fobs <- (dif^2) / ((g - 1) * se2)
    pval <- 1 - stats::pf(Fobs, df1, df_error)

    # Diferencia critica (umbral en escala de medias)
    diff_crit <- sqrt((g - 1) * Fcrit * se2)

    sig <- ifelse(pval < 0.001, "***",
                  ifelse(pval < 0.01,  "**",
                         ifelse(pval < 0.05,   "*", "ns")))

    Comparacion[i]   <- paste(sort(c(g1, g2)), collapse = " - ")
    Diferencia[i]    <- round(dif, 4)
    SE2[i]           <- round(se2, 6)
    F_obs[i]         <- round(Fobs, 4)
    Valor_Critico[i] <- round(diff_crit, 4)
    p_value[i]       <- round(pval, 4)
    Significancia[i] <- sig
  }

  resultados <- data.frame(
    Comparacion   = Comparacion,
    Diferencia    = Diferencia,
    SE2           = SE2,
    F_obs         = F_obs,
    Valor_Critico = Valor_Critico,
    p_value       = p_value,
    Significancia = Significancia,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados   = resultados,
    Promedios    = medias,
    Orden_Medias = names(sort(medias, decreasing = TRUE)),
    Metodo       = "Scheffe",
    Termino      = term_label,
    MSerror      = MSerror,
    df_error     = df_error,
    N            = n
  )
  class(out) <- c("comparaciones", "scheffe")
  return(out)
}
