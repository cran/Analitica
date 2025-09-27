#' Tukey HSD Test for Multiple Comparisons v2.0
#'
#' Performs Tukey's Honest Significant Difference (HSD) test for all pairwise
#' comparisons after fitting an ANOVA model. This post hoc method uses the
#' studentized range distribution and is appropriate when variances are equal
#' across groups and observations are independent.
#'
#' Tukey's test controls the family-wise error rate and is widely used when group
#' comparisons have not been planned in advance.
#'
#' Advantages:
#' - Strong control of Type I error rate.
#' - Ideal for balanced designs with equal variances.
#'
#' Disadvantages:
#' - Assumes equal variances and sample sizes.
#' - Less powerful with heteroscedasticity.
#'
#' @references Tukey, J. W. (1949). "Comparing individual means in the analysis of variance." \emph{Biometrics}, 5(2), 99–114. <https://doi.org/10.2307/3001913>
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{"tukey"} and \code{"comparaciones"} containing:
#' \itemize{
#'   \item \code{Resultados}: a data.frame with columns \code{Comparacion}, \code{Diferencia}, \code{SE}, \code{t_value},
#'         \code{p_value} (unadjusted), \code{p_ajustada} (Tukey), \code{Valor_Critico} (critical difference), and \code{Significancia}.
#'   \item \code{Promedios}: a named vector of group means as defined by \code{comparar}.
#'   \item \code{Orden_Medias}: group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: "Tukey test".
#'   \item \code{Termino}: the term being compared (e.g., "A", "B", or "A:B").
#'   \item \code{MSerror}, \code{df_error}, \code{N}: useful for plots with error bars.
#' }
#'
#' @export
#' @importFrom stats qtukey ptukey deviance qf
#' @importFrom utils combn
#'
#' @examples
#'
#' #Caso DCA
#' data(d_e, package = "Analitica")
#' mod1 <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' summary(mod1)
#' resultado <- TukeyTest(mod1)
#' summary(resultado)
#' plot(resultado)
#'
#' #Caso DBA
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' summary(mod2)
#' # Comparar niveles de 'tratamiento' (ajustando el error por el modelo con bloque)
#' res <- TukeyTest(mod2, comparar = "as.factor(labor)")
#' summary(res)
#' plot(res)
#'
#' #Caso DFA Two Ways
#' mod3 <- aov(Sueldo_actual ~ as.factor(labor) * Sexo, data = d_e)
#' summary(mod3)
#' # promedios de as.factor(labor) (promediando sobre B)
#' resA <- TukeyTest(mod3, comparar = "as.factor(labor)")
#' summary(resA)
#' plot(resA)
#'
#' # promedios de la interaccion entre factor A y factor B
#' resB <- TukeyTest(mod3, comparar = c("as.factor(labor)","Sexo"))
#' summary(resB)
#' plot(resB)
#'
TukeyTest <- function(modelo, comparar = NULL, alpha = 0.05) {

  if (is.null(modelo$model)) {
    stop("The 'model' object must contain the data (try aov/lm with embedded data).")
  }

  # Detecta factores disponibles en el modelo
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("The model has no factors in 'xlevels'. Make sure that the categorical variables are 'factor'.")
  }

  # Si no se especifica, toma el primer factor del modelo
  if (is.null(comparar)) {
    comparar <- names(xlv)[1]
  }

  # Acepta 1+ nombres: efecto principal o interaccion
  comparar <- as.character(comparar)
  mf <- modelo$model
  resp_name <- names(mf)[1]
  respuesta <- mf[[1]]

  # Verifica que existan y sean factores; de lo contrario, fuerza a factor
  for (nm in comparar) {
    if (!nm %in% names(mf)) {
      stop(sprintf("The term '%s' is not in the model data.", nm))
    }
    if (!is.factor(mf[[nm]])) {
      mf[[nm]] <- factor(mf[[nm]])
    }
  }

  # Construye el factor de grupos a comparar (principal o interaccion)
  if (length(comparar) == 1) {
    grupos <- mf[[comparar]]
    term_label <- comparar
  } else {
    grupos <- interaction(mf[, comparar, drop = FALSE], drop = TRUE)
    term_label <- paste(comparar, collapse = ":")
  }

  # Medias y tamaños por grupo
  medias <- tapply(respuesta, grupos, mean)
  n <- tapply(respuesta, grupos, length)
  nombres_grupos <- names(medias)

  # Orden para reporte
  orden_medias <- order(medias, decreasing = TRUE)
  etiquetas_ordenadas <- nombres_grupos[orden_medias]

  # MS de error desde el modelo completo (usa bloqueos y demas para reducir el error)
  df_error <- modelo$df.residual
  MSerror <- deviance(modelo) / df_error
  ng <- length(medias)
  if (ng < 2) stop("At least 2 levels/groups are needed in the term to be compared.")

  # Todas las comparaciones pareadas
  comparaciones <- combn(nombres_grupos, 2, simplify = FALSE)

  resultados <- data.frame(
    Comparacion   = character(),
    Diferencia    = numeric(),
    Valor_Critico = numeric(),
    p_value       = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  # Valor critico de Tukey (q)
  q_crit <- qtukey(1 - alpha, ng, df_error)

  for (par in comparaciones) {
    g1 <- par[1]; g2 <- par[2]

    dif <- abs(medias[g1] - medias[g2])
    # Tukey-Kramer para tamaños desiguales
    SE <- sqrt(MSerror * (1 / n[g1] + 1 / n[g2]) / 2) * sqrt(2)  # equivale a sqrt(MS*(1/ni+1/nj))

    # valor critico en escala de diferencias (HSD)
    valor_critico <- q_crit * sqrt(MSerror/2) * sqrt(1/n[g1] + 1/n[g2])

    # estadistico observado en escala q (studentized range)
    q_obs <- dif / (sqrt(MSerror/2) * sqrt(1/n[g1] + 1/n[g2]))
    p_val <- 1 - ptukey(q_obs, ng, df_error)

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01,  "**",
                         ifelse(p_val < 0.05,   "*", "ns")))

    comparacion <- paste(sort(c(g1, g2)), collapse = " - ")

    resultados <- rbind(resultados, data.frame(
      Comparacion   = comparacion,
      Diferencia    = round(dif, 4),
      Valor_Critico = round(valor_critico, 4),
      p_value       = round(p_val, 4),
      Significancia = sig,
      stringsAsFactors = FALSE
    ))
  }

  out <- list(
    Resultados   = resultados,
    Promedios    = medias,
    Orden_Medias = etiquetas_ordenadas,
    Metodo       = "Tukey",
    Termino      = term_label
  )
  class(out) <- c("comparaciones", "tukey")
  return(out)
}
