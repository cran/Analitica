#' Games-Howell Post Hoc Test v.2.2.0
#'
#' Performs the Games-Howell test for pairwise comparisons after ANOVA,
#' without assuming equal variances or sample sizes. It is suitable when
#' Levene or Bartlett test indicates heterogeneity of variances.
#'
#' Advantages:
#' - Excellent for heteroscedastic data.
#' - Controls Type I error across unequal group sizes.
#'
#' Disadvantages:
#' - Slightly conservative in small samples.
#' - More complex to compute than Tukey.
#'
#' @references Games, P. A., & Howell, J. F. (1976).
#' "Pairwise Multiple Comparison Procedures with Unequal N's and/or Variances:
#' A Monte Carlo Study". Journal of Educational Statistics, 1(2), 113–125.
#' <https://doi.org/10.1002/j.2162-6057.1976.tb00211.x>
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"gameshowell"} and \code{"comparaciones"},
#' which contains:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, including:
#'         \code{Comparacion}, \code{Diferencia}, \code{t_value}, \code{gl},
#'         \code{p_value}, and \code{Significancia}.
#'   \item \code{Promedios}: A named numeric vector of group means as defined by \code{comparar}.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string indicating the method used ("Games-Howell").
#'   \item \code{Termino}: The term being compared (e.g., "A", "B", or "A:B").
#' }
#'
#' @export
#'
#' @importFrom stats pt
#' @importFrom utils combn
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' # Comparación sobre el primer factor del modelo
#' resultado <- GHTest(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' # Con bloques, comparando solo el factor de interés
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' res2 <- GHTest(mod2, comparar = "as.factor(labor)")
#' summary(res2)
#' plot(res2)
#'
#' # Modelo con interacción
#' mod3 <- aov(Sueldo_actual ~ as.factor(labor) * Sexo, data = d_e)
#' # efecto principal
#' resA <- GHTest(mod3, comparar = "as.factor(labor)")
#' # interacción
#' resAB <- GHTest(mod3, comparar = c("as.factor(labor)", "Sexo"))
#' summary(resAB)
#' plot(resAB)
GHTest <- function(modelo, comparar = NULL, alpha = 0.05) {

  if (is.null(modelo$model)) {
    stop("The 'model' object must contain the data (try aov/lm with embedded data).")
  }

  # Factores disponibles en el modelo
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("The model has no factors in 'xlevels'. Make sure that the categorical variables are 'factor'.")
  }

  # Si no se especifica, toma el primer factor del modelo
  if (is.null(comparar)) {
    comparar <- names(xlv)[1]
  }

  comparar <- as.character(comparar)
  mf <- modelo$model
  resp_name <- names(mf)[1]
  respuesta <- mf[[1]]

  # Verifica que existan y sean factores; si no, los fuerza a factor
  for (nm in comparar) {
    if (!nm %in% names(mf)) {
      stop(sprintf("The term '%s' is not in the model data.", nm))
    }
    if (!is.factor(mf[[nm]])) {
      mf[[nm]] <- factor(mf[[nm]])
    }
  }

  # Construye el factor de grupos a comparar (principal o interacción)
  if (length(comparar) == 1) {
    grupos <- mf[[comparar]]
    term_label <- comparar
  } else {
    grupos <- interaction(mf[, comparar, drop = FALSE], drop = TRUE)
    term_label <- paste(comparar, collapse = ":")
  }

  # Medias, tamaños y varianzas por grupo
  medias <- tapply(respuesta, grupos, mean)
  n     <- tapply(respuesta, grupos, length)
  s2    <- tapply(respuesta, grupos, var)
  nombres_grupos <- names(medias)

  if (length(nombres_grupos) < 2) {
    stop("At least 2 levels/groups are needed in the term to be compared.")
  }

  # Ordenar nombres de los grupos por media (de mayor a menor)
  orden_medias <- order(medias, decreasing = TRUE)
  etiquetas_ordenadas <- nombres_grupos[orden_medias]

  comparaciones <- combn(nombres_grupos, 2, simplify = FALSE)

  resultados <- data.frame(
    Comparacion = character(),
    Diferencia  = numeric(),
    t_value     = numeric(),
    gl          = numeric(),
    p_value     = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  for (par in comparaciones) {
    g1 <- par[1]
    g2 <- par[2]

    dif   <- abs(medias[g1] - medias[g2])
    se_ij <- sqrt((s2[g1] / n[g1]) + (s2[g2] / n[g2]))

    # Welch-Satterthwaite degrees of freedom
    df_num <- (s2[g1] / n[g1] + s2[g2] / n[g2])^2
    df_den <- ((s2[g1]^2) / (n[g1]^2 * (n[g1] - 1))) +
      ((s2[g2]^2) / (n[g2]^2 * (n[g2] - 1)))
    gl <- df_num / df_den

    t_val <- dif / se_ij
    p_val <- 2 * stats::pt(-abs(t_val), df = gl)  # two-tailed

    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01,  "**",
                         ifelse(p_val < 0.05,   "*", "ns")))

    # Ordenar nombres para coincidir con plot.comparacion()
    nombres_ordenados <- sort(c(g1, g2))
    comparacion <- paste(nombres_ordenados, collapse = " - ")

    resultados <- rbind(
      resultados,
      data.frame(
        Comparacion = comparacion,
        Diferencia  = round(dif, 4),
        t_value     = round(t_val, 4),
        gl          = round(gl, 2),
        p_value     = round(p_val, 4),
        Significancia = sig,
        stringsAsFactors = FALSE
      )
    )
  }

  out <- list(
    Resultados   = resultados,
    Promedios    = medias,
    Orden_Medias = etiquetas_ordenadas,
    Metodo       = "Games-Howell",
    Termino      = term_label
  )
  class(out) <- c("comparaciones", "gameshowell")

  return(out)
}
