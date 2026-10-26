#' Dunnett's T3 Post Hoc Test v2.2.0
#'
#' Performs Dunnett's T3 test for pairwise comparisons after an ANOVA model.
#' This test is recommended when group variances are unequal and sample sizes differ.
#' It is based on the studentized range distribution and provides conservative control
#' over Type I error without assuming homoscedasticity.
#'
#' Advantages:
#' - More powerful than T2 when group sizes are small.
#' - Adjusted for unequal variances.
#'
#' Disadvantages:
#' - Complex critical value estimation.
#' - Less frequently used and harder to find in software.
#'
#' @references Dunnett, C. W. (1980). "Pairwise multiple comparisons in the unequal variance case."
#' \emph{Journal of the American Statistical Association}, 75(372), 796–800.
#' <https://doi.org/10.1080/01621459.1980.10477558>
#'
#' @param modelo An object of class \code{aov} or \code{lm} (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: interaction (e.g., \code{c("A","B")} for \code{A:B})
#'   If omitted, it uses the first factor in \code{modelo$xlevels}.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"dunnettt3"} and \code{"comparaciones"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with pairwise comparisons, mean differences,
#'         \code{q_value}, \code{gl}, \code{p_value}, and significance indicators.
#'   \item \code{Promedios}: A named numeric vector of group means as defined by \code{comparar}.
#'   \item \code{Orden_Medias}: A character vector of group names ordered from highest to lowest mean.
#'   \item \code{Metodo}: A character string with the test name ("Dunnett T3").
#'   \item \code{Termino}: The term being compared (e.g., "A", "B", or "A:B").
#' }
#'
#' @export
#' @importFrom stats ptukey qf
#' @importFrom utils combn
#'
#' @examples
#' data(d_e, package = "Analitica")
#' mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' resultado <- T3Test(mod)
#' summary(resultado)
#' plot(resultado)
#'
#' # Con bloques, comparando solo el factor de interés
#' mod2 <- aov(Sueldo_actual ~ as.factor(labor) + Sexo, data = d_e)
#' res2 <- T3Test(mod2, comparar = "as.factor(labor)")
#' summary(res2)
#' plot(res2)
#'
#' # Modelo con interacción
#' mod3 <- aov(Sueldo_actual ~ as.factor(labor) * Sexo, data = d_e)
#' res3 <- T3Test(mod3, comparar = c("as.factor(labor)", "Sexo"))
#' summary(res3)
#' plot(res3)
#'
T3Test <- function(modelo, comparar = NULL, alpha = 0.05) {

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

  comparaciones <- .pares_ordenados(medias)

  resultados <- data.frame(
    Comparacion   = character(),
    Diferencia    = numeric(),
    q_value       = numeric(),
    gl            = numeric(),
    p_value       = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  .g1 <- .g2 <- character(0); .p <- numeric(0)

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

    # q observado (escala de rango studentizado)
    q_val <- dif / se_ij * sqrt(2)

    # p-value desde la distribución de rango studentizado
    p_val <- 1 - stats::ptukey(q_val, nmeans = length(nombres_grupos), df = gl)

    sig <- .sig_stars(p_val, alpha)

    comparacion <- paste(g1, "-", g2)

    .g1 <- c(.g1, g1); .g2 <- c(.g2, g2); .p <- c(.p, unname(p_val))

    resultados <- rbind(
      resultados,
      data.frame(
        Comparacion   = comparacion,
        Diferencia    = round(dif, 4),
        q_value       = round(q_val, 4),
        gl            = round(gl, 2),
        p_value       = round(p_val, 4),
        Significancia = sig,
        stringsAsFactors = FALSE
      )
    )
  }

  out <- list(
    Resultados   = resultados,
    Promedios    = medias,
    Orden_Medias = etiquetas_ordenadas,
    Metodo       = "Dunnett T3",
    Termino      = term_label
  )
  out <- .agregar_info_pares(out, .g1, .g2, .p)
  class(out) <- c("dunnettt3", "comparaciones")

  return(out)
}
