#' Nemenyi Test for Multiple Comparisons (Non-Parametric)
#'
#' Performs the Nemenyi test after a significant Kruskal-Wallis or Friedman test.
#' Based on the studentized range distribution applied to mean ranks.
#'
#' Advantages:
#' - Easy to implement for equal-sized groups.
#' - Conservative control of family-wise error rate.
#'
#' Disadvantages:
#' - Only valid with equal group sizes.
#' - No p-values are directly calculated (based on critical differences only).
#'
#' @param formula A formula of the form \code{y ~ group}.
#' @param data A data frame containing the variables.
#' @param alpha Significance level (default is 0.05).
#'
#' @return An object of class \code{"nemenyi"} and \code{"comparaciones"}, including:
#' \itemize{
#'   \item \code{Resultados}: Data frame with group comparisons, rank differences, critical value, p-values, and significance codes.
#'   \item \code{Promedios}: Mean ranks of each group.
#'   \item \code{Orden_Medias}: Group names ordered from highest to lowest rank.
#'   \item \code{Metodo}: Name of the method ("Nemenyi (no paramétrico)").
#' }
#'
#' @references Nemenyi, P. (1963). Distribution-free Multiple Comparisons.
#'
#'@examples
#'set.seed(123)
#'datos <- data.frame(
#'  grupo = rep(c("A", "B", "C", "D"), each = 10),
#'  medida = c(
#'    rnorm(10, mean = 10),
#'    rnorm(10, mean = 12),
#'    rnorm(10, mean = 15),
#'    rnorm(10, mean = 11)
#'  )
#')
# Verificar que los tamaños de grupo sean iguales (requisito de Nemenyi)
#'table(datos$grupo)
#'#> A  B  C  D
#'#>10 10 10 10
#'# Aplicar el test de Nemenyi
#'resultado <- NemenyiTest(medida ~ grupo, data = datos)
#'# Ver los resultados
#'summary(resultado)
#'# O simplemente
#'resultado$Resultados
#'# Ver orden de medias (rangos)
#'resultado$Orden_Medias
#'
#'
#'
#' @export
#' @importFrom stats qtukey
#' @importFrom utils combn

NemenyiTest <- function(formula, data, alpha = 0.05) {
  mf <- model.frame(formula, data)
  respuesta <- mf[[1]]
  grupo <- as.factor(mf[[2]])

  niveles <- levels(grupo)
  k <- length(niveles)
  n_grupo <- table(grupo)

  if (length(unique(n_grupo)) != 1) {
    stop("El test de Nemenyi requiere n iguales por grupo.")
  }

  n <- unique(n_grupo)[1]
  N <- length(respuesta)
  ranks <- rank(respuesta)
  Rj <- tapply(ranks, grupo, mean)

  comparaciones <- combn(niveles, 2, simplify = FALSE)
  resultados <- data.frame(
    Comparacion = character(),
    Diferencia = numeric(),
    Valor_Critico = numeric(),
    p_value = numeric(),
    Significancia = character(),
    stringsAsFactors = FALSE
  )

  SE <- sqrt(k * (k + 1) / (6 * n))
  q_crit <- qtukey(1 - alpha, k, Inf) / sqrt(2)
  diff_crit <- q_crit * SE

  for (par in comparaciones) {
    g1 <- par[1]
    g2 <- par[2]

    dif <- abs(Rj[g1] - Rj[g2])
    z <- dif / SE
    pval <- 2 * (1 - pnorm(z))

    sig <- ifelse(dif > diff_crit,
                  ifelse(dif > diff_crit * 1.5, "**", "*"),
                  "ns")

    resultados <- rbind(resultados, data.frame(
      Comparacion = paste(sort(c(g1, g2)), collapse = " - "),
      Diferencia = round(dif, 4),
      Valor_Critico = round(diff_crit, 4),
      p_value = round(pval, 4),
      Significancia = sig
    ))
  }

  out <- list(
    Resultados = resultados,
    Promedios = Rj,
    Orden_Medias = names(sort(Rj, decreasing = TRUE)),
    Metodo = "Nemenyi (no parametrico)"
  )
  class(out) <- c("comparaciones", "nemenyi")

  return(out)
}
