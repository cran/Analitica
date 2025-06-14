#' Brunner-Munzel Test for Two Independent Samples
#'
#' Performs the Brunner-Munzel nonparametric test for two independent groups,
#' which estimates the probability that a randomly selected value from one group
#' is less than a randomly selected value from the other group.
#'
#' This test is suitable when group variances are unequal and/or sample sizes differ.
#' It does not assume equal variances and is often used as a more robust alternative to the Wilcoxon test.
#'
#' Advantages:
#' - Handles unequal variances and non-normality.
#' - Recommended when variance homogeneity is questionable.
#'
#' Disadvantages:
#' - Less well-known and supported.
#' - Requires large sample sizes for accurate inference.
#'
#' @param grupo1 Numeric vector of values from group 1.
#' @param grupo2 Numeric vector of values from group 2.
#' @param alpha Significance level (default = 0.05).
#' @param alternative Character string specifying the alternative hypothesis.
#'        One of \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#'
#' @return An object of class \code{"comparacion"} and \code{"brunnermunzel"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with test statistics, p-value, and estimated effect size.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: Group names ordered by their mean values (descending).
#'   \item \code{Metodo}: A character string describing the test and hypothesis.
#'   \item \code{p_hat}: Estimated probability that a value from grupo1 is less than a value from grupo2 (plus 0.5 * ties).
#' }
#'
#' @references Brunner, E., & Munzel, U. (2000). "The nonparametric Behrens-Fisher problem: Asymptotic theory and a small-sample approximation." \emph{Biometrical Journal}, 42(1), 17–25. <https://doi.org/10.1002/(SICI)1521-4036(200001)42:1%3C17::AID-BIMJ17%3E3.0.CO;2-U>
#'
#' @examples
#' data(d_e, package = "Analitica")
#' g1 <- d_e$Sueldo_actual[d_e$labor == 1]
#' g2 <- d_e$Sueldo_actual[d_e$labor == 2]
#' resultado <- BMTest(g1, g2, alternative = "greater")
#' summary(resultado)
#'
#' @export
#' @importFrom stats pt qt
BMTest <- function(grupo1, grupo2, alpha = 0.05, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)

  # Validations
  if (!is.numeric(grupo1) || !is.numeric(grupo2)) {
    stop("Both groups must be numeric vectors.")
  }
  if (length(grupo1) < 2 || length(grupo2) < 2) {
    stop("Each group must have at least two observations.")
  }

  m <- length(grupo1)
  n <- length(grupo2)
  datos <- c(grupo1, grupo2)
  rangos <- rank(datos)

  R1 <- rangos[1:m]
  R2 <- rangos[(m + 1):(m + n)]

  R1_mean <- mean(R1)
  R2_mean <- mean(R2)

  S1_sq <- sum((R1 - R1_mean)^2) / (m - 1)
  S2_sq <- sum((R2 - R2_mean)^2) / (n - 1)

  S_sq <- (S1_sq / m) + (S2_sq / n)
  SE <- sqrt(S_sq)
  T_stat <- (R1_mean - R2_mean) / SE

  df_num <- S_sq^2
  df_den <- ((S1_sq^2) / (m^2 * (m - 1))) + ((S2_sq^2) / (n^2 * (n - 1)))
  df <- df_num / df_den

  p_val <- switch(alternative,
                  "two.sided" = 2 * (1 - pt(abs(T_stat), df)),
                  "greater"   = 1 - pt(T_stat, df),
                  "less"      = pt(T_stat, df))

  t_crit <- qt(1 - alpha / 2, df)
  # Optional confidence interval (commented out for now)
  # CI_lower <- (R1_mean - R2_mean) - t_crit * SE
  # CI_upper <- (R1_mean - R2_mean) + t_crit * SE

  # Effect size: P(X < Y) + 0.5 * P(X = Y)
  count_less <- sum(outer(grupo1, grupo2, "<"))
  count_equal <- sum(outer(grupo1, grupo2, "=="))
  total_pairs <- m * n
  p_hat <- (count_less + 0.5 * count_equal) / total_pairs

  # Results and ordering
  nombres <- c("Grupo1", "Grupo2")
  medias <- c(mean(grupo1), mean(grupo2))
  names(medias) <- nombres
  orden <- order(medias, decreasing = TRUE)
  etiquetas_ordenadas <- names(medias)[orden]
  diferencia <- medias[1] - medias[2]
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))

  resultados <- data.frame(
    Comparacion = paste(nombres[1], nombres[2], sep = " - "),
    Diferencia = round(diferencia, 4),
    df = round(df, 2),
    SE = round(SE, 4),
    t_critical = round(t_crit, 4),
    p_value = round(p_val, 4),
    p_hat = round(p_hat, 4),
    Significancia = sig,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    df = df,
    Orden_Medias = etiquetas_ordenadas,
    Metodo = paste("Brunner-Munzel (", alternative, ")", sep = ""),
    p_hat = p_hat
  )
  class(out) <- c("comparacion", "brunnermunzel")

  return(out)
}
