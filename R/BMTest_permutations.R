#' Brunner-Munzel Test (Permutation Version) for Two Independent Groups
#'
#' Performs the Brunner-Munzel test using a permutation approach, suitable for
#' comparing two independent samples when the assumption of equal variances may not hold.
#'
#' This version computes an empirical p-value based on resampling, without relying
#' on the t-distribution approximation.
#'
#' @param grupo1 A numeric vector representing the first group.
#' @param grupo2 A numeric vector representing the second group.
#' @param alpha Significance level (default is 0.05).
#' @param alternative Character string specifying the alternative hypothesis:
#'        one of \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#' @param nperm Number of permutations to perform (default = 10000).
#' @param seed Optional random seed for reproducibility (default is NULL).
#'
#' @return An object of class \code{"comparacion"} and \code{"brunnermunzel_perm"}, containing:
#' \itemize{
#'   \item \code{Resultados}: A data frame with comparison name, mean difference, empirical p-value, and significance.
#'   \item \code{Promedios}: A named numeric vector of group means.
#'   \item \code{Orden_Medias}: Group names ordered by their mean.
#'   \item \code{Metodo}: Description of the method used.
#' }
#'
#' @references Brunner, E., & Munzel, U. (2000). "The nonparametric Behrens-Fisher problem:
#' Asymptotic theory and a small-sample approximation." \emph{Biometrical Journal}, 42(1), 17–25.
#'
#' @examples
#' data(d_e, package = "Analitica")
#' g1 <- d_e$Sueldo_actual[d_e$labor == 1]
#' g2 <- d_e$Sueldo_actual[d_e$labor == 2]
#' resultado <- BMpTest(g1, g2)
#' summary(resultado)
#'
#' @export
BMpTest <- function(grupo1, grupo2, alpha = 0.05,
                    alternative = c("two.sided", "less", "greater"),
                    nperm = 10000, seed=NULL) {
  if(!is.null(seed)){
    set.seed(seed)
  }

  alternative <- match.arg(alternative)

  if (!is.numeric(grupo1) || !is.numeric(grupo2)) {
    stop("Both groups must be numeric vectors.")
  }

  if (length(grupo1) < 2 || length(grupo2) < 2) {
    stop("Each group must contain at least two observations.")
  }

  m <- length(grupo1)
  n <- length(grupo2)
  datos <- c(grupo1, grupo2)
  total <- m + n

  # Internal function to compute the Brunner-Munzel statistic
  calc_stat <- function(x, y) {
    datos <- c(x, y)
    ranks <- rank(datos)
    R1 <- ranks[1:length(x)]
    R2 <- ranks[(length(x) + 1):length(datos)]

    R1_mean <- mean(R1)
    R2_mean <- mean(R2)

    S1_sq <- sum((R1 - R1_mean)^2) / (length(x) - 1)
    S2_sq <- sum((R2 - R2_mean)^2) / (length(y) - 1)
    S_sq <- (S1_sq / length(x)) + (S2_sq / length(y))
    T_stat <- (R1_mean - R2_mean) / sqrt(S_sq)

    return(T_stat)
  }

  # Observed statistic
  T_obs <- calc_stat(grupo1, grupo2)

  # Permutations
  T_perm <- numeric(nperm)
  for (i in seq_len(nperm)) {
    idx <- sample(total)
    x_perm <- datos[idx[1:m]]
    y_perm <- datos[idx[(m + 1):total]]
    T_perm[i] <- calc_stat(x_perm, y_perm)
  }

  # Empirical p-value
  p_val <- switch(alternative,
                  "two.sided" = mean(abs(T_perm) >= abs(T_obs)),
                  "greater"   = mean(T_perm >= T_obs),
                  "less"      = mean(T_perm <= T_obs)
  )

  # Calcular P(X<Y) + 0.5*P(X=Y)
  comparaciones <- outer(grupo1, grupo2, FUN = function(xi, yj) {
    ifelse(xi < yj, 1, ifelse(xi == yj, 0.5, 0))
  })
  p_hat <- mean(comparaciones)

  # Output formatting
  nombres <- c("Grupo1", "Grupo2")
  medias <- c(mean(grupo1), mean(grupo2))
  names(medias) <- nombres
  orden <- order(medias, decreasing = TRUE)
  etiquetas_ordenadas <- names(medias)[orden]
  diferencia <- abs(diff(medias))
  sig <- ifelse(p_val < alpha, "*", "ns")

  resultados <- data.frame(
    Comparacion = paste(nombres[1], nombres[2], sep = " - "),
    Diferencia = round(diferencia, 4),
    Valor_Critico = NA,
    p_value = round(p_val, 4),
    p_hat = round(p_hat, 4),
    Significancia = sig,
    stringsAsFactors = FALSE
  )

  out <- list(
    Resultados = resultados,
    Promedios = medias,
    Orden_Medias = etiquetas_ordenadas,
    Metodo = paste("Brunner-Munzel (perm, ", alternative, ")", sep = "")
  )
  class(out) <- c("comparacion", "brunnermunzel_perm")
  return(out)
}
