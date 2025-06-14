#' Summary Method for Objects of Class 'comparacion'
#'
#' Displays a formatted summary of the results from a pairwise comparison test
#' of two independent groups. Compatible with objects returned by functions like
#' \code{BMTest()} or \code{MWTest()}.
#'
#' @param object An object of class \code{"comparacion"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns a one-row data frame with the summary statistics.
#' @export
summary.comparacion <- function(object, ...) {
  if (!inherits(object, "comparacion")) {
    stop("The object must be of class 'comparacion'.")
  }

  if (is.null(object$Resultados) || is.null(object$Promedios) || is.null(object$Orden_Medias)) {
    stop("The 'comparacion' object is incomplete or malformed.")
  }

  resultados <- object$Resultados[1, , drop = FALSE]
  promedios <- object$Promedios
  orden <- object$Orden_Medias
  metodo <- if (!is.null(object$Metodo)) object$Metodo else "comparacion"

  cat("========================================\n")
  cat("  Summary of Pairwise Comparison Test\n")
  cat("========================================\n")
  cat("Method:", metodo, "\n\n")

  # Funcion auxiliar para mostrar columnas con nombres alternativos
  mostrar_campo <- function(nombres_posibles, etiqueta) {
    nombre_real <- intersect(nombres_posibles, names(resultados))
    if (length(nombre_real) > 0) {
      valor <- resultados[[nombre_real[1]]]
      cat(sprintf("%-20s: %s\n", etiqueta, valor))
    }
  }

  mostrar_campo(c("Comparacion", "Comparation"), "Comparison")
  mostrar_campo(c("Diferencia", "Difference"), "Mean difference")
  mostrar_campo(c("df"), "Degrees of freedom")
  mostrar_campo(c("SE"), "Standard error")
  mostrar_campo(c("t_critical"), "t critical value")

  if (all(c("CI_lower", "CI_upper") %in% names(resultados))) {
    cat(sprintf("%-20s: [%s, %s]\n", "Confidence interval",
                resultados$CI_lower, resultados$CI_upper))
  }

  mostrar_campo("p_value", "p-value")

  # Advertencia si p esta cerca de 0.05
  if ("p_value" %in% names(resultados)) {
    p <- resultados$p_value
    if (!is.na(p) && abs(p - 0.05) < 0.005) {
      cat("\n The p-value is near the significance threshold (0.05). Interpret with caution.\n")
    }
  }

  mostrar_campo("p_hat", "P(X < Y) + 0.5 P(X = Y)")
  mostrar_campo(c("Significancia", "Sig"), "Significance")
  mostrar_campo("Valor_Critico", "Critical value")

  cat("\nGroup means (ordered from highest to lowest):\n")
  ordenado <- round(promedios[orden], 3)
  for (g in names(ordenado)) {
    cat("  ", g, ": ", ordenado[g], "\n", sep = "")
  }

  resultados$Metodo <- metodo  # metodo a la salida
  invisible(resultados)
}

