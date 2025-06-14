#' @export
summary.comparaciones <- function(object, ...) {
  cat("=====================================\n")
  cat("  Multiple Comparison Method Summary\n")
  cat("=====================================\n")
  cat("Method used:", object$Metodo, "\n\n")

  cat(">> Group means:\n")
  print(round(object$Promedios, 4))

  cat("\n>> Order of means (from highest to lowest):\n")
  print(object$Orden_Medias)

  cat("\n>> Pairwise comparisons:\n")
  print(object$Resultados)

  invisible(object)
}

