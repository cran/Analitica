#' Summary Method for Objects of Class 'comparaciones'
#'
#' Prints the method, the group means (or mean ranks for rank-based procedures),
#' their order from highest to lowest, and the table of pairwise comparisons.
#'
#' @param object An object of class \code{"comparaciones"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{object}.
#' @export
summary.comparaciones <- function(object, ...) {
  es_rango <- inherits(object, "comparaciones_np")
  etiqueta <- if (es_rango) "Mean ranks" else "Group means"

  cat("=====================================\n")
  cat("  Multiple Comparison Method Summary\n")
  cat("=====================================\n")
  cat("Method used:", object$Metodo, "\n")
  if (!is.null(object$Ajuste)) cat("P-value adjustment:", object$Ajuste, "\n")
  cat("\n")

  if (!is.null(object$KruskalWallis)) {
    kw <- object$KruskalWallis
    cat(sprintf(">> Kruskal-Wallis: H = %.4f, df = %d, p = %.4g\n\n",
                unname(kw$statistic), as.integer(kw$parameter), kw$p.value))
  }

  cat(">> ", etiqueta, " (from highest to lowest):\n", sep = "")
  prom <- object$Promedios
  if (!is.null(object$Orden_Medias) && all(object$Orden_Medias %in% names(prom))) {
    prom <- prom[object$Orden_Medias]
  }
  print(round(prom, 4))

  cat("\n>> Pairwise comparisons:\n")
  print(object$Resultados, row.names = FALSE)

  invisible(object)
}
