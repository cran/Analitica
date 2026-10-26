#' Print Method for Objects of Class 'comparaciones'
#'
#' Short description of a multiple-comparison object. Use \code{summary()} to
#' see the full results.
#'
#' @param x An object of class \code{"comparaciones"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.comparaciones <- function(x, ...) {
  cat("Multiple Comparison Object\n")
  cat("Method used:", x$Metodo, "\n")
  cat("Number of groups:", length(x$Promedios), "\n")
  cat("Number of pairwise comparisons:", nrow(x$Resultados), "\n")
  cat("Use summary() to view full results.\n")

  invisible(x)
}
