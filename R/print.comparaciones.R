#' @export
print.comparaciones <- function(x, ...) {
  cat("Multiple Comparison Object\n")
  cat("Method used:", x$Metodo, "\n")
  cat("Number of groups:", length(x$Promedios), "\n")
  cat("Number of pairwise comparisons:", nrow(x$Resultados), "\n")
  cat("Use summary() to view full results.\n")

  invisible(x)
}
