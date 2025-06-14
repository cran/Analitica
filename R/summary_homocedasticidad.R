#' Summary Method for Homoscedasticity Test Results
#'
#' Displays a summary of variance homogeneity tests such as Bartlett, Fligner-Killeen,
#' or Levene, applied to a fitted formula using numeric data and groupings.
#'
#' @param object An object of class \code{"homocedasticidad"}.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns the input object (invisible). Printed output includes:
#' test name, statistic, degrees of freedom, p-value, and decision at the 0.05 level.
#' @export
summary.homocedasticidad <- function(object, ...) {
  cat("\n--- Homoscedasticity Test Summary ---\n\n")

  method <- if (!is.null(object$Method)) object$Method else "Unknown"
  stat   <- if (!is.null(object$Statistic)) object$Statistic else NA
  df     <- object$df
  pval   <- if (!is.null(object$p_value)) object$p_value else NA
  sig    <- if (!is.null(object$Significance)) object$Significance else ""
  decision <- if (!is.null(object$Decision)) object$Decision else "[Not available]"

  cat("Method applied         :", method, "\n")

  if (grepl("Levene", method, ignore.case = TRUE)) {
    cat("F Statistic            :", stat, "\n")
    if (is.numeric(df) && length(df) == 2 &&
        all(c("df_between", "df_within") %in% names(df))) {
      cat("Degrees of freedom     :", df["df_between"], "(between),",
          df["df_within"], "(within)\n")
    } else {
      cat("Degrees of freedom     : [Invalid or missing]\n")
    }
  } else if (grepl("Bartlett", method, ignore.case = TRUE) ||
             grepl("Fligner", method, ignore.case = TRUE)) {
    cat("Chi-squared Statistic  :", stat, "\n")
    cat("Degrees of freedom     :", if (!is.null(df)) df else "[Missing]", "\n")
  } else {
    cat("Test Statistic         :", stat, "\n")
    cat("Degrees of freedom     :", if (!is.null(df)) df else "[Missing]", "\n")
  }

  cat("p-value                :", pval, sig, "\n")  # <- AÑADIDA AQUÍ LA SIGNIFICANCIA
  cat("Decision (alpha = 0.05):", decision, "\n")
  cat("----------------------------------------\n")

  invisible(object)
}


