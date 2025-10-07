#' Summary Method for Homoscedasticity Test Results (extended) v2.0
#'
#' Displays a summary of variance homogeneity tests such as Bartlett, Fligner-Killeen,
#' or Levene (1-via global and if exists, factorial decomposition on |Y - centro_celda|).
#'
#' @param object An object of class \code{"homocedasticidad"}.
#' @param digits Number of digits for F; default 4.
#' @param ... Currently ignored.
#' @return Invisibly returns the input object (invisible).
#' @export
summary.homocedasticidad <- function(object, digits = 4, ...) {
  cat("\n--- Homoscedasticity Test Summary ---\n\n")

  method   <- if (!is.null(object$Method)) object$Method else "Unknown"
  stat     <- if (!is.null(object$Statistic)) object$Statistic else NA
  df       <- object$df
  pval     <- if (!is.null(object$p_value)) object$p_value else NA
  sig      <- if (!is.null(object$Significance)) object$Significance else ""
  decision <- if (!is.null(object$Decision)) object$Decision else "[Not available]"

  cat("Method applied         :", method, "\n")

  if (grepl("Levene", method, ignore.case = TRUE)) {
    # Levene global (1-via sobre celdas)
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

  cat("p-value                :", pval, sig, "\n")
  cat("Decision (alpha = 0.05):", decision, "\n")
  cat("----------------------------------------\n")

  ## ====== Descomposicion factorial (si existe) ======
  ## Espera que tu funcion Levene.Test (factorial) haya guardado:
  ## object$twoway$table (data.frame ANOVA), object$twoway$method, object$twoway$note (opc)
  ## object$twoway$cols (lista con nombres de columnas df/F/p) y object$twoway$interpretation
  if (!is.null(object$twoway) && is.list(object$twoway) && !is.null(object$twoway$table)) {
    cat("\n>> ANOVA sobre desviaciones absolutas (descomposicion factorial)\n")
    cat("Metodo:                 ", if (!is.null(object$twoway$method)) object$twoway$method else "[desconocido]", "\n", sep = "")
    if (!is.null(object$twoway$note)) cat("Nota:                    ", object$twoway$note, "\n", sep = "")
    cat("\n")

    tab <- object$twoway$table
    rn  <- rownames(tab)

    # Descubrir columnas (flexible a distintos paquetes)
    df_col <- if (!is.null(object$twoway$cols$df)) object$twoway$cols$df else "Df"
    f_col  <- if (!is.null(object$twoway$cols$F))  object$twoway$cols$F  else {
      cand <- c("F value", "F")
      cand[cand %in% colnames(tab)][1]
    }
    p_col  <- if (!is.null(object$twoway$cols$p))  object$twoway$cols$p  else {
      cand <- grep("^Pr\\(>F\\)", colnames(tab), value = TRUE)
      cand[1]
    }

    DF <- if (!is.null(df_col) && df_col %in% colnames(tab)) tab[[df_col]] else NA
    Fv <- if (!is.null(f_col)  && f_col  %in% colnames(tab)) tab[[f_col]]  else NA
    pv <- if (!is.null(p_col)  && p_col  %in% colnames(tab)) tab[[p_col]]  else NA

    sig_star <- function(p) {
      if (is.na(p)) return("")
      if (p < 0.001) return("***")
      if (p < 0.01)  return("**")
      if (p < 0.05)  return("*")
      "ns"
    }

    # Intento de interpretacion, si se proveyo
    interp_vec <- rep(NA_character_, length(rn))
    names(interp_vec) <- rn
    if (!is.null(object$twoway$interpretation)) {
      inames <- names(object$twoway$interpretation)
      interp_vec[inames] <- object$twoway$interpretation
    }

    # Mapa opcional de nombres (si quieres ver A, B, A:B)
    # Cambia este mapeo segun tus factores reales:
    # term_map <- c("Suelo" = "A", "Fertilizante" = "B", "Suelo:Fertilizante" = "A:B")
    # rn_clean <- ifelse(rn %in% names(term_map), term_map[rn], rn)
    rn_clean <- rn  # por defecto, muestra nombres reales de factores

    apa <- data.frame(
      Term           = rn_clean,
      DF             = DF,
      F              = ifelse(is.na(Fv), NA, round(Fv, digits)),
      p              = ifelse(is.na(pv), NA, signif(pv, 4)),
      Signif         = vapply(pv, sig_star, character(1)),
      Interpretacion = unname(interp_vec[rn]),
      row.names = NULL
    )

    print(apa, row.names = FALSE)
    cat("----------------------------------------\n")
  }

  invisible(object)
}


