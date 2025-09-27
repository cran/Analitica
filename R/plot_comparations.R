#' Generic plot for multiple-comparison tests (with multcompView letters)
#'
#' @param x An object of class \code{comparaciones} (includes "bonferroni" or "tukey").
#' @param alpha Significance threshold for the letters (default 0.05).
#' @param p_column Which p-value column to use: "auto" (detect), "p_ajustada", "p_value", or "p".
#' @param horizontal If TRUE, draw horizontal bars.
#' @param fill Bar fill color.
#' @param label_size Letter size.
#' @param angle_x Angle of x-axis labels (if \code{horizontal = FALSE}).
#' @param show_se If TRUE and \code{x$MSerror} and \code{x$N} exist, draws standard error (SE) bars.
#' @param ... Not used.
#' @return A \code{ggplot} object.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_errorbar labs theme_minimal coord_cartesian coord_flip theme element_text
#' @importFrom multcompView multcompLetters
#' @importFrom rlang .data
#'
#'
plot.comparaciones <- function(x, alpha = 0.05,
                               p_column = c("auto","p_ajustada","p_value","p"),
                               horizontal = FALSE,
                               fill = "steelblue",
                               label_size = 5,
                               angle_x = 45,
                               show_se = FALSE,
                               ...) {
  if (!inherits(x, "comparaciones")) stop("Object must be of class 'comparaciones'.")

  resultados <- x$Resultados
  promedios  <- x$Promedios
  orden      <- x$Orden_Medias
  metodo     <- if (!is.null(x$Metodo)) x$Metodo else "Comparaciones"
  termino    <- if (!is.null(x$Termino)) paste0(" - ", x$Termino) else ""

  pick_col <- function(df, choices) { nm <- intersect(choices, names(df)); if (length(nm)) nm[1] else NA_character_ }
  col_comp <- pick_col(resultados, c("Comparacion","Comparation","Comparison"))
  if (is.na(col_comp)) stop("No comparison column found (e.g., 'Comparacion').")

  p_column <- match.arg(p_column)
  if (p_column == "auto") {
    col_p <- pick_col(resultados, c("p_ajustada","padj","p.adjust","p.adjusted","p_value","p"))
  } else {
    col_p <- p_column
    if (!col_p %in% names(resultados)) stop(sprintf("Column '%s' does not exist in 'Resultados'.", col_p))
  }
  if (is.na(col_p)) stop("No p-value column found ('p_ajustada', 'p_value' or 'p').")

  pares <- resultados[[col_comp]]
  pvals <- resultados[[col_p]]
  if (anyNA(pares) || anyNA(pvals)) stop("There are NA values in comparisons or p-values.")
  names(pvals) <- gsub(" - ", "-", pares, fixed = TRUE)

  letras <- multcompView::multcompLetters(pvals, threshold = alpha)$Letters

  df_medias <- data.frame(
    Grupo = names(promedios),
    Media = as.numeric(promedios),
    Letra = letras[names(promedios)],
    stringsAsFactors = FALSE
  )
  if (!is.null(orden) && all(orden %in% df_medias$Grupo)) {
    df_medias <- df_medias[match(orden, df_medias$Grupo), ]
  } else {
    df_medias <- df_medias[order(df_medias$Media, decreasing = TRUE), ]
  }
  df_medias$Grupo <- factor(df_medias$Grupo, levels = df_medias$Grupo)

  rng  <- range(df_medias$Media, na.rm = TRUE); span <- diff(rng); if (span == 0) span <- max(1, abs(rng[2]))
  pad  <- 0.04 * span
  df_medias$y_lab <- df_medias$Media + ifelse(df_medias$Media >= 0, pad, -pad)

  # --- Robust SE bars ---
  have_se <- isTRUE(show_se) && !is.null(x$MSerror) && !is.null(x$N)
  if (have_se) {
    if (is.null(names(x$N)) && length(x$N) == length(promedios)) {
      names(x$N) <- names(promedios)
    }
    group_ids <- as.character(df_medias$Grupo)
    missing_n <- setdiff(group_ids, names(x$N))
    if (length(missing_n)) {
      warning("Missing sample sizes (N) for: ", paste(missing_n, collapse = ", "),
              ". Error bars will not be drawn.")
      have_se <- FALSE
    } else {
      se <- sqrt(x$MSerror / x$N[group_ids])
      df_medias$SE <- as.numeric(se)
    }
  }

  p <- ggplot2::ggplot(df_medias, ggplot2::aes(x = .data$Grupo, y = .data$Media)) +
    ggplot2::geom_col(fill = fill, width = 0.6)

  if (have_se) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$Media - .data$SE, ymax = .data$Media + .data$SE),
      width = 0.15
    )
  }

  p <- p +
    ggplot2::geom_text(ggplot2::aes(y = .data$y_lab, label = .data$Letra), size = label_size) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0("Group means (", metodo, termino, ")"),
      x = "Group", y = "Mean"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = if (horizontal) 0 else angle_x,
                                          hjust = if (horizontal) 0.5 else 1)
    )

  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  } else {
    ymax <- max(df_medias$y_lab, df_medias$Media, na.rm = TRUE)
    ymin <- min(0, min(df_medias$Media, df_medias$y_lab, na.rm = TRUE))
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, ymax + 0.02 * span))
  }

  p
}
