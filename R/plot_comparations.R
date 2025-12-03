#' Generic plot for multiple-comparison tests (with multcompView letters) v2.0.1
#'
#' @param x An object of class \code{comparaciones}.
#' @param alpha Significance threshold for the letters (default 0.05).
#' @param p_column Which p-value column: "auto","p_ajustada","p_value","p".
#' @param horizontal If TRUE, draw horizontal bars.
#' @param fill Bar fill color.
#' @param label_size Letter size.
#' @param label_color Letter color.
#' @param angle_x Angle of x-axis labels (if \code{horizontal = FALSE}).
#' @param show_se If TRUE and \code{x$MSerror}, \code{x$N} exist, draws error bars.
#' @param se_type "se" (default) or "ci95".
#' @param pad_frac Fraction of y-span used to place letters (default 0.35).
#' @param errorbar_width Width of errorbar whiskers.
#' @param ... Not used.
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_errorbar labs theme_minimal coord_cartesian coord_flip theme element_text
#' @importFrom multcompView multcompLetters
#' @importFrom rlang .data
#' @importFrom stats setNames
plot.comparaciones <- function(
    x, alpha = 0.05,
    p_column = c("auto","p_ajustada","p_value","p"),
    horizontal = FALSE,
    fill = "steelblue",
    label_size = 5,
    label_color = "black",
    angle_x = 45,
    show_se = FALSE,
    se_type = c("se","ci95"),
    pad_frac = 0.35,
    errorbar_width = 0.2,
    ...
) {
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
    col_p <- pick_col(resultados, c("p_ajustada","padj","p.adjust","p.adjusted","p_value","p","P.value"))
  } else {
    col_p <- p_column
    if (!col_p %in% names(resultados)) stop(sprintf("Column '%s' does not exist in 'Resultados'.", col_p))
  }
  if (is.na(col_p)) stop("No p-value column found ('p_ajustada', 'p_value', 'p', 'P.value').")

  pares <- resultados[[col_comp]]
  pvals <- resultados[[col_p]]
  if (anyNA(pares) || anyNA(pvals)) stop("There are NA values in comparisons or p-values.")
  names(pvals) <- gsub(" - ", "-", pares, fixed = TRUE)

  # Letras base de multcompView
  letras <- multcompView::multcompLetters(pvals, threshold = alpha)$Letters

  # Tabla de medias
  df_medias <- data.frame(
    Grupo = names(promedios),
    Media = as.numeric(promedios),
    Letra = letras[names(promedios)],
    stringsAsFactors = FALSE
  )

  # Letras faltantes -> placeholder
  if (anyNA(df_medias$Letra)) {
    faltan <- which(is.na(df_medias$Letra))
    df_medias$Letra[faltan] <- LETTERS[seq_along(faltan)]
    warning("Some groups did not receive letters from 'multcompLetters'; assigned stand-in letters.")
  }

  # Orden de barras
  if (!is.null(orden) && all(orden %in% df_medias$Grupo)) {
    df_medias <- df_medias[match(orden, df_medias$Grupo), ]
  } else {
    df_medias <- df_medias[order(df_medias$Media, decreasing = TRUE), ]
  }
  df_medias$Grupo <- factor(df_medias$Grupo, levels = df_medias$Grupo)

  # === Reordenar letras de acuerdo al orden de las medias ===
  # Idea: recorrer las etiquetas en el orden actual de df_medias,
  # registrar la PRIMERA aparición de cada símbolo de letra,
  # y remapear a a,b,c,... en ese mismo orden (respetando combinaciones como "ab").
  make_targets <- function(n) {
    if (n <= length(letters)) letters[seq_len(n)] else c(letters, paste0("a", seq_len(n - length(letters))))
  }
  # Símbolos únicos en orden de aparición por df_medias$Letra (uniendo por grupo en orden actual)
  symbols_in_order <- character(0)
  for (lab in df_medias$Letra) {
    chars <- unique(strsplit(lab, "", fixed = TRUE)[[1]])
    for (ch in chars) if (!ch %in% symbols_in_order) symbols_in_order <- c(symbols_in_order, ch)
  }
  if (length(symbols_in_order)) {
    target <- make_targets(length(symbols_in_order))
    remap  <- setNames(target, symbols_in_order)
    remap_label <- function(s) {
      chars <- unique(strsplit(s, "", fixed = TRUE)[[1]])
      paste0(remap[chars], collapse = "")
    }
    df_medias$Letra <- vapply(df_medias$Letra, remap_label, character(1))
  }
  # === fin reordenamiento de letras ===

  # Colocación de letras
  rng  <- range(df_medias$Media, na.rm = TRUE)
  span <- diff(rng); if (span == 0) span <- max(1, abs(rng[2]))
  pad  <- max(1e-8, pad_frac * span)
  df_medias$y_lab <- df_medias$Media + ifelse(df_medias$Media >= 0, pad, -pad)

  # Barras de error
  have_se <- isTRUE(show_se) && !is.null(x$MSerror) && !is.null(x$N)
  se_type <- match.arg(se_type)
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
    } else if (is.null(x$MSerror)) {
      warning("MSerror not found; error bars will not be drawn.")
      have_se <- FALSE
    } else {
      se_vec <- sqrt(x$MSerror / x$N[group_ids])
      if (se_type == "ci95") {
        df_err <- if (!is.null(x$gl_error)) as.numeric(x$gl_error) else Inf
        mult <- if (is.finite(df_err)) stats::qt(0.975, df = df_err) else 1.96
        df_medias$SE_low <- df_medias$Media - mult * se_vec
        df_medias$SE_up  <- df_medias$Media + mult * se_vec
      } else {
        df_medias$SE_low <- df_medias$Media - se_vec
        df_medias$SE_up  <- df_medias$Media + se_vec
      }
    }
  }

  p <- ggplot2::ggplot(df_medias, ggplot2::aes(x = .data$Grupo, y = .data$Media)) +
    ggplot2::geom_col(fill = fill, width = 0.6)

  if (have_se) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$SE_low, ymax = .data$SE_up),
      width = errorbar_width
    )
  }

  p <- p +
    ggplot2::geom_text(ggplot2::aes(y = .data$y_lab, label = .data$Letra),
                       size = label_size, color = label_color) +
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
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, ymax + 0.1 * pad))
  }

  p
}
