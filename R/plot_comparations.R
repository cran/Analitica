#' Plot for parametric multiple-comparison tests (with multcompView letters)
#'
#' Bar chart of group means with compact letter display. Groups sharing a letter
#' are not significantly different. For rank-based (non-parametric) procedures
#' see \code{\link{plot.comparaciones_np}}.
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
#' @seealso \code{\link{plot.comparaciones_np}}
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

  promedios  <- x$Promedios
  orden      <- x$Orden_Medias
  metodo     <- if (!is.null(x$Metodo)) x$Metodo else "Comparaciones"
  termino    <- if (!is.null(x$Termino)) paste0(" - ", x$Termino) else ""

  # Letras de grupos homogeneos (calculadas con p-valores sin redondear si existen)
  letras <- .letras_cld(x, alpha = alpha, p_column = match.arg(p_column))

  df_medias <- data.frame(
    Grupo = names(promedios),
    Media = as.numeric(promedios),
    Letra = unname(letras[names(promedios)]),
    stringsAsFactors = FALSE
  )
  if (anyNA(df_medias$Letra)) {
    stop("Some groups did not receive letters; check that 'Resultados' contains all pairs.")
  }

  # Orden de barras
  if (!is.null(orden) && all(orden %in% df_medias$Grupo)) {
    df_medias <- df_medias[match(orden, df_medias$Grupo), ]
  } else {
    df_medias <- df_medias[order(df_medias$Media, decreasing = TRUE), ]
  }
  df_medias$Grupo <- factor(df_medias$Grupo, levels = df_medias$Grupo)

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
