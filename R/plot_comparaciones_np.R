#' Plot for rank-based (non-parametric) multiple-comparison tests
#'
#' Rank-based procedures (Dunn, Conover-Iman, Nemenyi, DSCF, Bonferroni-Mann-Whitney)
#' compare distributions through ranks, not means, so a bar chart of means can
#' contradict the test. This method shows the data instead: a boxplot with the
#' individual observations (default), a violin plot, or the mean ranks used by
#' the test. Groups are ordered from highest to lowest mean rank and the sample
#' size of each group is shown on the x axis.
#'
#' @param x An object of class \code{comparaciones_np} (returned by
#'   \code{\link{DunnTest}}, \code{\link{ConoverTest}}, \code{\link{NemenyiTest}},
#'   \code{\link{DSCFTest}} or \code{\link{BonferroniNPTest}}).
#' @param tipo Type of plot: \code{"boxplot"} (boxplot + points, default),
#'   \code{"violin"} (violin + narrow boxplot) or \code{"rangos"} (mean rank per group).
#' @param letras If \code{TRUE}, adds the compact letter display above each group.
#'   Groups sharing a letter were not found to be different. Default \code{FALSE}.
#' @param alpha Significance threshold for the letters (default 0.05).
#' @param p_column Which p-value column to use for the letters: "auto", "p_ajustada", "p_value" or "p".
#' @param fill Fill color.
#' @param point_color Color of the individual points.
#' @param label_size Letter size.
#' @param angle_x Angle of the x-axis labels.
#' @param ... Not used.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(d_e, package = "Analitica")
#' res <- DunnTest(Sueldo_actual ~ labor, data = d_e)
#' plot(res)
#' plot(res, tipo = "violin")
#' plot(res, tipo = "rangos", letras = TRUE)
#'
#' @seealso \code{\link{plot.comparaciones}}
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_violin geom_point geom_segment geom_text scale_x_discrete labs theme_minimal theme element_text
#' @importFrom rlang .data
plot.comparaciones_np <- function(
    x,
    tipo = c("boxplot", "violin", "rangos"),
    letras = FALSE,
    alpha = 0.05,
    p_column = c("auto", "p_ajustada", "p_value", "p"),
    fill = "steelblue",
    point_color = "steelblue4",
    label_size = 5,
    angle_x = 0,
    ...
) {
  tipo <- match.arg(tipo)
  p_column <- match.arg(p_column)
  metodo <- if (!is.null(x$Metodo)) x$Metodo else "Comparaciones"
  orden <- x$Orden_Medias

  if (is.null(x$Datos) && tipo != "rangos") {
    warning("The object does not contain the data ('Datos'); showing mean ranks instead.")
    tipo <- "rangos"
  }

  if (tipo == "rangos") {
    n_g <- if (!is.null(x$Datos)) table(x$Datos$grupo)[orden] else NULL
  } else {
    d <- x$Datos
    d$grupo <- factor(as.character(d$grupo), levels = orden)
    n_g <- table(d$grupo)
  }
  etiquetas_x <- if (is.null(n_g)) orden else paste0(orden, "\nn = ", as.integer(n_g))

  if (tipo == "rangos") {
    r <- data.frame(grupo = factor(orden, levels = orden),
                    rango = as.numeric(x$Promedios[orden]))
    p <- ggplot2::ggplot(r, ggplot2::aes(x = .data$grupo, y = .data$rango)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$grupo, y = 0, yend = .data$rango),
                            colour = fill, linewidth = 1) +
      ggplot2::geom_point(size = 3.5, colour = point_color)
    y_top <- r$rango
    titulo <- paste0("Mean ranks (", metodo, ")")
    y_lab <- "Mean rank"
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$grupo, y = .data$y))
    if (tipo == "violin") {
      p <- p +
        ggplot2::geom_violin(fill = fill, alpha = 0.3, colour = NA) +
        ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white")
    } else {
      p <- p +
        ggplot2::geom_boxplot(fill = fill, alpha = 0.35, outlier.shape = NA, width = 0.6) +
        ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 1.2,
                             colour = point_color)
    }
    y_top <- as.numeric(tapply(d$y, d$grupo, max, na.rm = TRUE))
    titulo <- paste0("Distribution by group (", metodo, ")")
    y_lab <- "Response"
  }

  if (isTRUE(letras)) {
    let <- .letras_cld(x, alpha = alpha, p_column = p_column)
    rango_y <- diff(range(y_top, 0, na.rm = TRUE))
    if (!is.finite(rango_y) || rango_y == 0) rango_y <- 1
    df_let <- data.frame(grupo = factor(orden, levels = orden),
                         y = y_top + 0.06 * rango_y,
                         Letra = unname(let[orden]))
    p <- p + ggplot2::geom_text(data = df_let,
                                ggplot2::aes(x = .data$grupo, y = .data$y, label = .data$Letra),
                                size = label_size, fontface = "bold", inherit.aes = FALSE)
  }

  p +
    ggplot2::scale_x_discrete(labels = etiquetas_x) +
    ggplot2::labs(
      title = titulo, x = "Group", y = y_lab,
      caption = if (isTRUE(letras)) {
        "Groups sharing a letter were not found to be significantly different"
      } else NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle_x,
                                                       hjust = if (angle_x == 0) 0.5 else 1))
}
