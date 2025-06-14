#' Generic plot for multiple comparison tests (with multcompView letters)
#'
#' @param x An object of class \code{comparaciones}.
#' @param ... Additional arguments (currently not used).
#'
#' @return No return value. Called for side effects: displays a bar plot with significance letters.
#'
#' @description
#' This function generates a bar plot displaying group means along with significance letters
#' based on multiple comparisons. It uses \code{multcompView} to assign letters indicating
#' statistically different groups.
#'
#' @seealso \code{\link[multcompView]{multcompLetters}}, \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#' # Assuming you have an object of class 'comparaciones' named res
#' # plot(res)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal coord_cartesian theme element_text
#' @importFrom multcompView multcompLetters
plot.comparaciones <- function(x, ...) {
  if (!inherits(x, "comparaciones")) {
    stop("El objeto debe ser de clase 'comparaciones'")
  }

  resultados <- x$Resultados
  promedios <- x$Promedios
  orden <- x$Orden_Medias
  metodo <- if (!is.null(x$Metodo)) x$Metodo else "Comparaciones"

  # Preparar data.frame de comparaciones para multcompView
  # NOTA: Necesitamos una matriz de p-valores para multcompLetters
  grupos <- unique(unlist(strsplit(resultados$Comparacion, " - ")))
  matriz_p <- matrix(1, nrow = length(grupos), ncol = length(grupos),
                     dimnames = list(grupos, grupos))

  for (i in 1:nrow(resultados)) {
    par <- unlist(strsplit(resultados$Comparacion[i], " - "))
    pval <- resultados$p_value[i]
    matriz_p[par[1], par[2]] <- pval
    matriz_p[par[2], par[1]] <- pval
  }

  # Asignar letras usando multcompLetters
  letras <- multcompView::multcompLetters(matriz_p, threshold = 0.05)$Letters

  # Preparar data.frame final
  df_medias <- data.frame(
    Grupo = names(promedios),
    Media = as.numeric(promedios),
    Letra = letras[names(promedios)],
    stringsAsFactors = FALSE
  )
  df_medias <- df_medias[match(orden, df_medias$Grupo), ]
  df_medias$Grupo <- factor(df_medias$Grupo, levels = df_medias$Grupo)

  # Generar grafico
  ggplot(df_medias, aes(x = Grupo, y = Media)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    geom_text(aes(label = Letra), vjust = -0.5, size = 5) +
    theme_minimal() +
    labs(
      title = paste("Promedios por grupo ", metodo),
      x = "Grupo", y = "Media"
    ) +
    coord_cartesian(ylim = c(0, max(df_medias$Media, na.rm = TRUE) * 1.1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
