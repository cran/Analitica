# Funciones internas compartidas por los procedimientos de comparaciones multiples.
# No se exportan.

#' @noRd
.sig_stars <- function(p, alpha = 0.05) {
  # "ns" si p >= alpha; si es significativo, asteriscos segun 0.01 y 0.001
  out <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  out[p >= alpha] <- "ns"
  out[is.na(p)] <- NA_character_
  out
}

#' @noRd
.pares_ordenados <- function(valores) {
  # Todos los pares (g1, g2) con g1 de mayor valor que g2 (orden de mayor a menor)
  orden <- names(sort(valores, decreasing = TRUE))
  utils::combn(orden, 2, simplify = FALSE)
}

#' @noRd
.factor_empates <- function(y) {
  # C = 1 - sum(t^3 - t) / (N^3 - N), t = tamano de cada grupo de empates
  N <- length(y)
  t <- table(y)
  1 - sum(t^3 - t) / (N^3 - N)
}

#' @noRd
.agregar_info_pares <- function(out, g1, g2, p_value, p_ajustada = NULL) {
  # Guarda los pares y los p-valores sin redondear, para que plot() arme
  # las letras con los mismos valores que definen la columna Significancia.
  out$Pares <- data.frame(g1 = g1, g2 = g2, stringsAsFactors = FALSE)
  pe <- data.frame(p_value = as.numeric(p_value))
  if (!is.null(p_ajustada)) pe$p_ajustada <- as.numeric(p_ajustada)
  out$P_exactos <- pe
  out
}

#' @noRd
.letras_cld <- function(x, alpha = 0.05, p_column = "auto") {
  resultados <- x$Resultados
  pick_col <- function(df, choices) {
    nm <- intersect(choices, names(df)); if (length(nm)) nm[1] else NA_character_
  }

  if (p_column == "auto") {
    col_p <- pick_col(resultados, c("p_ajustada", "padj", "p.adjust", "p.adjusted",
                                    "p_value", "p", "P.value"))
  } else {
    col_p <- p_column
    if (!col_p %in% names(resultados)) {
      stop(sprintf("Column '%s' does not exist in 'Resultados'.", col_p))
    }
  }
  if (is.na(col_p)) stop("No p-value column found ('p_ajustada', 'p_value', 'p', 'P.value').")

  # Pares y p-valores: preferir los guardados sin redondear
  if (!is.null(x$Pares) && !is.null(x$P_exactos) && col_p %in% names(x$P_exactos)) {
    g1 <- x$Pares$g1; g2 <- x$Pares$g2
    pvals <- x$P_exactos[[col_p]]
  } else {
    col_comp <- pick_col(resultados, c("Comparacion", "Comparation", "Comparison"))
    if (is.na(col_comp)) stop("No comparison column found (e.g., 'Comparacion').")
    partes <- strsplit(as.character(resultados[[col_comp]]), " - ", fixed = TRUE)
    if (any(lengths(partes) != 2)) {
      stop("Could not split the comparison labels into two groups.")
    }
    g1 <- vapply(partes, `[`, character(1), 1)
    g2 <- vapply(partes, `[`, character(1), 2)
    pvals <- resultados[[col_p]]
  }
  if (anyNA(pvals)) stop("There are NA values in the p-values.")

  grupos <- names(x$Promedios)
  if (!is.null(x$Orden_Medias) && all(x$Orden_Medias %in% grupos)) grupos <- x$Orden_Medias

  # Matriz simetrica de p-valores (evita depender de nombres con guiones)
  M <- matrix(1, length(grupos), length(grupos), dimnames = list(grupos, grupos))
  for (i in seq_along(pvals)) {
    M[g1[i], g2[i]] <- pvals[i]
    M[g2[i], g1[i]] <- pvals[i]
  }
  letras <- multcompView::multcompLetters(M, threshold = alpha)$Letters
  letras <- letras[grupos]

  # Reasignar letras en orden de aparicion (a para el grupo mayor, etc.)
  simbolos <- character(0)
  for (lab in letras) {
    for (ch in unique(strsplit(lab, "", fixed = TRUE)[[1]])) {
      if (!ch %in% simbolos) simbolos <- c(simbolos, ch)
    }
  }
  if (length(simbolos)) {
    destino <- if (length(simbolos) <= length(letters)) {
      letters[seq_along(simbolos)]
    } else {
      c(letters, paste0("a", seq_len(length(simbolos) - length(letters))))
    }
    remap <- stats::setNames(destino, simbolos)
    letras <- vapply(letras, function(s) {
      ch <- unique(strsplit(s, "", fixed = TRUE)[[1]])
      paste0(sort(remap[ch]), collapse = "")
    }, character(1))
  }
  stats::setNames(letras, grupos)
}
