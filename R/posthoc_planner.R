#' Post Hoc Planner for FWER and Test Recommendation v1.5
#'
#' Computes the number of pairwise comparisons (\eqn{m}), expected FWER under independence
#' (no adjustment, Bonferroni, Sidak) for the term(s) you plan to compare, and proposes a
#' reasonable post hoc test (e.g., Tukey HSD, Holm, Gabriel, Duncan, LSD, Scheffe) depending
#' on design size, (un)balanced group sizes, and variance assumptions.
#'
#' Advadata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCCntages:
#' - One-shot planning: know \eqn{m}, FWER (informative) and a suggested post hoc before testing.
#' - Handles main effects or full-cell comparisons (interactions).
#' - Output formatted by rows (metrics as rows; each unit as a column).
#'
#' Disadvantages:
#' - Heuristic recommendation (not a formal decision rule).
#' - FWER values "under independence" are informative, not guaranteed when tests are dependent.
#'
#' @param modelo An \code{aov} or \code{lm} object (full model: includes blocks, factors, etc.).
#' @param comparar Character vector with the name(s) of the factor(s) to compare:
#'   - One name: main effect (e.g., "treatment" or "A")
#'   - Several names: if \code{scope="cells"} compares \code{A:B:...} cells;
#'     if \code{scope="factor"}, reports each factor separately.
#'   If omitted, uses all factors in \code{modelo$xlevels} when \code{scope="factor"},
#'   or the first factor when \code{scope="cells"} unless multiple are provided.
#' @param alpha Global significance level (FWER target), default \code{0.05}.
#' @param scope \code{"factor"} compares each factor separately (main effects);
#'   \code{"cells"} compares all interaction cells (full combination of \code{comparar}).
#' @param equal_var Logical; assume homoscedasticity (default \code{TRUE}).
#' @param unequal_n Logical; expect moderate group size imbalance (default \code{FALSE}).
#' @param independence Logical; if \code{TRUE} report FWER "under independence" for Bonf/Sidak (default \code{TRUE}).
#' @param liberal_ok Logical; allow more liberal suggestions (e.g., LSD/Duncan) when power is prioritized (default \code{FALSE}).
#' @param orientation Output table orientation: \code{"rows"} (metrics as rows, default) or \code{"cols"}.
#' @param digits Numeric formatting (non-count numeric columns), default \code{4}.
#' @param percent_digits Digits for percentage rendering of FWERs, default \code{1}.
#'
#' @return A \code{data.frame}. By default (rows orientation), the first column is \code{Metrica}
#'   and subsequent columns are "Unidades" (each factor or the full-cell unit).
#'   Metrics include: number of levels \code{g}, comparisons \code{m}, global \eqn{\alpha},
#'   Bonferroni/Sidak per-comparison \eqn{\alpha}, and FWERs (under independence).
#'   Also includes a \code{Sugerencia post hoc} string with the recommended test.
#'
#' @references
#' Hochberg, Y. & Tamhane, A. C. (1987). Multiple Comparison Procedures.
#'
#' Wiley. Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46, 561-584.
#'
#' @examples
#' # DCA (one-way):
#' # data(d_e, package = "Analitica")
#' # mod1 <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
#' # Posthoc_planner(mod1, comparar = "as.factor(labor)")
#'
#' # RCBD / DBA: y ~ tratamiento + bloque
#' # mod2 <- aov(y ~ as.factor(labor) + Sexo, data = d_e)
#' # Posthoc_planner(mod2, comparar = "as.factor(labor)", scope = "factor")
#'
#' # Factorial: y ~ A * B  (compare full cells A:B)
#' # mod3 <- aov(y ~ as.factor(labor) * Sexo, data = d_e)
#' # Posthoc_planner(mod3, comparar = c("as.factor(labor)","Sexo"), scope = "cells")
#'
#' @export
Posthoc_planner <- function(modelo,
                            comparar = NULL,
                            alpha = 0.05,
                            scope = c("factor","cells"),
                            equal_var = TRUE,
                            unequal_n = FALSE,
                            independence = TRUE,
                            liberal_ok = FALSE,
                            orientation = c("rows","cols"),
                            digits = 4,
                            percent_digits = 1) {

  scope <- match.arg(scope)
  orientation <- match.arg(orientation)

  if (is.null(modelo$model))
    stop("The 'modelo' object must contain the data (fit aov/lm with embedded data).")

  # -------- Helpers --------
  choose2     <- function(g) ifelse(g >= 2, g*(g-1)/2, 0)
  fwer_indep  <- function(alpha_star, m) ifelse(m > 0, 1 - (1 - alpha_star)^m, 0)
  alpha_bonf  <- function(alpha, m) ifelse(m > 0, alpha / m, NA_real_)
  alpha_sidak <- function(alpha, m) ifelse(m > 0, 1 - (1 - alpha)^(1/m), NA_real_)

  recommend_test <- function(g, m, equal_var, unequal_n, scope, liberal_ok){
    if (!equal_var) {
      return("Games-Howell (heterocedastico) o Holm (robusto)")
    }
    if (unequal_n) {
      if (g <= 12) return("Gabriel (desbalance moderado) o Tukey HSD (Kramer) / Holm")
      return("Holm (m grande) o Tukey HSD (Kramer)")
    }
    # balanceado y homocedastico
    if (scope == "cells" && m > 45) {
      return("Holm (potente) o Scheffe (contrastes generales)")
    }
    if (m <= 6) {
      return(if (liberal_ok) "LSD (gate F) o Tukey HSD / Holm"
             else            "Tukey HSD o Holm (LSD si aceptas mayor Tipo I)")
    }
    if (m <= 20) {
      return(if (liberal_ok) "Duncan (mas potente) o Tukey HSD / Holm"
             else            "Holm (step-down) o Tukey HSD")
    }
    if (m <= 45) return("Holm (step-down) o Tukey HSD")
    "Holm (step-down) o Scheffe (si hay mas que pares)"
  }

  # -------- Parse factors / groups --------
  xlv <- modelo$xlevels
  if (is.null(xlv) || length(xlv) == 0) {
    stop("Model has no factors in 'xlevels'. Convert categorical predictors to factor.")
  }

  mf <- modelo$model

  # If comparar is NULL: for scope="factor" , then  all factors; for "cells" then first factor unless multiple supplied
  if (is.null(comparar)) {
    if (scope == "factor") {
      comparar <- names(xlv)                      # all factors
    } else { # "cells"
      comparar <- names(xlv)[1]                   # first factor only (unless user provides several)
    }
  }
  comparar <- as.character(comparar)

  # Ensure the listed terms exist and are factors
  for (nm in comparar) {
    if (!nm %in% names(mf))
      stop(sprintf("Term '%s' is not found in model data.", nm))
    if (!is.factor(mf[[nm]])) mf[[nm]] <- factor(mf[[nm]])
  }

  # Count levels per requested units
  level_count <- function(varname) nlevels(mf[[varname]])

  # -------- Build per-unit summaries --------
  build_row <- function(label, g, alpha, equal_var, unequal_n, scope, independence, liberal_ok){
    m  <- choose2(g)
    aB <- alpha_bonf(alpha, m)
    aS <- alpha_sidak(alpha, m)
    data.frame(
      Unidad = label,
      `N niveles (g)` = g,
      `Comparaciones (m)` = m,
      `alpha global (FWER)` = alpha,
      `alpha_B = alpha/m (Bonf)` = aB,
      `FWER Bonf (indep.)` = if (independence) fwer_indep(aB, m) else NA_real_,
      `alpha* de Sidak` = aS,
      `FWER Sidak (indep.)` = if (independence) fwer_indep(aS, m) else NA_real_,
      `FWER sin ajuste (indep.)` = if (independence) fwer_indep(alpha, m) else NA_real_,
      `Sugerencia post hoc` = recommend_test(g, m, equal_var, unequal_n, scope, liberal_ok),
      check.names = FALSE
    )
  }

  res_list <- list()

  if (scope == "factor") {
    # One column per factor requested
    for (nm in comparar) {
      g <- level_count(nm)
      res_list[[length(res_list)+1]] <-
        build_row(paste0("Factor: ", nm), g, alpha, equal_var, unequal_n, scope, independence, liberal_ok)
    }
  } else { # "cells": full combination of provided factors
    g_vec <- vapply(comparar, level_count, FUN.VALUE = numeric(1))
    G <- prod(g_vec)
    lbl <- paste0("Celdas (", paste(comparar, collapse=":"), ")")
    res_list[[1]] <- build_row(lbl, G, alpha, equal_var, unequal_n, scope, independence, liberal_ok)
  }

  res_wide <- do.call(rbind, res_list)
  rownames(res_wide) <- NULL

  # Nice column order
  col_order <- c("Unidad","N niveles (g)","Comparaciones (m)","alpha global (FWER)",
                 "alpha_B = alpha/m (Bonf)","FWER Bonf (indep.)",
                 "alpha* de Sidak","FWER Sidak (indep.)",
                 "FWER sin ajuste (indep.)","Sugerencia post hoc")
  col_order <- col_order[col_order %in% names(res_wide)]
  res_wide <- res_wide[, col_order, drop = FALSE]

  if (orientation == "cols") return(res_wide)

  # -------- Row-oriented presentation (metrics as rows) --------
  is_num <- sapply(res_wide, is.numeric)
  fmt_num <- function(x) ifelse(is.na(x), NA, formatC(x, format = "f", digits = digits))
  fmt_pct <- function(x) ifelse(is.na(x), NA, paste0(formatC(100*x, format = "f", digits = percent_digits), "%"))

  pct_cols <- intersect(c("FWER Bonf (indep.)","FWER Sidak (indep.)","FWER sin ajuste (indep.)"),
                        names(res_wide))

  res_fmt <- res_wide
  for (cl in names(res_fmt)) {
    if (cl %in% pct_cols) {
      res_fmt[[cl]] <- fmt_pct(res_fmt[[cl]])
    } else if (is_num[cl] && cl != "N niveles (g)" && cl != "Comparaciones (m)") {
      res_fmt[[cl]] <- fmt_num(res_fmt[[cl]])
    }
  }

  unidades <- res_fmt$Unidad
  metricas <- setdiff(names(res_fmt), "Unidad")

  # Use 'optional=TRUE' to avoid name mangling; no duplicate check.names here
  tabla_core <- as.data.frame(t(as.matrix(res_fmt[, metricas, drop = FALSE])),
                              optional = TRUE, stringsAsFactors = FALSE)
  colnames(tabla_core) <- unidades

  tabla_rows <- data.frame(
    Metrica = rownames(tabla_core),
    tabla_core,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  metrica_order <- c("N niveles (g)", "Comparaciones (m)", "alpha global (FWER)",
                     "alpha_B = alpha/m (Bonf)", "FWER Bonf (indep.)",
                     "alpha* de Sidak", "FWER Sidak (indep.)",
                     "FWER sin ajuste (indep.)", "Sugerencia post hoc")
  metrica_order <- metrica_order[metrica_order %in% tabla_rows$Metrica]
  tabla_rows <- tabla_rows[match(metrica_order, tabla_rows$Metrica), , drop = FALSE]

  tabla_rows
}

