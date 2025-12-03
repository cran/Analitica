#' Post Hoc Planner for FWER and Test Recommendation v1.6
#'
#' @description
#' One-shot planner for factor or cell comparisons, reporting m, FWER,
#' suggested adjustments (Bonferroni/Sidak) and a post hoc recommendation
#' (Holm, Tukey, Duncan, Gabriel, Scheffe, SNK, etc.) before testing.
#'
#' @references
#' #' Bonferroni, C. (1936). *Teoria statistica delle classi e calcolo delle probabilità*.
#' Pubblicazioni del R. Istituto Superiore di Scienze Economiche e Commerciali di Firenze.
#'
#' Fisher, R. A. (1935). *The design of experiments*. Oliver & Boyd.
#'
#' Duncan, D. B. (1955). Multiple range and multiple F tests. *Biometrics, 11*(1), 1–42.
#'
#' Gabriel, K. R. (1978). A simple method of multiple comparisons of means.
#' *Journal of the American Statistical Association, 73*(364), 724–729.
#'
#' Games, P. A., & Howell, J. F. (1976). Pairwise multiple comparison procedures with
#' unequal n’s and/or variances: A Monte Carlo study. *Journal of Educational Statistics, 1*(2), 113–125.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics, 6*(2), 65–70.
#'
#' Newman, D. (1939). The distribution of range in samples from a normal population,
#' expressed in terms of an independent estimate of standard deviation.
#' *Biometrika, 31*(1/2), 20–36.
#'
#' Scheffé, H. (1953). A method for judging all contrasts in the analysis of variance.
#' *Biometrika, 40*(1–2), 87–104.
#'
#' Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate
#' normal distributions. *Journal of the American Statistical Association, 62*(318), 626–633.
#'
#' Tukey, J. W. (1949). Comparing individual means in the analysis of variance.
#' *Biometrics, 5*(2), 99–114.
#'
#' @param model aov or lm object (complete model). Data are reconstructed with model.frame().
#' @param compare Character with the name(s) of the factor(s) to compare:
#'   - One name: main effect.
#'   - Several names: if scope="cells" compares A:B:... cells; if scope="factor", reports each factor.
#'   If omitted, uses all factors when scope="factor", or the first factor when scope="cells".
#' @param alpha Overall significance level (FWER target), default 0.05.
#' @param scope "factor" compares each factor separately; "cells" compares interaction cells.
#' @param equal_var Logical; assume homoscedasticity (default TRUE).
#' @param unequal_n Logical; expect moderate imbalance of group sizes (default FALSE).
#' @param independence Logical; if TRUE reports FWER "under independence" (default TRUE).
#' @param liberal_ok Logical; allows more liberal suggestions (LSD/Duncan/SNK) (default FALSE).
#' @param orientation "rows" (metrics as rows, default) or "cols".
#' @param digits Decimal places for numeric output, default 4.
#' @param percent_digits Decimal places for percentages, default 1.
#' @param observed_cells Logical; in scope="cells", count only observed cells (drop NA). Default TRUE.
#'
#'
#' @return data.frame.
#'   - orientation="rows": first column "Metric", rest columns are units (factor/cells).
#'   - orientation="cols": one row per unit, metrics as columns.
#'   Includes: g levels, m comparisons, global alpha, Bonferroni/Sidak alphas,
#'   FWERs (under independence), "Suggested p-value adjustment" and "Post hoc suggestion".
#'
#' @examples
#' ## =========================
#' ## Ejemplo 1: One-way ANOVA
#' ## =========================
#' # Datos: PlantGrowth (3 grupos balanceados)
#' data(PlantGrowth)
#' m1 <- aov(weight ~ group, data = PlantGrowth)

#' # Comparar por factor (default scope="factor")
#' Posthoc_planner(m1)
#' # Variante: salida por columnas
#' Posthoc_planner(m1, orientation = "cols")
#' # Variante: alpha más estricto
#' Posthoc_planner(m1, alpha = 0.01)
#'
#'## ==============================================
#' ## Ejemplo 2: Dos factores y comparación de celdas
#' ## ==============================================
#' # Datos: ToothGrowth (suplemento x dosis)
#' data(ToothGrowth)
#' TG <- ToothGrowth
#' TG$dose <- factor(TG$dose)  # tratar "dose" como factor

#' m2 <- aov(len ~ supp * dose, data = TG)

#' # scope="cells" compara celdas de la interacción (supp:dose)
#' Posthoc_planner(
#'   m2,
#' compare = c("supp","dose"),
#' scope   = "cells",        # comparar celdas
#' observed_cells = TRUE     # contar solo celdas observadas (default)
#' )

#' # También puedes pedir el resumen por factor dentro del mismo modelo
#' Posthoc_planner(
#' m2,
#' compare = c("supp","dose"),
#' scope   = "factor"        # reporte por cada factor por separado
#' )
#'
#'
#' @export
Posthoc_planner <- function(model,
                            compare = NULL,
                            alpha = 0.05,
                            scope = c("factor","cells"),
                            equal_var = TRUE,
                            unequal_n = FALSE,
                            independence = TRUE,
                            liberal_ok = FALSE,
                            orientation = c("rows","cols"),
                            digits = 4,
                            percent_digits = 1,
                            observed_cells = TRUE) {

  scope <- match.arg(scope)
  orientation <- match.arg(orientation)

  # -------- Reconstruir datos del model --------
  mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (is.null(mf)) stop("The model data could not be reconstructed with model.frame(model).")

  # -------- Validaciones --------
  if (!is.numeric(alpha) || length(alpha)!=1 || is.na(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be a number in (0,1).")
  if (digits < 0 || percent_digits < 0) stop("digits y percent_digits must be >= 0.")

  # Detectar factores presentes
  factor_names <- names(Filter(is.factor, mf))
  if (length(factor_names) == 0)
    stop("The model contains no factor predictors in the reconstructed data (model.frame).")

  # Si compare es NULL
  if (is.null(compare) || length(compare) == 0) {
    if (scope == "factor") {
      compare <- factor_names
    } else { # "cells"
      compare <- factor_names[1]
    }
  }
  compare <- as.character(compare)

  # Validar existencia y asegurar factor + droplevels
  not_found <- setdiff(compare, names(mf))
  if (length(not_found))
    stop(sprintf(
      "The following terms are not in the model data: %s. Available: %s.",
      paste(not_found, collapse = ", "),
      paste(names(mf), collapse = ", ")
    ))
  for (nm in compare) {
    if (!is.factor(mf[[nm]])) mf[[nm]] <- factor(mf[[nm]])
    mf[[nm]] <- droplevels(mf[[nm]])
  }

  # -------- Helpers --------
  choose2     <- function(g) ifelse(g >= 2, g*(g-1)/2, 0)
  fwer_indep  <- function(alpha_star, m) ifelse(m > 0, 1 - (1 - alpha_star)^m, 0)
  alpha_bonf  <- function(alpha, m) ifelse(m > 0, alpha / m, NA_real_)
  alpha_sidak <- function(alpha, m) ifelse(m > 0, 1 - (1 - alpha)^(1/m), NA_real_)

  # Elegir ajuste p-valor sugerido entre Bonferroni y Sidak (segun FWER)
  choose_adjust <- function(alpha, m, independence) {
    if (is.na(m) || m <= 1) return("-")
    aB <- alpha_bonf(alpha, m)
    aS <- alpha_sidak(alpha, m)

    if (!independence) return(sprintf("Bonferroni (alpha* = %.6f)", aB))

    fwerB <- 1 - (1 - aB)^m
    fwerS <- 1 - (1 - aS)^m
    tol <- 1e-12
    okB <- (fwerB <= alpha + tol)
    okS <- (fwerS <= alpha + tol)

    if (okB && okS) {
      if (abs(alpha - fwerS) <= abs(alpha - fwerB)) {
        return(sprintf("Sidak (indep.) (alpha* = %.6f)", aS))
      } else {
        return(sprintf("Bonferroni (alpha* = %.6f)", aB))
      }
    } else if (okS) {
      return(sprintf("Sidak (indep.) (alpha* = %.6f)", aS))
    } else if (okB) {
      return(sprintf("Bonferroni (alpha* = %.6f)", aB))
    } else {
      return(sprintf("Bonferroni (alpha* = %.6f)", aB))
    }
  }

  # Recomendador de prueba post hoc (incluye Gabriel, Scheffe, SNK)
  recommend_test <- function(g, m, equal_var, unequal_n, scope, liberal_ok){
    if (!equal_var) {
      return("Games-Howell (heterocedastic) or Holm (robust)")
    }
    if (unequal_n) {
      if (g <= 12) return("Gabriel (moderate imbalance) or Tukey-Kramer / Holm")
      return("Holm (large m) or Tukey-Kramer")
    }
    # balanceado y homocedastico
    if (scope == "cells" && m > 45) {
      return("Holm (powerful) or Scheffe (many contrasts)")
    }
    if (m <= 6) {
      if (liberal_ok) return("SNK or LSD (with F for gate) | Alternate: Tukey HSD / Holm")
      return("Tukey HSD or Holm (LSD/SNK if you accept major Type I)")
    }
    if (m <= 20) {
      if (liberal_ok) return("SNK or Duncan (more powerful) | Alternative: Tukey HSD / Holm")
      return("Holm (step-down) o Tukey HSD")
    }
    if (m <= 45) return("Holm (step-down) or Tukey HSD")
    "Holm (step-down) or Scheffe (if there are more than pairs)"
  }

  # -------- Conteo de niveles / celdas --------
  level_count <- function(varname) nlevels(mf[[varname]])

  build_row <- function(label, g, alpha, equal_var, unequal_n, scope, independence, liberal_ok){
    m  <- choose2(g)
    aB <- alpha_bonf(alpha, m)
    aS <- alpha_sidak(alpha, m)
    data.frame(
      Unit = label,
      `N levels (g)` = g,
      `Comparasions (m)` = m,
      `alpha global (FWER)` = alpha,
      `alpha_B = alpha/m (Bonf)` = aB,
      `FWER Bonf (indep.)` = if (independence) fwer_indep(aB, m) else NA_real_,
      `alpha* Sidak` = aS,
      `FWER Sidak (indep.)` = if (independence) fwer_indep(aS, m) else NA_real_,
      `FWER without (indep.)` = if (independence) fwer_indep(alpha, m) else NA_real_,
      `Suggested p-value adj` = choose_adjust(alpha, m, independence),
      `Suggestion post hoc` = recommend_test(g, m, equal_var, unequal_n, scope, liberal_ok),
      check.names = FALSE
    )
  }

  res_list <- list()

  if (scope == "factor") {
    for (nm in compare) {
      g <- level_count(nm)
      res_list[[length(res_list)+1]] <-
        build_row(paste0("Factor: ", nm), g, alpha, equal_var, unequal_n, scope, independence, liberal_ok)
    }
  } else { # "cells": combinaciin completa de los factores provistos
    if (length(compare) < 2) {
      g_vec <- vapply(compare, level_count, FUN.VALUE = numeric(1))
      G <- prod(g_vec)
    } else {
      if (observed_cells) {
        inter <- interaction(mf[compare], drop = TRUE, sep = ":")
        G <- nlevels(inter)
      } else {
        g_vec <- vapply(compare, level_count, FUN.VALUE = numeric(1))
        G <- prod(g_vec)
      }
    }
    lbl <- paste0("Celdas (", paste(compare, collapse=":"), ")")
    res_list[[1]] <- build_row(lbl, G, alpha, equal_var, unequal_n, scope, independence, liberal_ok)
  }

  res_wide <- do.call(rbind, res_list)
  rownames(res_wide) <- NULL

  # Orden de columnas
  col_order <- c(
    "Unit","N levels (g)","Comparasions (m)","alpha global (FWER)",
    "alpha_B = alpha/m (Bonf)","FWER Bonf (indep.)",
    "alpha* Sidak","FWER Sidak (indep.)",
    "FWER without adj (indep.)","Suggested p-value adj","Suggestion post hoc"
  )
  col_order <- col_order[col_order %in% names(res_wide)]
  res_wide <- res_wide[, col_order, drop = FALSE]

  if (orientation == "cols") {
    attr(res_wide, "note") <- "FWER under independence is indicative; with actual dependence it may differ."
    return(res_wide)
  }

  # -------- Presentaciin por filas --------
  is_num <- vapply(res_wide, is.numeric, logical(1))
  fmt_num <- function(x) ifelse(is.na(x), NA, formatC(x, format = "f", digits = digits))
  fmt_pct <- function(x) ifelse(is.na(x), NA, paste0(formatC(100*x, format = "f", digits = percent_digits), "%"))

  pct_cols <- intersect(
    c("FWER Bonf (indep.)","FWER Sidak (indep.)","FWER without adj (indep.)"),
    names(res_wide)
  )

  res_fmt <- res_wide
  for (cl in names(res_fmt)) {
    if (cl %in% pct_cols) {
      res_fmt[[cl]] <- fmt_pct(res_fmt[[cl]])
    } else if (isTRUE(is_num[[cl]]) && !cl %in% c("N levels (g)", "Comparasions (m)")) {
      res_fmt[[cl]] <- fmt_num(res_fmt[[cl]])
    }
  }

  unidades <- res_fmt$Unidad
  metricas <- setdiff(names(res_fmt), "Unit")

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

  # -------- Reordenamiento robusto (sin perder filas) --------
  conocidas <- c("N levels (g)", "Comparasions (m)", "alpha global (FWER)",
                 "alpha_B = alpha/m (Bonf)", "FWER Bonf (indep.)",
                 "alpha* Sidak", "FWER Sidak (indep.)",
                 "FWER without adj (indep.)", "Suggested p-value adj", "Suggestion post hoc")
  extra <- setdiff(tabla_rows$Metrica, conocidas)
  orden_final <- c(conocidas, extra)
  idx <- match(orden_final, tabla_rows$Metrica)
  idx <- idx[!is.na(idx)]
  tabla_rows <- tabla_rows[idx, , drop = FALSE]

  attr(tabla_rows, "note") <- "FWER under independence is indicative; with dependent tests it may not match."
  tabla_rows
}
