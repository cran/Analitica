#' Levene's Test for Homogeneity of Variances (Manual Implementation) v2.1
#'
#' Performs Levene's test for equality of variances across groups using a formula interface.
#' This test evaluates the null hypothesis that the variances are equal across groups,
#' and is commonly used as a preliminary test before ANOVA or other parametric analyses.
#'
#' Levene's test is based on an analysis of variance (ANOVA) applied to the absolute deviations
#' from each group's center (either the mean or, more robustly, the median).
#' It is less sensitive to departures from normality than Bartlett's test.
#'
#' Advantages:
#' - Robust to non-normality, especially when using the median.
#' - Suitable for equal or unequal sample sizes across groups.
#' - Widely used in practice for checking homoscedasticity.
#'
#' Disadvantages:
#' - Less powerful than parametric alternatives under strict normality.
#'
#' @references
#' Levene, H. (1960). "Robust Tests for Equality of Variances." In Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling (pp.278-292). Stanford University Press.
#'
#' @param formula y ~ factors (e.g., y ~ A or y ~ A * B).
#' @param data data.frame with variables in the formula.
#' @param alpha Significance level (default 0.05).
#' @param center "median" (Brown-Forsythe, default) or "mean" (classical Levene).
#' @param decompose logical. If TRUE and there are >= 2 factors, run ANOVA on |Y - cell_center|.
#' @param anova_type "I", "II" (2-way only) or "III" (any number of factors, no 'car').
#'
#' @return An object of class \code{"homocedasticidad"}, containing:
#' \describe{
#'   \item{Statistic}{F statistic of the Levene test.}
#'   \item{df}{Degrees of freedom (between and within groups).}
#'   \item{p_value}{The p-value for the test.}
#'   \item{Decision}{\code{"Homoscedastic"} or \code{"Heteroscedastic"} depending on the test result.}
#'   \item{Method}{A string indicating the method used ("Levene").}
#' }
#'
#' @export
#' @importFrom stats model.frame na.omit terms delete.response model.response update ave pf anova lm drop1
#'
#' @examples
#' data(d_e, package = "Analitica")
#' res <- Levene.Test(Sueldo_actual ~ as.factor(labor), data = d_e)
#' summary(res)
#'
#' # RCBD
#' resB<-Levene.Test(Sueldo_actual ~ as.factor(labor)+Sexo, data = d_e)
#' summary(resB)
#'
#' # anova 2-ways
#' resC<-Levene.Test(Sueldo_actual ~ as.factor(labor)*Sexo, data = d_e)
#' summary(resC)
#'
Levene.Test <- function(formula, data, alpha = 0.05,
                        center = "median",
                        decompose = TRUE,
                        anova_type = c("I","II","III")) {

  anova_type <- match.arg(anova_type)
  if (missing(formula) || missing(data)) stop("Both 'formula' and 'data' must be provided.")
  if (!center %in% c("median","mean")) stop("Argument 'center' must be 'median' or 'mean'.")

  # 1) Model frame (no NA) and evaluated terms (supports as.factor(), etc.)
  mf    <- stats::model.frame(formula, data, na.action = stats::na.omit)
  tt_mf <- stats::terms(mf)
  if (attr(tt_mf, "response") != 1) stop("Formula must have a response (y ~ ...).")

  y <- stats::model.response(mf)
  if (!is.numeric(y)) stop("The response must be numeric.")

  # RHS terms as they appear in mf (expanded labels; may include colnames like `as.factor(A)`)
  rhs_terms_all <- attr(tt_mf, "term.labels")
  if (length(rhs_terms_all) < 1) stop("Provide at least one grouping factor on the RHS.")

  # Main effects only to build CELLS (strip ':' terms)
  base_terms <- rhs_terms_all[!grepl(":", rhs_terms_all)]
  if (length(base_terms) < 1) stop("Need at least one main-effect term on RHS to define cells.")

  rhs_df <- mf[, base_terms, drop = FALSE]
  for (j in seq_along(rhs_df)) rhs_df[[j]] <- as.factor(rhs_df[[j]])

  # CELLS = interaction of main effects (evaluated)
  cell <- interaction(rhs_df, drop = TRUE)

  # 2) Cell center and Z = |Y - cell_center|
  if (center == "mean") {
    cell_center <- stats::ave(y, cell, FUN = mean)
  } else {
    med_by_cell <- tapply(y, cell, median)
    cell_center <- med_by_cell[cell]
  }
  Z <- abs(y - cell_center)

  # 3) Global Levene: 1-way over CELLS
  meanZ_cell <- tapply(Z, cell, mean)
  n_cell     <- tapply(Z, cell, length)
  Zbar       <- mean(Z)

  SS_between <- sum(n_cell * (meanZ_cell - Zbar)^2)
  SS_within  <- sum((Z - meanZ_cell[cell])^2)

  k  <- length(meanZ_cell)
  N  <- length(Z)
  df_between <- k - 1
  df_within  <- N - k
  MS_between <- SS_between / df_between
  MS_within  <- SS_within  / df_within
  F_stat     <- MS_between / MS_within
  p_val      <- 1 - stats::pf(F_stat, df_between, df_within)

  decision <- ifelse(p_val < alpha, "Heteroscedastic (cells)", "Homoscedastic (cells)")
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01,  "**",
                       ifelse(p_val < 0.05,  "*", "ns")))

  # helper to backtick terms containing ':' or parentheses
  backtick_term <- function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    paste0("`", parts, "`", collapse = ":")
  }

  # 4) Factorial decomposition on Z (if >= 2 main effects)
  twoway <- NULL
  if (isTRUE(decompose) && length(base_terms) >= 2) {
    Zdf <- cbind.data.frame(Z = Z, rhs_df)

    # Z formula that respects RHS terms as in mf (includes interactions)
    rhs_bt    <- vapply(rhs_terms_all, backtick_term, character(1))
    z_formula <- stats::as.formula(paste("Z ~", paste(rhs_bt, collapse = " + ")))

    n_factors <- ncol(rhs_df)

    if (anova_type == "I") {
      fit    <- stats::lm(z_formula, data = Zdf)
      aov_tab <- stats::anova(fit)  # Type I (sequential)
      method_str <- sprintf("ANOVA on |Y - cell_center| (%s), SS Type I (sequential)", center)
      p_col <- grep("^Pr\\(>F\\)", colnames(aov_tab), value = TRUE)
      f_col <- grep("^F value$|^F$", colnames(aov_tab), value = TRUE)

    } else if (anova_type == "III") {
      # Type III without 'car': drop1 with F and sum-to-zero contrasts
      op <- options(contrasts = c("contr.sum","contr.poly"))
      on.exit(options(op), add = TRUE)
      fit  <- stats::lm(z_formula, data = Zdf)
      aod  <- stats::drop1(fit, test = "F")
      keep <- rownames(aod) != "<none>"
      aov_tab <- data.frame(
        Df       = aod$Df[keep],
        `Sum Sq` = NA_real_,
        `F value`= aod$`F value`[keep],
        `Pr(>F)` = aod$`Pr(>F)`[keep],
        row.names = rownames(aod)[keep],
        check.names = FALSE
      )
      method_str <- sprintf("ANOVA on |Y - cell_center| (%s), SS Type III (drop1, contr.sum)", center)
      p_col <- "Pr(>F)"; f_col <- "F value"

    } else if (anova_type == "II") {
      # Type II without external packages (implemented for 2-way)
      if (n_factors != 2) {
        # Fallback to Type III for m > 2
        op <- options(contrasts = c("contr.sum","contr.poly"))
        on.exit(options(op), add = TRUE)
        fit  <- stats::lm(z_formula, data = Zdf)
        aod  <- stats::drop1(fit, test = "F")
        keep <- rownames(aod) != "<none>"
        aov_tab <- data.frame(
          Df       = aod$Df[keep],
          `Sum Sq` = NA_real_,
          `F value`= aod$`F value`[keep],
          `Pr(>F)` = aod$`Pr(>F)`[keep],
          row.names = rownames(aod)[keep],
          check.names = FALSE
        )
        method_str <- sprintf("ANOVA on |Y - cell_center| (%s), SS Type III (drop1, contr.sum)", center)
        p_col <- "Pr(>F)"; f_col <- "F value"
      } else {
        # 2-way: nested models
        fA <- names(rhs_df)[1]; fB <- names(rhs_df)[2]
        fit_main <- stats::lm(stats::as.formula(paste("Z ~", fA, "+", fB)), data = Zdf)  # Z ~ A + B
        fit_full <- stats::lm(stats::as.formula(paste("Z ~", fA, "*", fB)), data = Zdf)  # Z ~ A * B

        MS_res_main <- stats::deviance(fit_main) / stats::df.residual(fit_main)
        MS_res_full <- stats::deviance(fit_full) / stats::df.residual(fit_full)

        # A: compare (B) vs (A + B)
        fit_B <- stats::lm(stats::as.formula(paste("Z ~", fB)), data = Zdf)
        SS_A  <- stats::deviance(fit_B) - stats::deviance(fit_main)
        df_A  <- stats::df.residual(fit_B) - stats::df.residual(fit_main)
        F_A   <- (SS_A / df_A) / MS_res_main
        p_A   <- stats::pf(F_A, df_A, stats::df.residual(fit_main), lower.tail = FALSE)

        # B: compare (A) vs (A + B)
        fit_A0 <- stats::lm(stats::as.formula(paste("Z ~", fA)), data = Zdf)
        SS_B   <- stats::deviance(fit_A0) - stats::deviance(fit_main)
        df_B   <- stats::df.residual(fit_A0) - stats::df.residual(fit_main)
        F_B    <- (SS_B / df_B) / MS_res_main
        p_B    <- stats::pf(F_B, df_B, stats::df.residual(fit_main), lower.tail = FALSE)

        # A:B: compare (A + B) vs (A + B + A:B)
        SS_AB <- stats::deviance(fit_main) - stats::deviance(fit_full)
        df_AB <- stats::df.residual(fit_main) - stats::df.residual(fit_full)
        F_AB  <- (SS_AB / df_AB) / MS_res_full
        p_AB  <- stats::pf(F_AB, df_AB, stats::df.residual(fit_full), lower.tail = FALSE)

        aov_tab <- data.frame(
          Df       = c(df_A, df_B, df_AB, stats::df.residual(fit_full)),
          `Sum Sq` = c(SS_A, SS_B, SS_AB, stats::deviance(fit_full)),
          `F value`= c(F_A, F_B, F_AB, NA),
          `Pr(>F)` = c(p_A, p_B, p_AB, NA),
          row.names = c(fA, fB, paste(fA, fB, sep=":"), "Residuals"),
          check.names = FALSE
        )
        method_str <- sprintf("ANOVA on |Y - cell_center| (%s), SS Type II (2-way, no external packages)", center)
        p_col <- "Pr(>F)"; f_col <- "F value"
      }
    }

    # Simple interpretation per term
    interp <- rep(NA_character_, nrow(aov_tab)); names(interp) <- rownames(aov_tab)
    if (!is.null(p_col) && p_col %in% colnames(aov_tab)) {
      pv <- aov_tab[[p_col]]
      interp <- ifelse(is.na(pv) | pv >= alpha, "no evidence", "variance differs")
      names(interp) <- rownames(aov_tab)
    }

    twoway <- list(
      table = aov_tab,
      method = method_str,
      note = NULL,
      interpretation = interp,
      cols = list(df = "Df", F = f_col, p = p_col)
    )
  }

  # 5) Output
  out <- list(
    Method      = paste0("Levene (", center, ") - global by cells"),
    Statistic   = F_stat,
    df          = c(df_between = df_between, df_within = df_within),
    p_value     = p_val,
    Significance= sig,
    Decision    = decision,
    Components  = list(
      SS_between = SS_between, SS_within = SS_within,
      MS_between = MS_between, MS_within = MS_within,
      k_groups = k, N = N, Zbar = Zbar
    ),
    decompose   = !is.null(twoway),
    twoway      = twoway
  )
  class(out) <- "homocedasticidad"
  return(out)
}
