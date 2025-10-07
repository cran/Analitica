#' Mauchly's Test for Sphericity (1 factor) v1.0
#'
#' Tests for sphericity in repeated measures designs. Uses an interface for
#' type formula \code{dv ~ within | id}, where:
#' \itemize{
#'  \item \code{dv}: variable numerical response,
#'  \item \code{within}: within-subjects factor (repeated levels),
#'  \item \code{id}: subject/sample identifier.
#' }
#'
#' Calculates Mauchly's statistic \eqn{W}, its approximation \eqn{W}
#' corrected, the p-value and the correction coefficients for lack of
#' sphericity (Greenhouse–Geisser y Huynh–Feldt).
#'
#' @references Mauchly, J. W. (1940). Significance test for sphericity of a normal n-variate distribution. The Annals of Mathematical Statistics, 11(2), 204–209. https://doi.org/10.1214/aoms/1177731915
#'
#' @param formula Formula \code{dv ~ within | id}.
#' @param data \code{data.frame} with the variables.
#' @param alpha Significance level (default 0.05).
#' @param digits Decimals for printing (default 4).
#' @param do_print if \code{TRUE}, print a friendly summary.
#'
#' @return Objeto de clase \code{"sphericity"} with:
#' \describe{
#'   \item{Method}{Cadena con el método.}
#'   \item{Statistic}{Lista con \code{W} y \code{Chi2}.}
#'   \item{df}{Grados de libertad.}
#'   \item{p_value}{Valor-p.}
#'   \item{Decision}{\code{"Sphericity"} o \code{"No sphericity"} según \code{alpha}.}
#'   \item{Epsilons}{\code{GG} and \code{HF}.}
#'   \item{Components}{List with \code{n} (subjects), \code{k} (levels), \code{S} (covariances), \code{eigen} (eigenvalues).}
#' }
#'
#' @importFrom stats model.frame terms model.response as.formula na.omit
#' @importFrom stats cov contr.helmert pchisq reshape
#'
#'
#' @examples
#' # Ejemplo mínimo (datos ficticios):
#' set.seed(1)
#' d <- data.frame(
#'   id = rep(1:10, each = 4),
#'   within = rep(paste0("t", 1:4), times = 10),
#'   y = as.numeric(rep(rnorm(10, 10, 2), each = 4)) +
#'       rep(c(0, .5, 1.2, .8), times = 10) + rnorm(40, 0, 1)
#' )
#' res <- MauchlyTest(y ~ within | id, data = d, do_print = TRUE)
#' summary(res)
#' @export
MauchlyTest <- function(formula, data, alpha = 0.05, digits = 4, do_print = TRUE) {

  # ---- helpers ----
  .parse_rm_formula <- function(f) {
    # Expect "y ~ within | id"
    tf <- deparse(formula)
    # split by "~"
    parts <- strsplit(tf, "~", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("Formula must be of the form: dv ~ within | id")
    lhs <- trimws(parts[1])
    rhs <- trimws(parts[2])
    # split RHS by "|"
    rhs_parts <- strsplit(rhs, "|", fixed = TRUE)[[1]]
    if (length(rhs_parts) != 2) stop("Right-hand side must be 'within | id'.")
    within_expr <- trimws(rhs_parts[1])
    id_expr     <- trimws(rhs_parts[2])
    list(lhs = lhs, within = within_expr, id = id_expr)
  }

  .stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)

  # ---- 0) Preparación ----
  .stop_if(missing(formula) || missing(data), "Provide 'formula' and 'data'.")
  pieces <- .parse_rm_formula(formula)

  mf <- stats::model.frame(
    stats::as.formula(paste(pieces$lhs, "~", pieces$within, "+", pieces$id)),
    data = data, na.action = stats::na.omit
  )

  y      <- stats::model.response(mf)
  within <- mf[[pieces$within]]
  id     <- mf[[pieces$id]]

  .stop_if(!is.numeric(y), "The response (dv) must be numeric.")

  within <- as.factor(within)
  id     <- as.factor(id)

  # Pasar a formato ancho: una columna por nivel de 'within'
  if (requireNamespace("tidyr", quietly = TRUE)) {
    wide <- tidyr::pivot_wider(
      data.frame(id = id, within = within, y = y),
      id_cols = "id", names_from = "within", values_from = "y"
    )
  } else {
    wide <- reshape(
      data.frame(id = id, within = within, y = y),
      idvar = "id", timevar = "within", direction = "wide"
    )
    # renombrar columnas y.*
    old <- setdiff(names(wide), "id")
    names(wide)[names(wide) %in% old] <- gsub("^y\\.", "", old)
  }

  mat <- as.matrix(wide[, setdiff(names(wide), "id"), drop = FALSE])
  .stop_if(anyNA(mat), "Missing data per subject/within level: Mauchly requires a complete design.")

  n <- nrow(mat)         # sujetos
  k <- ncol(mat)         # niveles (k >= 3)
  .stop_if(k < 3, "Mauchly requires at least 3 within levels (k >= 3).")

  # ---- 1) Matriz de covarianzas entre niveles ----
  S <- stats::cov(mat)

  # ---- 2) Contrastes ortonormales (Helmert) ----
  L <- stats::contr.helmert(k)  # k x (k-1)
  H <- t(L)                     # (k-1) x k
  H <- diag(1 / sqrt(rowSums(H^2))) %*% H  # ortonormalizar filas

  # ---- 3) Subespacio y autovalores ----
  M  <- H %*% S %*% t(H)                       # (k-1) x (k-1)
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values

  # ---- 4) Estadístico W ----
  W <- prod(ev) / ((sum(ev) / (k - 1))^(k - 1))

  # ---- 5) Chi-cuadrado corregido ----
  df <- (k - 1) * (k - 2) / 2
  u  <- (2 * k^2 + k + 2) / (6 * (k - 1) * (n - 1))
  X2 <- -(n - 1) * (1 - u) * log(W)
  p  <- stats::pchisq(X2, df = df, lower.tail = FALSE)

  # ---- 6) Epsilons GG y HF ----
  epsGG    <- (sum(ev)^2) / ((k - 1) * sum(ev^2))
  epsHFraw <- (n * (k - 1) * epsGG - 2) / ((k - 1) * (n - 1) - (k - 1) * epsGG)
  epsHF    <- max(min(epsHFraw, 1), epsGG)

  decision <- if (p < alpha) "No sphericity" else "Sphericity"

  out <- list(
    Method    = "Mauchly (sphericity) - manual",
    Statistic = list(W = as.numeric(W), Chi2 = as.numeric(X2)),
    df        = df,
    p_value   = as.numeric(p),
    alpha     = alpha,
    Decision  = decision,
    Epsilons  = c(GG = as.numeric(epsGG), HF = as.numeric(epsHF)),
    Components= list(n = n, k = k, S = S, eigen = ev, u = u)
  )
  class(out) <- "esfericidad"

  if (isTRUE(do_print)) {
    cat(
      sprintf("Mauchly (sphericity): dv ~ %s | %s\n", pieces$within, pieces$id),
      sprintf("Subjects: n = %d | Levels: k = %d\n", n, k),
      sprintf("W = %.*f | Chi^2 = %.*f | gl = %d | p = %.*f --> %s\n",
              digits, W, digits, X2, df, digits, p, decision),
      sprintf("Epsilons: GG = %.*f | HF = %.*f\n", digits, epsGG, digits, epsHF),
      if (p < alpha)
        "Recommendation: apply GG/HF correction for within-subject effects.\n"
      else
        "Recommendation: within-subject effects can be reported without sphericity correction.\n",
      sep = ""
    )
  }

  invisible(out)
}

# ---- Métodos S3 ----
#' @export
#' @method print esfericidad
print.esfericidad <- function(x, digits = 4, ...) {
  cat("Mauchly's Test for Sphericity (manual)\n")
  cat(sprintf("W = %.*f | Chi^2 = %.*f | gl = %d | p = %.*f\n",
              digits, x$Statistic$W, digits, x$Statistic$Chi2, x$df, digits, x$p_value))
  cat(sprintf("Decision (@ alpha=%.3f): %s\n", x$alpha, x$Decision))
  cat(sprintf("Epsilons: GG = %.*f | HF = %.*f\n",
              digits, x$Epsilons["GG"], digits, x$Epsilons["HF"]))
  invisible(x)
}

#' @export
#' @method summary esfericidad
summary.esfericidad <- function(object, digits = 4, ...) {
  cat("== Mauchly's Test for Sphericity ==\n")
  cat(sprintf("Method: %s\n", object$Method))
  cat(sprintf("n (subjects) = %d | k (levels) = %d\n", object$Components$n, object$Components$k))
  cat(sprintf("W = %.*f\n", digits, object$Statistic$W))
  cat(sprintf("Chi^2 = %.*f | df = %d | p = %.*f\n",
              digits, object$Statistic$Chi2, object$df, digits, object$p_value))
  cat(sprintf("Epsilons -> GG: %.*f | HF: %.*f\n",
              digits, object$Epsilons["GG"], digits, object$Epsilons["HF"]))
  cat("----\nEigenvalues (H S H^T):\n")
  print(round(object$Components$eigen, digits))
  invisible(object)
}

# USO con tus nombres:
# mauchly_manual(datos, id = "muestra", within = "orden", dv = "Pconc")
