#' Outlier Detection Using Grubbs' Test (Iterative)
#'
#' Detects one or more outliers in a numeric variable using the iterative Grubbs' test,
#' which assumes the data follow a normal distribution.
#'
#' @param dataSet A \code{data.frame} containing the data.
#' @param vD Unquoted name of the numeric variable to be tested for outliers.
#' @param alpha Significance level for the test (default is \code{0.05}).
#'
#' @return A \code{data.frame} identical to the input, with an added logical column \code{outL}
#' indicating which observations were identified as outliers (\code{TRUE} or \code{FALSE}).
#'
#' @details
#' The function applies Grubbs' test iteratively, removing the most extreme value
#' and retesting until no further significant outliers are found.
#' The test is valid only under the assumption of normality.
#'
#' @references Grubbs, F. E. (1969). "Procedures for Detecting Outlying Observations in Samples."
#' \emph{Technometrics}, 11(1), 1–21. \doi{10.1080/00401706.1969.10490657}
#'
#' @export
#' @importFrom stats qt
#' @importFrom rlang enquo as_name
#'
#' @examples
#' data(d_e, package = "Analitica")
#' d<-grubbs_outliers(d_e, Sueldo_actual)
grubbs_outliers <- function(dataSet, vD, alpha = 0.05) {
  vD_quo <- rlang::enquo(vD)
  vD_name <- rlang::as_name(vD_quo)

  if (!vD_name %in% names(dataSet)) {
    stop(paste("Variable", vD_name, "not found in the dataset."))
  }

  valores <- dataSet[[vD_name]]

  if (!is.numeric(valores)) {
    stop("The variable must be numeric.")
  }

  dataSet$outL <- FALSE
  indices <- seq_along(valores)

  grubbs_stat <- function(x, alpha) {
    n <- length(x)
    if (n < 3) return(list(G = 0, G_crit = Inf))  # Minimum required sample size
    G <- max(abs(x - mean(x))) / sd(x)
    t_crit <- qt(1 - alpha / (2 * n), df = n - 2)
    G_crit <- ((n - 1) / sqrt(n)) * sqrt(t_crit^2 / (n - 2 + t_crit^2))
    list(G = G, G_crit = G_crit)
  }

  repeat {
    g <- grubbs_stat(valores, alpha)
    if (g$G > g$G_crit) {
      idx_out <- which.max(abs(valores - mean(valores)))
      idx_original <- indices[idx_out]
      dataSet$outL[idx_original] <- TRUE
      valores <- valores[-idx_out]
      indices <- indices[-idx_out]
    } else {
      break
    }
  }

  return(dataSet)
}

