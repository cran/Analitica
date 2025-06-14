#' Analitica: Tools for Exploratory Data Analysis and Group Comparisons
#'
#' The \strong{Analitica} package provides tools for exploratory statistical analysis,
#' data visualization, and comparison of groups using both parametric and non-parametric methods.
#' It supports univariate and grouped descriptive summaries, outlier detection, homoscedasticity testing,
#' and multiple post hoc procedures.
#'
#' Designed for applied analysis workflows, this package includes intuitive plotting functions
#' and manual implementations of key statistical tests often needed in educational or research contexts.
#'
#' @section Main Features:
#' \itemize{
#'   \item \code{\link{descripYG}}: Descriptive statistics with visualizations (histograms, boxplots, density ridges).
#'   \item \code{\link{Levene.Test}}: Manual implementation of Levene’s test for homogeneity of variances.
#'   \item \code{\link{BartlettTest}}: Manual implementation of Bartlett’s test.
#'   \item \code{\link{FKTest}}: Manual implementation of the Fligner-Killeen test.
#'   \item \code{\link{grubbs_outliers}}: Outlier detection based on Grubbs' test.
#'   \item \code{\link{GHTest}}, \code{\link{DuncanTest}}, \code{\link{SNKTest}}, etc.: Post hoc comparison procedures.
#' }
#'
#' @author Carlos Jiménez-Gallardo
#' @keywords package
#' @name Analitica
NULL
