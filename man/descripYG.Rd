% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descripYG.R
\name{descripYG}
\alias{descripYG}
\title{Descriptive Analysis With Optional Grouping}
\usage{
descripYG(dataset, vd, vi = NULL)
}
\arguments{
\item{dataset}{A \code{data.frame} or \code{tibble} containing the variables.}

\item{vd}{A numeric variable to analyze (dependent variable).}

\item{vi}{An optional grouping variable (independent variable, categorical or numeric).}
}
\value{
A \code{data.frame} with descriptive statistics. Also prints plots to the graphics device.
}
\description{
Performs a descriptive analysis on a numeric dependent variable, either globally
or grouped by an independent variable. Displays summary statistics such as mean,
standard deviation, skewness, and kurtosis, and generates associated plots
(histogram, boxplot, or density ridges).
}
\examples{
data(d_e, package = "Analitica")
descripYG(d_e, vd = Sueldo_actual)
descripYG(d_e, vd = Sueldo_actual, vi = labor)
descripYG(d_e,Sueldo_actual,labor)


}
