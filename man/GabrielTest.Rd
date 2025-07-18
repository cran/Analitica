% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GabrielTest.R
\name{GabrielTest}
\alias{GabrielTest}
\title{Gabriel’s Post Hoc Test for Multiple Comparisons}
\usage{
GabrielTest(modelo, alpha = 0.05)
}
\arguments{
\item{modelo}{An object of class \code{aov} or \code{lm}.}

\item{alpha}{Significance level (default is 0.05).}
}
\value{
An object of class \code{"gabriel"} and \code{"comparaciones"}, containing:
\itemize{
  \item \code{Resultados}: Data frame with comparisons, mean differences, adjusted critical value, p-value, and significance level.
  \item \code{Promedios}: Named numeric vector of group means.
  \item \code{Orden_Medias}: Vector of group names ordered from highest to lowest mean.
  \item \code{Metodo}: Name of the method used ("Gabriel").
}
}
\description{
A modification of Tukey's test for use with moderately unequal sample sizes.
}
\details{
Advantages:
- More powerful than Tukey for unequal group sizes.
- Controls error rates effectively with moderate imbalance.

Disadvantages:
- Can be anti-conservative with large differences in group sizes.
- Less common in standard statistical software.
}
\examples{
data(d_e, package = "Analitica")
mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
resultado <- GabrielTest(mod)
summary(resultado)
plot(resultado)


}
\references{
Hochberg, Y., & Tamhane, A. C. (1987). Multiple Comparison Procedures.
}
