% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BrunnerMunzelTest.R
\name{BMTest}
\alias{BMTest}
\title{Brunner-Munzel Test for Two Independent Samples}
\usage{
BMTest(
  grupo1,
  grupo2,
  alpha = 0.05,
  alternative = c("two.sided", "less", "greater")
)
}
\arguments{
\item{grupo1}{Numeric vector of values from group 1.}

\item{grupo2}{Numeric vector of values from group 2.}

\item{alpha}{Significance level (default = 0.05).}

\item{alternative}{Character string specifying the alternative hypothesis.
One of \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.}
}
\value{
An object of class \code{"comparacion"} and \code{"brunnermunzel"}, containing:
\itemize{
  \item \code{Resultados}: A data frame with test statistics, p-value, and estimated effect size.
  \item \code{Promedios}: A named numeric vector of group means.
  \item \code{Orden_Medias}: Group names ordered by their mean values (descending).
  \item \code{Metodo}: A character string describing the test and hypothesis.
  \item \code{p_hat}: Estimated probability that a value from grupo1 is less than a value from grupo2 (plus 0.5 * ties).
}
}
\description{
Performs the Brunner-Munzel nonparametric test for two independent groups,
which estimates the probability that a randomly selected value from one group
is less than a randomly selected value from the other group.
}
\details{
This test is suitable when group variances are unequal and/or sample sizes differ.
It does not assume equal variances and is often used as a more robust alternative to the Wilcoxon test.

Advantages:
- Handles unequal variances and non-normality.
- Recommended when variance homogeneity is questionable.

Disadvantages:
- Less well-known and supported.
- Requires large sample sizes for accurate inference.
}
\examples{
data(d_e, package = "Analitica")
g1 <- d_e$Sueldo_actual[d_e$labor == 1]
g2 <- d_e$Sueldo_actual[d_e$labor == 2]
resultado <- BMTest(g1, g2, alternative = "greater")
summary(resultado)

}
\references{
Brunner, E., & Munzel, U. (2000). "The nonparametric Behrens-Fisher problem: Asymptotic theory and a small-sample approximation." \emph{Biometrical Journal}, 42(1), 17–25. <https://doi.org/10.1002/(SICI)1521-4036(200001)42:1%3C17::AID-BIMJ17%3E3.0.CO;2-U>
}
