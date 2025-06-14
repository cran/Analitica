## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
library(Analitica)
data(d_e, package = "Analitica")

## ----desc-example-------------------------------------------------------------
data(d_e, package = "Analitica")
descripYG(d_e, vd = Sueldo_actual)
descripYG(d_e, vd = Sueldo_actual, vi = labor)

## ----homo-tests---------------------------------------------------------------
Levene.Test(Sueldo_actual ~ labor, data = d_e)
BartlettTest(Sueldo_actual ~ labor, data = d_e)
FKTest(Sueldo_actual ~ labor, data = d_e)

## ----outliers-----------------------------------------------------------------
res <- grubbs_outliers(d_e, Sueldo_actual)
head(res[res$outL == TRUE, ])

## ----comparisons--------------------------------------------------------------
mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
resultado <- GHTest(mod)
summary(resultado)
plot(resultado)

## ----np-tests-----------------------------------------------------------------
g1 <- d_e$Sueldo_actual[d_e$labor == 1]
g2 <- d_e$Sueldo_actual[d_e$labor == 2]
MWTest(g1, g2)
BMTest(g1, g2)
BMpTest(g1, g2)

