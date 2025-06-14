---
title: "Analitica: Exploratory Data Analysis and Group Comparison Tools"
output: rmarkdown::github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analitica <img src="https://img.shields.io/badge/status-active-brightgreen" alt="Status" align="right"/>

`Analitica` is an R package that provides tools for descriptive statistics, exploratory visualization,
outlier detection, homogeneity of variance tests, and post hoc group comparisons—both parametric
and non-parametric.

It is especially useful for applied analysis, teaching, and reproducible research.

## 📦 Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("<your_github_username>/Analitica")
```

## ✨ Features

- `descripYG()`: Descriptive summaries with histograms, boxplots, or ridge plots
- `Levene.Test()`, `BartlettTest()`, `FKTest()`: Homoscedasticity tests
- `grubbs_outliers()`: Univariate outlier detection using Grubbs' test
- `GHTest()`, `TukeyTest()`, `ScheffeTest()`, `SNKTest()`: Parametric post hoc tests
- `MWTest()`, `BMTest()`, `BMTest_perm()`: Non-parametric group comparisons

## 📊 Example

```r
library(Analitica)
data(d_e, package = "Analitica")

descripYG(d_e, vd = Sueldo_actual, vi = labor)

mod <- aov(Sueldo_actual ~ as.factor(labor), data = d_e)
resultado <- GHTest(mod)
summary(resultado)
plot(resultado)
```

## 📄 License

MIT © [Carlos Jiménez-Gallardo](mailto:carlos.jimenez@ufrontera.cl)

## 📚 Citation

If you use this package, please cite it as:

```text
Jiménez-Gallardo, C. (2025). Analitica: Exploratory Data Analysis and Group Comparison Tools. R package version 1.6.0.

