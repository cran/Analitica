#' Descriptive Analysis With Optional Grouping
#'
#' Performs a descriptive analysis on a numeric dependent variable, either globally
#' or grouped by an independent variable. Displays summary statistics such as mean,
#' standard deviation, skewness, and kurtosis, and generates associated plots
#' (histogram, boxplot, or density ridges).
#'
#' @param dataset A \code{data.frame} or \code{tibble} containing the variables.
#' @param vd A numeric variable to analyze (dependent variable).
#' @param vi An optional grouping variable (independent variable, categorical or numeric).
#'
#' @return A \code{data.frame} with descriptive statistics. Also prints plots to the graphics device.
#' @export
#'
#' @importFrom dplyr group_by summarise rename n across tibble
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom stats IQR median quantile sd
#' @importFrom ggplot2 ggplot aes geom_histogram geom_boxplot
#'             theme_minimal theme_classic xlab ylab labs
#' @importFrom ggridges geom_density_ridges2
#' @importFrom moments kurtosis skewness
#' @importFrom rlang sym
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' data(d_e, package = "Analitica")
#' descripYG(d_e, vd = Sueldo_actual)
#' descripYG(d_e, vd = Sueldo_actual, vi = labor)
#' descripYG(d_e,Sueldo_actual,labor)
#'
#'
descripYG <- function(dataset, vd, vi = NULL) {

  vd_name <- as.character(substitute(vd))
  vi_name <- if (!missing(vi)) as.character(substitute(vi)) else NULL

  if (!vd_name %in% names(dataset)) stop(paste("Dependent variable not found:", vd_name))
  if (!is.null(vi_name) && !vi_name %in% names(dataset)) stop(paste("Independent variable not found:", vi_name))

  dataset <- dataset %>% tidyr::drop_na(tidyselect::all_of(c(vd_name, vi_name)))

  vd_sym <- rlang::sym(vd_name)

  if (is.null(vi_name)) {
    vd_vec <- dataset[[vd_name]]
    IQR_val <- IQR(vd_vec)
    summary_table <- dplyr::tibble(
      n = length(vd_vec),
      Mean = mean(vd_vec),
      Median = median(vd_vec),
      SD = sd(vd_vec),
      Kurtosis = moments::kurtosis(vd_vec),
      Skewness = moments::skewness(vd_vec),
      CV = sd(vd_vec) / abs(mean(vd_vec)),
      Min = min(vd_vec),
      Max = max(vd_vec),
      P25 = quantile(vd_vec, 0.25),
      P75 = quantile(vd_vec, 0.75),
      IQR = IQR_val,
      Fence_Low = quantile(vd_vec, 0.25) - 1.5 * IQR_val,
      Fence_High = quantile(vd_vec, 0.75) + 1.5 * IQR_val
    )

    bins <- trunc(3.322 * log10(length(vd_vec)) + 1)

    g_hist <- ggplot2::ggplot(dataset, ggplot2::aes(x = !!vd_sym)) +
      ggplot2::geom_histogram(color = 'white', fill = 'steelblue', bins = bins) +
      ggplot2::theme_classic()

    g_box <- ggplot2::ggplot(dataset, ggplot2::aes(x = !!vd_sym)) +
      ggplot2::geom_boxplot(fill = "lightgreen", width = 0.05)

    print(patchwork::wrap_plots(g_hist, g_box, ncol = 1, heights = c(3, 1)))

    return(as.data.frame(summary_table))

  } else {
    vi_sym <- rlang::sym(vi_name)

    summary_grouped <- dataset %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(vi_name))) %>%
      dplyr::summarise(
        n = dplyr::n(),
        Mean = mean(!!vd_sym),
        Median = median(!!vd_sym),
        SD = sd(!!vd_sym),
        Kurtosis = moments::kurtosis(!!vd_sym),
        Skewness = moments::skewness(!!vd_sym),
        CV = SD / abs(Mean),
        Min = min(!!vd_sym),
        Max = max(!!vd_sym),
        P25 = quantile(!!vd_sym, 0.25),
        P75 = quantile(!!vd_sym, 0.75),
        IQR = IQR(!!vd_sym),
        .groups = "drop"
      ) %>%
      dplyr::rename(Group = tidyselect::all_of(vi_name))

    g_ridge <- ggplot2::ggplot(dataset, ggplot2::aes(x = !!vd_sym, y = as.factor(!!vi_sym))) +
      ggridges::geom_density_ridges2(fill = "steelblue", alpha = 0.7) +
      ggplot2::theme_minimal()

    g_box <- ggplot2::ggplot(dataset, ggplot2::aes(x = factor(!!vi_sym), y = !!vd_sym)) +
      ggplot2::geom_boxplot(color = 'darkslategray', fill = 'steelblue') +
      ggplot2::theme_classic()

    print(g_ridge)
    print(g_box)

    return(as.data.frame(summary_grouped))
  }
}


