#' Bar Plot with Error Bars (Standard Deviation or Standard Error)
#'
#' Creates a bar plot of group means with error bars representing either
#' the standard deviation (SD) or the standard error (SE).
#'
#' @param dataSet A \code{data.frame} or \code{tibble} containing the data.
#' @param vD A string indicating the name of the numeric dependent variable.
#' @param vI A string indicating the name of the categorical independent variable (grouping variable).
#' @param variation Type of variation to display: \code{"sd"} for standard deviation
#'   or \code{"se"} for standard error. Default is \code{"sd"}.
#' @param title Title of the plot. Default is \code{"Bar plot with error bars"}.
#' @param label_y Label for the Y-axis. Default is \code{"Y Axis"}.
#' @param label_x Label for the X-axis. Default is \code{"X Axis"}.
#'
#' @return A \code{ggplot} object representing the plot.
#'
#' @examples
#' data(d_e, package = "Analitica")
#' bar_error(d_e, vD = Sueldo_actual, vI = labor, variation = "sd")
#'
#' @export
#' @importFrom dplyr group_by summarize n
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar geom_text labs theme_minimal theme element_text
#' @importFrom rlang ensym
bar_error <- function(dataSet, vD, vI, variation = "sd",
                      title = "Bar plot with error bars",
                      label_y = "Y Axis", label_x = "X Axis") {

  if (!variation %in% c("sd", "se")) {
    stop("The 'variation' argument must be either 'sd' or 'se'.")
  }

  vD_quo <- rlang::enquo(vD)
  vI_quo <- rlang::enquo(vI)

  vD_name <- rlang::as_name(vD_quo)
  vI_name <- rlang::as_name(vI_quo)

  if (!(vD_name %in% names(dataSet))) {
    stop(paste("Dependent variable not found:", vD_name))
  }
  if (!(vI_name %in% names(dataSet))) {
    stop(paste("Independent variable not found:", vI_name))
  }

  summary_data <- dataSet %>%
    dplyr::group_by(!!vI_quo) %>%
    dplyr::summarize(
      mean = mean(!!vD_quo, na.rm = TRUE),
      sd = sd(!!vD_quo, na.rm = TRUE),
      n = dplyr::n(),
      se = sd / sqrt(n),
      .groups = "drop"
    )

  summary_data$variation_value <- if (variation == "sd") summary_data$sd else summary_data$se

  ggplot2::ggplot(summary_data, ggplot2::aes(x = factor(!!vI_quo), y = mean)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - variation_value, ymax = mean + variation_value),
                           width = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = round(mean, 2)), vjust = -0.5, size = 3.5) +
    ggplot2::labs(title = title, x = label_x, y = label_y) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
}
