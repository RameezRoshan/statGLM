if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".data", "density"))
}

#' Modern Publication-Ready Theme for statGLM
#'
#' @description
#' Provides a clean, elegant, high-definition ggplot2 theme tailored for
#' statistical reports and academic publications (e.g. Nature, Science, JASA, Biometrika).
#' Features crisp vector borders, subtle gridlines, legible typography, and seamless
#' multi-panel alignment.
#'
#' @param base_size Base font size in points (default 11).
#' @param base_family Base font family (default \code{""}).
#' @param grid Character indicating grid lines: \code{"major"} (default), \code{"both"}, \code{"x"}, \code{"y"}, or \code{"none"}.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_statglm()
#'
#' @export
theme_statglm <- function(base_size = 11, base_family = "", grid = c("major", "both", "x", "y", "none")) {
  grid <- match.arg(grid)
  
  th <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Panel background and borders
      panel.background = ggplot2::element_rect(fill = "#ffffff", color = NA),
      plot.background = ggplot2::element_rect(fill = "#ffffff", color = NA),
      panel.border = ggplot2::element_rect(color = "#333333", fill = NA, linewidth = 0.6),
      
      # Titles & Typography
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1.15),
        face = "bold",
        color = "#111111",
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(0.92),
        color = "#555555",
        margin = ggplot2::margin(b = 8)
      ),
      plot.caption = ggplot2::element_text(
        size = ggplot2::rel(0.8),
        color = "#666666",
        hjust = 1,
        margin = ggplot2::margin(t = 6)
      ),
      
      # Axis formatting
      axis.title = ggplot2::element_text(size = ggplot2::rel(0.95), face = "bold", color = "#222222"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(0.85), color = "#333333"),
      axis.ticks = ggplot2::element_line(color = "#333333", linewidth = 0.5),
      
      # Legend formatting
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9), face = "bold", color = "#222222"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.82), color = "#333333"),
      legend.position = "right",
      
      # Facet strips
      strip.background = ggplot2::element_rect(fill = "#f2f4f8", color = "#cccccc", linewidth = 0.5),
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.9), face = "bold", color = "#222222",
                                         margin = ggplot2::margin(t = 4, b = 4))
    )
  
  # Configure grid lines
  if (grid == "none") {
    th <- th + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank())
  } else if (grid == "major") {
    th <- th + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "#f0f0f0", linewidth = 0.4),
                              panel.grid.minor = ggplot2::element_blank())
  } else if (grid == "x") {
    th <- th + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = "#f0f0f0", linewidth = 0.4),
                              panel.grid.major.y = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank())
  } else if (grid == "y") {
    th <- th + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "#f0f0f0", linewidth = 0.4),
                              panel.grid.major.x = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank())
  } else if (grid == "both") {
    th <- th + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "#f0f0f0", linewidth = 0.4),
                              panel.grid.minor = ggplot2::element_line(color = "#f8f8f8", linewidth = 0.25))
  }
  
  th
}

#' Save statGLM Plots in High Definition for Publication
#'
#' @description
#' Exports statGLM plots to vector (PDF, SVG) or high-resolution raster (PNG, TIFF)
#' formats conforming to publication standards (300 to 600 DPI).
#'
#' @param plot A ggplot or gtable / grob object to save.
#' @param filename File name with extension (e.g. \code{"figure1.png"} or \code{"figure1.pdf"}).
#' @param width Width in inches (default 7).
#' @param height Height in inches (default 5).
#' @param dpi Dots per inch for raster outputs (default 300, publication grade).
#' @param ... Additional arguments passed to \code{ggplot2::ggsave}.
#'
#' @return The filename string invisibly.
#'
#' @export
save_glm_plot <- function(plot, filename, width = 7, height = 5, dpi = 300, ...) {
  if (inherits(plot, "gtable") || inherits(plot, "grob")) {
    ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi, ...)
  } else {
    ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
  invisible(filename)
}
