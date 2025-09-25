utils::globalVariables("hbma_cell_type_palette")

#' ggplot2 fill scale for cell types using HBMA convention
#'
#' @param name Legend title (default = "Cell type")
#' @param ... Additional arguments passed to `scale_fill_manual()`
#'
#' @return A ggplot2 fill scale
#'
#' @importFrom ggplot2 scale_fill_manual
#'
#' @examples
#' library(ggplot2)
#' set.seed(2464)
#' df <- data.frame(
#'   cell_type = c(
#'     "immune", "vascular", "cortex GABAergic", "astrocyte", "ependymal", "OPC",
#'     "oligodendrocyte", "striatum MSN", "striatum non-MSN",
#'     "cortex glutamatergic"
#'   ), stringsAsFactors = FALSE
#' )
#' df$value <- runif(nrow(df))
#'
#' ggplot(df, aes(x = cell_type, y = value, fill = cell_type)) +
#'   geom_col() +
#'   hbma_cell_type_fill_scale()
#'
#' @export
hbma_cell_type_fill_scale <- function(name = "Cell type", ...) {
  scale_fill_manual(name = name, values = get_hbma_cell_type_palette(), ...)
}

#' ggplot2 color scale for cell types using HBMA convention
#'
#' @param name Legend title (default = "Cell type")
#' @param ... Additional arguments passed to `scale_color_manual()`
#'
#' @return A ggplot2 color scale
#'
#' @importFrom ggplot2 scale_color_manual
#'
#' @examples
#' library(ggplot2)
#' set.seed(2464)
#' df <- data.frame(
#'   cell_type = c(
#'     "immune", "vascular", "cortex GABAergic", "astrocyte", "ependymal", "OPC",
#'     "oligodendrocyte", "striatum MSN", "striatum non-MSN",
#'     "cortex glutamatergic"
#'   ), stringsAsFactors = FALSE
#' )
#' df$value <- runif(nrow(df))
#'
#' ggplot(df, aes(x = cell_type, y = value, color = cell_type)) +
#'   geom_point() +
#'   hbma_cell_type_color_scale()
#'
#' @export
hbma_cell_type_color_scale <- function(name = "Cell type", ...) {
  scale_color_manual(name = name, values = get_hbma_cell_type_palette(), ...)
}
