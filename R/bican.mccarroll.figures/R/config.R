#' Configure the data root directory
#'
#' Set and query the root directory of the released data tree used by
#' figure-generation functions in this package. When configured, plotting
#' functions can derive their default input and output locations from a single
#' data root directory, avoiding repeated hard-coded paths in function calls.
#'
#' The root directory is stored in the session option
#' `bican.mccarroll.figures.data_root_dir`. You may set this option directly
#' using `options()`, but `set_data_root_dir()` performs basic validation and
#' path normalization.
#'
#' @param root Path to the data root directory.
#' @param must_work Logical; if `TRUE`, require that `root` exists.
#'
#' @return `set_data_root_dir()` returns the normalized root path invisibly.
#'   `get_data_root_dir()` returns the configured root path, or `NULL` if unset.
#'
#' @examples
#' \dontrun{
#' set_data_root_dir("/path/to/CAP_freeze_3_analysis")
#' get_data_root_dir()
#' }
#'
#' @export
set_data_root_dir <- function(root, must_work = TRUE) {
    if (!is.character(root) || length(root) != 1L || is.na(root) || root == "") {
        stop("`root` must be a non-empty character scalar.", call. = FALSE)
    }

    root2 <- normalizePath(root, winslash = "/", mustWork = must_work)

    if (must_work && !dir.exists(root2)) {
        stop("Directory does not exist: ", root2, call. = FALSE)
    }

    options(bican.mccarroll.figures.data_root_dir = root2)

    invisible(root2)
}

#' @rdname set_data_root_dir
#' @export
get_data_root_dir <- function() {
    getOption("bican.mccarroll.figures.data_root_dir", default = NULL)
}

get_data_paths <- function() {
    root <- get_data_root_dir()
    if (is.null(root)) {
        stop(
            "Data root directory is not configured. Call set_data_root_dir('/path/to/data_root').",
            call. = FALSE
        )
    }

    fig_dir <- file.path(root, "figure_repository")
    cache_dir <- file.path(fig_dir, "data_cache")

    list(
        root = root,
        fig_dir = fig_dir,
        cache_dir = cache_dir
    )
}
