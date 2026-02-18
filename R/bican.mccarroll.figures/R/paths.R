## R/paths.R

#' Get the configured data root directory
#'
#' Returns the currently configured data root directory.
#' This directory should point to the root of the unpacked
#' data tarball containing all required relative input paths.
#'
#' The value is taken from the option
#' `bican.mccarroll.figures.data_root_dir`.
#'
#' @return Character scalar or NULL if not configured.
#'
#' @examples
#' options(bican.mccarroll.figures.data_root_dir = "/path/to/data_root")
#' get_data_root_dir()
#'
#' @export
get_data_root_dir <- function() {
    getOption("bican.mccarroll.figures.data_root_dir", default = NULL)
}


#' Get the configured figure output directory
#'
#' Returns the currently configured output directory where
#' figures will be written. This directory is independent
#' of the data root directory and is typically outside
#' the unpacked tarball.
#'
#' The value is taken from the option
#' `bican.mccarroll.figures.out_dir`.
#'
#' @return Character scalar or NULL if not configured.
#'
#' @examples
#' options(bican.mccarroll.figures.out_dir = "/path/to/output")
#' get_out_dir()
#'
#' @export
get_out_dir <- function() {
    getOption("bican.mccarroll.figures.out_dir", default = NULL)
}

#' Get the configured cache directory
#'
#' Returns the currently configured cache directory used
#' for intermediate results. This directory is independent
#' of the data root directory.
#'
#' The value is taken from the option
#' `bican.mccarroll.figures.cache_dir`.
#'
#' @return Character scalar or NULL if not configured.
#'
#' @examples
#' options(bican.mccarroll.figures.cache_dir = "/path/to/cache")
#' get_cache_dir()
#'
#' @export
get_cache_dir <- function() {
    getOption("bican.mccarroll.figures.cache_dir", default = NULL)
}

#' Get all configured figure-related paths
#'
#' Returns a named list containing the currently configured
#' data root directory, output directory, and cache directory.
#' Values may be NULL if not yet configured.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{data_root_dir}{Character scalar or NULL}
#'     \item{out_dir}{Character scalar or NULL}
#'     \item{cache_dir}{Character scalar or NULL}
#'   }
#'
#' @examples
#' options(
#'   bican.mccarroll.figures.data_root_dir = "/path/to/data_root",
#'   bican.mccarroll.figures.out_dir = "/path/to/output",
#'   bican.mccarroll.figures.cache_dir = "/path/to/cache"
#' )
#' get_figure_paths()
#'
#' @export
get_figure_paths <- function() {
    list(
        data_root_dir = get_data_root_dir(),
        out_dir = get_out_dir(),
        cache_dir = get_cache_dir()
    )
}

## -------------------------------------------------------------------------
## Internal helpers (not exported)
## -------------------------------------------------------------------------

.resolve_required_path <- function(arg_value, option_name, what) {
    if (!is.null(arg_value)) return(arg_value)

    opt <- getOption(option_name, default = NULL)
    if (is.null(opt) || !nzchar(opt)) {
        stop(
            what, " is not set. Either pass it explicitly or set options('",
            option_name, "' = '/path/to/...').",
            call. = FALSE
        )
    }
    opt
}

.resolve_data_root_dir <- function(data_root_dir = NULL) {
    .resolve_required_path(
        arg_value = data_root_dir,
        option_name = "bican.mccarroll.figures.data_root_dir",
        what = "Data root directory"
    )
}

.resolve_out_dir <- function(out_dir = NULL) {
    .resolve_required_path(
        arg_value = out_dir,
        option_name = "bican.mccarroll.figures.out_dir",
        what = "Output directory"
    )
}

.resolve_cache_dir <- function(cache_dir = NULL) {
    .resolve_required_path(
        arg_value = cache_dir,
        option_name = "bican.mccarroll.figures.cache_dir",
        what = "Cache directory"
    )
}

.ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(NULL)
}

.is_absolute_path <- function(p) {
    ## POSIX absolute path or Windows drive / UNC.
    if (!is.character(p) || length(p) != 1L || is.na(p) || p == "") return(FALSE)
    if (substr(p, 1L, 1L) == "/") return(TRUE)
    if (grepl("^[A-Za-z]:[\\\\/]", p)) return(TRUE)
    if (grepl("^\\\\\\\\", p)) return(TRUE)
    FALSE
}

.resolve_under_root <- function(root, p) {
    ## If p is absolute, return as-is; if relative, interpret under root.
    if (.is_absolute_path(p)) return(p)
    file.path(root, p)
}
