#' bican.mccarroll.figures
#'
#' Figure generation utilities for BICAN McCarroll analyses.
#'
#' This package provides functions that orchestrate downstream analysis
#' workflows and render manuscript-quality figures. Each figure function
#' is capable of reproducing the necessary analysis steps from processed
#' input data. For long-running computations, intermediate results may
#' be cached to avoid recomputation.
#'
#' # Path Configuration
#'
#' All figure functions follow a consistent rule for resolving file paths:
#'
#' 1. If a path argument is supplied explicitly, it is used as-is.
#' 2. If a path argument is NULL, the package attempts to resolve the path
#'    relative to a configured data root directory.
#' 3. If neither a path argument nor a configured option is available,
#'    the function stops with an informative error.
#'
#' The data root directory should point to the top-level directory of the
#' unpacked analysis data tarball. Default relative paths are interpreted
#' beneath this directory.
#'
#' # Required Options
#'
#' Users must configure the following options before running figure
#' functions, unless all required paths are supplied explicitly:
#'
#' - `bican.mccarroll.figures.data_root_dir`
#' - `bican.mccarroll.figures.out_dir`
#' - `bican.mccarroll.figures.cache_dir`
#'
#' These may be set using `options()`:
#'
#' ```
#' options(
#'   bican.mccarroll.figures.data_root_dir = "/path/to/data_root",
#'   bican.mccarroll.figures.out_dir = "/path/to/output",
#'   bican.mccarroll.figures.cache_dir = "/path/to/cache"
#' )
#' ```
#'
#' # Absolute vs Relative Paths
#'
#' For input arguments:
#'
#' - Absolute paths are used directly.
#' - Relative paths are interpreted relative to
#'   `bican.mccarroll.figures.data_root_dir`.
#'
#' Output and cache directories are never interpreted relative to the
#' data root directory. They must be either supplied explicitly or
#' configured via options.
#'
#' No default paths are set when the package is loaded. Users must
#' configure paths explicitly in their session or script.
#'
#' @keywords internal
"_PACKAGE"
