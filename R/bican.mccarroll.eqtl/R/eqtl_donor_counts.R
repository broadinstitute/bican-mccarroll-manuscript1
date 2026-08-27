# bican.mccarroll.eqtl::count_eqtl_donors(
#     "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/eqtls/results/LEVEL_6"
# )


#' Count donors per cell type / region from eQTL covariates files
#'
#' Scans every \code{<cell_type>__<region>} subdirectory of an eQTL results
#' directory and counts the number of donors used for that cell type/region
#' test, read from the first line of its tensorQTL covariates file
#' (\code{<cell_type>__<region>.covariates.txt}), which lists \code{id}
#' followed by one column per donor.
#'
#' Subdirectory names are split on the literal \code{"__"} separator (the
#' same convention used throughout this package, e.g.
#' \code{compare_eqtl_runs_ctr()}). Cell type names may themselves contain
#' region-like tokens (e.g. \code{GABA_CGE_DFC__DFC}), so subdirectories are
#' required to split into exactly two parts on \code{"__"}; anything else is
#' a naming-convention violation and stops execution.
#'
#' If a subdirectory's covariates file is missing or malformed, that row is
#' still included with \code{n_donor = NA} and a warning naming the
#' subdirectory, rather than being silently dropped.
#'
#' @param eqtl_dir Character scalar. Directory containing
#'   \code{<cell_type>__<region>} eQTL result subdirectories (e.g. one
#'   \code{LEVEL_*} directory).
#' @param output_path Character scalar or \code{NULL}. If non-NULL, the
#'   result table is written to this path as a tab-delimited file.
#'
#' @return A \code{data.frame} with columns \code{cell_type}, \code{region},
#'   and \code{n_donor}, one row per subdirectory, ordered by cell type then
#'   region.
#'
#' @importFrom logger log_info
#' @importFrom data.table fwrite
#' @export
count_eqtl_donors <- function(eqtl_dir, output_path = NULL) {
  file_separator <- "__"

  data_sets <- list.dirs(eqtl_dir, full.names = FALSE, recursive = FALSE)
  if (length(data_sets) == 0) {
    stop("No data sets found in eQTL directory: ", eqtl_dir, call. = FALSE)
  }

  z <- strsplit(data_sets, file_separator, fixed = TRUE)
  bad <- data_sets[lengths(z) != 2L]
  if (length(bad) > 0) {
    stop(
      "Expected subdirectory names to split into exactly 2 parts on '",
      file_separator, "': ", paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  cell_types <- vapply(z, `[`, character(1L), 1L)
  regions <- vapply(z, `[`, character(1L), 2L)

  n_donor <- vapply(
    data_sets,
    function(d) {
      covariate_file <- file.path(eqtl_dir, d, paste0(d, ".covariates.txt"))
      .count_covariates_donors(covariate_file)
    },
    integer(1L)
  )

  df <- data.frame(
    cell_type = cell_types,
    region = regions,
    n_donor = n_donor,
    stringsAsFactors = FALSE
  )

  df <- df[order(df$cell_type, df$region, method = "radix"), ]
  rownames(df) <- NULL

  logger::log_info(
    "Counted donors for {nrow(df)} cell type/region eQTL data sets in {eqtl_dir}"
  )

  if (!is.null(output_path)) {
    data.table::fwrite(df, output_path, sep = "\t")
    logger::log_info("Written to: {output_path}")
  }

  df
}


.count_covariates_donors <- function(path) {
  if (!file.exists(path)) {
    warning("Covariates file not found: ", path, call. = FALSE)
    return(NA_integer_)
  }

  header <- readLines(path, n = 1L)
  fields <- strsplit(header, "\t", fixed = TRUE)[[1L]]

  if (length(fields) == 0 || !grepl("^id$", fields[1L], ignore.case = TRUE)) {
    warning(
      "Covariates file does not start with an 'id' header row: ", path,
      call. = FALSE
    )
    return(NA_integer_)
  }

  as.integer(length(fields) - 1L)
}
