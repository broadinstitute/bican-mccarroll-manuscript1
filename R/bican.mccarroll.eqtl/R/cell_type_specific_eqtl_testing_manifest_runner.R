#' Run cell-type-specific eQTL analyses from a manifest
#'
#' Runs the eQTL data-gathering and linear-model testing steps for each row
#' of a tab-separated manifest. Each manifest row defines a cluster, one or
#' more in-group cell types, and one or more regions. Shared input files and
#' model covariates are supplied once as function arguments.
#'
#' The manifest must contain the following columns:
#'
#' \describe{
#'   \item{\code{cluster}}{Integer cluster identifier.}
#'   \item{\code{in_cell_types}}{
#'     One or more in-group cell types. Multiple values should be
#'     comma-separated.
#'   }
#'   \item{\code{regions}}{
#'     One or more brain regions. Multiple values should be comma-separated.
#'   }
#' }
#'
#' Output filenames are generated from the cluster, in-group cell types, and
#' regions. For example, cluster 10, astrocyte, and CaH produce the prefix
#' \code{cluster10_astrocyte__CaH}.
#'
#' @param manifest_file Path to a tab-separated manifest containing the
#'   columns \code{cluster}, \code{in_cell_types}, and \code{regions}.
#' @param output_dir Directory in which analysis outputs will be written.
#'   The directory is created recursively if it does not already exist.
#' @param cell_type_file Path to the region and cell-type metadata file passed
#'   to \code{\link{build_cluster_dataframe}}.
#' @param cluster_assignment_file Path to the cluster-assignment file passed
#'   to \code{\link{build_cluster_dataframe}}.
#' @param covariates Character vector of donor-level covariates to include in
#'   the gathered analysis dataset.
#' @param eqtl_root Root directory containing the cell-type-specific eQTL
#'   result files.
#' @param force_factor_covariates Character vector of covariates that should
#'   be encoded as factors during linear-model testing.
#'
#' @return Invisibly returns a list with one element per manifest row. Each
#'   element contains the paths of the files generated for that analysis.
#'
#' @details
#' For each manifest row, the function first calls
#' \code{\link{build_cluster_dataframe}} and then calls
#' \code{\link{run_cell_type_specific_eqtl_tests}} using the generated data
#' and covariate-mapping files.
#'
#' Multiple values in the \code{in_cell_types} and \code{regions} manifest
#' columns must be separated by commas without quoting the individual values.
#'
#' @examples
#' \dontrun{
#' run_cell_type_specific_eqtl_tests_from_manifest(
#'     manifest_file = "private_eqtl_manifest.tsv",
#'     output_dir = "private_eqtl_analysis",
#'     cell_type_file = "region_cell_type.tsv",
#'     cluster_assignment_file = "cluster_assignments.tsv",
#'     covariates = c(
#'         "biobank",
#'         "age",
#'         "imputed_sex",
#'         paste0("PC", 1:5),
#'         "pmi_hr",
#'         paste0("InferredCov", 1:10)
#'     ),
#'     eqtl_root = "results/LEVEL_6",
#'     force_factor_covariates = c(
#'         "biobank",
#'         "imputed_sex"
#'     )
#' )
#' }
#'
#' @importFrom data.table fread
#' @export
run_cell_type_specific_eqtl_tests_from_manifest <- function(
        manifest_file,
        output_dir,
        cell_type_file,
        cluster_assignment_file,
        covariates,
        eqtl_root,
        force_factor_covariates
) {
    manifest <- data.table::fread(manifest_file)

    dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    output_paths <- vector(
        "list",
        nrow(manifest)
    )

    for (i in seq_len(nrow(manifest))) {
        in_cell_types <- .parse_manifest_vector(
            manifest$in_cell_types[[i]]
        )

        regions <- .parse_manifest_vector(
            manifest$regions[[i]]
        )

        output_paths[[i]] <- .run_private_eqtl_analysis(
            cluster = manifest$cluster[[i]],
            in_cell_types = in_cell_types,
            regions = regions,
            output_dir = output_dir,
            cell_type_file = cell_type_file,
            cluster_assignment_file = cluster_assignment_file,
            covariates = covariates,
            eqtl_root = eqtl_root,
            force_factor_covariates = force_factor_covariates
        )
    }

    invisible(output_paths)
}
run_private_eqtl_manifest <- function(
        manifest_file,
        output_dir,
        cell_type_file,
        cluster_assignment_file,
        covariates,
        eqtl_root,
        force_factor_covariates
) {
    manifest <- data.table::fread(manifest_file)

    dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    output_paths <- vector(
        "list",
        nrow(manifest)
    )

    for (i in seq_len(nrow(manifest))) {
        in_cell_types <- .parse_manifest_vector(
            manifest$in_cell_types[[i]]
        )

        regions <- .parse_manifest_vector(
            manifest$regions[[i]]
        )

        output_paths[[i]] <- .run_private_eqtl_analysis(
            cluster = manifest$cluster[[i]],
            in_cell_types = in_cell_types,
            regions = regions,
            output_dir = output_dir,
            cell_type_file = cell_type_file,
            cluster_assignment_file = cluster_assignment_file,
            covariates = covariates,
            eqtl_root = eqtl_root,
            force_factor_covariates = force_factor_covariates
        )
    }

    invisible(output_paths)
}

.parse_manifest_vector <- function(x) {
    strsplit(x, ",", fixed = TRUE)[[1]]
}


.make_analysis_id <- function(cluster, in_cell_types, regions) {
    paste0(
        "cluster",
        cluster,
        "_",
        paste(in_cell_types, collapse = "-"),
        "__",
        paste(regions, collapse = "-")
    )
}


.make_output_paths <- function(output_dir, analysis_id) {
    prefix <- file.path(output_dir, analysis_id)

    list(
        gathered_file = paste0(prefix, ".data.txt"),
        covariate_mapping_file = paste0(prefix, ".covariate_list.txt"),
        report_outfile = paste0(prefix, ".report.txt"),
        results_outfile = paste0(prefix, ".linear_model_test.txt"),
        outPDF = paste0(prefix, ".linear_model_test.pdf")
    )
}


.run_private_eqtl_analysis <- function(
        cluster,
        in_cell_types,
        regions,
        output_dir,
        cell_type_file,
        cluster_assignment_file,
        covariates,
        eqtl_root,
        force_factor_covariates
) {
    analysis_id <- .make_analysis_id(
        cluster = cluster,
        in_cell_types = in_cell_types,
        regions = regions
    )

    paths <- .make_output_paths(
        output_dir = output_dir,
        analysis_id = analysis_id
    )

    message("Running ", analysis_id)

    bican.mccarroll.eqtl::build_cluster_dataframe(
        cell_type_file = cell_type_file,
        cluster_assignment_file = cluster_assignment_file,
        cluster = cluster,
        in_cell_types = in_cell_types,
        covariates = covariates,
        eqtl_root = eqtl_root,
        regions = regions,
        outfile = paths$gathered_file,
        covariate_list_outfile = paths$covariate_mapping_file,
        report_outfile = paths$report_outfile
    )

    bican.mccarroll.eqtl::run_cell_type_specific_eqtl_tests(
        gathered_file = paths$gathered_file,
        covariate_mapping_file = paths$covariate_mapping_file,
        force_factor_covariates = force_factor_covariates,
        results_outfile = paths$results_outfile,
        outPDF = paths$outPDF
    )

    invisible(paths)
}


