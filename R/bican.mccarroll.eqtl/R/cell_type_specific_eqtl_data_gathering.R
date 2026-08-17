#' Build a cluster-level eQTL analysis data frame
#'
#' Builds a long-form data frame for one eQTL cluster. Each row represents
#' one donor, cell-type/region dataset, and gene-variant pair. Both inverse-
#' normal-transformed and TPM expression values are gathered. Genotypes are
#' recovered across all in-group genotype files.
#'
#' @param cell_type_file Path to a TSV containing cell_type and region.
#' @param cluster_assignment_file Path to a TSV containing gene, variant_id,
#'   and cluster.
#' @param cluster Numeric cluster identifier.
#' @param in_cell_types Character vector defining the in-group cell types.
#' @param covariates Character vector of requested covariate names. Exact
#'   covariates are retained as numeric columns. Covariates represented by rows
#'   named `<covariate>_<level>` are automatically collapsed into a factor.
#' @param eqtl_root Root directory containing cell-type/region subdirectories.
#' @param regions Optional character vector of regions to retain.
#' @param outfile Optional path for the output TSV or TSV.GZ file.
#' @param covariate_list_outfile Optional path for the resolved covariate map.
#' @param report_outfile Optional path for the long-form processing report.
#'
#' @return A list with data, covariate_mapping, and report elements.
#' @export
build_cluster_dataframe <- function(
  cell_type_file,
  cluster_assignment_file,
  cluster,
  in_cell_types,
  covariates,
  eqtl_root,
  regions = NULL,
  outfile = NULL,
  covariate_list_outfile = NULL,
  report_outfile = NULL
) {
  .require_data_table()

  if (length(cluster) != 1L || is.na(cluster) || !is.numeric(cluster)) {
    stop("cluster must be one non-missing numeric value.", call. = FALSE)
  }
  if (!length(in_cell_types) || anyNA(in_cell_types)) {
    stop("in_cell_types must contain at least one non-missing value.", call. = FALSE)
  }
  if (anyNA(covariates)) {
    stop("covariates cannot contain missing values.", call. = FALSE)
  }

  datasets <- .read_dataset_table(cell_type_file, regions, in_cell_types)
  cluster_table <- .read_cluster_table(cluster_assignment_file, cluster)
  paths <- .construct_input_paths(datasets, in_cell_types, eqtl_root)

  .assert_files_exist(paths$expression_normalized, "normalized expression")
  .assert_files_exist(paths$expression_tpm, "TPM expression")
  .assert_files_exist(paths$genotype, "in-group genotype")
  .assert_files_exist(paths$covariate, "covariate")

  covariate_input <- .read_covariates(paths$covariate, covariates)
  covariate_mapping <- covariate_input$mapping

  normalized_headers <- lapply(paths$expression_normalized, .read_bed_header)
  normalized_donors <- lapply(normalized_headers, function(x) x$donors)
  tpm_headers <- lapply(paths$expression_tpm, .read_bed_header)
  tpm_donors <- lapply(tpm_headers, function(x) x$donors)
  common_donors <- Reduce(intersect, normalized_donors)

  if (!length(common_donors)) {
    stop("No donors are shared across all selected expression matrices.", call. = FALSE)
  }

  genotype_headers <- lapply(paths$genotype, .read_bed_header)
  genotype_donors <- unique(unlist(lapply(genotype_headers, function(x) x$donors)))
  missing_genotype_donors <- setdiff(common_donors, genotype_donors)
  if (length(missing_genotype_donors)) {
    warning(
      "Expression donors missing from all in-group genotype files: ",
      paste(missing_genotype_donors, collapse = ", "),
      call. = FALSE
    )
  }

  missing_covariate_donors <- setdiff(common_donors, covariate_input$donors)
  if (length(missing_covariate_donors)) {
    stop(
      "Expression donors missing from covariate file: ",
      paste(missing_covariate_donors, collapse = ", "),
      call. = FALSE
    )
  }

  covariate_long <- covariate_input$data[
    match(common_donors, covariate_input$data$donor), ,
    drop = FALSE
  ]

  genotype_result <- .read_genotypes(
    paths = paths$genotype,
    variants = cluster_table$variant_id,
    donors = common_donors,
    cell_type_regions = paths$genotype_cell_type_region
  )

  expression_results <- vector("list", nrow(datasets))
  per_dataset_reports <- vector("list", nrow(datasets))

  for (i in seq_len(nrow(datasets))) {
    result <- .read_expression_dataset(
      normalized_path = paths$expression_normalized[[i]],
      tpm_path = paths$expression_tpm[[i]],
      dataset = datasets[i, , drop = FALSE],
      genes = cluster_table$gene,
      donors = common_donors,
      normalized_donors_in_file = length(normalized_donors[[i]]),
      tpm_donors_in_file = length(tpm_donors[[i]])
    )
    expression_results[[i]] <- result$data
    per_dataset_reports[[i]] <- result$report
  }

  expression_long <- data.table::rbindlist(expression_results, use.names = TRUE)
  cluster_dt <- data.table::as.data.table(cluster_table)
  genotype_dt <- data.table::as.data.table(genotype_result$data)
  covariate_dt <- data.table::as.data.table(covariate_long)

  output <- merge(expression_long, cluster_dt, by = "gene", all.x = TRUE, sort = FALSE)
  output <- merge(
    output,
    genotype_dt,
    by = c("variant_id", "donor"),
    all.x = TRUE,
    sort = FALSE
  )
  output <- merge(output, covariate_dt, by = "donor", all.x = TRUE, sort = FALSE)

  output_columns <- c(
    "donor",
    "cell_type",
    "region",
    "cell_type_region",
    "gene",
    "variant_id",
    "cluster",
    "expression_normalized",
    "expression_tpm",
    "genotype",
    covariate_mapping$output_column,
    "in_group"
  )
  data.table::setcolorder(output, output_columns)
  data.table::setorderv(
    output,
    c("gene", "variant_id", "cell_type_region", "donor")
  )

  report <- .build_report(
    datasets = datasets,
    cluster = cluster,
    cluster_table = cluster_table,
    common_donors = common_donors,
    output = output,
    covariate_mapping = covariate_mapping,
    genotype_result = genotype_result,
    per_dataset_reports = per_dataset_reports
  )

  if (!is.null(outfile)) {
    data.table::fwrite(output, outfile, sep = "\t", na = "NA", compress = "auto")
  }

  if (!is.null(covariate_list_outfile)) {
    data.table::fwrite(
      covariate_mapping,
      covariate_list_outfile,
      sep = "\t",
      na = "NA",
      compress = "auto"
    )
  }

  if (!is.null(report_outfile)) {
    data.table::fwrite(
      report,
      report_outfile,
      sep = "\t",
      na = "NA",
      compress = "auto"
    )
  }

  list(
    data = as.data.frame(output),
    covariate_mapping = as.data.frame(covariate_mapping),
    report = as.data.frame(report)
  )
}

.require_data_table <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required.", call. = FALSE)
  }
}

.read_dataset_table <- function(cell_type_file, regions, in_cell_types) {
  x <- data.table::fread(cell_type_file, sep = "\t", data.table = FALSE)
  .require_columns(x, c("cell_type", "region"), "cell_type_file")

  if (!is.null(regions)) {
    x <- x[x$region %in% regions, , drop = FALSE]
  }
  if (!nrow(x)) {
    stop("No datasets remain after applying the region filter.", call. = FALSE)
  }

  missing_in <- setdiff(in_cell_types, unique(x$cell_type))
  if (length(missing_in)) {
    stop(
      "Requested in-group cell types absent after region filtering: ",
      paste(missing_in, collapse = ", "),
      call. = FALSE
    )
  }

  x$cell_type_region <- paste(x$cell_type, x$region, sep = "__")
  if (anyDuplicated(x$cell_type_region)) {
    duplicated_values <- unique(x$cell_type_region[duplicated(x$cell_type_region)])
    stop(
      "Duplicate cell-type/region datasets in cell_type_file: ",
      paste(duplicated_values, collapse = ", "),
      call. = FALSE
    )
  }
  x$in_group <- as.integer(x$cell_type %in% in_cell_types)
  x
}

.read_cluster_table <- function(cluster_assignment_file, cluster) {
  x <- data.table::fread(cluster_assignment_file, sep = "\t", data.table = FALSE)
  required <- c("gene", "variant_id", "cluster")
  .require_columns(x, required, "cluster_assignment_file")

  if (anyNA(x[, required, drop = FALSE])) {
    stop("cluster_assignment_file contains missing required values.", call. = FALSE)
  }

  selected <- x[x$cluster == cluster, required, drop = FALSE]
  if (!nrow(selected)) {
    stop("No gene-variant pairs found for cluster ", cluster, ".", call. = FALSE)
  }
  if (anyDuplicated(selected$gene)) {
    genes <- unique(selected$gene[duplicated(selected$gene)])
    stop(
      "Genes occur more than once in the selected cluster: ",
      paste(genes, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(selected[c("gene", "variant_id")])) {
    stop("Duplicate gene-variant rows occur in the selected cluster.", call. = FALSE)
  }

  selected
}

.construct_input_paths <- function(datasets, in_cell_types, eqtl_root) {
  dataset_dirs <- file.path(eqtl_root, datasets$cell_type_region)
  expression_normalized <- file.path(
    dataset_dirs,
    paste0(datasets$cell_type_region, ".gene_expression_normalized.bed.gz")
  )
  expression_tpm <- file.path(
    dataset_dirs,
    paste0(datasets$cell_type_region, ".gene_expression_tpm.bed.gz")
  )

  in_index <- which(datasets$cell_type %in% in_cell_types)
  in_names <- datasets$cell_type_region[in_index]
  in_dirs <- dataset_dirs[in_index]

  first_name <- in_names[[1L]]
  first_dir <- in_dirs[[1L]]

  list(
    expression_normalized = expression_normalized,
    expression_tpm = expression_tpm,
    genotype = file.path(in_dirs, paste0(in_names, ".genotype_matrix.bed.gz")),
    genotype_cell_type_region = in_names,
    covariate = file.path(first_dir, paste0(first_name, ".covariates_peer.txt"))
  )
}

.assert_files_exist <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Missing expected ", label, " file(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

.read_bed_header <- function(path) {
  header <- names(data.table::fread(path, nrows = 0L, sep = "\t"))
  if (length(header) < 5L) {
    stop("BED file has fewer than five columns: ", path, call. = FALSE)
  }
  expected <- c("#chr", "start", "end", "pid")
  if (!identical(header[seq_len(4L)], expected)) {
    stop(
      "Unexpected BED annotation columns in ", path,
      ". Expected: ", paste(expected, collapse = ", "),
      call. = FALSE
    )
  }
  list(columns = header, donors = header[-seq_len(4L)])
}

.read_covariates <- function(path, requested) {
  x <- data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  if (!ncol(x) || names(x)[1L] != "ID") {
    stop("Covariate file must have ID as its first column: ", path, call. = FALSE)
  }
  if (anyDuplicated(x$ID)) {
    duplicated_ids <- unique(x$ID[duplicated(x$ID)])
    stop(
      "Duplicate covariate IDs in covariate file: ",
      paste(duplicated_ids, collapse = ", "),
      call. = FALSE
    )
  }

  resolved <- .resolve_covariates(requested, x$ID)
  donors <- names(x)[-1L]
  covariate_data <- data.frame(
    donor = donors,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  mapping_rows <- vector("list", length(resolved))

  for (i in seq_along(resolved)) {
    item <- resolved[[i]]
    selected <- x[match(item$source_columns, x$ID), , drop = FALSE]
    values <- .covariate_numeric_matrix(selected, item$requested_covariate)

    if (item$encoding == "numeric") {
      value <- as.numeric(values[1L, ])
    } else {
      value <- .collapse_one_hot_covariate(
        values = values,
        source_columns = item$source_columns,
        covariate_name = item$requested_covariate,
        donors = donors
      )
    }

    covariate_data[[item$output_column]] <- value
    mapping_rows[[i]] <- data.frame(
      requested_covariate = item$requested_covariate,
      output_column = item$output_column,
      encoding = item$encoding,
      source_columns = paste(item$source_columns, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  if (anyNA(covariate_data[-1L])) {
    missing_by_column <- names(covariate_data)[-1L][
      vapply(covariate_data[-1L], anyNA, logical(1))
    ]
    stop(
      "Missing covariate values detected in: ",
      paste(missing_by_column, collapse = ", "),
      call. = FALSE
    )
  }

  mapping <- if (length(mapping_rows)) {
    do.call(rbind, mapping_rows)
  } else {
    data.frame(
      requested_covariate = character(),
      output_column = character(),
      encoding = character(),
      source_columns = character(),
      stringsAsFactors = FALSE
    )
  }

  list(data = covariate_data, mapping = mapping, donors = donors)
}

.resolve_covariates <- function(requested, available) {
  if (!length(requested)) {
    return(list())
  }
  if (anyDuplicated(requested)) {
    duplicated_requests <- unique(requested[duplicated(requested)])
    stop(
      "Covariates were requested more than once: ",
      paste(duplicated_requests, collapse = ", "),
      call. = FALSE
    )
  }

  lapply(requested, function(name) {
    exact <- which(available == name)
    prefix_pattern <- paste0("^", .escape_regex(name), "_")
    prefixed <- grep(prefix_pattern, available, perl = TRUE)

    if (length(exact) && length(prefixed)) {
      stop(
        "Ambiguous covariate request '", name,
        "': both an exact row and one-hot encoded rows exist.",
        call. = FALSE
      )
    }
    if (!length(exact) && !length(prefixed)) {
      stop("Requested covariate not found: ", name, call. = FALSE)
    }

    if (length(exact)) {
      list(
        requested_covariate = name,
        output_column = name,
        encoding = "numeric",
        source_columns = available[exact]
      )
    } else {
      list(
        requested_covariate = name,
        output_column = name,
        encoding = "factor",
        source_columns = available[prefixed]
      )
    }
  })
}

.covariate_numeric_matrix <- function(selected, covariate_name) {
  converted <- lapply(selected[-1L], function(column) {
    value <- suppressWarnings(as.numeric(column))
    invalid <- is.na(value) & !is.na(column)
    if (any(invalid)) {
      stop(
        "Covariate '", covariate_name,
        "' contains non-numeric values.",
        call. = FALSE
      )
    }
    value
  })

  values <- as.matrix(as.data.frame(converted, check.names = FALSE))
  rownames(values) <- selected$ID
  colnames(values) <- names(selected)[-1L]
  storage.mode(values) <- "numeric"
  values
}

.collapse_one_hot_covariate <- function(
  values,
  source_columns,
  covariate_name,
  donors
) {
  if (anyNA(values)) {
    bad_donors <- donors[colSums(is.na(values)) > 0L]
    stop(
      "Missing values in one-hot covariate '", covariate_name,
      "' for donors: ", paste(bad_donors, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(!values %in% c(0, 1))) {
    stop(
      "One-hot covariate '", covariate_name,
      "' contains values other than 0 and 1.",
      call. = FALSE
    )
  }

  n_selected <- colSums(values == 1)
  invalid <- n_selected != 1L
  if (any(invalid)) {
    details <- paste0(donors[invalid], "=", n_selected[invalid])
    stop(
      "One-hot covariate '", covariate_name,
      "' must have exactly one selected level per donor. Invalid donors: ",
      paste(details, collapse = ", "),
      call. = FALSE
    )
  }

  prefix_pattern <- paste0("^", .escape_regex(covariate_name), "_")
  levels <- sub(prefix_pattern, "", source_columns, perl = TRUE)
  if (any(!nzchar(levels)) || anyDuplicated(levels)) {
    stop(
      "Could not derive unique factor levels for one-hot covariate '",
      covariate_name, "'.",
      call. = FALSE
    )
  }

  selected_row <- max.col(t(values), ties.method = "first")
  factor(levels[selected_row], levels = levels)
}

.escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\?.])", "\\\\\\1", x)
}


.read_genotypes <- function(paths, variants, donors, cell_type_regions) {
  # Make R CMD CHECK happy
  pid <- .N <- . <- variant_id <- donor <- N <- NULL

  variants <- unique(as.character(variants))
  donors <- unique(as.character(donors))

  if (length(paths) != length(cell_type_regions)) {
    stop(
      "`paths` and `cell_type_regions` must have the same length.",
      call. = FALSE
    )
  }

  recovered <- matrix(
    NA_real_,
    nrow = length(variants),
    ncol = length(donors),
    dimnames = list(variants, donors)
  )

  sources <- matrix(
    NA_character_,
    nrow = length(variants),
    ncol = length(donors),
    dimnames = list(variants, donors)
  )

  variants_seen <- character()
  duplicate_variants <- character()
  duplicate_variant_rows <- 0L
  per_file <- vector("list", length(paths))

  for (i in seq_along(paths)) {
    header <- .read_bed_header(paths[[i]])
    donors_available <- intersect(donors, header$donors)
    select_columns <- c("pid", donors_available)

    x <- data.table::fread(
      paths[[i]],
      sep = "\t",
      select = select_columns,
      data.table = TRUE,
      showProgress = FALSE
    )

    x <- x[pid %in% variants]

    counts <- table(x$pid)
    duplicates_here <- names(counts[counts > 1L])

    if (length(duplicates_here)) {
      duplicate_variants <- union(
        duplicate_variants,
        duplicates_here
      )

      duplicate_variant_rows <- duplicate_variant_rows +
        sum(counts[counts > 1L] - 1L)

      x <- x[!duplicated(pid)]
    }

    variants_seen <- union(variants_seen, x$pid)
    present_variants <- intersect(variants, x$pid)

    if (length(present_variants) && length(donors_available)) {
      row_index <- match(present_variants, x$pid)

      values <- as.matrix(
        x[row_index, donors_available, with = F]
      )

      storage.mode(values) <- "numeric"
      rownames(values) <- present_variants
      colnames(values) <- donors_available

      target <- recovered[
        present_variants,
        donors_available,
        drop = FALSE
      ]

      conflict <- (
        !is.na(target) &
          !is.na(values) &
          target != values
      )

      if (any(conflict)) {
        positions <- which(conflict, arr.ind = TRUE)

        details <- apply(
          positions,
          1L,
          function(pos) {
            variant <- present_variants[[pos[[1L]]]]
            donor <- donors_available[[pos[[2L]]]]

            paste0(
              variant,
              "/",
              donor,
              " (",
              target[pos[[1L]], pos[[2L]]],
              " vs ",
              values[pos[[1L]], pos[[2L]]],
              ")"
            )
          }
        )

        stop(
          paste0(
            "Conflicting recovered genotypes across ",
            "in-group files: "
          ),
          paste(details, collapse = ", "),
          call. = FALSE
        )
      }

      fill <- is.na(target) & !is.na(values)
      target[fill] <- values[fill]

      recovered[
        present_variants,
        donors_available
      ] <- target

      source_target <- sources[
        present_variants,
        donors_available,
        drop = FALSE
      ]

      source_target[fill] <- cell_type_regions[[i]]

      sources[
        present_variants,
        donors_available
      ] <- source_target
    }

    missing_here <- setdiff(variants, x$pid)

    per_file[[i]] <- data.table::data.table(
      scope = "genotype_file",
      cell_type_region = cell_type_regions[[i]],
      metric = c(
        "donors_in_genotype",
        "requested_variants",
        "variants_found",
        "variants_missing",
        "duplicate_variant_rows"
      ),
      value = as.character(
        c(
          length(header$donors),
          length(variants),
          length(unique(x$pid)),
          length(missing_here),
          if (length(counts)) {
            sum(counts[counts > 1L] - 1L)
          } else {
            0L
          }
        )
      ),
      details = c(
        NA_character_,
        NA_character_,
        NA_character_,
        paste(missing_here, collapse = ","),
        paste(duplicates_here, collapse = ",")
      )
    )
  }

  missing_variants <- setdiff(variants, variants_seen)

  if (length(missing_variants)) {
    warning(
      paste0(
        "Requested variants missing from all in-group ",
        "genotype files: "
      ),
      paste(missing_variants, collapse = ", "),
      call. = FALSE
    )
  }

  genotype_dt <- data.table::as.data.table(
    recovered,
    keep.rownames = "variant_id"
  )

  genotype_long <- data.table::melt(
    genotype_dt,
    id.vars = "variant_id",
    measure.vars = donors,
    variable.name = "donor",
    value.name = "genotype",
    variable.factor = FALSE
  )

  duplicate_keys <- genotype_long[
    ,
    .N,
    by = .(variant_id, donor)
  ][N > 1L]

  if (nrow(duplicate_keys)) {
    stop(
      paste0(
        "Internal error: recovered genotype data contain ",
        "duplicate variant-donor pairs."
      ),
      call. = FALSE
    )
  }

  list(
    data = genotype_long,
    duplicate_variants = duplicate_variants,
    duplicate_variant_rows = duplicate_variant_rows,
    variants_found = length(variants_seen),
    recovered_values = sum(!is.na(recovered)),
    missing_values = sum(is.na(recovered)),
    per_file_reports = per_file
  )
}

.read_expression_dataset <- function(
  normalized_path,
  tpm_path,
  dataset,
  genes,
  donors,
  normalized_donors_in_file,
  tpm_donors_in_file
) {
  # Make R CMD CHECK happy
  cell_type <- region <- cell_type_region <- in_group <- NULL

  normalized <- .read_expression_matrix(normalized_path, genes, donors)
  tpm <- .read_expression_matrix(tpm_path, genes, donors)

  normalized_dt <- data.table::as.data.table(
    normalized$matrix,
    keep.rownames = "gene"
  )
  normalized_long <- data.table::melt(
    normalized_dt,
    id.vars = "gene",
    measure.vars = donors,
    variable.name = "donor",
    value.name = "expression_normalized",
    variable.factor = FALSE
  )

  tpm_dt <- data.table::as.data.table(tpm$matrix, keep.rownames = "gene")
  tpm_long <- data.table::melt(
    tpm_dt,
    id.vars = "gene",
    measure.vars = donors,
    variable.name = "donor",
    value.name = "expression_tpm",
    variable.factor = FALSE
  )

  expression_long <- merge(
    normalized_long,
    tpm_long,
    by = c("gene", "donor"),
    all = TRUE,
    sort = FALSE
  )
  expression_long[, cell_type := dataset$cell_type[[1L]]]
  expression_long[, region := dataset$region[[1L]]]
  expression_long[, cell_type_region := dataset$cell_type_region[[1L]]]
  expression_long[, in_group := dataset$in_group[[1L]]]

  report <- data.table::data.table(
    scope = "cell_type_region",
    cell_type_region = dataset$cell_type_region[[1L]],
    metric = c(
      "donors_in_normalized_expression",
      "donors_in_tpm_expression",
      "intersected_donors",
      "requested_genes",
      "normalized_genes_found",
      "normalized_genes_missing",
      "tpm_genes_found",
      "tpm_genes_missing",
      "rows_emitted",
      "missing_normalized_expression_values",
      "missing_tpm_expression_values",
      "in_group"
    ),
    value = as.character(c(
      normalized_donors_in_file,
      tpm_donors_in_file,
      length(donors),
      length(genes),
      normalized$genes_found,
      length(genes) - normalized$genes_found,
      tpm$genes_found,
      length(genes) - tpm$genes_found,
      nrow(expression_long),
      sum(is.na(expression_long$expression_normalized)),
      sum(is.na(expression_long$expression_tpm)),
      dataset$in_group[[1L]]
    )),
    details = NA_character_
  )

  list(data = expression_long, report = report)
}

.read_expression_matrix <- function(path, genes, donors) {
  # Make R CMD CHECK happy
  pid <- NULL

  header <- .read_bed_header(path)
  donors_available <- intersect(donors, header$donors)

  x <- data.table::fread(
    path,
    sep = "\t",
    select = c("pid", donors_available),
    data.table = TRUE,
    showProgress = FALSE
  )
  x <- x[pid %in% genes]

  expression_matrix <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(donors),
    dimnames = list(genes, donors)
  )
  row_index <- match(genes, x$pid)
  present <- !is.na(row_index)
  if (any(present) && length(donors_available)) {
    values <- as.matrix(x[row_index[present], donors_available, with = F])
    storage.mode(values) <- "numeric"
    expression_matrix[present, donors_available] <- values
  }

  list(
    matrix = expression_matrix,
    genes_found = sum(present)
  )
}

.build_report <- function(
  datasets,
  cluster,
  cluster_table,
  common_donors,
  output,
  covariate_mapping,
  genotype_result,
  per_dataset_reports
) {
  duplicate_details <- if (length(genotype_result$duplicate_variants)) {
    paste(genotype_result$duplicate_variants, collapse = ",")
  } else {
    NA_character_
  }

  overall <- data.table::data.table(
    scope = "overall",
    cell_type_region = NA_character_,
    metric = c(
      "selected_cluster",
      "selected_datasets",
      "in_group_datasets",
      "gene_variant_pairs",
      "variants_found_across_in_group",
      "variants_missing_across_in_group",
      "duplicate_variant_rows",
      "intersected_donors",
      "rows_emitted",
      "missing_normalized_expression_values",
      "missing_tpm_expression_values",
      "recovered_genotype_values",
      "missing_genotype_values",
      "resolved_covariate_columns"
    ),
    value = as.character(c(
      cluster,
      nrow(datasets),
      sum(datasets$in_group),
      nrow(cluster_table),
      genotype_result$variants_found,
      0L,
      genotype_result$duplicate_variant_rows,
      length(common_donors),
      nrow(output),
      sum(is.na(output$expression_normalized)),
      sum(is.na(output$expression_tpm)),
      genotype_result$recovered_values,
      sum(is.na(output$genotype)),
      nrow(covariate_mapping)
    )),
    details = c(
      rep(NA_character_, 6L),
      duplicate_details,
      rep(NA_character_, 7L)
    )
  )

  data.table::rbindlist(
    c(
      list(overall),
      genotype_result$per_file_reports,
      per_dataset_reports
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

.require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}
