# Load required libraries
# library(edgeR)
# library(data.table)


metacell_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
manifest_file="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metadata/DE_manifest.txt"
cell_metadata_file="/broad/bican_um1_mccarroll/RNAseq/analysis/cellarium_upload/CAP_freeze_2/CAP_cell_metadata.annotated.txt.gz"
cellTypeProportionsPCAFile=NULL
has_village=T;

# this doesn't include donor, region or village, as those are
# automatically added by the manifest processing.
metadata_columns=c("biobank", "cohort", "age", "imputed_sex", "pmi_hr", "PC1", "PC2", "PC3", "PC4", "PC5", "toxicology_group", "single_cell_assay", "pct_intronic", "frac_contamination", "hbcac_status")
outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_2_analysis/differential_expression/metacells"
outName="donor_rxn_DGEList"

# z=bican.mccarroll.differentialexpression::build_merged_dge(manifest_file, metacell_dir, cell_metadata_file, metadata_columns,has_village, outDir, outName, validate_round_trip=T)


# NOTE: For limma/voom runs, there is some metacells groups that are essentially duplicates.
# We'll want to filter those out for limma/voom runs.
# column retain_for_differential_expression should be true for all limma voom or partition runs.


#' Build a Merged DGEList from Metacell Files and Metadata
#'
#' This function processes a manifest and corresponding metacell count matrices to
#' build a merged \code{DGEList} object suitable for analysis with edgeR. The input
#' matrices are expected to be pseudobulked gene expression data with genes in rows
#' and donors in columns. Column names are renamed to include donor, cell type,
#' and region. Additional donor-level metadata is merged in from an external file.
#'
#' @param manifest_file Path to the manifest file, with required columns: \code{output_name}, \code{cell_type_name}, \code{region_name}, and \code{cell_type_filter_str}.
#' @param metacell_dir Path to the directory containing \code{*.metacells.txt.gz} files.
#' @param cell_metadata_file Path to the TSV file containing donor-level metadata.
#' @param metadata_columns Character vector of metadata column names to include from \code{cell_metadata_file}.  The
#' program will automatically capture donor, cell type, region, and village (if applicable).
#' @param has_village Logical; if \code{TRUE}, the donor names in the metacell files are expected to be in the format \code{donor_village}.
#' @param outDir Optional; directory to save the merged \code{DGEList} object as two gzipped TSV files.
#' @param outName Optional; prefix for the output files. If \code{NULL}, no files are saved.
#' @param validate_round_trip Logical; if \code{TRUE}, validates the saved DGEList by loading it back and comparing to the original.
#'
#' @return A merged \code{DGEList} object containing all samples with donor, cell type, region,
#'   and specified metadata columns in \code{dge$samples}.
#'
#' @export
#' @import data.table edgeR
build_merged_dge <- function(manifest_file, metacell_dir, cell_metadata_file, metadata_columns, has_village=T,
                             outDir=NULL, outName=NULL, validate_round_trip=F) {
    # Read cell metadata
    cell_metadata <- data.table::fread(cell_metadata_file, data.table=F)

    # handle mixed assay types
    cell_metadata<- addMixedAssayType(cell_metadata)

    # Add cell type proportions PCA scores to the cell metadata
    # if (!is.null(cellTypeProportionsPCAFile)) {
    #     cell_metadata <- addCellTypeProportionsPCA(cellTypeProportionsPCAFile, cell_metadata, numPCs=4, scalePCs=TRUE)
    # }

    # To enable the python->R filtering, we need to replace NA values in character/factor columns with "NA"
    cell_metadata <- replace_na_strings(cell_metadata)

    manifest <- data.table::fread(manifest_file)

    dge_list <- vector("list", length = nrow(manifest))
    for (i in seq_len(nrow(manifest))) {
        dge_list[[i]] <- process_metacell_file(manifest[i], metacell_dir, cell_metadata, metadata_columns, has_village=has_village)
    }

    dge_list <- Filter(Negate(is.null), dge_list)
    if (length(dge_list) == 0) stop("No valid DGELists generated.")

    # Get union of all gene names
    gene_union <- Reduce(union, lapply(dge_list, function(x) rownames(x$counts)))

    # Pad all DGELists to the union
    dge_list <- lapply(dge_list, pad_dge_to_union, all_genes = gene_union)

    # Combine using edgeR::cbind
    dge_merged <- do.call(cbind, dge_list)

    # write to disk if the outDir and outName are provided
    if (!is.null(outDir) && !is.null(outName)) {
        saveDGEList(dge_merged, dir = outDir, prefix = outName)
        if (validate_round_trip) {
            r=loadDGEList(dir = outDir, prefix = outName)
            compareDGEList(dge_merged, r)
        }

    }



    return(dge_merged)
}



#' Process a Single Metacell Count Matrix
#'
#' Reads a single pseudobulk count matrix, renames its columns using donor, cell type, and region
#' information from the manifest, and adds donor-level metadata from a cell metadata table.
#'
#' This function is intended for internal use in building a combined DGEList object from multiple files.
#'
#' @param manifest_row A single-row \code{data.frame} or \code{data.table} representing one row from the manifest.
#' @param metacell_dir Directory containing the \code{*.metacells.txt.gz} count matrix files.
#' @param cell_metadata A \code{data.frame} containing per-cell metadata used to extract donor-level annotations.
#' @param metadata_columns A character vector of column names to be pulled from \code{cell_metadata} for each donor.
#' @param has_village Logical; if \code{TRUE}, split the column names in the metacell into donor_village.
#'
#' @return A \code{DGEList} object with renamed sample columns and additional metadata in \code{dge$samples}.
#' @import data.table edgeR
#' @keywords internal

#manifest_row=manifest[1];
process_metacell_file <- function(manifest_row, metacell_dir, cell_metadata, metadata_columns, has_village=TRUE) {
    file_path <- file.path(metacell_dir, paste0(manifest_row$output_name, ".metacells.txt.gz"))

    if (!file.exists(file_path)) {
        warning("Missing file: ", file_path)
        return(NULL)
    }

    # Read counts matrix (genes x donors)
    mat <- data.table::fread(file_path)
    if (dim (mat)[2]<2) {
        warning("file didn't have any donors: ", file_path)
        return (NULL)
    }
    gene_names <- mat[[1]]
    mat <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat) <- gene_names

    # Rename columns: donor__celltype__region
    donors <- colnames(mat)
    villages=NULL
    if (has_village) {
        # Split donor names into donor and village
        donor_parts <- strsplit(donors, "_", fixed = TRUE)
        donors <- sapply(donor_parts, function(x) x[1])  # Take the first part as donor
        villages <- sapply(donor_parts, function(x) paste(x[-1], collapse = "_"))  # Join the rest as village
        new_names <- paste0(donors, "__", villages, "__", manifest_row$cell_type_name, "__", manifest_row$region_name)
    } else {
        new_names <- paste0(donors, "__", manifest_row$cell_type_name, "__", manifest_row$region_name)
    }

    colnames(mat) <- new_names

    # Create DGEList
    dge <- edgeR::DGEList(counts = mat)

    # Add sample metadata
    dge$samples$sample_name <- new_names
    dge$samples$donor <- donors
    dge$samples$cell_type <- manifest_row$cell_type_name
    dge$samples$region <- manifest_row$region_name

    if (has_village)
        dge$samples$village <- villages

    #add donor level metadata
    # first subset the cell metadata to the filtered cell type information based on the manifest
    # This is to handle complex cell type filters like "MSN_D1 matrix" which rely on multiple columns.
    cell_metadata_this <- filter_df_auto_df_name(cell_metadata, manifest_row$cell_type_filter_str)
    # Filter by region list
    cell_metadata_this<- filter_by_region_list(df=cell_metadata_this, region_str=manifest_row$region_list)

    #unique (cell_metadata_this$brain_region_abbreviation)
    # d=cell_metadata_this[cell_metadata_this$donor_external_id=="MD6434",]
    # table (d[,c("brain_region_abbreviation", "village")])

    if (nrow(cell_metadata_this) == 0) {
        warning("No matching cell metadata for filter: ", manifest_row$cell_type_filter_str)
        return(dge)
    }
    # Add donor annotations
    dge <- add_donor_annotations(dge, cell_metadata_this, metadata_columns)

    #add the "use" columns.
    useCols=c("MDS", "differential_expression", "eQTL_Analysis")
    for (col in useCols) {
        if (col %in% colnames(manifest_row)) {
            dge$samples[[col]] <- manifest_row[[col]]
        } else {
            dge$samples[[col]] <- NA
        }
    }

    return(dge)
}

#' Add Donor-Level Annotations to a DGEList
#'
#' Merges donor-level metadata into the \code{samples} data frame of a \code{DGEList} object
#' by matching donor identifiers to aggregated values in a cell-level metadata table.
#'
#' @param dge_merged A \code{DGEList} object with a \code{donor} column in \code{dge$samples}.
#' @param cell_metadata A data frame of per-cell metadata used to extract donor-level summaries.
#' @param metadata_columns Character vector of columns to extract and merge into the \code{samples} data frame.
#' @param add_counts Logical; if \code{TRUE}, adds a column \code{n_cells} with the number of nuclei per donor (num_nuclei)
#'
#' @return A \code{DGEList} object with additional donor-level columns in \code{dge$samples}.
#'
#' @keywords internal
add_donor_annotations<-function (dge_merged, cell_metadata, metadata_columns, add_counts=TRUE) {
    #aggregate donor-level metadata
    donor_metadata <- getAnnotations(cell_metadata, aggregationFeatures = "donor_external_id",  additionalColumns=metadata_columns, addCounts = add_counts)

    #if add_counts is true, add that to the metadata_columns
    if (add_counts && !("num_nuclei" %in% metadata_columns)) {
        metadata_columns <- c(metadata_columns, "num_nuclei")
    }

    # Merge with cell metadata
    idx=match(dge_merged$samples$donor, donor_metadata$donor_external_id)
    dge_merged$samples <- cbind(dge_merged$samples, donor_metadata[idx, metadata_columns, drop = FALSE])

    return(dge_merged)
}


#' Pad a DGEList to Include a Union of Genes
#'
#' Ensures that the count matrix in a \code{DGEList} contains all genes in a specified union,
#' adding rows of zeros for missing genes. The resulting matrix is reordered to match
#' \code{all_genes}.
#'
#' @param dge A \code{DGEList} object to be padded.
#' @param all_genes A character vector of gene names representing the union across all datasets.
#'
#' @return A \code{DGEList} object whose \code{counts} matrix includes all genes in \code{all_genes}.
#'
#' @keywords internal
pad_dge_to_union <- function(dge, all_genes) {
    cur_genes <- rownames(dge$counts)
    missing_genes <- setdiff(all_genes, cur_genes)

    if (length(missing_genes) > 0) {
        zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(dge$counts))
        rownames(zero_mat) <- missing_genes
        colnames(zero_mat) <- colnames(dge$counts)
        new_counts <- rbind(dge$counts, zero_mat)
        new_counts <- new_counts[all_genes, , drop = FALSE]
        dge$counts <- new_counts
    } else {
        dge$counts <- dge$counts[all_genes, , drop = FALSE]
    }

    return(dge)
}

#' Aggregate Cell-Level Metadata to Donor-Level Summaries
#'
#' Aggregates one or more columns from a cell-level metadata data frame into donor-level summaries.
#' Numeric columns are averaged, and categorical columns must be consistent within each group.
#'
#' @param metricsDF A data frame with per-cell metadata, including donor ID and columns to aggregate.
#' @param aggregationFeatures A character vector of column names to group by (e.g., \code{"donor_external_id"}).
#' @param additionalColumns A character vector of column names to summarize for each group.
#' @param addCounts Logical; if \code{TRUE}, adds a column \code{n_cells} with the number of cells per group.
#'
#' @return A data frame with one row per donor and columns for the aggregated metadata.
#'
#' @keywords internal
#' @importFrom stats na.omit
getAnnotations <- function(metricsDF,
                                 aggregationFeatures = c("donor_external_id"),
                                 additionalColumns = c("age", "imputed_sex"),
                                 addCounts = FALSE) {

    missing_grouping <- setdiff(aggregationFeatures, colnames(metricsDF))
    missing_additional <- setdiff(additionalColumns, colnames(metricsDF))
    if (length(missing_grouping) > 0) {
        stop("Missing aggregationFeatures in metricsDF: ", paste(missing_grouping, collapse = ", "))
    }
    if (length(missing_additional) > 0) {
        stop("Missing additionalColumns in metricsDF: ", paste(missing_additional, collapse = ", "))
    }

    allColumns=unique(c(aggregationFeatures, additionalColumns))

    # Subset relevant columns
    subsetDF <- metricsDF[, allColumns, drop = FALSE]

    # Split by group key
    split_groups <- split(subsetDF, f = subsetDF[, aggregationFeatures], drop = TRUE)

    # Compute per-group summaries
    result_list <- lapply(split_groups, function(df) {
        row <- df[1, aggregationFeatures, drop = FALSE]

        for (col in additionalColumns) {
            x <- df[[col]]
            if (is.numeric(x)) {
                row[[col]] <- mean(x, na.rm = TRUE)
            } else {
                unique_vals <- unique(stats::na.omit(x))
                if (length(unique_vals) == 1) {
                    row[[col]] <- unique_vals
                } else {
                    stop(sprintf("Non-unique categorical values in column '%s' for group: %s",
                                 col, paste(row[1, aggregationFeatures], collapse = ", ")))
                }
            }
        }

        if (addCounts) {
            row[["num_nuclei"]] <- nrow(df)
        }

        return(row)
    })

    # Combine summaries
    annotationDF <- do.call(rbind, result_list)

    return(annotationDF)
}


#' Filter a Data Frame Using a Python-Style Expression
#'
#' This function applies a Python-style filter string (e.g., one that uses
#' `adata.obs["column"] == "value"` syntax) to an R data frame. It automatically
#' detects the name of the data frame passed in and rewrites the filter expression
#' to use R syntax (e.g., `df$column == "value"`). Logical `NA` values in the
#' evaluated expression are excluded from the result.
#'
#' @param df_expr An unquoted R data frame symbol (e.g., `df` or `my_data`). This is
#'   captured using non-standard evaluation, so do not quote it.
#' @param python_filter_string A character string representing a Python-style filter,
#'   such as that used on `.obs` in `AnnData`. Column references must be of the form
#'   `adata.obs["column_name"]`. Logical operators `&`, `|`, and comparison operators
#'   `==`, `!=`, `>=`, `<=`, `>`, `<` are supported.
#'
#' @return A subset of the original data frame, filtered using the evaluated expression.
#'   Rows with `NA` in the logical index are excluded.
#'
#' @export
# df_expr= cell_metadata; python_filter_string=manifest[26]$cell_type_filter_str
filter_df_auto_df_name <- function(df_expr, python_filter_string) {
    # Capture the unevaluated expression (e.g., df2)
    df_name <- as.character(substitute(df_expr))

    # Evaluate the dataframe itself
    df <- eval(substitute(df_expr), envir = parent.frame())

    # Replace adata.obs["col"] with dfname$col
    expr <- gsub('adata\\.obs\\[["\']([^"\']+)["\']\\]',
                 paste0(df_name, '$\\1'), python_filter_string)

    # Evaluate the logical expression
    idx <- eval(parse(text = expr), envir = parent.frame())

    # Return the subset of the dataframe
    result <- df[which(idx), ]

    return (result)
}

#' Filter a Data Frame Using a Python-Style Expression, incl. .isin()
filter_df_auto_df_name <- function(df_expr, python_filter_string) {
    df_name <- as.character(substitute(df_expr))
    df <- eval(substitute(df_expr), envir = parent.frame())

    expr <- python_filter_string

    # .isin(['A','B']) or .isin(('A','B')) -> %in% c('A','B')
    expr <- gsub("\\.isin\\s*\\(\\s*\\[", " %in% c(", expr, perl = TRUE)
    expr <- gsub("\\.isin\\s*\\(\\s*\\(", " %in% c(", expr, perl = TRUE)
    expr <- gsub("\\]\\s*\\)", ")", expr, perl = TRUE)
    expr <- gsub("\\)\\s*\\)", ")", expr, perl = TRUE)

    # Optional: normalize Python literals
    expr <- gsub("\\bNone\\b", "NA", expr, perl = TRUE)
    expr <- gsub("\\bTrue\\b", "TRUE", expr, perl = TRUE)
    expr <- gsub("\\bFalse\\b", "FALSE", expr, perl = TRUE)

    # adata.obs["col"] -> df$col
    expr <- gsub('adata\\.obs\\[\\s*["\']([^"\']+)["\']\\s*\\]',
                 paste0(df_name, '$\\1'), expr, perl = TRUE)

    idx <- eval(parse(text = expr), envir = parent.frame())
    df[which(idx), , drop = FALSE]
}


#' Filter Cell Metadata by Brain Region Abbreviation
#'
#' Given a comma-separated list of region names (e.g., \code{"NAC,NACc,NACs"}), this function
#' filters the input data frame to include only rows where \code{brain_region_abbreviation}
#' matches any of the listed regions.
#'
#' @param df A data frame containing a column named \code{brain_region_abbreviation}.
#' @param region_str A character string of one or more region names separated by commas.
#'
#' @return A filtered data frame with only matching brain regions.
#'
#' @keywords internal
filter_by_region_list <- function(df, region_str) {
    regions <- strsplit(region_str, ",")[[1]]
    regions <- trimws(regions)
    result=df[df$brain_region_abbreviation %in% regions, , drop = FALSE]
    return (result)
}


#' Convert NA Values in Character or Factor Columns to the String "NA"
#'
#' This function mutates all character or factor columns in a data frame,
#' replacing \code{NA} values with the literal string \code{"NA"}, to match
#' Python-style behavior in filtering expressions.
#'
#' @param df A data frame whose character/factor columns may contain NA values.
#'
#' @return A modified data frame where NA entries in string columns are replaced with "NA".
#'
#' @keywords internal
replace_na_strings <- function(df) {
    for (col in names(df)) {
        if (is.character(df[[col]]) || is.factor(df[[col]])) {
            df[[col]][is.na(df[[col]])] <- "NA"
        }
    }
    return(df)
}

#' Flag donors with multiple assay types in the same region
#'
#' This function identifies donors that are associated with more than one unique
#' sequencing assay in the same brain region and sets their `sequencing_assay`
#' value to `"NextGEM:GEMX"` in the input data frame.
#'
#' @param metricsDF A data frame containing at least the columns
#'   `single_cell_assay`, `donor_external_id`, and `brain_region_abbreviation`.
#'  Each row represents a nuclei for a donor with a corresponding assay type / region
#'
#' @return A data frame identical to `metricsDF` except that for donors associated
#'   with more than one distinct `single_cell_assay` in a region, the `single_cell_assay`
#'   value is set to `"mixed"` in all rows for that donor-region combination.
addMixedAssayType <- function(metricsDF) {
    # Get unique combinations of assay, donor, and region
    distinct_rows <- unique(metricsDF[, c("single_cell_assay", "donor_external_id", "brain_region_abbreviation")])

    # Create a composite label for donor-region
    distinct_rows$donor_region <- paste(distinct_rows$donor_external_id,
                                        distinct_rows$brain_region_abbreviation,
                                        sep = "_")

    # Identify donor-region pairs with more than one assay
    donor_region_counts <- table(distinct_rows$donor_region)
    mixed_donor_regions <- names(donor_region_counts[donor_region_counts > 1])

    # Tag original rows for those mixed donor-regions
    donor_region_all <- paste(metricsDF$donor_external_id,
                              metricsDF$brain_region_abbreviation,
                              sep = "_")

    idx_mixed <- which(donor_region_all %in% mixed_donor_regions)
    if (length(idx_mixed) > 0) {
        metricsDF$single_cell_assay[idx_mixed] <- "NextGEM:GEMX"
    }

    return(metricsDF)
}

#' Parse the cell type proportions PCA file and add to the metrics DF.
#'
#' This function reads a PCA file containing cell type proportions,
#' scales the PCA scores if requested, and adds them to the metrics data frame.
#'
#' @param cellTypeProportionsPCAFile Path to the PCA file containing cell type proportions.
#' @param metricsDF A data frame containing cell-level metadata, which will be updated with PCA scores.
#' @param numPCs The number of principal components to retain from the PCA file (default is 4).
#' @param scalePCs Logical; if \code{TRUE}, scales the PCA scores to have mean 0 and standard deviation 1 (default is \code{FALSE}).
#' @return A data frame identical to \code{metricsDF} but with additional columns for the PCA scores of cell type proportions.
#' @import logger
#' @export
addCellTypeProportionsPCA<-function (cellTypeProportionsPCAFile, metricsDF, numPCs=4, scalePCs=FALSE) {
    logger::log_info("Adding cell type proportions PCA scores to cell level metadata")
    a=utils::read.table(cellTypeProportionsPCAFile, header=T, sep="\t")
    #TODO REMOVE THIS HACK
    #a=a[,c(1,2,4)]

    idxPC=grep("PC", colnames(a))
    idxNotPCs=setdiff(1:ncol(a), idxPC)
    colnames(a)[idxPC]=paste("cell_proportion", colnames(a)[idxPC], sep="_")

    if (scalePCs) {
        a[,idxPC] <- scale(a[,idxPC], center = TRUE, scale = TRUE)
    }

    #select the first X PCs.
    if (numPCs > length(idxPC)) {
        stop(paste("Requested", numPCs, "PCs but only found", length(idxPC), "in the file."))
    }

    idxPC=idxPC[1:numPCs]
    a=a[,c(idxNotPCs,idxPC), drop=FALSE]

    #a compound key to match the rownames in a.
    metricsDF$key=paste(metricsDF$donor_external_id, metricsDF$brain_region_abbreviation, sep="_")

    idx=match(metricsDF$key, rownames(a))
    if (any (is.na(idx))) {
        logger::log_warn(paste("Missing cell type proportions PCA for donors [", paste(metricsDF$donor_external_id[is.na(idx)], collapse=","), "]"))
    }
    #add the PCA scores to the metricsDF
    metricsDF=cbind(metricsDF, a[idx, idxPC, drop=FALSE])
    return (metricsDF)
}

# Functions to serialize and deserialize edgeR::DGEList objects
# using two human-readable TSV files (gzipped by default):
# one for counts (with optional gene annotation) and one for sample metadata,
# each prefixed by a user-defined string.

#' Save a DGEList to disk as two gzipped TSV files with a filename prefix
#'
#' @param dge A DGEList object from edgeR
#' @param dir Directory path where the files will be written (will be created if needed)
#' @param prefix Character prefix for the filenames (default "dge").
#'               Will produce files: prefix_counts.tsv.gz and prefix_samples.tsv.gz
#' @export
#' @importFrom utils write.table
saveDGEList <- function(dge, dir, prefix = "dge") {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    # Prepare counts with optional gene annotation
    counts_df <- as.data.frame(dge$counts, check.names = FALSE)
    if (!is.null(dge$genes)) {
        genes_df <- as.data.frame(dge$genes, stringsAsFactors = FALSE)
        counts_df <- cbind(genes_df[rownames(counts_df), , drop = FALSE], counts_df)
    }
    # Add explicit GeneID column
    counts_df <- cbind(GeneID = rownames(counts_df), counts_df)

    # Write gzipped counts TSV
    counts_file <- file.path(dir, paste0(prefix, "_counts.tsv.gz"))
    gz_counts_con <- gzfile(counts_file, "w")
    utils::write.table(counts_df,
                file = gz_counts_con,
                sep = "\t", quote = FALSE, row.names = FALSE)
    close(gz_counts_con)

    # Prepare sample metadata
    samples_df <- as.data.frame(dge$samples, stringsAsFactors = FALSE, check.names = FALSE)
    # Drop lib.size and norm.factors to avoid persisting old values
    samples_df <- samples_df[, setdiff(colnames(samples_df), c("lib.size", "norm.factors")), drop = FALSE]
    samples_df <- cbind(SampleID = rownames(dge$samples), samples_df)

    # Write gzipped samples TSV
    samples_file <- file.path(dir, paste0(prefix, "_samples.tsv.gz"))
    gz_samples_con <- gzfile(samples_file, "w")
    utils::write.table(samples_df,
                file = gz_samples_con,
                sep = "\t", quote = FALSE, row.names = FALSE)
    close(gz_samples_con)
}

#' Load a DGEList from disk saved by saveDGEList with a filename prefix
#'
#' @param dir Directory path containing the files
#' @param prefix Character prefix used when saving (default "dge").
#' @return A DGEList object
#' @export
loadDGEList <- function(dir, prefix = "dge") {
    # Construct filenames for gzipped TSVs
    counts_file  <- file.path(dir, paste0(prefix, "_counts.tsv.gz"))
    samples_file <- file.path(dir, paste0(prefix, "_samples.tsv.gz"))

    # Read in the tables via gzfile connections
    gz_counts_con  <- gzfile(counts_file, "r")
    counts_df      <- utils::read.delim(gz_counts_con, stringsAsFactors = FALSE, check.names = FALSE)
    close(gz_counts_con)

    gz_samples_con <- gzfile(samples_file, "r")
    samples_df     <- utils::read.delim(gz_samples_con, stringsAsFactors = FALSE, check.names = FALSE)
    close(gz_samples_con)

    # Extract sample IDs
    sample_ids <- samples_df$SampleID

    # Identify count columns by matching sample IDs
    count_cols <- intersect(names(counts_df), sample_ids)
    counts_mat <- as.matrix(counts_df[, count_cols, drop = FALSE])
    rownames(counts_mat) <- counts_df$GeneID

    # Extract gene annotation columns if present
    gene_cols <- setdiff(names(counts_df), c("GeneID", count_cols))
    genes_df <- if (length(gene_cols) > 0) counts_df[, gene_cols, drop = FALSE] else NULL
    if (!is.null(genes_df)) rownames(genes_df) <- counts_df$GeneID

    # Prepare samples data.frame: drop lib.size/norm.factors if present, then set rownames
    samples_df <- samples_df[, setdiff(colnames(samples_df), c("lib.size", "norm.factors", "SampleID")), drop = FALSE]
    rownames(samples_df) <- sample_ids

    # Construct and return DGEList; lib.size and norm.factors will be recomputed
    edgeR::DGEList(
        counts  = counts_mat,
        genes   = genes_df,
        samples = samples_df
    )
}

#' Compare two DGEList objects for near equality
#'
#' @param dge1 A DGEList object
#' @param dge2 A DGEList object (e.g., loaded)
#' @param tol Numeric tolerance for comparing numeric values (default 1e-8)
#' @return TRUE if equal within tolerance, otherwise throws an error describing the mismatch
#' @export
compareDGEList <- function(dge1, dge2, tol = 1e-8) {
    # Counts: dimensions and names
    if (!all(dim(dge1$counts) == dim(dge2$counts))) stop("Counts dimensions differ")
    if (!all(rownames(dge1$counts) == rownames(dge2$counts))) stop("Counts rownames differ")
    if (!all(colnames(dge1$counts) == colnames(dge2$counts))) stop("Counts colnames differ")
    # Numeric compare with tolerance, treat NA==NA as equal
    diff_counts <- abs(dge1$counts - dge2$counts)
    na_both <- is.na(dge1$counts) & is.na(dge2$counts)
    mismatch <- (diff_counts > tol) & !na_both
    if (any(mismatch, na.rm = TRUE)) stop("Counts differ beyond tolerance")

    # Genes slot
    if (!is.null(dge1$genes) || !is.null(dge2$genes)) {
        if (is.null(dge1$genes) || is.null(dge2$genes)) stop("One genes slot is NULL and the other is not")
        if (!all(dim(dge1$genes) == dim(dge2$genes))) stop("Genes dimensions differ")
        if (!all(rownames(dge1$genes) == rownames(dge2$genes))) stop("Genes rownames differ")
        if (!all(colnames(dge1$genes) == colnames(dge2$genes))) stop("Genes colnames differ")
        for (col in colnames(dge1$genes)) {
            v1 <- dge1$genes[[col]]
            v2 <- dge2$genes[[col]]
            if (is.numeric(v1) && is.numeric(v2)) {
                diff <- abs(v1 - v2)
                na_both_g <- is.na(v1) & is.na(v2)
                if (any((diff > tol) & !na_both_g, na.rm = TRUE)) stop(paste("Gene annotation", col, "differs"))
            } else {
                neq <- !(v1 == v2) & !(is.na(v1) & is.na(v2))
                if (any(neq, na.rm = TRUE)) stop(paste("Gene annotation", col, "differs"))
            }
        }
    }

    # Samples slot
    if (!all(dim(dge1$samples) == dim(dge2$samples))) stop("Samples dimensions differ")
    if (!all(rownames(dge1$samples) == rownames(dge2$samples))) stop("Samples rownames differ")
    if (!all(colnames(dge1$samples) == colnames(dge2$samples))) stop("Samples colnames differ")
    for (col in colnames(dge1$samples)) {
        v1 <- dge1$samples[[col]]
        v2 <- dge2$samples[[col]]
        if (is.numeric(v1) && is.numeric(v2)) {
            diff <- abs(v1 - v2)
            na_both_s <- is.na(v1) & is.na(v2)
            if (any((diff > tol) & !na_both_s, na.rm = TRUE)) stop(paste("Sample column", col, "differs"))
        } else {
            neq <- !(v1 == v2) & !(is.na(v1) & is.na(v2))
            if (any(neq, na.rm = TRUE)) stop(paste("Sample column", col, "differs"))
        }
    }
    # If we reach here, everything matches
    logger::log_info("DGELists are equal within tolerance of ", tol)
    return(TRUE)
}
