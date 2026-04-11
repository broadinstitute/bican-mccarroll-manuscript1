# MIT License
#
# Copyright 2026 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Make a less verbose version of coloc.abf, and get the result from $results
coloc.abf.quiet <- purrr::quietly(coloc::coloc.abf)

#' Run coloc for the gene-variant pairs in the tensorQTL pairs VCF file and the GWAS variants in the GWAS VCF file,
#' and write the results to the output files.
#'
#' @param eqtl_file Required path to the tensorQTL eQTL file.
#' @param pairs_vcf_file Required path to the pairs VCF file.
#' @param gwas_vcf_file Required path to the GWAS VCF file.
#' @param coloc_summary_file Required path to the output coloc summary file.
#' @param coloc_results_file Required path to the output coloc results file.
#' @param eqtl_qval_threshold Optional threshold for filtering eQTLs by q-value. Default is 0.05.
#' @param summary_threshold Optional threshold for filtering coloc results by PP.H4.abf. Default is 0.7.
#' @param variant_threshold Optional threshold for filtering coloc results by SNP.PP.H4. Default is 0.01.
#' @param cumulative_threshold Optional threshold for filtering coloc results by cumulative SNP.PP.H4. Default is 0.95.
#' @param chrs Optional vector of chromosomes to process. Default is NULL, which processes all chromosomes (chr1-chr22 and chrX).
#' @param bcftools_path Optional path to the bcftools executable. Default is "bcftools", which assumes it is in the system PATH.
#'
#' @returns Invisibly returns NULL.
#'
#' @importFrom coloc coloc.abf
#' @importFrom data.table fread
#' @importFrom dplyr arrange bind_rows desc distinct filter group_by if_all mutate n select slice_min ungroup
#' @importFrom purrr quietly
#' @importFrom readr col_character col_double cols_only read_tsv write_tsv
#' @importFrom rlang .data .env
#' @importFrom tibble tibble
#' @importFrom tidyr separate_wider_delim
#' @importFrom tidyselect everything
#'
#' @export
coloc_pairs_gwas <- function(
  eqtl_file,
  pairs_vcf_file,
  gwas_vcf_file,
  coloc_summary_file,
  coloc_results_file,
  eqtl_qval_threshold = 0.05,
  summary_threshold = 0.7,
  variant_threshold = 0.01,
  cumulative_threshold = 0.95,
  chrs = NULL,
  bcftools_path = "bcftools"
) {
  assert_file_is_readable(eqtl_file)
  assert_file_is_readable(pairs_vcf_file)
  assert_file_is_readable(gwas_vcf_file)
  assert_file_is_writable(coloc_summary_file)
  assert_file_is_writable(coloc_results_file)

  if (is.null(chrs) || is.na(chrs)) {
    chrs <- paste0("chr", c(1:22, "X"))
  }

  eqtl_sites_df <- read_eqtl_sites_df(eqtl_file, eqtl_qval_threshold)
  write_coloc_eqtl_pp_h4(
    eqtl_sites_df,
    pairs_vcf_file,
    gwas_vcf_file,
    coloc_summary_file,
    coloc_results_file,
    summary_threshold,
    variant_threshold,
    cumulative_threshold,
    chrs,
    bcftools_path
  )

  return(invisible())
}

assert_file_is_readable <- function(file) {
  if (!file.exists(file)) {
    stop(sprintf("Error: File '%s' does not exist", file))
  }
  if (!file.access(file, 4) == 0) {
    stop(sprintf("Error: File '%s' is not readable", file))
  }
}

assert_file_is_writable <- function(file) {
  if (file.exists(file) && !file.access(file, 2) == 0) {
    stop(sprintf("Error: File '%s' is not writable", file))
  }
  dir <- dirname(file)
  if (!dir.exists(dir)) {
    stop(sprintf("Error: Directory '%s' does not exist", dir))
  }
  if (!file.access(dir, 2) == 0) {
    stop(sprintf("Error: Directory '%s' is not writable", dir))
  }
}


# Read eQTL sites

read_eqtl_sites_df <- function(eqtl_file, eqtl_qval_threshold) {
  eqtl_file |>
    read_tsv(
      col_types = cols_only(
        phenotype_id = col_character(),
        variant_id = col_character(),
        qval = col_double()
      )
    ) |>
    # Remove phenotypes with multiple variants.
    group_by(.data$phenotype_id) |>
    dplyr::filter(n() == 1) |>
    ungroup() |>
    dplyr::filter(.data$qval < eqtl_qval_threshold) |>
    separate_wider_delim(
      "variant_id",
      names = c("chr", "pos", "allele1", "allele2"), delim = ":",
      cols_remove = FALSE
    ) |>
    dplyr::select("phenotype_id", "chr") |>
    distinct()
}

# Read the GWAS data

read_gwas_vcf <- function(gwas_vcf_file, chr, bcftools_path) {
  full_query <- "%CHROM:%POS:%REF:%ALT\t%POS\t[%ES]\t[%SE]\t[%LP]"
  col_names <- c("variant_id", "pos", "effect_size", "standard_error", "pvalue_nlt")
  # col_classes <- c('character', 'integer', 'double', 'double', 'double')
  cmd_line <- paste0(
    bcftools_path, " norm --regions ", chr, " --multiallelics -any --output-type u ", gwas_vcf_file, " | ",
    bcftools_path, " query --format '", full_query, "'",
    # Write a result that shouldn't parse, so we can see if the command failed vs. returning no results
    " || echo !!!failed!!!"
  )
  vcf_positions <- fread(
    cmd = cmd_line,
    header = FALSE,
    col.names = col_names,
    na.strings = ".",
    # Could not get this to work with colClasses, so using mutate below instead
    # colClasses = col_classes,
    showProgress = TRUE
  )
  if (length(vcf_positions) == 0) {
    return(
      tibble(
        variant_id = character(),
        pos = integer(),
        effect_size = double(),
        standard_error = double(),
        pvalue = double()
      )
    )
  }
  vcf_positions <-
    vcf_positions |>
    # remove rows with NA values in any of the columns
    dplyr::filter(if_all(everything(), ~ !is.na(.))) |>
    mutate(
      pvalue = 10^(-.data$pvalue_nlt),
      effect_size = as.double(.data$effect_size),
      standard_error = as.double(.data$standard_error),
      pvalue_nlt = NULL
    )
  vcf_positions
}

# Read the pairs data

read_pairs_vcf <- function(pairs_vcf_file, chr, bcftools_path) {
  full_query <- "%CHROM:%POS:%REF:%ALT\t%POS\t%INFO/phenotype_id\t%INFO/slope\t%INFO/slope_se\t%INFO/pval_nominal"
  col_names <- c("variant_id", "pos", "phenotype_id", "effect_size", "standard_error", "pvalue")
  # col_classes <- c('character', 'integer', 'double', 'double', 'double')
  bcftools_command <- paste0(
    bcftools_path, " query --regions ", chr, " --format '", full_query, "' ", pairs_vcf_file,
    # Write a result that shouldn't parse, so we can see if the command failed vs. returning no results
    " || echo !!!failed!!!"
  )
  vcf_positions <- fread(
    cmd = bcftools_command,
    header = FALSE,
    col.names = col_names,
    na.strings = ".",
    # Could not get this to work with colClasses, so using mutate below instead
    # colClasses = col_classes,
    showProgress = TRUE
  )
  if (length(vcf_positions) == 0) {
    return(
      tibble(
        variant_id = character(),
        pos = integer(),
        phenotype_id = character(),
        effect_size = double(),
        standard_error = double(),
        pvalue = double()
      )
    )
  }
  vcf_positions <-
    vcf_positions |>
    # remove rows with NA values in any of the columns
    dplyr::filter(if_all(everything(), ~ !is.na(.))) |>
    mutate(
      pvalue = as.double(.data$pvalue),
      effect_size = as.double(.data$effect_size),
      standard_error = as.double(.data$standard_error)
    )
  vcf_positions
}

# Setup the coloc dataset for the gene variant pairs

get_gwas_pairs_dfs <- function(pairs_df, gwas_df, phenotype_id) {
  pairs_window_df <-
    pairs_df |>
    dplyr::filter(.data$phenotype_id == .env$phenotype_id) |>
    # remove rows with NA values in any of the columns
    dplyr::filter(if_all(everything(), ~ !is.na(.))) |>
    # remove rows with no standard error that produce NA in coloc.abf
    dplyr::filter(.data$standard_error != 0)

  if (nrow(pairs_window_df) == 0) {
    return(NULL)
  }

  pos_list <- pairs_window_df$pos
  pos_min <- min(pos_list)
  pos_max <- max(pos_list)

  gwas_window_df <-
    gwas_df |>
    dplyr::filter(.data$pos >= pos_min & .data$pos <= pos_max) |>
    dplyr::filter(.data$variant_id %in% pairs_window_df$variant_id) |>
    # remove rows with no standard error that produce NA in coloc.abf
    dplyr::filter(.data$standard_error != 0) |>
    # Some GWAS datasets have multiple variants for the same position, so we take the one with the lowest p-value
    # Ex: PGC3_SCZ_wave3.primary.public.v3.hg38.bcf chr6:61307790
    group_by(.data$variant_id) |>
    slice_min(.data$pvalue, with_ties = FALSE) |>
    ungroup()

  return(list(
    pairs_window_df = pairs_window_df,
    gwas_window_df = gwas_window_df
  ))
}

get_coloc_phentotype <- function(pairs_df, gwas_df, phenotype_id) {
  gwas_pairs_dfs <- get_gwas_pairs_dfs(pairs_df, gwas_df, phenotype_id)
  if (is.null(gwas_pairs_dfs)) {
    return(NULL)
  }

  pairs_window_df <- gwas_pairs_dfs$pairs_window_df
  gwas_window_df <- gwas_pairs_dfs$gwas_window_df

  variant_count <- length(gwas_window_df$variant_id)
  message(sprintf("%d variants found for %s in GWAS data", variant_count, phenotype_id))
  if (length(gwas_window_df$variant_id) == 0) {
    return(NULL)
  }

  gwas_dataset <-
    list(
      beta = gwas_window_df$effect_size,
      varbeta = gwas_window_df$standard_error^2,
      snp = gwas_window_df$variant_id,
      position = gwas_window_df$pos,
      pvalues = gwas_window_df$pvalue,
      type = "cc"
    )

  pairs_dataset <-
    list(
      beta = pairs_window_df$effect_size,
      varbeta = pairs_window_df$standard_error^2,
      snp = pairs_window_df$variant_id,
      position = pairs_window_df$pos,
      pvalues = pairs_window_df$pvalue,
      # This youtube mentioned that tensorQTL is normalizing the MAF so it's ok to set sdY to 1?
      # https://youtu.be/tty5cdznLTE?t=1856
      sdY = 1, # NOTE: Set sdY to 1.
      type = "quant"
    )

  coloc_result <-
    coloc.abf.quiet(
      dataset1 = pairs_dataset,
      dataset2 = gwas_dataset
    )

  return(coloc_result$result)
}

get_coloc_gene_pp_h4 <- function(
  pairs_df,
  gwas_df,
  phenotype_id,
  summary_threshold,
  variant_threshold,
  cumulative_threshold
) {
  empty_summary <- tibble(phenotype_id = character(), pp_h4 = double())
  empty_results <- tibble(phenotype_id = phenotype_id, variant_id = character(), snp_pp_h4 = double())

  coloc_result <- get_coloc_phentotype(pairs_df, gwas_df, phenotype_id)
  if (is.null(coloc_result)) {
    return(
      list(
        summary = empty_summary,
        results = empty_results
      )
    )
  }

  pp_h4 <- coloc_result$summary[["PP.H4.abf"]]

  summary_df <- tibble(phenotype_id = phenotype_id, pp_h4 = pp_h4)

  if (pp_h4 < summary_threshold) {
    return(
      list(
        summary = summary_df,
        results = empty_results
      )
    )
  }

  results_df <- tibble(
    phenotype_id = phenotype_id,
    variant_id = coloc_result$results$snp,
    snp_pp_h4 = coloc_result$results$SNP.PP.H4
  ) |>
    arrange(desc(.data$snp_pp_h4)) |>
    mutate(cum_pp_h4 = cumsum(.data$snp_pp_h4)) |>
    dplyr::filter(.data$cum_pp_h4 <= cumulative_threshold | .data$snp_pp_h4 > variant_threshold) |>
    dplyr::select(-"cum_pp_h4")

  return(
    list(
      summary = summary_df,
      results = results_df
    )
  )
}

# loop through the genes in the eQTL sites and run coloc for each gene

get_coloc_pp_h4_chr <- function(
  eqtl_sites_df,
  pairs_vcf_file,
  gwas_vcf_file,
  summary_threshold,
  variant_threshold,
  cumulative_threshold,
  chr,
  bcftools_path
) {
  chr_sites_df <- eqtl_sites_df |>
    dplyr::filter(.data$chr == .env$chr) |>
    arrange(phenotype_id)
  pairs_df <- read_pairs_vcf(pairs_vcf_file, chr, bcftools_path)
  gwas_df <- read_gwas_vcf(gwas_vcf_file, chr, bcftools_path)

  summary_dfs <- list()
  results_dfs <- list()
  for (phenotype_id in chr_sites_df$phenotype_id) {
    coloc_result <- get_coloc_gene_pp_h4(
      pairs_df,
      gwas_df,
      phenotype_id,
      summary_threshold,
      variant_threshold,
      cumulative_threshold
    )
    summary_dfs[[phenotype_id]] <- coloc_result$summary
    results_dfs[[phenotype_id]] <- coloc_result$results
  }
  return(
    list(
      summary = bind_rows(summary_dfs),
      results = bind_rows(results_dfs)
    )
  )
}

get_coloc_eqtl_pp_h4 <- function(
  eqtl_sites_df,
  pairs_vcf_file,
  gwas_vcf_file,
  summary_threshold,
  variant_threshold,
  cumulative_threshold,
  chrs,
  bcftools_path
) {
  eqtl_chrs <- intersect(chrs, eqtl_sites_df$chr)
  summary_dfs <- list()
  results_dfs <- list()
  for (chr in eqtl_chrs) {
    message(sprintf("Processing %s", chr))
    chr_result <- get_coloc_pp_h4_chr(
      eqtl_sites_df,
      pairs_vcf_file,
      gwas_vcf_file,
      summary_threshold,
      variant_threshold,
      cumulative_threshold,
      chr,
      bcftools_path
    )
    summary_dfs[[chr]] <- chr_result$summary
    results_dfs[[chr]] <- chr_result$results
  }
  return(
    list(
      summary = bind_rows(summary_dfs),
      results = bind_rows(results_dfs)
    )
  )
}

write_coloc_eqtl_pp_h4 <- function(
  eqtl_sites_df,
  pairs_vcf_file,
  gwas_vcf_file,
  coloc_summary_file,
  coloc_results_file,
  summary_threshold,
  variant_threshold,
  cumulative_threshold,
  chrs,
  bcftools_path
) {
  coloc_result <- get_coloc_eqtl_pp_h4(
    eqtl_sites_df,
    pairs_vcf_file,
    gwas_vcf_file,
    summary_threshold,
    variant_threshold,
    cumulative_threshold,
    chrs,
    bcftools_path
  )
  coloc_result$summary |> write_tsv(file = coloc_summary_file)
  coloc_result$results |> write_tsv(file = coloc_results_file)
}
