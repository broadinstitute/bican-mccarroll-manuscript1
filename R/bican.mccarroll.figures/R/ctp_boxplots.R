
#' Regroup and Aggregate Cell Type Proportions
#'
#' @description
#' Aggregates cell type proportions based on a user-defined mapping list. 
#' This is particularly useful for collapsing fine-grained annotations 
#'
#' @param long_ctp_df A data frame in long format containing cell type proportions.
#'   Must contain sample_id, donor_external_id, brain_region_abbreviation_simple,
#'   and frac_nuclei column
#' @param cell_type_column Character string specifying the column name that 
#'   contains the cell type annotations.
#' @param groupings A named list where names are the new aggregate labels and 
#'   values are character vectors of the original cell types to be summed.
#'
#' @return A consolidated data frame (tibble) containing the aggregated 
#'   proportions for each sample and region.
#' 
#' @details 
#' The function uses \code{purrr::imap_dfr} to iterate over the mapping list. 
#' Note that any cell types present in the input data but \emph{not} included 
#' in the \code{groupings} list will be dropped from the resulting data frame.
#'
#' @export
#'
#' @importFrom dplyr filter group_by summarise bind_rows
#' @importFrom purrr imap_dfr
#' @importFrom rlang sym :=
#'
#' @examples
#' \dontrun{
#' mappings <- list(
#'   "Glia" = c("astrocyte", "microglia"),
#'   "Neurons" = c("excitatory", "inhibitory")
#' )
#' regroup_cell_types(my_data, "annotation", mappings)
#' }
regroup_cell_types <- function(long_ctp_df, cell_type_column, groupings) {
  
  grouped_df <- purrr::imap_dfr(
    groupings,
    ~ {
      # Identify the target column dynamically using rlang
      target_col <- rlang::sym(cell_type_column)
      
      long_ctp_df |>
        # Filter for the original subtypes mapped to this group (.x)
        dplyr::filter(!!target_col %in% .x) |>
        # Maintain experimental design metadata
        dplyr::group_by(
          sample_id, 
          donor_external_id, 
          brain_region_abbreviation_simple
        ) |>
        # Sum fractions and rename the annotation to the group label (.y)
        dplyr::summarise(
          !!target_col := .y,
          fraction_nuclei = sum(fraction_nuclei, na.rm = TRUE),
          .groups = "drop"
        )
    }
  )
  
  return(grouped_df)
}




#' Create Boxplots for Cell Type Proportions given a CTP Data Frame
#'
#' @description
#' Generates a grouped boxplot visualizing cell type fractions across different 
#' experimental conditions (e.g., brain regions). The function handles dynamic 
#' column mapping, optional factor reordering, and aesthetic customizations 
#' suitable for multi-panel figures.
#'
#' @param ctp_df A data frame in long format. Must contain a numeric column 
#'   \code{fraction_nuclei} and columns specified by \code{cell_type_column} 
#'   and \code{fill_feature}.
#' @param cell_type_column Character string. The column name in \code{ctp_df} 
#'   representing cell type identities (plotted on the x-axis).
#' @param fill_feature Character string. The column name used for boxplot 
#'   grouping and fill color. Defaults to \code{"brain_region_abbreviation_simple"}.
#' @param cell_type_label_map A named character vector for x-axis labels. 
#'   The names should match values in \code{cell_type_column} and define the 
#'   plotting order (factor levels).
#' @param region_label_map A named character vector for legend labels. 
#'   The names should match values in \code{fill_feature} and define the 
#'   fill order.
#' @param vline_int Numeric. Optional x-axis intercept(s) to draw dashed 
#'   vertical separator lines (e.g., to separate glia from neurons).
#' @param title Character string. The plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @details 
#' This function employs tidy evaluation via \code{rlang}. If label maps are 
#' provided, the respective variables are converted to factors using the 
#' names of the map as levels, ensuring a deterministic categorical order 
#' in the final visualization.
#'
#' @export
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot position_dodge geom_vline 
#'   scale_x_discrete scale_fill_discrete labs theme_bw theme element_text 
#'   element_blank ggtitle
#' @importFrom rlang sym
make_ctp_boxplot <- function(ctp_df, 
                             cell_type_column, 
                             fill_feature = "brain_region_abbreviation_simple",
                             cell_type_label_map = NULL,
                             region_label_map = NULL,
                             vline_int = NULL,
                             title = "Cell Type Proportions as a Fraction of All Nuclei") {
  
  # Ensure regions follow the provided map order
  if (!is.null(region_label_map)) {
    ctp_df[[fill_feature]] <- factor(
      ctp_df[[fill_feature]], 
      levels = names(region_label_map)
    )
  }
  
  # Ensure cell types follow the provided map order
  if (!is.null(cell_type_label_map)) {
    ctp_df[[cell_type_column]] <- factor(
      ctp_df[[cell_type_column]], 
      levels = names(cell_type_label_map)
    )
  }
  
  # Construct the plot
  ctp_boxplot <- ggplot2::ggplot(
    ctp_df,
    ggplot2::aes(
      x = !!rlang::sym(cell_type_column),
      y = fraction_nuclei,
      fill = !!rlang::sym(fill_feature)
    )
  ) + 
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.8),
      outlier.size = 0.5
    ) +
    # Conditional layers added via list evaluation
    list(
      if (!is.null(vline_int)) {
        ggplot2::geom_vline(xintercept = vline_int, linetype = "dashed", color = "black")
      },
      if (!is.null(cell_type_label_map)) {
        ggplot2::scale_x_discrete(labels = cell_type_label_map)
      },
      if (!is.null(region_label_map)) {
        ggplot2::scale_fill_discrete(labels = region_label_map)
      }
    ) +
    ggplot2::labs(
      x = "",
      y = "Fraction of nuclei",
      fill = "Region"
    ) +
    # base_size is set here to avoid the theme hierarchy warning
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(title) 
  
  return(ctp_boxplot)
}




#' Generate Standardized Cell Type Proportion Plots for manuscript figure
#'
#' @description
#' A high-level wrapper function that executes the full visualization pipeline:
#' loading raw proportion data, aggregating specific neuronal and glial classes,
#' and generating a publication-ready boxplot. This function is specifically
#' configured for the BICAN UM1 McCarroll analysis freeze.
#'
#' @return A \code{ggplot} object visualizing cell type fractions across 
#'   Striatal and DFC regions.
#' 
#' @details 
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads Level 1 cell type proportions from the Broad BICAN directory.
#'   \item Identifies and filters "Out of Domain" (OOD) neurons.
#'   \item Groups fine-grained annotations into broad classes: Interneurons 
#'         (GABAergic + Striatal Cholinergic) and Projection Neurons (MSNs + IT + L5/6 Glutamatergic).
#'   \item Retains Glial classes (Oligodendrocytes, Astrocytes, OPC, Microglia) 
#'         as individual categories.
#'   \item Applies custom region and cell type labels with pre-defined plotting order.
#' }
#' 
#' @note This function contains hard-coded file paths and specific cell-type 
#'   naming conventions relevant to the CAP freeze 3 analysis.
#'
#' @export
#' 
#' @importFrom data.table fread
#' @importFrom rlang set_names
make_region_ctp_plot_pretty <- function() {
  
  # --- Data Loading ---
  ctp_dir <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/cell_type_proportions/LEVEL_1"
  ctp_annotations <- data.table::fread(
    file.path(ctp_dir, "donor_region.annotation.cell_type_proportions.txt"),
  )
  
  cell_types <- unique(ctp_annotations$annotation)
  
  # --- Taxonomy Logic ---
  # Define neurons that do not fit the primary Interneuron/Projection groupings
  OOD_neurons <- c("VTR-HTH_glut", "SN-VTR-HTH_GABA", "other_GABA")
  
  # Interneurons: All GABAergic types plus striatal cholinergic, excluding OOD
  interneurons <- c(grep("GABA", cell_types, value = TRUE),
                    "striatal_cholinergic") |> setdiff(OOD_neurons)
  
  # Projection Neurons: MSNs, IT (Intratelencephalic), and specific Cortical Glutamatergic layers
  projection_neurons <- c(grep("MSN", cell_types, value = TRUE),
                          grep("IT", cell_types, value = TRUE),
                          "cortical_glutamatergic_L5ET",
                          "cortical_glutamatergic_L56NP",
                          "cortical_glutamatergic_L6")
  
  # --- Grouping Configuration ---
  neuron_groupings <- list(
    "interneurons" = interneurons,
    "projection_neurons" = projection_neurons
  )
  
  # Idiomatically convert glia vector to a named list for regrouping
  glia_of_interest <- c("oligodendrocyte", "astrocyte", "OPC", "microglia") |> 
    rlang::set_names() |> 
    as.list()
  
  groupings <- c(glia_of_interest, neuron_groupings)
  
  # --- Processing & Plotting ---
  # Aggregate raw counts/fractions into the defined groups
  ctp_df <- regroup_cell_types(
    long_ctp_df = ctp_annotations,
    cell_type_column = "annotation",
    groupings = groupings
  )
  
  # Define display labels and factor ordering for the legend
  region_label_map <- c(
    "CaH" = "CaH (n=XX)",
    "Pu"  = "Pu (n=XX)",
    "NAC" = "NAC (n=XX)",
    "ic"  = "ic (n=XX)",
    "DFC" = "DFC (n=XX)"
  )
  
  # Define display labels and factor ordering for the x-axis
  cell_type_labels <- c(
    "oligodendrocyte"    = "Oligodendrocyte",
    "astrocyte"          = "Astrocyte",
    "OPC"                = "OPC",
    "microglia"          = "Microglia",
    "projection_neurons" = "Projection Neurons",
    "interneurons"       = "Interneurons"
  )
  
  # Generate the final boxplot
  make_ctp_boxplot(
    ctp_df = ctp_df,
    cell_type_column = "annotation",
    fill_feature = "brain_region_abbreviation_simple",
    cell_type_label_map = cell_type_labels,
    region_label_map = region_label_map,
    vline_int = 4.5, # Vertical line to separate Glia from Neurons
    title = "Cell Type Proportions as a Fraction of All Nuclei in Striatum and DFC"
  )
}

make_region_ctp_plot_pretty()













