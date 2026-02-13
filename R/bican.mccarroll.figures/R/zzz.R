.onLoad <- function(libname, pkgname) {
    default_root <- "/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis"
    if (is.null(getOption("bican.mccarroll.figures.data_root_dir"))) {
        options(bican.mccarroll.figures.data_root_dir = default_root)
    }
}
