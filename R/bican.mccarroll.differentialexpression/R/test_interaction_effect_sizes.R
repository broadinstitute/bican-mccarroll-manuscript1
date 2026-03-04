plot_interaction_sanity <- function(interaction_path,
                                    subset_one_path,
                                    subset_two_path) {

    read_de <- function(path) {
        x <- read.table(path, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, row.names = 1)
        x$gene <- rownames(x)
        x
    }

    interaction <- read_de(interaction_path)
    subset_one  <- read_de(subset_one_path)
    subset_two  <- read_de(subset_two_path)

    merged <- merge(subset_one[, c("gene","logFC")],
                    subset_two[, c("gene","logFC")],
                    by="gene",
                    suffixes=c("_base","_comp"))

    merged <- merge(merged,
                    interaction[, c("gene","logFC")],
                    by="gene")

    names(merged)[names(merged)=="logFC"] <- "logFC_interaction"

    merged$expected <- merged$logFC_comp - merged$logFC_base

    r <- cor(merged$expected, merged$logFC_interaction, use="complete.obs")

    p <- ggplot2::ggplot(merged,
                         ggplot2::aes(x=expected, y=logFC_interaction)) +
        ggplot2::geom_point(alpha=0.3) +
        ggplot2::geom_abline(slope=1, intercept=0, color="red") +
        ggplot2::labs(
            x="Expected slope difference (comparison - baseline)",
            y="Interaction logFC",
            title=paste("Interaction sanity check (R =", round(r,3),")")
        )

    print(p)

    invisible(list(data=merged, correlation=r))
}

plot_interaction_sanity(
    interaction_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog/results/interaction/astrocyte_0__ic__vs__astrocyte_1__CaH__age_DE_results.txt",
    subset_one_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog/results/subset/astrocyte_0__ic__age_DE_results.txt",
    subset_two_path="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/emuratog/results/subset/astrocyte_1__CaH__age_DE_results.txt"

)
