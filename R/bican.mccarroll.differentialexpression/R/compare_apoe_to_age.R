# A little ad-hoc analysis of APOE status DE vs age DE

# age_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old//LEVEL_2_apoe_score/sex_age/cell_type"
# apoe_dir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3_analysis/differential_expression/results/old/LEVEL_2/sex_age/cell_type"
# cell_type="microglia"

compare_apoe_to_age_de<-function (age_dir, apoe_dir, cell_type) {

    age_file=file.path(age_dir, paste0(cell_type, "_age_DE_results.txt"))
    apoe_file=file.path(apoe_dir, paste0(cell_type, "_apoe_score_DE_results.txt"))
    if (!file.exists(age_file)) {
        stop(paste("Age DE results file does not exist:", age_file))
    }
    if (!file.exists(apoe_file)) {
        stop(paste("APOE DE results file does not exist:", apoe_file))
    }

    age=read.table(age_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    apoe=read.table(apoe_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)

    both=intersect (rownames (age), rownames (apoe))
    age=age[both, ]
    apoe=apoe[both, ]

    #subset to results significant in either
    sig_age=which(age$adj.P.Val<0.05 | apoe$adj.P.Val <0.05)
    age=age[sig_age, ]
    apoe=apoe[sig_age, ]

    plot (age$logFC, apoe$logFC, xlab="Age logFC", ylab="APOE logFC", main=paste0(" Age vs APOE DE in ", cell_type))

}
