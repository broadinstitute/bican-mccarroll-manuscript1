inDegTable="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/original_data//DEGs-Table 1.tsv"
inDegCellKey="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/original_data/key for cell types-Table 1.tsv"
ontologyMapFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/original_data/cell_ontology_map_file.txt"
outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_40903571/voom-like"

a=read.table(inDegCellKey, fill=T, sep="\t")
b=read.table(inDegTable, sep="\t")
colnames(b)[1:7]=c("gene", "logFC", "ratio_elderly", "ratio_adult", "wilcox_p", "q-value", "cell_type")
b$cell_type_full=a[match(b$cell_type, a$V1),]$V2

ontologymap=read.table(ontologyMapFile, header=T, stringsAsFactors = F, sep="\t")

idx=match(b$cell_type_full, ontologymap$paper)
b$bican_cell_type=ontologymap[idx,]$bican
b=b[!is.na(b$bican_cell_type),]

cell_types=unique (b$bican_cell_type)

# convert each cell type into a limma-voom like file
#ct=cell_types[1]
for (ct in cell_types) {
    # logFC	AveExpr	t	P.Value	adj.P.Val	B	z.std
    z=b[b$bican_cell_type==ct,]
    df=data.frame(logFC=z$logFC, AveExpr=NA, t=NA, P.Value=z$wilcox_p, adj.P.Val=z$`q-value`, B=NA, z.std=NA)
    rownames(df)=z$gene
    outFile=paste(outDir, "/", ct, ".age_DE_results.txt", sep="")
    write.table(df, outFile, row.names=T, col.names=T, quote=F, sep="\t")
}
