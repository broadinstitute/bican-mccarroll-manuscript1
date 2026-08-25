inDegTableS6="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_39227716/original_data/Table S6-Table 1.tsv"
inDegTableS9="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_39227716/original_data/Table S9-Table 1.tsv"

ontologyMapFile="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_39227716/original_data/cell_ontology_map_file.txt"
outDir="/broad/bican_um1_mccarroll/RNAseq/analysis/CAP_freeze_3.1_analysis/differential_expression/external_comparison_PMID_39227716/voom-like"

#skip the first 3 lines, which are a title/header in the sheet.
b1=read.table(inDegTableS6, sep="\t", header=T, skip=3)
b2=read.table(inDegTableS9, sep="\t", header=T, skip=3)

#before merging, remove assays from b1 that are in b2.
#this is because the FDR is somehow slightly different.
b1=b1[!b1$assay %in% b2$assay,]
b=rbind(b1, b2)

ontologymap=read.table(ontologyMapFile, header=T, stringsAsFactors = F, sep="\t")
all (ontologymap$paper %in% b$assay)

b=b[b$assay %in% ontologymap$paper,]

idx=match(b$assay, ontologymap$paper)
b$bican_cell_type=ontologymap[idx,]$bican
b=b[!is.na(b$bican_cell_type),]

table (b$assay, b$bican_cell_type)

cell_types=unique (b$bican_cell_type)
region="OFC"

# convert each cell type into a limma-voom like file
#ct=cell_types[1]
for (ct in cell_types) {
    # logFC	AveExpr	t	P.Value	adj.P.Val	B	z.std
    z=b[b$bican_cell_type==ct,]
    z=unique (z)
    if (length(unique(z$ID))!=dim(z)[1]) stop ("Mismatch between number of rows and number of unique gene IDs")
    #Note, we're going to scale the logFC to be in decades instead of years so their results match our scaling.
    df=data.frame(logFC=10*z$logFC, AveExpr=z$AveExpr, t=z$t, P.Value=z$P.Value, adj.P.Val=z$adj.P.Val, B=z$B, z.std=z$z.std)
    rownames(df)=z$ID
    #outFile=paste(outDir, "/", ct, ".age_DE_results.txt", sep="")
    outFile <- file.path( outDir, paste0(ct, "__", region, "__age_DE_results.txt"))
    write.table(df, outFile, row.names=T, col.names=T, quote=F, sep="\t")
}
