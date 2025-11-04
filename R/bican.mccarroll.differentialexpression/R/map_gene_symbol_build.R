# library (data.table)
# library(rtracklayer)

#gene_symbols_source=model$feature
# source_gtf_file="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.reduced.gtf"
# target_gtf_file="/broad/mccarroll/software/metadata/individual_reference/GRCh38_maskedAlt.89/GRCh38_maskedAlt.reduced.gtf"

gene_symbol_mapper<-function (gene_symbols_source, source_gtf_file, target_gtf_file) {
    source <- fread(source_gtf_file, header = T, sep = "\t", data.table = FALSE)
    target<- fread(target_gtf_file, header = T, sep = "\t", data.table = FALSE)

    source_genes <- source[source$annotationType == "gene", ]
    target_genes <- target[target$annotationType == "gene", ]

    #map from the source gene symbol to the ENSG ID
    idx=match(gene_symbols_source, source_genes$gene_name)
    len_missing=length(which(is.na(idx)))
    if (len_missing>0){
        warning(paste0(len_missing, " gene symbols could not be mapped to ENSG IDs in the source GTF file. They will be returned as NA."))
    }

    source_ensg_ids <- source_genes$gene_id[idx]

    df=data.frame(source_gene_symbol=gene_symbols_source, source_ensg_id=source_ensg_ids)

    #map from the ENSG ID to the target gene symbol
    idx2=match(df$source_ensg_id, target_genes$gene_id)
    target_gene_symbols <- target_genes$gene_name[idx2]

    df$target_gene_symbol=target_gene_symbols
    len_missing2=length(which(is.na(idx2)))
    if (len_missing2>0){
        warning(paste0(len_missing2, " ENSG IDs could not be mapped to gene symbols in the target GTF file. They will be returned as NA."))
    }
    return (df)

}



# returns data.frame with ENSG (no version), gene_type, gene_name
#gtf_path="/broad/mccarroll/software/metadata/individual_reference/GRCh38_ensembl_v43/GRCh38_ensembl_v43.gtf"
parse_gtf_genes_df <- function(gtf_path) {
    gr <- BiocIO::import(gtf_path)
    df <- as.data.frame(gr)

    # keep only gene-level rows
    df <- df[df$type == "gene", ]

    # normalize ENSG IDs (drop version suffix)
    df$ENSG <- sub("\\.\\d+$", "", as.character(df$gene_id))

    # handle type naming differences
    if ("gene_type" %in% names(df)) {
        df$gene_type <- as.character(df$gene_type)
    } else if ("gene_biotype" %in% names(df)) {
        df$gene_type <- as.character(df$gene_biotype)
    } else {
        df$gene_type <- NA_character_
    }

    df$gene_name <- as.character(df$gene_name)

    df=df[, c("seqnames", "gene_type", "gene_name", "ENSG")]
    return (df)
}

