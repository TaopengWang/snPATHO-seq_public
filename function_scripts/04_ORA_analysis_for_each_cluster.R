
# aim ----
    # run simple ORA analysis on the provided gene list 







# arguments
    # gene_list                 # named gene list, name are gene names or entrez ID, filtered already
    # entrez_conversion         # if names should be converted from HUGO ID to entrez ID
    # gmt_file                  # pathway file to be used
    # method                    # 'ORA' or 'GSEA' analysis



clusterprofiler_ora_single_list <- function(
    gene_list,
    entrez_conversion = TRUE,
    gmt_file,
    method = 'ORA'
) {
    # load required packages
    library(clusterProfiler)
    
    # convert data to entrez ID format if needed
    if (entrez_conversion == TRUE) {
        library(org.Hs.eg.db)
        entrez_ID_df <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID",
                     OrgDb="org.Hs.eg.db")
        # sometimes not all genes can be mapped properly
        gene_df <- as.data.frame(gene_list)
        entrez_ID_df <- base::merge(entrez_ID_df, gene_df, by.x = 'SYMBOL', by.y = 0, all.x = TRUE)
        gene_list_entrez <- entrez_ID_df$gene_list
        names(gene_list_entrez) <- entrez_ID_df$ENTREZID
    } else {
        gene_list_entrez <- gene_list
    }

    # order gene list from high to low
    gene_list_entrez <- gene_list_entrez[order(gene_list_entrez, decreasing = TRUE)]

    # load gmt file
    working_gmt <- read.gmt(gmt_file)

    # run the analysis
    if (method == 'ORA') {
        path_res <- enricher(names(gene_list_entrez), TERM2GENE = working_gmt)
    } else if (method == 'GSEA') {
        path_res <- GSEA(gene_list_entrez, 
            TERM2GENE = working_gmt,
            minGSSize = 15,
            maxGSSize = 500)
    }

    # ouptut
    return(path_res)
}











