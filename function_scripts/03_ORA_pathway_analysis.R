
# aim ----
    # script to do pathway analysis using DE results from single-cell data




# default 
    # use Padj < 0.05 as significance threshold



# arguments
    # de_res                    # default output from seurat analysis, 
                                    # with columns "gene", "p_val_adj", "avg_log2FC", "cluster"
    # entrez_conversion         # if input DE results shoudl be converted to entrez ID
    # FC_threshold              # FC threshold to select DE genes
    # gmt_file                  # full path to gmt file







# example
# path_res <- cluster_profiler_ora(de_res,
#                                 entrez_conversion = TRUE,
#                                 FC_threshold = 0.25,
#                                 gmt_file = './path/file')




cluster_profiler_ora <- function(
    de_res,
    entrez_conversion = TRUE,
    FC_threshold = 0.25,
    n_top_genes = NULL,
    gmt_file
) {
    # load required packages
    library(clusterProfiler)
    

    # convert data to entrez ID format if needed
    if (entrez_conversion == TRUE) {
        library(org.Hs.eg.db)
        entrez_ID_df <- bitr(de_res$gene, fromType="SYMBOL", toType="ENTREZID",
                     OrgDb="org.Hs.eg.db")
        de_res_entrez <- base::merge(de_res, entrez_ID_df,
                             by.x = 'gene', by.y = 'SYMBOL')           
    } else {
        de_res_entrez <- de_res
    }

    # filter DEGs by padj and FC
    if (!is.null(n_top_genes)) {
        de_res_entrez_filtered <- de_res_entrez %>%
            group_by(cluster) %>% 
            filter(p_val_adj < 0.05 & avg_log2FC > FC_threshold) %>% 
            slice_max(order_by = avg_log2FC, n = n_top_genes) %>% 
            ungroup()
    } else {
        de_res_entrez_filtered <- de_res_entrez %>%
            group_by(cluster) %>% 
            filter(p_val_adj < 0.05 & avg_log2FC > FC_threshold) %>% 
            ungroup()
    }
    

    # convert gene dataframe to gene lists for each cluster
    gene_list <- list()
    for (cl in unique(de_res_entrez_filtered$cluster)) {
        working_genes <- de_res_entrez_filtered[de_res_entrez_filtered$cluster == cl, 'ENTREZID'] %>% pull()
        output_name <- paste0('cluster_', cl)
        gene_list[[output_name]] <- working_genes
    }

    # load gmt file
    working_gmt <- read.gmt(gmt_file)

    # run pathway analysis 
    path_res <- compareCluster(geneCluster = gene_list,
                                fun = enricher,
                                TERM2GENE = working_gmt)

    return(path_res)
}



