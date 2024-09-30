
# aim ----
    # run ORA analysis for all modules derived




# load packages ----
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)




# arguments ----
gm_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/meta_programs.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis/pathway_analysis'
dir.create(out_dir, recursive = T)

reactome_gs <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/pathway_gmt/c2.cp.reactome.v2022.1.Hs.entrez.gmt'
gobp_gs <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/pathway_gmt/c5.go.bp.v7.4.entrez.gmt'






# load data ----
gm <- readRDS(file = gm_file)

# rename the gene modules to match the paper plot
colnames(gm) <- c('C6', 'C15', 'C10', 'C11', 'C2', 'C3', 'C4', 'C7', 'C8', 'C9', 'C13',
                    'C1', 'C5', 'exclude', 'C12', 'C14')


gene_list <- list()
for (cl in colnames(gm)[colnames(gm) != 'exclude']) {
    entrez_ID_df <- bitr(gm[,cl], fromType="SYMBOL", toType="ENTREZID",
                     OrgDb="org.Hs.eg.db")

    gene_list[[cl]] <- entrez_ID_df[['ENTREZID']]
}






# run ORA analysis
working_gmt <- read.gmt(gobp_gs)

res <- compareCluster(geneCluster = gene_list, fun = enricher, TERM2GENE = working_gmt, pvalueCutoff=0.05)

saveRDS(file = paste0(out_dir, '/', 'ORA_GOBP_analysis_results.rds'), res)








# plot top pathways for each 
top_n <- 5

pdf(file = paste0(out_dir, '/', 'GOBP_top_', top_n, '_pathway_dotplot.pdf'), width = 15, height = 20) 
    dotplot(res, showCategory = top_n)
                
dev.off()


scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
                    axis.text.y = element_text(size = 12))










