
# aim ----
    # map the MP to UMAP 4066


# define environment ----
library(Seurat)

library(RColorBrewer)
library(clusterProfiler)

library(AUCell)
library(tidyverse)
library(readxl)









# arguments ----
obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4411_cellbender_filtered_integration_final_annotation.rds'
cell_ranking_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4411_AUCell_cellrankings.rds'

mp_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/meta_programs.rds'
# hm_gs_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/pathway_gmt/h.all.v2022.1.Hs.symbols.gmt'
gavish_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/meta_programs.xlsx'


out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis/map_GM_4411_seurat'
dir.create(out_dir, recursive = T)











########################################################################################################################
########################################################################################################################
########################################################################################################################
# load the seurat object 
obj <- readRDS(file = obj_file)






# load MP genes 
mp <- readRDS(file = mp_file)

# rename MPs so that they are the same in the heatmap
colnames(mp) <- c('C6', 'C16', 'C11', 'C12', 'C2', 'C3', 'C4', 'C8', 'C9', 'C10', 'C14', 'C1', 'C5', 'C7', 'C13', 'C15')


gene_list <- list()
for (cl in colnames(mp)) {
    gene_list[[cl]] <- mp[, cl]
}


# setdiff(mp[, 'C11'], mp[, 'C10'])

# setdiff(mp[, 'C10'], mp[, 'C11'])



gavish_mp <- read_excel(gavish_file)

gene_list[['Cell Cycle - G2_M']] <- gavish_mp[['Cell Cycle - G2/M']]
gene_list[['Cell Cycle - G1_S']] <- gavish_mp[['Cell Cycle - G1/S']]
gene_list[['Hypoxia']] <- gavish_mp[['Hypoxia']]






# score using seurat method 
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj)

for (gm in names(gene_list)) {
    working_obj <- AddModuleScore(obj, features = gene_list[gm], nbin = 30, ctrl = 50)
    obj@meta.data[paste0(gm, '_module_score')] <- working_obj$Cluster1
}   








# plot AUCell scores on UMAP 
umap_dir <- paste0(out_dir, '/', 'Gene_module_socre_on_UMAP')
dir.create(umap_dir, recursive = T)

for (gs in names(gene_list)) {
    
    pdf(file = paste0(umap_dir, '/', gs, '_scaled.pdf'))
        print(
            FeaturePlot(obj, features = paste0(gs, '_module_score'),
                order = TRUE, slot = 'scale.data') +
                labs(title = gs) +
                scale_color_gradient2(low = "#1b1b1b",
                                    mid = "#e6e6e6",
                                    high = "red") +
                theme(text = element_text(size = 20))
        )
    dev.off()
}








