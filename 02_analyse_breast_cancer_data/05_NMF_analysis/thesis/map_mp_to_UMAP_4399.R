
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
obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4399_cellbender_filtered_integration_final_annotation.rds'
cell_ranking_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4399_AUCell_cellrankings.rds'

mp_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/meta_programs.rds'
# hm_gs_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/pathway_gmt/h.all.v2022.1.Hs.symbols.gmt'
gavish_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/meta_programs.xlsx'
nmf_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_programs.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis/map_GM_4399'
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



gavish_mp <- read_excel(gavish_file)

gene_list[['Cell Cycle - G2_M']] <- gavish_mp[['Cell Cycle - G2/M']]
gene_list[['Cell Cycle - G1_S']] <- gavish_mp[['Cell Cycle - G1/S']]
gene_list[['Hypoxia']] <- gavish_mp[['Hypoxia']]




nmf_programs <- readRDS(file = nmf_file)
gene_list[["4399_3'_9.8"]] <- nmf_programs[,'4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.9.8']
gene_list[["4399_3'_7.6"]] <- nmf_programs[, '4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.6']



# load defined ranking file 
cells_rankings <- readRDS(file = cell_ranking_file)








# calculate the AUCell scores
cells_AUC <- AUCell_calcAUC(geneSets = gene_list, cells_rankings, nCores = 4)
saveRDS(file = paste0(out_dir, '/', 'AUCell_scores.rds'), cells_AUC)







aucell_assay <- CreateAssayObject(data = cells_AUC@assays@data@listData$AUC)
obj[['AUCell']] <- aucell_assay
DefaultAssay(obj) <- 'AUCell'
obj <- ScaleData(obj)






# plot AUCell scores on UMAP 
umap_dir <- paste0(out_dir, '/', 'AUCell_on_UMAP')
dir.create(umap_dir, recursive = T)

for (gs in rownames(obj[['AUCell']])) {
    
    pdf(file = paste0(umap_dir, '/', gs, '_scaled.pdf'))
        print(
            FeaturePlot(obj, features = gs,
                order = TRUE, slot = 'scale.data') +
                labs(title = gs) +
                scale_color_gradient2(low = "blue",
                                    mid = "white",
                                    high = "red") +
                theme(text = element_text(size = 20))
        )
    dev.off()

    pdf(file = paste0(umap_dir, '/', gs, '.pdf'))
        print(
            FeaturePlot(obj, features = gs, cols = c('lightgrey', 'red'),
                order = TRUE) +
                labs(title = gs)+
                theme(text = element_text(size = 20))
        )
    dev.off()
}








