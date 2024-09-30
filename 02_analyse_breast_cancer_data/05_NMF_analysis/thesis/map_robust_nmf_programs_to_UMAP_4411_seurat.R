
# aim ----
    # get Seurat gene module scores of each robust NMF program on each sample




# define environment ----
library(Seurat)
library(tidyverse)
library(RColorBrewer)


# arguments ----
obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4411_cellbender_filtered_integration_final_annotation.rds'

robust_nmf_program_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_programs.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis/map_robust_nmf_programs_4411_seurat'
dir.create(out_dir, recursive = T)









########################################################################################################################
########################################################################################################################
########################################################################################################################
# load the seurat object 
obj <- readRDS(file = obj_file)






# load MP genes 
robust_nmf_program <- readRDS(file = robust_nmf_program_file)

gene_list <- list()
for (cl in colnames(robust_nmf_program)) {
    gene_list[[cl]] <- robust_nmf_program[, cl]
}






# score using seurat method 
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj)

for (gm in names(gene_list)) {
    working_obj <- AddModuleScore(obj, features = gene_list[gm], nbin = 30, ctrl = 50)
    obj@meta.data[paste0(gm, '_module_score')] <- working_obj$Cluster1
}   








# plot AUCell scores on UMAP 
for (gs in names(gene_list)) {
    
    pdf(file = paste0(out_dir, '/', gs, '_scaled.pdf'))
        print(
            FeaturePlot(obj, features = paste0(gs, '_module_score'),
                order = TRUE, slot = 'scale.data') +
                labs(title = NULL) +
                scale_color_gradient2(low = "#1b1b1b",
                                    mid = "#e6e6e6",
                                    high = "red") +
                theme(text = element_text(size = 20))
        )
    dev.off()
}





