
# aim ----
    # perform a simple dimension reduction





# load packages ----
library(Seurat)
library(tidyverse)







# arguments ----
data_obj_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/filtered/seurat_obj'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/no_integration/25k_downsampled/cellbender/all_samples'
dir.create(out_dir, recursive = T)







########################################################################################################################
########################################################################################################################
########################################################################################################################
# load data ----
obj_list <- list()
for (f in dir(data_obj_dir)) {
    obj_list[[f]] <- readRDS(file = paste0(data_obj_dir, '/', f))
}


# basic filtering
filtered_obj_list <- lapply(obj_list, function(x){
    x <- subset(x, subset = nFeature_RNA > 200)
    x <- subset(x, subset = nFeature_RNA < 8000)
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- subset(x, subset = percent.mt < 10)
})










# merge data 
merged_data <- merge(filtered_obj_list[[1]], y = filtered_obj_list[2: length(filtered_obj_list)])




# dimension reduction
merged_data <-  merged_data %>% 
                NormalizeData() %>% 
                FindVariableFeatures() %>% 
                ScaleData() %>% 
                RunPCA() %>% 
                FindNeighbors(dims = 1:30) %>% 
                RunUMAP(dims = 1:30)




# modify annotation 
merged_data$workflow <- ifelse(merged_data$sample_id %in% c('4066FFPE', '4399FFPE', '4411FFPE'), 'FFPE-snPATHO-Seq',
                        ifelse(merged_data$sample_id %in% c('4066GEX', '4399GEX', '4411GEX'), "Frozen-3'", 'Frozen-Flex'))

merged_data$patient <- ifelse(merged_data$sample_id %in% c("4066FFPE", "4066GEX", "4066SNAPFix"), '4066',
                        ifelse(merged_data$sample_id %in% c("4399FFPE", "4399GEX", "4399SNAPFix"), "4399", '4411'))




source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')


# plot to check results 
pdf(file = paste0(out_dir, '/', 'UMAP_workflow.pdf'))
    print(
        DimPlot(merged_data, group.by = 'workflow', cols = workflow_colors)
    )
dev.off()

pdf(file = paste0(out_dir, '/', 'UMAP_workflow.pdf'), width = 21)
    print(
        DimPlot(merged_data, group.by = 'workflow', split.by = 'workflow', cols = workflow_colors, ncol = 3) + 
            theme(text = element_text(size = 24))
    )
dev.off()

pdf(file = paste0(out_dir, '/', 'UMAP_samples.pdf'))
    print(
        DimPlot(merged_data, group.by = 'patient') + 
            theme(text = element_text(size = 24))
    )
dev.off()

