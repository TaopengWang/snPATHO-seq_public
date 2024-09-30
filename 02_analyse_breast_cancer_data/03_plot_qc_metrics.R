
# aim ----
    # simply plot QC metrics of all data 





# define environment ----
library(here)
setwd(here())

# make sure renv is loaded
source(paste0(here(), '/', '.Rprofile'))
print(.libPaths())

# load other packages
library(Seurat)
library(tidyverse)
library(ggpubr)









# arguments ----
data_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/integration/25k_downsampled/cellbender/all_data_QC'
dir.create(out_dir, recursive = T)





# load data ----
obj_4066 <- readRDS(file = paste0(data_dir, '/', '4066_cellbender_filtered_integration.rds'))
obj_4399 <- readRDS(file = paste0(data_dir, '/', '4399_cellbender_filtered_integration.rds'))
obj_4411 <- readRDS(file = paste0(data_dir, '/', '4411_cellbender_filtered_integration.rds'))




# merge all data 
merged_obj <- merge(obj_4066, list(obj_4399, obj_4411))

merged_obj$patient_id <- ifelse(merged_obj$sample_id %in% c("4066FFPE", "4066SNAPFix", "4066GEX"), '4066',
    ifelse(merged_obj$sample_id %in% c("4411FFPE", "4411SNAPFix", "4411GEX"), '4411', '4399'))




# plot numbers of UMIs per nucleus
pdf(file = paste0(out_dir, '/', 'nCounts_all_data.pdf'), width = 10, height = 4)
    print(VlnPlot(merged_obj, features = 'nCount_RNA', 
        group.by = 'sample_id', split.by = 'patient_id',
        pt.size = 0, log = TRUE) +
        labs(x = "Samples",
            y = 'Numbers of counts per nuclei',
            title = NULL))
dev.off()







# plot numbers of genes per nucleus
pdf(file = paste0(out_dir, '/', 'nFeature_all_data.pdf'), width = 10, height = 4)
    print(VlnPlot(merged_obj, features = 'nFeature_RNA', 
        group.by = 'sample_id', split.by = 'patient_id',
        pt.size = 0) +
        labs(x = "Samples",
            y = 'Genes detected per nuclei',
            title = NULL))
dev.off()




