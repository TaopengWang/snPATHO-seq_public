
# aim ----
    # merge the PBMC data together and then do dimension reduction
    # don't do data integration




# define environment ----
library(Seurat)
library(tidyverse)
library(RColorBrewer)




# arguments ----
data_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/external_data/FHIL_Data/Aligned_30K_downsampled'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/snRNA_external_data/FHIL_PBMC/no_integration'
dir.create(out_dir, recursive = T)

out_data_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/external_data/FHIL_Data/processed_data'

working_dir <- out_dir

flex_probe_ref_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv'






# load data ----
obj_list <- list()

mtx <- Read10X(paste0(data_dir, '/FRP_PBMC1/outs/per_sample_outs/FRP_PBMC1/count/sample_filtered_feature_bc_matrix'))
obj <- CreateSeuratObject(counts = mtx)
obj$sample_id <- 'FLEX_exp_1'
obj_list[['FLEX_exp_1']] <- obj


mtx <- Read10X(paste0(data_dir, '/FRP_PBMC2/outs/per_sample_outs/FRP_PBMC2/count/sample_filtered_feature_bc_matrix'))
obj <- CreateSeuratObject(counts = mtx)
obj$sample_id <- 'FLEX_exp_2'
obj_list[['FLEX_exp_2']] <- obj


mtx <- Read10X(paste0(data_dir, '/V31_PBMC1/outs/filtered_feature_bc_matrix'))
obj <- CreateSeuratObject(counts = mtx)
obj$sample_id <- '3p_exp_1'
obj_list[['3p_exp_1']] <- obj


mtx <- Read10X(paste0(data_dir, '/V31_PBMC2/outs/filtered_feature_bc_matrix'))
obj <- CreateSeuratObject(counts = mtx)
obj$sample_id <- '3p_exp_2'
obj_list[['3p_exp_2']] <- obj







#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# filtering
filtered_obj_list <- lapply(obj_list, function(x){
    x <- subset(x, subset = nFeature_RNA > 200)
    x <- subset(x, subset = nFeature_RNA < 8000)
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- subset(x, subset = percent.mt < 10)
})








########################################################################################################################
########################################################################################################################
########################################################################################################################
# identify and remove doublet from the dataset ----
db_removal_dir <- paste0(working_dir, '/', 'doublet_removal')
dir.create(db_removal_dir, recursive = T)

source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/05_doublet_finder.R')

# identify doublet
doublet_identified_obj_list <- lapply(filtered_obj_list, function(x){
    output_label <- unique(x$sample_id)
    x <- DB_finder_removal(x,
                        expected_doublet_rate = 0.08,
                        out_dir = db_removal_dir,
                        output_name = output_label,
                        output_col_name = 'doublet_status')
    return(x)
})


# remove the cells
doublet_filtered_obj_list <- lapply(doublet_identified_obj_list, function(x){
    Idents(x) <- 'doublet_status'
    x <- subset(x, idents = 'Singlet')
})








########################################################################################################################
########################################################################################################################
########################################################################################################################
# merge data objects and dimension reduction 

merged_obj <- merge(doublet_filtered_obj_list[[1]], doublet_filtered_obj_list[2:length(doublet_filtered_obj_list)])

merged_obj$processing_method <- 
    ifelse(merged_obj$sample_id %in% c('FLEX_exp_1', 'FLEX_exp_2'), 'Frozen-Flex',
    ifelse(merged_obj$sample_id %in% c('3p_exp_1', '3p_exp_2'), "Frozen-3'", 'unlabelled'))

merged_obj$sample_id <- 
    ifelse(merged_obj$sample_id == 'FLEX_exp_1', "PBMC_replicate_1_Frozen-Flex",
    ifelse(merged_obj$sample_id == 'FLEX_exp_2', "PBMC_replicate_2_Frozen-Flex",
    ifelse(merged_obj$sample_id == '3p_exp_1', "PBMC_replicate_1_Frozen-3'",
    ifelse(merged_obj$sample_id == '3p_exp_2', "PBMC_replicate_2_Frozen-3'", 'unlabelled'))))


# clean up annotations 
for (cl in grep('pANN*', colnames(merged_obj@meta.data), value = T)) {
    merged_obj@meta.data[[cl]] <- NULL
}

for (cl in grep('DF*', colnames(merged_obj@meta.data), value = T)) {
    merged_obj@meta.data[[cl]] <- NULL
}

# merged_obj$doublet_status <- NULL


# dimension reduction
merged_obj <- merged_obj %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:30) %>% 
    RunUMAP(dims = 1:30) 





# load predefined colors
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')
workflow_colors <- workflow_colors[order(names(workflow_colors))]



# plot check batch factors and QC values
pdf(paste0(working_dir, '/', 'umap_batch_factors.pdf'))
    print(
        DimPlot(merged_obj, reduction = "umap", group.by = 'processing_method',
            cols = workflow_colors) +
            theme(text = element_text(size = 30))
    )
dev.off()







# split the plot
pdf(file = paste0(working_dir, '/', 'umap_split_by_processing_method.pdf'), width = 10)
    print(
        DimPlot(merged_obj, reduction = "umap", group.by = 'processing_method', 
            split.by = 'processing_method', cols = workflow_colors) +
            theme(text = element_text(size = 30))
    )
dev.off()



pdf(file = paste0(working_dir, '/', 'umap_split_by_sample_id.pdf'), width = 20)
    print(
        DimPlot(merged_obj, reduction = "umap", group.by = 'processing_method', 
            split.by = 'sample_id', cols = workflow_colors) +
                 theme(text = element_text(size = 20))
    )
dev.off()












########################################################################################################################
########################################################################################################################
########################################################################################################################
# explore PCA plots 

pdf(file = paste0(working_dir, '/', 'pca_batch_factors.pdf'))
    print(
        DimPlot(merged_obj, reduction = 'pca', group.by = 'processing_method', 
            cols = workflow_colors) +
            theme(text = element_text(size = 30)) +
            labs(title = 'Processing method')
    )
dev.off()


pdf(file = paste0(working_dir, '/', 'pca_split_by_processing_methods.pdf'), width = 14)
    print(
        DimPlot(merged_obj, reduction = 'pca', group.by = 'processing_method',
            split.by = 'processing_method', cols = workflow_colors) 
    )
dev.off()


pdf(file = paste0(working_dir, '/', 'pca_split_by_sample_id.pdf'), width = 16)
    print(
        DimPlot(merged_obj, reduction = "pca", group.by = 'processing_method', 
            split.by = 'sample_id', cols = workflow_colors) +
            theme(text = element_text(size = 18)) +
            labs(title = NULL)
    )
dev.off()




# plot drivers for top PCs
pdf(file = paste0(working_dir, '/', 'top_pca_drivers_heatmap.pdf'), width = 14, height = 14)
    print(
        DimHeatmap(merged_obj, dims = 1:9, cells = 500, balanced = TRUE)
    )
dev.off()




# plot a few markers on PCA 
markers <- c('CD3D', 'CD14', 'MS4A1')
for (m in markers) {
    pdf(file = paste0(working_dir, '/', 'PCA_', m, '.pdf'))
        print(
            FeaturePlot(merged_obj, reduction = 'pca', features = m, order = T) +
                theme(text = element_text(size = 30))
        )
    dev.off()
}







# plot pca loadings and color genes by panel 
flex_probe_ref <- read.csv(file = flex_probe_ref_file, skip = 5)
flex_probe_ref <- flex_probe_ref[flex_probe_ref$included == 'TRUE', ]
flex_genes <- unique(unlist(lapply(str_split(flex_probe_ref$probe_id, '\\|'), '[[', 2)))

for (dim in c(1, 2, 3, 4, 5)) {
    # get top 30 genes by pca loadings 
    working_pc_name <- paste0('PC_', dim)
    loadings <- merged_obj@reductions$pca@feature.loadings[, working_pc_name]
    loadings <- loadings[order(abs(loadings), decreasing = T)]
    plot_genes <- names(loadings[1:30])
    plot_loadings <- loadings[plot_genes]
    plot_loadings <- plot_loadings[order(plot_loadings, decreasing = T)]
    reordered_plot_genes <- names(plot_loadings)

    plot_colors <- ifelse(reordered_plot_genes %in% flex_genes, 'black', 'red')
    names(plot_colors) <- reordered_plot_genes

    pdf(file = paste0(working_dir, '/', 'top_pca_loadings_', working_pc_name, '.pdf'), width = 4, height = 7)
        print(
            VizDimLoadings(merged_obj, dims = dim, reduction = "pca", ncol = 1) +
                theme(axis.text.y = element_text(color = rev(plot_colors)),
                    text = element_text(size = 24))
        )
    dev.off()

}








########################################################################################################################
########################################################################################################################
########################################################################################################################
# plot the expression of common cell type markers 
marker_dir <- paste0(working_dir, '/', 'markers')
dir.create(marker_dir, recursive = T)

markers <- c('CD3D', 'CD14', 'CD4', 'CD8A', 'MS4A1', 'CD79A')

for (g in markers) {
    png(file = paste0(marker_dir, '/', g, '.png'), units = 'in', width = 7, height = 7, res = 300)
        print(
            FeaturePlot(merged_obj, features = g, order = T, raster = F) +
                theme(text = element_text(size = 30))
        )
    dev.off()
}







########################################################################################################################
########################################################################################################################
########################################################################################################################
# save processed data 
saveRDS(file = paste0(out_data_dir, '/', 'PBMC_no_integration_obj.rds'), merged_obj)

