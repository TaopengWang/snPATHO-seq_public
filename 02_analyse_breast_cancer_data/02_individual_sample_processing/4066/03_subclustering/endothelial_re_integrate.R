

# aim ----
    # reintegrate and cluster cancer epithelial cells in the data







# define environment ----
library(here)
setwd(here())

# make sure renv is loaded
source(paste0(here(), '/', '.Rprofile'))
print(.libPaths())

# load other packages
library(Seurat)
library(tidyverse)
library(RColorBrewer)

library(ComplexHeatmap)
library(circlize)




# arguments ----
obj_file_path <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4066_cellbender_filtered_integration_initial_annotation.rds'

sample_id <- '4066'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/integration/25k_downsampled/cellbender'



















####################################################################################################
####################################################################################################
####################################################################################################
# load data ----
obj <- readRDS(file = obj_file_path)

subset_obj <- subset(obj, subset = initial_annotation %in% c('Endohelial'))



# split data by sample id
obj_list <- list()
for (ds in unique(subset_obj$sample_id)) {
    obj_list[[ds]] <- subset(subset_obj, subset = sample_id == ds)
}
















####################################################################################################
####################################################################################################
####################################################################################################
# redo integration
# load function
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/01_Seurat_integration.R')


# integrate data 
integrated_obj <- seurat_integration(obj_list,
                                    normalisation_method = 'RNA',
                                    nfeatue_each_dataset = 500,
                                    nfeatures_common = 500,
                                    vars_to_regress = NULL,
                                    integration_method = 'cca',
                                    nthread = NULL)
















####################################################################################################
####################################################################################################
####################################################################################################
# check integration quality and annotate clusters
working_dir <- paste0(out_dir, '/', sample_id, '/', 'subclustering', '/', 'Endothelial_reintegrated')
dir.create(working_dir, recursive = T)

# plot check batch factors and QC values
pdf(paste0(working_dir, '/', 'umap_batch_factors.pdf'))
    print(
        DimPlot(integrated_obj, reduction = "umap", group.by = 'sample_id')
    )
    print(
        FeaturePlot(integrated_obj, features = "nCount_RNA", reduction = 'umap') +
            scale_color_continuous(trans='log2', type = 'viridis')
    )
    print(
        FeaturePlot(integrated_obj, features = "nFeature_RNA", reduction = 'umap') +
            scale_color_continuous(trans='log2', type = 'viridis')
    )
    print(
        FeaturePlot(integrated_obj, features = "percent.mt", reduction = 'umap')
    )
dev.off()









# split umap by sample
pdf(file = paste0(working_dir, '/', 'umap_split_by_sample_id.pdf'),
    width = 5 * length(unique(integrated_obj$sample_id)),
    height = 5)
    print(
        DimPlot(integrated_obj, reduction = "umap", group.by = 'sample_id', 
            split.by = 'sample_id', ncol = 3)
    )
dev.off()





















####################################################################################################
####################################################################################################
####################################################################################################
# find optimal clustering resolution ----
clust_dir <- paste0(working_dir, '/', 'cluster_optimisation')
dir.create(clust_dir, recursive = T)

# load function
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/02_Seurat_clustering_optimisation.R')

# make sure the default assay is the integrated assay
DefaultAssay(integrated_obj) <- 'integrated'

# get distance matrix for silhouette score calcualtion
dist.matrix <- dist(x = Embeddings(integrated_obj, reduction = "pca"))

# cluster and find optimal clustering resolution
integrated_obj <- clustering_resolution_optimisation(integrated_obj,
                                        resolution_range = seq(0.2, 2, by = 0.1),
                                        plot_dir = clust_dir,
                                        batch_factor = "sample_id",
                                        run_silhouette = TRUE,
                                        dist.matrix = dist.matrix)






















####################################################################################################
####################################################################################################
####################################################################################################
# annotate subclusters
optimal_clustering_dir <- paste0(working_dir, '/', 'optimal_cluster')
dir.create(optimal_clustering_dir, recursive = T)


opt_res <- 0.4
cluster_name <- grep(paste0('res.', opt_res), colnames(integrated_obj@meta.data), value = T)[1]



# make a plot with label at the optimal cluster resolution ----
pdf(file = paste0(optimal_clustering_dir, '/', 'umap_', opt_res, '.pdf'))
    print(
        DimPlot(integrated_obj, reduction = "umap", group.by = cluster_name, label = T, repel = T)
    )
dev.off()





pdf(file = paste0(optimal_clustering_dir, '/', 'nCounts_violin.pdf'), width = 8, height = 4)
    print(
        VlnPlot(integrated_obj, features = 'nCount_RNA', group.by = cluster_name)
    )
dev.off()

pdf(file = paste0(optimal_clustering_dir, '/', 'nFeatures_violin.pdf'), width = 8, height = 4)
    print(
        VlnPlot(integrated_obj, features = 'nFeature_RNA', group.by = cluster_name)
    )
dev.off()

















####################################################################################################
####################################################################################################
####################################################################################################
# find top markers in one of the clusters
de_dir <- paste0(optimal_clustering_dir, '/', 'de_res')
dir.create(de_dir, recursive = T)

Idents(integrated_obj) <- cluster_name
de_res <- FindAllMarkers(integrated_obj, assay = 'RNA', test.use = 't', only.pos = TRUE)
write.csv(file = paste0(out_dir, '/', 'de_res.csv'), de_res)



























##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# check top pathways ----
reactome_gs <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/pathway_gmt/c2.cp.reactome.v2022.1.Hs.entrez.gmt'
output_label <- 'reactome_gs'


source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/03_ORA_pathway_analysis.R')

path_res <- cluster_profiler_ora(de_res, 
    entrez_conversion = TRUE,
    FC_threshold = 0.25,
    gmt_file = eval(sym(output_label)))
saveRDS(file = paste0(optimal_clustering_dir, '/', output_label, '_analysis_results.rds'), path_res)

plot <- dotplot(path_res, showCategory = 5) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

pdf(file = paste0(optimal_clustering_dir, '/', output_label, '_results.pdf'),
    width = 10, height = 16)
    print(plot)
dev.off()
















##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# plot heatmap to check top markers in each cluster ----
reactome_ecm_organisation_genes <- c('A2M','ACAN','ACTN1','ADAM10','ADAM12','ADAM15','ADAM17','ADAM19','ADAM8','ADAM9','ADAMTS1','ADAMTS14','ADAMTS16','ADAMTS18','ADAMTS2','ADAMTS3','ADAMTS4','ADAMTS5','ADAMTS8','ADAMTS9','AGRN','APP','ASPN','BCAN','BGN','BMP1','BMP10','BMP2','BMP4','BMP7','BSG','CAPN1','CAPN10','CAPN11','CAPN12','CAPN13','CAPN14','CAPN15','CAPN2','CAPN3','CAPN5','CAPN6','CAPN7','CAPN8','CAPN9','CAPNS1','CAPNS2','CASK','CASP3','CAST','CD151','CD44','CD47','CDH1','CEACAM1','CEACAM6','CEACAM8','CMA1','COL10A1','COL11A1','COL11A2','COL12A1','COL13A1','COL14A1','COL15A1','COL16A1','COL17A1','COL18A1','COL19A1','COL1A1','COL1A2','COL20A1','COL21A1','COL22A1','COL23A1','COL24A1','COL25A1','COL26A1','COL27A1','COL28A1','COL2A1','COL3A1','COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6','COL5A1','COL5A2','COL5A3','COL6A1','COL6A2','COL6A3','COL6A5','COL6A6','COL7A1','COL8A1','COL8A2','COL9A1','COL9A2','COL9A3','COLGALT1','COLGALT2','COMP','CRTAP','CTRB1','CTRB2','CTSB','CTSD','CTSG','CTSK','CTSL','CTSS','CTSV','DAG1','DCN','DDR1','DDR2','DMD','DMP1','DSPP','DST','EFEMP1','EFEMP2','ELANE','ELN','EMILIN1','EMILIN2','EMILIN3','F11R','FBLN1','FBLN2','FBLN5','FBN1','FBN2','FBN3','FGA','FGB','FGF2','FGG','FMOD','FN1','FURIN','GDF5','HAPLN1','HSPG2','HTRA1','IBSP','ICAM1','ICAM2','ICAM3','ICAM4','ICAM5','ITGA1','ITGA10','ITGA11','ITGA2','ITGA2B','ITGA3','ITGA4','ITGA5','ITGA6','ITGA7','ITGA8','ITGA9','ITGAD','ITGAE','ITGAL','ITGAM','ITGAV','ITGAX','ITGB1','ITGB2','ITGB3','ITGB4','ITGB5','ITGB6','ITGB7','ITGB8','JAM2','JAM3','KDR','KLK2','KLK7','KLKB1','LAMA1','LAMA2','LAMA3','LAMA4','LAMA5','LAMB1','LAMB2','LAMB3','LAMC1','LAMC2','LAMC3','LOX','LOXL1','LOXL2','LOXL3','LOXL4','LRP4','LTBP1','LTBP2','LTBP3','LTBP4','LUM','MADCAM1','MATN1','MATN3','MATN4','MFAP2','MFAP3','MFAP4','MFAP5','MMP1','MMP10','MMP11','MMP12','MMP13','MMP14','MMP15','MMP16','MMP17','MMP19','MMP2','MMP20','MMP24','MMP25','MMP3','MMP7','MMP8','MMP9','MUSK','NCAM1','NCAN','NCSTN','NID1','NID2','NRXN1','NTN4','OPTC','P3H1','P3H2','P3H3','P4HA1','P4HA2','P4HA3','P4HB','PCOLCE','PCOLCE2','PDGFA','PDGFB','PECAM1','PHYKPL','PLEC','PLG','PLOD1','PLOD2','PLOD3','PPIB','PRKCA','PRSS1','PRSS2','PSEN1','PTPRS','PXDN','SCUBE1','SCUBE3','SDC1','SDC2','SDC3','SDC4','SERPINE1','SERPINH1','SH3PXD2A','SPARC','SPOCK3','SPP1','TGFB1','TGFB2','TGFB3','THBS1','TIMP1','TIMP2','TLL1','TLL2','TMPRSS6','TNC','TNN','TNR','TNXB','TPSAB1','TRAPPC4','TTR','VCAM1','VCAN','VTN','VWF')




# complexheatmap
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/08_prepare_data_for_complex_heatmap.R')

grouping_annotations <- c(cluster_name)


plot_data_list <- prepare_data_for_heatmap(
    seurat_obj = integrated_obj,
    seurat_assay = 'RNA', 
    serat_slot = 'data',
    genes_to_plot = reactome_ecm_organisation_genes,
    grouping_annotations = grouping_annotations,
    additional_annotations = NULL)




# put the heatmap together ----
# scale the expression matrix
plot_data_list[[1]] <- plot_data_list[[1]][rowSums(plot_data_list[[1]]) > 0, ]
plot_data_list[[1]] <- scale(t(plot_data_list[[1]]))
plot_data_list[[1]] <- t(plot_data_list[[1]])

# set expression data colors 
col_fun <- colorRamp2(c(floor(min(plot_data_list[[1]])),
                        mean(plot_data_list[[1]]),
                        ceiling(max(plot_data_list[[1]]))), 
                      c("#2166AC", "white", "#B2182B"))









# set annotation colors & make annotation ----
cluster_colors <- set_colors(plot_data_list[[2]], cluster_name, palette = 'Accent')



# make annotation
right_annot <- HeatmapAnnotation(
    Clusters = plot_data_list[[2]][[cluster_name]],
    which = 'row',
    col = list(Clusters = cluster_colors), 
    annotation_name_gp = grid::gpar(fontsize = 12),
    annotation_legend_param = list(
        Clusters = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12))))







# plot heatmap ----
plot <- Heatmap(t(plot_data_list[[1]]),
                name = 'Scaled exp',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold")),
                col = col_fun,
                right_annotation = right_annot,
                show_column_names= TRUE,
                show_row_names = TRUE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                column_names_gp = gpar(fontsize = 6),
                row_names_gp = gpar(fontsize = 8))

pdf(file = paste0(optimal_clustering_dir, '/', 'reactome_ecm_organisation_genes_heatmap', '.pdf'),
    width = 18, height = 7)
    draw(plot)
dev.off()


















##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# plot expression of some of these ECM genes
marker_dir <- paste0(working_dir, '/', 'markers')
dir.create(marker_dir, recursive = T)

markers <- c('COL6A2', 'CALD1', 'IL7R', 'TRAC', 'TRBC2', 'KLRK1', 'CD96', 'CD2', 'CELF2', 'CXCR4',
    'B2M', 'FTL', 'MT−ND6', 'MGP', 'COL1A1', 'COL1A2', 'LUM', 'SPARC', 'COL3A1', 'ACTB', 'CD9', 'CLU',
    'LAMP1', 'MYADM', 'ATP6AP1', 'HSPA2', 'SPINT2', 'DHCR24', 'DSP', 'RGS1', 'ZEB2', 'SGK1', 'MS4A6A', 'SPI1', 
    'STAB1', 'ITGAX', 'CIITA', 'CSF2RA', 'MS4A7', 
    'SNAIL1', 'SNAIL2', 'ZEB1', 'ZEB2', 'CDH2', 'VIM', 'FN1', 'CDH1', 'CTNNB1',
    'KRT8', 'ERBB2',
    'FOXC2', 'JUN', 'PRRX1', 'ID2', 'SNAI2', 'MSX1')




DefaultAssay(integrated_obj) <- 'RNA'
to_plot <- markers[markers %in% rownames(integrated_obj)]
for (g in to_plot) {
    pdf(paste0(marker_dir, '/', g, '.pdf'))
        print(FeaturePlot(integrated_obj, features = g, reduction = 'umap', order = T))
    dev.off()
}

targeted_markers <- c(
    'KRT8', 'ERBB2', 'CDH1', 'FN1', 'COL1A1', 
    'ZEB1', 'ZEB2', 'FOXC2', 'JUN', 'PRRX1',  
    'ID2', 'SNAI2', 'MSX1', 'IL7R', 'ITGAX')
pdf(paste0(working_dir, '/', 'targeted_markers.pdf'), width = 20, height = 12)
    print(FeaturePlot(integrated_obj, features = targeted_markers, reduction = 'umap', order = T, ncol = 5))
dev.off()














####################################################################################################
####################################################################################################
####################################################################################################
# export cancer annotation
new_annotation <- integrated_obj@meta.data[, cluster_name, drop = FALSE]
colnames(new_annotation) <- 'new_annotation'
new_annotation$barcode <- rownames(new_annotation)
new_annotation$new_annotation <- paste0('Ca_', new_annotation$new_annotation)

write.csv(file = paste0(optimal_clustering_dir, '/', 'new_cancer_annotation.csv'), new_annotation)






# # modify the annotation in the main object
# old_annotation <- obj@meta.data[,'initial_annotation', drop = FALSE]
# colnames(old_annotation) <- 'old_annotation'
# old_annotation$barcode <- rownames(old_annotation)

# new_annotation <- integrated_obj@meta.data[, cluster_name, drop = FALSE]
# colnames(new_annotation) <- 'new_annotation'
# new_annotation$barcode <- rownames(new_annotation)
# new_annotation$new_annotation <- paste0('Ca_', new_annotation$new_annotation)

# merged_annotation <- base::merge(old_annotation, new_annotation, by = 'barcode', all = TRUE)
# merged_annotation$modified_annotation <- 
#     ifelse(is.na(as.character(merged_annotation$new_annotation)), as.character(merged_annotation$old_annotation), as.character(merged_annotation$new_annotation))
# rownames(merged_annotation) <- merged_annotation$barcode
# merged_annotation_ordered <- merged_annotation[rownames(old_annotation), ]

# obj@meta.data[['modified_cancer_annotation']] <- as.character(merged_annotation_ordered$modified_annotation)


# # plot umap to check 
# pdf(file = paste0(working_dir, '/', 'all_data_modified_annotation.pdf'), width = 10, height = 7)
#     print(
#         DimPlot(obj, group.by = 'modified_cancer_annotation', label = TRUE, repel = TRUE)
#     )
# dev.off()










# # save whole dataset with modified annotation
# saveRDS(file = "/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4066_cancer_reclustered.rds", obj)



