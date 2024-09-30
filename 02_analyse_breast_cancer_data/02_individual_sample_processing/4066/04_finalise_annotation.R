
# aim ----
    # update annotation based on subclustering




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

library(psych)
library(corrplot)















# arguments ----
data_obj_file_path <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4066_cellbender_filtered_integration_initial_annotation.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/integration/25k_downsampled/cellbender/4066/final_annotation'
dir.create(out_dir, recursive = TRUE)













#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# load data ----
obj <- readRDS(file = data_obj_file_path)



new_cancer_annot <- read.csv(file = '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/integration/25k_downsampled/cellbender/4066/subclustering/cancer_epithelial_reintegrated/optimal_cluster/new_cancer_annotation.csv',
    row.names = 1)
new_cancer_annot <- new_cancer_annot[, 'new_annotation', drop = FALSE]
new_DC_B_annot <- read.csv(file = '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/integration/25k_downsampled/cellbender/4066/subclustering/myeloid_reintegrate/optimal_cluster/new_DC_B_annotation.csv',
    row.names = 1)
colnames(new_DC_B_annot) <- c('new_annotation')

new_annotations <- rbind(new_cancer_annot, new_DC_B_annot)

new_meta_data <- base::merge(obj@meta.data, new_annotations, by = 0, all.x = TRUE)
new_meta_data$new_annotation <- ifelse(is.na(new_meta_data$new_annotation), as.character(new_meta_data$initial_annotation), new_meta_data$new_annotation)
new_meta_data <- column_to_rownames(new_meta_data, var = 'Row.names')
new_meta_data <- new_meta_data[rownames(obj@meta.data), ]


obj$new_annotation <- new_meta_data$new_annotation







# make major cell type annotations
obj$major_annotation <- 
    ifelse(as.character(obj$new_annotation) %in% c('Ca_5', 'Ca_4',
                                                    'Ca_6', 'Ca_2',
                                                    'Ca_0', 'Ca_7', 'Ca_1', 'Ca_3', 'Epithelial_proliferative'), 'Epithelial_cancer',
    ifelse(as.character(obj$new_annotation) %in% c('CAF_CXCL12', 'CAF_ECM'), 'CAF', 
    ifelse(as.character(obj$new_annotation) %in% c("T_CD8", "T_CD4"), 'T_cells', as.character(obj$new_annotation))))




# plot UMAP
pdf(file = paste0(out_dir, '/', 'umap_high_level_annotation.pdf'), width = 8)
    print(DimPlot(obj, group.by = 'major_annotation', label = TRUE, repel = TRUE,
        cols = 'Paired'))
dev.off()



pdf(file = paste0(out_dir, '/', 'umap_split_by_sample.pdf'), width = 18, height = 7)
    print(DimPlot(obj, group.by = 'sample_id', split.by = 'sample_id', ncol = 3))
dev.off()



obj$major_annotation <- ifelse(obj$major_annotation == 'Endohelial', 'Endothelial', obj$major_annotation)




























#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# DE between major cell types
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/08_prepare_data_for_complex_heatmap.R')

grouping_annotations <- c('major_annotation', 'sample_id')


plot_data_list <- prepare_data_for_heatmap(
    seurat_obj = obj,
    seurat_assay = 'RNA', 
    serat_slot = 'data',
    genes_to_plot = unique(top_markers),
    grouping_annotations = grouping_annotations,
    additional_annotations = NULL)





# put the heatmap together ----
# scale the expression matrix
plot_data_list[[1]] <- scale(t(plot_data_list[[1]]))
plot_data_list[[1]] <- t(plot_data_list[[1]])

# set expression data colors 
col_fun <- colorRamp2(c(floor(min(plot_data_list[[1]])),
                        mean(plot_data_list[[1]]),
                        ceiling(max(plot_data_list[[1]]))), 
                      c("#2166AC", "white", "#B2182B"))



# set annotation colors & make annotation ----
sample_id_colors <- set_colors(plot_data_list[[2]], 'sample_id', palette = 'Dark2')
cell_type_colors <- set_colors(plot_data_list[[2]], 'major_annotation', palette = 'Paired')




# make annotation
right_annot <- HeatmapAnnotation(
    sample_id = plot_data_list[[2]]$sample_id,
    cell_type = plot_data_list[[2]]$major_annotation,
    which = 'row',
    col = list(sample_id = sample_id_colors,
                cell_type = cell_type_colors), 
    annotation_name_gp = grid::gpar(fontsize = 12),
    annotation_legend_param = list(
        sample_id = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)),
        cell_type = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12))))







# plot top 10 genes in each cluster
top_markers_to_label <- de_res %>% 
    group_by(cluster) %>%
    filter(p_val_adj < 0.05) %>%
    slice_max(avg_log2FC, n = 5) %>% 
    pull(gene)

all_genes_to_label <- unique(top_markers_to_label)




# make sure the genes to annotate are in the heatmap matrix 
all_genes_to_label <- all_genes_to_label[all_genes_to_label %in% rownames(plot_data_list[[1]])]

# get location of the genes to be annotated
mtx_index <- seq(1:length(plot_data_list[[1]]))
names(mtx_index) <- rownames(plot_data_list[[1]])
position_in_matrix <- mtx_index[names(mtx_index) %in% all_genes_to_label]

bottom_annotation <- HeatmapAnnotation(markers = anno_mark(at = position_in_matrix, 
    labels = names(position_in_matrix),
    which = "column", side = "bottom"))



# plot heatmap ----
plot <- Heatmap(t(plot_data_list[[1]]),
                name = 'Scaled mean expression',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold")),
                col = col_fun,
                right_annotation = right_annot,
                bottom_annotation = bottom_annotation,
                show_column_names= FALSE,
                show_row_names = TRUE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                column_names_gp = gpar(fontsize = 6),
                row_names_gp = gpar(fontsize = 8))

pdf(file = paste0(out_dir, '/', 'top_DE_genes_heatmap', '.pdf'),
    width = 18, height = 7)
    draw(plot)
dev.off()






















#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# plot cell type proportion bar plots 

source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')

# modify processing method annotation 
obj$processing_method <- 
    ifelse(obj$processing_method == 'FFPE_snPATHO', 'FFPE-snPATHO-Seq',
    ifelse(obj$processing_method == 'SNAP_snPATHO', 'Frozen-Flex',
    ifelse(obj$processing_method == 'SNAP_3p', "Frozen-3'", 'unlabelled')))










cell_type_prop_summary_df <- as.data.frame(table(obj$major_annotation, obj$processing_method))
cell_type_prop_summary_df <- cell_type_prop_summary_df %>% 
    group_by(Var2) %>% 
    mutate(prop = Freq / sum(Freq)) %>% 
    as.data.frame()


# cell_type_prop_summary_df$prop <- cell_type_prop_summary_df$Freq / sum(cell_type_prop_summary_df$Freq)

cell_type_colors <- cell_type_colors[names(cell_type_colors) %in% unique(obj$major_annotation)]


plot <- ggplot(cell_type_prop_summary_df, aes(x = Var2, y = prop, color = Var1, fill = Var1)) +
            geom_col() +
            scale_color_manual(values = cell_type_colors) +
            scale_fill_manual(values = cell_type_colors) +
            theme_bw() +
            labs(x = 'processing_method',
                y = 'Proportion',
                fill = 'Cell types',
                color = 'Cell types') +
            theme(text = element_text(size = 20),
                  panel.grid = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)
                )
ggsave(file = paste0(out_dir, '/', 'cell_type_proportion_barplot.pdf'), plot, width = 7, height = 12)

















# proportion of each sample in each cluster ----
sample_contribution_df <- as.data.frame(table(obj$processing_method, obj$major_annotation))

# plot samples contribution to each cluster
total <- sample_contribution_df %>% 
    group_by(Var1) %>% 
    summarise(sum = sum(Freq))

sample_contribution_df_summarised <- base::merge(sample_contribution_df, total, by = 'Var1') 
sample_contribution_df_summarised$prop_of_total <- sample_contribution_df_summarised$Freq / sample_contribution_df_summarised$sum

source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')

plot <- ggplot(sample_contribution_df_summarised, 
    aes(x = Var1,
        y = Freq,
        fill = Var2)) +
    geom_col() +
    scale_fill_manual(values = cell_type_colors) +
    labs(x = 'processing_method',
        y = 'Numbers of cells',
        fill = 'Cell types') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(size = 20),
        panel.grid = element_blank())
ggsave(file = paste0(out_dir, '/', 'sample_contribution_cell_numbers.pdf'),
        plot, width = 7,
        height = 12, limitsize = FALSE)
















#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# plot targeted markers
cancer_epithelial <- c('ERBB2', 'GRB7', 'MKI67')
basal_epithelial <- c('KRT5', 'KRT14', 'KRT6B')
Luminal <- c('KRT23', 'LUM', 'ESR1')
myoepithelial <- c('MYLK', 'MYH11', 'MYLK')

T <- c('CD3D', 'CD8A', 'CD4', 'TRAC')
B <- c('MS4A1', 'CD19', 'CD79A')
MAST <- c('HDC', 'GATA2', 'CPA3')

DC <- c('CLEC9A', 'CD1C', 'ITGAX', 'XCR1', 'IRF7')
macrophage <- c('CD68', 'S100A9', 'S100A8', 'CSF1R')

CAF <- c('PDGFRB', 'PDGFRA', 'ACTA2')
PVL <- c('MCAM', 'RGS5')
endothelial <- c('PECAM1', 'VWF', 'EGFL7')

all_markers <- c(cancer_epithelial, basal_epithelial, Luminal, myoepithelial, T, B, MAST,
    DC, macrophage, CAF, PVL, endothelial)


# complexheatmap
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/08_prepare_data_for_complex_heatmap.R')

grouping_annotations <- c('major_annotation', 'sample_id')


plot_data_list <- prepare_data_for_heatmap(
    seurat_obj = obj,
    seurat_assay = 'RNA', 
    serat_slot = 'data',
    genes_to_plot = unique(all_markers),
    grouping_annotations = grouping_annotations,
    additional_annotations = NULL)





# put the heatmap together ----
# scale the expression matrix
plot_data_list[[1]] <- scale(t(plot_data_list[[1]]))
plot_data_list[[1]] <- t(plot_data_list[[1]])

# set expression data colors 
col_fun <- colorRamp2(c(floor(min(plot_data_list[[1]])),
                        mean(plot_data_list[[1]]),
                        ceiling(max(plot_data_list[[1]]))), 
                      c("#2166AC", "white", "#B2182B"))



# set annotation colors & make annotation ----
sample_id_colors <- set_colors(plot_data_list[[2]], 'sample_id', palette = 'Dark2')
cell_type_colors <- set_colors(plot_data_list[[2]], 'major_annotation', palette = 'Set3')




# make annotation
right_annot <- HeatmapAnnotation(
    sample_id = plot_data_list[[2]]$sample_id,
    cell_type = plot_data_list[[2]]$major_annotation,
    which = 'row',
    col = list(sample_id = sample_id_colors,
                cell_type = cell_type_colors), 
    annotation_name_gp = grid::gpar(fontsize = 12),
    annotation_legend_param = list(
        sample_id = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)),
        cell_type = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12))))




# plot heatmap ----
plot <- Heatmap(t(plot_data_list[[1]]),
                name = 'Scaled mean expression',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold")),
                col = col_fun,
                right_annotation = right_annot,
                show_column_names= TRUE,
                show_row_names = TRUE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                column_names_gp = gpar(fontsize = 6),
                row_names_gp = gpar(fontsize = 8))

pdf(file = paste0(out_dir, '/', 'selected_markers_heatmap', '.pdf'),
    width = 16, height = 7)
    draw(plot)
dev.off()


















#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# save object with final annotation
saveRDS(file = '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4066_cellbender_filtered_integration_final_annotation.rds', obj)








