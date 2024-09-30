
# aim ----
    # make plots for paper




# load other packages
library(Seurat)
library(tidyverse)
library(RColorBrewer)

library(ComplexHeatmap)
library(circlize)







# arguments ----
data_obj_file_path <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4399_cellbender_filtered_integration_final_annotation.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/integration/25k_downsampled/cellbender/4399/final_annotation'
dir.create(out_dir, recursive = T)




# load data ----
obj <- readRDS(file = data_obj_file_path)



# modify processing method annotation 
obj$processing_method <- 
    ifelse(obj$sample_id == '4399FFPE_run1', 'FFPE-snPATHO-Seq',
    ifelse(obj$sample_id == '4399SNAPFix', 'Frozen-Flex',
    ifelse(obj$sample_id == '4399GEX', "Frozen-3'", 'unlabelled')))




# # rename processing method to make the naming more consistent
# obj$processing_method <- ifelse(obj$sample_id == '4399FFPE_run1', 'snPATHO',
#     ifelse(obj$sample_id == '4399SNAPFix', 'FLEX', '3p'))


# obj$major_annotation <- ifelse(obj$major_annotation == 'Epithelial', 'Epithelial_cancer', obj$major_annotation)


# load predefined colors
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')

# only keep colors for cell types present
cell_type_colors <- cell_type_colors[names(cell_type_colors) %in% unique(obj$major_annotation)]
cell_type_colors <- cell_type_colors[order(names(cell_type_colors))]

# snRNA-seq workflow colors
workflow_colors <- workflow_colors[order(names(workflow_colors))]









########################################################################################################################
########################################################################################################################
########################################################################################################################
# plot umap with modified cell type colors
plot <- DimPlot(obj, group.by = 'major_annotation', cols = cell_type_colors) +
            labs(title = NULL)
plot <- LabelClusters(plot = plot, id = "major_annotation", repel = TRUE, box = TRUE, force = 1, color = 'white')

pdf(file = paste0(out_dir, '/', 'umap_high_level_annotation_paper_figure.pdf'), width = 8)
    print(
        plot
    )     
dev.off()



pdf(file = paste0(out_dir, '/', 'umap_split_by_sample_paper_figure.pdf'),
    width = 15,
    height = 5)
    print(
        DimPlot(obj, reduction = "umap", group.by = 'processing_method', 
            split.by = 'processing_method', ncol = 3, cols = workflow_colors) +
            labs(title = NULL)
    )
dev.off()





















########################################################################################################################
########################################################################################################################
########################################################################################################################
# run DE between annotated cell types
Idents(obj) <- 'major_annotation'
de_res <- FindAllMarkers(obj, assay = 'RNA', test.use = 't', only.pos = TRUE)
write.csv(file = paste0(out_dir, '/', 'de_between_major_cell_types.csv'), de_res)







# plot expression of top DE genes between clusters - annotate canonical cell type markers
# get top 200 DE genes from each cluster to plot
top_markers <- de_res %>% 
    group_by(cluster) %>%
    filter(p_val_adj < 0.05) %>%
    slice_max(avg_log2FC, n = 200) %>% 
    pull(gene)



# get mean expresion of the selected genes for each processing method 
plot_mtx <- as.data.frame(AverageExpression(obj,assays = 'RNA', slot = 'data', features = top_markers, group.by = c('sample_id', 'major_annotation'))[[1]])
plot_mtx$gene_names <- rownames(plot_mtx)
plot_mtx_long <- gather(plot_mtx, key = 'annotation', value = 'mean_expression', -gene_names)

# annotate the expression dataframe
plot_mtx_long$processing_method <- 
    ifelse(plot_mtx_long$annotation %in% grep('FFPE', plot_mtx_long$annotation, value = TRUE, ignore.case = TRUE), 'snPATHO', 
    ifelse(plot_mtx_long$annotation %in% grep('GEX', plot_mtx_long$annotation, value = TRUE, ignore.case = TRUE), '3p',
    ifelse(plot_mtx_long$annotation %in% grep('SNAPfix', plot_mtx_long$annotation, value = TRUE, ignore.case = TRUE), 'FLEX', 'unlabelled')))
plot_mtx_long$cell_type <- unlist(lapply(str_split(plot_mtx_long$annotation, '\\_'), function(x) paste(x[-1], collapse = '_')))
plot_mtx_long$cell_type <- gsub('run1_', '', plot_mtx_long$cell_type)

# scale the gene expression within each processing method
plot_mtx_long <- plot_mtx_long %>% 
        group_by(gene_names, processing_method) %>%
        mutate(scale(mean_expression)) %>% 
        as.data.frame() 
  
plot_order <- plot_mtx_long %>% 
    arrange(factor(cell_type, levels = c('RBCs', 'Epithelial_cancer', 'Endothelial',
        'LSEC', 'Lymphatic_endothelial', 'Cholangiocyte', 'Macrophage', 'Mixed_lymphocytes',
        'Hepatocyte', 'PVL', 'CAF'))) %>% 
    pull(annotation) %>% 
    unique()


scaled_mtx <- spread(plot_mtx_long[, c('gene_names', 'annotation', 'scale(mean_expression)')],
        value = 'scale(mean_expression)', key = 'annotation')
scaled_mtx <- column_to_rownames(scaled_mtx, var = 'gene_names')
# cap the scaled expression at -5 to 5
scaled_mtx[scaled_mtx > 5] <- 5
scaled_mtx[scaled_mtx < -5] <- -5
scaled_mtx <- scaled_mtx[rowSums(is.na(scaled_mtx)) == 0, ] 
scaled_mtx <- scaled_mtx[, plot_order]

# make heatmap annotation
annot_df <- plot_mtx_long[, c('annotation', 'processing_method', 'cell_type')]
annot_df <- distinct(annot_df)
annot_df <- annot_df[match(colnames(scaled_mtx), annot_df$annotation), ]

















# set colors
col_fun <- colorRamp2(c(-5, 0, 5), 
                      c("#2166AC", "white", "#B2182B"))

# set annotation colors & make annotation ----
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')
workflow_colors <- workflow_colors[names(workflow_colors) %in% unique(annot_df$processing_method)]
cell_type_colors <- cell_type_colors[names(cell_type_colors) %in% unique(annot_df$cell_type)]

# make annotation
left_annot <- HeatmapAnnotation(
        Workflows = annot_df$processing_method,
        Cell_types = annot_df$cell_type,
        which = 'row',
        col = list(Workflows = workflow_colors,
                Cell_types = cell_type_colors), 
        annotation_name_gp = grid::gpar(fontsize = 12),
        annotation_legend_param = list(
            Workflows = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)),
            Cell_types = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12))))













# annotate selected genes 
genes_to_plot <- c(
    'EPCAM', 'ESR1', 'KRT8', 'MKI67', 
    'ANXA4', 'FXYD2', 'KRT7', 'KRT19',
    'COL1A1', 'PDGFRA', 'PDGFRB',
    'RGS5', 'MCAM', 'ACTA2',
    'PECAM1', 'VWF', 'FLT1',
    'CLEC4G', 'CLEC4M', 'STAB2',
    'LYVE1', 'PROX1', 'CCL21', 'TBX1',
    'ITGAX', 'CD68', 'CD14', 'FCGR3A', 'CD163', 'SPP1',
    'ALB', 'APOA1', 'ASGR1', 'APOC3', 'HPX',
    'TRAC', 'CD3D', 'MS4A1', 'JCHAIN', 'IL7R', 'NKG7', 'PRF1',
    'SPTA1', 'HBB', 'HBA1'
)

# make sure the selected markers are in the matrix
all_markers <- genes_to_plot[genes_to_plot %in% rownames(scaled_mtx)]

# get location of the genes to be annotated
mtx_index <- seq(1:nrow(scaled_mtx))
names(mtx_index) <- rownames(scaled_mtx)
position_in_matrix <- mtx_index[names(mtx_index) %in% all_markers]

bottom_annotation <- HeatmapAnnotation(markers = anno_mark(at = position_in_matrix, 
        labels = names(position_in_matrix),
        which = "column", side = "bottom"))

# plot heatmap ----
plot <- Heatmap(t(scaled_mtx),
                name = 'Scaled mean expression',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold")),
                col = col_fun,
                left_annotation = left_annot,
                bottom_annotation = bottom_annotation,
                show_column_names= FALSE,
                show_row_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = TRUE,
                split = annot_df$processing_method,
                column_names_gp = gpar(fontsize = 6),
                row_names_gp = gpar(fontsize = 8))


pdf(file = paste0(out_dir, '/', 'top_DE_gene_canonical_marker_heatmap', '.pdf'),
    width = 14, height = 7)
    draw(plot)
dev.off()











#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# plot cell type proportion bar plots 

source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')

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
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
                text = element_text(size = 20),
                panel.grid = element_blank())
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




