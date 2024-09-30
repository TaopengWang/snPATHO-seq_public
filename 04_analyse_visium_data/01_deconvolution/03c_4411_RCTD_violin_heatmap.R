
# aim ----
    # plot decon proportion violin based on the regions they are located







# define environment ----
library(tidyverse)
library(Matrix)


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)





# arguments ----
decon_res_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/deconvolution/RCTD'

path_annotation_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/Visium/spot_annotation'

sample_id <- '4411'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/deconvolution/RCTD'












#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# load decon results ----
res_to_load <- grep(sample_id, list.dirs(decon_res_dir, full.names = F, recursive = F), value = T)

summarised_df_list <- list()
for (r in res_to_load) {
    # load the results
    working_res <- readRDS(file = paste0(decon_res_dir, '/', r, '/', r, '_RCTD_results.rds'))
    # scale the raw weights between 0 and 1
    raw_weights <- as.matrix(working_res@results$weights)
    scaled_weights <- apply(raw_weights, 1, function(x){
        x <- x / sum(x)
    })
    scaled_weights <- as.data.frame(scaled_weights)
    scaled_weights$cell_types <- rownames(scaled_weights)
    scaled_weights_long <- gather(scaled_weights, key = 'barcodes', value = 'scaled_proportion', -cell_types)
    scaled_weights_long$dataset <- r

    summarised_df_list[[r]] <- scaled_weights_long
}

summarised_df <- Reduce(rbind, summarised_df_list)
summarised_df$processing_method <- ifelse(summarised_df$dataset == grep('FFPE', unique(summarised_df$dataset), value = TRUE), 'snPATHO',
    ifelse(summarised_df$dataset == grep('GEX', unique(summarised_df$dataset), value = TRUE), '3p', 'FLEX'))

















#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# load pathology annotation ----
path_annotation_file <- grep(sample_id, dir(path_annotation_dir), value = T)
path_annotation <- read.csv(file = paste0(path_annotation_dir, '/', path_annotation_file), row.names = 1)
colnames(path_annotation) <- 'path_annotation'
rownames(path_annotation) <- paste0(rownames(path_annotation), '_1')






# add pathology annotation to the deconvolution results ----
annotated_summarised_df <- base::merge(summarised_df, path_annotation, by.x = 'barcodes', by.y = 0, all.x = TRUE)

# remove regions to be excluded due to processing artefacts
annotated_summarised_df <- annotated_summarised_df[annotated_summarised_df$path_annotation != 'exclude', ]
















#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
plot_df <- annotated_summarised_df

# get predefined colors ----
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')
workflow_colors <- workflow_colors[names(workflow_colors) %in% unique(plot_df$processing_method)]


# control order of the plot 
annotated_summarised_df$path_annotation <- factor(annotated_summarised_df$path_annotation, 
    levels = c('blood_vessel', 'cancer_adjacent_liver', 'Mixed_cancer_and_liver', 'cancer', 'Mixed_cancer_and_stroma', 
        'Mixed_cancer_and_immune', 'stroma', 'immune'))

annotated_summarised_df$cell_types <- factor(annotated_summarised_df$cell_types,
    levels = c('Cholangiocyte', 'LSEC', 'Hepatocyte', 'Epithelial_cancer', 
        'CAF', 'Endothelial', 'Macrophage', 'Mixed_lymphocytes'))



# plot violin plot for all cells types and based on pathology annotation ----
plot <- ggplot(plot_df, aes(x = cell_types , y = scaled_proportion, fill = processing_method)) +
    geom_boxplot(outlier.size = 0) +
    facet_grid(rows = vars(path_annotation)) +
    scale_fill_manual(values = workflow_colors) +
    theme_bw() +
    labs(x = 'Cell types',
        y = 'Scaled deconvolution proportion',
        fill = 'snRNA-seq workflows')

ggsave(file = paste0(out_dir, '/', sample_id, '_decon_proportion_violin.pdf'),
    plot, device = 'pdf', width = 14, height = 14)



















#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# plot results together as heatmap ----

# modify processing method names for the paper
annotated_summarised_df$processing_method <- 
    ifelse(annotated_summarised_df$processing_method == 'snPATHO', 'FFPE-snPATHO-Seq',
    ifelse(annotated_summarised_df$processing_method == 'FLEX', 'Frozen-Flex', "Frozen-3'"))

plot_df <- annotated_summarised_df
plot_df <- plot_df %>% 
    group_by(cell_types, processing_method, path_annotation) %>% 
    summarise(median_proportion = median(scaled_proportion)) 

plot_df$unique_label <- paste0(plot_df$processing_method, '_', plot_df$cell_types)


plot_mtx <- spread(plot_df[, c('unique_label', 'median_proportion', 'path_annotation')], 
    value = 'median_proportion', key = 'path_annotation')
plot_mtx <- column_to_rownames(plot_mtx, var = 'unique_label')
plot_mtx <- as.matrix(plot_mtx)

annotation_df <- data.frame(dataset_label = plot_df$unique_label, 
    processing_method = plot_df$processing_method, 
    cell_types = plot_df$cell_types)

annotation_df <- distinct(annotation_df)
rownames(annotation_df) <- annotation_df$dataset_label
annotation_df <- annotation_df[rownames(plot_mtx), ]


# control plot order
annotation_df <- arrange(annotation_df, factor(annotation_df$cell_types,
    levels = c('Epithelial_cancer', 'CAF', 'Mixed_lymphocytes', 'Endothelial', 
        'Macrophage', 'RBCs', 'Hepatocyte', 'Cholangiocyte', 'LSEC')))
plot_mtx <- plot_mtx[rownames(annotation_df), ]




# set heatmap colors
col_fun <- colorRamp2(c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), 
                      brewer.pal(name = 'YlOrRd', n = 9))




# make annotation
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')
workflow_colors <- workflow_colors[names(workflow_colors) %in% unique(annotation_df$processing_method)]
# cell_type_colors <- cell_type_colors[names(cell_type_colors) %in% unique(annotation_df$cell_types)]

top_annot <- HeatmapAnnotation(
    Workflows = annotation_df$processing_method,
    which = 'column',
    col = list(Workflows = workflow_colors), 
    annotation_name_gp = grid::gpar(fontsize = 16),
    annotation_legend_param = list(
        Workflows = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16))))





bottom_annot <- columnAnnotation(Cell_types = anno_text(annotation_df$cell_types), 
                                                    gp = gpar(fontsize = 8))




plot <- Heatmap(t(plot_mtx),
                name = 'Scaled proportion',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold")),
                col = col_fun,
                column_split = annotation_df$processing_method,
                top_annotation = top_annot,
                bottom_annotation = bottom_annot,
                show_column_names= FALSE,
                show_row_names = TRUE,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 14),
                column_names_gp = gpar(fontsize = 8))

pdf(file = paste0(out_dir, '/', sample_id, '_decon_proportion_heatmap', '.pdf'),
    width = 10, height = 6)
    draw(plot, merge_legend = TRUE)
dev.off()






