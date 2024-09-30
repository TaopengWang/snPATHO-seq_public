
# aim ----
    # map robust NMF programs on Visium data


# define environment ----
library(Seurat)
library(STutility)

library(RColorBrewer)
library(tidyverse)







# arguments ----
visium_obj_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/Visium/processed/STutility'

robust_nmf_program_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_programs.rds'

samples <- c('4066', '4399', '4411')


out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis/map_robust_nmf_programs_spatial'
dir.create(out_dir, recursive = T)






# load object ----
obj_list <- list()
for (s in samples) {
    obj <- readRDS(file = paste0(visium_obj_dir, '/', s, 'FFPE.rds'))
    obj <- NormalizeData(obj)
    obj_list[[s]] <- obj
}


# load NMF gene lists ----
robust_nmf_program <- readRDS(file = robust_nmf_program_file)

gene_list <- list()
for (cl in colnames(robust_nmf_program)) {
    gene_list[[cl]] <- robust_nmf_program[, cl]
}










# get module scores
for (o in names(obj_list)) {
    working_obj <- obj_list[[o]]

    working_out <- paste0(out_dir, '/', o, '_seurat_module_scores')
    dir.create(working_out, recursive = T)

    image <- GetStaffli(working_obj)
    x_dim <- image@meta.data$x
    y_dim <- image@meta.data$y



    for (n in names(gene_list)) {
        annotated_obj <- AddModuleScore(
            object = working_obj,
            features = gene_list[n],
            assay = 'Spatial')
        working_obj@meta.data[n] <- annotated_obj$Cluster1

        plot_df <- data.frame(x = x_dim,
                            y = y_dim,
                            value = working_obj@meta.data[[n]])

        plot <- ggplot(plot_df, aes(x = x, y = y, color = value, fill = value)) +
            geom_point(shape = 21) +
            scale_fill_gradient2(low = '#1b1b1b', mid = '#e6e6e6', high = 'red') +
            scale_color_gradient2(low = '#1b1b1b', mid = '#e6e6e6', high = 'red') +
            theme_bw() +
            theme(panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank(),
                aspect.ratio = (max(y_dim)*1.7) / max(x_dim),
                text = element_text(size = 20)) +
            labs(title = NULL,
                fill = 'Gene module score',
                x = NULL,
                y = NULL) +
            guides(color = 'none')
        ggsave(file = paste0(working_out, '/', n, '.pdf'),
            device = 'pdf', plot, width = 8, height = 7)
    } 
}
















# plot 4411 again as the image is more elongated
for (o in names(obj_list)[[3]]) {
    working_obj <- obj_list[[o]]

    working_out <- paste0(out_dir, '/', o, '_seurat_module_scores')
    dir.create(working_out, recursive = T)

    image <- GetStaffli(working_obj)
    x_dim <- image@meta.data$x
    y_dim <- image@meta.data$y



    for (n in names(all_gene_list)) {
        annotated_obj <- AddModuleScore(
            object = working_obj,
            features = all_gene_list[n],
            assay = 'Spatial',
            name = n)
        annotated_obj@meta.data[[n]] <- annotated_obj@meta.data[[paste0(n, '1')]] 

        plot_df <- data.frame(x = x_dim,
                            y = y_dim,
                            value = annotated_obj@meta.data[[n]])

        plot <- ggplot(plot_df, aes(x = x, y = y, color = value, fill = value)) +
            geom_point(shape = 21) +
            scale_fill_gradient2(low = '#1878b2', mid = '#e5e5e5', high = '#B2182B') +
            scale_color_gradient2(low = '#1878b2', mid = '#e5e5e5', high = '#B2182B') +
            theme_bw() +
            theme(panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank(),
                aspect.ratio = (max(y_dim)*1.3) / max(x_dim)) +
            labs(title = n,
                fill = 'Gene module score',
                x = NULL,
                y = NULL) +
            guides(color = 'none')
        ggsave(file = paste0(working_out, '/', n, '.pdf'),
            device = 'pdf', plot, width = 8, height = 7)
    } 
}


    
