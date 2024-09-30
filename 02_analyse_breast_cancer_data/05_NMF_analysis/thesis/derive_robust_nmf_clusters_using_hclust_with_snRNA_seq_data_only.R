
# aim ----
    # derive robust NMF clusters using the derived robust NMF programs 




# define environment ----
library(tidyverse)

library(viridis)
library(ComplexHeatmap)
library(circlize)

library(ape)
library(phylogram)

# library("readxl")




# arguments ----
robust_nmf_program_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_programs.rds'

# gavish_program_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/meta_programs.xlsx'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis'
dir.create(out_dir, recursive = T)







# load data ----
robust_nmf_programs <- readRDS(file = robust_nmf_program_file)
# gavish_program <- read_excel(gavish_program_file)
# colnames(gavish_program) <- paste0('Gavish - ', colnames(gavish_program))

# merged_programs <- cbind(robust_nmf_programs, gavish_program)






# cluster programs 
nmf_intersect <- apply(robust_nmf_programs , 2, function(x) apply(robust_nmf_programs , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# order programs based on their similarity
# Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=10))-1)  ) , decreasing = TRUE)




# save plot order 
# phylo_tree <- as.phylo(as.dendrogram(nmf_intersect_hc))
# write.csv(file = paste0(out_dir, '/', 'clustered_robust_NMF_programs_heamtap_hclustered_plot_order.csv'), as.data.frame(phylo_tree$tip.label))
plot_order <- colnames(nmf_intersect)
write.csv(file = paste0(out_dir, '/', 'clustered_robust_NMF_programs_heamtap_hclustered_plot_order.csv'), as.data.frame(plot_order))














###################################################################################################
###################################################################################################
###################################################################################################
# plot a heatmap with the filtered NMF programs 

# prepare heatmap annotation ----
annotation_df <- data.frame(programs = colnames(nmf_intersect))
annotation_df$sample <- ifelse(annotation_df$programs %in% grep('4066', annotation_df$programs, value = T), '4066-pBC-HER2',
    ifelse(annotation_df$programs %in% grep('4399', annotation_df$programs, value = T), '4399-mBC-TNBC',
    ifelse(annotation_df$programs %in% grep('4411', annotation_df$programs, value = T), '4411-mBC-Lum', 'Gavish_et_al_2023')))

annotation_df$processing_method <- ifelse(annotation_df$programs %in% grep('FFPE', annotation_df$programs, value = T), 'FFPE-snPATHO-Seq',
    ifelse(annotation_df$programs %in% grep('SNAP', annotation_df$programs, value = T), 'Frozen-Flex',
    ifelse(annotation_df$programs %in% grep('GEX', annotation_df$programs, value = T), "Frozen-3'", 'Gavish_et_al_2023')))




# plot heatmap ----
col_fun <- colorRamp2(c(0, 5, 10, 15, 20, 25), 
                      c('white', rev(inferno(n = 6))[2:length(rev(inferno(n = 6)))]))



# load preset colors for snRNA-seq workflows
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')
workflow_colors <- workflow_colors[unique(annotation_df$processing_method)]
sample_colors <- sample_colors[unique(annotation_df$sample)]

# make annotation
top_annot <- HeatmapAnnotation(
    Samples = annotation_df$sample,
    Workflows = annotation_df$processing_method,
    which = 'column',
    col = list(Samples = sample_colors,
                Workflows = workflow_colors), 
    annotation_name_gp = grid::gpar(fontsize = 16),
    annotation_legend_param = list(
        Samples = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16)),
        Workflows = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16))))



plot <- Heatmap(nmf_intersect,
                name = 'Overlapping',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold")),
                col = col_fun,
                top_annotation = top_annot,
                show_column_names= FALSE,
                show_row_names = FALSE,
                cluster_rows = F,
                cluster_columns = F,
                # cluster_rows = as.dendrogram(nmf_intersect_hc),
                # cluster_columns = as.dendrogram(nmf_intersect_hc),
                column_names_gp = gpar(fontsize = 14),
                row_names_gp = gpar(fontsize = 14))

pdf(file = paste0(out_dir, '/', 'robust_NMF_cluster_heatmap', '.pdf'),
    width = 10, height = 8)
    draw(plot, merge_legend = TRUE)
dev.off()





















###################################################################################################
###################################################################################################
###################################################################################################
# remove these undetected NMF programs and plot again 

to_remove <- c(
    "Gavish - Respiration",
    "Gavish - Translation initiation",
    "Gavish - MYC",
    "Gavish - Cell Cycle HMG-rich",

)











# prepare heatmap annotation ----
annotation_df <- data.frame(programs = colnames(nmf_intersect))
annotation_df$sample <- ifelse(annotation_df$programs %in% grep('4066', annotation_df$programs, value = T), '4066-pBC-HER2',
    ifelse(annotation_df$programs %in% grep('4399', annotation_df$programs, value = T), '4399-mBC-TNBC',
    ifelse(annotation_df$programs %in% grep('4411', annotation_df$programs, value = T), '4411-mBC-Lum', 'Gavish_et_al_2023')))

annotation_df$processing_method <- ifelse(annotation_df$programs %in% grep('FFPE', annotation_df$programs, value = T), 'FFPE-snPATHO-Seq',
    ifelse(annotation_df$programs %in% grep('SNAP', annotation_df$programs, value = T), 'Frozen-Flex',
    ifelse(annotation_df$programs %in% grep('GEX', annotation_df$programs, value = T), "Frozen-3'", 'Gavish_et_al_2023')))




# plot heatmap ----
col_fun <- colorRamp2(c(0, 5, 10, 15, 20, 25), 
                      c('white', rev(inferno(n = 6))[2:length(rev(inferno(n = 6)))]))



# load preset colors for snRNA-seq workflows
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')
workflow_colors <- workflow_colors[unique(annotation_df$processing_method)]
sample_colors <- sample_colors[unique(annotation_df$sample)]

# make annotation
top_annot <- HeatmapAnnotation(
    Samples = annotation_df$sample,
    Workflows = annotation_df$processing_method,
    which = 'column',
    col = list(Samples = sample_colors,
                Workflows = workflow_colors), 
    annotation_name_gp = grid::gpar(fontsize = 16),
    annotation_legend_param = list(
        Samples = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16)),
        Workflows = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16))))



plot <- Heatmap(nmf_intersect,
                name = 'Overlapping',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold")),
                col = col_fun,
                top_annotation = top_annot,
                show_column_names= FALSE,
                show_row_names = FALSE,
                cluster_rows = F,
                cluster_columns = F,
                # cluster_rows = as.dendrogram(nmf_intersect_hc),
                # cluster_columns = as.dendrogram(nmf_intersect_hc),
                column_names_gp = gpar(fontsize = 14),
                row_names_gp = gpar(fontsize = 14))

pdf(file = paste0(out_dir, '/', 'robust_NMF_cluster_heatmap', '.pdf'),
    width = 10, height = 8)
    draw(plot, merge_legend = TRUE)
dev.off()




