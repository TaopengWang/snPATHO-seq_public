
# aim ----
    # plot a heatmap for clustered robust NMF programs 





# define environment ----
library(tidyverse)
library(readxl)

library(viridis)
library(ComplexHeatmap)
library(circlize)

library(ape)
library(phylogram)





# arguments ----
sorted_NMF_similarity_matrix_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_programs_clustering_plot_order.csv'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis'
dir.create(out_dir, recursive = T)






# load data ----
sorted_NMF_similarity_matrix <- read.csv(file = sorted_NMF_similarity_matrix_file, row.names = 1)
rownames(sorted_NMF_similarity_matrix) <- colnames(sorted_NMF_similarity_matrix)







# prepare heatmap annotation ----
annotation_df <- data.frame(programs = colnames(sorted_NMF_similarity_matrix))
annotation_df$sample <- ifelse(annotation_df$programs %in% grep('4066', annotation_df$programs, value = T), '4066',
    ifelse(annotation_df$programs %in% grep('4399', annotation_df$programs, value = T), '4399',
    ifelse(annotation_df$programs %in% grep('4411', annotation_df$programs, value = T), '4411', 'unlabelled')))

annotation_df$processing_method <- ifelse(annotation_df$programs %in% grep('FFPE', annotation_df$programs, value = T), 'snPATHO',
    ifelse(annotation_df$programs %in% grep('SNAP', annotation_df$programs, value = T), 'FLEX',
    ifelse(annotation_df$programs %in% grep('GEX', annotation_df$programs, value = T), "3'", 'unlabelled')))

annotation_df$specificity <- 
    ifelse(annotation_df$programs %in% c(
        "X4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.4.2",
        "X4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
        "X4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.6",
        "X4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.4.1",
        "X4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.3",
        "X4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.5",
        "X4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.5.3",
        "X4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.2",
        "X4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
        "X4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
        "X4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.1",
        "X4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
        "X4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.2",
        "X4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.4.4",
        "X4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.5",
        "X4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
        "X4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.4.4",
        "X4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.6",
        "X4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
        "X4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.5.1",
        "X4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.1",
        "X4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.1",
        "X4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.1",
        "X4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.6",
        "X4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.9.3",
        "X4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.6.4",
        "X4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.7",
        "X4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.4",
        "X4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.9.3",
        "X4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.6.1",
        "X4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.6.6",
        "X4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.5",
        "X4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2"
    ), 'cross-sample, cross-workflow',
    ifelse(annotation_df$programs %in% c(
        "X4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.6",
        "X4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.5",
        "X4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.6.2",
        "X4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.9.8",
        "X4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2",
        "X4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.9.5",
        "X4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.5",
        "X4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.6",
        "X4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2",
        "X4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.7",
        "X4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.6.1",
        "X4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.1",
        "X4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.2",
        "X4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.7.5",
        "X4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
        "X4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.6.2",
        "X4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.3",
        "X4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.3"
    ), 'tumour-specific, cross-workflow', ''))




# subset data to only plot NMF programs contributed to clusters with at least 3 members
annotation_df <- annotation_df[annotation_df$specificity %in% c("cross-sample, cross-workflow", "tumour-specific, cross-workflow"), ]

sorted_NMF_similarity_matrix <- sorted_NMF_similarity_matrix[annotation_df$programs, annotation_df$programs]

# hcluster on similarity matrix 
nmf_intersect_flex_filtered_hc <- hclust(as.dist(50-sorted_NMF_similarity_matrix), method="average") 

# save plot order 
phylo_tree <- as.phylo(as.dendrogram(nmf_intersect_flex_filtered_hc))
write.csv(file = paste0(out_dir, '/', 'clustered_robust_NMF_programs_heamtap_hclustered_plot_order.csv'), as.data.frame(phylo_tree$tip.label))




# plot heatmap ----
col_fun <- colorRamp2(c(0, 5, 10, 15, 20, 25), 
                      c('white', rev(inferno(n = 6))[2:length(rev(inferno(n = 6)))]))



# load preset colors for snRNA-seq workflows
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/colors.R')
workflow_colors <- workflow_colors[unique(annotation_df$processing_method)]
sample_colors <- sample_colors[unique(annotation_df$sample)]
specificity_colors <- c('#ff68b4', '#68b4ff')
names(specificity_colors) <- c('tumour-specific, cross-workflow', 'cross-sample, cross-workflow')

# make annotation
top_annot <- HeatmapAnnotation(
    Samples = annotation_df$sample,
    Workflows = annotation_df$processing_method,
    Specificity = annotation_df$specificity,
    which = 'column',
    col = list(Samples = sample_colors,
                Workflows = workflow_colors,
                Specificity = specificity_colors), 
    annotation_name_gp = grid::gpar(fontsize = 16),
    annotation_legend_param = list(
        Samples = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16)),
        Workflows = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16)),
        Specificity = list(title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 16))))



plot <- Heatmap(sorted_NMF_similarity_matrix,
                name = 'Overlapping',
                heatmap_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold")),
                col = col_fun,
                top_annotation = top_annot,
                show_column_names= FALSE,
                show_row_names = FALSE,
                cluster_rows = as.dendrogram(nmf_intersect_flex_filtered_hc),
                cluster_columns = as.dendrogram(nmf_intersect_flex_filtered_hc),
                column_names_gp = gpar(fontsize = 14),
                row_names_gp = gpar(fontsize = 14))

pdf(file = paste0(out_dir, '/', 'clustered_robust_NMF_programs_heamtap_hclustered', '.pdf'),
    width = 10.5, height = 8)
    draw(plot, merge_legend = TRUE)
dev.off()













