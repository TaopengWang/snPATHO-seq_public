
# aim ----
    # instead of mapping the whole GM
    # check the expression of individual gene on UMAP and spatial data to refine the genes to look at




# define environment ----
library(Seurat)
library(STutility)




# arguments ----
snRNA_obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4411_cellbender_filtered_integration_final_annotation.rds'

visium_obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/Visium/processed/STutility/4411FFPE.rds'

core_gene_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/core_genes/nmf_4411_1.csv'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/4411_unique_GMs'







# get the genes to plot 
core_genes <- read.csv(file = core_gene_file, row.names = 1)[[1]]

marker_dir <- paste0(out_dir, '/', 'GM_1_markers')
dir.create(marker_dir, recursive = T)













# check the expression of these genes on UMAP
snRNA_obj <- readRDS(file = snRNA_obj_file)
DefaultAssay(snRNA_obj) <- 'RNA'

for (m in core_genes) {
    pdf(file = paste0(marker_dir, '/', m, '_snRNA.pdf'), width = 8)
        print(
            FeaturePlot(snRNA_obj, features = m, order = T, cols = c('#e0e0e0', '#ff0000'))
        )
    dev.off()
}







# check the expression of these genes on spatial data
visium_obj <- readRDS(file = visium_obj_file)
DefaultAssay(visium_obj) <- 'Spatial'
visium_obj <- NormalizeData(visium_obj)

image <- GetStaffli(visium_obj)
x_dim <- image@meta.data$x
y_dim <- image@meta.data$y


for (m in core_genes) {
    plot_df <- data.frame(x = x_dim,
                        y = y_dim,
                        value = as.numeric(GetAssayData(visium_obj, assay = 'Spatial', slot = 'data')[m,]))

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
            labs(title = m,
                fill = 'Expression',
                x = NULL,
                y = NULL) +
            guides(color = 'none')
        ggsave(file = paste0(marker_dir, '/', m, '_spatial.pdf'),
            device = 'pdf', plot, width = 8, height = 7)

}



