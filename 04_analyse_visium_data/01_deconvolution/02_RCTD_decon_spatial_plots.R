
# aim ----
    # plot deconvolution results spatially





# define environment ----
library(STutility)
library(Seurat)
library(tidyverse)
library(RColorBrewer)

library(spacexr)
library(Matrix)









# arguments ----
# visium_obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/Visium/processed/STutility/4066FFPE.rds'
visium_obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/Visium/processed/STutility/4399FFPE.rds'
# visium_obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/Visium/processed/STutility/4411FFPE.rds'

# decon_res_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/deconvolution/RCTD/4411FFPE/4411FFPE_RCTD_results.rds'
# decon_res_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/deconvolution/RCTD/4399FFPE_run1/4399FFPE_run1_RCTD_results.rds'
# decon_res_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/deconvolution/RCTD/4399GEX/4399GEX_RCTD_results.rds'
decon_res_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/deconvolution/RCTD/4399SNAPFix/4399SNAPFix_RCTD_results.rds'




# 4399FFPE_run1  4399SNAPFix 4399GEX
# 4411FFPE 4411GEX 4411SNAPFix

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/deconvolution/RCTD/4399SNAPFix'
dir.create(out_dir, recursive = T)










# load visium data
visium_obj <- readRDS(file = visium_obj_file)


# load decon results
decon_res <- readRDS(file = decon_res_file)



# some visium sptos can be removed due to low counts, make sure the spots are matched
visium_obj_subset <- subset(visium_obj, cells = rownames(decon_res@results$weights))




## add raw and scaled weights to visium object ----
raw_weights <- decon_res@results$weights
# scaale raw eights
scaled_weights <- apply(raw_weights, 1, function(x){
    x <- x / sum(x)
})
scaled_weights_assay <- CreateAssayObject(data = as.matrix(scaled_weights))
visium_obj_subset[['RCTD']] <- scaled_weights_assay












# get x and y coordinates of the visium spots
image <- GetStaffli(visium_obj_subset)
image_subset <- image@meta.data[colnames(visium_obj_subset), ]
x_dim <- image_subset$x 
y_dim <- image_subset$y






# make spatial plot for results
for (ct in rownames(visium_obj_subset[['RCTD']])) {
    plot_df <- data.frame(x = x_dim,
                        y = y_dim,
                        value = visium_obj_subset@assays$RCTD@data[ct, ])
    plot <- ggplot(plot_df, aes(x = x, y = y, color = value, fill = value)) +
            geom_point(shape = 21) +
            scale_fill_gradient(low = '#e5e5e5', high = 'darkred') +
            scale_color_gradient(low = '#e5e5e5', high = 'darkred') +
            theme_bw() +
            theme(panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank(),
                aspect.ratio = (max(y_dim)*1.3) / max(x_dim),
                text = element_text(size = 16),
                plot.title = element_text(size = 30)) +
            labs(title = ct,
                fill = 'Fraction per spot',
                x = NULL,
                y = NULL) +
            guides(color = 'none')
    ggsave(file = paste0(out_dir, '/', ct, '_sclaed_porportion.pdf'),
        device = 'pdf', plot, width = 8, height = 7)
}








