
# aim ----
    # plot dotplot reflecting the expression of canonical cell type markers



# define environment ----
library(Seurat)
library(tidyverse)
library(RColorBrewer)






# arguments ----
data_obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/external_data/FHIL_Data/processed_data/PBMC_integrated_annotation_modifed_by_subclustering.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/snRNA_external_data/FHIL_PBMC/dotplot_cell_markers'
dir.create(out_dir, recursive = T)













# load data ----
data_obj <- readRDS(file = data_obj_file)


t_cells <- c('CD3D', 'CD4', 'CD8A', 'TRAC')
memory <- c('ITGB1', 'ITGA4')
naive_cm_t <- c('CCR7', 'SELL')
cytotoxic <- c('GZMA', 'PRF1', 'GZMB')

Treg <- c('FOXP3', 'LAG3', 'TIGIT')

b_cells <- c('MS4A1', 'CD79A')
naive_B <- c('IGHM', 'IGHD')
memory_b <- c('BANK1', 'AIM2', 'SSPN', 'IGHG1', 'IGHA1')

nk <- c('NCAM1', 'NCR1', 'NKG7', 'CX3CR1')

myeloid <- c('CD14', 'ITGAX')
CD14_monocytes <- c('VCAN', 'S100A9')
CD16_monocytes <- c('FCGR3A')
cDC <- c('CLEC9A', 'XCR1', 'CD1C', 'CLEC10A')
pDC <- c('JCHAIN', 'IL3RA')

hspc <- c('GATA2', 'CD34', 'SOX4')
ILC <- c('KIT', 'SOX4')
MAIT <- c('KLRB1', 'TRAJ33')
platelet <- c('PPBP', "PF4")

markers_to_plot <- unique(c(
    b_cells, naive_B, memory_b, 
    myeloid, CD14_monocytes, CD16_monocytes, cDC, pDC,
    t_cells, memory, naive_cm_t, cytotoxic, Treg, 
    nk, hspc, ILC, MAIT,
    platelet  
))






 





data_obj$assay_type <- 
    ifelse(data_obj$sample_id %in% c('FLEX_exp_1', 'FLEX_exp_2'), 'Frozen - Flex', "Frozen - 3'")
data_obj$plot_annot <- paste0(as.character(data_obj$new_annotation), ' - ', as.character(data_obj$assay_type))

working_df <- data_obj@meta.data[, c('plot_annot', 'new_annotation', 'assay_type')]
working_df <- working_df %>% 
    arrange(new_annotation, assay_type)

data_obj$plot_annot <- factor(data_obj$plot_annot, levels = unique(working_df$plot_annot))

Idents(data_obj) <- 'plot_annot'

pdf(file = paste0(out_dir, '/', 'Dotplot_canonical_pbmc_markers.pdf'), width = 18, height = 12)
    print(
        DotPlot(data_obj, features = markers_to_plot, assay = 'RNA',
            cols = c("lightgrey", "red"),
            scale = F) +
            theme(text = element_text(size = 24),
                axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
    )
dev.off()













