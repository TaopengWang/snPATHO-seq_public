
# aim ----
    # use gene module scores to select cells
        # only keep cancer cells
        # define cells using GM score > 1
        # 
    






# define environment ----
library(Seurat)
library(tidyverse)
library(RColorBrewer)



# arguments ----
obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4411_cellbender_filtered_integration_final_annotation.rds'

nmf_gene_list_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/core_genes_no_spatial'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/4411_unique_GMs'
dir.create(out_dir, recursive = T)

min_score <- 0.5
















################################################################################################################
################################################################################################################
################################################################################################################
# get genes from the robust NMF programs
programs <- c('nmf_4411_1', 'nmf_4411_2', 'nmf_4411_3', 'nmf_hypoxia', 'nmf_cell_cycle', 'nmf_stress', 'nmf_mixed_1', 'nmf_mixed_3')

NMF_list <- list()
for (p in programs) {
    genes <- read.csv(file = paste0(nmf_gene_list_dir, '/', p, '.csv'), row.names = 1)
    NMF_list[[p]] <- genes[[1]]
}


# NMF_cell_cycle <- c('AC091057.6','ANLN','ARHGAP11A','ARHGAP11B','ARHGAP33','ASPM','ATAD2','ATAD5','AURKB','BRCA1','BRCA2','BRIP1','BUB1','BUB1B','C21orf58','CCNB2','CCNF','CDC45','CDK1','CENPE','CENPF','CENPI','CENPK','CENPM','CENPU','CEP128','CEP152','CIT','DIAPH3','DNA2','DSCC1','E2F1','ECT2','ESPL1','EZH2','FANCA','FANCD2','FANCI','FOXM1','GMNN','HELLS','HIST1H1B','HIST1H1D','HMGB2','INCENP','IQGAP3','KIF11','KIF15','KIF18B','KIF22','KIF23','KIF2C','KIFC1','KNL1','KNTC1','LMNB1','LMNB2','MASTL','MCM3','MDC1','MELK','MIS18BP1','MKI67','MXD3','MYBL2','NCAPD2','NCAPG','NCAPG2','NDE1','NEURL1B','NUSAP1','ORC6','PKMYT1','PLK1','POLA2','POLQ','PRC1','PSRC1','RAD54B','RAD54L','RECQL4','RRM2','SMC4','SPAG5','SPC24','SPC25','STIL','TACC3','TCF19','TICRR','TMPO','TONSL','TOP2A','TPX2','TRIP13','TROAP','TYMS','UBE2C','WDR62','WDR76','XRCC2','ZNF519','ZWINT')
# NMF_hypoxia <- c('ADM','AL117329.1','ANGPTL4','ARID3A','BNIP3L','CAMK1D','CAMK2N1','CASP14','CSF3R','DNAH11','EGFR','EGLN3','ELL2','ENO1','ENO2','ERO1A','ERRFI1','ESYT2','ETS2','FAM13A','FAM83D','FUT11','FZD8','GBE1','GPI','HERC3','HILPDA','HK2','IL1RAP','ILVBL','INHBA','JPH2','LIMCH1','LOXL2','LPCAT1','LRP2','MAFK','MXD1','N4BP3','NDRG1','NECTIN4','P4HA1','PDGFB','PDK1','PFKFB3','PFKFB4','PFKP','PIAS2','PKM','PLIN2','PLOD2','PPFIA4','PPP1R13L','PRELID2','PRKAA2','PROM1','RBCK1','RUNX1','SCARB1','SEMA4B','SEMA6A','SFXN3','SH3D21','SLC16A3','SLC2A1','SLC2A3','SPAG4','STC1','SYDE1','TGFBI','TNFAIP8','TNFRSF11A','TTYH3','UBC','VEGFA','YEATS2','ZMYND8')
# NMF_stress <- c('AARS','ABHD3','ADM2','AHR','ALDH1L2','ANXA1','BRIX1','CARS','CEBPG','CLDN1','CPA4','DDIT3','DDX21','GARS','GPNMB','HSPA9','IARS','ID4','LSMEM1','MAFK','MAL2','MTHFD1L','MTHFD2','NOP58','PSAT1','RAB11FIP1','SDCBP2','SESN2','SLC38A1','SLC38A2','SLC6A9','SLC7A1','STC2','TARS','TRIB3','TUBE1','WDR43','YARS')

# NMF_4411_unique_1 <- c('ABCC3','AHR','ARHGEF38','ARHGEF6','ARNT2','ASPH','BMPR1B','CASD1','CHRD','CLK1','CRISP3','CYBRD1','ELP2','ENPP1','ERN1','FRY','GLI3','ITPRID2','KLRC3','KLRC4','LDLRAD4','MAP3K20','MAPK10','MBNL2','ME3','MGAT4A','MTCL1','NAV2','NR4A2','OSMR','PGGHG','PIEZO2','PPP1R9A','PRRT2','RAPH1','SGCE','SH3BGRL','SLC39A6','SNX22','SPAG16','SPAG17','TAT','TCIM','TLK1','TPM1','TRERF1','UGDH','ZNF697')
# NMF_4411_unique_2 <- c('ABCC5','AKAP12','ARHGAP39','ARHGAP8','ATG16L2','CACNA1H','CCDC9B','DACH1','EME2','EPS15L1','ERBB4','ESR1','EVL','FMN1','GRB14','GRHL2','HIVEP3','HOXB3','HOXB4','HOXC8','HOXD4','IDUA','IL17REL','KIF12','LRRC75B','MPPED2','NBEAL2','NDST1','NECAB3','NOS1AP','OPLAH','PIGQ','PISD','PLXNB1','RAP2C','RCN3','RHPN1','SCN2A','SCUBE2','SSH3','SYTL4','THSD4','TNFRSF25','TRAPPC9','UNC5C','ZNF552','ZNF587B')
# NMF_4411_unique_3 <- c('BRF2','BRIX1','C1orf52','CCNL1','CLK4','CNNM4','COPS2','DDX3X','DNAJA4','EIF1B','EIF5','FEM1C','GNL3','HSPA4','HSPA9','HSPH1','IFRD1','KAT6A','KDM3A','KLHL21','MALSU1','NR1D1','PAK1IP1','POLR2B','POLR2C','PSMD6','RBM22','RINT1','RIOK3','SLC20A1','TAF7','THAP3','TMEM87A','TUFT1','TXNL1','USP15','UTP3','WDR43','WEE1','YTHDF2','ZFAND2A','ZNF14','ZNF622')

# NMF_list <- list()
# for (t in c('NMF_cell_cycle', 'NMF_hypoxia', 'NMF_stress', 'NMF_4411_unique_1', 'NMF_4411_unique_2', 'NMF_4411_unique_3')) {
#     NMF_list[[t]] <- eval(parse(text = t))
# }









# get module scores for each gene modules
obj <- readRDS(file = obj_file)

for (t in names(NMF_list)) {
    obj <- AddModuleScore(
        object = obj,
        features = NMF_list[t],
        assay = 'RNA',
        name = t)
    obj@meta.data[[t]] <- obj@meta.data[[paste0(t, '1')]] 
    obj@meta.data <- obj@meta.data[, !(colnames(obj@meta.data) == paste0(t, '1'))]
}














# check assignmetn ----
# only keep cancer cells 
subset_obj <- subset(obj, subset = major_annotation == 'Epithelial_cancer')


# get maximum score each cell has been assigned
score_columns <- grep('NMF', colnames(subset_obj@meta.data), value = T, ignore.case = T)
max_scores <- apply(subset_obj@meta.data[, score_columns], 1, function(x) max(x))
subset_obj$max_scores <- max_scores
# remove cells with no score over minimal requirement
subset_obj <- subset(subset_obj, subset = max_scores > min_score)
# assign cells based on maximum scores
assignment <- apply(subset_obj@meta.data[, score_columns], 1, function(x){
    x <- as.numeric(x)
    score_columns[which(x == max(x))]
})
subset_obj$assigment <- assignment







# plot to check assignment ----
pdf(file = paste0(out_dir, '/', 'GM_assignment.pdf'), width = 8, height = 7)
    print(
        DimPlot(subset_obj, group.by = 'assigment')
    )
dev.off()




saveRDS(file = paste0(out_dir, '/', '4411_GM_assigned.rds'), subset_obj)









# plot assignment on all cells - all other cells label as unassigned
meta_df <- obj@meta.data
subset_meta_df <- subset_obj@meta.data
merged_meta_df <- base::merge(meta_df, subset_meta_df[, 'assigment', drop = FALSE], by = 0, all.x = T)
merged_meta_df$assigment <- ifelse(is.na(merged_meta_df$assigment), 'no_assignment', merged_meta_df$assigment)

merged_meta_df <- column_to_rownames(merged_meta_df, var = 'Row.names')
merged_meta_df <- merged_meta_df[match(rownames(obj@meta.data), rownames(merged_meta_df)),]
obj@meta.data <- merged_meta_df

colors <- brewer.pal(name = 'Set3', n = 6)
colors <- c(colors, '#e5e5e5')
names(colors) <- c('nmf_cell_cycle', 'nmf_hypoxia', 'nmf_stress',
    'nmf_4411_1', 'nmf_4411_2', 'nmf_4411_3',
    'no_assignment')


pdf(file = paste0(out_dir, '/', 'GM_assignment_paper.pdf'), width = 8, height = 7)
    print(
        DimPlot(obj, group.by = 'assigment', cols = colors, shuffle = TRUE)
    )
dev.off()















################################################################################################################
################################################################################################################
################################################################################################################
# DE between cancer populations 
de_obj <- subset(obj, subset = major_annotation == 'Epithelial_cancer')

Idents(de_obj) <- 'assigment'
de_res <- FindAllMarkers(de_obj, assay = 'RNA', slot = "data", logfc.threshold = 0, only.pos = FALSE, return.thresh = Inf)
saveRDS(file = paste0(out_dir, '/', 'DE_between_GM_cells.rds'), de_res)



# ################################################################################################################
# ################################################################################################################
# ################################################################################################################
# # DE between cancer populations 

# # subset further to only keep cells with unique gene modules 
# subset_obj <- subset(subset_obj, subset = assigment %in% c('NMF_4411_unique_1', 'NMF_4411_unique_2', 'NMF_4411_unique_3'))

# # only keep genes in the FLEX panel
# flex_probe_set <- read.csv(file = flex_probe_file, skip = 5)
# all_flex_genes <- flex_probe_set[flex_probe_set$include == TRUE, 'probe_id']
# all_flex_genes <- unique(unlist(lapply(str_split(all_flex_genes, '\\|'), '[[', 2)))

# Idents(subset_obj) <- 'assigment'
# de_res <- FindAllMarkers(subset_obj, features = all_flex_genes,
#     min.pct = 0, assay = 'RNA', logfc.threshold = 0, 
#     only.pos = FALSE, return.thresh = Inf)
# saveRDS(file = paste0(out_dir, '/', 'DE_between_GM_cells.rds'), de_res)























