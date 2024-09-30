
# aim ----
    # map gene modules on umap
    # map MsigDB and Itay's signatures on umap for comparison






# define environment ----
.libPaths(
    c(
        '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/renv/library/R-4.1/x86_64-conda-linux-gnu',
        .libPaths()
    )
)



library(Seurat)

library(RColorBrewer)
library(clusterProfiler)

library(AUCell)
library(tidyverse)





# arguments ----
obj_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4066_cellbender_filtered_integration_final_annotation.rds'
cell_ranking_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/data/snPATHO_seq/25k_downsampled_processed_objs/integrated/4066_AUCell_cellrankings.rds'

hm_gs_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/resources/pathway_gmt/h.all.v2022.1.Hs.symbols.gmt'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/map_GM_4066'
dir.create(out_dir, recursive = T)

ncores <- 12











########################################################################################################################
########################################################################################################################
########################################################################################################################
# load core genes from the robust NMF programs
NMF_cell_cycle <- c('AC091057-6','ANLN','ARHGAP11A','ARHGAP11B','ARHGAP33','ASPM','ATAD2','ATAD5','AURKB','BRCA1','BRCA2','BRIP1','BUB1','BUB1B','C21orf58','CCNB2','CCNF','CDC45','CDK1','CENPE','CENPF','CENPI','CENPK','CENPM','CENPU','CEP128','CEP152','CIT','DIAPH3','DNA2','DSCC1','E2F1','ECT2','ESPL1','EZH2','FANCA','FANCD2','FANCI','FOXM1','GMNN','HELLS','HIST1H1B','HIST1H1D','HMGB2','INCENP','IQGAP3','KIF11','KIF15','KIF18B','KIF22','KIF23','KIF2C','KIFC1','KNL1','KNTC1','LMNB1','LMNB2','MASTL','MCM3','MDC1','MELK','MIS18BP1','MKI67','MXD3','MYBL2','NCAPD2','NCAPG','NCAPG2','NDE1','NEURL1B','NUSAP1','ORC6','PKMYT1','PLK1','POLA2','POLQ','PRC1','PSRC1','RAD54B','RAD54L','RECQL4','RRM2','SMC4','SPAG5','SPC24','SPC25','STIL','TACC3','TCF19','TICRR','TMPO','TONSL','TOP2A','TPX2','TRIP13','TROAP','TYMS','UBE2C','WDR62','WDR76','XRCC2','ZNF519','ZWINT')
NMF_hypoxia <- c('ADM','AL117329-1','ANGPTL4','ARID3A','BNIP3L','CAMK1D','CAMK2N1','CASP14','CSF3R','DNAH11','EGFR','EGLN3','ELL2','ENO1','ENO2','ERO1A','ERRFI1','ESYT2','ETS2','FAM13A','FAM83D','FUT11','FZD8','GBE1','GPI','HERC3','HILPDA','HK2','IL1RAP','ILVBL','INHBA','JPH2','LIMCH1','LOXL2','LPCAT1','LRP2','MAFK','MXD1','N4BP3','NDRG1','NECTIN4','P4HA1','PDGFB','PDK1','PFKFB3','PFKFB4','PFKP','PIAS2','PKM','PLIN2','PLOD2','PPFIA4','PPP1R13L','PRELID2','PRKAA2','PROM1','RBCK1','RUNX1','SCARB1','SEMA4B','SEMA6A','SFXN3','SH3D21','SLC16A3','SLC2A1','SLC2A3','SPAG4','STC1','SYDE1','TGFBI','TNFAIP8','TNFRSF11A','TTYH3','UBC','VEGFA','YEATS2','ZMYND8')
NMF_stress <- c('AARS','ABHD3','ADM2','AHR','ALDH1L2','ANXA1','BRIX1','CARS','CEBPG','CLDN1','CPA4','DDIT3','DDX21','GARS','GPNMB','HSPA9','IARS','ID4','LSMEM1','MAFK','MAL2','MTHFD1L','MTHFD2','NOP58','PSAT1','RAB11FIP1','SDCBP2','SESN2','SLC38A1','SLC38A2','SLC6A9','SLC7A1','STC2','TARS','TRIB3','TUBE1','WDR43','YARS')
NMF_4066_unique_1 <- c('ACACB','ADGRV1','ANKRD30B','CADPS2','CAPN13','CLPSL1','CPED1','CYP21A2','DNAH14','DPY19L3','EIF4G3','GCNT2','GLI3','HOXC9','INTS6','ITPR2','KRTAP2-4','LURAP1L','MAST4','MGAT4A','MPV17','MYEOV','NAV1','NFATC4','OGT','PCNX2','PLA2R1','PLK2','PRKG1','RHOBTB3','SLC4A7','SLC4A8','TBC1D30','TTC6','ZNF512B','ZNF580')
NMF_4066_unique_2 <- c('ABCA1','ADGRF1','CD55','CEACAM5','DAAM1','DAPK2','DTNA','DTX4','DUSP10','EFCAB13','EGLN3','EXPH5','FAM214A','FGD6','GDF15','GSDMC','HILPDA','KDM3A','KLF5','KLHL24','KRT15','KRT80','MACC1','MALT1','MBOAT1','NABP1','NAMPT','NDRG1','NFKBIZ','NT5C2','OCLN','OSBPL10','PDLIM5','RAB11FIP1','RND3','RTKN2','SPNS2','STK38','TCP11L2','TLE1','TMEM44','USP53','WSB1')

NMF_list <- list()
for (t in c('NMF_cell_cycle', 'NMF_hypoxia', 'NMF_stress', 'NMF_4066_unique_1', 'NMF_4066_unique_2')) {
    NMF_list[[t]] <- eval(parse(text = t))
}




# load genes from Gavish_et_al_2023 meta programs
Gavish_G2M <- c('TOP2A','UBE2C','HMGB2','NUSAP1','CENPF','CCNB1','TPX2','CKS2','BIRC5','PRC1','PTTG1','KPNA2','MKI67','CDC20','CDK1','CCNB2','CDKN3','SMC4','NUF2','ARL6IP1','CKAP2','ASPM','PLK1','CKS1B','CCNA2','AURKA','MAD2L1','GTSE1','HMMR','UBE2T','CENPE','CENPA','KIF20B','AURKB','CDCA3','CDCA8','UBE2S','KNSTRN','KIF2C','PBK','TUBA1B','DLGAP5','TACC3','STMN1','DEPDC1','ECT2','CENPW','ZWINT','HIST1H4C','KIF23')
Gavish_G1S <- c('PCNA','RRM2','FEN1','GINS2','TYMS','MCM3','GMNN','HIST1H4C','CLSPN','ATAD2','TK1','KIAA0101','DUT','HELLS','MCM7','UBE2T','MCM4','CENPU','DHFR','ZWINT','ASF1B','MCM5','DNAJC9','RFC4','HMGB2','CDC6','RRM1','ORC6','CDK1','RAD51AP1','RNASEH2A','CHAF1A','CENPK','CDCA5','SLBP','MCM6','TMEM106C','CENPM','MYBL2','E2F1','USP1','DNMT1','PKMYT1','MAD2L1','PSMC3IP','CDCA4','RFC2','CDC45','UHRF1','MCM2')
Gavish_hypoxia <- c('NDRG1','BNIP3','SLC2A1','P4HA1','VEGFA','ENO2','ADM','HILPDA','PGK1','BNIP3L','FAM162A','PLOD2','EGLN3','LDHA','MT1X','NRN1','DDIT4','DDIT3','ERO1A','ANGPTL4','IGFBP3','PFKP','INSIG2','NDUFA4L2','PDK1','CA9','AKAP12','C4ORF3','IGFBP5','GBE1','LGALS3','TMEM45A','WSB1','GPI','S100A10','SLC16A3','EPB41L4A-AS1','ALDOA','ANKRD37','CAV1','AK4','DNAJB9','ERRFI1','MT2A','SLC3A2','VIM','ENO1','PFKFB3','HK2','CA12')
Gavish_stress <- c('ATF3','EGR1','FOS','PPP1R15A','JUN','DUSP1','JUNB','FOSB','IER2','GADD45B','BTG2','ZFP36','DNAJB1','NR4A1','SERTAD1','DDIT3','KLF6','NFKBIA','HSPA1A','HSPA1B','CYR61','MAFF','NFKBIZ','IER3','CCNL1','MCL1','RHOB','SOCS3','IRF1','CDKN1A','HES1','ID2','KLF4','KLF10','CXCL2','ID1','TOB1','TRIB1','HSPH1','RASD1','SAT1','DDIT4','EGR2','HERPUD1','PLK2','DNAJA1','HSPA6','PMAIP1','CITED2','ID3')
Gavish_stress_invitro <- c('DDIT3','SLC3A2','TRIB3','ASNS','HERPUD1','PPP1R15A','EPB41L4A-AS1','SNHG12','GARS','GADD45A','MTHFD2','XBP1','CEBPG','TAF1D','NUPR1','SARS','SNHG8','ZFAS1','DDIT4','PHGDH','SHMT2','TXNIP','ATF4','HSPA5','GDF15','CARS','MAP1LC3B','SQSTM1','YARS','DNAJB9','PSAT1','GADD45B','ATF3','SNHG15','RSRC2','CEBPB','ARF4','EIF1B','OSER1','BTG1','SNHG7','WARS','ZFAND2A','C6ORF48','STC2','ZNF622','HSPA9','PDRG1','UPP1','CCDC174')

Gavish_list <- list()
for (t in c('Gavish_G2M', 'Gavish_G1S', 'Gavish_hypoxia', 'Gavish_stress', 'Gavish_stress_invitro')) {
    Gavish_list[[t]] <- eval(parse(text = t))
}






# load genes from MsigDB
hm_gs <- read.gmt(hm_gs_file)

hm_gs_list <- list()
for (t in unique(hm_gs$term)) {
    hm_gs_list[[t]] <- hm_gs[hm_gs$term == t, 'gene']
}








# merge all sets together as a list 
all_gene_list <- c(NMF_list, Gavish_list, hm_gs_list)












########################################################################################################################
########################################################################################################################
########################################################################################################################
# map all gene sets on umap using AUCell
obj <- readRDS(file = obj_file)



# build and save cell rankings
if (file.exists(cell_ranking_file)) {
    cells_rankings <- readRDS(file = cell_ranking_file)
} else {
    # get expression matrix 
    mtx <- GetAssayData(obj, assay = 'RNA', slot = 'counts')

    # split the matrix to minimise memory usage
    max_cell <- 3000
    breaks <- seq(from = 1, to = ncol(mtx), by = max_cell)
    # build rankings for each sub matrix
    cellranking_list <- list()
    for (b in breaks) {
        next_break <- b + 3000
        if (next_break > ncol(mtx)) {
            next_break <- ncol(mtx)
            mtx_subset <- mtx[, b:next_break]
        } else {
            mtx_subset <- mtx[, b:(next_break-1)]
        }
        
        cellranking_subset <- AUCell_buildRankings(mtx_subset, nCores = ncores)
        cellranking_list[[as.character(b)]] <- cellranking_subset
    }
    # bind all sub rankings together
    cells_rankings <- do.call(AUCell::cbind, cellranking_list)
    # cells_rankings <- AUCell_buildRankings(mtx, plotStats=FALSE, nCores = ncores)
    saveRDS(file = cell_ranking_file, cells_rankings)
}






# calculate the AUCell scores
cells_AUC <- AUCell_calcAUC(geneSets = all_gene_list, cells_rankings, nCores = ncores)
saveRDS(file = paste0(out_dir, '/', 'AUCell_scores.rds'), cells_AUC)












# assign cells based on AUCell results 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, nCores = ncores, plotHist=FALSE, assign=TRUE) 
saveRDS(file = paste0(out_dir, '/', 'AUcell_assignment.rds'), cells_assignment)





# check AUCell score histograms
histogram_dir <- paste0(out_dir, '/', 'AUCell_histograms_test')
dir.create(histogram_dir, recursive = T)

for (gs in rownames(cells_AUC)) {
    threshold <- cells_assignment[[gs]]$aucThr$selected
    pdf(file = paste0(histogram_dir, '/', gs, '.pdf'))
        print(
            AUCell_plotHist(cells_AUC[gs, ], aucThr = threshold)
        )
    dev.off()
}








# plot AUCell scores on UMAP 
umap_dir <- paste0(out_dir, '/', 'AUCell_on_UMAP_test')
dir.create(umap_dir, recursive = T)

for (gs in rownames(cells_AUC)) {
    cell_assigned <- cells_assignment[[gs]]$assignment
    obj$AUCell_assignment <- ifelse(rownames(obj@meta.data) %in% cell_assigned, 'assigned', 'not_assigned')
    pdf(file = paste0(umap_dir, '/', gs, '.pdf'))
        print(
            DimPlot(obj, group.by = 'AUCell_assignment', cols = c('red', 'lightgrey'),
                order = TRUE) +
                labs(title = gs)
        )
    dev.off()
}









