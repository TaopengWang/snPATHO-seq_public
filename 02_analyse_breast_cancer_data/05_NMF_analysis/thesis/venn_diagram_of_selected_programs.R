
# aim ----
    # plot venndiagram of selected programs





# load packages ----
library(tidyverse)
library(ggVennDiagram)



# arguments ----
roubst_nmf_program_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_paper_format.csv'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis/venn_diagram'
dir.create(out_dir, recursive = T)








# load data ----
roubst_nmf_program <- read.csv(roubst_nmf_program_file, header = F)
colnames(roubst_nmf_program) <- roubst_nmf_program[1, ]
roubst_nmf_program <- roubst_nmf_program[2:nrow(roubst_nmf_program),]














####################################################################################
####################################################################################
####################################################################################
# make venn diagram for selected programs - C3
plot_list <- list(
    roubst_nmf_program[['4066_3p_Epithelial_cancer_rank4_9_nruns10.RDS.8.7']],
    roubst_nmf_program[['4066_FLEX_Epithelial_cancer_rank4_9_nruns10.RDS.6.1']],
    roubst_nmf_program[['4066_FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.6.1']]
)

names(plot_list) <- c("Frozen-3'", "Frozen-Flex", "FFPE-snPATHO-seq")

plot <- ggVennDiagram(plot_list) + 
        scale_fill_gradient(low="grey90",high = "red")
ggsave(file = paste0(out_dir, '/', 'venn_diagram_C3.pdf'), plot, width = 8)


# common
intersect(intersect(plot_list[[1]], plot_list[[2]]), plot_list[[3]])
# "MYEOV"   "PRKG1"   "ACACB"   "LURAP1L" "ITPR2"   "GCNT2"   "CAPN13" 
# [8] "MAST4"   "PLK2"

# 3' specific
setdiff(setdiff(plot_list[[1]], plot_list[[2]]), plot_list[[3]])
# "HPSE2"       "VSNL1"       "RHBDL3"      "CLSTN2"      "KCNIP4"     
#  [6] "ECE1"        "AC091646.1"  "EDARADD"     "AC018816.1"  "SCD"        
# [11] "CADM2"       "GALNT14"     "RERG"        "TMTC1"       "PDE10A"     
# [16] "TBC1D4"      "WNK2"        "TMEM178B"    "FIRRE"       "ROBO1"      
# [21] "KCNMA1"      "CASC11"      "TLN2"        "ADCY2"       "AC012501.2" 
# [26] "GALNT7"      "KCNQ5"       "RERE"        "ERBB4"       "AC078923.1" 
# [31] "ACADSB"      "ASCC3"       "GHR"         "PDCD4"       "LURAP1L-AS1"
# [36] "RORA"        "DEGS1"       "GTSCR1"      "GUCY1A1"     "COBL"

# Flex specific
setdiff(setdiff(plot_list[[2]], plot_list[[1]]), plot_list[[3]])
# "KLF8"     "ZNF283"   "NICN1"    "MPHOSPH6" "ZNF680"   "PCCA"    
#  [7] "CCDC7"    "ZNF605"   "TARBP1"   "ELL2"     "MNAT1"    "UGGT2"   
# [13] "GLYCTK"   "ANKRD30A" "ZNF587B"

# snPATHO-seq
setdiff(setdiff(plot_list[[3]], plot_list[[1]]), plot_list[[2]])
# "FBXL8"   "CAB39L"  "LYPD6B"  "SMYD3"   "C5orf30" "MOGS"    "ELP2"   
#  [8] "BHLHE40" "ZNF493"  "SLC7A6"  "TMEM267" "SYNE4"   "ACACA"   "COQ4"

# Flex and snPATHO-seq common
setdiff(intersect(plot_list[[2]], plot_list[[3]]), plot_list[[1]])
# "TTC6"     "ADGRV1"   "ANKRD30B" "CYP21A2"  "ZNF580"   "RHOBTB3" 
#  [7] "SLC4A8"   "INTS6"    "PLA2R1"   "DPY19L3"  "EIF4G3"   "SLC4A7"  
# [13] "PCNX2"    "KRTAP2-4" "OGT"      "CPED1"    "DNAH14"   "GLI3"    
# [19] "HOXC9"    "CLPSL1"   "TBC1D30"  "MPV17"    "CADPS2"   "MGAT4A"  
# [25] "ZNF512B"  "NFATC4"









####################################################################################
####################################################################################
####################################################################################
# make venn diagram for selected programs - C5
plot_list <- list(
    roubst_nmf_program[['4399_3p_Epithelial_cancer_rank4_9_nruns10.RDS.9.8']],
    roubst_nmf_program[['4399_FLEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.2']],
    roubst_nmf_program[['4399_FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.5']]
)

names(plot_list) <- c("Frozen-3'", "Frozen-Flex", "FFPE-snPATHO-seq")

plot <- ggVennDiagram(plot_list) + 
        scale_fill_gradient(low="grey90",high = "red")
ggsave(file = paste0(out_dir, '/', 'venn_diagram_C5.pdf'), plot, width = 8)





# common
intersect(intersect(plot_list[[1]], plot_list[[2]]), plot_list[[3]])
# "TNC"     "P3H2"    "LDLRAD4" "ANTXR1"  "TNS3"    "DZIP1L"

# 3' specific
setdiff(setdiff(plot_list[[1]], plot_list[[2]]), plot_list[[3]])
# "MIR924HG"   "UNC5C"      "DLGAP1"     "LINC01122"  "AC011447.3"
#  [6] "DIAPH3"     "CEP128"     "ZNF519"     "CACNA1C"    "LINC01572" 
# [11] "MIR181A1HG" "ZNF90"      "CYYR1"      "PDE4B"      "AP001347.1"
# [16] "FANCI"      "MELK"       "SESN3"      "SPATA5"     "CENPP"     
# [21] "PDE7B"      "PLCB4"      "MIR205HG"   "MEGF9"      "MMS22L"    
# [26] "PDE3A"      "MIR99AHG"   "SLC24A3"    "ADGRA3"     "PBX3"      
# [31] "MIR100HG"   "ZNF618"     "CENPF"      "CCDC88A"    "NRCAM"     
# [36] "SNHG14"     "PRIM2"      "MIS18BP1"   "TMEM108"    "ATAD2"     
# [41] "ZNF608"

# Flex specific
setdiff(setdiff(plot_list[[2]], plot_list[[1]]), plot_list[[3]])
# "MDFI"     "NFATC1"   "TNS4"     "ATAD3C"   "ZNF362"   "SLC12A4" 
#  [7] "KALRN"    "TNFRSF25" "KCNC3"    "NCKAP5L"  "LCAT"     "UNC13D"  
# [13] "FOXO6"    "SHFL"     "PLXNA1"   "GNAQ"     "LRRC75B"  "PLEKHG5" 
# [19] "PML"      "TNNI2"    "PDE9A"

# snPATHO-seq
setdiff(setdiff(plot_list[[3]], plot_list[[1]]), plot_list[[2]])
# "STEAP4"  "GPRC5A"  "FLNA"    "SKIL"    "TPM2"    "PARP9"   "SULF2"  
#  [8] "EFNA1"   "ABCA1"   "JUNB"    "ZFP36L2" "STOM"    "NNMT"    "COL4A2" 
# [15] "PLK2"    "AACS"    "TUBA4A"  "OSMR"    "SMAD7"   "SMAD6"   "FN1"    
# [22] "FHL2"




# Flex and snPATHO-seq common
setdiff(intersect(plot_list[[2]], plot_list[[3]]), plot_list[[1]])
# "ADAMTS10"    "PHLDB1"      "CDH3"        "GABRE"       "IL18BP"     
#  [6] "RND1"        "RNF121"      "VWA1"        "COL5A1"      "BGN"        
# [11] "MFAP2"       "FAM189A2"    "CSF1R"       "WNK2"        "PMEPA1"     
# [16] "MCF2L"       "INPP5B"      "CADPS2"      "PALM2-AKAP2" "ANKRD10"    
# [21] "PRKCD" 


