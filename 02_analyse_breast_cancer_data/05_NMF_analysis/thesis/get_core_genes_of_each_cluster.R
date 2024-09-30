
# aim ----
    # get genes detected in at least 2 NMF programs in each cluster



# define environment ----
library(tidyverse)
library(readxl)




# arguments ----
robust_NMF_programs_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_programs.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/thesis_analysis/core_genes_summary'
dir.create(out_dir, recursive = T)





# load data ----
robust_NMF_programs <- readRDS(file = robust_NMF_programs_file)
robust_NMF_programs <- as.data.frame(robust_NMF_programs)









########################################################################################################
########################################################################################################
########################################################################################################
cluster_1 <- c(
    "4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.4.2",
    "4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
    "4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.6",
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.4.1",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.3",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.5",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.5.3"
)

cluster_2 <- c(
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
    '4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.9.5'
)

cluster_3 <- c(
    "4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.6",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.5",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.6.2"
)

cluster_4 <- c(
    "4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.9.8",
    "4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2",
    "4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.9.5"
)

cluster_5 <- c(
    "4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.5",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.6",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2"
)

cluster_6 <- c(
    "4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.7",
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.6.1",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.1"
)

cluster_7 <- c(
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.4",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.7"
)


cluster_8 <- c(
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.2",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.1",
    "4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
    "4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.2",
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.4.4",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.5",
    "4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
    "4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.4.4"
)

cluster_9 <- c(
)




