
# aim ----
    # after clustering with Gavish_et_al_2023 meta programs, get core genes from different clusters
        # any gene present in more than 2 robust NMF programs with in the same cluster is considered as a core gene






# define environment ----
library(tidyverse)
library(readxl)





# arguments ----
robust_NMF_programs_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/robust_NMF_programs.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/core_genes_no_spatial'
dir.create(out_dir, recursive = T)














# load data ----
robust_NMF_programs <- readRDS(file = robust_NMF_programs_file)
robust_NMF_programs <- as.data.frame(robust_NMF_programs)





# # not summarised - only 2 programs
# "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
# "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.9.5"

# # not summarised - only 1 program
    # potentially related to estrogen response
        # ESR1, BCL2, GFRA1, RETREG1, MED13L
# "4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.9.2"

# remove potential liver gene contamination
# nmf_mixed_2 <- c(
#     "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.4",
#     "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.7"
# )



########################################################################################################
########################################################################################################
########################################################################################################
nmf_cell_cycle <- c(
    "4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.4.2",
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.4.1",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.3",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.5",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.5.3",
    "4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
    "4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.6"
)


nmf_4411_2 <- c(
    '4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.6',
    '4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.6.2',
    '4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.5')


nmf_4399 <- c(
    '4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.9.8',
    '4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.9.5',
    '4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2'
)


nmf_4411_1 <- c(
    '4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.5',
    '4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.6',
    '4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2'
)



nmf_4066 <- c(
    '4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.7',
    '4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.6.1',
    '4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.6.1'
)


nmf_mixed_3 <- c(
    "4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
    "4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.4.4",
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.4.4",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.5.5",
    "4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.4.3",
    "4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.2",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.1",
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.2",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.8"
)





nmf_mixed_1 <- c(
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.1",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.1",
    "4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.8.8",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.9.1",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.6",
    "4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.5.1",
    "4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.6"
)





nmf_hypoxia <- c(
    "4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.9.3",
    "4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.9.3",
    "4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.6.4",
    "4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.7",
    "4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.4",
    "4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.6.1",
    "4066GEX_Epithelial_cancer_rank4_9_nruns10.RDS.6.6",
    "4066FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.5",
    "4066SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.2"
)






nmf_stress <- c('4399GEX_Epithelial_cancer_rank4_9_nruns10.RDS.7.2',
    '4399FFPE_run1_Epithelial_cancer_rank4_9_nruns10.RDS.7.5',
    '4399SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.8.8')





nmf_4411_3 <- c('4411GEX_Epithelial_cancer_rank4_9_nruns10.RDS.6.2',
    '4411FFPE_Epithelial_cancer_rank4_9_nruns10.RDS.7.3',
    '4411SNAPFix_Epithelial_cancer_rank4_9_nruns10.RDS.7.3')

















all_clusters <- c("nmf_cell_cycle", "nmf_4411_2", "nmf_4399", "nmf_4411_1", "nmf_4066",
    'nmf_mixed_3', "nmf_mixed_1", "nmf_hypoxia", 'nmf_stress', "nmf_4411_3")



for (c in all_clusters) {
    working_nmf <- robust_NMF_programs[, eval(parse(text = c))]
    working_nmf_long <- gather(working_nmf, key = 'module_names', value = 'genes')
    detection_freq <- table(working_nmf_long$genes)
    core_genes <- names(detection_freq[detection_freq >= 2])
    write.csv(file = paste0(out_dir, '/', c, '.csv'), core_genes)
}















