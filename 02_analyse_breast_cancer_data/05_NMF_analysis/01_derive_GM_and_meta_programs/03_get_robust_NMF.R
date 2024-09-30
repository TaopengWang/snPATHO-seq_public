
# aim ----
    # define robust NMF programs from the results




# define environment ----
library(Seurat)
library(tidyverse)

library(viridis)
library(scales)



# arguments ----
NMF_list_file <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis/merged_Epithelial_cancer_module_list.rds'

out_dir <- '/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/results/gene_module_analysis'

Min_intersect_initial <- 10    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 10    # the minimal intersection cutoff for adding a new NMF to the forming cluster 
Min_group_size        <- 0     # the minimal group size to consider for defining the first NMF program in a cluster 


















####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# load NMF results
NMF_list <- readRDS(file = NMF_list_file)
# only keep top 50 genes of each program
nmf_programs <- lapply(NMF_list, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))





# select robust NMF within each sample ----
# load predefined function
source('/share/ScratchGeneral/tonwan/projects/snPATHO-Seq/scripts/function_scripts/Itay_2023_gene_module_paper/robust_nmf_programs.R')

robust_NMF_names <- robust_nmf_programs(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)

filtered_nmf_programs <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), robust_NMF_names),drop=F])
filtered_nmf_programs <- do.call(cbind, filtered_nmf_programs)

# save robust NMF programs
saveRDS(file = paste0(out_dir, '/', 'robust_NMF_programs.rds'), filtered_nmf_programs)












####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# derive meta programs ----
# calculate similarity between robust nmf programs 
nmf_intersect <- apply(filtered_nmf_programs , 2, function(x) apply(filtered_nmf_programs , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# order programs based on their similarity
Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)






# derive metaprograms
Cluster_list              <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1] > Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster, names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  # Genes in the forming MP are first chosen to be those in the first NMF. 
  # Genes_MP always has only 50 genes and evolves during the formation of the cluster
  Genes_MP <- filtered_nmf_programs[ ,names(Sorted_intersection[1])] 
  # remove selected NMF
  filtered_nmf_programs <- filtered_nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(filtered_nmf_programs))]  
  # intersection between all other NMFs and Genes_MP
  Intersection_with_Genes_MP  <- sort(apply(filtered_nmf_programs, 2, function(x) length(intersect(Genes_MP,x))), decreasing = TRUE)  
  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF 
  NMF_history <- Genes_MP   
  
  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.    
  while (Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {  
    
    Curr_cluster <- c(Curr_cluster, names(Intersection_with_Genes_MP)[1])

    ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_MP_temp <- sort(table(c(NMF_history, filtered_nmf_programs[ ,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)
    ### genes with overlap equal to the 50th gene 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   
    
    if (length(Genes_at_border) > 1) {
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study <- paste(strsplit(i, "[.]")[[1]][1 : which(strsplit(i, "[.]")[[1]] == "RDS")], collapse = ".")
        Q <- nmf_programs[[curr_study]][match(names(Genes_at_border), toupper(rownames(nmf_programs[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(nmf_programs[[curr_study]]))))], i] 
        ### sometimes when adding genes the names do not appear
        names(Q) <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(nmf_programs[[curr_study]]))))])   
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score, Q)
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score, decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      
      Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]), names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history <- c(NMF_history, filtered_nmf_programs[ ,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP <- Genes_MP_temp[1:50]

    # remove selected NMF
    filtered_nmf_programs <- filtered_nmf_programs[ ,-match(names(Intersection_with_Genes_MP)[1], colnames(filtered_nmf_programs))]  
    # intersection between all other NMFs and Genes_MP 
    Intersection_with_Genes_MP <- sort(apply(filtered_nmf_programs, 2, function(x) length(intersect(Genes_MP,x))), decreasing = TRUE)  
    
  }
  
  Cluster_list[[paste0("Cluster_", k)]] <- Curr_cluster
  MP_list[[paste0("MP_", k)]] <- Genes_MP
  k <- k+1
  
  # Remove current chosen cluster
  nmf_intersect <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect)), -match(Curr_cluster, colnames(nmf_intersect))]  
  # Sort intersection of remaining NMFs not included in any of the previous clusters
  Sorted_intersection <- sort(apply(nmf_intersect, 2, function(x) (length(which(x>=Min_intersect_initial))-1)), decreasing = TRUE)   
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}


MP_list <-  do.call(cbind, MP_list)

saveRDS(file = paste0(out_dir, '/', 'meta_programs', '.rds'), MP_list)
write.csv(file = paste0(out_dir, '/', 'meta_programs', '.csv'), as.data.frame(MP_list))



saveRDS(file = paste0(out_dir, '/', 'robust_nmf_program_clusters', '.rds'), Cluster_list)
Cluster_list <- do.call(cbind, Cluster_list)
write.csv(file = paste0(out_dir, '/', 'robust_nmf_program_clusters', '.csv'), Cluster_list)



###  Sort Jaccard similarity plot according to new clusters:

inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
  
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered

nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new]) 

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))


plot <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

ggsave(file = paste0(out_dir, '/', 'robust_NMF_programs_clustering.pdf'), width = 10, height = 10, plot, device = 'pdf')



write.csv(file = paste0(out_dir, '/', 'robust_NMF_programs_clustering_plot_order.csv'), nmf_intersect_original[inds_new,inds_new])







