rm(list=ls())
source("utils.R")
library(stringr)
library(reshape2)
library(scales)
library(scater)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
pathway_file <- "KEGG_metabolism.gmt"

#1. Loading the data
selected_impute_sce <- readRDS("selected_sce.rds")
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
all_locations <- as.vector(selected_impute_sce$Location_General)
locations <- unique(all_locations)

#some genes occur in multiple pathways.
gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_impute_sce)[rowData(selected_impute_sce)$metabolic])
set.seed(123)
norm_tpm <- readRDS("Normalization/Deconvolution_tpm.rds")

##Calculate the pathway activities mean ratio of genes in each pathway for each cell type ----------------------------------------
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(locations),dimnames = list(pathway_names,locations))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(locations),dimnames = list(pathway_names,locations))

###calculate the pvalues using shuffle method ------------------------
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(locations),dimnames = (list(pathway_names, locations)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next
  
  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_locations, mean))
  
  #remove genes which are zeros in any location to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]
  
  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
  
  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_locations, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_locations, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[locations]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[locations]
  
  ##shuffle 1000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(locations,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_locations_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:1000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_locations_list <- lapply(times,function(x) sample(all_locations)) 
  names(shuffle_locations_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(locations),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- locations
  for(c in locations){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 1000 
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 1000
    }
    if(pval>0.05) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}

all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]

#### Draw the heatmap
dat <- mean_expression_shuffle
colnames(dat) <- gsub("_", " ", colnames(dat))
sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1
pdf(file.path(outDir, "KEGGpathway_activity_heatmap_group.pdf"),onefile=T,width=6,height=9)
mybreaks <- c(
  seq(0, 0.75, length.out=33),
  seq(0.76, 1.2, length.out=33),
  seq(1.21, 2,length.out=34)
) 
color <- colorRampPalette(c("blue","white","red"))(100)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,breaks=mybreaks)
dev.off()

#saveRDS(dat, file.path(outDir,"heatmap_location.rds"))
#write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")