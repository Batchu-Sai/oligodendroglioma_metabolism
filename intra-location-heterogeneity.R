rm(list = ls())
library(scater)
library(stringr)
library(pheatmap)
library(gtools)
library(scran)
source("utils.R")
source("runGSEA_preRank.R")
options(stringsAsFactors=FALSE)
pathway_file <- "KEGG_metabolism.gmt"

#1. Loading the data
selected_tumor_sce <- readRDS("selected_sce.rds")
selected_tumor_metabolic_sce <- selected_tumor_sce[rowData(selected_tumor_sce)$metabolic,]

#=========================================================================
location <- unique(selected_tumor_sce$Location_General)
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)
pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),tumor=character())

for (t in location){
  each_metabolic_sce <- selected_tumor_metabolic_sce[,selected_tumor_metabolic_sce$Location_General==t]
  each_metabolic_tpm <- assay(each_metabolic_sce,"exprs")
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- as.matrix(each_metabolic_tpm)
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  
  ###select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>=0.8)[1]
  
  ###plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             tumor=rep(t,length(percentVar)))
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
  
  ####
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  
  runGSEA_preRank(pre_rank_matrix,pathway_file,t)
  #get the result
  result_dir <- list.files(path="preRankResults",pattern = paste0("^",t,".GseaPreranked(.*)"),full.names=T)
  result_file <- list.files(path=result_dir,pattern="gsea_report_for_na_pos_(.*).xls",full.names=T)
  gsea_result <- read.table(result_file,header = T,sep="\t",row.names=1)
  gsea_pathways <- str_to_title(rownames(gsea_result))
  gsea_pathways <- str_replace(gsea_pathways,"Tca","TCA")
  gsea_pathways <- str_replace(gsea_pathways,"Gpi","GPI")
  enrich_data_df <- rbind(enrich_data_df,data.frame(x=t,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$NOM.p.val))
}

#write.csv(enrich_data_df, quote = F, "enriched_data_location.csv", row.names = F)

#remove pvalue <0.05 pathways
min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN=min)
select_pathways <- names(min_pval)[(min_pval<=0.05)]
select_enrich_data_df <- enrich_data_df[enrich_data_df$y%in% select_pathways,]

#converto pvalue to -log10
pvals <- select_enrich_data_df$PVAL
pvals[pvals<=0] = 1e-10
select_enrich_data_df$PVAL <- -log10(pvals)
select_enrich_data_df$x <- gsub("_", " ", select_enrich_data_df$x)

#sort
pathway_pv_sum <- by(select_enrich_data_df$PVAL,select_enrich_data_df$y,FUN=sum)
names(pathway_pv_sum)[names(pathway_pv_sum) == 'Glycolysis  Gluconeogenesis'] <- 'Glycolysis / Gluconeogenesis'
pathway_order <- names(pathway_pv_sum)[order(pathway_pv_sum,decreasing = T)]

#top 10
pathway_order <- pathway_order[1:10]
select_enrich_data_df <- select_enrich_data_df[select_enrich_data_df$y %in% pathway_order,]
select_enrich_data_df$x <- factor(select_enrich_data_df$x, levels = mixedsort(location))
select_enrich_data_df$y <- factor(select_enrich_data_df$y,levels = pathway_order)

##bubble plot
p <- ggplot(select_enrich_data_df, aes(x = x, y = y))+
  geom_point(aes(size = PVAL, fill = NES), alpha = 0.75, shape = 21) +
  labs(x = NULL, y = NULL,size = "-log10 pvalue", fill = "NES") +
  scale_size(range = c(0, 10)) +
  scale_fill_gradient( low = "white", high = "blue") +
  theme(legend.key=element_blank(),
        legend.position = "right", 
        legend.key.size = unit(1.0, "cm"),
        legend.text = element_text(size = 10, face ="bold", colour ="black"),
        legend.title = element_text(size = 10, face = "bold"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold", angle = 0, vjust = 0, hjust = 0.5), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        panel.border = element_blank(), 
        panel.background = element_blank())
p
ggsave(file.path(outDir,"group_enriched_pathway.pdf"),p,width = 3.8,height=2.3,units="in",device="pdf",useDingbats=FALSE)

##plot variance
pc_plotdata$tumor <- factor(pc_plotdata$tumor,levels=mixedsort(location))
p <- ggplot(pc_plotdata) + geom_point(aes(x,y,colour=factor(sel)),size=0.5) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~tumor,scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(colour="black", size = 6),
        strip.background = element_rect(fill="white",size=0.2,colour = NULL),
        strip.text=element_text(size=6))
p
ggsave(file.path(outDir,"malignant_PC_variance_plot.pdf"),p,device="pdf",useDingbats=FALSE)
unlink("preRankResults",recursive=T)
unlink("prerank.rnk")
date_string <- Sys.Date()
date_split <- strsplit(as.character(date_string),"-")[[1]]
unlink(paste0(tolower(month.abb[as.numeric(date_split[2])]),date_split[3]),recursive=T)
