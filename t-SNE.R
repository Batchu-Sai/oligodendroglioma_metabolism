rm(list = ls())
library(scater)
library(reshape2)
options(stringsAsFactors=FALSE)

# Loading the tumor data ------------------------------------------------------------------
selected_sce <- readRDS("selected_sce.rds")
selected_metabolic_sce <- selected_sce[rowData(selected_sce)$metabolic,]

# T-SNE ===========================================================================
library("Rtsne")
set.seed(12345)
tsne_metabolic <- Rtsne(t(assay(selected_metabolic_sce,"exprs")),initial_dims=20,theta=0.0,perplexity = 30)
tsne_metabolic_out <- data.frame(x=tsne_metabolic$Y[,1],y=tsne_metabolic$Y[,2],location = colData(selected_metabolic_sce)$Location_General, location_specific = colData(selected_metabolic_sce)$Location)

ggplot(tsne_metabolic_out) + geom_point(aes(x, y, colour = location), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Using metabolic genes")

