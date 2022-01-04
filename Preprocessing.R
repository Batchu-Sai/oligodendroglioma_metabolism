rm(list = ls())
library(Seurat)
library(RCurl)
library(dplyr)
library(fgsea)
library(scater)
source("utils.R")
options(stringsAsFactors = F)

# Load Data and change colnames ----------------------------------------------------------------------------------------------------------------------
# Data downloaded from https://singlecell.broadinstitute.org/single_cell/study/SCP12/oligodendroglioma-intra-tumor-heterogeneity#/

counts <- read.table("Data/SCP12/expression/OG_processed_data_portal.txt", sep = "\t", header = 1, row.names = 1)
cell_metadata <- read.table("Data/SCP12/metadata/cell_type_assignment_portal.txt", sep = "\t", header = 1)
patient_metadata <- read.table("Patient_MetaData.csv", sep = ",", header=1)

# change some cell names that don't have "MGH' in front of them
cell_metadata[grep("^93_", cell_metadata$NAME),]$NAME <- paste0("MGH", cell_metadata[grep("^93_", cell_metadata$NAME),]$NAME)
cell_metadata[grep("^97_", cell_metadata$NAME),]$NAME <- paste0("MGH", cell_metadata[grep("^97_", cell_metadata$NAME),]$NAME)

cell_metadata[grep("^93_", cell_metadata$NAME),]$NAME <- paste0("MGH", cell_metadata[grep("^93_", cell_metadata$NAME),]$NAME)
cell_metadata[grep("^97_", cell_metadata$NAME),]$NAME <- paste0("MGH", cell_metadata[grep("^97_", cell_metadata$NAME),]$NAME)

names <- colnames(counts)
names[grep("^X97_", names)] <- gsub("X97", "MGH97", names[grep("^X97_", names)])
names[grep("^X93_", names)] <- gsub("X93", "MGH93", names[grep("^X93_", names)])
colnames(counts) <- names

# keep only cells that are tumor cells
cell_metadata <- cell_metadata[cell_metadata$CLUSTER == 'malignant',]
counts <- counts[,colnames(counts) %in% cell_metadata$NAME]

# counts are in log2(TPM/10+1), so convert them to TPM only
counts <- (2^counts -1)*10
all_data <- counts

# Get col_data for scater object and duplicate rows to match cell numbers -------------------------------------------------------------------------------------------------------
cell_metadata$Patient <- gsub("_.*", "", cell_metadata$NAME)
table(cell_metadata$Patient)
col_data <- merge(patient_metadata, cell_metadata, by.x = "Designation", by.y = "Patient")
rownames(col_data) <- col_data$NAME

# mark the metabolic genes --------------------------------------------------------------------------------------------------------------------------------------------------------------------
dim(all_data)
all_data <- all_data[,colnames(all_data) %in% col_data$NAME,]
dim(all_data)
all_data <- data.matrix(all_data)
pathways <- gmtPathways("KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE

# Build scater object --------------------------------------------------------------------------------------------------------------------------------------------------------------------

sce <- SingleCellExperiment(
  assays = list(tpm=all_data,exprs=all_data),
  colData = col_data,
  rowData = row_data
)

# Save final scater object --------------------------------------------------------------------------------
saveRDS(sce,"selected_sce.rds")

