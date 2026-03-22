## Script 4: 10X single cell data - differential gene expression

library(tidyverse)
library(scuttle)
library(scran)

## Step 1: load sce object, add metadata and get chromosome annotations ####
sce <- readRDS("../../../../TST_blisters/NEW_n=31/sce_post-QC_norm_dimred_int_clust3_VDJ_annot.rds")
meta <- read.csv("../../../../TST_blisters/NEW_n=31/metadata.csv") %>% dplyr::select(Sample,Gender)

# add sex to sce object
colData(sce) <- colData(sce) %>%
  as.data.frame() %>%
  left_join(meta) %>%
  DataFrame

# chromosome annotation of genes
chr.annot <- as.data.frame(rowData(sce))
chr.annot <- chr.annot %>% select(Symbol, Chromosome)
colnames(chr.annot)[1] <- "gene"

## Step 2: create pseudobulk matrices per sample/CellType combination
# remove erythrocyte and undefined cluster
sce_filt <- sce[,!sce$CellType %in% c("undefined","erythrocytes")]

# raw counts are summed 
summed <- aggregateAcrossCells(sce_filt,
                               id = DataFrame(
                                 label = sce_filt$CellType,
                                 sample = sce_filt$Sample))
# filter out sample-label combination with low cell count (e.g. 20)
summed.filt <- summed[,summed$ncells >= 20]
#check <- as.data.frame(colData(summed.filt))

# convert character to factor columns
summed.filt$Gender <- factor(summed.filt$Gender)
summed.filt$Sample <- factor(summed.filt$Sample)

## Step 3: run the differential expression analysis
# set significance thresholds
p.thr <- 0.05
padj.thr <- 0.05
lfc.thr <- 0

# using default settings for pseudoBulkDGE (robust=TRUE)
de.results <- pseudoBulkDGE(summed.filt, 
                            label = summed.filt$CellType, # specify cluster column
                            design = ~Gender, 
                            coef = "GenderMale", # this means pos.FC associated with male
                            condition = summed.filt$Gender, # specify experimental condition column; only used for abundance-based gene filtering
                            lfc = lfc.thr, # specify the log fold change filter here
                            sorted = TRUE,
                            robust = TRUE
)

# filter the output based on significance thresholds
de.filt.padj <- lapply(de.results, function(x) subset(x, FDR < padj.thr))

# add a column to indicate in which group genes are expressed higher
de.filt.padj <- lapply(de.filt.padj, function(x) 
  cbind(x, Gender_marker = ifelse(x$logFC > 0, "up in male", "up in female")))

for (i in 1:length(de.filt.padj)){
  de.filt.padj[[i]]$cluster <- rep(names(de.filt.padj[i]),times = nrow(de.filt.padj[[i]]))
}
de.filt.padj <- do.call(rbind, de.filt.padj)
de.filt.padj$gene <- row.names(de.filt.padj)
de.filt.padj <- as.data.frame(de.filt.padj)
de.filt.padj <- left_join(de.filt.padj,chr.annot)

write.csv(de.filt.padj,paste0("../../../data/SourceData_10X_DEbyCellType_NoUndefinedOrEC_PseudobulkFilter20.csv"),
          row.names = FALSE)
