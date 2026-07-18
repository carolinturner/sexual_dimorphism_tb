## Script 7: Differential gene expression analysis - Latent TB: Day 2 TST (single cell)

library(HDF5Array)
library(scater)
library(tidyverse)
library(scran)
library(openxlsx)

## Step 1: load sce object and get chromosome annotations ####
sce <- loadHDF5SummarizedExperiment("../../../data/Processed_sce_TSTD2_LatentTB_h5")

# chromosome annotation of genes
chr.annot <- as.data.frame(rowData(sce))
chr.annot <- chr.annot %>% dplyr::select(Symbol, Chromosome)
colnames(chr.annot)[1] <- "gene"

## Step 2: create pseudobulk matrices per sample/CellType combination
# remove erythrocyte and undefined cluster
sce_filt <- sce[,!sce$CellType %in% c("undefined","erythrocytes")]

# raw counts are summed 
summed <- scuttle::aggregateAcrossCells(sce_filt,
                               id = DataFrame(
                                 label = sce_filt$CellType,
                                 sample = sce_filt$Sample))
# filter out sample-label combination with low cell count (e.g. 20)
summed.filt <- summed[,summed$ncells >= 20]
# save breakdown of pseudobulk datasets
check <- as.data.frame(colData(summed.filt))
check_table <- as.data.frame(table(check$CellType,check$Sex)) %>%
  pivot_wider(names_from = Var2, values_from = Freq)
colnames(check_table) <- c("CellType","Female_n","Male_n")
write.csv(check_table,"../../../data/TableS2_DE_TST_NumberPseudobulksBySex.csv", row.names = F)

# convert character to factor columns
summed.filt$Sex <- factor(summed.filt$Sex)
summed.filt$Sample <- factor(summed.filt$Sample)

## Step 3: run the differential expression analysis
# set significance thresholds
p.thr <- 0.05
padj.thr <- 0.05
lfc.thr <- 0

# using default settings for pseudoBulkDGE (robust=TRUE)
de.results <- pseudoBulkDGE(summed.filt, 
                            label = summed.filt$CellType, # specify cluster column
                            design = ~Sex, 
                            coef = "SexMale", # this means pos.FC associated with male
                            condition = summed.filt$Sex, # specify experimental condition column; only used for abundance-based gene filtering
                            lfc = lfc.thr, # specify the log fold change filter here
                            sorted = TRUE,
                            robust = TRUE
)

# filter the output based on significance thresholds
de.filt.padj <- lapply(de.results, function(x) subset(x, FDR < padj.thr))

# add a column to indicate in which group genes are expressed higher
de.filt.padj <- lapply(de.filt.padj, function(x) 
  cbind(x, Sex_marker = ifelse(x$logFC > 0, "up in male", "up in female")))

for (i in 1:length(de.filt.padj)){
  de.filt.padj[[i]]$cluster <- rep(names(de.filt.padj[i]),times = nrow(de.filt.padj[[i]]))
}
de.filt.padj <- do.call(rbind, de.filt.padj)
de.filt.padj$gene <- row.names(de.filt.padj)
de.filt.padj <- as.data.frame(de.filt.padj)
results <- de.filt.padj %>% left_join(chr.annot) %>% dplyr::select(cluster,gene,Chromosome,Sex_marker,logFC,logCPM,F,PValue,FDR)

## Step 4: Add to Source Data
wb <- loadWorkbook("../../../SourceData.xlsx")

addWorksheet(wb,"Fig3B")
writeData(wb, "Fig3B", results)

saveWorkbook(wb, "../../../SourceData.xlsx", overwrite = TRUE)
