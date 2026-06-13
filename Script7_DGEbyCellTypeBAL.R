## Script 7: BAL single cell data - differential gene expression

library(tidyverse)
library(AnnotationHub)
library(ensembldb)
library(scater)
library(scran)

# load data
sce <- readRDS("GSE326212_TB.rds")
coldat_df <- colData(sce) %>% as.data.frame()

# Step 1: get chromosome annotations
ens.hs <- query(AnnotationHub(),c("Homo sapiens","EnsDb", 98))[[1]]

chr_annot <- AnnotationDbi::select(ens.hs,
                                   keys = rownames(sce), 
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL","SEQNAME")) %>%
  set_names(c("gene","Chromosome"))

# Step 2: create pseudobulk matrices per sample/CellType combination
# raw counts are summed
summed <- aggregateAcrossCells(sce,
                               id = DataFrame(
                                 label = sce$celltypes,
                                 sample = sce$sample
                               ))
# filter out sample-label combinations with low cell count (e.g. 20)
summed.filt <- summed[,summed$ncells >=20]
# save breakdown of pseudobulk datasets
check <- as.data.frame(colData(summed.filt))
check_table <- as.data.frame(table(check$celltypes,check$sex)) %>%
  pivot_wider(names_from = Var2, values_from = Freq)
colnames(check_table) <- c("CellType","Female_n","Male_n")
write.csv(check_table,"../../../data/DE_BAL_NumberPseudobulksBySex.csv", row.names = F)

# convert character to factor columns
summed.filt$sex <- factor(summed.filt$sex)
summed.filt$sample <- factor(summed.filt$sample)

# Step 3: run the differential expression analysis
# set significance thresholds
p.thr <- 0.05
padj.thr <- 0.05
lfc.thr <- 0

# using default settings for pseudoBulkDGE (robust = TRUE)
de.results <- pseudoBulkDGE(summed.filt,
                            label = summed.filt$celltypes, # specify cluster column
                            design = ~sex,
                            coef = "sexM", # this means pos. FC associated with male
                            condition = summed.filt$sex, # specify experimental condition column; only used for abundance-based gene filtering
                            lfc = lfc.thr, # specify the log fold change filter here
                            sorted = TRUE,
                            robust = TRUE)

# filter the output based on significance thresholds
de.filt.padj <- lapply(de.results, function(x) subset(x, FDR < padj.thr))

# add a column to indicate in which group genes are expressed higher
de.filt.padj <- lapply(de.filt.padj, function(x)
  cbind(x, sex_marker = ifelse(x$logFC > 0, "up in male", "up in female")))

for (i in 1:length(de.filt.padj)){
  de.filt.padj[[i]]$cluster <- rep(names(de.filt.padj[i]),times =nrow(de.filt.padj[[i]]))
}
de.filt.padj <- do.call(rbind, de.filt.padj)
de.filt.padj$gene <- row.names(de.filt.padj)
df <- as.data.frame(de.filt.padj)
df <- left_join(df,chr_annot)

write.csv(df, "../../../data/SourceData_BAL_DEbyCellType_PseudobulkFilter20.csv", row.names = F)

