## Script 3: 10X single cell TST data module analyses

library(tidyverse)
library(biomaRt)
library(HGNChelper)

## Step 0: update the gene symbol map ####
currentHumanMap <- getCurrentHumanMap()

## Step 1: load modules and update gene symbols ####
modules <- read.csv("../../../data/modules_TNF_IFN.csv") %>%
  filter(!Source.sample == "Blood" & !Module.name %in% c("IL17_KC","IFNB1_TST.URA")) %>%
  dplyr::select(Module.name,Module.genes)

# correct old gene symbols
module_genes_old <- strsplit(modules$Module.genes, ",")
names(module_genes_old) <- modules$Module.name
all_genes <- unique(unlist(module_genes_old))
updates <- checkGeneSymbols(all_genes,
                            unmapped.as.na = FALSE,
                            map = currentHumanMap)
gene_updates <- setNames(updates$Suggested.Symbol, updates$x) # named vector of corrected symbols

# apply corrections to each module
module_genes <- lapply(module_genes_old, function(g) {
  g <- trimws(g)
  g_updated <- ifelse(!is.na(gene_updates[g]), gene_updates[g], g)
  unname(g_updated)
})

## Step 2: load sce object, add metadata and update gene symbols ####
sce <- readRDS("../../../../TST_blisters/NEW_n=31/sce_post-QC_norm_dimred_int_clust3_VDJ_annot.rds")
meta <- read.csv("../../../../TST_blisters/NEW_n=31/metadata.csv") %>% dplyr::select(Sample,Gender)

# add sex to sce object
colData(sce) <- colData(sce) %>%
  as.data.frame() %>%
  left_join(meta) %>%
  DataFrame

# correct old gene symbols
data_genes <- rownames(rowData(sce))
updates <- checkGeneSymbols(data_genes,
                            unmapped.as.na = FALSE,
                            map = currentHumanMap)
all(rownames(rowData(sce)) == updates$x) # must be TRUE
rownames(rowData(sce)) <- updates$Suggested.Symbol

## Step 3: check module genes missing from dataset ####
for (mod in names(module_genes)) {
  genes_in_module <- module_genes[[mod]]
  missing_genes <- setdiff(genes_in_module, rownames(rowData(sce)))
  
  if (length(missing_genes) > 0) {
    cat("\nMissing genes in module", mod, ":\n")
    print(missing_genes)
  } else {
    cat("\nAll genes in module", mod, "are present in the data.\n")
  }
}

## load annotation
#annot <- read.csv("../../../../TST_blisters/NEW_n=31/cluster_order_ont_fun.csv")


## Step 4: calculate module scores ####
# calculate module scores
# summary data frame for module gene coverage in dataset
summary <- data.frame()

# loop through all modules, and add module scores to sce  
for (i in 1:length(module_genes)){
  module.genes <- module_genes[[i]]
  rowSubset(sce, field = names(module_genes)[i]) <- module.genes
  # determine number and percentage of module genes found in dataset
  module.size <- length(module.genes)
  detected.n <- length(which(rowData(sce)[[names(module_genes)[i]]] == TRUE))
  detected.pc <- detected.n/module.size*100
  row <- c(names(module_genes)[i], module.size, detected.n, detected.pc)
  summary <- rbind(summary,row)
  # calculate column means of module genes
  sub <- sce[which(rowData(sce)[[names(module_genes)[i]]] == TRUE),]
  means <- colMeans(logcounts(sub))
  # add column means as module score to colData to full sce object
  # make sure that the cell names are in the same order
  all(colnames(sce) == names(means))
  colData(sce)[[paste0(names(module_genes)[i],"_module-score")]] <- means
}
colnames(summary) <- c("module","module_size","module_coverage_n","module_coverage_pct")
write.csv(summary,"../../../data/ModuleCoverage_10x.csv",row.names = F)

# extract module scores and cluster label
scores <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::select("Sample","CellType",paste0(names(module_genes),"_module.score"))

# add sex 
scores <- scores %>%
  left_join(meta) %>%
  dplyr::select(Sample,Gender,CellType,everything())

# convert from wide to long
scores_long <- scores %>% rownames_to_column("cell")
scores_long <- gather(scores_long, module, score, colnames(scores_long[5:ncol(scores_long)]))

# calculate Z scores per cell
zscores <- scores_long %>% group_by(module) %>% mutate(zscore = scale(score)) %>% ungroup()
# average by group
avg <- zscores %>% group_by(module,CellType,Sample,Gender) %>% summarise(average = mean(zscore)) %>% ungroup()

# remove suffix from module names
avg$module <- gsub("_module.score","",avg$module)

write.csv(avg,"../../../data/SourceData_10X_ModuleAnalysis.csv",row.names = F)

