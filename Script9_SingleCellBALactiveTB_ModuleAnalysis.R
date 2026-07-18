## Script 9: Module analyses - Active TB: BAL (single cell)

library(openxlsx)
library(tidyverse)
library(biomaRt)
library(HGNChelper)
library(HDF5Array)
library(scater)

### Step 1: update gene symbols from modules and data matrix ####
currentHumanMap <- getCurrentHumanMap()

## modules ###
modules <- read.xlsx("../../../SupplementaryTables.xlsx",
                     sheet = "Supplementary Table 1",
                     startRow = 3) %>%
  pivot_longer(cols = everything(), names_to = "module", values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  group_by(module) %>%
  summarise(genes = paste(gene, collapse = ","), .groups = "drop") %>%
  filter(module %in% c("IFN1_MDM","IFN2_MDM","TNF_MDM")) 

# correct old gene symbols
module_genes_old <- strsplit(modules$genes, ",")
names(module_genes_old) <- modules$module
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

## data ###
sce <- readRDS("../../../data/GSE326212_TB.rds")

data_genes_old <- rownames(rowData(sce))
updates <- checkGeneSymbols(data_genes_old,
                            unmapped.as.na = FALSE,
                            map = currentHumanMap)
all(rownames(rowData(sce)) == updates$x) # must be TRUE
rownames(rowData(sce)) <- updates$Suggested.Symbol

## Step 2: check module genes missing from dataset ####
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

## Step 3: calculate module scores ####
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
write.csv(summary,"../../../data/ModuleCoverage_SingleCell_BAL_ActiveTB.csv",row.names = F)

# extract module scores and cluster label
scores <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::select("sample","sex","celltypes",paste0(names(module_genes),"_module.score"))

# convert from wide to long
scores_long <- scores %>% 
  rownames_to_column("cell") %>%
  pivot_longer(cols = -c(cell,sample,sex,celltypes), names_to = "module", values_to = "score")

# calculate Z scores per cell
zscores <- scores_long %>% group_by(module) %>% mutate(zscore = scale(score)) %>% ungroup()
# average by group
avg <- zscores %>% group_by(module,celltypes,sample,sex) %>% summarise(Zscore_average = mean(zscore)) %>% ungroup()

# remove suffix from module names
avg$module <- gsub("_module.score","",avg$module)

## Step 4: add to Source Data ####
wb <- loadWorkbook("../../../SourceData.xlsx")

addWorksheet(wb,"Fig4B")
writeData(wb, "Fig4B", avg)

saveWorkbook(wb, "../../../SourceData.xlsx", overwrite = TRUE)
