## Script 3: Module analyses Latent TB TST (D2+D7)

library(openxlsx)
library(tidyverse)
library(HGNChelper)

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
  filter(!module %in% c("STAT1_Blood","TNF_Blood"))

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
dat <- read.csv("../../../data/tpm_PC0.001_log2_genesymbol_dedup.csv")

# correct relevant old gene symbols
data_genes_old <- dat %>% pull(external_gene_name)
updates <- checkGeneSymbols(data_genes_old,
                            unmapped.as.na = FALSE,
                            map = currentHumanMap)
updates_relevant <- updates %>%
  filter(Suggested.Symbol %in% gene_updates) %>%
  dplyr::select(-Approved) %>%
  dplyr::rename(external_gene_name =x, updated_symbol = Suggested.Symbol)

dat_clean <- dat %>%
  filter(external_gene_name %in% updates_relevant$external_gene_name) %>%
  left_join(updates_relevant) %>%
  dplyr::select(-external_gene_name) %>%
  column_to_rownames("updated_symbol")

## Step 2: check module genes missing from dataset ####
for (mod in names(module_genes)) {
  genes_in_module <- module_genes[[mod]]
  missing_genes <- setdiff(genes_in_module, rownames(dat_clean))
  
  if (length(missing_genes) > 0) {
    cat("\nMissing genes in module", mod, ":\n")
    print(missing_genes)
  } else {
    cat("\nAll genes in module", mod, "are present in the data.\n")
  }
}

## Step 3: calculate module score per sample and add sex metadata ####
summary <- data.frame()

for (mod in names(module_genes)){
  genes <- module_genes[[mod]]
  score <- dat_clean %>%
    filter(row.names(dat_clean) %in% genes) %>%
    summarise(across(all_of(names(dat_clean)),mean))
  rownames(score) <- mod
  summary <- rbind(summary,score)
}

summary.t <- as.data.frame(t(summary)) %>% rownames_to_column("sample")
meta <- read.table("../../../data/E-MTAB-14687.sdrf.txt",header=T,sep="\t") %>%
  dplyr::select(Source.Name,Characteristics.stimulus.,Characteristics.sex.)%>%
  dplyr::rename(sample = "Source.Name",
                stimulant = "Characteristics.stimulus.",
                sex = "Characteristics.sex.")

results <- summary.t %>% left_join(meta)

## Step 4: Save source data file ####
wb <- loadWorkbook("../../../SourceData.xlsx")

addWorksheet(wb,"Fig2AB,S2")
writeData(wb, "Fig2AB,S2", results)

saveWorkbook(wb, "../../../SourceData.xlsx", overwrite = TRUE)
