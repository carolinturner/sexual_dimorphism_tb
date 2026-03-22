## Script 1: BAR-TB module analyses

library(tidyverse)
library(biomaRt)
library(HGNChelper)

## Step 0: update the gene symbol map ####
currentHumanMap <- getCurrentHumanMap()

## Step 1: load modules and update gene symbols ####
modules <- read.csv("../../../data/modules_TNF_IFN.csv") %>%
  dplyr::select(Module.name,Module.genes) %>%
  filter(Module.name %in% c("STAT1_Covidsortium","TNF.Blood.RachelByngMaddick",
                            "IFN1_MDM","IFN2_MDM","IFNG_MDM","TNF_MDM"))

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


## Step 2: load BARTB data, annotate and update gene symbols ####
tpm_id <- read.csv("../../../../BAR-TB/Re-analysis_Feb24/BARTB_post_sva_tb-hiv-gender-age-race_13var.csv") %>%
  dplyr::rename(ensembl_gene_id = "X")

# annotate matrix with gene annotations and remove rows without gene name annotation
ensembl <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl")

biomart <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart = ensembl
)
tpm_annot <- left_join(tpm_id,biomart) %>%
  dplyr::select(-ensembl_gene_id) %>%
  dplyr::select(external_gene_name,everything()) %>%
  mutate_if(is.character, ~na_if(., '')) %>%
  na.omit()

# remove duplicate genes (keeping duplicate with highest average expression)
tpm_old <- tpm_annot %>%
  mutate(mean = rowMeans(across(where(is.numeric)))) %>%
  arrange(desc(mean)) %>%
  distinct(external_gene_name, .keep_all = T) %>%
  dplyr::select(-mean) %>%
  column_to_rownames("external_gene_name")

# correct old gene symbols
data_genes <- row.names(tpm_old)
updates <- checkGeneSymbols(data_genes,
                            unmapped.as.na = FALSE,
                            map = currentHumanMap)
all(rownames(tpm_old) == updates$x) # must be TRUE
tpm <- tpm_old
rownames(tpm) <- updates$Suggested.Symbol


## Step 3: check module genes missing from dataset
for (mod in names(module_genes)) {
  genes_in_module <- module_genes[[mod]]
  missing_genes <- setdiff(genes_in_module, rownames(tpm_old))
  
  if (length(missing_genes) > 0) {
    cat("\nMissing genes in module", mod, ":\n")
    print(missing_genes)
  } else {
    cat("\nAll genes in module", mod, "are present in the data.\n")
  }
}

## Step 4: calculate module score per sample and append metadata
summary <- data.frame()

for (mod in names(module_genes)){
  genes <- module_genes[[mod]]
  score <- tpm %>%
    filter(row.names(tpm) %in% genes) %>%
    summarise(across(all_of(names(tpm)),mean))
  rownames(score) <- mod
  summary <- rbind(summary,score)
}

summary.t <- as.data.frame(t(summary)) %>% rownames_to_column("sample")

# append metadata and drop sample id
meta <- read.csv("../../../../BAR-TB/Re-analysis_Feb24/BARTB_pheno_updated_Jan24.csv") %>%
  dplyr::select(UCL_BAR,tb,gender) %>%
  dplyr::rename(sample = "UCL_BAR")

dat.meta <- meta %>%
  left_join(dat) %>%
  dplyr::select(-sample) %>%
  filter(tb == "TB")

# convert into long format
d <- dat.meta %>%
  pivot_longer(3:ncol(dat.meta),names_to = "module_name", values_to = "moduleTPM") %>% # CHECK ncol!!
  mutate_at(c("moduleTPM"), as.numeric)

write.csv(d,"../../../data/SourceData_BARTB_ModuleAnalysis.csv",row.names = F)
