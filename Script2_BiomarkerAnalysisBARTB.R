## Script 2: BAR-TB biomarker analyses

library(tidyverse)
library(ggpubr)

# load data
meta <- read.csv("../../../../BAR-TB/Re-analysis_Feb24/BARTB_pheno_updated_Jan24.csv") %>%
  dplyr::select(UCL_BAR,tb,gender) %>%
  dplyr::rename(sample = "UCL_BAR")
scores <- read.csv("../../../../BAR-TB/Re-analysis_Feb24/BARTB_master_SVA_newTBdef_Apr24.csv") %>%
  dplyr::rename(sample = "X")

# select signatures
sig <- c("Sweeney3","Roe3","Kaforou25","BATF2")
dat <- scores %>% select(sample,all_of(sig))

# select samples
dat1 <- dat %>% left_join(meta) %>% filter(tb == "TB")
dat_long <- dat1 %>% pivot_longer(cols = all_of(sig), names_to = "Module")

write.csv(dat_long,"../../../data/SourceData_BARTB_BiomarkerAnalysis.csv", row.names = F)
