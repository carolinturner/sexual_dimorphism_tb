## Table 1: Signature performance Subclinical TB Blood

library(openxlsx)
library(tidyverse)
library(pROC)
library(flextable)

## Load data
d1 <- read.xlsx("../../../SourceData.xlsx",sheet="Table1")

## Stats
# number of patients
length(unique(d1$study_id)) 
# sex breakdown all
d1 %>% group_by(sex) %>% summarise(n_participants = n_distinct(study_id))
# sex breakdown subclinical tb
d1 %>% filter(subclinical_tb0to12 == 1) %>% group_by(sex) %>% summarise(n_patients = n_distinct(study_id)) 

## Define functions
# calculate AUCs
auc <- function(gene, dataset, outcome, by_sex){ 
  dataset <- dataset %>% 
    filter(sex == by_sex) %>%
    filter(!is.na(.data[[gene]]) & !is.na(.data[[outcome]]))
  
  # ROC and AUC
  roc1 <- roc(dataset[[outcome]] ~ dataset[[gene]], quiet = TRUE)
  auc_ci <- ci.auc(roc1)
  
  auc_t <- data.frame(
    gene = gene,
    sex = by_sex,
    auc.ci = paste0(round(auc_ci[2], 2), " (", round(auc_ci[1], 2), " - ", round(auc_ci[3], 2), ")")
  )
  return(auc_t)
}

# calculate p-values
sex_comparison_p <- function(gene, outcome, dataset) {
  male_data <- dataset %>% 
    filter(sex == "Male", !is.na(.data[[gene]]), !is.na(.data[[outcome]]))
  female_data <- dataset %>% 
    filter(sex == "Female", !is.na(.data[[gene]]), !is.na(.data[[outcome]]))
  roc_male   <- roc(male_data[[outcome]], male_data[[gene]], quiet = TRUE)
  roc_female <- roc(female_data[[outcome]], female_data[[gene]], quiet = TRUE)
  roc_test <- roc.test(roc_male, roc_female, method = "d", paired = FALSE)
  return(roc_test$p.value)
}

## Run AUC and p-value functions for top signatures
auc_results <- data.frame()
for (i in top_signatures) {
  for (k in c("Female", "Male")) {
    temp <- auc(gene=i, outcome="subclinical_tb0to12", dataset=combined_data, by_sex = k)
    auc_results <- rbind(auc_results, temp)
  }
}

sex_p_val <- data.frame(
  gene = top_signatures,
  p_sex_comparison = sapply(top_signatures, function(x) {
    sex_comparison_p(x, "subclinical_tb0to12", combined_data)
  })
) %>% 
  mutate(adjusted_p_sex = p.adjust(p_sex_comparison, method = "fdr"))

## Export summary table
summary <- auc_results %>%
  pivot_wider(names_from = sex, values_from = auc.ci, names_prefix = "AUC ") %>%
  left_join(sex_p_val, by = "gene") %>%
  select(
    Signature = gene,
    `AUC Female` = `AUC Female`,
    `AUC Male` = `AUC Male`,
    `Adjusted p-values` = adjusted_p_sex
  ) %>%
  mutate(`Adjusted p-values` = round(`Adjusted p-values`, 3))

# export
table1 <- flextable(summary) %>%
  set_header_labels(
    Signature = "Gene Signature",
    `AUC Female` = "AUC Female (95% CI)",
    `AUC Male` = "AUC Male (95% CI)",
    `adjusted p-values` = "Adjusted p-value*"
  ) %>%
  autofit() %>%
  theme_vanilla() 

# save as Word Document
save_as_docx(table1, path = "../../../figures/Table1.docx")
