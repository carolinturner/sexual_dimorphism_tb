## Script 8: Differential abundance analysis - Active TB: BAL (single cell)

library(scater)
library(tidyverse)
library(openxlsx)
library(edgeR)

## Step 1: load data and add cell annotations to Source Data
sce <- readRDS("../../../data/GSE326212_TB.rds")
coldat_df <- colData(sce) %>% as.data.frame() %>% select(cell_id,sample,sex,celltypes)

wb <- loadWorkbook("../../../SourceData.xlsx")
addWorksheet(wb,"Fig4A")
writeData(wb, "Fig4A", coldat_df)
saveWorkbook(wb, "../../../SourceData.xlsx", overwrite = TRUE)

## Step 2: Differential abundance analysis
# Quantify number of cells assigned to each cluster
abundances <- table(sce$celltypes,sce$sample)
abundances <- unclass(abundances)
head(abundances)

# attaching some column metadata
extra.info <- colData(sce)[match(colnames(abundances),sce$sample),]
y.ab <- DGEList(abundances, samples=extra.info)
# convert character to factor columns
y.ab[["samples"]][["sex"]] <- factor(y.ab[["samples"]][["sex"]])
y.ab[["samples"]][["sample"]] <- factor(y.ab[["samples"]][["sample"]])

# Filter out low-abundance lables
keep <- filterByExpr(y.ab, group = y.ab$samples$sex)
summary(keep) # removes three celltypes
keep # removes 20_Migratory DC, 24_pDC_Plasma cell, 26_Mo-derived AM_TPP1
y.ab <- y.ab[keep,]

design <- model.matrix(~sex, y.ab$samples)
design # display design matrix to inform value for coef below

y.ab <- estimateDisp(y.ab, design, tren="none")
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# test for differences in abundance between groups
res <- glmQLFTest(fit.ab, coef="sexM") # coef as defined in design matrix; this is the group for wich fold change is calculated
summary(decideTests(res)) # 1 down in male, 23 not sig, 0 up

topTags(res, n=24)

# save results to dataframe
df <- as.data.frame(topTags(res, n=24))
df <- tibble::rownames_to_column(df, "CellType")

# add a column to indicate which group contributes more to a cluster
df$FC_direction <- ifelse(df$logFC >0, "up in male", "up in female")

# write to file
write.csv(df,"../../../data/TableS4_DA_BAL_CellType.csv",row.names = FALSE)
