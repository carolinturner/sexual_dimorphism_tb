## Script 5: Differential abundance analysis - Latent TB: Day 2 TST (single cell)

library(HDF5Array)
library(scater)
library(edgeR)
library(tidyverse)

## load sce object 
sce <- loadHDF5SummarizedExperiment("../../../data/Processed_sce_TSTD2_LatentTB_h5")


#### CellType ####
# Quantify number of cells assigned to each cluster
abundances <- table(sce$CellType, sce$Sample) 
abundances <- unclass(abundances) 
head(abundances)

# Attaching some column metadata.
extra.info <- colData(sce)[match(colnames(abundances), sce$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)
# convert character to factor columns
y.ab[["samples"]][["Sex"]] <- factor(y.ab[["samples"]][["Sex"]])
y.ab[["samples"]][["Sample"]] <- factor(y.ab[["samples"]][["Sample"]])

# Filter out low-abundance labels (= very rare subpopulations that contain only a handful of cells).
# Most clusters will not be of low-abundance (otherwise there would not have been enough evidence to define the cluster in the first place).
keep <- filterByExpr(y.ab, group=y.ab$samples$Sex) # erythrocytes removed 
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~Sex, y.ab$samples)
design # display design matrix to inform value for coef below

y.ab <- estimateDisp(y.ab, design, trend="none") # turn off trend as not enough points for its stable estimation.
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# test for differences in abundance between sample groups
res <- glmQLFTest(fit.ab, coef="SexMale") # coef as defined in design matrix; this is the group for which fold change is calculated
summary(decideTests(res))

# display results for all cell types
# the log-fold change here refers to the change in cell abundance between sample groups
topTags(res, n=13) # coefficient is displayed at top of table, this is the group for which fold change is calculated

# save results to dataframe
df <- as.data.frame(topTags(res, n=12))
df <- tibble::rownames_to_column(df, "CellType")

# add a column to indicate which group contributes more to a cluster
df$FC_direction <- ifelse(df$logFC > 0, "up in male", "up in female")

# write to file
write.csv(df, "../../../data/TableS3_DA_TST_CellType.csv", row.names = FALSE)

#### Cluster ####
# Quantify number of cells assigned to each cluster
abundances <- table(sce$label.3, sce$Sample) 
abundances <- unclass(abundances) 
head(abundances)

# Attaching some column metadata.
extra.info <- colData(sce)[match(colnames(abundances), sce$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)
# convert character to factor columns
y.ab[["samples"]][["Sex"]] <- factor(y.ab[["samples"]][["Sex"]])
y.ab[["samples"]][["Sample"]] <- factor(y.ab[["samples"]][["Sample"]])

# Filter out low-abundance labels (= very rare subpopulations that contain only a handful of cells).
# Most clusters will not be of low-abundance (otherwise there would not have been enough evidence to define the cluster in the first place).
keep <- filterByExpr(y.ab, group=y.ab$samples$Sex) 
y.ab <- y.ab[keep,]
summary(keep) # 69 of 101 clusters retained

design <- model.matrix(~Sex, y.ab$samples)
design # display design matrix to inform value for coef below

y.ab <- estimateDisp(y.ab, design, trend="none") # turn off trend as not enough points for its stable estimation.
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# test for differences in abundance between sample groups
res <- glmQLFTest(fit.ab, coef="SexMale") # coef as defined in design matrix; this is the group for which fold change is calculated
summary(decideTests(res))

# display results for all 69 clusters
# the log-fold change here refers to the change in cell abundance between sample groups
topTags(res, n=101) # coefficient is displayed at top of table, this is the group for which fold change is calculated

# save results to dataframe
df <- as.data.frame(topTags(res, n=101))
df <- tibble::rownames_to_column(df, "label.3")

# add a column to indicate which group contributes more to a cluster
df$FC_direction <- ifelse(df$logFC > 0, "up in male", "up in female")

# add annotations
annot <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::select(label.3,CellType) %>%
  unique()
df <- left_join(df, annot)
colnames(df)[1] <- "Cluster"
colnames(df)[8] <- "Annotation"

# write to file
write.csv(df, "../../../data/TableS3_DA_TST_Cluster.csv", row.names = FALSE)
