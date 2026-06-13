## Script 5: 10X single cell TST data - differential abundance

library(scater)
library(scran)
library(tidyverse)
library(edgeR)
library(ggh4x)

## load sce object and metadata
sce <- readRDS("../../../../TST_blisters/NEW_n=31/sce_post-QC_norm_dimred_int_clust3_VDJ_annot.rds")
meta <- read.csv("../../../../TST_blisters/NEW_n=31/metadata.csv") %>% dplyr::select(Sample,Gender)

## add sex to sce object
colData(sce) <- colData(sce) %>%
  as.data.frame() %>%
  left_join(meta) %>%
  DataFrame

#### CellType ####
# Quantify number of cells assigned to each cluster
abundances <- table(sce$CellType, sce$Sample) 
abundances <- unclass(abundances) 
head(abundances)

# Attaching some column metadata.
extra.info <- colData(sce)[match(colnames(abundances), sce$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)
# convert character to factor columns
y.ab[["samples"]][["Gender"]] <- factor(y.ab[["samples"]][["Gender"]])
y.ab[["samples"]][["Sample"]] <- factor(y.ab[["samples"]][["Sample"]])

# Filter out low-abundance labels (= very rare subpopulations that contain only a handful of cells).
# Most clusters will not be of low-abundance (otherwise there would not have been enough evidence to define the cluster in the first place).
keep <- filterByExpr(y.ab, group=y.ab$samples$Gender) # erythrocytes removed 
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~Gender, y.ab$samples)
design # display design matrix to inform value for coef below

y.ab <- estimateDisp(y.ab, design, trend="none") # turn off trend as not enough points for its stable estimation.
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# test for differences in abundance between sample groups
res <- glmQLFTest(fit.ab, coef="GenderMale") # coef as defined in design matrix; this is the group for which fold change is calculated
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
write.csv(df, "../../../data/DA_edgeR_TST_CellType.csv", row.names = FALSE)

#### Cluster ####
# Quantify number of cells assigned to each cluster
abundances <- table(sce$label.3, sce$Sample) 
abundances <- unclass(abundances) 
head(abundances)

# Attaching some column metadata.
extra.info <- colData(sce)[match(colnames(abundances), sce$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)
# convert character to factor columns
y.ab[["samples"]][["Gender"]] <- factor(y.ab[["samples"]][["Gender"]])
y.ab[["samples"]][["Sample"]] <- factor(y.ab[["samples"]][["Sample"]])

# Filter out low-abundance labels (= very rare subpopulations that contain only a handful of cells).
# Most clusters will not be of low-abundance (otherwise there would not have been enough evidence to define the cluster in the first place).
keep <- filterByExpr(y.ab, group=y.ab$samples$Gender) 
y.ab <- y.ab[keep,]
summary(keep) # 69 of 101 clusters retained

design <- model.matrix(~Gender, y.ab$samples)
design # display design matrix to inform value for coef below

y.ab <- estimateDisp(y.ab, design, trend="none") # turn off trend as not enough points for its stable estimation.
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)

# test for differences in abundance between sample groups
res <- glmQLFTest(fit.ab, coef="GenderMale") # coef as defined in design matrix; this is the group for which fold change is calculated
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
annot <- read.csv("../../../../TST_blisters/NEW_n=31/cluster_order_ont_fun.csv") %>%
  select(label.3,CellType)
df <- left_join(df, annot)
colnames(df)[1] <- "Cluster"

# write to file
write.csv(df, "../../../data/DA_edgeR_TST_cluster.csv", row.names = FALSE)


