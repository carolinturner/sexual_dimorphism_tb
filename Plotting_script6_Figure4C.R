library(pheatmap)
library(tidyverse)

# load data
dat <- read.csv("../../../data/SourceData_BAL_DEbyCellType_PseudobulkFilter20.csv") %>%
  select(logFC,cluster,gene,Chromosome)

# reformat data
d <- dat %>%
  pivot_wider(names_from = cluster, values_from = logFC) %>%
  column_to_rownames("gene") %>%
  mutate(Chromosome = ifelse(Chromosome %in% c("X","Y"),Chromosome,"autosomal")) %>%
  arrange(desc(Chromosome))
annot <- d %>% select(Chromosome)
mat <- d %>% select(-Chromosome) %>%
  as.matrix()

# define annotation colors
mycolors <- list(
  Chromosome = c(autosomal="grey30",
                 X = "lightcoral",
                 Y = "turquoise")
)

# define heatmap colour range to center symmetrically around 0
rg <- max(abs(mat), na.rm = TRUE) # absolute maximum fold change
breaks <- seq(-rg, rg, length.out = 100) # 100 breakpoints 
my_palette <- colorRampPalette(c("red", "white", "blue"))(length(breaks) - 1) # 99 colors


p <- pheatmap(mat,
              cluster_rows = F,
              cluster_cols = F,
              scale = "none",
              color = my_palette,
              breaks = seq(-rg, rg, length.out = 100),
              show_rownames = T,
              annotation_row = annot,
              annotation_colors = mycolors)
ggsave("../../../figures/Figure4C.svg",p,width = 15,height = 12,units = "cm",dpi = 300)
