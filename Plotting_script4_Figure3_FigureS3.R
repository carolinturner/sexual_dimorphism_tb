library(openxlsx)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(patchwork)
library(pheatmap)
library(gtable)

#My_Theme
t = 8 #size of text
m = 4 #size of margin around text
tc = "black" #colour of text
My_Theme = theme(
  axis.title.x = element_text(size = t, face = "bold", margin = margin(t = m)),
  axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = 0, hjust = 0.5),
  axis.title.y.left = element_text(size = t, face = "bold", margin = margin(r = 10)),
  axis.title.y.right = element_text(size = t, face = "bold", margin = margin(l = m)),
  axis.text.y = element_text(size = t, face = "bold", colour = tc),
  legend.title = element_text(size=t, face = "bold", colour = tc),
  legend.text = element_text(size=t, face = "bold", colour = tc),
  plot.title = element_text(size=t, face = "bold", colour = tc),
  strip.text = element_text(size=t, face = "bold", colour = tc),
  strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.5),
  panel.border = element_rect(fill = NA, linewidth = 0.5, colour = tc),
  panel.background = element_rect(fill = "grey95"),
  legend.position = "right", legend.justification = "top",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0,0,0,0),
  legend.box.spacing = unit(c(0,0,0,0),"cm"),
  panel.spacing = unit(0.05,"cm"),
  plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm")
)

# Annotation colors
mycolors <- list(
  Chromosome = c(autosomal="grey30",
                 X = "#F8766D",
                 Y = "#00BFC4")
)

### load module data for Figure 3A and S3 ####
d1 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig3A,S3") %>%
  filter(!CellType %in% c("undefined","erythrocytes")) %>%
  mutate(Sex = recode(Sex,
                      "Female" = "F",
                      "Male" = "M"),
         CellType = recode(CellType,
                           CD4_T = "CD4",
                           CD8_T = "CD8",
                           atypical_CD8 = "CD8_atyp.",
                           gd_T_V1 = "gd_V1",
                           gd_T_V2 = "gd_V2",
                           keratinocytes_basal = "KC_b",
                           keratinocytes_suprabasal = "KC_sb",
                           melanocytes = "MC",
                           neutrophils = "neutro"))
d1$CellType <- factor(d1$CellType, levels = c("CD4","CD8","CD8_atyp.","gd_V1","gd_V2","NKT","NK","myeloid","neutro","KC_b","KC_sb","MC"))

# subset data to keep only selected modules
dat1 <- subset(d1, module %in% c("IFNA2_TST","IFNG_TST","TNF_TST")) # TST modules
dat2 <- subset(d1, module %in% c("IFNA_KC","IFNG_KC","TNF_KC")) # KC modules
dat3 <- subset(d1, module %in% c("IFN1_MDM","IFN2_MDM","TNF_MDM")) # MDM modules

# stats (multiple testing correction)
stats1 <- dat1 %>%
  group_by(module, CellType) %>%
  wilcox_test(Zscore_average ~ Sex) %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup() %>%
  mutate(y.position = max(dat1$Zscore_average, na.rm = TRUE) + 0.5,
         p.adj.label = ifelse(p.adj > 0.05, "ns", p.adj))
stats2 <- dat2 %>%
  group_by(module, CellType) %>%
  wilcox_test(Zscore_average ~ Sex) %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup() %>%
  mutate(y.position = max(dat2$Zscore_average, na.rm = TRUE) + 0.5,
         p.adj.label = ifelse(p.adj > 0.05, "ns", p.adj))
stats3 <- dat3 %>%
  group_by(module, CellType) %>%
  wilcox_test(Zscore_average ~ Sex) %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup() %>%
  mutate(y.position = max(dat3$Zscore_average, na.rm = TRUE) + 0.5,
         p.adj.label = ifelse(p.adj > 0.05, "ns", p.adj))

### load differential gene expression data for Figure 3B ####
d2 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig3B") %>%
  select(logFC,cluster,gene,Chromosome)

# reformat data
dat <- d2 %>%
  pivot_wider(names_from = cluster, values_from = logFC) %>%
  column_to_rownames("gene") %>%
  mutate(Chromosome = ifelse(Chromosome %in% c("X","Y"),Chromosome,"autosomal")) %>%
  arrange(desc(Chromosome)) %>%
  dplyr::rename(CD4 = "CD4_T",
                CD8 = "CD8_T",
                "CD8_atyp." = "atypical_CD8",
                gd_V1 = "gd_T_V1",
                gd_V2 = "gd_T_V2",
                KC_sb = "keratinocytes_suprabasal")
annot <- dat %>% select(Chromosome)
mat <- dat %>% select(-Chromosome) %>% select(CD4,CD8,CD8_atyp.,gd_V1,gd_V2,NKT,NK,myeloid,KC_sb) %>% as.matrix()

### Figure 3A ####
fig3a <- ggplot(data = dat1, aes(x=Sex, y=Zscore_average, fill=Sex)) + 
  geom_jitter(size = 0.6, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_grid(module~CellType) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_pvalue_manual(stats1, label = "p.adj.label", tip.length = 0.01) + 
  ylab("Module Z-score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())

### Figure 3B ####
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
              annotation_colors = mycolors,
              fontsize = 8
)

### Combine Figure 3 ####
p_padded <- gtable_add_padding(p$gtable, padding = unit(c(1, 0, 1, 1.5), "cm"))

bottom_row <- (wrap_elements(full = p_padded) | plot_spacer()) +
  plot_layout(widths = c(3, 1))

combined_fig3 <- fig3a / bottom_row +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

# save svg
ggsave("../../../figures/Figure3.svg",
       plot = combined_fig3,
       unit = "cm",
       width = 26,
       height = 35,
       dpi = 300)

### Supplementary Figure 3 ####
figS3a <- ggplot(data = dat2, aes(x=Sex, y=Zscore_average, fill=Sex)) + 
  geom_jitter(size = 0.6, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA) +
  facet_grid(module~CellType) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_pvalue_manual(stats2, label = "p.adj.label", tip.length = 0.01) + 
  ylab("Module Z-score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())

figS3b <- ggplot(data = dat3, aes(x=Sex, y=Zscore_average, fill=Sex)) + 
  geom_jitter(size = 0.6, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape=NA) +
  facet_grid(module~CellType) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_pvalue_manual(stats3, label = "p.adj.label", tip.length = 0.01) +
  ylab("Module Z-score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())

combined_figS3 <- (figS3a / figS3b) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

# save svg
ggsave("../../../figures/FigureS3.svg",
       plot = combined_figS3,
       unit = "cm",
       width = 26,
       height = 30,
       dpi = 300)
