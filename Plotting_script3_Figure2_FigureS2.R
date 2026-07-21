library(openxlsx)
library(tidyverse)
library(ggpubr)
library(biomaRt)
library(ggrepel)
library(ggbreak)
library(patchwork)

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

### Figure 2A ####
d1 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig2AB,S2") %>%
  pivot_longer(cols=-c(sample,stimulant,sex),names_to = "module_name",values_to = "moduleTPM") %>%
  filter(module_name %in% c("IFNA_KC","IFNG_KC","TNF_KC")) 

# calculate Z-scores, using saline samples as controls
ctrl <- d1 %>%
  filter(stimulant == "saline") %>%
  group_by(module_name) %>%
  summarise(mean.ctrl = mean(moduleTPM),
            sd.ctrl = sd(moduleTPM),
            .groups = "drop")
zdat <- d1 %>%
  filter(!stimulant == "saline") %>%
  left_join(ctrl) %>%
  mutate(stimulant = recode(stimulant,
                            "TST_D2" = "Day 2 TST",
                            "TST_D7" = "Day 7 TST"),
         module_Zscore = (moduleTPM - mean.ctrl)/sd.ctrl)

fig2a <- ggplot(data = zdat, aes(x = sex, y = module_Zscore, fill = sex)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(
    comparisons = list(c("Female", "Male")),
    label = "p.signif",
    size = 2.5,
    bracket.size = 0.2
  ) + 
  scale_y_continuous(limits = c(-3, 15), expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = c("Female" = "F", "Male" = "M")) +
  facet_grid(stimulant ~ module_name) +
  xlab("") +
  ylab("Module Z-score (log2)") +
  My_Theme +
  theme(legend.position = "none")

### Figure 2B ####
d2 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig2AB,S2") %>%
  pivot_longer(cols=-c(sample,stimulant,sex),names_to = "module_name",values_to = "moduleTPM") %>%
  filter(module_name %in% c("IFNA2_TST","IFNG_TST","TNF_TST")) 

# calculate Z-scores, using saline samples as controls
ctrl <- d2 %>%
  filter(stimulant == "saline") %>%
  group_by(module_name) %>%
  summarise(mean.ctrl = mean(moduleTPM),
            sd.ctrl = sd(moduleTPM),
            .groups = "drop")
zdat <- d2 %>%
  filter(!stimulant == "saline") %>%
  left_join(ctrl) %>%
  mutate(stimulant = recode(stimulant,
                            "TST_D2" = "Day 2 TST",
                            "TST_D7" = "Day 7 TST"),
         module_Zscore = (moduleTPM - mean.ctrl)/sd.ctrl)

fig2b <- ggplot(data = zdat, aes(x = sex, y = module_Zscore, fill = sex)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(
    comparisons = list(c("Female", "Male")),
    label = "p.signif",
    size = 2.5,
    bracket.size = 0.2
  ) + 
  scale_y_continuous(limits = c(-3, 15), expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = c("Female" = "F", "Male" = "M")) +
  facet_grid(stimulant ~ module_name) +
  xlab("") +
  ylab("Module Z-score (log2)") +
  My_Theme +
  theme(legend.position = "none")

### Figure 2C ####
d3 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig2C") %>%
  mutate(diffexpressed = case_when(
    log2FoldChange >= 1 & padj < 0.05 ~ "Enriched in Male",
    log2FoldChange <= -1 & padj < 0.05 ~ "Enriched in Female",
    TRUE ~ "No differential expression"))

# biomaRt to create label with gene name and chromosome
mart <- biomaRt::useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "https://may2025.archive.ensembl.org"
)

mapping <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = d3$ensembl_gene_id,
  mart = mart
)

d3_annot <- d3 %>%
  left_join(mapping) %>% 
  mutate(label = if_else(diffexpressed != "No differential expression" & !is.na(external_gene_name),
                         paste0(external_gene_name, " (", chromosome_name, ")"),
                         NA_character_),
         panel_label = "Day 2 TST")

# make a custom colour vector
mycolour <- c("Enriched in Female" = "#F8766D",
              "Enriched in Male" = "#00BFC4",
              "No differential expression" = "grey50")


# make "volcano plot" (= scatter plot with log2 fold change on x-axis and adj p value on y-axis)
# shared base 
text_theme <- element_text(face = "bold", size = t, colour =tc)

base <- ggplot(d3_annot, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = c(-1, 1), col = "grey") +
  geom_hline(yintercept = -log10(0.05), col = "grey") +
  scale_x_continuous(limits = c(-2, 2), breaks = c(-1, 0, 1)) +
  scale_color_manual(values = mycolour) +
  labs(col = "Genes") +
  theme(plot.background  = element_blank(),
        axis.title.x     = text_theme,
        axis.title.y     = element_blank(),
        axis.text.x      = text_theme,
        axis.text.y      = text_theme,
        legend.title     = text_theme,
        legend.text      = text_theme,
        panel.background = element_rect(fill = "grey95"),
        panel.border     = element_rect(fill = NA, linewidth = 0.5, colour = tc))

# top panel
top <- base +
  facet_wrap(~ panel_label) +
  coord_cartesian(ylim = c(50, 55)) +
  scale_y_continuous(breaks = 50) +
  theme(strip.text       = text_theme,
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.5),
        axis.title.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        legend.position  = "none")

# bottom panel
bottom <- base +
  geom_text_repel(data = subset(d3_annot, !is.na(label) & -log10(padj) < 14),
                  aes(label = label), colour = "black", size = 3,
                  max.overlaps = 30, box.padding = 0.3,
                  segment.size = 0.2, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, 14)) +
  scale_y_continuous(breaks = c(0, 4, 8, 12))

# assemble + single shared y-axis label 
ylab <- wrap_elements(grid::textGrob("-log10(padj)", rot = 90,
                                     gp = grid::gpar(fontface = "bold", fontsize = 8)))

# fig2c 
stack_c <- (top / bottom) + plot_layout(heights = c(1, 5))
fig2c <- wrap_plots(ylab, stack_c, widths = c(0.04, 1), guides = "collect") & theme(legend.position = "bottom")

### Figure 2D ####
d4 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig2D") %>%
  mutate(diffexpressed = case_when(
    log2FoldChange >= 1 & padj < 0.05 ~ "Enriched in Male",
    log2FoldChange <= -1 & padj < 0.05 ~ "Enriched in Female",
    TRUE ~ "No differential expression"))

d4_annot <- d4 %>%
  left_join(mapping) %>% 
  mutate(label = if_else(diffexpressed != "No differential expression" & !is.na(external_gene_name),
                         paste0(external_gene_name, " (", chromosome_name, ")"),
                         NA_character_),
         panel_label = "Day 7 TST")

# make "volcano plot" (= scatter plot with log2 fold change on x-axis and adj p value on y-axis)
# shared base 
base <- ggplot(d4_annot, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = c(-1, 1), col = "grey") +
  geom_hline(yintercept = -log10(0.05), col = "grey") +
  scale_x_continuous(limits = c(-2, 2), breaks = c(-1, 0, 1)) +
  scale_color_manual(values = mycolour) +
  labs(col = "Genes") +
  theme(plot.background  = element_blank(),
        axis.title.x     = text_theme,
        axis.title.y     = element_blank(),
        axis.text.x      = text_theme,
        axis.text.y      = text_theme,
        legend.title     = text_theme,
        legend.text      = text_theme,
        panel.background = element_rect(fill = "grey95"),
        panel.border     = element_rect(fill = NA, linewidth = 0.5, colour = tc))

# top panel
top <- base +
  facet_wrap(~ panel_label) +
  coord_cartesian(ylim = c(48, 50)) +
  scale_y_continuous(breaks = 48) +
  theme(strip.text       = text_theme,
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.5),
        axis.title.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        legend.position  = "none")

# bottom panel
bottom <- base +
  geom_text_repel(data = subset(d4_annot, !is.na(label) & -log10(padj) < 14),
                  aes(label = label), colour = "black", size = 3,
                  max.overlaps = 30, box.padding = 0.3,
                  segment.size = 0.2, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(breaks = c(0, 4, 8))

# assemble + single shared y-axis label 
ylab <- wrap_elements(grid::textGrob("-log10(padj)", rot = 90,
                                     gp = grid::gpar(fontface = "bold", fontsize = 8)))

# fig2d 
stack_d <- (top / bottom) + plot_layout(heights = c(1, 5))
fig2d <- wrap_plots(ylab, stack_d, widths = c(0.04, 1), guides = "collect") & theme(legend.position = "bottom")

#### Combine Figure 2 ####
combined_fig2 <-
  wrap_elements(fig2a) +
  wrap_elements(fig2b) +
  wrap_elements(fig2c) +
  wrap_elements(fig2d) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

combined_fig2


# save svg
ggsave("../../../figures/Figure2.svg",
       plot = combined_fig2,
       unit = "cm",
       width = 26,
       height = 24,
       dpi = 300)


### Figure S2 ####
d5 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig2AB,S2") %>%
  pivot_longer(cols=-c(sample,stimulant,sex),names_to = "module_name",values_to = "moduleTPM") %>%
  filter(module_name %in% c("IFN1_MDM","IFN2_MDM","TNF_MDM")) 

# calculate Z-scores, using saline samples as controls
ctrl <- d5 %>%
  filter(stimulant == "saline") %>%
  group_by(module_name) %>%
  summarise(mean.ctrl = mean(moduleTPM),
            sd.ctrl = sd(moduleTPM),
            .groups = "drop")
zdat <- d5 %>%
  filter(!stimulant == "saline") %>%
  left_join(ctrl) %>%
  mutate(stimulant = recode(stimulant,
                            "TST_D2" = "Day 2 TST",
                            "TST_D7" = "Day 7 TST"),
         module_Zscore = (moduleTPM - mean.ctrl)/sd.ctrl)

figS2 <- ggplot(data = zdat, aes(x = sex, y = module_Zscore, fill = sex)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(
    comparisons = list(c("Female", "Male")),
    label = "p.signif",
    size = 2.5,
    bracket.size = 0.2
  ) + 
  scale_y_continuous(limits = c(-3, 15), expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = c("Female" = "F", "Male" = "M")) +
  facet_grid(stimulant ~ module_name) +
  xlab("") +
  ylab("Module Z-score (log2)") +
  My_Theme +
  theme(legend.position = "none")

# save svg
ggsave("../../../figures/FigureS2.svg", plot = figS2, unit = "cm", width = 20, height = 16, dpi = 300)
