library(openxlsx)
library(tidyverse)
library(ggpubr)
library(biomaRt)
library(ggrepel)
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

### Figure 1A ####
d1 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig1A,S1") %>%
  pivot_longer(cols=-c(sample,sex),names_to = "module_name",values_to = "moduleTPM") %>%
  filter(module_name %in% c("STAT1_Blood","TNF_Blood")) 

fig1a <- ggplot(data = d1, aes(x=sex, y=moduleTPM, fill=sex)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA) +
  facet_wrap(~module_name,
             scales = "free",
             ncol=2) +
  scale_x_discrete(labels = c("Female" = "F", "Male" = "M")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_compare_means(
    comparisons = list(c("Female", "Male")),
    method = "wilcox",
    label = "p.signif",
    size = 2.5,
    bracket.size = 0.2
  ) + 
  ylab("Module score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())

### Figure 1B ####
d2 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig1B") %>%
  pivot_longer(cols=-c(sample,sex),names_to = "signature",values_to = "value")

fig1b <- ggplot(data = d2, aes(x=sex, y=value, fill=sex)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA) +
  facet_wrap(~signature,
             scales = "free",
             ncol=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = c("Female" = "F", "Male" = "M")) +
  stat_compare_means(
    comparisons = list(c("Female", "Male")),
    method = "wilcox",
    label = "p.signif",
    size = 2.5,
    bracket.size = 0.2) + 
  ylab("Signature score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())

### Figure 1C ####
d3 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig1C") %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Enriched in Male",
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
                      NA_character_))

# make a custom colour vector
mycolour <- c("Enriched in Female" = "#F8766D",
              "Enriched in Male" = "#00BFC4",
              "No differential expression" = "grey50")


# make "volcano plot" (= scatter plot with log2 fold change on x-axis and adj p value on y-axis)
fig1c <- ggplot(d3_annot, aes(x=log2FoldChange,y=-log10(padj),col=diffexpressed, size=diffexpressed)) +
  geom_point(alpha=1) +
  geom_vline(xintercept = c(-1,1), col="grey") +
  scale_x_continuous(limits = c(-2, 2), breaks = c(-1, 0, 1)) +
  geom_hline(yintercept = -log10(0.05), col="grey") +
  scale_color_manual(values = mycolour) +
  scale_size_manual(values = c("No differential expression" = 1.5,
                               "Enriched in Female" = 1.5,
                               "Enriched in Male" = 1.5)) +
  geom_text_repel(data = subset(d3_annot, !is.na(label)),
                  aes(label = label),
                  colour = "black",
                  size = 3,
                  max.overlaps = 30,
                  box.padding = 0.3,
                  segment.size = 0.2,
                  show.legend = FALSE) +
  labs(col = "Genes", size = "Genes") +
  theme(legend.position = "right",
        axis.title.x = element_text(face = "bold", size = 8, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 8, colour = "black", margin = margin(r = 10)),
        axis.text.x = element_text(face = "bold", size = 8, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 8, colour = "black"),
        legend.title = element_text(face = "bold", size = 8, colour = "black"),
        legend.text = element_text(face = "bold", size = 8, colour = "black"),
        panel.background = element_rect(fill = "grey95"),
        panel.border = element_rect(fill = NA, linewidth = 0.5, colour = "black"))

### Combine Figure 1 ####
combined_fig1 <- ((fig1a | fig1c) + plot_layout(widths = c(1, 1))) /
  (fig1b) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("A", "C", "B"))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("../../../figures/Figure1.svg",
       plot = combined_fig1,
       unit = "cm",
       width = 28,
       height = 18,
       dpi = 300)


### Supplementary Figure 1 ####
d4 <- read.xlsx("../../../SourceData.xlsx",sheet="Fig1A,S1") %>%
  pivot_longer(cols=-c(sample,sex),names_to = "module_name",values_to = "moduleTPM") %>%
  filter(module_name %in% c("IFN1_MDM","IFN2_MDM","TNF_MDM")) 

figS1 <- ggplot(data = d4, aes(x=sex, y=moduleTPM, fill=sex)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~module_name,
             scales = "free",
             ncol=3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = c("Female" = "F", "Male" = "M")) +
  stat_compare_means(
    comparisons = list(c("Female", "Male")),
    method = "wilcox",
    label = "p.signif",
    size = 2,
    bracket.size = 0.2) + 
  ylab("Module score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("../../../figures/FigureS1.svg",plot=figS1,unit="cm",width=20,height=8,dpi=300)
