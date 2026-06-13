library(tidyverse)
library(rstatix)
library(ggpubr)

#My_Theme
t = 8 #size of text
m = 4 #size of margin around text
tc = "black" #colour of text
My_Theme = theme(
  axis.title.x = element_text(size = t, face = "bold", margin = margin(t = m)),
  axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = 0, hjust = 0.5),
  axis.title.y.left = element_text(size = t, face = "bold", margin = margin(r = m)),
  axis.title.y.right = element_text(size = t, face = "bold", margin = margin(l = m)),
  axis.text.y = element_text(size = t, face = "bold", colour = tc),
  legend.title = element_text(size=t, face = "bold", colour = tc),
  legend.text = element_text(size=t, face = "bold", colour = tc),
  plot.title = element_text(size=t, face = "bold", colour = tc),
  strip.text = element_text(size=t, face = "bold", colour = tc),
  strip.background = element_rect(fill = "gray90", colour = "black", linewidth = 0.5),
  panel.border = element_rect(fill = NA, linewidth = 0.5, colour = tc),
  panel.background = element_rect(fill = "gray97"),
  legend.position = "right", legend.justification = "top",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0,0,0,0),
  legend.box.spacing = unit(c(0,0,0,0),"cm"),
  panel.spacing = unit(0.05,"cm"),
  plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm")
)

# load data
d <- read.csv("../../../data/SourceData_10X_ModuleAnalysis.csv") %>%
  filter(!CellType %in% c("undefined","erythrocytes")) %>%
  mutate(Gender = recode(Gender,
                         "Female" = "F",
                         "Male" = "M"),
         module = recode(module,
                         IFNA2_TST.URA = "IFNA2_TST",
                         IFNG_TST.URA = "IFNG_TST",
                         TNF_TST.URA = "TNF_TST",
                         IFNG_KC_v2 = "IFNG_KC",
                         TNF_KC_v2 = "TNF_KC"),
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
d$CellType <- factor(d$CellType, levels = c("CD4","CD8","CD8_atyp.","gd_V1","gd_V2","NKT","NK","myeloid","neutro","KC_b","KC_sb","MC"))

# subset data to keep only selected modules
dat1 <- subset(d, module %in% c("IFNA2_TST","IFNG_TST","TNF_TST")) # TST modules
dat2 <- subset(d, module %in% c("IFNA_KC","IFNG_KC","TNF_KC")) # KC modules
dat3 <- subset(d, module %in% c("IFN1_MDM","IFN2_MDM","IFNG_MDM","TNF_MDM")) # MDM modules

# stats (multiple testing correction)
stats1 <- dat1 %>%
  group_by(module, CellType) %>%
  wilcox_test(average ~ Gender) %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup() %>%
  mutate(y.position = max(dat1$average, na.rm = TRUE) + 0.5,
         p.adj.label = ifelse(p.adj > 0.05, "ns", p.adj))
stats2 <- dat2 %>%
  group_by(module, CellType) %>%
  wilcox_test(average ~ Gender) %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup() %>%
  mutate(y.position = max(dat2$average, na.rm = TRUE) + 0.5,
         p.adj.label = ifelse(p.adj > 0.05, "ns", p.adj))
stats3 <- dat3 %>%
  group_by(module, CellType) %>%
  wilcox_test(average ~ Gender) %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup() %>%
  mutate(y.position = max(dat3$average, na.rm = TRUE) + 0.5,
         p.adj.label = ifelse(p.adj > 0.05, "ns", p.adj))

### Figure 3A ####
fig3a <- ggplot(data = dat1, aes(x=Gender, y=average, fill=Gender)) + 
  geom_jitter(size = 0.6, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_grid(module~CellType) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_pvalue_manual(stats1, label = "p.adj.label", tip.length = 0.01) +
  ylab("module Z-score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
fig3a
ggsave("../../../figures/Figure3A.svg",plot=fig3a,units="cm",width=20,height=10,dpi=300)

### Supplementary Figure 3 ####
figS3a <- ggplot(data = dat2, aes(x=Gender, y=average, fill=Gender)) + 
  geom_jitter(size = 0.6, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA) +
  facet_grid(module~CellType) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_pvalue_manual(stats2, label = "p.adj.label", tip.length = 0.01) +
  ylab("module Z-score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
figS3a
ggsave("../../../figures/FigureS3A.svg",plot=figS3a,units="cm",width=20,height=10,dpi=300)

figS3b <- ggplot(data = dat3, aes(x=Gender, y=average, fill=Gender)) + 
  geom_jitter(size = 0.6, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape=NA) +
  facet_grid(module~CellType) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_pvalue_manual(stats3, label = "p.adj.label", tip.length = 0.01) +
  ylab("module Z-score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
figS3b
ggsave("../../../figures/FigureS3B.svg",plot=figS3b,units="cm",width=20,height=12,dpi=300)
