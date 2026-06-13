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
d <- read.csv("../../../data/SourceData_BAL_ModuleAnalysis.csv") %>%
  filter(module %in% c("IFN1_MDM","IFN2_MDM","IFNG_MDM","TNF_MDM")) %>% # MDM modules
  mutate(cluster = recode(celltypes,
                           "0_AM_APOC2" = "0",
                           "1_AM_IFI27" = "1",
                           "2_Neutrophil_IFIT2_IL1RN" = "2",
                           "3_Neutrophil_CXCR4-high" = "3",
                           "4_CD4 T cell_Effector" = "4",
                           "5_CD8A T cell" = "5",
                           "6_CD4 T cell_Naive_Stem" = "6",
                           "7_gd T cell_NK_ILC" = "7",
                           "8_Treg_FOXP3" = "8",
                           "9_T cell_proliferating" = "9",
                           "10_AM" = "10",
                           "11_AM" = "11",
                           "12_Mo-derived AM_APOE" = "12",
                           "13_AM_FN1_CCL18" = "13",
                           "14_AM_Mito-high" = "14",
                           "15_Mo-derived AM_FPR3" = "15",
                           "16_AM_IFN-responsive" = "16",
                           "17_AM_stressed" = "17",
                           "18_Monocyte" = "18",
                           "19_cDC1_Monocyte" = "19",
                           "20_Migratory DC" = "20",
                           "21_AM_proliferating" = "21",
                           "22_B cell" = "22",
                           "23_Mast cell" = "23",
                           "24_pDC_Plasma cell" = "24",
                           "25_AM_RSRP1_NRAMP1" = "25",
                           "26_Mo-derived AM_TPP1" = "26"),
         celltype_broad = ifelse(cluster %in% c("0","1","10","11","12","13","14","15","16","17","21","25","26"),"AM",
                                 ifelse(cluster %in% c("2","3"),"Neut",
                                        ifelse(cluster %in% c("4","5","6","7","8","9"),"T",
                                               ifelse(cluster == "22","B",
                                                      ifelse(cluster %in% c("18","19","20","24"),"MoDC",
                                                             ifelse(cluster == "23","Mast",NA)))))))
d$celltype_broad <- factor(d$celltype_broad,levels = c("AM","MoDC","Neut","T","B","Mast"))

# stats (multiple testing correction)
stats <- d %>%
  group_by(module, cluster,celltype_broad) %>%
  wilcox_test(average ~ sex) %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup() %>%
  mutate(y.position = max(d$average, na.rm = TRUE) + 0.1,
         p.adj.label = ifelse(p.adj > 0.05, "ns", p.adj))

# Figure 4b
fig4b <- ggplot(data = d, aes(x=sex, y=average, fill=sex)) + 
  geom_jitter(size = 0.6, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape=NA) +
  facet_grid(module~celltype_broad + cluster) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  stat_pvalue_manual(stats, label = "p.adj.label", tip.length = 0.01) +
  ylab("module Z-score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
fig4b
ggsave("../../../figures/Figure4B.svg",plot=fig4b,units="cm",width=25,height=12,dpi=300)
