library(tidyverse)
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
  panel.spacing = unit(0.1,"cm"),
  plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")
)

# load data
d1 <- read.csv("../../../data/SourceData_BARTB_ModuleAnalysis.csv") # module analysis data
d2 <- read.csv("../../../data/SourceData_BARTB_BiomarkerAnalysis.csv") # biomarker analysis data

# subset data to keep only selected modules
dat1 <- subset(d1, module_name %in% c("STAT1_Covidsortium","TNF.Blood.RachelByngMaddick")) # Blood modules
dat2 <- subset(d1, module_name %in% c("IFN1_MDM","IFN2_MDM","IFNG_MDM","TNF_MDM")) # MDM modules

### Figure 1A-B ####
fig1a <- ggplot(data = dat1, aes(x=gender, y=moduleTPM, fill=gender)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA) +
  facet_wrap(~module_name,
             labeller = as_labeller(c(
               TNF.Blood.RachelByngMaddick="TNF_Blood",
               STAT1_Covidsortium="STAT1_Blood")),
             scales = "free",
             ncol=3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_compare_means(method = "wilcox",label = "p.signif",label.x=1.4) + 
  ylab("module score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("../../../figures/Figure1A.svg",plot=fig1a,unit="cm",width=7,height=6,dpi=300)

fig1b <- ggplot(data = d2, aes(x=gender, y=value, fill=gender)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA) +
  facet_wrap(~Module,
             scales = "free",
             ncol=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_compare_means(method = "wilcox",label = "p.signif",label.x = 1.4) + 
  ylab("signature score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("../../../figures/Figure1B.svg",plot=fig1b,unit="cm",width=14,height=6,dpi=300)


### Supplementary Figure 1 ####
figS1 <- ggplot(data = dat2, aes(x=gender, y=moduleTPM, fill=gender)) + 
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0)+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~module_name,
             scales = "free",
             ncol=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_compare_means(method = "wilcox",label = "p.signif",label.x = 1.4) + 
  ylab("module score (log2)") +
  My_Theme+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("../../../figures/FigureS1.svg",plot=figS1,unit="cm",width=14,height=6,dpi=300)
