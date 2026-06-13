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
  panel.spacing = unit(0.05,"cm"),
  plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm")
)

# load data
dat <- read.csv("../../../data/SourceData_BAL_DA.csv")

# calculate abundance of each cell type as percentage of all BAL cells
d1 <- dat %>%
  group_by(sample,sex) %>%
  mutate(n_total = n()) %>%
  group_by(sample,sex,celltypes) %>%
  summarise(n=n(),
         pct=n/n_total*100) %>%
  ungroup() %>%
  unique()

# select only B cells
d2 <- d1 %>%
  filter(celltypes == "22_B cell")

# copy stats from edgeR differential abundance analysis
stats <- tibble::tribble(
  ~group1, ~group2, ~p.adj,
  "F","M",0.015
)

# plot figure
fig4A <- ggplot(d2, aes(x=sex,y=pct,colour=sex))+
  geom_jitter(width = 0.2, alpha = 1, size = 2, show.legend = F) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black") +
  labs(title = "22_B cell",
       y = "% of all BAL cells") +
  stat_pvalue_manual(stats, y.position = 2, label = "p.adj") +
  My_Theme
fig4A
ggsave("../../../figures/Figure4A.svg",fig4A,width = 6,height = 6,units = "cm",dpi = 300)
