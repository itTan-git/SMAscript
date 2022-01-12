# Visualyzation of Serratia spp. phylogenic tree
# Display ANI matrix
# Display 16S percent identity matrix

# library ====
library(tidyverse)
library(ggtree)
library(ggtext)

# graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# data import ====
S_spp_tree <- read.tree(file = "DATA/Serratia_spp/Serratia_spp_20210127_tree_reroot.nwk") 

## phylogenic tree ====
### change tip label ====
S_spp_tree$tip.label <- c("*Y. enterocolitica*<sup>T</sup>", "*S. microhaemolytica*<sup>T</sup>", 
                          "*S. fonticola*<sup>T</sup>", "*S. oryzae*<sup>T</sup>",
                          "*S. rubidaea*<sup>T</sup>", "*S. odorifera*<sup>T</sup>",
                          "*S. plymuthica*<sup>T</sup>", "*S. liquefaciens*<sup>T</sup>",
                          "*S. grimesii*<sup>T</sup>", "*S. proteamaculans*<sup>T</sup>",
                          "*S. quinivorans*<sup>T</sup>", "*S. marcescens*<sup>T</sup> (Sma<sup>T</sup>)",
                          "*S. marcescens* subsp. *sakuensis*<sup>T</sup> (Sma_sak<sup>T</sup>)", 
                          "*S. nematodiphila*<sup>T</sup> (Sne<sup>T</sup>)", 
                          "*S. ureilytica*<sup>T</sup> (Sur<sup>T</sup>)",  "Sma Db11", "Sma SM39",
                          "*S. surfactantfaciens*<sup>T</sup> (Ssu<sup>T</sup>)", "*S. ficaria*<sup>T</sup>",
                          "*S. symbiotica*<sup>T</sup>")

### draw tree ==== 
gg_S_spp_tree <- ggtree(tr = S_spp_tree, ladderize = T, size = .2,
                        layout = "rectangular")
data_tree <- gg_S_spp_tree$data %>% 
  filter(isTip)
gg_S_spp_tree <- gg_S_spp_tree +
  geom_richtext(data = data_tree, mapping = aes(x = x, y = y, label = label),
                fill = NA, label.colour = NA, size = fontsize, hjust = 0) +
  geom_treescale(x = 0, y = 4.5, offset = 0.5, fontsize = fontsize) +
  xlim(0, 1.1) +
  theme(plot.margin = unit(c(0,0,0,0), units = "line")) 

ggsave(filename = "F01a_serratia_spp_phylo.pdf", plot = gg_S_spp_tree, device = "pdf", 
       width = fig_width, height = fig_height / 3, units = "mm", path = "RESULTS/Serratia_spp/", dpi = 1000)

## ANI matrix ====
### get tree order ====
tree_order <- c("Sma_sak<sup>T</sup>", 
                "Sma<sup>T</sup>",
                "Sne<sup>T</sup>",
                "Sma Db11", 
                "Sur<sup>T</sup>",
                "Sma SM39",
                "Ssu<sup>T</sup>",
                "*S. symbiotica*<sup>T</sup>",
                "*S. ficaria*<sup>T</sup>")

### ANI matrix manipulation ====
ANI_score <- c("96-100%", "95-96%", "94-95%", "< 94%")
S_spp_ani <- read_tsv(file = "DATA/Serratia_spp/Serratia_spp_SMcomplex_ANIb.tsv") %>% 
  mutate(tree_order = tree_order) %>% 
  select(tree_order, starts_with(match = "Serratia")) %>% 
  `colnames<-`(c("name", tree_order)) %>% 
  pivot_longer(cols = -name, names_to = "to", values_to = "ANI") %>% 
  mutate(name = factor(name, levels = tree_order),
         to = factor(to, levels = tree_order),
         ANI = as.numeric(sprintf("%.2f", ANI*100))) %>% 
  mutate(ANI_category = case_when(96 <= ANI ~ ANI_score[1],
                                  95 <= ANI & ANI < 96 ~ ANI_score[2],
                                  94 <= ANI & ANI < 95 ~ ANI_score[3],
                                  TRUE ~ ANI_score[4]))

ANI_color <- c("#ff9366", "#68cfaf", "#669cff", "#ffffff")

### draw ANI matrix ====
g_S_spp_ani <- ggplot(data = S_spp_ani, mapping = aes(x = name, y = to, fill = ANI_category)) +
  geom_tile(colour = "#000000") +
  geom_text(mapping = aes(label = ifelse(ANI == 100, ANI, sprintf(fmt = "%.2f", ANI))), size = fontsize) +
  scale_x_discrete(name = "", position = "top") +
  scale_y_discrete(name = "", limits = rev) +
  scale_fill_manual(name = "", values = ANI_color, breaks = ANI_score, guide = guide_legend(nrow = 1)) +
  theme_minimal() +
  theme(axis.text.x.top = element_markdown(hjust = 0, vjust = .5, angle = 90, size = fontsize_theme, colour = "#000000"),
        axis.text.y = element_markdown(size = fontsize_theme, hjust = 1, vjust = .5, colour = "#000000"), 
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0, 1.5, 0), units = "line"),
        legend.title = element_blank(),
        legend.position = c(.4, -.05),
        legend.key.size = unit(x = .25, units = "cm"),
        legend.text = element_text(size = fontsize_legend))

ggsave(filename = "F01b_serratia_spp_ANI.pdf", plot = g_S_spp_ani, device = "pdf", 
       width = fig_width, height = fig_height / 3, units = "mm", path = "RESULTS/Serratia_spp/", dpi = 1000)

# combine F01a and F01b, and add tag with Illustrator

## 16S sequence blastn identity matrix ====
### data manipulation ====
S_spp_16s <- read_tsv(file = "DATA/Serratia_spp/16SrRNA/16SrRNA_PIM_matrix.tsv") %>% 
  head(9) %>% 
  mutate(tree_order = tree_order) %>% 
  select(tree_order, starts_with(match = "Serratia")) %>% 
  select(1:10) %>% 
  `colnames<-`(c("name", tree_order)) %>% 
  pivot_longer(cols = -name, names_to = "to", values_to = "ANI") %>% 
  mutate(name = factor(name, levels = tree_order),
         to = factor(to, levels = tree_order),
         ANI = as.numeric(sprintf("%.2f", ANI)))

### draw matrix ====
g_S_spp_16s <- ggplot(data = S_spp_16s, mapping = aes(x = name, y = to, fill = ANI)) +
  geom_tile(colour = "#000000") +
  geom_text(mapping = aes(label = ifelse(ANI == 100, ANI, sprintf(fmt = "%.2f", ANI))), size = fontsize) +
  scale_x_discrete(name = "", position = "top") +
  scale_y_discrete(name = "", limits = rev) +
  scale_fill_gradient2(name = "", 
                       low = "#005aff", mid = "#ffffff", high = "#ff4b00", 
                       midpoint = mean(c(max(S_spp_16s$ANI), min(S_spp_16s$ANI))),
                       guide = guide_colorbar(frame.colour = "#000000", ticks.colour = "#000000")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_markdown(hjust = 0, angle = 90, size = fontsize_theme, colour = "#000000"),
        axis.text.y = element_markdown(size = fontsize_theme, colour = "#000000"),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = fontsize * 4))

ggsave(filename = "SF03_Serratia_spp_16S_pim.pdf", plot = g_S_spp_16s, device = "pdf",
       width = fig_width * 1.8, height = fig_height / 2.3, units = "mm", path = "RESULTS/Serratia_spp/", dpi = 1000)
