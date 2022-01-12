# compare phylogenic tree and pangenome clustering dendrogram

# library ====
library(tidyverse)
library(ggtree)

# graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# data import ====
ono_tree_ori <- read.tree(file = "DATA/pangenome_analysis/775strains_tree_20200227.nwk") 
pangenome_tree_ori <- read.tree(file = "RESULTS/pangenome_analysis/sample_pangenome_clustering.nwk") 
sample_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>%  
  select(sample_ID = "Strain ID in this paper", clade = "ANI clade")

# ANI clade parameter ====
ANIbclade <- as.character(c(1:14))
ANIbcolor <- c("#dc143c", "#ff8c00", "#ff00ff", "#ffca80",	
               "#000000", "#ffb6c1", "#ffff00", "#005aff",	
               "#00bfff", "#00ffff", "#bfe4ff", "#03af7a", 
               "#00ff7f", "#a9a9a9")
clade_color <- data.frame(clade=ANIbclade, color=ANIbcolor, stringsAsFactors = F)


# change phylogenic tree tip label to dendrogram ====
otu <- ono_tree_ori$tip.label 
otu_sep <- stringr::str_split(otu, "_", simplify = T) 
ono_tree_ori$tip.label <- otu_sep[,1] 

# adjust tree size ====
ono_tree <- ono_tree_ori ; pangenome_tree <- pangenome_tree_ori
relative_value <- mean(ono_tree$edge.length)/mean(pangenome_tree$edge.length) 
pangenome_tree$edge.length <- pangenome_tree$edge.length * relative_value

# add colour info to tree order ====
onotree_label <- data.frame(label=ono_tree$tip.label, stringsAsFactors = F) %>%  
  inner_join(sample_info[,c("sample_ID", "clade")], by = c("label" = "sample_ID")) %>%  
  mutate(clade = factor(clade)) %>% 
  inner_join(clade_color, by = "clade") 

# display phylogenic tree ====
gg_onotree <- ggtree(tr = ono_tree,
                     ladderize = T,
                     layout = "rectangular") +
  layout_dendrogram() +
  scale_y_reverse() 

rgg_onotree <- ggtree::rotate(gg_onotree,1318) %>%  ggtree::rotate(1319) 

# display pangenome clustering dendrogram ====
gg_pangenome_tree <- ggtree(tr = pangenome_tree,
                            ladderize = T,
                            layout = "rectangular") +
  geom_tiplab(align = F) + 
  scale_y_reverse() 

rgg_pangenome_tree <- ggtree::rotate(gg_pangenome_tree,776) %>% 
  ggtree::rotate(913) %>% ggtree::rotate(1259) %>% ggtree::rotate(1480) %>% 
  ggtree::rotate(1481) %>% ggtree::rotate(1510) %>%  ggtree::rotate(1260) %>% 
  ggtree::rotate(1390) %>% ggtree::rotate(1402) %>% ggtree::rotate(914) %>% 
  ggtree::rotate(915) %>% ggtree::rotate(996) %>% ggtree::rotate(1135) %>% 
  ggtree::rotate(1351) %>% ggtree::rotate(1434) %>% ggtree::rotate(777) 

# compare phylogenic tree and dendrogram ====
data_ono_tree <- rgg_onotree$data 
data_pangenome_tree <- rgg_pangenome_tree$data %>%  
  mutate(x = max(x) - x + max(data_ono_tree$x) + .3) 

data_pangenome_tree_dd <- data_pangenome_tree %>% 
  select(label, pan_x = x, pan_y = y)

dd <- data_ono_tree %>% 
  mutate(x = max(x)) %>% 
  filter(isTip) %>% 
  select(label, phy_x = x, phy_y = y) %>% 
  inner_join(data_pangenome_tree_dd, by = "label") %>% 
  inner_join(onotree_label, by = "label") %>% 
  mutate(Clade = factor(clade, levels = 1:14))

# for output ====
letter_tree <- ggtree(tr = ono_tree, ladderize = T,
                      layout = "rectangular", size = .2)+
  geom_treescale(x = .4, y = -200, offset = 20, fontsize = fontsize)+ 
  layout_dendrogram()+
  scale_y_reverse()

## add dendrogram ====
tree_plot <- letter_tree %>% ggtree::rotate(1318) %>% ggtree::rotate(1319) +
  geom_tree(data = data_pangenome_tree, size = .2) + 
  geom_segment(data = dd, mapping = aes(x = phy_x, xend = pan_x, y = phy_y, yend = pan_y, color = Clade), 
               size = .2, show.legend = TRUE, key_glyph = draw_key_rect) +
  scale_color_manual(breaks = ANIbclade, values = ANIbcolor, 
                     guide = guide_legend(ncol = 2)) +  
  annotate("text", x = -.45, y = 400, label = "core gene phylogeny", 
           fontface = "bold", size = fontsize * 1.5)+ 
  annotate("text", x = .55, y = 400, label = "pangenome", fontface = "bold", size = fontsize * 1.5) 

final_plot <- tree_plot +
  theme(legend.title = element_text(size = fontsize_theme * 1.5, face = "bold"),
        legend.text = element_text(size = fontsize_theme),
        legend.key.size = unit(x = .3, units = "cm"),
        legend.position = c(0.05, 0.97),
        legend.justification = c(0, 1))

ggsave(filename = "F05_phylo_vs_dendro.pdf", path = "RESULTS/pangenome_analysis/", plot = final_plot, 
       device = "pdf", width = fig_width * 2, height = fig_height / 2, units = "mm",) 

