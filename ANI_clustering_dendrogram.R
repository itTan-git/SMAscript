# visualize dendrogram of clustering based on ANI

# library ====
library(tidyverse)
library(ggtree)
library(ggnewscale)

# graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# import data set ====
ANI_tree <- read.tree(file = "RESULTS/ANIm_clustering/ANI_clustering_tree.nwk")
node_minANI <- read_tsv(file = "RESULTS/ANIm_clustering/minimumANI_eachnode.tsv")
sample_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv")

# visualise ANI clustering tree ====
gg_ANI_tree <- ggtree(tr = ANI_tree, size = .1) +
  layout_dendrogram()
gg_ANI_tree_rotate <- gg_ANI_tree %>%  
  ggtree::rotate(776) %>% 
  ggtree::rotate(778) %>% 
  ggtree::rotate(899) %>% 
  ggtree::rotate(779) 

## make data for adding TS tip label ====
TS_info <- sample_info %>% 
  select(sample_ID = "Strain ID in this paper",
         species = Species,
         type_strain = "Type strain",
         clade = "ANI clade",
         strain_name = "strain name") %>% 
  inner_join(gg_ANI_tree_rotate$data, by = c("sample_ID" = "label")) %>%
  filter(!is.na(type_strain) | sample_ID == "DL020" | sample_ID == "DL536") %>% 
  arrange(y) %>% 
  mutate(tip_lab = ifelse(is.na(type_strain), strain_name, 
                          str_c(species, "^T", sep = ""))) %>%  
  select(label = sample_ID, ts_y = y, tip_lab)   

## make data for adding nodepoint ====
low97_node <- node_minANI %>% 
  mutate(minANI_bin = case_when(minimum_ANI < 94 ~ " - 94.0",
                                minimum_ANI > 94.0 & minimum_ANI < 94.5 ~ "94.0 - 94.5",
                                minimum_ANI > 94.5 & minimum_ANI < 95.0 ~ "94.5 - 95.0",
                                minimum_ANI > 95.0 & minimum_ANI < 95.5 ~ "95.0 - 95.5",
                                minimum_ANI > 95.5 & minimum_ANI < 96.0 ~ "95.5 - 96.0",
                                minimum_ANI > 96.0 & minimum_ANI < 96.5 ~ "96.0 - 96.5",
                                minimum_ANI > 96.5 & minimum_ANI < 97.0 ~ "96.5 - 97.0",
                                minimum_ANI > 97.0 & minimum_ANI < 97.5 ~ "97.0 - 97.5",
                                TRUE ~ "97.5 - ")) %>%  
  filter(minimum_ANI < 97.5) %>% 
  select(-minimum_ANI)

## node point parameter ====
node_fill <- c("#000000", "#83919e", "#ffffff", "#ffffff", "#83919e", "#000000")
node_colour <- c("#000000", "#83919e", "#000000", "#000000", "#83919e", "#000000")
node_shape <- rep(c(21, 25), each = 3)

## visualize tree part ====
gg_ANI_tree_ts_node <- gg_ANI_tree_rotate %<+% TS_info %<+% low97_node +
  geom_tiplab(mapping = aes(y = ts_y, label = tip_lab), parse = TRUE,
              align = TRUE, size = 1.5, angle = 270, vjust = .5,
              linetype = NULL, offset = -.25, colour = "black") + 
  geom_tippoint(mapping = aes(x = x+.25, y = ts_y), size = .5, shape = 23, 
                fill = "red", colour = "black", stroke = .1, na.rm = T) +
  geom_nodepoint(mapping = aes(subset = !is.na(minANI_bin), 
                               colour = minANI_bin, 
                               shape = minANI_bin, 
                               fill = minANI_bin),
                 stroke = .2) +
  scale_color_manual(values = node_colour,
                     guide = guide_legend(title = "minimum ANI\nwithin node", title.vjust = .5)) +
  scale_fill_manual(values = node_fill, 
                     guide = guide_legend(title = "minimum ANI\nwithin node", title.vjust = .5)) +
  scale_shape_manual(values = node_shape,
                     guide = guide_legend(title = "minimum ANI\nwithin node", title.vjust = .5, override.aes = list(size = 3))) +
  new_scale_fill()

# for adding clade info to tree ====
## make clade info data ====
clade_data <- sample_info %>% 
  select(label = `Strain ID in this paper` ,clade = `ANI clade`) %>% 
  inner_join(gg_ANI_tree_ts_node$data, by = "label") %>% 
  select(label, clade, y)

clade_number_plot <- clade_data %>% 
  group_by(clade) %>% 
  summarise(pos_x = 9, pos_y = mean(y))
clade_number_plot_keep <- clade_number_plot %>% 
  filter(clade != 8 & clade != 14)
clade_number_plot_adjust <- clade_number_plot %>% 
  filter(clade == 8 | clade == 14) %>% 
  mutate(pos_x = 11)

clade_data_hm <- clade_data %>% 
  mutate(Clade = factor(clade)) %>% 
  select(Clade) %>% 
  data.frame()
rownames(clade_data_hm) <- sample_info$`Strain ID in this paper`

## color setting ====
ANIbclade <- as.character(c(1:14))
ANIbcolor <- c("#dc143c", "#ff8c00", "#ff00ff", "#ffca80",	
               "#000000", "#ffb6c1", "#ffff00", "#005aff",	
               "#00bfff", "#00ffff", "#bfe4ff", "#03af7a", 
               "#00ff7f", "#a9a9a9") 

gg_ANI_tree_ts_node_clade <- gheatmap(p = gg_ANI_tree_ts_node, data = clade_data_hm, 
                                      offset = 4.5, width = .05, color = NULL, 
                                      colnames = T, colnames_position = "bottom",
                                      colnames_offset_y = -16, font.size = fontsize) + 
  annotate(geom = "text", x = clade_number_plot_keep$pos_x, y = clade_number_plot_keep$pos_y, 
           label = clade_number_plot_keep$clade, hjust = .5, vjust = 1, size = fontsize) + 
  annotate(geom = "text", x = clade_number_plot_adjust$pos_x, y = clade_number_plot_adjust$pos_y, 
           label = clade_number_plot_adjust$clade, hjust = .5, vjust = 1, size = fontsize) + 
  annotate(geom = "segment", x = 8, xend = clade_number_plot_adjust$pos_x, 
           y = clade_number_plot_adjust$pos_y, yend = clade_number_plot_adjust$pos_y, size = .1) + 
  scale_fill_manual(values = ANIbcolor, breaks = ANIbclade, guide = "none") +
  theme(legend.position = c(0.02, 1), 
        legend.justification = c(0,1), 
        legend.key.size = unit(x = .4, units = "cm"),
        legend.text = element_text(size = fontsize_theme),
        legend.title = element_text(size = fontsize_theme * 1.2, face = "bold"))

ggsave(filename = "F03_ANI_clustering_tree.pdf", plot = gg_ANI_tree_ts_node_clade, 
       device = "pdf", width = fig_width * 2, height = fig_height / 2.5, units = "mm", dpi = 1000,
       path = "RESULTS/ANIm_clustering/") 
