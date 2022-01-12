# variations of genome size and GC within cluster

# library ====
library(tidyverse)
library(scales)
library(ggrepel)
library(grid)
library(gridExtra)
library(gtable)
library(patchwork)
library(ggtext)

# graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# import data ====
sample_info <- read_tsv(file = "RESULTS/snp_distance_clustering/sample_info_SNP10cluster.tsv") %>%  
  mutate(clade = as.factor(clade), cluster = as.factor(cluster)) 

ANIbcolor <- c("#dc143c", "#ff8c00", "#ff00ff", 
               "#005aff", "#00bfff", "#00ffff", 
               "#03af7a", "#00ff7f" ) 

# get number of member of each cluster ====
cluster_sample_num <- sample_info %>% 
  select(sample_ID, cluster) %>% 
  group_by(cluster) %>% 
  summarise(sample_num = n()) %>% 
  group_by(sample_num) %>% 
  summarise(cluster_num = n()) %>% 
  ungroup() %>% 
  data.table::transpose() %>% 
  `rownames<-`(c("Cluster size", "Numbers of SNP10 clusters"))

## cluster size summary table ====
summary_table_theme <- ttheme_default(core = list(fg_params = list(fontsize = fontsize_theme), 
                                                  bg_params = list(fill = "#ffffff", col = "#000000", lwd = 1)), 
                                      rowhead = list(bg_params = list(fill = rep("#c8c8cb", nrow(cluster_sample_num)), col = "#000000", lwd = 1),
                                                     fg_params = list(fontsize = fontsize_theme * 1.5, fontface = 2))) 
summary_table <- tableGrob(d = cluster_sample_num, cols = NULL, theme = summary_table_theme) 
summary_table$widths[-1] <- rep(unit(x = 7.532, units = "mm"), 7)
table_plot <- wrap_plots(summary_table) +
  labs(tag = "A")

ggsave(filename = "SF04a_SNP_cluster_summary.pdf", plot = table_plot, device = "pdf",
       width = 85 * 1.5, height = 225 / 4, units = "mm", path = "RESULTS/snp_distance_clustering/", dpi = 1000)
# change SNP10 lettering with Illustrator

# variation in SNP10 cluster ====
## calculate var ====
cluster_summary <- sample_info %>% 
  select(sample_ID, clade, size, GC, cluster) %>% 
  group_by(clade, cluster) %>% 
  summarise(sample_num = n(), 
            size_ave = mean(size), size_diff = max(size) - min(size), 
            GC_ave = mean(GC * 100), GC_diff = max(GC * 100) - min(GC * 100)) %>%
  filter(sample_num != 1) %>% 
  ungroup() %>% 
  arrange(sample_num) %>%  mutate(sample_num = factor(sample_num)) %>% 
  arrange(clade, cluster) %>% 
  mutate(x_order = factor(c(1:n()))) 

size_diff_significant <- cluster_summary %>% 
  filter(size_diff > 150000) %>% 
  mutate(size_label = str_c(cluster, " (", comma(size_diff, accuracy = 1), " bp)", sep = ""), 
         significant = "+") 

GC_diff_significant <- cluster_summary %>% 
  filter(GC_diff > 0.2) %>% 
  mutate(GC_label = str_c(cluster, " (", comma(GC_diff), "%)", sep = ""), 
         significant = "+") 

cluster_data <- cluster_summary %>% 
  left_join(size_diff_significant[,c("cluster", "size_label")], by = "cluster") %>% 
  left_join(GC_diff_significant[,c("cluster", "GC_label")], by = "cluster") 

## plot genome size variation ====
gg_size_diff <- cluster_data %>% 
  ggplot(data = ., mapping = aes(x = x_order, y = size_diff, 
                                 label = size_label,  color = clade,
                                 shape = sample_num)) + 
  geom_point(size = .8, stroke = .2, stat = "identity", show.legend = FALSE) + 
  geom_text_repel(segment.size = .05, size = fontsize, color = "#000000", na.rm = TRUE, nudge_x = 1) + 
  scale_x_discrete(name = "", label = NULL) + 
  scale_y_continuous(name = "Genome size\nmax - minimum (bp)", labels = comma) + 
  scale_color_manual(name = "Clade", values = ANIbcolor) + 
  labs(tag = "B", shape = "Cluster\nsize") + 
  theme(panel.grid = element_line(colour = "grey92"), 
        panel.background = element_rect(fill = "#ffffff", colour = NULL, size = .25), 
        legend.key = element_blank(), 
        legend.title = element_text(size = fontsize_legend, face = "bold"), 
        legend.text = element_text(size = fontsize_legend), 
        axis.title = element_text(size = fontsize_legend), 
        axis.text = element_text(size = fontsize_legend, colour = "#000000"), 
        axis.ticks.x = element_blank()) 

## plot GC contents variation ==== 
gg_GC_diff <- cluster_data %>%  
  ggplot(data = ., mapping = aes(x = x_order, y = GC_diff, 
                                 label = GC_label, color = clade,
                                 shape = sample_num)) + 
  geom_point(size = .8, stroke = .2, stat = "identity") + 
  geom_text_repel(segment.size = .05, size = fontsize, color = "#000000", na.rm = TRUE, nudge_x = 1) + 
  scale_x_discrete(name = "Multimember SNP<sub>10</sub> clusters", label = NULL) + 
  scale_y_continuous(name = "GC content\nmax - minimum (%)") + 
  scale_color_manual(name = "Clade", values = ANIbcolor) + 
  labs(tag = "C", shape = "Cluster size") + 
  theme(panel.grid = element_line(colour = "grey92"), 
        panel.background = element_rect(fill = "#ffffff", colour = NULL, size = .25), 
        legend.key = element_blank(), 
        legend.title = element_text(size = fontsize_legend, face = "bold"), 
        legend.text = element_text(size = fontsize_legend), 
        axis.title = element_text(size = fontsize_legend), 
        axis.title.x = element_markdown(size = fontsize_legend),
        axis.text = element_text(size = fontsize_legend, colour = "#000000"), 
        axis.ticks.x = element_blank()) 

# output ====
final_plot <- gg_size_diff / gg_GC_diff + 
  guides(shape = guide_legend(override.aes = list(size = 3)),
         color = guide_legend(override.aes = list(size = 3))) +
  plot_layout(guides = "collect")

ggsave(filename = "SF04b_SNP_cluster_variation_plot.pdf", plot = final_plot, device = "pdf",
       width = 85 * 2, height = 225 / 2, units = "mm", path = "RESULTS/snp_distance_clustering/", dpi = 300)

# combine Figure A, B, C with Illustrator
