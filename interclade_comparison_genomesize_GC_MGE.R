# Interclade difference

# library ====
library(tidyverse)
library(scales)
library(ggrepel)
library(grid)
library(gridExtra)
library(gtable)
library(patchwork)
library(ggpmisc)
library(ggbeeswarm)

# graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

## genome size and GC contents ====
## data import
sample_data <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper", picked = "Representative strain in each SNP10 cluster", 
         clade = "ANI clade", size = "Total length (bp)", GC = "GC (%)") %>%
  mutate(clade = factor(clade), GC = GC / 100)

## calculate average in each clade
sample_clade_mean <- sample_data %>% 
  filter(picked == "Rep") %>% 
  select(clade, size, GC) %>% 
  group_by(clade) %>% 
  summarise(size_mean = mean(size), GC_mean = mean(GC))

## get strain number of each clade
N_clade <- sample_data %>% 
  select(clade, picked) %>% 
  group_by(clade) %>% 
  nest() %>% 
  arrange(clade) %>% 
  mutate(SNP10 = map_dbl(data, function(x){nrow(na.omit(x))}), all = map_dbl(data, nrow)) %>% 
  select(-data) %>% 
  data.table::transpose()

make_info <- function(x){
  str_c(str_c(x[1:2], collapse = "\n"),"\n(",x[3],")", sep = "")
}
info_combine <- apply(N_clade, 2, make_info); names(info_combine) <- NULL

# ANOVA
size_anova_result <- summary(aov(size~clade, d=filter(sample_data, picked == "Rep")))
GC_anova_result <- summary(aov(GC~clade, d=filter(sample_data, picked == "Rep")))

sink("RESULTS/Gsize_GC/Gsize_GC_ANOVA_result.txt")
cat("Gsize ANOVA\n"); size_anova_result
cat("\nGC ANOVA\n"); GC_anova_result
sink()

## Tukey-Kramer
## significant threshold: P < 0.05
size_tukey <- TukeyHSD(aov(size~clade, d=filter(sample_data, picked == "Rep")))
GC_tukey <- TukeyHSD(aov(GC~clade, d=filter(sample_data, picked == "Rep")))

sink("RESULTS/Gsize_GC/Gsize_GC_Tukey-Kramer_result.txt")
cat("Gsize Tukey-Kramer\n"); size_tukey
cat("\nGC Tukey-Kramer\n"); GC_tukey
sink()

## make significant point data ====
### genome size
size_tk_result <- size_tukey$clade %>% 
  as_tibble(rownames = "clade") %>% 
  mutate(from = str_split(string = .$clade, pattern = "-", simplify = TRUE)[,1], 
         to = str_split(string = .$clade, pattern = "-", simplify = TRUE)[,2]) %>% 
  filter(`p adj` < 0.01)

size_sig_clade2 <- size_tk_result %>% 
  filter(from == 2 | to == 2) %>% 
  mutate(xpos = as.numeric(case_when(from == 2 ~ to, to == 2 ~ from)),
         ypos_max = 6500000, ypos_min = 6450000) %>% 
  select(xpos, ypos_max, ypos_min)
size_sig_clade11 <- size_tk_result %>% 
  filter(from == 11 | to == 11) %>% 
  mutate(xpos = as.numeric(case_when(from == 11 ~ to, to == 11 ~ from)),
         ypos_max = 6400000, ypos_min = 6350000) %>% 
  select(xpos, ypos_max, ypos_min)

### GC content
gc_tk_result <- GC_tukey$clade %>% 
  as_tibble(rownames = "clade") %>% 
  mutate(from = str_split(string = .$clade, pattern = "-", simplify = TRUE)[,1], 
         to = str_split(string = .$clade, pattern = "-", simplify = TRUE)[,2]) %>% 
  filter(`p adj` < 0.01) 

gc_sig_clade1 <- gc_tk_result %>% 
  filter(from == 1 | to == 1) %>% 
  mutate(xpos = as.numeric(case_when(from == 1 ~ to, to == 1 ~ from)),
         ypos_max = 0.581, ypos_min = 0.58) %>% 
  select(xpos, ypos_max, ypos_min)
gc_sig_clade2 <- gc_tk_result %>% 
  filter(from == 2 | to == 2) %>% 
  mutate(xpos = as.numeric(case_when(from == 2 ~ to, to == 2 ~ from)),
         ypos_max = 0.578, ypos_min = 0.577) %>% 
  select(xpos, ypos_max, ypos_min)
gc_sig_clade3 <- gc_tk_result %>% 
  filter(from == 3 | to == 3) %>% 
  mutate(xpos = as.numeric(case_when(from == 3 ~ to, to == 3 ~ from)),
         ypos_max = 0.575, ypos_min = 0.574) %>% 
  select(xpos, ypos_max, ypos_min)
gc_sig_clade5 <- gc_tk_result %>% 
  filter(from == 5 | to == 5) %>% 
  mutate(xpos = as.numeric(case_when(from == 5 ~ to, to == 5 ~ from)),
         ypos_max = 0.572, ypos_min = 0.571) %>% 
  select(xpos, ypos_max, ypos_min)
gc_sig_clade11 <- gc_tk_result %>% 
  filter(from == 11 | to == 11) %>% 
  mutate(xpos = as.numeric(case_when(from == 11 ~ to, to == 11 ~ from)),
         ypos_max = 0.569, ypos_min = 0.568) %>% 
  select(xpos, ypos_max, ypos_min)

## ANI clade color
ANIbclade <- as.character(c(1:14))
ANIbcolor <- c("#dc143c", "#ff8c00", "#ff00ff", "#ffca80",	
               "#000000", "#ffb6c1", "#ffff00", "#005aff",	
               "#00bfff", "#00ffff", "#bfe4ff", "#03af7a", 
               "#00ff7f", "#a9a9a9") 

## genome size boxplot ====
p_Gsize <- ggplot(data = filter(sample_data, picked == "Rep"), aes(x = clade, y = size, fill = clade)) +
  stat_boxplot(geom = "errorbar", width = 0.3) + 
  geom_boxplot(outlier.shape = 1, outlier.size = .5, outlier.stroke = .2, 
               show.legend = F, size = .25) + 
  stat_summary(data = filter(sample_data, picked == "Rep", clade != 5, clade != 14),
               fun = "mean", geom = "point", colour = "#000000", fill = "#ffffff",
               shape = 23, size = 1.25, stroke = .2, show.legend = F) + 
  scale_fill_manual(values = ANIbcolor) + 
  scale_y_continuous(name = "Genome size (bp)", limits = c(4800000, 6500000), 
                     minor_breaks = NULL, breaks = c(5000000, 5500000, 6000000, 6500000),
                     labels = c("5.0M", "5.5M", "6.0M", "6.5M"))+ 
  scale_x_discrete(name = "", labels = NULL) + 
  ### clade 2 significance
  geom_segment(data = size_sig_clade2, 
               mapping = aes(x = xpos, xend = xpos, y = ypos_max, yend = ypos_min), 
               colour = ANIbcolor[2], lineend = "square", inherit.aes = FALSE) + 
  geom_segment(data = size_sig_clade2, 
               mapping = aes(x = min(xpos), xend = max(xpos), y = ypos_max, yend = ypos_max), 
               colour = ANIbcolor[2], lineend = "square", inherit.aes = FALSE) +
  geom_point(mapping = aes(x = 2, y = 6500000), shape = 23, 
             fill = ANIbcolor[2], colour = ANIbcolor[2], 
             show.legend = FALSE) + 
  ### clade 11 significance
  geom_segment(data = size_sig_clade11, 
               mapping = aes(x = xpos, xend = xpos, y = ypos_max, yend = ypos_min), 
               colour = ANIbcolor[11], lineend = "square", inherit.aes = FALSE) + 
  geom_segment(data = size_sig_clade11, 
               mapping = aes(x = min(xpos), xend = max(xpos), y = ypos_max, yend = ypos_max), 
               colour = ANIbcolor[11], lineend = "square", inherit.aes = FALSE) + 
  geom_point(mapping = aes(x = 11, y = 6400000), shape = 23, 
             fill = ANIbcolor[11], colour = ANIbcolor[11], 
             show.legend = FALSE) + 
  labs(tag = "A") +
  theme(panel.grid = element_blank(), plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
        panel.background = element_rect(fill = "#ffffff", colour = "#000000", size = .25), 
        plot.margin = unit(c(0,0,0,0), units = "line"),
        plot.tag = element_text(vjust = 1.2),
        plot.tag.position = c(0,.9),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = fontsize_theme * 1.2, colour = "#000000"), 
        axis.text.y = element_text(size = fontsize_theme, colour = "#000000")) 

## GC content boxplot ====
p_GC <- ggplot(data = filter(sample_data, picked == "Rep"), aes(x = clade, y = GC, fill = clade)) +
  stat_boxplot(geom = "errorbar", width = 0.3) + 
  geom_boxplot(outlier.shape = 1, outlier.size = .5, outlier.stroke = .2, 
               show.legend = F, size = .25) + 
  stat_summary(data = filter(sample_data, picked == "Rep", clade != 5, clade != 14),
               fun = "mean", geom = "point", colour = "#000000", fill = "#ffffff",
               shape = 23, size = 1.25, stroke = .2, show.legend = F) + 
  scale_fill_manual(values = ANIbcolor) + 
  scale_y_continuous(name = "GC content (%)", limits = c(0.565, 0.605), 
                     minor_breaks = NULL, breaks = c(0.57, 0.58, 0.59, 0.60),
                     labels = c("57", "58", "59", "60"))+ 
  scale_x_discrete(name = "", labels = NULL)+ 
  ### clade 1 significance
  geom_segment(data = gc_sig_clade1, 
               mapping = aes(x = xpos, xend = xpos, y = ypos_max, yend = ypos_min), 
               colour = ANIbcolor[1], lineend = "square", inherit.aes = FALSE) + 
  geom_segment(data = gc_sig_clade1, 
               mapping = aes(x = 1, xend = max(xpos), y = ypos_min, yend = ypos_min), 
               colour = ANIbcolor[1], lineend = "square", inherit.aes = FALSE) +
  geom_point(mapping = aes(x = 1, y = 0.58), shape = 23, 
             fill = ANIbcolor[1], colour = ANIbcolor[1], 
             show.legend = FALSE) + 
  ### clade 2 significance
  geom_segment(data = gc_sig_clade2, 
               mapping = aes(x = xpos, xend = xpos, y = ypos_max, yend = ypos_min), 
               colour = ANIbcolor[2], lineend = "square", inherit.aes = FALSE) + 
  geom_segment(data = gc_sig_clade2, 
               mapping = aes(x = min(xpos), xend = max(xpos), y = ypos_min, yend = ypos_min), 
               colour = ANIbcolor[2], lineend = "square", inherit.aes = FALSE) +
  geom_point(mapping = aes(x = 2, y = 0.577), shape = 23, 
             fill = ANIbcolor[2], colour = ANIbcolor[2], 
             show.legend = FALSE) + 
  ### clade 3 significance
  geom_segment(data = gc_sig_clade3, 
               mapping = aes(x = xpos, xend = xpos, y = ypos_max, yend = ypos_min), 
               colour = ANIbcolor[3], lineend = "square", inherit.aes = FALSE) + 
  geom_segment(data = gc_sig_clade3, 
               mapping = aes(x = min(xpos), xend = max(xpos), y = ypos_min, yend = ypos_min), 
               colour = ANIbcolor[3], lineend = "square", inherit.aes = FALSE) +
  geom_point(mapping = aes(x = 3, y = 0.574), shape = 23, 
             fill = ANIbcolor[3], colour = ANIbcolor[3], 
             show.legend = FALSE) + 
  ### clade 5 significance
  geom_segment(data = gc_sig_clade5, 
               mapping = aes(x = xpos, xend = xpos, y = ypos_max, yend = ypos_min), 
               colour = ANIbcolor[5], lineend = "square", inherit.aes = FALSE) + 
  geom_segment(data = gc_sig_clade5, 
               mapping = aes(x = min(xpos), xend = max(xpos), y = ypos_min, yend = ypos_min), 
               colour = ANIbcolor[5], lineend = "square", inherit.aes = FALSE) +
  geom_point(mapping = aes(x = 5, y = 0.571), shape = 23, 
             fill = ANIbcolor[5], colour = ANIbcolor[5], 
             show.legend = FALSE) + 
  ### clade 11 significance
  geom_segment(data = gc_sig_clade11, 
               mapping = aes(x = xpos, xend = xpos, y = ypos_max, yend = ypos_min), 
               colour = ANIbcolor[11], lineend = "square", inherit.aes = FALSE) + 
  geom_segment(data = gc_sig_clade11, 
               mapping = aes(x = min(xpos), xend = max(xpos), y = ypos_min, yend = ypos_min), 
               colour = ANIbcolor[11], lineend = "square", inherit.aes = FALSE) +
  geom_point(mapping = aes(x = 11, y = 0.568), shape = 23, 
             fill = ANIbcolor[11], colour = ANIbcolor[11], 
             show.legend = FALSE) + 
  labs(tag = "B") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#ffffff", colour = "#000000", size = .25), 
        plot.tag.position = c(0,.91),
        plot.margin = unit(c(.2,0,0,0), units = "line"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = fontsize_theme * 1.2, vjust = -2), 
        axis.text.y = element_text(size = fontsize_theme, colour = "#000000")) 

## plasmid number beeswarm plot ====
## add plasmid info to sample_data
sample_plasmid <- read_tsv(file = "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/plasmid_finder/Serratia_plasmidfinder_plasmid_num.tsv") %>% 
  left_join(x = sample_data, y = ., by = "sample_ID") 

# beeswarm plot
gg_plasmid_num <- sample_plasmid %>% 
  filter(picked == "Rep") %>% 
  mutate(clade = factor(x = clade, levels = 1:14)) %>% 
  ggplot(mapping = aes(x = clade, y = rep_num, fill = clade)) +
  geom_quasirandom(size = 1.5, shape = 21, colour = "#000000", stroke = .2, groupOnX = TRUE) + 
  stat_summary(data = filter(sample_plasmid, picked == "Rep", clade != 5, clade != 14),
               fun = "mean", geom = "point", colour = "#000000", fill = "#ffffff",
               shape = 23, size = 1.25, stroke = .2, show.legend = F) + 
  scale_fill_manual(values = ANIbcolor) + 
  scale_x_discrete(name = "", labels = NULL) + 
  scale_y_continuous(name = "Number of plasmid replicon", breaks = seq(0,max(sample_plasmid$rep_num),1),
                     minor_breaks = NULL) + 
  guides(fill = "none") + 
  labs(tag = "C") +
  theme_bw() +
  theme(plot.tag.position = c(0,.91),
        plot.margin = unit(c(.2,0,0,0), units = "line"),
        axis.title.y = element_text(size = fontsize_theme * 1.2, vjust = -2), 
        axis.text = element_text(size = fontsize_theme, colour = "#000000"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

## integrase number beeswarm plot ====
sample_int <- read_tsv(file = "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/prophage_integrase/integrase_pangenome_num.tsv") %>% 
  left_join(x = sample_data, y = ., by = c("sample_ID" = "strain")) 

# beeswarm plot
gg_int_num <- sample_int %>% 
  filter(picked == "Rep") %>% 
  mutate(clade = factor(x = clade, levels = 1:14)) %>% 
  ggplot(mapping = aes(x = clade, y = int_num, fill = clade)) +
  geom_quasirandom(size = 1.5, shape = 21, colour = "#000000", stroke = .2, groupOnX = TRUE) + 
  stat_summary(data = filter(sample_int, picked == "Rep", clade != 5, clade != 14),
               fun = "mean", geom = "point", colour = "#000000", fill = "#ffffff",
               shape = 23, size = 1.25, stroke = .2, show.legend = F) + 
  scale_fill_manual(values = ANIbcolor) + 
  scale_x_discrete(name = "", labels = info_combine) + 
  scale_y_continuous(name = "Number of integrase", breaks = seq(0,max(sample_int$int_num),1),
                     minor_breaks = NULL) +
  guides(fill = "none") + 
  labs(tag = "D") +
  theme_bw() +
  theme(plot.tag.position = c(0,.91),
        plot.margin = unit(c(.2,0,0,0), units = "line"),
        axis.title.y = element_text(size = fontsize_theme * 1.2, vjust = -2), 
        axis.text = element_text(size = fontsize_theme, colour = "#000000"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(lineheight = 1.2))

## genome size and GC correlation ====

## make data
GC_genomesize <- sample_data %>% 
  filter(!is.na(picked)) %>% 
  group_by(clade) %>% 
  mutate(N = n()) %>% 
  filter(N > 10) %>% 
  ungroup() %>% 
  arrange(clade) %>% 
  select(ID = sample_ID, Gsize = size, GC, clade)

clade_label <- tibble(clade = unique(GC_genomesize$clade),
                      label = str_c("clade ", unique(GC_genomesize$clade)),
                      x = rep(Inf, length(unique(GC_genomesize$clade))),
                      y = rep(Inf, length(unique(GC_genomesize$clade))))

## make scatter plot
p_all <- ggplot(data = GC_genomesize, mapping = aes(x = GC, y = Gsize)) + 
  geom_point(mapping = aes(fill = clade), shape = 21, colour = "#000000", size = .8, stroke = .2) + 
  stat_smooth(mapping = aes(x = GC, y = Gsize), 
              formula = y ~ x, method = "lm", size = .5, 
              se = FALSE, colour = "#000000", show.legend = FALSE) + 
  stat_poly_eq(formula = y ~ x, parse = TRUE, size = fontsize, 
               label.x = .05, label.y = .1) + 
  geom_text(data = clade_label, mapping = aes(x = x, y = y, label = label), 
            size = fontsize, hjust = 1.1, vjust = 1.5, fontface = "bold") +
  facet_wrap(~ clade, ncol = 1) +
  scale_fill_manual(name = "", values = ANIbcolor[unique(GC_genomesize$clade)]) + 
  scale_colour_manual(name = "", values = ANIbcolor[unique(GC_genomesize$clade)]) + 
  scale_x_continuous(name = "GC content (%)", limits = c(0.580, 0.605), 
                     minor_breaks = NULL, breaks = c(0.58, 0.59, 0.60),
                     labels = c("58", "59", "60")) + 
  scale_y_continuous(name = "Genome size (bp)", limits = c(4800000, 7000000), 
                     minor_breaks = NULL, breaks = c(5000000, 5500000, 6000000, 6500000),
                     labels = c("5.0M", "5.5M", "6.0M", "6.5M")) + 
  guides(fill = "none", colour = "none") +
  labs(tag = "E") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#ffffff", colour = "#000000", size = .25),
        plot.margin = unit(c(0, 0, 0, .5), units = "line"),
        legend.key = element_rect(fill = "#ffffff"),
        legend.text = element_text(size = fontsize_legend),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.tag = element_text(vjust = -4),
        axis.title.x = element_text(size = fontsize_theme * 1.2),
        axis.title.y = element_text(size = fontsize_theme * 1.2),
        axis.text.y = element_text(size = fontsize_theme, colour = "#000000"),
        axis.text.x = element_text(size = fontsize_theme, colour = "#000000")) 


## plasmid number and genome size correlation ====
plasmid_genomesize <- sample_plasmid %>% 
  filter(!is.na(picked)) %>% 
  group_by(clade) %>% 
  mutate(N = n()) %>% 
  filter(N > 10) %>% 
  ungroup() %>% 
  arrange(clade) %>% 
  select(clade, size, rep_num)

p_plasmid_genomesize <- ggplot(data = plasmid_genomesize, mapping = aes(x = rep_num, y = size)) + 
  geom_point(mapping = aes(fill = clade), shape = 21, colour = "#000000", size = .8, stroke = .2) + 
  stat_smooth(mapping = aes(x = rep_num, y = size), 
              formula = y ~ x, method = "lm", size = .5, 
              se = FALSE, colour = "#000000", show.legend = FALSE) + 
  stat_poly_eq(formula = y ~ x, parse = TRUE, size = fontsize, 
               label.x = .05, label.y = .9) + 
  geom_text(data = clade_label, mapping = aes(x = x, y = y, label = label), 
            size = fontsize, hjust = 1.1, vjust = 1.5, fontface = "bold") +
  facet_wrap(~ clade, ncol = 1) +
  scale_fill_manual(name = "", values = ANIbcolor[unique(plasmid_genomesize$clade)]) + 
  scale_colour_manual(name = "", values = ANIbcolor[unique(plasmid_genomesize$clade)]) + 
  scale_x_continuous(name = "Number of plasmid replicon") +
  scale_y_continuous(name = "Genome size (bp)", limits = c(4800000, 7000000), 
                     minor_breaks = NULL, breaks = c(5000000, 5500000, 6000000, 6500000),
                     labels = c("5.0M", "5.5M", "6.0M", "6.5M")) + 
  guides(fill = "none", colour = "none") +
  labs(tag = "") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#ffffff", colour = "#000000", size = .25),
        plot.margin = unit(c(0, 0, 0, 0), units = "line"),
        legend.key = element_rect(fill = "#ffffff"),
        legend.text = element_text(size = fontsize_legend),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.tag = element_blank(),
        axis.title.x = element_text(size = fontsize_theme * 1.2),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = fontsize_theme, colour = "#000000")) 

## integrase number and genome size correlation ====
integrase_genomesize <- sample_int %>% 
  filter(!is.na(picked)) %>% 
  group_by(clade) %>% 
  mutate(N = n()) %>% 
  filter(N > 10) %>% 
  ungroup() %>% 
  arrange(clade) %>% 
  select(clade, size, int_num)

p_integrase_genomesize <- ggplot(data = integrase_genomesize, mapping = aes(x = int_num, y = size)) + 
  geom_point(mapping = aes(fill = clade), shape = 21, colour = "#000000", size = .8, stroke = .2) + 
  stat_smooth(mapping = aes(x = int_num, y = size), 
              formula = y ~ x, method = "lm", size = .5, 
              se = FALSE, colour = "#000000", show.legend = FALSE) + 
  stat_poly_eq(formula = y ~ x, parse = TRUE, size = fontsize, 
               label.x = .05, label.y = .9) + 
  geom_text(data = clade_label, mapping = aes(x = x, y = y, label = label), 
            size = fontsize, hjust = 1.1, vjust = 1.5, fontface = "bold") +
  facet_wrap(~ clade, ncol = 1) +
  scale_fill_manual(name = "", values = ANIbcolor[unique(integrase_genomesize$clade)]) + 
  scale_colour_manual(name = "", values = ANIbcolor[unique(integrase_genomesize$clade)]) + 
  scale_x_continuous(name = "Number of integrase") +
  scale_y_continuous(name = "Genome size (bp)", limits = c(4800000, 7000000), 
                     minor_breaks = NULL, breaks = c(5000000, 5500000, 6000000, 6500000),
                     labels = c("5.0M", "5.5M", "6.0M", "6.5M")) + 
  guides(fill = "none", colour = "none") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#ffffff", colour = "#000000", size = .25),
        plot.margin = unit(c(0, 0, 0, 0), units = "line"),
        legend.key = element_rect(fill = "#ffffff"),
        legend.text = element_text(size = fontsize_legend),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.tag = element_blank(),
        axis.title.x = element_text(size = fontsize_theme * 1.2),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = fontsize_theme, colour = "#000000")) 

# output plot ====
letter_plot_abcd <- p_Gsize / p_GC / gg_plasmid_num / gg_int_num
letter_plot_e <- p_all | p_plasmid_genomesize | p_integrase_genomesize
letter_plot <- letter_plot_abcd | letter_plot_e *
  theme(axis.title.x = element_text(vjust = 5))
ggsave(filename = "RESULTS/Gsize_GC/F06a_final_output.pdf", plot = letter_plot, 
       device = "pdf", width = fig_height, height = fig_width * 2, units = "mm")

# add clade, SNP10, (all) in left of x axis of plot D with Illustrator
