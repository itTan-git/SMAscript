# Distribution of plasmid replicon

# script option ====
## library import ====
library(tidyverse)
library(data.table)
library(ggtree)
library(ggnewscale)

## file path ====
file_path <- "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/plasmid_finder/"

## graphic parameter ====
fontsize = 2
fontsize_AMR = 1.5
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# data import ====
## sample info ====
strain_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper",
         clade = "ANI clade",
         source = "clinical or hospital associated or nonclinical")

## plasmid number ====
Plasmid <- read_tsv(file = str_c(file_path, "Serratia_plasmidfinder_count.tsv")) %>%
  left_join(x = select(strain_info, sample_ID), y = ., by = c("sample_ID" = "strain")) %>% 
  mutate_at(.vars = -c(1,2), .funs = function(x){case_when(is.na(x) ~ 0, TRUE ~ 1)}) %>% 
  pivot_longer(cols = -sample_ID, names_to = "rep", values_to = "pre_ab") %>% 
  group_by(rep) %>% 
  nest() %>% 
  mutate(strain_num = map_dbl(.x = data, .f = function(x){sum(x$pre_ab)})) %>% 
  filter(strain_num != 0, !str_detect(string = rep, pattern = "Gram+")) %>% 
  arrange(-strain_num) %>% 
  select(-strain_num) %>% 
  unnest(cols = data) %>% 
  pivot_wider(names_from = "rep", values_from = "pre_ab")

### output table for supplementary table ====
write_tsv(Plasmid, file = str_c(file_path, "TableS7_plasmid_replicon.tsv"), na = "")

### get strain list which have plasmid replicon ====
have_plasmid <- Plasmid %>% 
  pivot_longer(cols = -sample_ID, names_to = "rep", values_to = "pre_ab") %>% 
  group_by(sample_ID) %>% 
  summarise(rep_num = sum(pre_ab)) %>% 
  mutate(rep_num = case_when(rep_num == 0 ~ "No", TRUE ~ "Yes"))

## phylogenic tree ====
ono_tree <- read.tree(file = "DATA/pangenome_analysis/775strains_tree_20200227.nwk")

### change tree label to sample info ====
otu <- ono_tree$tip.label
otu_sep <- stringr::str_split(otu, "_", simplify = T)
ono_tree$tip.label <- otu_sep[,1]

### make tree data ====
gg_onotree <- ggtree(tr = ono_tree, ladderize = T, size = .2,
                     layout = "rectangular") +
  layout_dendrogram() +
  scale_y_reverse() 
rgg_onotree <- ggtree::rotate(gg_onotree,1318) %>%  ggtree::rotate(1319)

## make label for rep strain ====
rep_data <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = `Strain ID in this paper`, rep = `show in tip`, 
         Species, strain_name = `strain name`, TS = `Type strain`) %>% 
  mutate(rep_label = ifelse(test = !is.na(rep), 
                            yes = ifelse(test = is.na(TS), 
                                         yes = strain_name, 
                                         no = str_c(Species, "^T", sep = "")),
                            no = NA)) %>%
  inner_join(x = rgg_onotree$data, y = ., by = c("label" = "sample_ID")) %>% 
  arrange(y) %>%
  mutate(rep_y = ifelse(is.na(rep), NA, y),
         norep_y = ifelse(!is.na(rep), NA, y)) %>%
  select(label, rep_label, rep_y, norep_y)
rep_tree <- rgg_onotree %<+% rep_data

## add complete strain data ====
complete_data <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = `Strain ID in this paper`, Level) %>% 
  inner_join(x = rgg_onotree$data, y = ., by = c("label" = "sample_ID")) %>%
  mutate(comp_y = ifelse(test = Level == "complete", y, NA)) %>%
  select(label, comp_y)
complete_tree <- rep_tree %<+% complete_data 

# draw graph ====
## main_tree ====
main_tree <- complete_tree +   
  geom_tiplab(mapping = aes(y = rep_y, label = NA), 
                align = T, linetype = 5, linesize = .1, offset = -.005,
              colour = "#ff4b00", na.rm = T) +
  geom_tiplab(mapping = aes(y = rep_y, label = rep_label),
              align = T, size = fontsize_AMR, parse = TRUE,
              hjust = 0, linetype = NULL, offset = -.005,
              geom = "text", colour = "#000000", na.rm = TRUE) + 
  geom_tiplab(mapping = aes(y = norep_y, label = NA), geom = "label", align = TRUE,
              linetype = 5, linesize = .05, colour = "#000000", na.rm = TRUE) + 
  geom_tippoint(mapping = aes(x = max(rgg_onotree$data$x) + .0015, y = comp_y), 
                size = .5, shape = 23, fill = "#005aff", 
                colour = "#000000", stroke = .1, na.rm = TRUE) + 
  geom_treescale(x = .42, y = -130, offset = -20, fontsize = fontsize) +
  geom_rootedge(rootedge = .005, size = .2) +
  new_scale_fill()

## clade info ====
### color data ====
ANIbclade <- as.character(c(1:14))
ANIbcolor <- c("#dc143c", "#ff8c00", "#ff00ff", "#ffca80",	
               "#000000", "#ffb6c1", "#ffff00", "#005aff",	
               "#00bfff", "#00ffff", "#bfe4ff", "#03af7a", 
               "#00ff7f", "#a9a9a9")

### make clade information dataframe ====
clade_info <- strain_info %>% 
  select(Clade = clade) %>% 
  mutate(Clade = factor(Clade)) %>%
  data.frame(row.names = strain_info$sample_ID)

### add clade information heatmap ====
hm_offset_1 = .02
clade_tree <- gheatmap(p = main_tree, data = clade_info, 
                       offset = hm_offset_1, width = .02, color = NULL, 
                       colnames = T, colnames_position = "bottom", hjust = 0,
                       colnames_offset_y = 5, font.size = fontsize) +
  scale_fill_manual(name = "Clade", values = ANIbcolor, breaks = ANIbclade, guide = guide_legend(order = 1, ncol = 2)) +
  new_scale_fill()

## source info ====
### source colour ====
source_type <- c("clinical", "hospital associated", "nonclinical")
source_color <- c("#ff4b00", "#ff8082", "#005aff")

### make source information dataframe ====
source_info <- strain_info %>% 
  select(Source = source) %>%
  mutate(Source = factor(ifelse(Source == "n.a.", "nonclinical", Source))) %>% 
  data.frame(row.names = strain_info$sample_ID)

### add source information heatmap to tree ====
hm_offset_2 = hm_offset_1 + (0.009333835 * 1.2)
source_tree <- gheatmap(p = clade_tree, data = source_info, 
                        offset = hm_offset_2, width = .02, color = NULL, 
                        colnames = T, colnames_position = "bottom", hjust = 0,
                        colnames_offset_y = 5, font.size = fontsize) +
  scale_fill_manual(name = "Source", values = source_color, breaks = source_type, 
                    label = c("clinical", "hospital associated", "others"),
                    guide = guide_legend(order = 2)) + 
  new_scale_fill()

## have plasmid info ====
### config ====
width_rep = .02

### color ====
have_type <- c("Yes", "No")
have_color <- c("#005aff", "#c8c8cb")

### make data for heatmap ====
have_plasmid_data <- have_plasmid %>% 
  select(Plasmid = rep_num) %>% 
  data.frame(row.names = have_plasmid$sample_ID, check.names = FALSE, stringsAsFactors = FALSE)

#### no plasmid ====
hm_offset_3 = hm_offset_2 + (0.009333835 * 1.2)
have_plasmid_tree <- gheatmap(p = source_tree, data = have_plasmid_data, 
                        offset = hm_offset_3, width = width_rep * ncol(have_plasmid_data), color = NULL, 
                        colnames = T, colnames_position = "bottom", hjust = 0, 
                        colnames_offset_y = 5, font.size = fontsize) +
  scale_fill_manual(name = "Plasmid", values = have_color , 
                    breaks = have_type, guide = guide_legend(order = 3)) +
  new_scale_fill()

## plasmid rep distribution ====
### color ====
plasmid_type <- c("+", "-")
plasmid_color <- c("#005aff", "#c8c8cb")

### make data for heatmap ====
Plasmid_data <- Plasmid %>% 
  select(-sample_ID) %>% 
  mutate_all(.funs = function(x){case_when(x == 0 ~ "-", TRUE ~ "+")}) %>% 
  data.frame(row.names = Plasmid$sample_ID, check.names = FALSE, stringsAsFactors = FALSE)

#### plasmid ====
hm_offset_4 = hm_offset_3 + (0.009333835 * 1.2)
plasmid_tree <- gheatmap(p = have_plasmid_tree, data = Plasmid_data, 
                        offset = hm_offset_4, width = width_rep * 2 * ncol(Plasmid_data), color = NULL, 
                        colnames = T, colnames_position = "bottom", hjust = 0,
                        colnames_offset_y = 5, font.size = fontsize_AMR) + 
  scale_fill_manual(name = "Plasmid distribution", values = plasmid_color, 
                    breaks = plasmid_type, guide = guide_legend(order = 4)) + 
  new_scale_fill()

# output ====
letter_plot <- plasmid_tree +
  theme(plot.margin = unit(x = c(0, 2, 0, 0), units = "line"),
        legend.title = element_text(size = fontsize_theme * 1.2, face = "bold"),
        legend.text = element_text(size = fontsize_theme),
        legend.key.size = unit(x = .4, units = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.justification = c(0, 1),
        legend.position = c(.05, .95), 
        legend.box = "horizontal"
  )

ggsave(filename = str_c(file_path, "SF05_plasmid_tree.pdf") , plot = letter_plot, 
       device = "pdf", height = fig_height, width = fig_width * 2, units = "mm", dpi = 1000) 

