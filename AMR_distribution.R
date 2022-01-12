# script option ====
## library import ====
library(tidyverse)
library(data.table)
library(ggtree)
library(ggnewscale)
library(ggtext)

## graphic parameter ====
fontsize = 2
fontsize_AMR = 1.5
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225
a4_width = 210
a4_height = 297
mBIO_width = 174
mBIO_height = 230

## function ====
get_class_data <- function(data, filt_class){
  class_data <- data %>% 
    filter(class == filt_class) %>% 
    select(-class) %>% 
    transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
    data.frame(row.names = .$sample_ID, check.names = FALSE) %>% 
    select(-sample_ID)
}

# data import ====
## sample info ====
strain_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper",
         clade = "ANI clade",
         source = "clinical or hospital associated or nonclinical")

## AMR number ====
### import data ====
AMR <- read_tsv("DATA/final_dataset/AMR_info.tsv") %>% 
  arrange(class, `gene name`) %>% 
  mutate(`gene name` = case_when(for_MDR_dicision == "carBla" ~ str_c("*", `gene name`, sep = " "),
                                 for_MDR_dicision == "ESBL" ~ str_c("**", `gene name`, sep = " "),
                                 TRUE ~ `gene name`)) %>%
  select(-c(category, for_MDR_dicision))

### category of AMR
category_AMR <- c("beta-lactamase Class A", "beta-lactamase Class B", "beta-lactamase Class C", "beta-lactamase Class D",
                  "Aminoglycoside", "M-L-S", "Phenicol", "Tet", "S-T", "Rif", "Quinolone")

### make gene distribution data of each category ====
beta_A <- get_class_data(data = AMR, filt_class = category_AMR[1])
beta_B <- get_class_data(data = AMR, filt_class = category_AMR[2])
beta_C <- get_class_data(data = AMR, filt_class = category_AMR[3])
beta_D <- get_class_data(data = AMR, filt_class = category_AMR[4])
Agly <- get_class_data(data = AMR, filt_class = category_AMR[5])
MLS <- get_class_data(data = AMR, filt_class = category_AMR[6])
Phe <- get_class_data(data = AMR, filt_class = category_AMR[7])
Tet <- get_class_data(data = AMR, filt_class = category_AMR[8])
ST <- get_class_data(data = AMR, filt_class = category_AMR[9])
Rif <- get_class_data(data = AMR, filt_class = category_AMR[10])
Qnr <- get_class_data(data = AMR, filt_class = category_AMR[11])
Gyr <- AMR %>% 
  filter(str_starts(string = `gene name`, pattern = "gyrA")) %>% 
  select(-class) %>% 
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  data.frame(row.names = .$sample_ID, check.names = FALSE) %>% 
  select(-sample_ID) %>% 
  mutate_all(function(x){ifelse(x != "-", "+", "-")}) %>% 
  mutate_all(function(x){factor(x, levels = c("+", "-"))})
parC <- AMR %>% 
  filter(`gene name` == "parC") %>% 
  select(-class) %>% 
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  data.frame(row.names = .$sample_ID, check.names = FALSE) %>% 
  select(-sample_ID) %>% 
  mutate(parC = case_when(parC == "His75Gln" ~ "H75Q", 
                          parC == "Ser80Arg" | parC == "Ser80Ile" ~ "S80R or S80I",
                          parC == "Ala81Pro" ~ "A81P",
                          parC == "Glu84Gly" | parC == "Glu84Lys" ~ "E84G or E84K",
                          parC == "Ala108Thr" ~ "A108T",
                          TRUE ~ "-")) %>%  
  mutate(parC = factor(x = parC, levels = c("H75Q", "S80R or S80I", "A81P", "E84G or E84K", "A108T", "-")))
parE <- AMR %>% 
  filter(`gene name` == "parE") %>% 
  select(-class) %>% 
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  data.frame(row.names = .$sample_ID, check.names = FALSE) %>% 
  select(-sample_ID) %>% 
  mutate(parE = case_when(parE == "Ile444Phe" ~ "I444F",
                          parE == "Ser458Ala" | parE == "Ser458Trp" ~ "S458A or S458W",
                          TRUE ~ "-")) %>% 
  mutate(parE = factor(x = parE, levels = c("I444F", "S458A or S458W", "-")))
  
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

### make clade data for heatmap ====
clade_info <- strain_info %>% 
  select(Clade = clade) %>%  
  mutate(Clade = factor(Clade)) %>%  
  data.frame(row.names = strain_info$sample_ID) 

### add clade info heatmap ====
hm_offset_1 = .02
clade_tree <- gheatmap(p = main_tree, data = clade_info, 
                       offset = hm_offset_1, width = .02, color = NULL, 
                       colnames = T, colnames_position = "bottom", hjust = 0,
                       colnames_offset_y = 5, font.size = fontsize) + 
  scale_fill_manual(name = "Clade", values = ANIbcolor, breaks = ANIbclade, guide = guide_legend(order = 1, ncol = 2)) + 
  new_scale_fill()

## source info ====
### source colour ====
source_type <- c("clinical", "hospital associated", "others")
source_color <- c("#ff4b00", "#ff8082", "#005aff")

### make isolation source data ====
source_info <- strain_info %>% 
  select(Source = source) %>% 
  mutate(Source = factor(ifelse(Source == "n.a.", "nonclinical", Source))) %>%  
  mutate(Source = str_replace(string = Source, pattern = "nonclinical", replacement = "others")) %>% 
  data.frame(row.names = strain_info$sample_ID) 

### add isolation source heatmap ====
hm_offset_2 = hm_offset_1 + (.01 * 1.2)
source_tree <- gheatmap(p = clade_tree, data = source_info, 
                        offset = hm_offset_2, width = .02, color = NULL, 
                        colnames = T, colnames_position = "bottom", hjust = 0,
                        colnames_offset_y = 5, font.size = fontsize) + 
  scale_fill_manual(name = "Source", values = source_color, breaks = source_type, guide = guide_legend(order = 2)) +
  new_scale_fill()

## AMR distribution ====
### color ====
AMR_type <- c("+", "-")
AMR_color <- c("#005aff", "#c8c8cb")

### config ====
width_gene = .02
category_title_pos <- 0.009333835

#### beta lactamase A ====
hm_offset_3 = hm_offset_2 + (category_title_pos * 2)
beta_A_tree <- gheatmap(p = source_tree, data = beta_A, 
                        offset = hm_offset_3, width = width_gene * ncol(beta_A), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = guide_legend(order = 3)) + 
  new_scale_fill()
column_beta_A_tree <- get_heatmap_column_position(beta_A_tree, by = "bottom")
category_title_beta_A_tree <- tibble(
  label = "β-lactamase Class A",
  x = column_beta_A_tree$x[1] - category_title_pos,
  y = 0)
beta_A_tree <- beta_A_tree +
  geom_text(data = column_beta_A_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_beta_A_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### beta lactamase B ====
hm_offset_4 = hm_offset_3 + (category_title_pos * (ncol(beta_A) + 1))
beta_B_tree <- gheatmap(p = beta_A_tree, data = beta_B, 
                        offset = hm_offset_4, width = width_gene * ncol(beta_B), color = NULL, 
                        colnames = F) + 
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") + 
  new_scale_fill()
column_beta_B_tree <- get_heatmap_column_position(beta_B_tree, by = "bottom")
category_title_beta_B_tree <- tibble(
  label = "β-lactamase Class B",
  x = column_beta_B_tree$x[1] - category_title_pos,
  y = 0)
beta_B_tree <- beta_B_tree +
  geom_text(data = column_beta_B_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_beta_B_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### beta lactamase C ====
hm_offset_5 = hm_offset_4 + (category_title_pos * (ncol(beta_B) + 1))
beta_C_tree <- gheatmap(p = beta_B_tree, data = beta_C, 
                        offset = hm_offset_5, width = width_gene * ncol(beta_C), color = NULL, 
                        colnames = F) + 
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") +
  new_scale_fill()
column_beta_C_tree <- get_heatmap_column_position(beta_C_tree, by = "bottom")
category_title_beta_C_tree <- tibble(
  label = "β-lactamase Class C",
  x = column_beta_C_tree$x[1] - category_title_pos,
  y = 0)
beta_C_tree <- beta_C_tree +
  geom_text(data = column_beta_C_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_beta_C_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### beta lactamase D ====
hm_offset_6 = hm_offset_5 + (category_title_pos * (ncol(beta_C) + 1))
beta_D_tree <- gheatmap(p = beta_C_tree, data = beta_D, 
                        offset = hm_offset_6, width = width_gene * ncol(beta_D), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") +
  new_scale_fill()
column_beta_D_tree <- get_heatmap_column_position(beta_D_tree, by = "bottom")
category_title_beta_D_tree <- tibble(
  label = "β-lactamase Class D",
  x = column_beta_D_tree$x[1] - category_title_pos,
  y = 0)
beta_D_tree <- beta_D_tree +
  geom_text(data = column_beta_D_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_beta_D_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### Aminoglycoside ====
hm_offset_7 = hm_offset_6 + (category_title_pos * (ncol(beta_D) + 1))
Agly_tree <- gheatmap(p = beta_D_tree, data = Agly, 
                        offset = hm_offset_7, width = width_gene * ncol(Agly), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") + 
  new_scale_fill()
column_Agly_tree <- get_heatmap_column_position(Agly_tree, by = "bottom")
category_title_Agly_tree <- tibble(
  label = category_AMR[5],
  x = column_Agly_tree$x[1] - category_title_pos,
  y = 0)
Agly_tree <- Agly_tree +
  geom_text(data = column_Agly_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_Agly_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### M-L-S ====
hm_offset_8 = hm_offset_7 + (category_title_pos * (ncol(Agly) + 1))
MLS_tree <- gheatmap(p = Agly_tree, data = MLS, 
                        offset = hm_offset_8, width = width_gene * ncol(MLS), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") + 
  new_scale_fill()
column_MLS_tree <- get_heatmap_column_position(MLS_tree, by = "bottom")
category_title_MLS_tree <- tibble(
  label = category_AMR[6],
  x = column_MLS_tree$x[1] - category_title_pos,
  y = 0)
MLS_tree <- MLS_tree +
  geom_text(data = column_MLS_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_MLS_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### Phenicol ====
hm_offset_9 = hm_offset_8 + (category_title_pos * (ncol(MLS) + 1))
Phe_tree <- gheatmap(p = MLS_tree, data = Phe, 
                        offset = hm_offset_9, width = width_gene * ncol(Phe), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") + 
  new_scale_fill()
column_Phe_tree <- get_heatmap_column_position(Phe_tree, by = "bottom")
category_title_Phe_tree <- tibble(
  label = category_AMR[7],
  x = column_Phe_tree$x[1] - category_title_pos,
  y = 0)
Phe_tree <- Phe_tree +
  geom_text(data = column_Phe_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_Phe_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### Tet ====
hm_offset_10 = hm_offset_9 + (category_title_pos * (ncol(Phe) + 1))
Tet_tree <- gheatmap(p = Phe_tree, data = Tet, 
                        offset = hm_offset_10, width = width_gene * ncol(Tet), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") + 
  new_scale_fill()
column_Tet_tree <- get_heatmap_column_position(Tet_tree, by = "bottom")
category_title_Tet_tree <- tibble(
  label = category_AMR[8],
  x = column_Tet_tree$x[1] - category_title_pos,
  y = 0)
Tet_tree <- Tet_tree +
  geom_text(data = column_Tet_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_Tet_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### ST ====
hm_offset_11 = hm_offset_10 + (category_title_pos * (ncol(Tet) + 1))
ST_tree <- gheatmap(p = Tet_tree, data = ST, 
                        offset = hm_offset_11, width = width_gene * ncol(ST), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") + 
  new_scale_fill()
column_ST_tree <- get_heatmap_column_position(ST_tree, by = "bottom")
category_title_ST_tree <- tibble(
  label = category_AMR[9],
  x = column_ST_tree$x[1] - category_title_pos,
  y = 0)
ST_tree <- ST_tree +
  geom_text(data = column_ST_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_ST_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")

#### Rif ====
hm_offset_12 = hm_offset_11 + (category_title_pos * (ncol(ST) + 1))
Rif_tree <- gheatmap(p = ST_tree, data = Rif, 
                        offset = hm_offset_12, width = width_gene * ncol(Rif), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") + 
  new_scale_fill()
column_Rif_tree <- get_heatmap_column_position(Rif_tree, by = "bottom")
category_title_Rif_tree <- tibble(
  label = category_AMR[10],
  x = column_Rif_tree$x[1] - category_title_pos,
  y = 0)
Rif_tree <- Rif_tree +
  geom_text(data = column_Rif_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_Rif_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")


#### Quinolone ====
hm_offset_13 = hm_offset_12 + (category_title_pos * (ncol(Rif) + 1))
Qnr_tree <- gheatmap(p = Rif_tree, data = Qnr, 
                     offset = hm_offset_13, width = width_gene * ncol(Qnr), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "AMR", values = AMR_color, 
                    breaks = AMR_type, guide = "none") +
  new_scale_fill()
column_Qnr_tree <- get_heatmap_column_position(Qnr_tree, by = "bottom")
category_title_Qnr_tree <- tibble(
  label = category_AMR[11],
  x = column_Qnr_tree$x[1] - category_title_pos,
  y = 0)
Qnr_tree <- Qnr_tree +
  geom_text(data = column_Qnr_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_Qnr_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")


## point mutation distribution ====
### gyrA ====
#### color ====
gyrA_PM_type <- c("+", "-")
gyrA_PM_color <- c("#005aff", "#c8c8cb")

#### heatmap ====
hm_offset_14 = hm_offset_13 + (category_title_pos * (ncol(Qnr) + 1))
gyrA_tree <- gheatmap(p = Qnr_tree, data = Gyr, 
                      offset = hm_offset_14, width = width_gene * ncol(Gyr), color = NULL, 
                      colnames = F) +
  scale_fill_manual(name = "*gyrA*", values = gyrA_PM_color, 
                    breaks = gyrA_PM_type, guide = guide_legend(order = 4)) + 
  new_scale_fill()
column_gyrA_tree <- get_heatmap_column_position(gyrA_tree, by = "bottom") %>% 
  mutate(label_gene = str_split(string = label, pattern = "_", simplify = TRUE)[,1],
         label_position = str_split(string = label, pattern = "_", simplify = TRUE)[,2]) %>% 
  mutate(label = str_c("*", label_gene, "*_", label_position, sep = ""))
category_title_gyrA_tree <- tibble(
  label = "FQR-mutation",
  x = column_gyrA_tree$x[1] - category_title_pos,
  y = 0)
gyrA_tree <- gyrA_tree +
  geom_richtext(data = column_gyrA_tree, mapping = aes(x = x, y = y, label = label),
                hjust = 0, nudge_y = 10, label.padding = unit(c(0, 0, 0, 0), "lines"), 
                size = fontsize_AMR, fill = NA, label.colour = NA) +
  geom_text(data = category_title_gyrA_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold")
  

### parC ====
#### color ====
parC_PM_type <- c("H75Q", "S80R or S80I", "A81P", "E84G or E84K", "A108T", "-")
parC_PM_color <- c("#fff100", "#ff4b00", "#03af7a", "#f6aa00", "#990099", "#c8c8cb") 

#### heatmap ====
hm_offset_15 = hm_offset_14 + (category_title_pos * (ncol(Gyr)))
parC_tree <- gheatmap(p = gyrA_tree, data = parC, 
                      offset = hm_offset_15, width = width_gene * ncol(parC), color = NULL, 
                      colnames = F) +
  scale_fill_manual(name = "*parC*", values = parC_PM_color, 
                    breaks = parC_PM_type, guide = guide_legend(order = 5)) + 
  new_scale_fill()
column_parC_tree <- get_heatmap_column_position(parC_tree, by = "bottom") %>% 
  mutate(label = str_c("*", label, "*", sep = ""))
parC_tree <- parC_tree +
  geom_richtext(data = column_parC_tree, mapping = aes(x = x, y = y, label = label),
                hjust = 0, nudge_y = 10, label.padding = unit(c(0, 0, 0, 0), "lines"), 
                size = fontsize_AMR, fill = NA, label.colour = NA)


### parE ====
#### color ====
parE_PM_type <- c("I444F", "S458A or S458W", "-")
parE_PM_color <- c("#fff100", "#ff4b00", "#c8c8cb")

#### heatmap ====
hm_offset_16 = hm_offset_15 + (category_title_pos * (ncol(parC)))
parE_tree <- gheatmap(p = parC_tree, data = parE, 
                      offset = hm_offset_16, width = width_gene * ncol(parE), color = NULL, 
                      colnames = F) +
  scale_fill_manual(name = "*parE*", values = parE_PM_color, 
                    breaks = parE_PM_type, guide = guide_legend(order = 6)) + 
  new_scale_fill()
column_parE_tree <- get_heatmap_column_position(parE_tree, by = "bottom") %>% 
  mutate(label = str_c("*", label, "*", sep = ""))
parE_tree <- parE_tree +
  geom_richtext(data = column_parE_tree, mapping = aes(x = x, y = y, label = label),
                hjust = 0, nudge_y = 10, label.padding = unit(c(0, 0, 0, 0), "lines"), 
                size = fontsize_AMR, fill = NA, label.colour = NA)

# output ====
letter_plot <- parE_tree +
  scale_y_reverse() +
  theme(plot.margin = unit(c(0,2,0,0), units = "line"),
        legend.title = element_markdown(size = fontsize_theme * 1.2, face = "bold"),
        legend.text = element_text(size = fontsize_theme),
        legend.key.size = unit(x = .4, units = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.justification = c(0, 1),
        legend.position = c(.05, .95), 
        legend.box = "horizontal" 
  )

ggsave(filename = "RESULTS/AMR/F08_AMR_tree.pdf", plot = letter_plot,
       device = cairo_pdf, height = fig_height, width = fig_width * 2, units = "mm", dpi = 1000) 

