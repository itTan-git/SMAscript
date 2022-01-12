# script option ====
## library import ====
library(tidyverse)
library(data.table)
library(ggtree)
library(ggnewscale)
library(ggtext)

## graphic parameter ====
fontsize = 1.5
fontsize_AMR = 1.2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
a4_width = 210
a4_height = 297

## function ====
get_category_data <- function(data, filt_category){
  category_data <- data %>% 
    filter(category == filt_category) %>% 
    select(-c(category, cluster_type)) %>% 
    transpose(keep.names = "sample_ID", make.names = "Gene") %>% 
    data.frame(row.names = .$sample_ID, check.names = FALSE) %>% 
    select(-sample_ID)
}

## file path ====
file_path <- "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/virulence_gene/"

# data import ====
## sample info ====
strain_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper",
         clade = "ANI clade",
         source = "clinical or hospital associated or nonclinical")

# virulence gene data import ====
vir_original <- read_tsv(file = str_c(file_path, "vir_gene_list.tsv")) 

## virulence gene data for analysis ====
vir <- vir_original %>% 
### convert presense absence info by SM39/Db11 cluster
  mutate(cluster_type = case_when(!is.na(SM39) & !is.na(Db11) ~ "SM39 and Db11", 
                                  !is.na(SM39) & is.na(Db11) ~ "SM39", 
                                  TRUE ~ "Db11")) %>% 
  select(cluster_type, Gene, category, starts_with(match = c("DL", "SE"))) %>% 
  pivot_longer(cols = -c(1:3), names_to = "strain", values_to = "pre_ab") %>% 
  mutate(pre_ab = case_when(pre_ab == "+" ~ cluster_type, TRUE ~ pre_ab)) %>% 
  pivot_wider(names_from = "strain", values_from = "pre_ab")

category_list <- unique(vir$category)

### each category data ====
sig_pep_p <- get_category_data(data = vir, filt_category = category_list[1])
sig_pep_n <- get_category_data(data = vir, filt_category = category_list[2])
TypeIV <- get_category_data(data = vir, filt_category = category_list[3])
flagellar <- get_category_data(data = vir, filt_category = category_list[4])
chemotaxis <- get_category_data(data = vir, filt_category = category_list[5])
methyl_chemotaxis <- get_category_data(data = vir, filt_category = category_list[6])
Fimbriae <- get_category_data(data = vir, filt_category = category_list[7])
TypeIV_pili <- get_category_data(data = vir, filt_category = category_list[8])
TypeII <- get_category_data(data = vir, filt_category = category_list[9])
hemolysin <- get_category_data(data = vir, filt_category = category_list[10])
ent <- get_category_data(data = vir, filt_category = category_list[11])
siderophore <- get_category_data(data = vir, filt_category = category_list[12])
fhu <- get_category_data(data = vir, filt_category = category_list[13])
fep <- get_category_data(data = vir, filt_category = category_list[14])
fec <- get_category_data(data = vir, filt_category = category_list[15])
has <- get_category_data(data = vir, filt_category = category_list[16])
hem <- get_category_data(data = vir, filt_category = category_list[17])
hms <- get_category_data(data = vir, filt_category = category_list[18])
TonB <- get_category_data(data = vir, filt_category = category_list[19])
serrawettin <- get_category_data(data = vir, filt_category = category_list[20])
vir_mutants <- get_category_data(data = vir, filt_category = category_list[21])

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
              align = T, size = fontsize, parse = TRUE,
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

### clade data info ====
clade_info <- strain_info %>% 
  select(Clade = clade) %>% 
  mutate(Clade = factor(Clade)) %>%
  data.frame(row.names = strain_info$sample_ID)

### add clade info to tree ====
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

### make isolation source info ====
source_info <- strain_info %>% 
  select(Source = source) %>% 
  mutate(Source = factor(ifelse(Source == "n.a.", "nonclinical", Source))) %>%  
  mutate(Source = str_replace(string = Source, pattern = "nonclinical", replacement = "others")) %>% 
  data.frame(row.names = strain_info$sample_ID) 

### add isolation source info to tree ====
hm_offset_2 = hm_offset_1 + (.01 * 1.2)
source_tree <- gheatmap(p = clade_tree, data = source_info, 
                        offset = hm_offset_2, width = .02, color = NULL, 
                        colnames = T, colnames_position = "bottom", hjust = 0,
                        colnames_offset_y = 5, font.size = fontsize) + 
  scale_fill_manual(name = "Source", values = source_color, breaks = source_type, guide = guide_legend(order = 2)) + 
  new_scale_fill()

## virulence gene distribution ====
### color ====
vir_type <- c("SM39", "SM39 and Db11", "Db11", "-")
vir_color <- c(ANIbcolor[1],"#005aff", ANIbcolor[9], "#c8c8cb")

### config ====
width_gene = .02

### category title position ===
category_title_pos <- 0.009333835

### locus ID of SM39 and Db11 ====
vir_ID <- vir_original %>% 
  select(Gene, SM39, Db11) %>% 
  mutate_at(.vars = -1, function(x){case_when(is.na(x) ~ "-", TRUE ~ x)})

#### Predictable signal peptide (+) ====
hm_offset_3 = hm_offset_2 + (category_title_pos * 2)
sig_pep_p_tree <- gheatmap(p = source_tree, data = sig_pep_p, 
                           offset = hm_offset_3, width = width_gene * ncol(sig_pep_p), color = NULL, 
                           colnames = F) +
  scale_fill_manual(name = "Potentially virulence-related genes/operons<br>in SM39 and Db11 (Iguchi, *et al.*, 2014)", 
                    values = vir_color, 
                    labels = c("SM39-specific genes/operons", "genes/operons shared by SM39 and Db11", 
                               "Db11-specific genes/operons", "absent"),
                    breaks = vir_type, guide = guide_legend(order = 3)) + 
  new_scale_fill()
column_sig_pep_p_tree <- get_heatmap_column_position(sig_pep_p_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_sig_pep_p_tree <- tibble(
  label = category_list[1],
  x = column_sig_pep_p_tree$x[1] - category_title_pos,
  y = 0)
sig_pep_p_tree <- sig_pep_p_tree +
  geom_text(data = column_sig_pep_p_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_sig_pep_p_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### Predictable signal peptide (-) ====
hm_offset_4 = hm_offset_3 + (category_title_pos * (ncol(sig_pep_p) + 1))
sig_pep_n_tree <- gheatmap(p = sig_pep_p_tree, data = sig_pep_n, 
                           offset = hm_offset_4, width = width_gene * ncol(sig_pep_n), color = NULL, 
                           colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_sig_pep_n_tree <- get_heatmap_column_position(sig_pep_n_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_sig_pep_n_tree <- tibble(
  label = category_list[2],
  x = column_sig_pep_n_tree$x[1] - category_title_pos,
  y = 0)
sig_pep_n_tree <- sig_pep_n_tree +
  geom_text(data = column_sig_pep_n_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_sig_pep_n_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### Type VI secretion system-associated ====
hm_offset_5 = hm_offset_4 + (category_title_pos * (ncol(sig_pep_n) + 1))
TypeIV_tree <- gheatmap(p = sig_pep_n_tree, data = TypeIV, 
                        offset = hm_offset_5, width = width_gene * ncol(TypeIV), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_TypeIV_tree <- get_heatmap_column_position(TypeIV_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_TypeIV_tree <- tibble(
  label = "T6SS-related",
  x = column_TypeIV_tree$x[1] - category_title_pos,
  y = 0)
TypeIV_tree <- TypeIV_tree +
  geom_text(data = column_TypeIV_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_TypeIV_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### flagellar biosynthesis ====
hm_offset_6 = hm_offset_5 + (category_title_pos * (ncol(TypeIV) + 1))
flagellar_tree <- gheatmap(p = TypeIV_tree, data = flagellar, 
                           offset = hm_offset_6, width = width_gene * ncol(flagellar), color = NULL, 
                           colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_flagellar_tree <- get_heatmap_column_position(flagellar_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_flagellar_tree <- tibble(
  label = "flagellar biosynthesis-related",
  x = column_flagellar_tree$x[1] - category_title_pos,
  y = 0)
flagellar_tree <- flagellar_tree +
  geom_text(data = column_flagellar_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_flagellar_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### chemotaxis-related ====
hm_offset_7 = hm_offset_6 + (category_title_pos * (ncol(flagellar) + 1))
chemotaxis_tree <- gheatmap(p = flagellar_tree, data = chemotaxis, 
                            offset = hm_offset_7, width = width_gene * ncol(chemotaxis), color = NULL, 
                            colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_chemotaxis_tree <- get_heatmap_column_position(chemotaxis_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_chemotaxis_tree <- tibble(
  label = category_list[5],
  x = column_chemotaxis_tree$x[1] - category_title_pos,
  y = 0)
chemotaxis_tree <- chemotaxis_tree +
  geom_text(data = column_chemotaxis_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_chemotaxis_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### methyl-accepting chemotaxis proteins ====
hm_offset_8 = hm_offset_7 + (category_title_pos * (ncol(chemotaxis) + 1))
methyl_chemotaxis_tree <- gheatmap(p = chemotaxis_tree, data = methyl_chemotaxis, 
                                   offset = hm_offset_8, width = width_gene * ncol(methyl_chemotaxis), color = NULL, 
                                   colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_methyl_chemotaxis_tree <- get_heatmap_column_position(methyl_chemotaxis_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_methyl_chemotaxis_tree <- tibble(
  label = category_list[6],
  x = column_methyl_chemotaxis_tree$x[1] - category_title_pos,
  y = 0)
methyl_chemotaxis_tree <- methyl_chemotaxis_tree +
  geom_text(data = column_methyl_chemotaxis_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_methyl_chemotaxis_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### Fimbriae operons ====
hm_offset_9 = hm_offset_8 + (category_title_pos * (ncol(methyl_chemotaxis) + 1))
Fimbriae_tree <- gheatmap(p = methyl_chemotaxis_tree, data = Fimbriae, 
                          offset = hm_offset_9, width = width_gene * ncol(Fimbriae), color = NULL, 
                          colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_Fimbriae_tree <- get_heatmap_column_position(Fimbriae_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_Fimbriae_tree <- tibble(
  label = category_list[7],
  x = column_Fimbriae_tree$x[1] - category_title_pos,
  y = 0)
Fimbriae_tree <- Fimbriae_tree +
  geom_text(data = column_Fimbriae_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_Fimbriae_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### putative type IV pili-related genes ====
hm_offset_10 = hm_offset_9 + (category_title_pos * (ncol(Fimbriae) + 1))
TypeIV_pili_tree <- gheatmap(p = Fimbriae_tree, data = TypeIV_pili, 
                             offset = hm_offset_10, width = width_gene * ncol(TypeIV_pili), color = NULL, 
                             colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_TypeIV_pili_tree <- get_heatmap_column_position(TypeIV_pili_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_TypeIV_pili_tree <- tibble(
  label = "Type IV pili-related",
  x = column_TypeIV_pili_tree$x[1] - category_title_pos,
  y = 0)
TypeIV_pili_tree <- TypeIV_pili_tree +
  geom_text(data = column_TypeIV_pili_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_TypeIV_pili_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### type II secretion pathway related protein ====
hm_offset_11 = hm_offset_10 + (category_title_pos * (ncol(TypeIV_pili) + 1))
TypeII_tree <- gheatmap(p = TypeIV_pili_tree, data = TypeII, 
                        offset = hm_offset_11, width = width_gene * ncol(TypeII), color = NULL, 
                        colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_TypeII_tree <- get_heatmap_column_position(TypeII_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_TypeII_tree <- tibble(
  label = "T2SS-related",
  x = column_TypeII_tree$x[1] - category_title_pos,
  y = 0)
TypeII_tree <- TypeII_tree +
  geom_text(data = column_TypeII_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_TypeII_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### hemolysin/hemagglutin-like two-partner Type V secretion systems ====
hm_offset_12 = hm_offset_11 + (category_title_pos * (ncol(TypeII) + 1))
hemolysin_tree <- gheatmap(p = TypeII_tree, data = hemolysin, 
                           offset = hm_offset_12, width = width_gene * ncol(hemolysin), color = NULL, 
                           colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_hemolysin_tree <- get_heatmap_column_position(hemolysin_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_hemolysin_tree <- tibble(
  label = "Two-partner T5SS-related",
  x = column_hemolysin_tree$x[1] - category_title_pos,
  y = 0)
hemolysin_tree <- hemolysin_tree +
  geom_text(data = column_hemolysin_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_hemolysin_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### ent operon and entA- and entD-homologs ====
hm_offset_13 = hm_offset_12 + (category_title_pos * (ncol(hemolysin) + 1))
ent_tree <- gheatmap(p = hemolysin_tree, data = ent, 
                     offset = hm_offset_13, width = width_gene * ncol(ent), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_ent_tree <- get_heatmap_column_position(ent_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = c("ent operon (entF)", "entA-homolog", "entD-homolog")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_ent_tree <- tibble(
  label = "Iron uptake systems",
  x = column_ent_tree$x[1] - category_title_pos,
  y = 0)
ent_tree <- ent_tree +
  geom_text(data = column_ent_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_ent_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### siderophore synthesis ====
hm_offset_14 = hm_offset_13 + (category_title_pos * (ncol(ent)))
siderophore_tree <- gheatmap(p = ent_tree, data = siderophore, 
                             offset = hm_offset_14, width = width_gene * ncol(siderophore), color = NULL, 
                             colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_siderophore_tree <- get_heatmap_column_position(siderophore_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = c("siderophore synthesis operon (fes_2)")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
siderophore_tree <- siderophore_tree +
  geom_text(data = column_siderophore_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### fhu operon ====
hm_offset_15 = hm_offset_14 + (category_title_pos * (ncol(siderophore)))
fhu_tree <- gheatmap(p = siderophore_tree, data = fhu, 
                     offset = hm_offset_15, width = width_gene * ncol(fhu), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_fhu_tree <- get_heatmap_column_position(fhu_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(category_list[13], " (", label, ")", sep = "")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
fhu_tree <- fhu_tree +
  geom_text(data = column_fhu_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### fep operon ====
hm_offset_16 = hm_offset_15 + (category_title_pos * (ncol(fhu)))
fep_tree <- gheatmap(p = fhu_tree, data = fep, 
                     offset = hm_offset_16, width = width_gene * ncol(fep), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_fep_tree <- get_heatmap_column_position(fep_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(category_list[14], " (", label, ")", sep = "")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
fep_tree <- fep_tree +
  geom_text(data = column_fep_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### fec operon ====
hm_offset_17 = hm_offset_16 + (category_title_pos * (ncol(fep)))
fec_tree <- gheatmap(p = fep_tree, data = fec, 
                     offset = hm_offset_17, width = width_gene * ncol(fec), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_fec_tree <- get_heatmap_column_position(fec_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(category_list[15], " (", label, ")", sep = "")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
fec_tree <- fec_tree +
  geom_text(data = column_fec_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### has operon ====
hm_offset_18 = hm_offset_17 + (category_title_pos * (ncol(fec)))
has_tree <- gheatmap(p = fec_tree, data = has, 
                     offset = hm_offset_18, width = width_gene * ncol(has), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_has_tree <- get_heatmap_column_position(has_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(category_list[16], " (", label, ")", sep = "")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
has_tree <- has_tree +
  geom_text(data = column_has_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### hem operon ====
hm_offset_19 = hm_offset_18 + (category_title_pos * (ncol(has)))
hem_tree <- gheatmap(p = has_tree, data = hem, 
                     offset = hm_offset_19, width = width_gene * ncol(hem), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_hem_tree <- get_heatmap_column_position(hem_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(category_list[17], " (", label, ")", sep = "")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
hem_tree <- hem_tree +
  geom_text(data = column_hem_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### hms operon ====
hm_offset_20 = hm_offset_19 + (category_title_pos * (ncol(hem)))
hms_tree <- gheatmap(p = hem_tree, data = hms, 
                     offset = hm_offset_20, width = width_gene * ncol(hms), color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_hms_tree <- get_heatmap_column_position(hms_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(category_list[18], " (", label, ")", sep = "")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
hms_tree <- hms_tree +
  geom_text(data = column_hms_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### Other TonB-dependent receptors/transporters ====
hm_offset_21 = hm_offset_20 + (category_title_pos * (ncol(hms) + 1))
TonB_tree <- gheatmap(p = hms_tree, data = TonB, 
                      offset = hm_offset_21, width = width_gene * ncol(TonB), color = NULL, 
                      colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_TonB_tree <- get_heatmap_column_position(TonB_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_TonB_tree <- tibble(
  label = category_list[19],
  x = column_TonB_tree$x[1] - category_title_pos,
  y = 0)
TonB_tree <- TonB_tree +
  geom_text(data = column_TonB_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_TonB_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### Serrawettin W2 ====
hm_offset_22 = hm_offset_21 + (category_title_pos * (ncol(TonB) + 1))
serrawettin_tree <- gheatmap(p = TonB_tree, data = serrawettin, 
                      offset = hm_offset_22, width = width_gene * ncol(serrawettin), color = NULL, 
                      colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_serrawettin_tree <- get_heatmap_column_position(serrawettin_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_serrawettin_tree <- tibble(
  label = category_list[20],
  x = column_serrawettin_tree$x[1] - category_title_pos,
  y = 0)
serrawettin_tree <- serrawettin_tree +
  geom_text(data = column_serrawettin_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_serrawettin_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

#### virulence-attenuated mutants ====
hm_offset_23 = hm_offset_22 + (category_title_pos * (ncol(serrawettin) + 1))
vir_mutants_tree <- gheatmap(p = serrawettin_tree, data = vir_mutants, 
                             offset = hm_offset_23, width = width_gene * ncol(vir_mutants), color = NULL, 
                             colnames = F) +
  scale_fill_manual(name = "Virulence associated gene", values = vir_color, 
                    breaks = vir_type, guide = "none") + 
  new_scale_fill()
column_vir_mutants_tree <- get_heatmap_column_position(vir_mutants_tree, by = "bottom") %>% 
  left_join(select(vir, Gene, cluster_type), by = c("label" = "Gene")) %>% 
  left_join(vir_ID, by = c("label" = "Gene")) %>% 
  mutate(label = str_c(SM39, "; ", Db11, " (", label, ")"))
category_title_vir_mutants_tree <- tibble(
  label = "virulence-attenuated mutations in Db11",
  x = column_vir_mutants_tree$x[1] - category_title_pos,
  y = 0)
vir_mutants_tree <- vir_mutants_tree +
  geom_text(data = column_vir_mutants_tree, mapping = aes(x = x, y = y, label = label, color = cluster_type),
            hjust = 0, nudge_y = 10, size = fontsize_AMR) +
  geom_text(data = category_title_vir_mutants_tree, mapping = aes(x = x, y = y, label = label),
            hjust = 0, nudge_y = 2, size = fontsize_AMR, fontface = "bold") +
  scale_color_manual(values = vir_color, breaks = vir_type, guide = "none") +
  new_scale_color()

# output ====
letter_plot <- vir_mutants_tree +
  scale_y_reverse() +
  theme(plot.margin = unit(x = c(0, 5, 0, 0), units = "line"), 
        legend.title = element_markdown(size = fontsize_theme * 1.2, face = "bold"),
        legend.text = element_text(size = fontsize_theme),
        legend.key.size = unit(x = .4, units = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.justification = c(0, 1),
        legend.position = c(.05, .95), 
        legend.box = "horizontal" 
  )

ggsave(filename = str_c(file_path, "SF08a_virulence_gene_distribution.pdf"), plot = letter_plot, 
       device = "pdf", height = a4_height, width = a4_width, units = "mm", dpi = 1000) 

# adjust legend title with Illustrator