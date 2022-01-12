# distribution of the genes specifically conserved in clade 1 and clade 2

# script option ====
## PATH ====
file_path <- "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/core_gene/"

## library import ====
library(tidyverse)
library(data.table)
library(ggtree)
library(ggnewscale)

## graphic parameter ====
fontsize = 1.5
fontsize_AMR = 1.2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
a4_width = 210
a4_height = 297

# data import ====
## sample info ====
strain_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper",
         clade = "ANI clade",
         source = "clinical or hospital associated or nonclinical")

## pangenome concentrated gene data import ====
pan_clade1and2 <- read_tsv(file = str_c(file_path, "scoary_1and2_vs_others.tsv")) %>% 
  arrange(`p value*`) %>% 
  head(100) %>% 
  select(Gene, SM39_order, SM39, Db11)

## pangenome gene distribution data ====
pangenome_preab <- read_csv(file = "DATA/pangenome_analysis/gene_presence_absence.csv")
pan_clade1and2_data <- pangenome_preab %>% 
  select(-c("Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", 
            "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", 
            "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc")) %>% 
  right_join(select(pan_clade1and2, Gene, SM39_order), by = "Gene") %>% 
  arrange(SM39_order) %>% 
  pivot_longer(cols = -Gene, names_to = "strain", values_to = "locusID") %>% 
  mutate(strain = str_split(string = strain, pattern = "_", simplify = TRUE)[,1]) %>% 
  filter(strain %in% strain_info$sample_ID) %>% 
  mutate(pre_ab = case_when(is.na(locusID) ~ "-", TRUE ~ "+")) %>% 
  select(-locusID) %>% 
  pivot_wider(names_from = "Gene", values_from = "pre_ab") %>% 
  data.frame(row.names = .$strain, check.names = FALSE) %>% 
  select(-strain)
  
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

### add clade heatmap ====
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

## concentrated gene ====
### color ====
pre_ab_type <- c("+", "-")
pre_ab_color <- c("#005aff", "#c8c8cb")

### config ====
width_gene = .02 * 1.5

#### genes enriched in clade 1 and clade 2 ====
hm_offset_3 = hm_offset_2 + (0.009333835 * 1.5)
enriched_gene_tree <- gheatmap(p = source_tree, data = pan_clade1and2_data, 
                                   offset = hm_offset_3, width = width_gene * ncol(pan_clade1and2_data), 
                                   color = NULL, colnames = F) +
  scale_fill_manual(name = "Accessory genes specifically\nconserved in both clades 1 and 2", values = pre_ab_color, 
                    breaks = pre_ab_type, guide = guide_legend(order = 3)) + 
  new_scale_fill()
column_enriched_gene_tree <- get_heatmap_column_position(enriched_gene_tree, by = "bottom") %>% 
  left_join(pan_clade1and2, by = c("label" = "Gene")) %>% 
  mutate(SM39_order_ID = str_split(string = SM39_order, pattern = "_", simplify = TRUE)[,2],
         Db11_ID = case_when(Db11 == "-" ~ "-", TRUE ~ str_split(string = Db11, pattern = "_", simplify = TRUE)[,2])) %>% 
  mutate(SM39_order_ID = case_when(str_starts(string = SM39, pattern = "between") ~ str_c(SM39_order_ID, "*", sep = ""),
                                is.na(SM39_order) ~ "-",
                                TRUE ~ SM39_order_ID),
         gene = str_c("(", label, ")")) 
column_description <- tibble(
  label = c("ID in\nSM39", ":", "ID in\nDb11"),
  x = tail(column_enriched_gene_tree$x, 1) + (1.5 * diff(column_enriched_gene_tree$x)[1]),
  y = c(-12, -23.5, -34)
)
enriched_gene_tree <- enriched_gene_tree +
  geom_text(data = column_enriched_gene_tree, mapping = aes(x = x, y = y, label = SM39_order_ID),
            hjust = 0, nudge_y = 5, size = fontsize) +
  geom_text(data = column_enriched_gene_tree, mapping = aes(x = x, y = y - 23.5), label = ":",
            hjust = .5, size = fontsize) +
  geom_text(data = column_enriched_gene_tree, mapping = aes(x = x, y = y - 27, label = Db11_ID),
            hjust = 0, size = fontsize) +
  geom_text(data = column_enriched_gene_tree, mapping = aes(x = x, y = y - 44, label = gene),
            hjust = 0, size = fontsize) +
  geom_text(data = column_description, mapping = aes(x = x, y = y, label = label),
            hjust = .5, size = fontsize, vjust = .5, lineheight = .8) +
  new_scale_color()

# output ====
letter_plot <- enriched_gene_tree +
  scale_y_reverse() +
  theme(plot.margin = unit(c(0, .5, 0, 0), units = "lines"),
        legend.title = element_text(size = fontsize_theme * 1.2, face = "bold"),
        legend.text = element_text(size = fontsize_theme),
        legend.key.size = unit(x = .4, units = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.justification = c(0, 1),
        legend.position = c(.05, .95), 
        legend.box = "horizontal" 
  )

ggsave(filename = str_c(file_path, "SF09a_clade1and2_conserved_gene.pdf"), plot = letter_plot, 
       device = "pdf", height = a4_height, width = a4_width, units = "mm", dpi = 1000) # A4用のサイズで保存 

# adjust legend title with Illustrator