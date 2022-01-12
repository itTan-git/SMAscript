# visualize phylogenic tree

# library ====
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggtext)

# graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# data import ====
ono_tree <- read.tree(file = "~/git/Ono_data/DATA/pangenome_analysis/775strains_tree_20200227.nwk") 
sample_info <- read_tsv(file = "~/git/Ono_data/DATA/final_dataset/sample_info.tsv") 

# tree part ====
## change tree label to sample info ====
otu <- ono_tree$tip.label 
otu_sep <- stringr::str_split(otu, "_", simplify = T)
ono_tree$tip.label <- otu_sep[,1]

## get tree order ====
gg_onotree <- ggtree(tr = ono_tree, ladderize = T, size = .2,
                     layout = "rectangular") +
  scale_y_reverse() 
rgg_onotree <- ggtree::rotate(gg_onotree,1318) %>%  ggtree::rotate(1319)

## make label for rep strain ====
rep_data <- sample_info %>% 
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
complete_data <- sample_info %>% 
  select(sample_ID = `Strain ID in this paper`, Level) %>% 
  inner_join(x = rgg_onotree$data, y = ., by = c("label" = "sample_ID")) %>% 
  mutate(comp_y = ifelse(test = Level == "complete", y, NA)) %>% 
  select(label, comp_y) 
complete_tree <- rep_tree %<+% complete_data 

## draw tree ====
main_tree <- complete_tree +   
  geom_tiplab(mapping = aes(y = rep_y, label = NA), 
              align = T, linetype = 5, linesize = .1, offset = .003,
              colour = "#ff4b00", na.rm = T) + 
  geom_tiplab(mapping = aes(y = rep_y, label = rep_label),  
              align = T, size = fontsize, parse = TRUE,
              hjust = 0, linetype = NULL, offset = .001,
              geom = "text", colour = "#000000", na.rm = TRUE) + 
  geom_tiplab(mapping = aes(y = norep_y, label = NA), geom = "label", align = TRUE,
              linetype = 5, linesize = .05, colour = "#000000", na.rm = TRUE) + 
  geom_tippoint(mapping = aes(x = max(rgg_onotree$data$x) + .0015, y = comp_y), 
                size = .5, shape = 23, fill = "#005aff", 
                colour = "#000000", stroke = .1, na.rm = TRUE) + 
  geom_treescale(x = 0, y = -130, offset = 4, fontsize = fontsize) +
  geom_rootedge(rootedge = .005, size = .2) +
  new_scale_fill()

# region ====
## region category and color ====
region <- c("Africa", "Asia/Oceania", "Europe w/o UK", "UK",
            "Japan", "North America", "South America", "missing")
region_color <- c("#fff100", "#990099", "#f6aa00", "#ff4b00", 
                  "#005aff", "#d8f255", "#03af7a", "#c8c8cb")

## make region dataframe ====
region_info <- sample_info %>% 
  select(Region) %>%  
  mutate_at("Region", as.factor) %>%  
  data.frame() %>%  
  `rownames<-`(sample_info$`Strain ID in this paper`) 

## add region data to tree ====
region_tree <- gheatmap(p = main_tree, data = region_info, 
                        offset = .035, width = .02, color = NULL, 
                        colnames = T, colnames_position = "bottom", hjust = 0,
                        colnames_offset_y = 5, colnames_angle = 90, font.size = fontsize) + 
  scale_fill_manual(name = "Region", values = region_color, breaks = region, 
                    label = c(region[1], "Asia/Oceania w/o Japan", region[3:8]),
                    guide = guide_legend(order = 1)) + 
  new_scale_fill()

# Source info ====
## source colour ====
source_type <- c("clinical", "hospital associated", "nonclinical")
source_color <- c("#ff4b00", "#ff8082", "#005aff")

## make source dataframe ====
source_info <- sample_info %>% 
  select(Source = `clinical or hospital associated or nonclinical`) %>% 
  mutate(Source = factor(ifelse(Source == "n.a.", "nonclinical", Source))) %>%  
  data.frame() %>% 
  `rownames<-`(sample_info$`Strain ID in this paper`) 

## add source data to tree ====
source_tree <- gheatmap(p = region_tree, data = source_info, 
                        offset = .047, width = .02, color = NULL, 
                        colnames = T, colnames_position = "bottom", hjust = 0,
                        colnames_offset_y = 5, colnames_angle = 90, font.size = fontsize) + 
  scale_fill_manual(name = "Source", values = source_color, breaks = source_type,
                    label = c("clinical", "hospial associated", "others"), 
                    guide = guide_legend(order = 2)) +  
  new_scale_fill()

# clade ====
## clade color ====
ANIbclade <- as.character(c(1:14))
ANIbcolor <- c("#dc143c", "#ff8c00", "#ff00ff", "#ffca80",	
               "#000000", "#ffb6c1", "#ffff00", "#005aff",	
               "#00bfff", "#00ffff", "#bfe4ff", "#03af7a", 
               "#00ff7f", "#a9a9a9")

## make clade dataframe ====
clade_info <- sample_info %>% 
  select(Clade = `ANI clade`) %>%  
  mutate(Clade = factor(Clade)) %>%
  data.frame() %>%  
  `rownames<-`(sample_info$`Strain ID in this paper`) 

## add clade data to tree ====
clade_tree <- gheatmap(p = source_tree, data = clade_info, 
                       offset = .059, width = .02, color = NULL, 
                       colnames = T, colnames_position = "bottom", hjust = 0,
                       colnames_offset_y = 5, colnames_angle = 90, font.size = fontsize) + 
  scale_fill_manual(name = "Clade", values = ANIbcolor, breaks = ANIbclade, guide = guide_legend(order = 3, ncol = 2)) + 
  new_scale_fill()

# pig gene cluster ====
## pig gene cluster color ====
pig_type <- c("+", "-")
pig_color <- c("#ff4b00", "#c8c8cb")

## make pig gene cluster dataframe ====
pig_info <- sample_info %>% 
  select(pig_cluster = `pig gene cluster`) %>%  
  mutate(pig_cluster = factor(ifelse(is.na(pig_cluster), "-", "+"))) %>%  
  select("*pig* cluster" = pig_cluster) %>% 
  data.frame(check.names = F) %>%  
  `rownames<-`(sample_info$`Strain ID in this paper`) 

## pig cluster gain/loss events ====
node_list <- clade_tree$data %>% 
  filter(!isTip) %>% 
  select(node)

node_pig_info <- tibble(node = node_list$node)
### each node pig info
for (N in 1:length(node_list$node)){
  temp_clade <- clade_tree %>% 
    groupClade(.node = node_list$node[N]) 
  temp_clade_tip <- temp_clade$data %>% 
    filter(group == 1 & isTip == T) %>% 
    pull(label) 
  temp_clade_pig <- pig_info %>% 
    rownames_to_column(var = "ID") %>% 
    mutate(pig_cluster = as.character(`*pig* cluster`)) %>% 
    filter(ID %in% temp_clade_tip) %>% 
    arrange(pig_cluster) %>% 
    pull(pig_cluster) %>% 
    unique() %>% 
    str_c(collapse = "")
  node_pig_info$pig[N] <- temp_clade_pig
}

pig_event <- node_pig_info %>% 
  filter(pig == "-+") %>% 
  pull(node)

pig_node_position <- clade_tree$data %>% 
  filter(node %in% pig_event)

## add pig gene cluster data to tree ====
pig_tree <- gheatmap(p = clade_tree, data = pig_info, 
                     offset = .071, width = .02, color = NULL, 
                     colnames = F) +
  scale_fill_manual(name = "*pig* cluster", values = pig_color, breaks = pig_type, guide = guide_legend(order = 4)) + 
  geom_point(data = pig_node_position, mapping = aes(x = x, y = y), color = "#ff4b00", size = .5, shape = 19, fill = "#ff4b00", stroke = .2) +
  new_scale_fill()
column_pig_tree <- get_heatmap_column_position(pig_tree)
pig_tree <- pig_tree +
  geom_richtext(data = column_pig_tree, mapping = aes(x = x, y = y, label = label),
                angle = 90, hjust = 0, size = fontsize, fill = NA, label.colour = NA)


# genome size ====
## genome size dataframe ====
genomeSize_info <- sample_info %>% 
  mutate(sizeMb = `Total length (bp)`/1000000) %>% 
  select(sample_ID = `Strain ID in this paper`, size = sizeMb) 

## add genome size data to tree ====
genomeSize_tree <- pig_tree + 
  geom_fruit(data = genomeSize_info, 
             geom = geom_point,
             mapping = aes(y = sample_ID, x = size),
             offset = .21, pwidth = .15,
             grid.params = list(vline = TRUE, size = .02),
             axis.params = list(axis = "x", text.angle = 90, text.size = fontsize, hjust = 0, line.size = .03, line.color = "#ffffff"),
             shape = 16, size = .28, color = "#ff4b00") + 
  annotate(geom = "text", x = .6, y = -25, label = "Genome size\n(Mbp)", size = 2.5, vjust = 0) 

# GC content ====
## GC content dataframe ====
GC_info <- sample_info %>% 
  select(sample_ID = `Strain ID in this paper`, GC = `GC (%)`) 

## add GC content data to tree ====
GC_tree <- genomeSize_tree + 
  geom_fruit(data = GC_info, 
             geom = geom_point,
             mapping = aes(y = sample_ID, x = GC),
             offset = .03, pwidth = .15,
             grid.params = list(vline = TRUE, size = .02),
             axis.params = list(axis = "x", text.angle = 90, text.size = fontsize, hjust = 0, line.size = .03, line.color = "#ffffff"),
             shape = 16, size = .28, color = "#005aff") + 
  annotate(geom = "text", x = .685, y = -25, label = "GC content\n(%)", size = 2.5, vjust = 0) 

# output ====
letter_plot <- GC_tree + 
  scale_y_reverse() +
  theme(legend.title = element_markdown(size = fontsize_theme * 1.5, face = "bold"),
        legend.text = element_text(size = fontsize_theme),
        legend.key.size = unit(x = .4, units = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.justification = c(0, 0),
        legend.position = c(.05, .05))

ggsave(filename = "F04a_phylogenic_tree.pdf", plot = letter_plot, 
       device = "pdf", width = fig_width * 2, height = fig_height, units = "mm", dpi = 1000,
       path = "RESULTS/pangenome_analysis/") 

# arrow show pig gene gain/loss event is added with Adobe Illustrator based red dot in tree.
