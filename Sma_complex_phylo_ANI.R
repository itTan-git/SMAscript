# phylogenic tree of representative strain of clade and Serratia spp. type strain
# and ANI analysis

# library
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggtext)

# graphic parameter
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# import tree data ====
SMc_tree <- read.tree(file = "RESULTS/Smacomplex_strains/RAxML_bipartitions.SF02_strain_SmaTS_reroot.nwk") 

# tip label data ====
tip_label <- 
c("Sma MSU97 (DL526)",  "Sma SM03 (DL542)",  "Sma QE07 (SE191)", 
  "Sma QE08 (SE192)",  "Sma QE10 (SE194)", "Sma E46 (SE084)",
  "Ssu YD25<sup>T</sup>",  "Sma QE30 (SE213)",  "Sma W2.3 (DL601)",
  "Sma PH1a (DL622)",  "Sma H1q (DL623)",  "Sma Db11 (DL020)",
  "Sma QE03 (SE187)",  "Sur Lr5/4 (DL659)",  "Sur CCUG:50595<sup>T</sup>",
  "Sma 12ES (DL474)",  "Sma_sak K27 (DL654)", "Sne DSM21420<sup>T</sup> (DL658)", 
  "Sma NCTC10211<sup>T</sup> (DL661)", "Sma_sak DSM17174<sup>T</sup>", "Sma SM39 (DL536)",
  "Sma A1 (DL241)",  "Sne CRK0003 (DL655)",  "Sma CAVp305 (DL465)",
  "Sma 332 (DL339)",  "Sma E29 (SE068)",  "Sma PWN146 (DL006)",
  "Sma E19 (SE061)",  "Sma 092713_C_TSB (DL545)",  "*S. ficaria*<sup>T</sup>")

# display tip label ====
SMc_tree$tip.label <- tip_label 
gg_SMc_tree_ori <- ggtree(tr = SMc_tree, ladderize = T, size = .2,
                          layout = "rectangular") +
  geom_tiplab(align = TRUE, mapping = aes(label = NA), linetype = 5, linesize = .05)

# clade data
clade_data <- read_tsv(file = "RESULTS/Smacomplex_strains/clade_data.tsv")

TS_data <- gg_SMc_tree_ori$data %>% 
  left_join(clade_data, by = c("label" = "tip_label"))
gg_SMc_tree_ori$data <- TS_data  

## ANI clade colour parameter ====
ANI_score <- c("< 93", "93-94", "94-95", "95-96", "96-97", "97-98", "98-100")
ANI_color <- c("#ffffff", "#7ff8c1", "#cafcea", "#fdf0ce", "#f6d072", "#f1996f", "#EB5929")

## data for ANI heatmap ====
ANI_table <- read_tsv(file = "RESULTS/Smacomplex_strains/ANIb_percentage_identity_mod.tsv")
SMc_ANI <- ANI_table %>% 
  select(-sample) %>% 
  data.frame(row.names = ANI_table$sample, check.names = FALSE) %>% 
  mutate_all(function(x){case_when(x < 0.93 ~ ANI_score[1],
                                   0.93 <= x & x < 0.94 ~ ANI_score[2],
                                   0.94 <= x & x < 0.95 ~ ANI_score[3],
                                   0.95 <= x & x < 0.96 ~ ANI_score[4],
                                   0.96 <= x & x < 0.97 ~ ANI_score[5],
                                   0.97 <= x & x < 0.98 ~ ANI_score[6],
                                   TRUE ~ ANI_score[7])})

# make data for display tiplabel and clade info ====
gg_SMc_tree <- rotate(gg_SMc_tree_ori, 31) %>% rotate(45) %>% rotate(51) %>% 
  rotate(52) %>% rotate(53) %>% rotate(59) %>% rotate(56)

display_data <- gg_SMc_tree$data %>% 
  filter(isTip) %>% 
  select(label, x, y, ANI_clade) %>% 
  mutate(tip_x = max(x),
         clade_x = tip_x + .05)

gg_SMc_tree_mod <- gg_SMc_tree +
  geom_richtext(data = display_data, mapping = aes(x = tip_x, y = y, label = label),
                size = fontsize, hjust = 0, fill = NA, label.colour = NA) +
  geom_treescale(x = 0, y = 20, offset = 0.5, fontsize = fontsize)

## add ANI matrix ====
ANI_tree <- gheatmap(p = gg_SMc_tree_mod, data = SMc_ANI, 
                     offset = .11, width = 1,
                     color = "#000000", colnames = F) +
  scale_fill_manual(name = "ANI", values = ANI_color, 
                    breaks = ANI_score, guide = guide_legend(nrow = 1, direction = "vertical", byrow = TRUE)) +
  ylim(c(0,37)) 
column_ANI_tree <- get_heatmap_column_position(ANI_tree, by = "top")
display_data$clade_x <- min(column_ANI_tree$x) - .015
letter_plot <- ANI_tree +
  geom_text(data = display_data, mapping = aes(x = clade_x, y = y, label = ANI_clade),
            size = fontsize, na.rm = TRUE) +
  geom_richtext(data = column_ANI_tree, mapping = aes(x = x, y = y, label = label),
                size = fontsize * .8, angle = 90, fill = NA, label.colour = NA, hjust = 0) +
  annotate(geom = "text", x = display_data$clade_x[1], y = 31, label = "Clade", size = fontsize) +
  theme(plot.margin = unit(c(0,0,1.5,0), units = "line") ,
        legend.position = c(.5,0),
        legend.title = element_text(size = fontsize_theme * 1.5, face = "bold"),
        legend.text = element_text(size = fontsize_theme),
        legend.key.size = unit(x = .4, units = "cm"))

ggsave(filename = "SF02_SMcomplex_phylo.pdf", plot = letter_plot,
       device = "pdf", height = fig_height / 1.5, width = fig_width * 2, units = "mm", 
       path = "RESULTS/Smacomplex_strains/", dpi = 1000) 
