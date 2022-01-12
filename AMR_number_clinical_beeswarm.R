# script option ====
## library import ====
library(tidyverse)
library(data.table)
library(ggbeeswarm)

## graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# data import ====
## sample info ====
strain <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper",
         cluster = "SNP10 cluster",
         picked = "Representative strain in each SNP10 cluster",
         source = "clinical or hospital associated or nonclinical") %>% 
  mutate(source = ifelse(source == "hospital associated", "clinical", 
                         ifelse(source == "n.a.", "nonclinical", 
                                source)))

## AMR number ====
AMR <- read_tsv("DATA/final_dataset/AMR_info.tsv") %>% 
  select(-c(class, category, for_MDR_dicision)) %>%  
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  mutate_at(colnames(.)[-1], function(x){as.numeric(ifelse(x == "-", "0", "1"))}) %>% 
  left_join(x = strain, y = .) %>% 
  select(-c(sample_ID, picked, source)) %>% 
  pivot_longer(cols = -cluster, names_to = "gene", values_to = "PA") %>% 
  group_by(cluster, gene) %>% 
  summarise(P_cluster_ratio = sum(PA)/n()) %>% 
  ungroup(gene) %>% 
  summarise(AMR_number = sum(P_cluster_ratio))
AMR_num <- strain %>% 
  filter(picked == "Rep") %>% 
  left_join(AMR, by = "cluster") 

# draw graph ====
source_AMRnum <- ggplot(data = AMR_num, mapping = aes(x = source, y = AMR_number, fill = source)) +
  geom_quasirandom(size = 1, shape = 21, colour = "#000000", stroke = .2) + 
  scale_fill_manual(breaks = c("clinical", "nonclinical"), 
                    values = c("#ff4b00", "#005aff")) + 
  scale_x_discrete(name = "", label = c("clinical/hospital\nenvironment", "others")) + 
  scale_y_continuous(name = "AMR gene number", breaks = c(seq(0, max(AMR_num$AMR_number), 2),26)) + 
  guides(fill = "none") + 
  theme_bw() +
  theme(axis.title = element_text(size = fontsize_theme * 1.2, colour = "#000000"), 
        axis.text = element_text(size = fontsize_theme * 1.2, colour = "#000000"))

# output ====
ggsave(filename = "RESULTS/clinical/SF07_clinical_AMR_plot.pdf", plot = source_AMRnum,
       device = "pdf", width = fig_width, height = fig_height / 2, units = "mm")
