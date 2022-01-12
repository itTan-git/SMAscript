# integrase search from pangenome

# library ====
library(tidyverse)
library(ggbeeswarm)

# file path ====
file_path <- "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/prophage_integrase/"


## graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# data import
pangenome <- read_csv(file = "DATA/pangenome_analysis/gene_presence_absence.csv")

## search "integrase" word from pangenome
integrase <- pangenome %>% 
  filter(str_detect(string = Annotation, pattern = "integrase")) %>% 
  write_tsv(file = str_c(file_path, "integrase_pangenome.tsv")) %>% 
  select(-c("Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", 
           "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", 
           "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc")) %>% 
  pivot_longer(cols = -Gene, names_to = "strain", values_to = "pre_ab") %>% 
  mutate(strain = str_split(string = strain, pattern = "_", simplify = TRUE)[,1],
         pre_ab = case_when(is.na(pre_ab) ~ 0, TRUE ~ 1)) %>% 
  group_by(strain) %>% 
  summarise(int_num = sum(pre_ab)) %>% 
  write_tsv(file = str_c(file_path, "integrase_pangenome_num.tsv"))
